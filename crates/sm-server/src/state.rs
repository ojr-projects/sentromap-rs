//! Server state: the index, annotations, calibrated cost model, query jobs and the result cache
//! (design §12.1).

use std::collections::{BTreeMap, HashMap};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Arc, Mutex};
use std::time::Instant;

use serde::Serialize;
use sm_agg::{OrderMode, ResultSet};
use sm_annot::Annotations;
use sm_index::Index;
use sm_search::{CostModel, Plan};
use tokio::sync::Semaphore;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Status {
    Queued,
    Searching,
    Ordering,
    Ready,
    Failed,
    Cancelled,
}

pub struct JobState {
    pub status: Status,
    pub result: Option<Arc<ResultSet>>,
    pub error: Option<String>,
    pub timings_ms: BTreeMap<&'static str, f64>,
}

pub struct Job {
    pub id: u64,
    /// The query as the user gave it.
    pub kmer: u64,
    /// True when the query is the reverse complement of its canonical form: results are
    /// computed for the canonical form and mirrored on output.
    pub flipped: bool,
    pub max_n: u8,
    pub order: OrderMode,
    pub plan: Plan,
    pub cached: bool,
    pub created: Instant,
    pub cancelled: AtomicBool,
    pub state: Mutex<JobState>,
}

impl Job {
    pub fn set_status(&self, s: Status) {
        self.state.lock().unwrap().status = s;
    }

    pub fn time(&self, stage: &'static str, ms: f64) {
        self.state.lock().unwrap().timings_ms.insert(stage, ms);
    }

    pub fn result(&self) -> Option<Arc<ResultSet>> {
        self.state.lock().unwrap().result.clone()
    }

    pub fn is_cancelled(&self) -> bool {
        self.cancelled.load(Ordering::Relaxed)
    }

    /// Orient a canonical-frame k-mer to the user's query orientation.
    pub fn orient(&self, k: u64) -> u64 {
        if self.flipped { sm_core::revcomp(k) } else { k }
    }
}

struct CacheEntry {
    rs: Arc<ResultSet>,
    last_used: u64,
}

/// LRU over result sets keyed by `(canonical query, max n)`, with a byte budget.
pub struct Cache {
    entries: HashMap<(u64, u8), CacheEntry>,
    clock: u64,
    pub budget: usize,
}

impl Cache {
    pub fn new(budget: usize) -> Self {
        Self { entries: HashMap::new(), clock: 0, budget }
    }

    /// A cached result for `canonical` whose max n covers `max_n`.
    pub fn get(&mut self, canonical: u64, max_n: u8) -> Option<Arc<ResultSet>> {
        self.clock += 1;
        let clock = self.clock;
        let e = self
            .entries
            .iter_mut()
            .filter(|((q, n), _)| *q == canonical && *n >= max_n)
            .min_by_key(|((_, n), _)| *n)?;
        e.1.last_used = clock;
        Some(e.1.rs.clone())
    }

    pub fn insert(&mut self, canonical: u64, rs: Arc<ResultSet>) {
        self.clock += 1;
        self.entries.insert((canonical, rs.max_n), CacheEntry { rs, last_used: self.clock });
        self.evict();
    }

    pub fn bytes(&self) -> usize {
        self.entries.values().map(|e| e.rs.bytes()).sum()
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    fn evict(&mut self) {
        while self.entries.len() > 1 && self.bytes() > self.budget {
            let oldest = *self.entries.iter().min_by_key(|(_, e)| e.last_used).unwrap().0;
            self.entries.remove(&oldest);
        }
    }
}

pub struct App {
    pub idx: Index,
    pub annots: Annotations,
    pub model: CostModel,
    pub name: String,
    pub jobs: Mutex<HashMap<u64, Arc<Job>>>,
    pub cache: Mutex<Cache>,
    /// Admission control: jobs expected to saturate memory bandwidth run one at a time.
    pub heavy: Semaphore,
    pub next_id: AtomicU64,
    pub started: Instant,
}

/// Jobs estimated under this run concurrently in the fast lane.
pub const FAST_LANE_MS: f64 = 50.0;
/// Jobs kept for polling.
pub const MAX_JOBS: usize = 10_000;

impl App {
    pub fn new(idx: Index, annots: Annotations, model: CostModel, cache_bytes: usize) -> Self {
        let name = if idx.manifest.source.name.is_empty() {
            idx.dir.file_name().map(|s| s.to_string_lossy().into_owned()).unwrap_or_default()
        } else {
            idx.manifest.source.name.clone()
        };
        Self {
            idx,
            annots,
            model,
            name,
            jobs: Mutex::new(HashMap::new()),
            cache: Mutex::new(Cache::new(cache_bytes)),
            heavy: Semaphore::new(1),
            next_id: AtomicU64::new(1),
            started: Instant::now(),
        }
    }

    pub fn job(&self, id: u64) -> Option<Arc<Job>> {
        self.jobs.lock().unwrap().get(&id).cloned()
    }

    pub fn add_job(&self, job: Arc<Job>) {
        let mut jobs = self.jobs.lock().unwrap();
        if jobs.len() >= MAX_JOBS {
            let oldest = jobs.values().min_by_key(|j| j.created).map(|j| j.id);
            if let Some(id) = oldest {
                jobs.remove(&id);
            }
        }
        jobs.insert(job.id, job);
    }

    /// Display name of a contig.
    pub fn display(&self, ci: usize) -> &str {
        self.annots.registry.contig_display.get(ci).map_or(self.idx.contigs.get(ci).name.as_str(), |s| s.as_str())
    }
}

/// Run a query job: search, assemble, order at max n; cache the result.
pub async fn run_job(app: Arc<App>, job: Arc<Job>) {
    let heavy = job.plan.estimate_ms > FAST_LANE_MS;
    let _permit = if heavy { Some(app.heavy.acquire().await.expect("semaphore closed")) } else { None };
    if job.is_cancelled() {
        job.set_status(Status::Cancelled);
        return;
    }
    job.set_status(Status::Searching);
    let queued_ms = job.created.elapsed().as_secs_f64() * 1e3;
    job.time("queued", queued_ms);
    let (a, j) = (app.clone(), job.clone());
    let out = tokio::task::spawn_blocking(move || -> Result<Arc<ResultSet>, String> {
        let canonical = sm_core::canonical(j.kmer);
        let t = Instant::now();
        let (v, _) = sm_search::search(&a.idx, canonical, j.max_n as u32, &a.model, Some(j.plan.engine));
        j.time("search", t.elapsed().as_secs_f64() * 1e3);
        if j.is_cancelled() {
            return Err("cancelled".into());
        }
        let t = Instant::now();
        let rs = Arc::new(ResultSet::new(&a.idx, v));
        j.time("assemble", t.elapsed().as_secs_f64() * 1e3);
        if j.is_cancelled() {
            return Err("cancelled".into());
        }
        j.set_status(Status::Ordering);
        let t = Instant::now();
        rs.order(j.order, j.max_n);
        j.time("order", t.elapsed().as_secs_f64() * 1e3);
        Ok(rs)
    })
    .await;
    let mut st = job.state.lock().unwrap();
    match out {
        Ok(Ok(rs)) => {
            app.cache.lock().unwrap().insert(sm_core::canonical(job.kmer), rs.clone());
            st.result = Some(rs);
            st.status = Status::Ready;
        }
        Ok(Err(_)) if job.is_cancelled() => st.status = Status::Cancelled,
        Ok(Err(e)) => {
            st.status = Status::Failed;
            st.error = Some(e);
        }
        Err(e) => {
            st.status = Status::Failed;
            st.error = Some(format!("job panicked: {e}"));
        }
    }
    st.timings_ms.insert("total", job.created.elapsed().as_secs_f64() * 1e3);
}

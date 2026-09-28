//! HTTP endpoints (design §12.2). Control and small responses are JSON; large arrays (rasters,
//! histograms, track summaries) are little-endian binary with dimensions in headers.

use std::sync::Arc;
use std::sync::atomic::Ordering;
use std::time::Instant;

use axum::Json;
use axum::body::Body;
use axum::extract::{Path, Query, State};
use axum::http::{HeaderMap, HeaderValue, StatusCode, header};
use axum::response::{IntoResponse, Response};
use serde::Deserialize;
use serde_json::{Value, json};
use sm_agg::{OrderMode, OrderView, ResultSet, View};
use sm_core::kmer::{self, K};

use crate::state::{App, Job, JobState, Status, run_job};

pub type AppState = Arc<App>;

/// API error: status plus message, rendered as `{"error": ...}`.
pub struct ApiError(StatusCode, String);

impl IntoResponse for ApiError {
    fn into_response(self) -> Response {
        (self.0, Json(json!({ "error": self.1 }))).into_response()
    }
}

fn bad(msg: impl Into<String>) -> ApiError {
    ApiError(StatusCode::BAD_REQUEST, msg.into())
}

fn not_found(msg: impl Into<String>) -> ApiError {
    ApiError(StatusCode::NOT_FOUND, msg.into())
}

type ApiResult<T> = Result<T, ApiError>;

/// Run CPU-heavy work off the async runtime.
async fn blocking<T: Send + 'static>(f: impl FnOnce() -> T + Send + 'static) -> ApiResult<T> {
    tokio::task::spawn_blocking(f).await.map_err(|e| ApiError(StatusCode::INTERNAL_SERVER_ERROR, e.to_string()))
}

pub const MAX_N_LIMIT: u8 = 17;
const MAX_CELLS: u64 = 16_000_000;
const MAX_WIDTH: u32 = 16_384;

fn genome_len(app: &App) -> u32 {
    app.idx.manifest.genome_len as u32
}

/// Clamp a requested global region to the genome; defaults to all of it.
fn region(app: &App, start: Option<u32>, end: Option<u32>) -> ApiResult<(u32, u32)> {
    let len = genome_len(app);
    let s = start.unwrap_or(0).min(len);
    let e = end.unwrap_or(len).min(len);
    if s >= e {
        return Err(bad("empty region: need start < end"));
    }
    Ok((s, e))
}

fn binary(bytes: Vec<u8>, headers: &[(&'static str, String)]) -> Response {
    let mut h = HeaderMap::new();
    h.insert(header::CONTENT_TYPE, HeaderValue::from_static("application/octet-stream"));
    for (k, v) in headers {
        h.insert(*k, HeaderValue::from_str(v).unwrap());
    }
    (h, bytes).into_response()
}

// ---------------------------------------------------------------------------------------------
// Genome, tracks, features, sequence

pub async fn health(State(app): State<AppState>) -> Json<Value> {
    let cache = app.cache.lock().unwrap();
    Json(json!({
        "ok": true,
        "uptime_s": app.started.elapsed().as_secs_f64(),
        "index": app.idx.dir.display().to_string(),
        "files": app.idx.manifest.files.iter().map(|(k, v)| (k.clone(), json!({"bytes": v.bytes, "crc32": v.crc32}))).collect::<serde_json::Map<_, _>>(),
        "cache": { "entries": cache.len(), "bytes": cache.bytes(), "budget": cache.budget },
        "jobs": app.jobs.lock().unwrap().len(),
        "cost_model": app.model,
    }))
}

pub async fn genome(State(app): State<AppState>) -> Json<Value> {
    let m = &app.idx.manifest;
    let contigs: Vec<Value> = app
        .idx
        .contigs
        .iter()
        .enumerate()
        .map(|(i, c)| json!({ "name": c.name, "display": app.display(i), "offset": c.offset, "length": c.len }))
        .collect();
    Json(json!({
        "name": app.name,
        "k": K,
        "genome_len": m.genome_len,
        "distinct_kmers": m.distinct_kmers,
        "sites": m.sites,
        "prefix_len": m.prefix_len,
        "has_copy_b": m.has_copy_b,
        "max_n_limit": MAX_N_LIMIT,
        "contigs": contigs,
    }))
}

pub async fn tracks(State(app): State<AppState>) -> Json<Value> {
    let r = &app.annots.registry;
    Json(json!({
        "tracks": r.tracks,
        "feature_sets": r.feature_sets.iter().map(|f| json!({
            "name": f.name, "description": f.description, "source": f.source, "count": f.count,
            "kinds": f.kinds, "kind_counts": f.kind_counts,
        })).collect::<Vec<_>>(),
    }))
}

#[derive(Deserialize)]
pub struct RegionQ {
    start: Option<u32>,
    end: Option<u32>,
    width: Option<u32>,
}

/// `(mean, max)` f32 pairs per column.
pub async fn track_summary(
    State(app): State<AppState>,
    Path(name): Path<String>,
    Query(q): Query<RegionQ>,
) -> ApiResult<Response> {
    let (s, e) = region(&app, q.start, q.end)?;
    let w = q.width.unwrap_or(1000).clamp(1, MAX_WIDTH) as usize;
    let a = app.clone();
    let vals = blocking(move || a.annots.track(&name).map(|t| t.summary(s, e, w)))
        .await?
        .ok_or_else(|| not_found("no such track"))?;
    let bytes: Vec<u8> = bytemuck_f32_pairs(&vals);
    Ok(binary(bytes, &[("x-width", w.to_string())]))
}

fn bytemuck_f32_pairs(v: &[[f32; 2]]) -> Vec<u8> {
    v.iter().flat_map(|p| p[0].to_le_bytes().into_iter().chain(p[1].to_le_bytes())).collect()
}

#[derive(Deserialize)]
pub struct FeaturesQ {
    start: Option<u32>,
    end: Option<u32>,
    width: Option<u32>,
    limit: Option<usize>,
    /// Comma-separated kind names to include.
    kinds: Option<String>,
}

/// Features overlapping a region, or a density per column when there are too many.
pub async fn features(
    State(app): State<AppState>,
    Path(name): Path<String>,
    Query(q): Query<FeaturesQ>,
) -> ApiResult<Json<Value>> {
    let (s, e) = region(&app, q.start, q.end)?;
    let limit = q.limit.unwrap_or(2000).min(50_000);
    let w = q.width.unwrap_or(1000).clamp(1, MAX_WIDTH) as usize;
    let a = app.clone();
    blocking(move || {
        let fs = a.annots.feature_set(&name).ok_or_else(|| not_found("no such feature set"))?;
        let kinds: Option<Vec<u16>> = q.kinds.as_ref().map(|k| {
            k.split(',').filter_map(|name| fs.info.kinds.iter().position(|x| x == name).map(|i| i as u16)).collect()
        });
        let mut out = Vec::new();
        let mut total = 0usize;
        for i in fs.overlapping(s, e) {
            let f = fs.get(i);
            if kinds.as_ref().is_some_and(|k| !k.contains(&f.kind)) {
                continue;
            }
            total += 1;
            if out.len() < limit {
                out.push(f);
            }
        }
        if total <= limit {
            let feats: Vec<Value> = out
                .iter()
                .map(|f| json!({ "start": f.start, "end": f.end, "kind": f.kind, "strand": f.strand, "name": f.name }))
                .collect();
            Ok(Json(json!({ "mode": "features", "kinds": fs.info.kinds, "features": feats })))
        } else {
            let d = fs.density(s, e, w, kinds.as_deref());
            Ok(Json(json!({ "mode": "density", "kinds": fs.info.kinds, "total": total, "density": d })))
        }
    })
    .await?
}

#[derive(Deserialize)]
pub struct SeqQ {
    start: u32,
    end: u32,
}

pub async fn sequence(State(app): State<AppState>, Query(q): Query<SeqQ>) -> ApiResult<Response> {
    let (s, e) = region(&app, Some(q.start), Some(q.end))?;
    if e - s > 1_000_000 {
        return Err(bad("sequence requests are limited to 1 Mb"));
    }
    let bytes = app.idx.seq.fetch(s, e);
    Ok(([(header::CONTENT_TYPE, "text/plain")], bytes).into_response())
}

#[derive(Deserialize)]
pub struct KmerQ {
    pos: u32,
}

/// The 31-mer starting at a global position (for click-to-query).
pub async fn kmer_at(State(app): State<AppState>, Query(q): Query<KmerQ>) -> ApiResult<Json<Value>> {
    let (ci, off) = app.idx.contigs.to_local(q.pos).ok_or_else(|| bad("position outside the genome"))?;
    let k = app.idx.kmer_at(ci, off);
    Ok(Json(json!({
        "pos": q.pos,
        "contig": app.idx.contigs.get(ci).name,
        "display": app.display(ci),
        "start": off + 1,
        "kmer": k.map(kmer::to_string),
        "valid": k.is_some(),
    })))
}

// ---------------------------------------------------------------------------------------------
// Queries

#[derive(Deserialize)]
pub struct NewQuery {
    /// 31 bases.
    kmer: Option<String>,
    /// Or a global 0-based position.
    pos: Option<u32>,
    /// Or a contig name (or display name) with a 1-based start.
    contig: Option<String>,
    start: Option<u32>,
    max_n: Option<u8>,
    order: Option<OrderMode>,
}

pub async fn create_query(State(app): State<AppState>, Json(q): Json<NewQuery>) -> ApiResult<Json<Value>> {
    let kmer = if let Some(s) = &q.kmer {
        kmer::parse(s.trim()).map_err(|e| bad(e.to_string()))?
    } else if let Some(p) = q.pos {
        let (ci, off) = app.idx.contigs.to_local(p).ok_or_else(|| bad("position outside the genome"))?;
        app.idx.kmer_at(ci, off).ok_or_else(|| bad("no valid 31-mer at that position"))?
    } else if let (Some(c), Some(s)) = (&q.contig, q.start) {
        let ci = app
            .idx
            .contigs
            .by_name(c)
            .or_else(|| (0..app.idx.contigs.len()).find(|&i| app.display(i) == c))
            .ok_or_else(|| bad(format!("no contig {c:?}")))?;
        app.idx.kmer_at(ci, s.saturating_sub(1)).ok_or_else(|| bad("no valid 31-mer at that position"))?
    } else {
        return Err(bad("give kmer, pos, or contig + start"));
    };
    let max_n = q.max_n.unwrap_or(10);
    if max_n > MAX_N_LIMIT {
        return Err(bad(format!("max_n is limited to {MAX_N_LIMIT}")));
    }
    let order = q.order.unwrap_or(OrderMode::PerN);
    let canonical = sm_core::canonical(kmer);
    let plan = app.model.plan(&app.idx, max_n as u32);
    let cached = app.cache.lock().unwrap().get(canonical, max_n);
    let job = Arc::new(Job {
        id: app.next_id.fetch_add(1, Ordering::Relaxed),
        kmer,
        flipped: kmer != canonical,
        max_n,
        order,
        plan,
        cached: cached.is_some(),
        created: Instant::now(),
        cancelled: Default::default(),
        state: std::sync::Mutex::new(JobState {
            status: if cached.is_some() { Status::Ready } else { Status::Queued },
            result: cached,
            error: None,
            timings_ms: Default::default(),
        }),
    });
    app.add_job(job.clone());
    if !job.cached {
        tokio::spawn(run_job(app.clone(), job.clone()));
    }
    Ok(Json(json!({
        "id": job.id,
        "kmer": kmer::to_string(kmer),
        "max_n": max_n,
        "order": order,
        "engine": plan.engine,
        "estimate_ms": if job.cached { 0.0 } else { plan.estimate_ms },
        "expected_chance_hits": plan.expected_chance_hits,
        "cached": job.cached,
    })))
}

fn job_or_404(app: &App, id: u64) -> ApiResult<Arc<Job>> {
    app.job(id).ok_or_else(|| not_found("no such query (it may have expired)"))
}

fn ready(job: &Job) -> ApiResult<Arc<ResultSet>> {
    let st = job.state.lock().unwrap();
    match (st.status, &st.result) {
        (Status::Ready, Some(rs)) => Ok(rs.clone()),
        (s, _) => Err(ApiError(StatusCode::CONFLICT, format!("query is {s:?}, not ready"))),
    }
}

pub async fn query_status(State(app): State<AppState>, Path(id): Path<u64>) -> ApiResult<Json<Value>> {
    let job = job_or_404(&app, id)?;
    let st = job.state.lock().unwrap();
    let distinct = app.idx.manifest.distinct_kmers;
    let counts: Option<Vec<Value>> = st.result.as_ref().map(|rs| {
        (0..=job.max_n)
            .map(|n| {
                json!({
                    "n": n,
                    "variants": rs.count(n),
                    "sites": rs.sites_within[n as usize],
                    "expected_chance_hits": sm_core::ball::expected_chance_hits(distinct, n as u32),
                })
            })
            .collect()
    });
    Ok(Json(json!({
        "id": job.id,
        "status": st.status,
        "kmer": kmer::to_string(job.kmer),
        "max_n": job.max_n,
        "order": job.order,
        "engine": job.plan.engine,
        "estimate_ms": job.plan.estimate_ms,
        "cached": job.cached,
        "elapsed_ms": job.created.elapsed().as_secs_f64() * 1e3,
        "timings_ms": st.timings_ms,
        "error": st.error,
        "counts": counts,
    })))
}

pub async fn cancel_query(State(app): State<AppState>, Path(id): Path<u64>) -> ApiResult<Json<Value>> {
    let job = job_or_404(&app, id)?;
    job.cancelled.store(true, Ordering::Relaxed);
    let mut st = job.state.lock().unwrap();
    if matches!(st.status, Status::Queued | Status::Searching | Status::Ordering) {
        st.status = Status::Cancelled;
    }
    Ok(Json(json!({ "id": id, "status": st.status })))
}

#[derive(Deserialize)]
pub struct ViewQ {
    n: Option<u8>,
    order: Option<OrderMode>,
    start: Option<u32>,
    end: Option<u32>,
    row0: Option<u32>,
    row1: Option<u32>,
    width: Option<u32>,
    height: Option<u32>,
    limit: Option<usize>,
    offset: Option<usize>,
    min_rows: Option<u32>,
    format: Option<String>,
}

struct Ctx {
    job: Arc<Job>,
    rs: Arc<ResultSet>,
    n: u8,
    mode: OrderMode,
}

fn ctx(app: &App, id: u64, q: &ViewQ) -> ApiResult<Ctx> {
    let job = job_or_404(app, id)?;
    let rs = ready(&job)?;
    let n = q.n.unwrap_or(job.max_n).min(job.max_n);
    Ok(Ctx { mode: q.order.unwrap_or(job.order), job, rs, n })
}

impl Ctx {
    fn view(&self) -> Arc<OrderView> {
        self.rs.order(self.mode, self.n)
    }
}

/// The variant-by-position raster (4 bytes per cell; see `sm_agg::Raster::to_bytes`).
pub async fn raster(State(app): State<AppState>, Path(id): Path<u64>, Query(q): Query<ViewQ>) -> ApiResult<Response> {
    let c = ctx(&app, id, &q)?;
    let (s, e) = region(&app, q.start, q.end)?;
    let (w, h) = (q.width.unwrap_or(1000).clamp(1, MAX_WIDTH), q.height.unwrap_or(500).clamp(1, MAX_WIDTH));
    if w as u64 * h as u64 > MAX_CELLS {
        return Err(bad("raster too large"));
    }
    blocking(move || {
        let t = Instant::now();
        let ov = c.view();
        let rows = ov.rows() as u32;
        let (r0, r1) = (q.row0.unwrap_or(0).min(rows), q.row1.unwrap_or(rows).min(rows));
        let view = View { n: c.n, start: s, end: e, row0: r0, row1: r1.max(r0), width: w, height: h };
        let r = sm_agg::raster(&c.rs, &ov, &view);
        let mut bytes = r.to_bytes();
        if c.job.flipped {
            // Strands are relative to the variant; mirror them into the query's orientation.
            for cell in bytes.as_chunks_mut::<4>().0 {
                cell[3] = (cell[3] & 1) << 1 | (cell[3] & 2) >> 1;
            }
        }
        binary(
            bytes,
            &[
                ("x-width", w.to_string()),
                ("x-height", h.to_string()),
                ("x-rows", rows.to_string()),
                ("x-row0", r0.to_string()),
                ("x-row1", r1.max(r0).to_string()),
                ("x-sites-drawn", r.sites_drawn.to_string()),
                ("x-render-ms", format!("{:.2}", t.elapsed().as_secs_f64() * 1e3)),
            ],
        )
    })
    .await
}

/// Sites per column by mismatch count: `width × (n + 1)` u32, column-major.
pub async fn histogram(
    State(app): State<AppState>,
    Path(id): Path<u64>,
    Query(q): Query<ViewQ>,
) -> ApiResult<Response> {
    let c = ctx(&app, id, &q)?;
    let (s, e) = region(&app, q.start, q.end)?;
    let w = q.width.unwrap_or(1000).clamp(1, MAX_WIDTH);
    blocking(move || {
        let h = sm_agg::histogram(&c.rs, c.n, s, e, w);
        let bytes: Vec<u8> = h.iter().flat_map(|x| x.to_le_bytes()).collect();
        binary(bytes, &[("x-width", w.to_string()), ("x-stride", (c.n as u32 + 1).to_string())])
    })
    .await
}

fn label(job: &Job, chars: &[(u8, u8)]) -> String {
    chars
        .iter()
        .map(|&(p, b)| {
            let (p, b) = if job.flipped { (K as u8 - 1 - p, 3 - b) } else { (p, b) };
            format!("{}{}", p + 1, kmer::BASE_CHAR[b as usize] as char)
        })
        .collect::<Vec<_>>()
        .join(" ")
}

/// Clades as nested row intervals, restricted to a row window and a minimum size.
pub async fn tree(State(app): State<AppState>, Path(id): Path<u64>, Query(q): Query<ViewQ>) -> ApiResult<Json<Value>> {
    let c = ctx(&app, id, &q)?;
    blocking(move || {
        let ov = c.view();
        let rows = ov.rows() as u32;
        let (r0, r1) = (q.row0.unwrap_or(0), q.row1.unwrap_or(rows).min(rows));
        let min = q.min_rows.unwrap_or(2).max(2);
        let limit = q.limit.unwrap_or(20_000).min(200_000);
        let clades: Vec<Value> = ov
            .clades
            .iter()
            .filter(|cl| cl.end > r0 && cl.start < r1 && cl.end - cl.start >= min && cl.depth > 0)
            .take(limit)
            .map(|cl| json!({ "depth": cl.depth, "start": cl.start, "end": cl.end, "label": label(&c.job, &ov.clade_label(cl)) }))
            .collect();
        Json(json!({
            "n": c.n,
            "order": ov.mode,
            "rows": rows,
            "tree_length": ov.tree_length,
            "star_length": ov.star_length,
            "clades": clades,
        }))
    })
    .await
}

/// Ordered variants: sequence (query orientation), mismatches and their positions, site count.
pub async fn variants(
    State(app): State<AppState>,
    Path(id): Path<u64>,
    Query(q): Query<ViewQ>,
) -> ApiResult<Json<Value>> {
    let c = ctx(&app, id, &q)?;
    blocking(move || {
        let ov = c.view();
        let off = q.offset.unwrap_or(0).min(ov.rows());
        let lim = q.limit.unwrap_or(100).min(10_000);
        let query = c.job.kmer;
        let rows: Vec<Value> = ov.order[off..(off + lim).min(ov.rows())]
            .iter()
            .enumerate()
            .map(|(i, &vi)| {
                let v = c.job.orient(c.rs.variants.kmer[vi as usize]);
                json!({
                    "row": off + i,
                    "kmer": kmer::to_string(v),
                    "mismatches": c.rs.variants.mismatches[vi as usize],
                    "positions": kmer::mismatch_positions(query, v).map(|p| p + 1).collect::<Vec<_>>(),
                    "sites": c.rs.site_count[vi as usize],
                })
            })
            .collect();
        Json(
            json!({ "n": c.n, "order": ov.mode, "rows": ov.rows(), "query": kmer::to_string(query), "variants": rows }),
        )
    })
    .await
}

/// Concrete sites in a region (and row window), for zoomed-in views.
pub async fn sites(State(app): State<AppState>, Path(id): Path<u64>, Query(q): Query<ViewQ>) -> ApiResult<Json<Value>> {
    let c = ctx(&app, id, &q)?;
    let (s, e) = region(&app, q.start, q.end)?;
    blocking(move || {
        let ov = c.view();
        let m = c.rs.count(c.n) as u32;
        let (r0, r1) = (q.row0.unwrap_or(0), q.row1.unwrap_or(u32::MAX));
        let lim = q.limit.unwrap_or(5000).min(100_000);
        let mut out = Vec::new();
        let mut total = 0usize;
        for i in c.rs.sites.range(s, e) {
            let var = c.rs.sites.var[i];
            let vi = var & !sm_agg::result::REVERSE;
            if vi >= m {
                continue;
            }
            let row = ov.row_of[vi as usize];
            if row < r0 || row >= r1 {
                continue;
            }
            total += 1;
            if out.len() >= lim {
                continue;
            }
            let pos = c.rs.sites.pos[i];
            let (ci, off) = app.idx.contigs.to_local(pos).unwrap();
            let rev = (var & sm_agg::result::REVERSE != 0) != c.job.flipped;
            out.push(json!({
                "pos": pos,
                "contig": app.idx.contigs.get(ci).name,
                "display": app.display(ci),
                "start": off + 1,
                "strand": if rev { "-" } else { "+" },
                "row": row,
                "mismatches": c.rs.variants.mismatches[vi as usize],
                "kmer": kmer::to_string(c.job.orient(c.rs.variants.kmer[vi as usize])),
            }));
        }
        Json(json!({ "n": c.n, "total": total, "sites": out }))
    })
    .await
}

/// Streamed export of every site within n: TSV or BED.
pub async fn export(State(app): State<AppState>, Path(id): Path<u64>, Query(q): Query<ViewQ>) -> ApiResult<Response> {
    let c = ctx(&app, id, &q)?;
    let bed = match q.format.as_deref().unwrap_or("tsv") {
        "tsv" => false,
        "bed" => true,
        f => return Err(bad(format!("unknown format {f:?}"))),
    };
    let (tx, rx) = tokio::sync::mpsc::channel::<Result<bytes::Bytes, std::io::Error>>(8);
    let a = app.clone();
    tokio::task::spawn_blocking(move || {
        let ov = c.view();
        let m = c.rs.count(c.n) as u32;
        let mut buf = String::with_capacity(1 << 20);
        if !bed {
            buf.push_str("contig\tstart\tend\tstrand\tkmer\tmismatches\trow\n");
        }
        for i in 0..c.rs.sites.pos.len() {
            let var = c.rs.sites.var[i];
            let vi = var & !sm_agg::result::REVERSE;
            if vi >= m {
                continue;
            }
            let pos = c.rs.sites.pos[i];
            let (ci, off) = a.idx.contigs.to_local(pos).unwrap();
            let name = &a.idx.contigs.get(ci).name;
            let rev = (var & sm_agg::result::REVERSE != 0) != c.job.flipped;
            let strand = if rev { '-' } else { '+' };
            let km = kmer::to_string(c.job.orient(c.rs.variants.kmer[vi as usize]));
            let mm = c.rs.variants.mismatches[vi as usize];
            let row = ov.row_of[vi as usize];
            use std::fmt::Write;
            if bed {
                let _ = writeln!(buf, "{name}\t{off}\t{}\t{km}\t{mm}\t{strand}", off + K);
            } else {
                let _ = writeln!(buf, "{name}\t{}\t{}\t{strand}\t{km}\t{mm}\t{row}", off + 1, off + K);
            }
            if buf.len() > 1 << 20 && tx.blocking_send(Ok(bytes::Bytes::from(std::mem::take(&mut buf)))).is_err() {
                return; // client went away
            }
        }
        let _ = tx.blocking_send(Ok(bytes::Bytes::from(buf)));
    });
    let body = Body::from_stream(tokio_stream::wrappers::ReceiverStream::new(rx));
    let fname = format!(
        "sentromap-{}-n{}.{}",
        kmer::to_string(c_kmer(&app, id)),
        q.n.unwrap_or(0),
        if bed { "bed" } else { "tsv" }
    );
    Ok((
        [
            (header::CONTENT_TYPE, "text/plain".to_string()),
            (header::CONTENT_DISPOSITION, format!("attachment; filename=\"{fname}\"")),
        ],
        body,
    )
        .into_response())
}

fn c_kmer(app: &App, id: u64) -> u64 {
    app.job(id).map_or(0, |j| j.kmer)
}

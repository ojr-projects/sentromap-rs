//! Engine selection by a calibrated cost model (design §7.5).

use std::time::Instant;

use sm_core::synth::Rng;
use sm_index::Index;

use crate::pigeonhole::{pigeonhole, work_estimate};
use crate::scan::scan;
use crate::variants::Variants;

#[derive(Clone, Copy, Debug, PartialEq, Eq, serde::Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Engine {
    Scan,
    Pigeonhole,
}

/// Per-unit costs in nanoseconds of wall time on the whole worker pool.
#[derive(Clone, Copy, Debug, serde::Serialize)]
pub struct CostModel {
    /// Per k-mer in a full scan.
    pub scan_ns: f64,
    /// Fixed overhead per search (setup, merge).
    pub fixed_ns: f64,
    /// Extra fixed overhead when a pigeonhole walk fans out to the pool.
    pub parallel_ns: f64,
    /// Per pigeonhole walk step (bucket reached), including scanning its entries. Entries per
    /// bucket are fixed for a given index (N / 4^P), so they are not costed separately.
    pub walk_ns: f64,
    /// Per variant emitted.
    pub hit_ns: f64,
}

impl Default for CostModel {
    /// Rough values for a 12-core desktop; replace with [`CostModel::calibrate`].
    fn default() -> Self {
        Self { scan_ns: 0.25, fixed_ns: 5_000.0, parallel_ns: 500_000.0, walk_ns: 2.0, hit_ns: 20.0 }
    }
}

/// A plan for one query: the chosen engine and its estimated time.
#[derive(Clone, Copy, Debug, serde::Serialize)]
pub struct Plan {
    pub engine: Engine,
    pub estimate_ms: f64,
    pub scan_ms: f64,
    pub pigeonhole_ms: Option<f64>,
    pub expected_chance_hits: f64,
}

impl CostModel {
    pub fn scan_ms(&self, idx: &Index) -> f64 {
        (self.fixed_ns + idx.manifest.distinct_kmers as f64 * self.scan_ns) / 1e6
    }

    pub fn pigeonhole_ms(&self, idx: &Index, n: u32) -> Option<f64> {
        idx.b.as_ref()?;
        let (steps, entries) = work_estimate(idx.manifest.prefix_len, n, idx.manifest.distinct_kmers);
        let parallel = if steps + entries < crate::pigeonhole::SEQUENTIAL_WORK { 0.0 } else { self.parallel_ns };
        Some((self.fixed_ns + parallel + steps * self.walk_ns) / 1e6)
    }

    pub fn plan(&self, idx: &Index, n: u32) -> Plan {
        let scan_ms = self.scan_ms(idx);
        let ph = self.pigeonhole_ms(idx, n);
        let engine = match ph {
            Some(t) if t < scan_ms => Engine::Pigeonhole,
            _ => Engine::Scan,
        };
        let expected = sm_core::ball::expected_chance_hits(idx.manifest.distinct_kmers, n);
        let search = if engine == Engine::Pigeonhole { ph.unwrap() } else { scan_ms };
        Plan {
            engine,
            estimate_ms: search + expected * self.hit_ns / 1e6,
            scan_ms,
            pigeonhole_ms: ph,
            expected_chance_hits: expected,
        }
    }

    /// Measure the constants on a loaded index with a short microbenchmark (a second or two).
    pub fn calibrate(idx: &Index) -> Self {
        let mut m = Self::default();
        let mut rng = Rng::new(0xC0FFEE);
        let queries: Vec<u64> = (0..3).map(|_| rng.kmer()).collect();
        // Warm the mapping, then time the scan at n = 0 (no output).
        scan(idx, queries[0], 0);
        let t = median((0..3).map(|i| time_ns(|| scan(idx, queries[i], 0))).collect());
        m.scan_ns = (t / idx.manifest.distinct_kmers.max(1) as f64).max(0.01);

        if idx.b.is_some() {
            // Fixed cost: an exact lookup (n = 0) does almost no walking.
            m.fixed_ns = median((0..3).map(|i| time_ns(|| pigeonhole(idx, queries[i], 0))).collect());
            // Fit `t = parallel + walk · steps` over parallel-regime walks (random queries, so
            // output is small), minimising relative error so every n counts equally.
            let p = idx.manifest.prefix_len;
            let mut rows = Vec::new();
            for n in (4..=17).step_by(2) {
                let (steps, entries) = work_estimate(p, n, idx.manifest.distinct_kmers);
                if steps + entries < crate::pigeonhole::SEQUENTIAL_WORK {
                    continue;
                }
                if steps * m.walk_ns > 20.0 * m.scan_ns * idx.manifest.distinct_kmers as f64 {
                    break; // far past the crossover; don't waste calibration time
                }
                let t = median((0..3).map(|i| time_ns(|| pigeonhole(idx, queries[i], n))).collect());
                rows.push((steps, (t - m.fixed_ns).max(1.0)));
            }
            if let Some((a, b)) = fit_relative(&rows) {
                m.parallel_ns = a;
                m.walk_ns = b;
            }
        }
        m
    }
}

fn time_ns<T>(f: impl FnOnce() -> T) -> f64 {
    let t = Instant::now();
    std::hint::black_box(f());
    t.elapsed().as_nanos() as f64
}

fn median(mut v: Vec<f64>) -> f64 {
    v.sort_by(f64::total_cmp);
    v[v.len() / 2]
}

/// Fit `t ≈ a + b·x` with `a, b ≥ 0`, minimising `Σ ((a + b·x − t) / t)²`.
fn fit_relative(rows: &[(f64, f64)]) -> Option<(f64, f64)> {
    if rows.len() < 2 {
        return rows.first().map(|&(x, t)| (0.0, t / x));
    }
    let (mut s00, mut s01, mut s11, mut s0t, mut s1t) = (0.0, 0.0, 0.0, 0.0, 0.0);
    for &(x, t) in rows {
        let w = 1.0 / (t * t);
        s00 += w;
        s01 += w * x;
        s11 += w * x * x;
        s0t += w * t;
        s1t += w * x * t;
    }
    let det = s00 * s11 - s01 * s01;
    let (a, b) = if det > 0.0 { ((s0t * s11 - s1t * s01) / det, (s1t * s00 - s0t * s01) / det) } else { (-1.0, -1.0) };
    if a >= 0.0 && b >= 0.0 {
        return Some((a, b));
    }
    // Clamp: slope only, through the origin.
    Some((0.0, (s1t / s11).max(0.0)))
}

/// Search with the engine the cost model prefers (or a forced one).
pub fn search(idx: &Index, query: u64, max_n: u32, model: &CostModel, force: Option<Engine>) -> (Variants, Plan) {
    let mut plan = model.plan(idx, max_n);
    if let Some(e) = force {
        plan.engine = e;
    }
    let v = match plan.engine {
        Engine::Scan => scan(idx, query, max_n),
        Engine::Pigeonhole => pigeonhole(idx, query, max_n),
    };
    (v, plan)
}

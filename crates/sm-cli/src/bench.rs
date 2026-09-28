//! `sentromap bench search`: engine comparison over n (port of the prototype's searchprof,
//! design §16.2).

use std::io::Write;
use std::time::Instant;

use anyhow::{Result, ensure};
use rayon::prelude::*;
use sm_core::kmer;
use sm_core::synth::Rng;
use sm_index::Index;
use sm_search::{CostModel, Engine};

pub struct QuerySet {
    pub class: &'static str,
    pub queries: Vec<u64>,
}

/// Unique (single-copy) windows, the highest-copy k-mer, and random (usually absent) k-mers.
pub fn query_sets(idx: &Index, per_class: usize, seed: u64) -> Vec<QuerySet> {
    let mut rng = Rng::new(seed);
    let n = idx.manifest.distinct_kmers as u32;
    let mut unique = Vec::new();
    while unique.len() < per_class {
        let leaf = (rng.next_u64() % n as u64) as u32;
        if idx.row_len(leaf) == 1 {
            unique.push(idx.a.kmer_at(leaf as usize));
        }
    }
    let top = (0..n).into_par_iter().max_by_key(|&l| (idx.row_len(l), std::cmp::Reverse(l))).unwrap();
    let random = (0..per_class).map(|_| rng.kmer()).collect();
    vec![
        QuerySet { class: "unique", queries: unique },
        QuerySet { class: "top-repeat", queries: vec![idx.a.kmer_at(top as usize)] },
        QuerySet { class: "random", queries: random },
    ]
}

fn time_ms<T>(f: impl FnOnce() -> T) -> (T, f64) {
    let t = Instant::now();
    let r = f();
    (r, t.elapsed().as_secs_f64() * 1e3)
}

pub fn search_bench(idx: &Index, max_n: u32, per_class: usize, out: &mut impl Write) -> Result<()> {
    let t = Instant::now();
    let model = CostModel::calibrate(idx);
    writeln!(out, "# calibrated in {:.2} s: {model:?}", t.elapsed().as_secs_f64())?;
    let sets = query_sets(idx, per_class, 42);
    for s in &sets {
        if s.class == "top-repeat" {
            let leaf = idx.a.lookup(kmer::canonical(s.queries[0])).unwrap();
            writeln!(out, "# top-repeat {} occurs {} times", kmer::to_string(s.queries[0]), idx.row_len(leaf))?;
        }
    }
    writeln!(
        out,
        "class\tn\tvariants\tsites\tscan_ms\tpigeonhole_ms\tchosen\tpredicted_ms\tscan_max_ms\tpigeonhole_max_ms"
    )?;
    for s in &sets {
        for n in 0..=max_n {
            let (mut ts, mut tp) = (Vec::new(), Vec::new());
            let (mut variants, mut sites) = (0usize, 0u64);
            for &q in &s.queries {
                let (mut a, t1) = time_ms(|| sm_search::scan(idx, q, n));
                ts.push(t1);
                if idx.b.is_some() {
                    let (mut b, t2) = time_ms(|| sm_search::pigeonhole(idx, q, n));
                    tp.push(t2);
                    a.sort_by_leaf();
                    b.sort_by_leaf();
                    ensure!(a == b, "engines disagree: {} n={n}", kmer::to_string(q));
                }
                variants += a.len();
                sites += a.site_count(idx);
            }
            let med = |v: &mut Vec<f64>| {
                v.sort_by(f64::total_cmp);
                if v.is_empty() { f64::NAN } else { v[v.len() / 2] }
            };
            let max = |v: &Vec<f64>| v.iter().copied().fold(f64::NAN, f64::max);
            let plan = model.plan(idx, n);
            let k = s.queries.len();
            writeln!(
                out,
                "{}\t{n}\t{:.0}\t{:.0}\t{:.2}\t{:.3}\t{}\t{:.3}\t{:.2}\t{:.3}",
                s.class,
                variants as f64 / k as f64,
                sites as f64 / k as f64,
                med(&mut ts),
                med(&mut tp),
                if plan.engine == Engine::Scan { "scan" } else { "pigeonhole" },
                plan.estimate_ms,
                max(&ts),
                max(&tp),
            )?;
            out.flush()?;
        }
    }
    Ok(())
}

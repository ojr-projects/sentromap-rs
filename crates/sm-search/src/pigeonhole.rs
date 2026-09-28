//! Engine 2: the pigeonhole two-copy walk (design §7.4).
//!
//! If a k-mer is within n of the query, its lead 15 bases or its back 16 bases are within
//! `r = ⌊n/2⌋`. Walk A enumerates copy A's prefixes within r of the query's lead bases; walk B
//! does the same over copy B (rotated, back bases leading). A hit is accepted by walk A when
//! its lead-15 distance is ≤ r, and by walk B only when it is > r, so every hit is found
//! exactly once. Both walks run for the query and for its reverse complement.

use rayon::prelude::*;
use sm_core::kmer::LEAD15;
use sm_core::{hamming, revcomp, rotate, unrotate};
use sm_index::Index;
use sm_index::kmers::{AnyKmers, SortedKmers, Word};

use crate::variants::Variants;

/// Low 30 bits of a rotated k-mer: the original lead 15 bases.
const ROT_LEAD15: u64 = (1 << 30) - 1;

/// Prefix bases enumerated up front to create parallel tasks (4^3 = 64 subtrees per walk).
const FAN_OUT: u32 = 3;

/// Below this estimated work (walk steps + entries) the walk runs on the calling thread.
/// Each bucket visit costs a couple of cache misses (~60 ns), so this is ~0.5 ms, about
/// what waking the pool costs.
pub const SEQUENTIAL_WORK: f64 = 10_000.0;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Copy {
    A,
    B,
}

/// Enumerate every bucket whose P-base prefix is within `r` of `qp` (the query's prefix),
/// starting from a fixed first `d0` bases (`pre0`, costing `c0`). Buckets come out ascending.
#[inline]
fn walk(qp: u64, p: u32, d0: u32, pre0: u64, c0: u32, r: u32, mut visit: impl FnMut(usize)) {
    let mut stack: Vec<(u32, u64, u32)> = Vec::with_capacity(4 * p as usize);
    stack.push((d0, pre0, c0));
    while let Some((d, pre, c)) = stack.pop() {
        if d == p {
            visit(pre as usize);
            continue;
        }
        let rest = p - d;
        if c == r {
            // Budget spent: the remaining bases must match the query.
            visit(((pre << (2 * rest)) | (qp & ((1 << (2 * rest)) - 1))) as usize);
            continue;
        }
        let qb = (qp >> (2 * (rest - 1))) & 3;
        for b in (0..4u64).rev() {
            stack.push((d + 1, (pre << 2) | b, c + (b != qb) as u32));
        }
    }
}

/// Walk tasks: `(first-bases prefix, cost)` for every subtree within budget.
fn subtrees(qp: u64, p: u32, r: u32) -> Vec<(u32, u64, u32)> {
    let f = FAN_OUT.min(p);
    let qf = qp >> (2 * (p - f));
    (0..1u64 << (2 * f)).map(|pre| (f, pre, hamming(pre, qf))).filter(|&(_, _, c)| c <= r).collect()
}

struct Task {
    copy: Copy,
    flipped: bool,
    /// The (possibly reverse-complemented) query this walk searches for.
    q: u64,
    d0: u32,
    pre0: u64,
    c0: u32,
}

/// Run the pigeonhole search. Panics if the index has no copy B.
pub fn pigeonhole(idx: &Index, query: u64, max_n: u32) -> Variants {
    let b = idx.b.as_ref().expect("pigeonhole search needs copy B");
    match (&idx.a, b) {
        (AnyKmers::U64(a), AnyKmers::U64(b)) => run(a, b, query, max_n),
        (AnyKmers::U32(a), AnyKmers::U32(b)) => run(a, b, query, max_n),
        _ => unreachable!("copies A and B share a word width"),
    }
}

fn run<W: Word>(a: &SortedKmers<W>, b: &SortedKmers<W>, query: u64, n: u32) -> Variants {
    let r = n / 2;
    let p = a.p();
    let shift = a.shift();
    let mut tasks = Vec::new();
    for (flipped, q) in [(false, query), (true, revcomp(query))] {
        for copy in [Copy::A, Copy::B] {
            let key = if copy == Copy::A { q } else { rotate(q) };
            for (d0, pre0, c0) in subtrees(key >> shift, p, r) {
                tasks.push(Task { copy, flipped, q, d0, pre0, c0 });
            }
        }
    }
    let run_task = |t: &Task| {
        let mut v = Variants::new(query, n as u8);
        let emit = |v: &mut Variants, k: u64, d: u32, leaf: u32| {
            let variant = if t.flipped { revcomp(k) } else { k };
            v.push(variant, d as u8, leaf, t.flipped);
        };
        match t.copy {
            Copy::A => {
                let q = t.q;
                let q_lead = q & LEAD15;
                walk(q >> shift, p, t.d0, t.pre0, t.c0, r, |bucket| {
                    let range = a.bucket_range(bucket);
                    let prefix = (bucket as u64) << shift;
                    let start = range.start;
                    for (i, w) in a.words()[range].iter().enumerate() {
                        let k = w.to_kmer(prefix);
                        let d = hamming(q, k);
                        if d <= n && hamming(q_lead, k & LEAD15) <= r {
                            emit(&mut v, k, d, (start + i) as u32);
                        }
                    }
                });
            }
            Copy::B => {
                let qr = rotate(t.q);
                let q_lead = qr & ROT_LEAD15;
                walk(qr >> shift, p, t.d0, t.pre0, t.c0, r, |bucket| {
                    let range = b.bucket_range(bucket);
                    let prefix = (bucket as u64) << shift;
                    for w in &b.words()[range] {
                        let kr = w.to_kmer(prefix);
                        let d = hamming(qr, kr);
                        if d <= n && hamming(q_lead, kr & ROT_LEAD15) > r {
                            let k = unrotate(kr);
                            let leaf = a.lookup(k).expect("copy B k-mer missing from copy A");
                            emit(&mut v, k, d, leaf);
                        }
                    }
                });
            }
        }
        v
    };
    // Small walks finish in microseconds; spawning tasks would dominate.
    let (steps, entries) = work_estimate(p, n, a.len() as u64);
    let parts: Vec<Variants> = if steps + entries < SEQUENTIAL_WORK {
        tasks.iter().map(run_task).collect()
    } else {
        tasks.par_iter().map(run_task).collect()
    };
    let mut out = Variants::new(query, n as u8);
    for mut part in parts {
        out.append(&mut part);
    }
    out
}

/// Work estimate for the cost model: `(walk steps, entries scanned)`, both strands.
pub fn work_estimate(p: u32, n: u32, distinct: u64) -> (f64, f64) {
    let r = n / 2;
    let ball = sm_core::ball::ball(p, r) as f64;
    let buckets = 4f64.powi(p as i32);
    let steps = 4.0 * ball.min(buckets);
    let entries = 4.0 * distinct as f64 * ball.min(buckets) / buckets;
    (steps, entries)
}

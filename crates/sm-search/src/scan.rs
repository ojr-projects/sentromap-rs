//! Engine 1: full scan of copy A, both strands in one pass (design §7.3).

use rayon::prelude::*;
use sm_core::{hamming, revcomp};
use sm_index::Index;
use sm_index::kmers::{AnyKmers, SortedKmers};

use crate::variants::{Variants, prefix_mask};

/// Leaves per parallel task.
const CHUNK: usize = 1 << 18;

/// Scan every k-mer in copy A for `query` and its reverse complement.
pub fn scan(idx: &Index, query: u64, max_n: u32) -> Variants {
    let rc = revcomp(query);
    let mut parts: Vec<Variants> = match &idx.a {
        AnyKmers::U64(c) => scan_flat(c, query, rc, max_n),
        AnyKmers::U32(c) => scan_buckets(c, query, rc, max_n),
    };
    let mut out = Variants::new(query, max_n as u8);
    for p in &mut parts {
        out.append(p);
    }
    out
}

/// Whole k-mers stored: stream the word array, the leaf id is the index.
fn scan_flat(c: &SortedKmers<u64>, q: u64, rc: u64, n: u32) -> Vec<Variants> {
    c.words()
        .par_chunks(CHUNK)
        .enumerate()
        .map(|(ci, words)| {
            let base = ci * CHUNK;
            let mut v = Variants::new(q, n as u8);
            for (i, &k) in words.iter().enumerate() {
                let dq = hamming(q, k);
                let dr = hamming(rc, k);
                if dq <= n {
                    v.push(k, dq as u8, (base + i) as u32, false);
                }
                if dr <= n {
                    v.push(revcomp(k), dr as u8, (base + i) as u32, true);
                }
            }
            v
        })
        .collect()
}

/// Suffix words stored: walk buckets, rebuilding k-mers from the implicit prefix, with
/// bucket-level early rejection.
fn scan_buckets(c: &SortedKmers<u32>, q: u64, rc: u64, n: u32) -> Vec<Variants> {
    let shift = c.shift();
    let pm = prefix_mask(c.p());
    let table = c.prefix_table();
    // Split the bucket space into ranges holding roughly CHUNK leaves each.
    let mut bounds = vec![0usize];
    let mut next = CHUNK as u32;
    while (next as usize) < c.len() {
        let b = table.partition_point(|&s| s < next);
        if b > *bounds.last().unwrap() {
            bounds.push(b);
        }
        next = next.saturating_add(CHUNK as u32);
    }
    bounds.push(c.buckets());
    bounds
        .par_windows(2)
        .map(|w| {
            let mut v = Variants::new(q, n as u8);
            for b in w[0]..w[1] {
                let r = c.bucket_range(b);
                if r.is_empty() {
                    continue;
                }
                let prefix = (b as u64) << shift;
                if hamming(q & pm, prefix) > n && hamming(rc & pm, prefix) > n {
                    continue;
                }
                let start = r.start;
                for (i, &s) in c.words()[r].iter().enumerate() {
                    let k = prefix | s as u64;
                    let dq = hamming(q, k);
                    let dr = hamming(rc, k);
                    if dq <= n {
                        v.push(k, dq as u8, (start + i) as u32, false);
                    }
                    if dr <= n {
                        v.push(revcomp(k), dr as u8, (start + i) as u32, true);
                    }
                }
            }
            v
        })
        .collect()
}

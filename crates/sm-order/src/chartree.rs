//! The shared-derived-character tree (design §9.2).
//!
//! Every variant differs from the query by up to n substitutions; each is a character
//! `(position, base)`. Characters are ranked by how many variants share them. A variant's key
//! is its characters' ranks, ascending (most shared first), packed 7 bits each into a `u128`.
//! Sorting the keys is a pre-order traversal of the trie over rank lists: ancestors come before
//! descendants and every clade is contiguous.

use rayon::prelude::*;
use sm_core::K;
use sm_core::kmer::{base_at, diff_mask};

/// Character ids: `position * 4 + base`.
pub const CHARACTERS: usize = 4 * K as usize;
/// Ranks packed per key; limits the usable n.
pub const MAX_RANKS: usize = 18;
const FIELD: u32 = 7;
/// Below this many variants, work runs on the calling thread (pool wake-up costs ~0.5 ms).
const PARALLEL_MIN: usize = 50_000;
const TOP_SHIFT: u32 = 128 - FIELD;

/// The characters of `variant` relative to `query`, as ids.
#[inline]
pub fn characters(query: u64, variant: u64) -> impl Iterator<Item = u8> {
    let mut m = diff_mask(query, variant);
    std::iter::from_fn(move || {
        if m == 0 {
            return None;
        }
        let bit = 63 - m.leading_zeros();
        m &= !(1u64 << bit);
        let pos = K - 1 - bit / 2;
        Some((pos * 4 + base_at(variant, pos) as u32) as u8)
    })
}

/// A ranking of characters: rank 1 is the most shared; 0 means absent.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Ranking {
    pub rank_of: [u8; CHARACTERS],
    /// Character id for each rank (index 0 unused).
    pub char_of: Vec<u8>,
    pub counts: [u64; CHARACTERS],
}

impl Ranking {
    /// Rank characters by frequency across the variants with `mismatches ≤ n`, ties by id.
    pub fn new(query: u64, kmers: &[u64], mismatches: &[u8], n: u8) -> Self {
        let count_chunk = |ks: &[u64], ms: &[u8]| {
            let mut c = [0u64; CHARACTERS];
            for (&k, &m) in ks.iter().zip(ms) {
                if m <= n {
                    for id in characters(query, k) {
                        c[id as usize] += 1;
                    }
                }
            }
            c
        };
        let counts = if kmers.len() < PARALLEL_MIN {
            count_chunk(kmers, mismatches)
        } else {
            kmers.par_chunks(1 << 16).zip(mismatches.par_chunks(1 << 16)).map(|(ks, ms)| count_chunk(ks, ms)).reduce(
                || [0u64; CHARACTERS],
                |mut a, b| {
                    for i in 0..CHARACTERS {
                        a[i] += b[i];
                    }
                    a
                },
            )
        };
        let mut ids: Vec<u8> = (0..CHARACTERS as u8).filter(|&i| counts[i as usize] > 0).collect();
        ids.sort_by_key(|&i| (std::cmp::Reverse(counts[i as usize]), i));
        let mut rank_of = [0u8; CHARACTERS];
        let mut char_of = vec![0u8];
        for (r, &id) in ids.iter().enumerate() {
            rank_of[id as usize] = (r + 1) as u8;
            char_of.push(id);
        }
        Self { rank_of, char_of, counts }
    }

    /// Packed key of a variant: its ranks ascending, 7 bits each from the top, zero-terminated.
    #[inline]
    pub fn key(&self, query: u64, variant: u64) -> u128 {
        let mut ranks = [0u8; 32];
        let mut len = 0;
        for id in characters(query, variant) {
            ranks[len] = self.rank_of[id as usize];
            len += 1;
        }
        debug_assert!(len <= MAX_RANKS, "variant beyond the key-packing limit");
        debug_assert!(ranks[..len].iter().all(|&r| r > 0), "character missing from the ranking");
        let ranks = &mut ranks[..len];
        ranks.sort_unstable();
        let mut key = 0u128;
        for (i, &r) in ranks.iter().enumerate() {
            key |= (r as u128) << (TOP_SHIFT - FIELD * i as u32);
        }
        key
    }

    /// The characters `(position, base)` on the path of the first `depth` ranks of `key`.
    pub fn path(&self, key: u128, depth: usize) -> Vec<(u8, u8)> {
        (0..depth)
            .map(|i| {
                let r = (key >> (TOP_SHIFT - FIELD * i as u32)) as u8 & 0x7f;
                let id = self.char_of[r as usize];
                (id / 4, id % 4)
            })
            .collect()
    }
}

/// Shared rank fields between two keys.
#[inline]
pub fn lcp(a: u128, b: u128) -> u8 {
    if a == b { MAX_RANKS as u8 } else { ((a ^ b).leading_zeros() / FIELD) as u8 }
}

/// A nested interval of the ordered variants: every member shares the first `depth` ranks.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Clade {
    pub depth: u8,
    /// Half-open range of positions in the order.
    pub start: u32,
    pub end: u32,
}

/// Variants in cladogram order, with the tree's structure.
#[derive(Clone, Debug, Default)]
pub struct Ordered {
    /// Indices into the input variants, in order.
    pub order: Vec<u32>,
    /// Keys in order (sorted ascending).
    pub keys: Vec<u128>,
    /// `lcp[i]`: shared ranks between ordered items `i − 1` and `i` (`lcp[0] = 0`: the root).
    pub lcp: Vec<u8>,
    /// Parsimony score: substitutions along the tree's edges.
    pub tree_length: u64,
}

impl Ordered {
    pub fn len(&self) -> usize {
        self.order.len()
    }

    pub fn is_empty(&self) -> bool {
        self.order.is_empty()
    }

    /// Clades as nested intervals (only those with two or more members, plus the root).
    /// Sorted by start, then outer before inner.
    pub fn clades(&self) -> Vec<Clade> {
        let m = self.order.len() as u32;
        let mut out = Vec::new();
        let mut stack: Vec<(u8, u32)> = vec![(0, 0)];
        for i in 1..=m {
            let l = if i < m { self.lcp[i as usize] } else { 0 };
            let mut start = i - 1;
            while stack.last().unwrap().0 > l {
                let (d, s) = stack.pop().unwrap();
                out.push(Clade { depth: d, start: s, end: i });
                start = s;
            }
            if stack.last().unwrap().0 < l {
                stack.push((l, start));
            }
        }
        out.push(Clade { depth: 0, start: 0, end: m });
        out.sort_by_key(|c| (c.start, std::cmp::Reverse(c.end), c.depth));
        out
    }
}

/// Order the variants with `mismatches ≤ n` by the character tree under `ranking`.
pub fn order(query: u64, kmers: &[u64], mismatches: &[u8], n: u8, ranking: &Ranking) -> Ordered {
    assert!(n as usize <= MAX_RANKS, "n = {n} exceeds the key-packing limit of {MAX_RANKS}");
    let keyed = if kmers.len() < PARALLEL_MIN {
        let mut keyed: Vec<(u128, u32)> = (0..kmers.len())
            .filter(|&i| mismatches[i] <= n)
            .map(|i| (ranking.key(query, kmers[i]), i as u32))
            .collect();
        keyed.sort_unstable_by_key(|&(k, _)| k);
        keyed
    } else {
        let mut keyed: Vec<(u128, u32)> = kmers
            .par_iter()
            .zip(mismatches.par_iter())
            .enumerate()
            .filter(|(_, (_, m))| **m <= n)
            .map(|(i, (&k, _))| (ranking.key(query, k), i as u32))
            .collect();
        keyed.par_sort_unstable_by_key(|&(k, _)| k);
        keyed
    };
    finish(keyed, mismatches)
}

fn finish(keyed: Vec<(u128, u32)>, mismatches: &[u8]) -> Ordered {
    let m = keyed.len();
    let mut lcp = vec![0u8; m];
    let mut tree_length = 0u64;
    for i in 0..m {
        if i > 0 {
            lcp[i] = self::lcp(keyed[i - 1].0, keyed[i].0);
        }
        tree_length += (mismatches[keyed[i].1 as usize] - lcp[i]) as u64;
    }
    let (keys, order) = keyed.into_iter().unzip();
    Ordered { order, keys, lcp, tree_length }
}

/// Per-n ordering (best trees): rank on the variants within n, then order them.
pub fn order_per_n(query: u64, kmers: &[u64], mismatches: &[u8], n: u8) -> (Ranking, Ordered) {
    let ranking = Ranking::new(query, kmers, mismatches, n);
    let o = order(query, kmers, mismatches, n, &ranking);
    (ranking, o)
}

/// Frozen ordering at `n` from the max-n order: a subsequence of it (design §9.6).
pub fn filter_frozen(max_order: &Ordered, mismatches: &[u8], n: u8) -> Ordered {
    let keyed: Vec<(u128, u32)> = max_order
        .keys
        .iter()
        .zip(&max_order.order)
        .filter(|&(_, &i)| mismatches[i as usize] <= n)
        .map(|(&k, &i)| (k, i))
        .collect();
    finish(keyed, mismatches)
}

/// The star tree's length: every variant hangs off the root.
pub fn star_length(mismatches: &[u8], n: u8) -> u64 {
    mismatches.iter().filter(|&&m| m <= n).map(|&m| m as u64).sum()
}

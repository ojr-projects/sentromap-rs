//! A query's result, held by the server (design §8): variants sorted by mismatch count (so
//! "within n" is a prefix), per-variant site counts, and every site sorted by genome position.

use std::collections::HashMap;
use std::sync::{Arc, Mutex};

use rayon::prelude::*;
use sm_index::Index;
use sm_search::Variants;

use crate::order::{OrderMode, OrderView};

/// Top bit of [`SiteTable::var`]: the site is on the reverse strand relative to the variant.
pub const REVERSE: u32 = 1 << 31;

/// Every site of a result, sorted by position.
#[derive(Default)]
pub struct SiteTable {
    pub pos: Vec<u32>,
    /// Variant index, with [`REVERSE`] set for '-' sites.
    pub var: Vec<u32>,
}

impl SiteTable {
    /// Index range of sites whose 31-mer overlaps global `[start, end)`.
    pub fn range(&self, start: u32, end: u32) -> std::ops::Range<usize> {
        let a = self.pos.partition_point(|&p| p + sm_core::K <= start);
        let b = self.pos.partition_point(|&p| p < end);
        a..b.max(a)
    }
}

pub struct ResultSet {
    /// The canonical form of the query; results are in its orientation.
    pub query: u64,
    pub max_n: u8,
    /// Sorted by `(mismatches, leaf, flipped)`.
    pub variants: Variants,
    /// `n_end[n]`: number of variants with at most n mismatches.
    pub n_end: Vec<u32>,
    pub site_count: Vec<u32>,
    pub sites: SiteTable,
    /// Sites with at most n mismatches, per n.
    pub sites_within: Vec<u64>,
    orders: Mutex<HashMap<(OrderMode, u8), Arc<OrderView>>>,
}

impl ResultSet {
    /// Assemble from a search for the canonical query.
    pub fn new(idx: &Index, mut v: Variants) -> Self {
        let max_n = v.max_n;
        let mut perm: Vec<usize> = (0..v.len()).collect();
        perm.par_sort_unstable_by_key(|&i| (v.mismatches[i], v.leaf[i], v.flipped[i]));
        v.permute(&perm);
        let mut n_end = vec![0u32; max_n as usize + 1];
        for &m in &v.mismatches {
            n_end[m as usize] += 1;
        }
        for n in 1..n_end.len() {
            n_end[n] += n_end[n - 1];
        }
        let site_count: Vec<u32> = v.leaf.par_iter().map(|&l| idx.row_len(l) as u32).collect();
        let total: usize = site_count.iter().map(|&c| c as usize).sum();
        // Gather sites in parallel chunks of variants, then sort by position.
        let chunks: Vec<Vec<(u32, u32)>> = (0..v.len())
            .collect::<Vec<_>>()
            .par_chunks(4096)
            .map(|ix| {
                let mut out = Vec::new();
                for &i in ix {
                    for (pos, strand) in idx.sites(v.leaf[i], v.flipped[i]) {
                        let rev = if strand == sm_core::Strand::Reverse { REVERSE } else { 0 };
                        out.push((pos, i as u32 | rev));
                    }
                }
                out
            })
            .collect();
        let mut pairs: Vec<(u32, u32)> = Vec::with_capacity(total);
        for c in chunks {
            pairs.extend(c);
        }
        pairs.par_sort_unstable();
        let (pos, var) = pairs.into_iter().unzip();
        let mut sites_within = vec![0u64; max_n as usize + 1];
        for (i, &c) in site_count.iter().enumerate() {
            sites_within[v.mismatches[i] as usize] += c as u64;
        }
        for n in 1..sites_within.len() {
            sites_within[n] += sites_within[n - 1];
        }
        Self {
            query: v.query,
            max_n,
            variants: v,
            n_end,
            site_count,
            sites: SiteTable { pos, var },
            sites_within,
            orders: Mutex::default(),
        }
    }

    /// Number of variants within `n`.
    pub fn count(&self, n: u8) -> usize {
        self.n_end[n.min(self.max_n) as usize] as usize
    }

    /// The ordering for `(mode, n)`, computed on first use and cached.
    pub fn order(&self, mode: OrderMode, n: u8) -> Arc<OrderView> {
        let n = n.min(self.max_n);
        let mode = mode.resolve(self.count(n));
        if let Some(o) = self.orders.lock().unwrap().get(&(mode, n)) {
            return o.clone();
        }
        let base = if mode == OrderMode::Frozen && n < self.max_n {
            Some(self.order(OrderMode::Frozen, self.max_n))
        } else {
            None
        };
        let o = Arc::new(OrderView::build(self, mode, n, base.as_deref()));
        self.orders.lock().unwrap().insert((mode, n), o.clone());
        o
    }

    /// Approximate heap size, for the cache budget.
    pub fn bytes(&self) -> usize {
        let v = self.variants.len() * (8 + 1 + 4 + 1 + 4);
        let s = self.sites.pos.len() * 8;
        let o: usize = self.orders.lock().unwrap().values().map(|o| o.bytes()).sum();
        v + s + o
    }
}

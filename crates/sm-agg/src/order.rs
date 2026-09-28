//! Orderings of a result at a given n (design §9): per-n character tree, frozen (subsequence of
//! the max-n order), or the exact MST for small sets.

use sm_order::{Clade, Ranking};

use crate::result::ResultSet;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum OrderMode {
    /// Re-rank characters at each n: best trees, order changes as n changes.
    PerN,
    /// Rank once at max n: every smaller-n order is a subsequence of the max-n order.
    Frozen,
    /// Exact MST when the set is small enough, else per-n.
    Mst,
}

impl OrderMode {
    pub fn resolve(self, count: usize) -> Self {
        if self == OrderMode::Mst && count > sm_order::MST_LIMIT { OrderMode::PerN } else { self }
    }
}

pub struct OrderView {
    pub mode: OrderMode,
    pub n: u8,
    /// Variant indices in row order.
    pub order: Vec<u32>,
    /// Row of each variant within n (indexed by variant).
    pub row_of: Vec<u32>,
    /// Nested clades over rows, sorted by start then outer first.
    pub clades: Vec<Clade>,
    /// For character trees: the ranking and each row's key, to label clades.
    pub ranking: Option<Ranking>,
    pub keys: Vec<u128>,
    pub tree_length: u64,
    pub star_length: u64,
}

impl OrderView {
    pub fn build(rs: &ResultSet, mode: OrderMode, n: u8, frozen_base: Option<&OrderView>) -> Self {
        let m = rs.count(n);
        let kmers = &rs.variants.kmer[..m];
        let mism = &rs.variants.mismatches[..m];
        let star_length = sm_order::star_length(mism, n);
        let (order, clades, ranking, keys, tree_length) = match mode {
            OrderMode::PerN | OrderMode::Frozen => {
                let (ranking, o) = match (mode, frozen_base) {
                    (OrderMode::Frozen, Some(base)) => {
                        let full = sm_order::Ordered {
                            order: base.order.clone(),
                            keys: base.keys.clone(),
                            lcp: Vec::new(),
                            tree_length: 0,
                        };
                        (base.ranking.clone().unwrap(), sm_order::filter_frozen(&full, &rs.variants.mismatches, n))
                    }
                    _ => sm_order::order_per_n(rs.query, kmers, mism, n),
                };
                let clades = o.clades();
                (o.order, clades, Some(ranking), o.keys, o.tree_length)
            }
            OrderMode::Mst => {
                let t = sm_order::mst(rs.query, kmers, mism, n);
                let clades = mst_clades(&t);
                (t.order, clades, None, Vec::new(), t.tree_length)
            }
        };
        let mut row_of = vec![u32::MAX; m];
        for (r, &v) in order.iter().enumerate() {
            row_of[v as usize] = r as u32;
        }
        Self { mode, n, order, row_of, clades, ranking, keys, tree_length, star_length }
    }

    pub fn rows(&self) -> usize {
        self.order.len()
    }

    /// The characters `(position, base)` shared by a clade (character-tree modes only).
    pub fn clade_label(&self, c: &Clade) -> Vec<(u8, u8)> {
        match &self.ranking {
            Some(r) if (c.start as usize) < self.keys.len() => r.path(self.keys[c.start as usize], c.depth as usize),
            _ => Vec::new(),
        }
    }

    pub fn bytes(&self) -> usize {
        self.order.len() * 8 + self.keys.len() * 16 + self.clades.len() * 12
    }
}

/// Subtrees of an MST in pre-order are contiguous: each node with children is a clade.
fn mst_clades(t: &sm_order::Mst) -> Vec<Clade> {
    let n = t.order.len();
    let mut pos = std::collections::HashMap::with_capacity(n);
    for (p, &v) in t.order.iter().enumerate() {
        pos.insert(v, p);
    }
    let mut size = vec![1u32; n];
    let mut depth = vec![0u8; n];
    for p in 0..n {
        let par = t.parent[t.order[p] as usize];
        if par != sm_order::mst::ROOT {
            depth[p] = depth[pos[&par]].saturating_add(1);
        }
    }
    for p in (0..n).rev() {
        let par = t.parent[t.order[p] as usize];
        if par != sm_order::mst::ROOT {
            size[pos[&par]] += size[p];
        }
    }
    let mut out: Vec<Clade> = (0..n)
        .filter(|&p| size[p] > 1)
        .map(|p| Clade { depth: depth[p] + 1, start: p as u32, end: p as u32 + size[p] })
        .collect();
    out.push(Clade { depth: 0, start: 0, end: n as u32 });
    out.sort_by_key(|c| (c.start, std::cmp::Reverse(c.end), c.depth));
    out
}

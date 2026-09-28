//! Columnar search output (design §7.1, §8): one row per variant, sites resolved lazily.

use sm_core::{Hit, K};
use sm_index::Index;

/// Every genome k-mer within `max_n` of the query, in the query's orientation.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Variants {
    pub query: u64,
    pub max_n: u8,
    pub kmer: Vec<u64>,
    pub mismatches: Vec<u8>,
    pub leaf: Vec<u32>,
    /// Came from the reverse-complement search: the variant is `revcomp(leaf k-mer)`.
    pub flipped: Vec<bool>,
}

impl Variants {
    pub fn new(query: u64, max_n: u8) -> Self {
        Self { query, max_n, ..Default::default() }
    }

    pub fn len(&self) -> usize {
        self.kmer.len()
    }

    pub fn is_empty(&self) -> bool {
        self.kmer.is_empty()
    }

    #[inline]
    pub fn push(&mut self, kmer: u64, mismatches: u8, leaf: u32, flipped: bool) {
        self.kmer.push(kmer);
        self.mismatches.push(mismatches);
        self.leaf.push(leaf);
        self.flipped.push(flipped);
    }

    pub fn append(&mut self, other: &mut Variants) {
        self.kmer.append(&mut other.kmer);
        self.mismatches.append(&mut other.mismatches);
        self.leaf.append(&mut other.leaf);
        self.flipped.append(&mut other.flipped);
    }

    /// Sort rows by `(leaf, flipped)` so engines can be compared directly.
    pub fn sort_by_leaf(&mut self) {
        let mut idx: Vec<usize> = (0..self.len()).collect();
        idx.sort_unstable_by_key(|&i| (self.leaf[i], self.flipped[i]));
        self.permute(&idx);
    }

    pub fn permute(&mut self, idx: &[usize]) {
        self.kmer = idx.iter().map(|&i| self.kmer[i]).collect();
        self.mismatches = idx.iter().map(|&i| self.mismatches[i]).collect();
        self.leaf = idx.iter().map(|&i| self.leaf[i]).collect();
        self.flipped = idx.iter().map(|&i| self.flipped[i]).collect();
    }

    /// Count of variants per mismatch count `0..=max_n`.
    pub fn histogram(&self) -> Vec<u64> {
        let mut h = vec![0u64; self.max_n as usize + 1];
        for &m in &self.mismatches {
            h[m as usize] += 1;
        }
        h
    }

    /// Expand to concrete sites, sorted like the oracle's output.
    pub fn sites(&self, idx: &Index) -> Vec<Hit> {
        let mut out = Vec::new();
        for i in 0..self.len() {
            for (pos, strand) in idx.sites(self.leaf[i], self.flipped[i]) {
                out.push(Hit { pos, strand, variant: self.kmer[i], mismatches: self.mismatches[i] });
            }
        }
        out.sort_unstable();
        out
    }

    pub fn site_count(&self, idx: &Index) -> u64 {
        self.leaf.iter().map(|&l| idx.row_len(l)).sum()
    }
}

/// Mask of the first `p` bases of a k-mer.
#[inline]
pub fn prefix_mask(p: u32) -> u64 {
    sm_core::kmer::KMER_MASK & !((1u64 << (2 * (K - p))) - 1)
}

//! Brute-force reference search (design §15.1). Every engine must agree with it.

use crate::kmer::{hamming, revcomp, windows};

/// Strand of a site relative to its variant (design §5.2).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Strand {
    /// The forward strand at the site reads the variant.
    Forward,
    /// The forward strand holds the variant's reverse complement.
    Reverse,
}

impl Strand {
    pub fn symbol(self) -> char {
        match self {
            Strand::Forward => '+',
            Strand::Reverse => '-',
        }
    }
}

/// One genome site within `n` of the query, with its variant in query orientation.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Hit {
    /// Global coordinate of the site's first base on the forward strand.
    pub pos: u32,
    pub strand: Strand,
    pub variant: u64,
    pub mismatches: u8,
}

/// Scan `seq` (a whole genome in global coordinates) window by window, on both strands.
/// `contig_bounds` are `(offset, len)` pairs; k-mers never span two contigs.
/// Returns hits sorted by `(pos, strand, variant)`.
pub fn search(seq: &[u8], contig_bounds: impl IntoIterator<Item = (u32, u32)>, query: u64, max_n: u32) -> Vec<Hit> {
    let mut hits = Vec::new();
    for (offset, len) in contig_bounds {
        let s = &seq[offset as usize..(offset + len) as usize];
        for (i, fwd) in windows(s) {
            let pos = offset + i as u32;
            let d = hamming(query, fwd);
            if d <= max_n {
                hits.push(Hit { pos, strand: Strand::Forward, variant: fwd, mismatches: d as u8 });
            }
            let rc = revcomp(fwd);
            let d = hamming(query, rc);
            if d <= max_n {
                hits.push(Hit { pos, strand: Strand::Reverse, variant: rc, mismatches: d as u8 });
            }
        }
    }
    hits.sort_unstable();
    hits
}

/// [`search`] over an in-memory [`Genome`](crate::fasta::Genome).
pub fn search_genome(g: &crate::fasta::Genome, query: u64, max_n: u32) -> Vec<Hit> {
    search(&g.seq, g.contigs.iter().map(|c| (c.offset, c.len)), query, max_n)
}

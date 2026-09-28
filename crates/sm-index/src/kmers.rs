//! A sorted set of k-mers with a prefix table: copy A, and later copy B (design §6.2–6.3).

use std::ops::Range;
use std::path::Path;

use anyhow::{Result, ensure};
use bytemuck::Pod;
use sm_core::K;

use crate::array::{Array, MapOptions};

/// Storage word for sorted k-mers. With P = 15 the 16-base suffix fits a `u32`;
/// for smaller P the whole k-mer is kept as a `u64`, so no reconstruction is needed.
pub trait Word: Pod + Ord + Send + Sync + 'static {
    const BITS: u32;
    fn from_kmer(k: u64) -> Self;
    /// Rebuild the k-mer from its bucket prefix (already shifted into place) and this word.
    fn to_kmer(self, prefix_bits: u64) -> u64;
}

impl Word for u64 {
    const BITS: u32 = 64;
    #[inline]
    fn from_kmer(k: u64) -> Self {
        k
    }
    #[inline]
    fn to_kmer(self, _prefix_bits: u64) -> u64 {
        self
    }
}

impl Word for u32 {
    const BITS: u32 = 32;
    #[inline]
    fn from_kmer(k: u64) -> Self {
        k as u32
    }
    #[inline]
    fn to_kmer(self, prefix_bits: u64) -> u64 {
        prefix_bits | self as u64
    }
}

pub struct SortedKmers<W: Word> {
    p: u32,
    prefix: Array<u32>,
    words: Array<W>,
}

impl<W: Word> SortedKmers<W> {
    pub fn open(dir: &Path, name: &str, p: u32, opts: MapOptions) -> Result<Self> {
        let prefix: Array<u32> = Array::open(&dir.join(format!("{name}.prefix")), "prefix", opts)?;
        let words: Array<W> = Array::open(&dir.join(format!("{name}.words")), "words", opts)?;
        ensure!(prefix.len() == (1usize << (2 * p)) + 1, "{name}: prefix table size does not match P = {p}");
        ensure!(*prefix.last().unwrap() as usize == words.len(), "{name}: prefix table does not end at N");
        if W::BITS == 32 {
            ensure!(p == 15, "u32 words need P = 15");
        }
        Ok(Self { p, prefix, words })
    }

    /// Prefix length in bases.
    pub fn p(&self) -> u32 {
        self.p
    }

    /// Bit shift that extracts the prefix: `kmer >> shift()`.
    pub fn shift(&self) -> u32 {
        2 * (K - self.p)
    }

    pub fn len(&self) -> usize {
        self.words.len()
    }

    pub fn is_empty(&self) -> bool {
        self.words.is_empty()
    }

    pub fn buckets(&self) -> usize {
        self.prefix.len() - 1
    }

    pub fn prefix_table(&self) -> &[u32] {
        &self.prefix
    }

    pub fn words(&self) -> &[W] {
        &self.words
    }

    #[inline]
    pub fn bucket_range(&self, b: usize) -> Range<usize> {
        self.prefix[b] as usize..self.prefix[b + 1] as usize
    }

    /// Index of `kmer` (the leaf id, for copy A), if present.
    #[inline]
    pub fn lookup(&self, kmer: u64) -> Option<u32> {
        let b = (kmer >> self.shift()) as usize;
        let r = self.bucket_range(b);
        let start = r.start;
        self.words[r].binary_search(&W::from_kmer(kmer)).ok().map(|i| (start + i) as u32)
    }

    /// The k-mer at index `i`.
    pub fn kmer_at(&self, i: usize) -> u64 {
        let b = self.prefix.partition_point(|&s| s as usize <= i) - 1;
        self.words[i].to_kmer((b as u64) << self.shift())
    }

    pub fn advise_sequential(&self) {
        self.words.advise(memmap2::Advice::Sequential);
    }
}

/// Copy A or B with its word width resolved at open time.
pub enum AnyKmers {
    U32(SortedKmers<u32>),
    U64(SortedKmers<u64>),
}

impl AnyKmers {
    pub fn open(dir: &Path, name: &str, p: u32, word_bits: u32, opts: MapOptions) -> Result<Self> {
        Ok(match word_bits {
            32 => AnyKmers::U32(SortedKmers::open(dir, name, p, opts)?),
            64 => AnyKmers::U64(SortedKmers::open(dir, name, p, opts)?),
            w => anyhow::bail!("unsupported word width {w}"),
        })
    }

    pub fn len(&self) -> usize {
        match self {
            AnyKmers::U32(c) => c.len(),
            AnyKmers::U64(c) => c.len(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn lookup(&self, kmer: u64) -> Option<u32> {
        match self {
            AnyKmers::U32(c) => c.lookup(kmer),
            AnyKmers::U64(c) => c.lookup(kmer),
        }
    }

    pub fn kmer_at(&self, i: usize) -> u64 {
        match self {
            AnyKmers::U32(c) => c.kmer_at(i),
            AnyKmers::U64(c) => c.kmer_at(i),
        }
    }
}

/// Dispatch a generic body over the concrete word type.
#[macro_export]
macro_rules! with_kmers {
    ($any:expr, $c:ident => $body:expr) => {
        match $any {
            $crate::kmers::AnyKmers::U32($c) => $body,
            $crate::kmers::AnyKmers::U64($c) => $body,
        }
    };
}

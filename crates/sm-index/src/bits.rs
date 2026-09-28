//! Bit vectors with a sampled select index (design §6.4: row starts, strand bits).

use std::path::Path;

use anyhow::Result;

use crate::array::{Array, ArrayWriter, MapOptions};

/// Streams bits into a `u64`-word array file; the bit length goes in the header's aux field.
pub struct BitWriter {
    words: ArrayWriter<u64>,
    cur: u64,
    len: u64,
}

impl BitWriter {
    pub fn create(path: &Path, kind: &str) -> Result<Self> {
        Ok(Self { words: ArrayWriter::create(path, kind)?, cur: 0, len: 0 })
    }

    #[inline]
    pub fn push(&mut self, bit: bool) -> Result<()> {
        self.cur |= (bit as u64) << (self.len % 64);
        self.len += 1;
        if self.len.is_multiple_of(64) {
            self.words.push(self.cur)?;
            self.cur = 0;
        }
        Ok(())
    }

    pub fn len(&self) -> u64 {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    pub fn finish(mut self) -> Result<()> {
        if !self.len.is_multiple_of(64) {
            self.words.push(self.cur)?;
        }
        self.words.set_aux(self.len);
        self.words.finish()
    }
}

/// A memory-mapped bit vector.
pub struct BitVec {
    words: Array<u64>,
    len: u64,
}

impl BitVec {
    pub fn open(path: &Path, kind: &str, opts: MapOptions) -> Result<Self> {
        let words = Array::open(path, kind, opts)?;
        let len = words.aux();
        anyhow::ensure!(len.div_ceil(64) == words.len() as u64, "{}: bit length mismatch", path.display());
        Ok(Self { words, len })
    }

    pub fn len(&self) -> u64 {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    #[inline]
    pub fn get(&self, i: u64) -> bool {
        debug_assert!(i < self.len);
        (self.words[(i / 64) as usize] >> (i % 64)) & 1 == 1
    }

    pub fn words(&self) -> &[u64] {
        &self.words
    }
}

/// One sample per this many set bits.
pub const SELECT_SAMPLE: u64 = 64;

/// Builds select samples while bits are streamed: `samples[j]` is the position of the
/// `(j * SELECT_SAMPLE)`-th set bit.
#[derive(Default)]
pub struct SelectBuilder {
    pub samples: Vec<u64>,
    ones: u64,
}

impl SelectBuilder {
    #[inline]
    pub fn observe(&mut self, pos: u64, bit: bool) {
        if bit {
            if self.ones.is_multiple_of(SELECT_SAMPLE) {
                self.samples.push(pos);
            }
            self.ones += 1;
        }
    }

    pub fn ones(&self) -> u64 {
        self.ones
    }
}

/// `select1` over a [`BitVec`] using stored samples.
pub struct Select {
    samples: Array<u64>,
    ones: u64,
}

impl Select {
    pub fn open(path: &Path, kind: &str, opts: MapOptions) -> Result<Self> {
        let samples = Array::open(path, kind, opts)?;
        let ones = samples.aux();
        Ok(Self { samples, ones })
    }

    pub fn ones(&self) -> u64 {
        self.ones
    }

    /// Position of the `j`-th set bit (0-based).
    #[inline]
    pub fn select1(&self, bits: &BitVec, j: u64) -> u64 {
        debug_assert!(j < self.ones);
        let start = self.samples[(j / SELECT_SAMPLE) as usize];
        let mut remaining = j % SELECT_SAMPLE;
        let words = bits.words();
        let mut w = (start / 64) as usize;
        let mut word = words[w] & (!0u64 << (start % 64));
        loop {
            let c = word.count_ones() as u64;
            if remaining < c {
                for _ in 0..remaining {
                    word &= word - 1;
                }
                return w as u64 * 64 + word.trailing_zeros() as u64;
            }
            remaining -= c;
            w += 1;
            word = words[w];
        }
    }
}

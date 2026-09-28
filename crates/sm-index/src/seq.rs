//! The genome sequence, 2-bit packed, with runs of non-ACGT bases and soft-masked (lowercase)
//! runs, so the server can serve sequence and resolve "the 31-mer at position x".

use std::path::Path;

use anyhow::Result;
use sm_core::K;
use sm_core::kmer::{BASE_CHAR, BASE_CODE, INVALID, KMER_MASK};

use crate::array::{Array, ArrayWriter, MapOptions};

/// Streams bases: 32 per `u64`, first base in the high bits. Non-ACGT bases are stored as A
/// and recorded in `other` runs; lowercase runs in `soft`.
pub struct SeqWriter {
    words: ArrayWriter<u64>,
    cur: u64,
    len: u64,
    other: RunWriter,
    soft: RunWriter,
}

struct RunWriter {
    out: ArrayWriter<u32>,
    open: Option<u32>,
}

impl RunWriter {
    fn create(path: &Path, kind: &str) -> Result<Self> {
        Ok(Self { out: ArrayWriter::create(path, kind)?, open: None })
    }
    #[inline]
    fn step(&mut self, pos: u32, inside: bool) -> Result<()> {
        match (self.open, inside) {
            (None, true) => self.open = Some(pos),
            (Some(s), false) => {
                self.out.extend(&[s, pos])?;
                self.open = None;
            }
            _ => {}
        }
        Ok(())
    }
    fn finish(mut self, end: u32) -> Result<()> {
        self.step(end, false)?;
        self.out.finish()
    }
}

impl SeqWriter {
    pub fn create(dir: &Path) -> Result<Self> {
        Ok(Self {
            words: ArrayWriter::create(&dir.join("seq.packed"), "seq")?,
            cur: 0,
            len: 0,
            other: RunWriter::create(&dir.join("seq.other"), "runs")?,
            soft: RunWriter::create(&dir.join("seq.soft"), "runs")?,
        })
    }

    pub fn extend(&mut self, seq: &[u8]) -> Result<()> {
        for &c in seq {
            let code = BASE_CODE[c as usize];
            let pos = self.len as u32;
            self.other.step(pos, code == INVALID)?;
            self.soft.step(pos, c.is_ascii_lowercase())?;
            let code = if code == INVALID { 0 } else { code as u64 };
            self.cur |= code << (62 - 2 * (self.len % 32));
            self.len += 1;
            if self.len.is_multiple_of(32) {
                self.words.push(self.cur)?;
                self.cur = 0;
            }
        }
        Ok(())
    }

    pub fn finish(mut self) -> Result<()> {
        if !self.len.is_multiple_of(32) {
            self.words.push(self.cur)?;
        }
        self.words.set_aux(self.len);
        self.words.finish()?;
        let end = self.len as u32;
        self.other.finish(end)?;
        self.soft.finish(end)
    }
}

pub struct Sequence {
    words: Array<u64>,
    len: u64,
    /// Flattened `[start, end)` pairs, ascending.
    other: Array<u32>,
    soft: Array<u32>,
}

impl Sequence {
    pub fn open(dir: &Path, opts: MapOptions) -> Result<Self> {
        let words: Array<u64> = Array::open(&dir.join("seq.packed"), "seq", opts)?;
        let len = words.aux();
        Ok(Self {
            words,
            len,
            other: Array::open(&dir.join("seq.other"), "runs", opts)?,
            soft: Array::open(&dir.join("seq.soft"), "runs", opts)?,
        })
    }

    pub fn len(&self) -> u64 {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    #[inline]
    fn code(&self, pos: u64) -> u8 {
        ((self.words[(pos / 32) as usize] >> (62 - 2 * (pos % 32))) & 3) as u8
    }

    fn in_runs(runs: &[u32], pos: u32) -> bool {
        Self::runs_in(runs, pos, pos + 1).next().is_some()
    }

    /// Runs of non-ACGT bases overlapping `[start, end)`.
    fn runs_in(runs: &[u32], start: u32, end: u32) -> impl Iterator<Item = (u32, u32)> + '_ {
        let pairs: &[[u32; 2]] = bytemuck::cast_slice(runs);
        let first = pairs.partition_point(|r| r[1] <= start);
        pairs[first..].iter().take_while(move |r| r[0] < end).map(|r| (r[0], r[1]))
    }

    /// Bases `[start, end)` as ASCII, with non-ACGT shown as `N` and soft-masking preserved.
    pub fn fetch(&self, start: u32, end: u32) -> Vec<u8> {
        let end = end.min(self.len as u32);
        let mut out: Vec<u8> = (start as u64..end as u64).map(|p| BASE_CHAR[self.code(p) as usize]).collect();
        for (s, e) in Self::runs_in(&self.soft, start, end) {
            for c in &mut out[(s.max(start) - start) as usize..(e.min(end) - start) as usize] {
                c.make_ascii_lowercase();
            }
        }
        for (s, e) in Self::runs_in(&self.other, start, end) {
            for c in &mut out[(s.max(start) - start) as usize..(e.min(end) - start) as usize] {
                *c = b'N';
            }
        }
        out
    }

    /// The forward k-mer starting at global `pos`, if all 31 bases are ACGT.
    /// (Callers must also check it does not cross a contig boundary.)
    pub fn kmer_at(&self, pos: u32) -> Option<u64> {
        if pos as u64 + K as u64 > self.len {
            return None;
        }
        if Self::runs_in(&self.other, pos, pos + K).next().is_some() {
            return None;
        }
        let mut k = 0u64;
        for p in pos as u64..pos as u64 + K as u64 {
            k = (k << 2) | self.code(p) as u64;
        }
        Some(k & KMER_MASK)
    }

    pub fn is_other(&self, pos: u32) -> bool {
        Self::in_runs(&self.other, pos)
    }
}

//! The positions "tether" (design §6.4): row `l` lists every genome position of leaf `l`.

use std::ops::Range;
use std::path::Path;

use anyhow::{Result, ensure};

use crate::array::{Array, MapOptions};
use crate::bits::{BitVec, Select};

pub struct Positions {
    sites: Array<u32>,
    /// 1 where the forward strand at the site holds the canonical k-mer.
    fwd_canonical: BitVec,
    /// 1 at the first site of each row, plus a sentinel 1 at `sites.len()`.
    row_start: BitVec,
    select: Select,
}

impl Positions {
    pub fn open(dir: &Path, opts: MapOptions) -> Result<Self> {
        let sites = Array::open(&dir.join("pos.sites"), "sites", opts)?;
        let fwd_canonical = BitVec::open(&dir.join("pos.strand"), "strand", opts)?;
        let row_start = BitVec::open(&dir.join("pos.rowstart"), "rowstart", opts)?;
        let select = Select::open(&dir.join("pos.select"), "select", opts)?;
        ensure!(fwd_canonical.len() == sites.len() as u64, "strand bits do not match site count");
        ensure!(row_start.len() == sites.len() as u64 + 1, "row-start bits do not match site count");
        Ok(Self { sites, fwd_canonical, row_start, select })
    }

    pub fn site_count(&self) -> u64 {
        self.sites.len() as u64
    }

    /// Number of rows (distinct k-mers).
    pub fn rows(&self) -> u64 {
        self.select.ones() - 1
    }

    /// Site index range of row `leaf`.
    #[inline]
    pub fn row(&self, leaf: u32) -> Range<u64> {
        let s = self.select.select1(&self.row_start, leaf as u64);
        let e = self.select.select1(&self.row_start, leaf as u64 + 1);
        s..e
    }

    /// `(global position, forward strand holds the canonical k-mer)` of site `i`.
    #[inline]
    pub fn site(&self, i: u64) -> (u32, bool) {
        (self.sites[i as usize], self.fwd_canonical.get(i))
    }

    pub fn sites(&self) -> &[u32] {
        &self.sites
    }

    pub fn row_start_bits(&self) -> &BitVec {
        &self.row_start
    }

    pub fn fwd_canonical_bits(&self) -> &BitVec {
        &self.fwd_canonical
    }
}

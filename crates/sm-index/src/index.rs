//! An opened, memory-mapped index.

use std::path::{Path, PathBuf};

use anyhow::{Result, ensure};
use sm_core::{Contigs, Strand};

use crate::array::MapOptions;
use crate::kmers::AnyKmers;
use crate::manifest::{self, Manifest};
use crate::positions::Positions;
use crate::seq::Sequence;

pub struct Index {
    pub dir: PathBuf,
    pub manifest: Manifest,
    pub contigs: Contigs,
    /// Copy A: all distinct canonical k-mers, sorted; a k-mer's index is its leaf id.
    pub a: AnyKmers,
    pub positions: Positions,
    pub seq: Sequence,
}

impl Index {
    pub fn open(dir: &Path) -> Result<Self> {
        Self::open_with(dir, MapOptions::default())
    }

    pub fn open_with(dir: &Path, opts: MapOptions) -> Result<Self> {
        let manifest = manifest::read_manifest(dir)?;
        ensure!(manifest.k == sm_core::K, "index built for k = {}", manifest.k);
        let contigs = manifest::read_contigs(&dir.join(manifest::CONTIGS))?;
        let a = AnyKmers::open(dir, "a", manifest.prefix_len, manifest.word_bits, opts)?;
        let positions = Positions::open(dir, opts)?;
        let seq = Sequence::open(dir, opts)?;
        ensure!(seq.len() == manifest.genome_len, "sequence length does not match the manifest");
        ensure!(a.len() as u64 == manifest.distinct_kmers, "copy A size does not match the manifest");
        ensure!(positions.rows() == manifest.distinct_kmers, "positions rows do not match copy A");
        ensure!(positions.site_count() == manifest.sites, "site count does not match the manifest");
        Ok(Self { dir: dir.to_path_buf(), manifest, contigs, a, positions, seq })
    }

    /// Sites of a variant from leaf `leaf`, `flipped` if it came from the reverse-complement
    /// search, as `(global position, strand relative to the variant)` (design §7.6).
    pub fn sites(&self, leaf: u32, flipped: bool) -> impl Iterator<Item = (u32, Strand)> + '_ {
        self.positions.row(leaf).map(move |i| {
            let (pos, fwd_canonical) = self.positions.site(i);
            (pos, if fwd_canonical != flipped { Strand::Forward } else { Strand::Reverse })
        })
    }

    /// The forward k-mer at 0-based `offset` in contig `contig`, if valid and inside the contig.
    pub fn kmer_at(&self, contig: usize, offset: u32) -> Option<u64> {
        let c = self.contigs.get(contig);
        if offset as u64 + sm_core::K as u64 > c.len as u64 {
            return None;
        }
        self.seq.kmer_at(c.offset + offset)
    }

    pub fn row_len(&self, leaf: u32) -> u64 {
        let r = self.positions.row(leaf);
        r.end - r.start
    }
}

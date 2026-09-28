//! Quantitative tracks (ChIP signal, conservation, GC%) stored as binned pyramids in global
//! coordinates: level 0 has `BIN0`-base bins, each level above merges `FACTOR` bins. Every bin
//! holds the coverage-weighted mean and the maximum (NaN where there is no data).

use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use rayon::prelude::*;
use sm_core::Contigs;
use sm_index::array::{Array, MapOptions, write_array};
use sm_index::seq::Sequence;

use crate::alias::Aliases;

pub const BIN0: u32 = 32;
pub const FACTOR: u32 = 4;

#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct TrackInfo {
    pub name: String,
    pub kind: String,
    pub description: String,
    pub source: String,
    pub bin0: u32,
    pub factor: u32,
    pub levels: u32,
    /// Robust display range: 1st and 99th percentiles of level-0 means.
    pub p01: f32,
    pub p99: f32,
    pub min: f32,
    pub max: f32,
    /// Fraction of level-0 bins with data.
    pub coverage: f32,
    pub unknown_chroms: Vec<String>,
}

/// Level-0 accumulators while importing.
struct Acc {
    sum: Vec<f64>,
    cov: Vec<u32>,
    max: Vec<f32>,
}

impl Acc {
    fn new(genome_len: u64) -> Self {
        let n = genome_len.div_ceil(BIN0 as u64) as usize;
        Self { sum: vec![0.0; n], cov: vec![0; n], max: vec![f32::NEG_INFINITY; n] }
    }

    /// Add `value` over global `[start, end)`.
    #[inline]
    fn add(&mut self, start: u32, end: u32, value: f32) {
        if !value.is_finite() {
            return;
        }
        let mut s = start;
        while s < end {
            let b = (s / BIN0) as usize;
            let e = end.min((b as u32 + 1) * BIN0);
            self.sum[b] += value as f64 * (e - s) as f64;
            self.cov[b] += e - s;
            if value > self.max[b] {
                self.max[b] = value;
            }
            s = e;
        }
    }

    fn level0(&self) -> Vec<[f32; 2]> {
        self.sum
            .par_iter()
            .zip(&self.cov)
            .zip(&self.max)
            .map(|((&s, &c), &m)| if c == 0 { [f32::NAN, f32::NAN] } else { [(s / c as f64) as f32, m] })
            .collect()
    }
}

fn level_path(dir: &Path, name: &str, level: u32) -> PathBuf {
    dir.join(format!("track.{name}.L{level}"))
}

/// Write the pyramid for level-0 accumulators; returns (levels, stats).
fn write_pyramid(dir: &Path, name: &str, acc: &Acc) -> Result<(u32, [f32; 5])> {
    let l0 = acc.level0();
    // Stats from level-0 means.
    let mut vals: Vec<f32> = l0.iter().map(|v| v[0]).filter(|v| !v.is_nan()).collect();
    let coverage = vals.len() as f32 / l0.len().max(1) as f32;
    let stats = if vals.is_empty() {
        [0.0, 0.0, 0.0, 0.0, 0.0]
    } else {
        vals.par_sort_unstable_by(f32::total_cmp);
        let q = |f: f64| vals[((vals.len() - 1) as f64 * f) as usize];
        [q(0.01), q(0.99), vals[0], *vals.last().unwrap(), coverage]
    };
    // Carry sums and coverage up the pyramid so means stay coverage-weighted.
    let mut sum = acc.sum.clone();
    let mut cov: Vec<u64> = acc.cov.iter().map(|&c| c as u64).collect();
    let mut max = acc.max.clone();
    let mut level = 0;
    write_array(&level_path(dir, name, 0), "track", &l0, BIN0 as u64)?;
    while sum.len() > 1024 {
        let n = sum.len().div_ceil(FACTOR as usize);
        let (mut s2, mut c2, mut m2) = (vec![0.0; n], vec![0u64; n], vec![f32::NEG_INFINITY; n]);
        for i in 0..sum.len() {
            let j = i / FACTOR as usize;
            s2[j] += sum[i];
            c2[j] += cov[i];
            m2[j] = m2[j].max(max[i]);
        }
        (sum, cov, max) = (s2, c2, m2);
        level += 1;
        let bin = BIN0 as u64 * (FACTOR as u64).pow(level);
        let lv: Vec<[f32; 2]> = (0..sum.len())
            .map(|i| if cov[i] == 0 { [f32::NAN, f32::NAN] } else { [(sum[i] / cov[i] as f64) as f32, max[i]] })
            .collect();
        write_array(&level_path(dir, name, level), "track", &lv, bin)?;
    }
    Ok((level + 1, stats))
}

fn info_from(
    name: &str,
    kind: &str,
    description: &str,
    source: &str,
    levels: u32,
    st: [f32; 5],
    unknown: Vec<String>,
) -> TrackInfo {
    TrackInfo {
        name: name.to_string(),
        kind: kind.to_string(),
        description: description.to_string(),
        source: source.to_string(),
        bin0: BIN0,
        factor: FACTOR,
        levels,
        p01: st[0],
        p99: st[1],
        min: st[2],
        max: st[3],
        coverage: st[4],
        unknown_chroms: unknown,
    }
}

/// Descriptive fields of a track.
pub struct TrackMeta<'a> {
    pub name: &'a str,
    pub kind: &'a str,
    pub description: &'a str,
    pub source: &'a str,
}

/// Import a bigWig file.
pub fn import_bigwig(
    dir: &Path,
    path: &Path,
    meta: &TrackMeta,
    contigs: &Contigs,
    aliases: &Aliases,
) -> Result<TrackInfo> {
    let TrackMeta { name, kind, description, source } = *meta;
    let mut bw = bigtools::BigWigRead::open_file(path).map_err(|e| anyhow::anyhow!("{}: {e}", path.display()))?;
    let chroms: Vec<(String, u32)> = bw.chroms().iter().map(|c| (c.name.clone(), c.length)).collect();
    let mut acc = Acc::new(contigs.total_len());
    let mut unknown = Vec::new();
    for (chrom, len) in chroms {
        let Some(ci) = aliases.contig(&chrom) else {
            unknown.push(chrom);
            continue;
        };
        let c = contigs.get(ci);
        let iter = bw.get_interval(&chrom, 0, len).map_err(|e| anyhow::anyhow!("{chrom}: {e}"))?;
        for v in iter {
            let v = v.map_err(|e| anyhow::anyhow!("{chrom}: {e}"))?;
            let end = v.end.min(c.len);
            if v.start < end {
                acc.add(c.offset + v.start, c.offset + end, v.value);
            }
        }
    }
    let (levels, st) = write_pyramid(dir, name, &acc)?;
    Ok(info_from(name, kind, description, source, levels, st, unknown))
}

/// GC fraction of the ACGT bases, computed from the stored sequence.
pub fn import_gc(dir: &Path, seq: &Sequence, contigs: &Contigs) -> Result<TrackInfo> {
    let mut acc = Acc::new(contigs.total_len());
    let chunk = 1u32 << 20;
    let parts: Vec<(u32, Vec<u8>)> = (0..contigs.total_len() as u32)
        .step_by(chunk as usize)
        .collect::<Vec<_>>()
        .into_par_iter()
        .map(|s| (s, seq.fetch(s, s.saturating_add(chunk))))
        .collect();
    for (s, bases) in parts {
        for (b, win) in bases.chunks(BIN0 as usize).enumerate() {
            let (mut gc, mut acgt) = (0u32, 0u32);
            for &c in win {
                match c.to_ascii_uppercase() {
                    b'G' | b'C' => {
                        gc += 1;
                        acgt += 1
                    }
                    b'A' | b'T' => acgt += 1,
                    _ => {}
                }
            }
            if acgt > 0 {
                let start = s + (b as u32) * BIN0;
                // Weight by ACGT bases so N runs don't dilute the mean.
                let bin = (start / BIN0) as usize;
                acc.sum[bin] += gc as f64;
                acc.cov[bin] += acgt;
                acc.max[bin] = acc.max[bin].max(gc as f32 / acgt as f32);
            }
        }
    }
    let (levels, st) = write_pyramid(dir, "GC", &acc)?;
    Ok(info_from("GC", "sequence", "GC fraction of ACGT bases", "computed from the genome", levels, st, Vec::new()))
}

pub struct Track {
    pub info: TrackInfo,
    levels: Vec<Array<[f32; 2]>>,
}

impl Track {
    pub fn open(dir: &Path, info: TrackInfo) -> Result<Self> {
        let levels = (0..info.levels)
            .map(|l| Array::open(&level_path(dir, &info.name, l), "track", MapOptions::default()))
            .collect::<Result<_>>()
            .with_context(|| format!("track {}", info.name))?;
        Ok(Self { info, levels })
    }

    /// `(mean, max)` per pixel column over global `[start, end)` split into `width` columns,
    /// from the coarsest level whose bins are no wider than a column.
    pub fn summary(&self, start: u32, end: u32, width: usize) -> Vec<[f32; 2]> {
        let span = (end - start) as f64;
        let per_col = span / width as f64;
        let mut level = 0usize;
        while level + 1 < self.levels.len() && self.levels[level + 1].aux() as f64 <= per_col {
            level += 1;
        }
        let bins = &self.levels[level];
        let bin = bins.aux() as u32;
        (0..width)
            .map(|c| {
                let s = start as f64 + c as f64 * per_col;
                let e = s + per_col;
                let b0 = (s as u32 / bin) as usize;
                let b1 = ((e.ceil() as u32).div_ceil(bin) as usize).max(b0 + 1).min(bins.len());
                let (mut sum, mut n, mut mx) = (0.0f64, 0u32, f32::NAN);
                for v in &bins[b0.min(bins.len())..b1] {
                    if !v[0].is_nan() {
                        sum += v[0] as f64;
                        n += 1;
                        mx = if mx.is_nan() { v[1] } else { mx.max(v[1]) };
                    }
                }
                if n == 0 { [f32::NAN, f32::NAN] } else { [(sum / n as f64) as f32, mx] }
            })
            .collect()
    }

    /// Level-0 mean at a global position.
    pub fn value_at(&self, pos: u32) -> f32 {
        self.levels[0].get((pos / BIN0) as usize).map_or(f32::NAN, |v| v[0])
    }
}

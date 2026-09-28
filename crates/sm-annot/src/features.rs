//! Interval features (genes, repeats, …) in global coordinates, stored as flat arrays sorted by
//! start with a running maximum of ends for overlap queries (design §11).

use std::collections::HashMap;
use std::io::BufRead;
use std::path::Path;

use anyhow::{Context, Result};
use sm_core::Contigs;
use sm_index::array::{Array, MapOptions, write_array};

use crate::alias::Aliases;

/// One feature, as parsed.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct Feature {
    /// Global coordinates, half-open.
    pub start: u32,
    pub end: u32,
    pub kind: String,
    pub name: String,
    /// b'+', b'-' or b'.'.
    pub strand: u8,
}

/// Counts from an import, for warnings.
#[derive(Clone, Debug, Default, serde::Serialize, serde::Deserialize)]
pub struct ImportStats {
    pub features: u64,
    pub unknown_seqid: u64,
    pub unknown_seqids: Vec<String>,
    pub skipped: u64,
}

struct Sink<'a> {
    aliases: &'a Aliases,
    contigs: &'a Contigs,
    out: Vec<Feature>,
    stats: ImportStats,
}

impl Sink<'_> {
    fn push(&mut self, seqid: &str, start0: u64, end: u64, kind: &str, name: &str, strand: u8) {
        let Some(ci) = self.aliases.contig(seqid) else {
            self.stats.unknown_seqid += 1;
            if self.stats.unknown_seqids.len() < 20 && !self.stats.unknown_seqids.iter().any(|s| s == seqid) {
                self.stats.unknown_seqids.push(seqid.to_string());
            }
            return;
        };
        let c = self.contigs.get(ci);
        let end = end.min(c.len as u64);
        if start0 >= end {
            self.stats.skipped += 1;
            return;
        }
        self.out.push(Feature {
            start: c.offset + start0 as u32,
            end: c.offset + end as u32,
            kind: kind.to_string(),
            name: name.to_string(),
            strand,
        });
        self.stats.features += 1;
    }
}

fn open_lines(path: &Path) -> Result<Box<dyn BufRead + Send>> {
    sm_core::fasta::open(path).with_context(|| format!("opening {}", path.display()))
}

fn gff_attr<'a>(attrs: &'a str, key: &str) -> Option<&'a str> {
    attrs.split(';').find_map(|kv| kv.strip_prefix(key).and_then(|v| v.strip_prefix('=')))
}

/// GFF3: every feature except whole-sequence `region` records. Name from `Name`, else `ID`.
pub fn parse_gff(path: &Path, contigs: &Contigs, aliases: &Aliases) -> Result<(Vec<Feature>, ImportStats)> {
    let mut sink = Sink { aliases, contigs, out: Vec::new(), stats: ImportStats::default() };
    for line in open_lines(path)?.lines() {
        let line = line?;
        if line.starts_with('#') || line.is_empty() {
            if line.starts_with("##FASTA") {
                break;
            }
            continue;
        }
        let f: Vec<&str> = line.splitn(9, '\t').collect();
        if f.len() < 9 || f[2] == "region" {
            sink.stats.skipped += 1;
            continue;
        }
        let (Ok(s), Ok(e)) = (f[3].parse::<u64>(), f[4].parse::<u64>()) else {
            sink.stats.skipped += 1;
            continue;
        };
        let name = gff_attr(f[8], "Name").or_else(|| gff_attr(f[8], "ID")).unwrap_or("");
        let name = percent_decode(name);
        sink.push(f[0], s.saturating_sub(1), e, f[2], &name, f[6].bytes().next().unwrap_or(b'.'));
    }
    Ok((sink.out, sink.stats))
}

fn percent_decode(s: &str) -> String {
    if !s.contains('%') {
        return s.to_string();
    }
    let b = s.as_bytes();
    let mut out = Vec::with_capacity(b.len());
    let mut i = 0;
    while i < b.len() {
        if b[i] == b'%'
            && i + 2 < b.len()
            && let Ok(v) = u8::from_str_radix(&s[i + 1..i + 3], 16)
        {
            out.push(v);
            i += 3;
            continue;
        }
        out.push(b[i]);
        i += 1;
    }
    String::from_utf8_lossy(&out).into_owned()
}

/// RepeatMasker `.out`: kind is the class/family (e.g. `LTR/Gypsy`), name the repeat.
pub fn parse_rmsk(path: &Path, contigs: &Contigs, aliases: &Aliases) -> Result<(Vec<Feature>, ImportStats)> {
    let mut sink = Sink { aliases, contigs, out: Vec::new(), stats: ImportStats::default() };
    for line in open_lines(path)?.lines() {
        let line = line?;
        let f: Vec<&str> = line.split_whitespace().collect();
        if f.len() < 11 || f[0].parse::<f64>().is_err() {
            continue;
        }
        let (Ok(s), Ok(e)) = (f[5].parse::<u64>(), f[6].parse::<u64>()) else { continue };
        let strand = if f[8] == "C" { b'-' } else { b'+' };
        sink.push(f[4], s.saturating_sub(1), e, f[10], f[9], strand);
    }
    Ok((sink.out, sink.stats))
}

/// BED (0-based, half-open). Kind is fixed; name from column 4 (TRF: the repeat unit, col 16).
pub fn parse_bed(path: &Path, kind: &str, contigs: &Contigs, aliases: &Aliases) -> Result<(Vec<Feature>, ImportStats)> {
    let mut sink = Sink { aliases, contigs, out: Vec::new(), stats: ImportStats::default() };
    for line in open_lines(path)?.lines() {
        let line = line?;
        if line.starts_with('#') || line.starts_with("track") || line.starts_with("browser") || line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 3 {
            continue;
        }
        let (Ok(s), Ok(e)) = (f[1].parse::<u64>(), f[2].parse::<u64>()) else { continue };
        let name = if f.len() >= 16 && f[3] == "trf" { f[15] } else { f.get(3).copied().unwrap_or("") };
        let strand = f.get(5).and_then(|s| s.bytes().next()).filter(|c| *c == b'+' || *c == b'-').unwrap_or(b'.');
        sink.push(f[0], s, e, kind, name, strand);
    }
    Ok((sink.out, sink.stats))
}

/// Metadata of a stored feature set.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct FeatureSetInfo {
    pub name: String,
    pub description: String,
    pub source: String,
    pub count: u64,
    /// Kind names, indexed by kind id; with per-kind counts.
    pub kinds: Vec<String>,
    pub kind_counts: Vec<u64>,
    pub stats: ImportStats,
}

fn file(dir: &Path, set: &str, part: &str) -> std::path::PathBuf {
    dir.join(format!("feat.{set}.{part}"))
}

/// Sort and write a feature set into the index directory.
pub fn write_feature_set(
    dir: &Path,
    name: &str,
    description: &str,
    source: &str,
    mut features: Vec<Feature>,
    stats: ImportStats,
) -> Result<FeatureSetInfo> {
    features.sort_unstable();
    let mut kind_ids: HashMap<String, u16> = HashMap::new();
    let mut kinds: Vec<String> = Vec::new();
    let mut kind_counts: Vec<u64> = Vec::new();
    let mut starts = Vec::with_capacity(features.len());
    let mut ends = Vec::with_capacity(features.len());
    let mut maxend = Vec::with_capacity(features.len());
    let mut meta = Vec::with_capacity(features.len()); // kind id | strand << 16
    let mut name_off: Vec<u32> = Vec::with_capacity(features.len() + 1);
    let mut names = Vec::new();
    let mut running = 0u32;
    for f in &features {
        let id = *kind_ids.entry(f.kind.clone()).or_insert_with(|| {
            kinds.push(f.kind.clone());
            kind_counts.push(0);
            (kinds.len() - 1) as u16
        });
        kind_counts[id as usize] += 1;
        starts.push(f.start);
        ends.push(f.end);
        running = running.max(f.end);
        maxend.push(running);
        let strand = match f.strand {
            b'+' => 1u32,
            b'-' => 2,
            _ => 0,
        };
        meta.push(id as u32 | strand << 16);
        name_off.push(names.len() as u32);
        names.extend_from_slice(f.name.as_bytes());
    }
    name_off.push(names.len() as u32);
    write_array(&file(dir, name, "start"), "fstart", &starts, 0)?;
    write_array(&file(dir, name, "end"), "fend", &ends, 0)?;
    write_array(&file(dir, name, "maxend"), "fmaxend", &maxend, 0)?;
    write_array(&file(dir, name, "meta"), "fmeta", &meta, 0)?;
    write_array(&file(dir, name, "nameoff"), "fnameoff", &name_off, 0)?;
    write_array(&file(dir, name, "names"), "fnames", &names, 0)?;
    Ok(FeatureSetInfo {
        name: name.to_string(),
        description: description.to_string(),
        source: source.to_string(),
        count: features.len() as u64,
        kinds,
        kind_counts,
        stats,
    })
}

/// A stored feature, borrowed from the mapped arrays.
#[derive(Clone, Copy, Debug, serde::Serialize)]
pub struct FeatureRef<'a> {
    pub start: u32,
    pub end: u32,
    pub kind: u16,
    pub strand: char,
    pub name: &'a str,
}

pub struct FeatureSet {
    pub info: FeatureSetInfo,
    starts: Array<u32>,
    ends: Array<u32>,
    maxend: Array<u32>,
    meta: Array<u32>,
    name_off: Array<u32>,
    names: Array<u8>,
}

impl FeatureSet {
    pub fn open(dir: &Path, info: FeatureSetInfo) -> Result<Self> {
        let o = MapOptions::default();
        let n = &info.name.clone();
        Ok(Self {
            starts: Array::open(&file(dir, n, "start"), "fstart", o)?,
            ends: Array::open(&file(dir, n, "end"), "fend", o)?,
            maxend: Array::open(&file(dir, n, "maxend"), "fmaxend", o)?,
            meta: Array::open(&file(dir, n, "meta"), "fmeta", o)?,
            name_off: Array::open(&file(dir, n, "nameoff"), "fnameoff", o)?,
            names: Array::open(&file(dir, n, "names"), "fnames", o)?,
            info,
        })
    }

    pub fn len(&self) -> usize {
        self.starts.len()
    }

    pub fn is_empty(&self) -> bool {
        self.starts.is_empty()
    }

    pub fn get(&self, i: usize) -> FeatureRef<'_> {
        let m = self.meta[i];
        let (a, b) = (self.name_off[i] as usize, self.name_off[i + 1] as usize);
        FeatureRef {
            start: self.starts[i],
            end: self.ends[i],
            kind: m as u16,
            strand: ['.', '+', '-'][(m >> 16) as usize & 3],
            name: std::str::from_utf8(&self.names[a..b]).unwrap_or(""),
        }
    }

    /// Indices of features overlapping `[start, end)`, in start order.
    pub fn overlapping(&self, start: u32, end: u32) -> impl Iterator<Item = usize> + '_ {
        let first = self.maxend.partition_point(|&m| m <= start);
        let last = self.starts.partition_point(|&s| s < end);
        (first..last.max(first)).filter(move |&i| self.ends[i] > start)
    }

    /// Features per pixel column over `[start, end)` split into `width` columns, optionally
    /// only of the given kinds. A feature counts in every column it overlaps.
    pub fn density(&self, start: u32, end: u32, width: usize, kinds: Option<&[u16]>) -> Vec<u32> {
        let mut out = vec![0u32; width];
        let span = (end - start) as f64;
        for i in self.overlapping(start, end) {
            let m = self.meta[i] as u16;
            if kinds.is_some_and(|k| !k.contains(&m)) {
                continue;
            }
            let s = self.starts[i].max(start) - start;
            let e = self.ends[i].min(end) - start;
            let c0 = (s as f64 / span * width as f64) as usize;
            let c1 = (((e as f64) / span * width as f64).ceil() as usize).clamp(c0 + 1, width);
            for c in &mut out[c0.min(width - 1)..c1] {
                *c += 1;
            }
        }
        out
    }

    /// Indices of features containing global position `pos`.
    pub fn at(&self, pos: u32) -> impl Iterator<Item = usize> + '_ {
        self.overlapping(pos, pos + 1)
    }
}

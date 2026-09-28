//! Index metadata: `manifest.json` and `contigs.json` (design §6.5).

use std::collections::BTreeMap;
use std::path::Path;

use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use sm_core::Contigs;

pub const MANIFEST: &str = "manifest.json";
pub const CONTIGS: &str = "contigs.json";

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Manifest {
    pub format_version: u32,
    pub k: u32,
    /// Prefix length P in bases (design §6.1).
    pub prefix_len: u32,
    /// Width of the stored words: 32 (suffix only, P = 15) or 64 (whole k-mer).
    pub word_bits: u32,
    /// Distinct canonical k-mers, N.
    pub distinct_kmers: u64,
    /// Occurrences (valid windows).
    pub sites: u64,
    pub genome_len: u64,
    pub contigs: usize,
    pub has_copy_b: bool,
    pub source: Source,
    pub build: BuildInfo,
    /// Per-file size and CRC-32.
    pub files: BTreeMap<String, FileInfo>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Source {
    pub fasta: String,
    pub fasta_bytes: u64,
    #[serde(default)]
    pub name: String,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct BuildInfo {
    pub tool_version: String,
    pub started_unix: u64,
    pub seconds: f64,
    pub partition_bits: u32,
    pub threads: usize,
    #[serde(default)]
    pub stage_seconds: BTreeMap<String, f64>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FileInfo {
    pub bytes: u64,
    pub crc32: u32,
}

#[derive(Serialize, Deserialize)]
struct ContigRow {
    name: String,
    description: String,
    offset: u32,
    len: u32,
}

pub fn write_json<T: Serialize>(path: &Path, v: &T) -> Result<()> {
    let s = serde_json::to_string_pretty(v)?;
    std::fs::write(path, s).with_context(|| format!("writing {}", path.display()))
}

pub fn read_manifest(dir: &Path) -> Result<Manifest> {
    let p = dir.join(MANIFEST);
    let s = std::fs::read_to_string(&p).with_context(|| format!("reading {}", p.display()))?;
    serde_json::from_str(&s).with_context(|| format!("parsing {}", p.display()))
}

pub fn write_contigs(path: &Path, contigs: &Contigs) -> Result<()> {
    let rows: Vec<ContigRow> = contigs
        .iter()
        .map(|c| ContigRow { name: c.name.clone(), description: c.description.clone(), offset: c.offset, len: c.len })
        .collect();
    write_json(path, &rows)
}

pub fn read_contigs(path: &Path) -> Result<Contigs> {
    let s = std::fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
    let rows: Vec<ContigRow> = serde_json::from_str(&s)?;
    let mut c = Contigs::new();
    for r in rows {
        let added = c.push(r.name, r.description, r.len as u64)?;
        anyhow::ensure!(added.offset == r.offset, "contig {} offset mismatch", added.name);
    }
    Ok(c)
}

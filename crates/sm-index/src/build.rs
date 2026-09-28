//! Partitioned external-memory index builder (design §13).
//!
//! Pass 1 streams the FASTA once and appends every occurrence `(canonical k-mer, position,
//! strand)` to one of 2^b partition files chosen by the k-mer's top bits. Pass 2 loads each
//! partition (in k-mer order), sorts it, and appends copy A's k-mers and the positions tether.
//! Finalisation picks P from the exact distinct count, writes the prefix table, and records
//! checksums. Peak RAM is about `threads × partition size × 2`.

use std::collections::BTreeMap;
use std::fs::{self, File, OpenOptions};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::time::{Instant, SystemTime, UNIX_EPOCH};

use anyhow::{Context, Result, ensure};
use rayon::prelude::*;
use sm_core::K;
use sm_core::fasta::{FastaReader, open};
use sm_core::kmer::windows_both;
use tracing::info;

use crate::array::{Array, ArrayWriter, MapOptions, write_array};
use crate::bits::{BitWriter, SelectBuilder};
use crate::manifest::{self, BuildInfo, FileInfo, Manifest, Source};
use crate::seq::SeqWriter;

#[derive(Clone, Debug)]
pub struct BuildOptions {
    /// Prefix length P; chosen from N when `None` (design §6.1).
    pub prefix_len: Option<u32>,
    /// Where partition files go (default: inside the output directory).
    pub tmp_dir: Option<PathBuf>,
    /// Target occurrences per partition, used to choose the partition count.
    pub partition_records: u64,
    /// Force the number of partition bits (even, 2..=12).
    pub partition_bits: Option<u32>,
    /// Display name of the genome.
    pub name: String,
}

impl Default for BuildOptions {
    fn default() -> Self {
        Self {
            prefix_len: None,
            tmp_dir: None,
            partition_records: 2_000_000,
            partition_bits: None,
            name: String::new(),
        }
    }
}

/// `P = clamp(round(log₄ N), 8, 15)`.
pub fn choose_prefix_len(distinct: u64) -> u32 {
    let p = ((distinct.max(1)) as f64).log(4.0).round() as u32;
    p.clamp(8, 15)
}

fn choose_partition_bits(est_sites: u64, target: u64) -> u32 {
    let mut b = 2;
    while b < 12 && est_sites >> b > target {
        b += 2;
    }
    b
}

const REC: usize = 12;
const FLUSH: usize = REC * 8192;

/// Pass-1 sink: buffers records per partition and appends them to files.
struct Partitioner {
    dir: PathBuf,
    bits: u32,
    bufs: Vec<Vec<u8>>,
    sizes: Vec<u64>,
}

impl Partitioner {
    fn new(dir: &Path, bits: u32) -> Self {
        let n = 1usize << bits;
        Self { dir: dir.to_path_buf(), bits, bufs: vec![Vec::new(); n], sizes: vec![0; n] }
    }

    fn path(&self, pid: usize) -> PathBuf {
        self.dir.join(format!("part{pid:05}"))
    }

    #[inline]
    fn push(&mut self, canon: u64, pos: u32, fwd_is_canonical: bool) -> Result<()> {
        let pid = (canon >> (2 * K - self.bits)) as usize;
        let key = canon | ((fwd_is_canonical as u64) << 62);
        let buf = &mut self.bufs[pid];
        buf.extend_from_slice(&key.to_le_bytes());
        buf.extend_from_slice(&pos.to_le_bytes());
        if buf.len() >= FLUSH {
            self.flush(pid)?;
        }
        Ok(())
    }

    fn flush(&mut self, pid: usize) -> Result<()> {
        if self.bufs[pid].is_empty() {
            return Ok(());
        }
        let path = self.path(pid);
        let mut f = OpenOptions::new().create(true).append(true).open(&path)?;
        f.write_all(&self.bufs[pid])?;
        self.sizes[pid] += (self.bufs[pid].len() / REC) as u64;
        self.bufs[pid].clear();
        Ok(())
    }

    fn finish(mut self) -> Result<Vec<u64>> {
        for pid in 0..self.bufs.len() {
            self.flush(pid)?;
        }
        Ok(self.sizes)
    }
}

/// Sorted output of one partition.
struct PartOut {
    kmers: Vec<u64>,
    sites: Vec<u32>,
    /// Bit 0: forward strand holds the canonical k-mer. Bit 1: first site of its row.
    flags: Vec<u8>,
}

fn process_partition(path: &Path, expected: u64) -> Result<PartOut> {
    let mut bytes = Vec::with_capacity(expected as usize * REC);
    if expected > 0 {
        File::open(path)?.read_to_end(&mut bytes)?;
        fs::remove_file(path).ok();
    }
    ensure!(bytes.len() as u64 == expected * REC as u64, "{}: size mismatch", path.display());
    // Sort key: k-mer, then position; the strand bit rides along in bit 0.
    let mut recs: Vec<u128> = bytes
        .as_chunks::<REC>()
        .0
        .iter()
        .map(|r| {
            let key = u64::from_le_bytes(r[..8].try_into().unwrap());
            let pos = u32::from_le_bytes(r[8..].try_into().unwrap());
            let kmer = key & sm_core::kmer::KMER_MASK;
            let strand = key >> 62 & 1;
            ((kmer as u128) << 33) | ((pos as u128) << 1) | strand as u128
        })
        .collect();
    drop(bytes);
    recs.sort_unstable();

    let mut out =
        PartOut { kmers: Vec::new(), sites: Vec::with_capacity(recs.len()), flags: Vec::with_capacity(recs.len()) };
    let mut prev = u64::MAX;
    for r in recs {
        let kmer = (r >> 33) as u64;
        let pos = (r >> 1) as u32;
        let new_row = kmer != prev;
        if new_row {
            out.kmers.push(kmer);
            prev = kmer;
        }
        out.sites.push(pos);
        out.flags.push((r & 1) as u8 | (new_row as u8) << 1);
    }
    Ok(out)
}

fn crc32_file(path: &Path) -> Result<FileInfo> {
    let mut f = File::open(path)?;
    let mut h = crc32fast::Hasher::new();
    let mut buf = vec![0u8; 1 << 22];
    let mut bytes = 0u64;
    loop {
        let n = f.read(&mut buf)?;
        if n == 0 {
            break;
        }
        h.update(&buf[..n]);
        bytes += n as u64;
    }
    Ok(FileInfo { bytes, crc32: h.finalize() })
}

/// Build an index for `fasta` into directory `out`.
pub fn build(fasta: &Path, out: &Path, opts: &BuildOptions) -> Result<Manifest> {
    let t0 = Instant::now();
    let started_unix = SystemTime::now().duration_since(UNIX_EPOCH).map(|d| d.as_secs()).unwrap_or(0);
    let mut stages = BTreeMap::new();
    fs::create_dir_all(out).with_context(|| format!("creating {}", out.display()))?;
    let fasta_bytes = fs::metadata(fasta).with_context(|| format!("reading {}", fasta.display()))?.len();
    let gz = fasta.extension().is_some_and(|e| e == "gz");
    let est_sites = if gz { fasta_bytes * 4 } else { fasta_bytes };
    let bits = opts.partition_bits.unwrap_or_else(|| choose_partition_bits(est_sites, opts.partition_records));
    ensure!(bits.is_multiple_of(2) && (2..=12).contains(&bits), "partition bits must be even, 2..=12");
    let tmp_root = opts.tmp_dir.clone().unwrap_or_else(|| out.to_path_buf());
    fs::create_dir_all(&tmp_root)?;
    let tmp = tempfile::Builder::new().prefix(".sentromap-build-").tempdir_in(&tmp_root)?;

    // Pass 1: partition occurrences.
    let t = Instant::now();
    let mut parts = Partitioner::new(tmp.path(), bits);
    let mut contigs = sm_core::Contigs::new();
    let mut reader = FastaReader::new(open(fasta)?);
    let mut seq = Vec::new();
    let mut windows = 0u64;
    let mut seq_w = SeqWriter::create(out)?;
    while let Some(h) = reader.next_record(&mut seq)? {
        let offset = contigs.push(h.name, h.description, seq.len() as u64)?.offset;
        seq_w.extend(&seq)?;
        for (i, fwd, rc) in windows_both(&seq) {
            let (canon, fwd_is_canonical) = if fwd < rc { (fwd, true) } else { (rc, false) };
            parts.push(canon, offset + i as u32, fwd_is_canonical)?;
            windows += 1;
        }
    }
    let sizes = parts.finish()?;
    seq_w.finish()?;
    ensure!(sizes.iter().sum::<u64>() == windows, "partition sizes do not add up to the window count");
    stages.insert("pass1_partition".into(), t.elapsed().as_secs_f64());
    info!(windows, contigs = contigs.len(), partitions = sizes.len(), secs = t.elapsed().as_secs_f64(), "pass 1 done");

    // Pass 2: sort partitions in parallel batches, append outputs in order.
    let t = Instant::now();
    let mut kmers_w: ArrayWriter<u64> = ArrayWriter::create(&out.join("a.words.tmp"), "words")?;
    let mut sites_w: ArrayWriter<u32> = ArrayWriter::create(&out.join("pos.sites"), "sites")?;
    let mut strand_w = BitWriter::create(&out.join("pos.strand"), "strand")?;
    let mut row_w = BitWriter::create(&out.join("pos.rowstart"), "rowstart")?;
    let mut select = SelectBuilder::default();
    let batch = rayon::current_num_threads().max(1);
    let pids: Vec<usize> = (0..sizes.len()).collect();
    for chunk in pids.chunks(batch) {
        let outs: Vec<PartOut> = chunk
            .par_iter()
            .map(|&pid| process_partition(&parts_path(tmp.path(), pid), sizes[pid]))
            .collect::<Result<_>>()?;
        for o in outs {
            kmers_w.extend(&o.kmers)?;
            sites_w.extend(&o.sites)?;
            for f in o.flags {
                strand_w.push(f & 1 == 1)?;
                let start = f & 2 == 2;
                select.observe(row_w.len(), start);
                row_w.push(start)?;
            }
        }
    }
    // Sentinel so that row l ends at select1(l + 1).
    select.observe(row_w.len(), true);
    row_w.push(true)?;
    let distinct = kmers_w.len();
    let sites = sites_w.len();
    kmers_w.finish()?;
    sites_w.finish()?;
    strand_w.finish()?;
    row_w.finish()?;
    write_array(&out.join("pos.select"), "select", &select.samples, select.ones())?;
    ensure!(select.ones() == distinct + 1, "row count mismatch");
    drop(tmp);
    stages.insert("pass2_sort".into(), t.elapsed().as_secs_f64());
    info!(distinct, sites, secs = t.elapsed().as_secs_f64(), "pass 2 done");

    // Finalise copy A: prefix table and word width from the exact N.
    let t = Instant::now();
    let p = opts.prefix_len.unwrap_or_else(|| choose_prefix_len(distinct));
    ensure!((1..=15).contains(&p), "prefix length must be 1..=15");
    let word_bits = if p == 15 { 32 } else { 64 };
    finalise_sorted(out, "a", p, word_bits)?;
    stages.insert("finalise_a".into(), t.elapsed().as_secs_f64());

    manifest::write_contigs(&out.join(manifest::CONTIGS), &contigs)?;

    let t = Instant::now();
    let names = [
        "a.prefix",
        "a.words",
        "pos.sites",
        "pos.strand",
        "pos.rowstart",
        "pos.select",
        "seq.packed",
        "seq.other",
        "seq.soft",
        manifest::CONTIGS,
    ];
    let files: BTreeMap<String, FileInfo> =
        names.par_iter().map(|n| Ok((n.to_string(), crc32_file(&out.join(n))?))).collect::<Result<_>>()?;
    stages.insert("checksums".into(), t.elapsed().as_secs_f64());

    let m = Manifest {
        format_version: crate::array::FORMAT_VERSION,
        k: K,
        prefix_len: p,
        word_bits,
        distinct_kmers: distinct,
        sites,
        genome_len: contigs.total_len(),
        contigs: contigs.len(),
        has_copy_b: false,
        source: Source {
            fasta: fasta.canonicalize().unwrap_or(fasta.to_path_buf()).display().to_string(),
            fasta_bytes,
            name: opts.name.clone(),
        },
        build: BuildInfo {
            tool_version: env!("CARGO_PKG_VERSION").to_string(),
            started_unix,
            seconds: t0.elapsed().as_secs_f64(),
            partition_bits: bits,
            threads: rayon::current_num_threads(),
            stage_seconds: stages,
        },
        files,
    };
    manifest::write_json(&out.join(manifest::MANIFEST), &m)?;
    info!(secs = m.build.seconds, p, distinct, sites, "build done");
    Ok(m)
}

fn parts_path(dir: &Path, pid: usize) -> PathBuf {
    dir.join(format!("part{pid:05}"))
}

/// Turn `{name}.words.tmp` (sorted `u64` k-mers) into `{name}.prefix` and `{name}.words`.
pub(crate) fn finalise_sorted(out: &Path, name: &str, p: u32, word_bits: u32) -> Result<()> {
    let tmp_path = out.join(format!("{name}.words.tmp"));
    let words_path = out.join(format!("{name}.words"));
    let kmers: Array<u64> = Array::open(&tmp_path, "words", MapOptions::default())?;
    kmers.advise(memmap2::Advice::Sequential);
    let shift = 2 * (K - p);
    let buckets = 1u64 << (2 * p);
    let mut table: ArrayWriter<u32> = ArrayWriter::create(&out.join(format!("{name}.prefix")), "prefix")?;
    let mut next = 0u64;
    ensure!(kmers.len() < u32::MAX as usize, "too many distinct k-mers for a u32 prefix table");
    for (i, &k) in kmers.iter().enumerate() {
        let b = k >> shift;
        while next <= b {
            table.push(i as u32)?;
            next += 1;
        }
    }
    while next <= buckets {
        table.push(kmers.len() as u32)?;
        next += 1;
    }
    table.finish()?;
    if word_bits == 32 {
        let mut w: ArrayWriter<u32> = ArrayWriter::create(&words_path, "words")?;
        for &k in kmers.iter() {
            w.push(k as u32)?;
        }
        w.finish()?;
        drop(kmers);
        fs::remove_file(&tmp_path)?;
    } else {
        drop(kmers);
        fs::rename(&tmp_path, &words_path)?;
    }
    Ok(())
}

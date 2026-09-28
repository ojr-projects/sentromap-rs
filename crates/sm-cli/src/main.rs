use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::time::Instant;

use anyhow::{Context, Result, bail};
use clap::{Parser, Subcommand};
use sm_core::kmer::{self, K};
use sm_index::{BuildOptions, Index};

#[derive(Parser)]
#[command(name = "sentromap", version, about = "Genome-wide Hamming-neighbourhood search for 31-mers")]
struct Cli {
    /// Worker threads (default: all cores).
    #[arg(long, global = true)]
    threads: Option<usize>,
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Build an index from a FASTA file (plain or gzip).
    Build {
        fasta: PathBuf,
        /// Output index directory.
        #[arg(short, long)]
        out: PathBuf,
        /// Prefix length P (default: round(log4 N) clamped to 8..=15).
        #[arg(long)]
        prefix_len: Option<u32>,
        /// Directory for temporary partition files (default: inside the output directory).
        #[arg(long)]
        tmp_dir: Option<PathBuf>,
        /// Display name of the genome.
        #[arg(long, default_value = "")]
        name: String,
        /// Spot-check this many random windows after building.
        #[arg(long, default_value_t = 10_000)]
        spot_checks: usize,
    },
    /// Search for every k-mer within n substitutions of a query, on both strands.
    Search {
        index: PathBuf,
        /// 31 bases, or CONTIG:POS (1-based) to take the k-mer from the genome.
        query: String,
        /// Maximum substitutions.
        #[arg(short, default_value_t = 2)]
        n: u32,
        /// Output: summary, variants (one row per variant) or sites (one row per site).
        #[arg(long, default_value = "summary")]
        output: String,
    },
    /// Validate an index: structure, checksums, and (with --fasta) spot checks.
    Verify {
        index: PathBuf,
        #[arg(long)]
        fasta: Option<PathBuf>,
        #[arg(long, default_value_t = 10_000)]
        spot_checks: usize,
    },
    /// Print an index's manifest.
    Info { index: PathBuf },
}

fn peak_rss_mb() -> Option<f64> {
    let s = std::fs::read_to_string("/proc/self/status").ok()?;
    let line = s.lines().find(|l| l.starts_with("VmHWM:"))?;
    let kb: f64 = line.split_whitespace().nth(1)?.parse().ok()?;
    Some(kb / 1024.0)
}

fn parse_query(idx: &Index, q: &str) -> Result<u64> {
    if let Some((contig, pos)) = q.rsplit_once(':') {
        let ci = idx.contigs.by_name(contig).with_context(|| format!("no contig named {contig:?}"))?;
        let pos: u32 = pos.replace(',', "").parse().context("position must be a number")?;
        if pos == 0 {
            bail!("positions are 1-based");
        }
        return idx.kmer_at(ci, pos - 1).with_context(|| format!("no valid {K}-mer at {q}"));
    }
    Ok(kmer::parse(q)?)
}

fn main() -> Result<()> {
    tracing_subscriber::fmt()
        .with_env_filter(tracing_subscriber::EnvFilter::try_from_default_env().unwrap_or_else(|_| "info".into()))
        .with_writer(std::io::stderr)
        .init();
    let cli = Cli::parse();
    if let Some(t) = cli.threads {
        rayon::ThreadPoolBuilder::new().num_threads(t).build_global()?;
    }
    match cli.cmd {
        Cmd::Build { fasta, out, prefix_len, tmp_dir, name, spot_checks } => {
            let opts = BuildOptions { prefix_len, tmp_dir, name, ..Default::default() };
            let m = sm_index::build(&fasta, &out, &opts)?;
            let idx = Index::open(&out)?;
            let t = Instant::now();
            sm_index::verify::verify_structure(&idx)?;
            let checked = sm_index::verify::spot_check(&idx, &fasta, spot_checks, 1)?;
            let bytes: u64 = m.files.values().map(|f| f.bytes).sum();
            eprintln!(
                "built {}: {} contigs, {} bases, {} sites, N = {} distinct, P = {}, {:.2} GB on disk, {:.1} s; verified ({} spot checks, {:.1} s); peak RSS {:.0} MB",
                out.display(),
                m.contigs,
                m.genome_len,
                m.sites,
                m.distinct_kmers,
                m.prefix_len,
                bytes as f64 / 1e9,
                m.build.seconds,
                checked,
                t.elapsed().as_secs_f64(),
                peak_rss_mb().unwrap_or(0.0)
            );
        }
        Cmd::Search { index, query, n, output } => {
            let idx = Index::open(&index)?;
            let q = parse_query(&idx, &query)?;
            let t = Instant::now();
            let v = sm_search::scan(&idx, q, n);
            let secs = t.elapsed().as_secs_f64();
            let mut out = BufWriter::new(std::io::stdout().lock());
            match output.as_str() {
                "summary" => {
                    let sites = v.site_count(&idx);
                    writeln!(out, "query\t{}", kmer::to_string(q))?;
                    writeln!(out, "variants\t{}\nsites\t{sites}\nsearch_ms\t{:.2}", v.len(), secs * 1e3)?;
                    let hist = v.histogram();
                    let mut cum_v = 0;
                    let mut cum_s = 0u64;
                    let mut per_m_sites = vec![0u64; n as usize + 1];
                    for i in 0..v.len() {
                        per_m_sites[v.mismatches[i] as usize] += idx.row_len(v.leaf[i]);
                    }
                    writeln!(out, "n\tvariants\tsites\tcum_variants\tcum_sites\texpected_chance_hits")?;
                    for m in 0..=n as usize {
                        cum_v += hist[m];
                        cum_s += per_m_sites[m];
                        let e = sm_core::ball::expected_chance_hits(idx.manifest.distinct_kmers, m as u32);
                        writeln!(out, "{m}\t{}\t{}\t{cum_v}\t{cum_s}\t{e:.3e}", hist[m], per_m_sites[m])?;
                    }
                }
                "variants" => {
                    writeln!(out, "variant\tmismatches\tpositions\tsites")?;
                    for i in 0..v.len() {
                        let pos: Vec<String> = kmer::mismatch_positions(q, v.kmer[i]).map(|p| p.to_string()).collect();
                        writeln!(
                            out,
                            "{}\t{}\t{}\t{}",
                            kmer::to_string(v.kmer[i]),
                            v.mismatches[i],
                            pos.join(","),
                            idx.row_len(v.leaf[i])
                        )?;
                    }
                }
                "sites" => {
                    writeln!(out, "contig\tstart\tend\tstrand\tvariant\tmismatches")?;
                    for h in v.sites(&idx) {
                        let (ci, off) = idx.contigs.to_local(h.pos).unwrap();
                        writeln!(
                            out,
                            "{}\t{}\t{}\t{}\t{}\t{}",
                            idx.contigs.get(ci).name,
                            off + 1,
                            off + K,
                            h.strand.symbol(),
                            kmer::to_string(h.variant),
                            h.mismatches
                        )?;
                    }
                }
                other => bail!("unknown output {other:?}"),
            }
        }
        Cmd::Verify { index, fasta, spot_checks } => {
            let idx = Index::open(&index)?;
            sm_index::verify::verify_structure(&idx)?;
            sm_index::verify::verify_checksums(&idx)?;
            if let Some(fa) = fasta {
                let n = sm_index::verify::spot_check(&idx, &fa, spot_checks, 1)?;
                eprintln!("{n} spot checks passed");
            }
            eprintln!("ok");
        }
        Cmd::Info { index } => {
            let m = sm_index::manifest::read_manifest(&index)?;
            println!("{}", serde_json::to_string_pretty(&m)?);
        }
    }
    Ok(())
}

use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::time::Instant;

use anyhow::{Context, Result, bail};
use clap::{Parser, Subcommand};
use sm_core::kmer::{self, K};
use sm_index::{BuildOptions, Index};
use sm_search::{CostModel, Engine};

mod bench;

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
        /// Skip copy B (full-scan-only index, about half the size).
        #[arg(long)]
        no_copy_b: bool,
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
        /// Engine: auto (cost model), scan or pigeonhole.
        #[arg(long, default_value = "auto")]
        engine: String,
    },
    /// Benchmarks.
    Bench {
        #[command(subcommand)]
        what: BenchCmd,
    },
    /// Validate an index: structure, checksums, and (with --fasta) spot checks.
    Verify {
        index: PathBuf,
        #[arg(long)]
        fasta: Option<PathBuf>,
        #[arg(long, default_value_t = 10_000)]
        spot_checks: usize,
    },
    /// Import annotations into an index: GFF3 features, a tracks.tsv manifest (bigWig signal,
    /// RepeatMasker .out, BED), chromosome aliases, and a GC% track from the genome.
    Annotate {
        index: PathBuf,
        /// GFF3 file(s), imported as the feature set "genes" (or NAME=PATH).
        #[arg(long)]
        gff: Vec<String>,
        /// Track manifest: name, kind, format (bigwig|rmsk-out|bed), file, source, description.
        #[arg(long)]
        tracks: Option<PathBuf>,
        /// Chromosome alias table (UCSC chromAlias.txt style).
        #[arg(long)]
        alias: Vec<PathBuf>,
        /// Skip the GC% track.
        #[arg(long)]
        no_gc: bool,
    },
    /// Serve the HTTP API (and the web front end) for an index.
    Serve {
        index: PathBuf,
        #[arg(long, default_value = "127.0.0.1:8080")]
        addr: std::net::SocketAddr,
        /// Static front-end directory served at / (default: web/dist if present).
        #[arg(long)]
        web: Option<PathBuf>,
        /// Result cache budget in MB.
        #[arg(long, default_value_t = 4096)]
        cache_mb: usize,
        /// Skip pre-faulting the index into memory.
        #[arg(long)]
        no_populate: bool,
        /// Skip the cost-model calibration at start-up.
        #[arg(long)]
        no_calibrate: bool,
        /// Allow cross-origin requests (for a front end served elsewhere).
        #[arg(long)]
        cors: bool,
    },
    /// Print an index's manifest.
    Info { index: PathBuf },
}

#[derive(Subcommand)]
enum BenchCmd {
    /// Compare the engines over n for unique, top-repeat and random queries (checks agreement).
    Search {
        index: PathBuf,
        #[arg(long, default_value_t = 15)]
        max_n: u32,
        /// Queries per class.
        #[arg(long, default_value_t = 5)]
        queries: usize,
    },
    /// Ordering cost and quality (character tree per-n and frozen, plain sort, exact MST).
    Order {
        index: PathBuf,
        /// Queries: 31 bases or CONTIG:POS. Default: the top repeat and a unique k-mer.
        #[arg(long)]
        query: Vec<String>,
        #[arg(long, value_delimiter = ',', default_value = "4,6,8,10,12,15")]
        n: Vec<u32>,
    },
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
        Cmd::Build { fasta, out, prefix_len, tmp_dir, name, no_copy_b, spot_checks } => {
            let opts = BuildOptions { prefix_len, tmp_dir, name, copy_b: !no_copy_b, ..Default::default() };
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
        Cmd::Search { index, query, n, output, engine } => {
            let idx = Index::open(&index)?;
            let q = parse_query(&idx, &query)?;
            let force = match engine.as_str() {
                "auto" => None,
                "scan" => Some(Engine::Scan),
                "pigeonhole" => Some(Engine::Pigeonhole),
                e => bail!("unknown engine {e:?}"),
            };
            if force == Some(Engine::Pigeonhole) && idx.b.is_none() {
                bail!("this index has no copy B; rebuild without --no-copy-b");
            }
            let t = Instant::now();
            let (v, plan) = sm_search::search(&idx, q, n, &CostModel::default(), force);
            let secs = t.elapsed().as_secs_f64();
            let mut out = BufWriter::new(std::io::stdout().lock());
            match output.as_str() {
                "summary" => {
                    let sites = v.site_count(&idx);
                    writeln!(out, "query\t{}", kmer::to_string(q))?;
                    writeln!(
                        out,
                        "variants\t{}\nsites\t{sites}\nengine\t{:?}\nsearch_ms\t{:.2}",
                        v.len(),
                        plan.engine,
                        secs * 1e3
                    )?;
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
        Cmd::Bench { what: BenchCmd::Search { index, max_n, queries } } => {
            let idx = Index::open_with(&index, sm_index::array::MapOptions { populate: true })?;
            bench::search_bench(&idx, max_n, queries, &mut std::io::stdout().lock())?;
        }
        Cmd::Bench { what: BenchCmd::Order { index, query, n } } => {
            let idx = Index::open_with(&index, sm_index::array::MapOptions { populate: true })?;
            let queries: Vec<(String, u64)> = if query.is_empty() {
                let sets = bench::query_sets(&idx, 1, 7);
                sets.iter().filter(|s| s.class != "random").map(|s| (s.class.to_string(), s.queries[0])).collect()
            } else {
                query.iter().map(|q| Ok((q.clone(), parse_query(&idx, q)?))).collect::<Result<_>>()?
            };
            for (label, q) in &queries {
                eprintln!("{label}: {}", kmer::to_string(*q));
            }
            bench::order_bench(&idx, &queries, &n, &mut std::io::stdout().lock())?;
        }
        Cmd::Annotate { index, gff, tracks, alias, no_gc } => {
            annotate(&index, &gff, tracks.as_deref(), &alias, !no_gc)?
        }
        Cmd::Serve { index, addr, web, cache_mb, no_populate, no_calibrate, cors } => {
            let web = web.or_else(|| {
                let d = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../web/dist");
                d.is_dir().then_some(d)
            });
            let cfg = sm_server::Config {
                index,
                addr,
                web,
                cache_bytes: cache_mb << 20,
                populate: !no_populate,
                cors,
                calibrate: !no_calibrate,
            };
            tokio::runtime::Builder::new_multi_thread().enable_all().build()?.block_on(sm_server::serve(cfg))?;
        }
        Cmd::Info { index } => {
            let m = sm_index::manifest::read_manifest(&index)?;
            println!("{}", serde_json::to_string_pretty(&m)?);
        }
    }
    Ok(())
}

fn annotate(
    dir: &std::path::Path,
    gffs: &[String],
    tracks: Option<&std::path::Path>,
    aliases: &[PathBuf],
    gc: bool,
) -> Result<()> {
    use rayon::prelude::*;
    use sm_annot::{Registry, features, tracks as tr};
    let idx = Index::open(dir)?;
    let mut al = sm_annot::Aliases::new(&idx.contigs);
    for a in aliases {
        let n = al.load(a)?;
        eprintln!("aliases: {n} rows from {}", a.display());
    }
    let mut reg = Registry::load(dir)?;
    reg.contig_display = al.display.clone();
    let report = |what: &str, s: &features::ImportStats| {
        eprintln!(
            "{what}: {} features, {} skipped, {} with unknown seqids {:?}",
            s.features, s.skipped, s.unknown_seqid, s.unknown_seqids
        );
    };
    for g in gffs {
        let (name, path) = g.split_once('=').map_or(("genes", g.as_str()), |(n, p)| (n, p));
        let (f, st) = features::parse_gff(std::path::Path::new(path), &idx.contigs, &al)?;
        report(name, &st);
        reg.upsert_features(features::write_feature_set(dir, name, "GFF3 features", path, f, st)?);
    }
    if let Some(tsv) = tracks {
        let base = tsv.parent().unwrap_or(std::path::Path::new("."));
        let text = std::fs::read_to_string(tsv)?;
        let rows: Vec<Vec<String>> = text
            .lines()
            .skip(1)
            .filter(|l| !l.trim().is_empty())
            .map(|l| l.split('\t').map(String::from).collect())
            .collect();
        let mut signal = Vec::new();
        for r in &rows {
            let [name, kind, format, file, source, desc] = &r[..] else { bail!("bad tracks.tsv row: {r:?}") };
            let path = base.join(file);
            match format.as_str() {
                "bigwig" => signal.push((name, kind, path, source, desc)),
                "rmsk-out" => {
                    let (f, st) = features::parse_rmsk(&path, &idx.contigs, &al)?;
                    report(name, &st);
                    reg.upsert_features(features::write_feature_set(dir, name, desc, source, f, st)?);
                }
                "bed" => {
                    let (f, st) = features::parse_bed(&path, name, &idx.contigs, &al)?;
                    report(name, &st);
                    reg.upsert_features(features::write_feature_set(dir, name, desc, source, f, st)?);
                }
                other => eprintln!("skipping {name}: unsupported format {other:?}"),
            }
        }
        let t = Instant::now();
        let infos: Vec<sm_annot::TrackInfo> = signal
            .par_iter()
            .map(|(name, kind, path, source, desc)| {
                let meta = tr::TrackMeta { name, kind, description: desc, source };
                tr::import_bigwig(dir, path, &meta, &idx.contigs, &al)
            })
            .collect::<Result<_>>()?;
        for i in infos {
            eprintln!(
                "track {}: coverage {:.1}%, p01 {:.3}, p99 {:.3}, unknown chroms {}",
                i.name,
                i.coverage * 100.0,
                i.p01,
                i.p99,
                i.unknown_chroms.len()
            );
            reg.upsert_track(i);
        }
        eprintln!("signal tracks imported in {:.1} s", t.elapsed().as_secs_f64());
    }
    if gc {
        reg.upsert_track(tr::import_gc(dir, &idx.seq, &idx.contigs)?);
    }
    reg.save(dir)?;
    eprintln!("registry: {} feature sets, {} tracks", reg.feature_sets.len(), reg.tracks.len());
    Ok(())
}

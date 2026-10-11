//! Doom checks (docs/doom.md): cheap tests designed to kill the companion-track idea.
//!
//! Kill test zero: do crease-aligned profiles of n-neighbourhoods beat a circular-shift null,
//! for any query class other than TE copies?

use std::io::Write;

use anyhow::{Context, Result, bail};
use rayon::prelude::*;
use sm_annot::{Annotations, Track};
use sm_core::kmer::{self, K};
use sm_core::synth::Rng;
use sm_core::{Hit, Strand};
use sm_index::Index;
use sm_search::CostModel;

/// Half-width of the profile window and its sampling step (the 32 bp base level of the tracks).
const HALF: i64 = 2048;
const STEP: i64 = 32;
/// Central samples: |offset| ≤ 96 bp from the k-mer centre. Flanks: |offset| ≥ 1024 bp.
const CENTRE: i64 = 96;
const FLANK: i64 = 1024;
const MAX_SITES: usize = 500;
const SHIFTS: usize = 100;
const NS: [u8; 4] = [0, 4, 8, 10];
const TRACKS: [&str; 4] = ["CTCF", "H3K9me3", "phastCons124way", "Input"];
/// The D. melanogaster CTCF core motif (CCAGCAGGGGGCG, as in vertebrate CTCF).
const CTCF_MOTIF: &[u8] = b"CCAGCAGGGGGCG";

#[derive(Clone)]
struct Query {
    class: &'static str,
    kmer: u64,
    /// The query's own genome position (none for absent k-mers).
    own: Option<u32>,
}

fn genomic_query(idx: &Index, class: &'static str, pos: u32) -> Option<Query> {
    let (ci, off) = idx.contigs.to_local(pos)?;
    let k = idx.kmer_at(ci, off)?;
    Some(Query { class, kmer: k, own: Some(pos) })
}

fn queries(idx: &Index, ann: &Annotations, per_class: usize, seed: u64, variant: u8) -> Result<Vec<Query>> {
    let mut rng = Rng::new(seed);
    let len = idx.manifest.genome_len;
    let major = |pos: u32| idx.contigs.locate(pos).is_some_and(|c| idx.contigs.get(c).len >= 1_000_000);
    let mut out = Vec::new();

    // Unique: single-copy 31-mers on the major chromosome arms.
    let mut n = 0;
    while n < per_class {
        let p = (rng.next_u64() % len) as u32;
        if !major(p) {
            continue;
        }
        if let Some(q) = genomic_query(idx, "unique", p)
            && idx.a.lookup(sm_core::canonical(q.kmer)).is_some_and(|l| idx.row_len(l) == 1)
        {
            out.push(q);
            n += 1;
        }
    }

    // TE: positions inside RepeatMasker LTR / LINE / DNA elements of at least 300 bp.
    let rm = ann.feature_set("RepeatMasker").context("RepeatMasker feature set missing; run `sentromap annotate`")?;
    let tes: Vec<usize> = (0..rm.len())
        .filter(|&i| {
            let f = rm.get(i);
            let class = rm.info.kinds[f.kind as usize].split('/').next().unwrap_or("");
            matches!(class, "LTR" | "LINE" | "DNA") && f.end - f.start >= 300
        })
        .collect();
    n = 0;
    while n < per_class {
        let f = rm.get(tes[rng.range(0, tes.len())]);
        let p = f.start + rng.range(0, (f.end - f.start - K) as usize) as u32;
        if let Some(q) = genomic_query(idx, "te", p) {
            out.push(q);
            n += 1;
        }
    }

    // CTCF motif: 31-mers centred on exact motif occurrences (either strand) on major contigs.
    let rc_motif = sm_core::synth::revcomp_seq(CTCF_MOTIF);
    let mut hits = Vec::new();
    for c in idx.contigs.iter().filter(|c| c.len >= 1_000_000) {
        let seq = idx.seq.fetch(c.offset, c.end()).to_ascii_uppercase();
        for (i, w) in seq.windows(CTCF_MOTIF.len()).enumerate() {
            let mm = |m: &[u8]| w.iter().zip(m).filter(|(a, b)| a != b).count();
            if (variant == 0 && (w == CTCF_MOTIF || w == rc_motif.as_slice()))
                || (variant == 1 && (mm(CTCF_MOTIF) <= 1 || mm(&rc_motif) <= 1))
            {
                hits.push(c.offset + i as u32);
            }
        }
    }
    if hits.is_empty() {
        bail!("no CTCF motif occurrences found");
    }
    n = 0;
    let mut tries = 0;
    while n < per_class.min(hits.len()) && tries < 100 * per_class {
        tries += 1;
        let m = hits[rng.range(0, hits.len())];
        let p = (m + CTCF_MOTIF.len() as u32 / 2).saturating_sub(K / 2);
        if let Some(q) = genomic_query(idx, "ctcf_motif", p)
            && !out.iter().any(|o| o.own == q.own)
        {
            out.push(q);
            n += 1;
        }
    }

    if variant == 1 {
        // CTCF summits: strongest local maxima (±512 bp) of the CTCF track on major contigs.
        let track = ann.track("CTCF").context("CTCF track missing")?;
        let mut peaks: Vec<(f32, u32)> = Vec::new();
        for c in idx.contigs.iter().filter(|c| c.len >= 1_000_000) {
            let (b0, b1) = (c.offset / 32 + 16, (c.end() / 32).saturating_sub(16));
            for b in b0..b1 {
                let v = track.value_at(b * 32);
                if v.is_nan() {
                    continue;
                }
                if (b - 16..=b + 16)
                    .all(|o| o == b || track.value_at(o * 32).partial_cmp(&v) != Some(std::cmp::Ordering::Greater))
                {
                    peaks.push((v, b * 32 + 16));
                }
            }
        }
        peaks.sort_by(|a, b| b.0.total_cmp(&a.0));
        peaks.truncate(1000);
        n = 0;
        let mut tries = 0;
        while n < per_class && tries < 100 * per_class && !peaks.is_empty() {
            tries += 1;
            let (_, centre) = peaks[rng.range(0, peaks.len())];
            if let Some(q) = genomic_query(idx, "ctcf_summit", centre.saturating_sub(K / 2))
                && !out.iter().any(|o| o.own == q.own)
            {
                out.push(q);
                n += 1;
            }
        }
        // Composition-matched twins: each genomic query's bases shuffled, absent from the genome.
        let genomic: Vec<Query> = out.clone();
        for q in genomic {
            let mut bases: Vec<u8> = (0..K).map(|i| kmer::base_at(q.kmer, i)).collect();
            for _ in 0..1000 {
                for i in (1..bases.len()).rev() {
                    let j = rng.range(0, i + 1);
                    bases.swap(i, j);
                }
                let k = bases.iter().fold(0u64, |a, &b| (a << 2) | b as u64);
                if idx.a.lookup(sm_core::canonical(k)).is_none() {
                    out.push(Query { class: shuffled_class(q.class), kmer: k, own: None });
                    break;
                }
            }
        }
    }

    // Absent: random 31-mers not in the genome (negative control: their neighbours are chance sites).
    n = 0;
    while n < per_class {
        let k = rng.kmer();
        if idx.a.lookup(sm_core::canonical(k)).is_none() {
            out.push(Query { class: "absent", kmer: k, own: None });
            n += 1;
        }
    }
    eprintln!("kt0: {} CTCF motif occurrences genome-wide; {} queries", hits.len(), out.len());
    Ok(out)
}

fn shuffled_class(c: &str) -> &'static str {
    match c {
        "unique" => "unique_shuf",
        "te" => "te_shuf",
        "ctcf_motif" => "ctcf_motif_shuf",
        "ctcf_summit" => "ctcf_summit_shuf",
        _ => "other_shuf",
    }
}

/// A site prepared for profiling: k-mer centre, orientation, and its contig's span.
#[derive(Clone, Copy)]
struct Site {
    centre: i64,
    reverse: bool,
    c_off: i64,
    c_len: i64,
}

/// `(shape, level)` for a set of sites shifted by `delta` (wrapping within each contig).
/// Shape: mean over sites of (centre mean − flank mean). Level: mean over sites of window mean.
fn stats(track: &Track, sites: &[Site], delta: i64) -> (f64, f64) {
    let (mut shape, mut ns, mut level, mut nl) = (0.0, 0usize, 0.0, 0usize);
    for s in sites {
        let c = s.c_off + (s.centre - s.c_off + delta).rem_euclid(s.c_len);
        let (mut cs, mut cn, mut fs, mut fnum, mut all, mut an) = (0.0, 0, 0.0, 0, 0.0, 0);
        let mut o = -HALF;
        while o <= HALF {
            let x = if s.reverse { c - o } else { c + o };
            if x >= s.c_off && x < s.c_off + s.c_len {
                let v = track.value_at(x as u32);
                if !v.is_nan() {
                    let v = v as f64;
                    all += v;
                    an += 1;
                    if o.abs() <= CENTRE {
                        cs += v;
                        cn += 1;
                    } else if o.abs() >= FLANK {
                        fs += v;
                        fnum += 1;
                    }
                }
            }
            o += STEP;
        }
        if cn > 0 && fnum > 0 {
            shape += cs / cn as f64 - fs / fnum as f64;
            ns += 1;
        }
        if an > 0 {
            level += all / an as f64;
            nl += 1;
        }
    }
    (if ns > 0 { shape / ns as f64 } else { f64::NAN }, if nl > 0 { level / nl as f64 } else { f64::NAN })
}

/// z-score and two-sided empirical p of `obs` against null draws.
fn test(obs: f64, null: &[f64]) -> (f64, f64) {
    let null: Vec<f64> = null.iter().copied().filter(|x| x.is_finite()).collect();
    if !obs.is_finite() || null.len() < 10 {
        return (f64::NAN, f64::NAN);
    }
    let m = null.iter().sum::<f64>() / null.len() as f64;
    let sd = (null.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (null.len() - 1) as f64).sqrt();
    let extreme = null.iter().filter(|x| (*x - m).abs() >= (obs - m).abs()).count();
    let z = if sd > 0.0 { (obs - m) / sd } else { 0.0 };
    (z, (extreme + 1) as f64 / (null.len() + 1) as f64)
}

struct Row {
    class: &'static str,
    kmer: u64,
    track: &'static str,
    n: u8,
    sites: usize,
    effective: usize,
    shape: f64,
    shape_z: f64,
    shape_p: f64,
    level_z: f64,
    level_p: f64,
}

fn evaluate(idx: &Index, ann: &Annotations, q: &Query, seed: u64) -> Vec<Row> {
    let (v, _) = sm_search::search(idx, q.kmer, 10, &CostModel::default(), None);
    let all: Vec<Hit> = v.sites(idx);
    let mut rows = Vec::new();
    for n in NS {
        // Leave the query's own site (and anything within the window of it) out.
        let mut hits: Vec<&Hit> = all
            .iter()
            .filter(|h| h.mismatches <= n && q.own.is_none_or(|p| (h.pos as i64 - p as i64).abs() > HALF + K as i64))
            .collect();
        let raw = hits.len();
        // One site per window: tandem neighbours share flanks (pseudoreplication).
        hits.sort_by_key(|h| h.pos);
        let mut eff: Vec<&Hit> = Vec::new();
        for h in hits {
            if eff.last().is_none_or(|l| h.pos as i64 - l.pos as i64 > 2 * HALF) {
                eff.push(h);
            }
        }
        let effective = eff.len();
        let mut rng = Rng::new(seed ^ q.kmer ^ n as u64);
        if eff.len() > MAX_SITES {
            for i in 0..MAX_SITES {
                let j = rng.range(i, eff.len());
                eff.swap(i, j);
            }
            eff.truncate(MAX_SITES);
        }
        let sites: Vec<Site> = eff
            .iter()
            .filter_map(|h| {
                let c = idx.contigs.get(idx.contigs.locate(h.pos)?);
                Some(Site {
                    centre: h.pos as i64 + K as i64 / 2,
                    reverse: h.strand == Strand::Reverse,
                    c_off: c.offset as i64,
                    c_len: c.len as i64,
                })
            })
            .collect();
        let deltas: Vec<i64> = (0..SHIFTS).map(|_| 10_000 + (rng.next_u64() % 10_000_000) as i64).collect();
        for name in TRACKS {
            let Some(track) = ann.track(name) else { continue };
            let mut row = Row {
                class: q.class,
                kmer: q.kmer,
                track: name,
                n,
                sites: raw,
                effective,
                shape: f64::NAN,
                shape_z: f64::NAN,
                shape_p: f64::NAN,
                level_z: f64::NAN,
                level_p: f64::NAN,
            };
            if sites.len() >= 3 {
                let (shape, level) = stats(track, &sites, 0);
                let null: Vec<(f64, f64)> = deltas.iter().map(|&d| stats(track, &sites, d)).collect();
                (row.shape_z, row.shape_p) = test(shape, &null.iter().map(|x| x.0).collect::<Vec<_>>());
                (row.level_z, row.level_p) = test(level, &null.iter().map(|x| x.1).collect::<Vec<_>>());
                row.shape = shape;
            }
            rows.push(row);
        }
    }
    rows
}

/// `variant` 0 is KT0 as first registered; 1 is KT0b (summit class, ≤1-mismatch motif,
/// composition-matched shuffled twins as each class's baseline).
pub fn kill_test_zero(
    idx: &Index,
    ann: &Annotations,
    per_class: usize,
    variant: u8,
    tsv: &std::path::Path,
    out: &mut impl Write,
) -> Result<()> {
    let seed = if variant == 0 { 0x5EED_0000 } else { 0x5EED_000B };
    let qs = queries(idx, ann, per_class, seed, variant)?;
    let t = std::time::Instant::now();
    let rows: Vec<Row> = qs.par_iter().flat_map(|q| evaluate(idx, ann, q, seed)).collect();
    eprintln!("kt0: evaluated in {:.1} s", t.elapsed().as_secs_f64());

    let mut f = std::io::BufWriter::new(std::fs::File::create(tsv)?);
    writeln!(f, "class\tkmer\ttrack\tn\tsites\teffective\tshape\tshape_z\tshape_p\tlevel_z\tlevel_p")?;
    for r in &rows {
        writeln!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.3}\t{:.4}\t{:.3}\t{:.4}",
            r.class,
            kmer::to_string(r.kmer),
            r.track,
            r.n,
            r.sites,
            r.effective,
            r.shape,
            r.shape_z,
            r.shape_p,
            r.level_z,
            r.level_p
        )?;
    }

    // Summary: fraction of testable queries significant (p ≤ 0.01; the null floor is 1/101).
    const ALPHA: f64 = 0.0101;
    writeln!(
        out,
        "class\ttrack\tn\ttestable\tmedian_eff_sites\tshape_sig\tshape_sig_frac\tmedian_shape_z\tlevel_sig_frac"
    )?;
    let classes: &[&str] = if variant == 0 {
        &["unique", "te", "ctcf_motif", "absent"]
    } else {
        &[
            "unique",
            "unique_shuf",
            "te",
            "te_shuf",
            "ctcf_motif",
            "ctcf_motif_shuf",
            "ctcf_summit",
            "ctcf_summit_shuf",
            "absent",
        ]
    };
    let mut baseline = std::collections::HashMap::new();
    let mut summary = Vec::new();
    for &class in classes {
        for track in TRACKS {
            for n in NS {
                let sel: Vec<&Row> = rows
                    .iter()
                    .filter(|r| r.class == class && r.track == track && r.n == n && r.shape_p.is_finite())
                    .collect();
                let testable = sel.len();
                let sig = sel.iter().filter(|r| r.shape_p <= ALPHA).count();
                let lsig = sel.iter().filter(|r| r.level_p <= ALPHA).count();
                let frac = |k: usize| if testable > 0 { k as f64 / testable as f64 } else { f64::NAN };
                let med = |mut v: Vec<f64>| {
                    v.sort_by(f64::total_cmp);
                    if v.is_empty() { f64::NAN } else { v[v.len() / 2] }
                };
                let med_eff = med(sel.iter().map(|r| r.effective as f64).collect());
                let med_z = med(sel.iter().map(|r| r.shape_z).collect());
                writeln!(
                    out,
                    "{class}\t{track}\t{n}\t{testable}\t{med_eff:.0}\t{sig}\t{:.3}\t{med_z:.2}\t{:.3}",
                    frac(sig),
                    frac(lsig)
                )?;
                baseline.insert((class, track, n), frac(sig));
                summary.push((class, track, n, testable, frac(sig)));
            }
        }
    }
    // Verdict (pre-registered in docs/doom.md): a non-TE genomic class beats the shifted null on
    // shape at > 3x the absent-k-mer hit rate (floored at alpha) with at least 5 significant queries.
    let mut survivors = Vec::new();
    for &(class, track, n, testable, frac) in &summary {
        if class == "te" || class == "absent" || class.ends_with("_shuf") || testable == 0 {
            continue;
        }
        if variant == 1 && track == "Input" {
            continue; // diagnostic only in KT0b
        }
        let control = if variant == 0 { "absent".to_string() } else { format!("{class}_shuf") };
        let base =
            baseline.get(&(control.as_str(), track, n)).copied().filter(|b| b.is_finite()).unwrap_or(0.0).max(ALPHA);
        let sig = (frac * testable as f64).round() as usize;
        if frac > 3.0 * base && sig >= 5 {
            let input = baseline.get(&(class, "Input", n)).copied().unwrap_or(f64::NAN);
            let flag = if variant == 1 && input >= frac { " [FLAG: Input as high: artefact?]" } else { "" };
            survivors.push(format!(
                "{class}/{track}/n={n} ({sig}/{testable}, {:.1}x {control}; Input rate {input:.2}){flag}",
                frac / base
            ));
        }
    }
    if survivors.is_empty() {
        writeln!(out, "\nVERDICT: KILL: no non-TE class beats the circular-shift null on shape")?;
    } else {
        writeln!(out, "\nVERDICT: CONTINUE: {}", survivors.join("; "))?;
    }
    Ok(())
}

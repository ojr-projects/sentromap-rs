//! Engines against the brute-force oracle on synthetic genomes (design §15.2).

use sm_core::kmer::{self, K};
use sm_core::oracle;
use sm_core::synth::{self, Rng, SynthSpec};
use sm_index::{BuildOptions, Index, build};

fn build_synth(
    spec: &SynthSpec,
    opts: &BuildOptions,
) -> (tempfile::TempDir, sm_core::fasta::Genome, Vec<Vec<u8>>, Index) {
    let dir = tempfile::tempdir().unwrap();
    let s = synth::generate(spec);
    let fa = dir.path().join("g.fa");
    s.genome.write_fasta(&mut std::fs::File::create(&fa).unwrap()).unwrap();
    build(&fa, &dir.path().join("idx"), opts).unwrap();
    let idx = Index::open(&dir.path().join("idx")).unwrap();
    (dir, s.genome, s.family_sources, idx)
}

/// Queries: family members (both orientations), genome windows, and random k-mers.
fn queries(g: &sm_core::fasta::Genome, sources: &[Vec<u8>], rng: &mut Rng) -> Vec<u64> {
    let mut qs = Vec::new();
    for src in sources {
        let at = rng.range(0, src.len() - K as usize);
        let q = kmer::parse(&String::from_utf8_lossy(&src[at..at + K as usize])).unwrap();
        qs.push(q);
        qs.push(kmer::revcomp(q));
    }
    let wins: Vec<u64> = kmer::windows(&g.seq).map(|(_, k)| k).collect();
    for _ in 0..6 {
        qs.push(wins[rng.range(0, wins.len())]);
    }
    for _ in 0..4 {
        qs.push(rng.kmer());
    }
    qs
}

fn check_engine(opts: BuildOptions, seeds: std::ops::RangeInclusive<u64>) {
    for seed in seeds {
        let spec = SynthSpec { seed, ..Default::default() };
        let (_dir, g, sources, idx) = build_synth(&spec, &opts);
        let mut rng = Rng::new(seed * 7919);
        for q in queries(&g, &sources, &mut rng) {
            for n in [0, 1, 2, 3, 4, 6] {
                let got = sm_search::scan(&idx, q, n).sites(&idx);
                let want = oracle::search_genome(&g, q, n);
                assert_eq!(got, want, "seed {seed} query {} n {n}", kmer::to_string(q));
            }
        }
    }
}

#[test]
fn scan_matches_oracle() {
    check_engine(BuildOptions::default(), 1..=4);
}

#[test]
fn scan_matches_oracle_many_partitions() {
    check_engine(BuildOptions { partition_bits: Some(8), prefix_len: Some(9), ..Default::default() }, 5..=6);
}

#[test]
fn index_round_trip_and_spot_checks() {
    let dir = tempfile::tempdir().unwrap();
    let s = synth::generate(&SynthSpec { seed: 11, ..Default::default() });
    let fa = dir.path().join("g.fa");
    s.genome.write_fasta(&mut std::fs::File::create(&fa).unwrap()).unwrap();
    let m =
        build(&fa, &dir.path().join("idx"), &BuildOptions { partition_bits: Some(4), ..Default::default() }).unwrap();
    let idx = Index::open(&dir.path().join("idx")).unwrap();
    sm_index::verify::verify_structure(&idx).unwrap();
    sm_index::verify::verify_checksums(&idx).unwrap();
    let checked = sm_index::verify::spot_check(&idx, &fa, 5000, 3).unwrap();
    assert_eq!(checked, 5000);
    // Every distinct canonical window is in copy A exactly once.
    let mut distinct: Vec<u64> = kmer::windows(&s.genome.seq).map(|(_, k)| kmer::canonical(k)).collect();
    // windows() over the concatenation may span contigs; restrict to per-contig windows.
    distinct.clear();
    for (i, _) in s.genome.contigs.iter().enumerate() {
        distinct.extend(kmer::windows(s.genome.contig_seq(i)).map(|(_, k)| kmer::canonical(k)));
    }
    let sites = distinct.len() as u64;
    distinct.sort_unstable();
    distinct.dedup();
    assert_eq!(m.distinct_kmers, distinct.len() as u64);
    assert_eq!(m.sites, sites);
    for (leaf, &k) in distinct.iter().enumerate() {
        assert_eq!(idx.a.lookup(k), Some(leaf as u32));
        assert_eq!(idx.a.kmer_at(leaf), k);
    }
    // Stored sequence round-trips, including N and soft-masking.
    let g = &s.genome;
    let back = idx.seq.fetch(0, g.seq.len() as u32);
    let norm: Vec<u8> = g
        .seq
        .iter()
        .map(|&c| {
            if kmer::BASE_CODE[c as usize] == kmer::INVALID {
                if c.is_ascii_lowercase() { b'n' } else { b'N' }
            } else {
                c
            }
        })
        .collect();
    let back_norm: Vec<u8> = back.iter().zip(&norm).map(|(&b, &n)| if n == b'n' { b'n' } else { b }).collect();
    assert_eq!(back_norm, norm);
}

/// The P = 15 layout with `u32` suffix words (human scale). Writes a 4.3 GB prefix table.
#[test]
#[ignore = "writes a 4.3 GB prefix table; run with --ignored"]
fn scan_matches_oracle_u32_words() {
    check_engine(BuildOptions { prefix_len: Some(15), ..Default::default() }, 1..=1);
}

//! Property tests (design §15.3) for the encoding primitives and the oracle.

use proptest::prelude::*;
use sm_core::kmer::*;
use sm_core::oracle;
use sm_core::synth::{self, SynthSpec};

fn kmer() -> impl Strategy<Value = u64> {
    any::<u64>().prop_map(|x| x & KMER_MASK)
}

/// Loop-based reverse complement as an independent reference.
fn revcomp_slow(k: u64) -> u64 {
    (0..K).fold(0, |acc, i| (acc << 2) | (3 - base_at(k, K - 1 - i)) as u64)
}

proptest! {
    #[test]
    fn revcomp_involution(k in kmer()) {
        prop_assert_eq!(revcomp(revcomp(k)), k);
        prop_assert_eq!(revcomp(k), revcomp_slow(k));
        prop_assert!(revcomp(k) <= KMER_MASK);
        prop_assert_ne!(revcomp(k), k, "odd k: never self-complementary");
    }

    #[test]
    fn canonical_idempotent(k in kmer()) {
        let c = canonical(k);
        prop_assert_eq!(canonical(c), c);
        prop_assert_eq!(canonical(revcomp(k)), c);
    }

    #[test]
    fn hamming_properties(a in kmer(), b in kmer()) {
        prop_assert_eq!(hamming(a, b), hamming(b, a));
        prop_assert_eq!(hamming(a, a), 0);
        prop_assert_eq!(hamming(revcomp(a), revcomp(b)), hamming(a, b));
        let slow = (0..K).filter(|&i| base_at(a, i) != base_at(b, i)).count() as u32;
        prop_assert_eq!(hamming(a, b), slow);
        prop_assert_eq!(mismatch_positions(a, b).count() as u32, slow);
    }

    #[test]
    fn rotate_bijection(k in kmer()) {
        let r = rotate(k);
        prop_assert!(r <= KMER_MASK);
        prop_assert_eq!(unrotate(r), k);
        // Rotation preserves distances.
        prop_assert_eq!(hamming(rotate(k), rotate(revcomp(k))), hamming(k, revcomp(k)));
    }

    #[test]
    fn pigeonhole_split(a in kmer(), b in kmer()) {
        let lead = hamming(a & LEAD15, b & LEAD15);
        let back = hamming(a & BACK16, b & BACK16);
        prop_assert_eq!(lead + back, hamming(a, b));
    }
}

/// Oracle invariants on synthetic genomes: every reported site really reads the variant
/// (or its reverse complement) and query/revcomp(query) results mirror each other.
#[test]
fn oracle_self_consistent() {
    for seed in 1..=4 {
        let s = synth::generate(&SynthSpec { seed, ..Default::default() });
        let g = &s.genome;
        let src = &s.family_sources[0];
        let q = parse(std::str::from_utf8(&src[10..10 + K as usize]).unwrap().to_ascii_uppercase().as_str()).unwrap();
        let hits = oracle::search_genome(g, q, 4);
        assert!(hits.len() >= 2, "family query should hit several copies (seed {seed}): {}", hits.len());
        for h in &hits {
            let window = &g.seq[h.pos as usize..h.pos as usize + K as usize];
            let fwd = windows(window).next().unwrap().1;
            let expect = match h.strand {
                oracle::Strand::Forward => h.variant,
                oracle::Strand::Reverse => revcomp(h.variant),
            };
            assert_eq!(fwd, expect);
            assert_eq!(hamming(q, h.variant), h.mismatches as u32);
        }
        let mirrored = oracle::search_genome(g, revcomp(q), 4);
        let mut flipped: Vec<_> = hits
            .iter()
            .map(|h| oracle::Hit {
                variant: revcomp(h.variant),
                strand: if h.strand == oracle::Strand::Forward {
                    oracle::Strand::Reverse
                } else {
                    oracle::Strand::Forward
                },
                ..*h
            })
            .collect();
        flipped.sort_unstable();
        assert_eq!(mirrored, flipped);
    }
}

#[test]
fn synth_has_inverted_hits_and_boundaries() {
    let s = synth::generate(&SynthSpec::default());
    let g = &s.genome;
    assert!(g.seq.iter().any(|c| c.is_ascii_lowercase()));
    assert!(g.seq.contains(&b'N'));
    let mut minus = 0;
    for src in &s.family_sources {
        let q = parse(std::str::from_utf8(&src[..K as usize]).unwrap()).unwrap();
        minus += oracle::search_genome(g, q, 6).iter().filter(|h| h.strand == oracle::Strand::Reverse).count();
    }
    assert!(minus > 0, "inverted repeats should produce '-' hits");
}

//! Seeded synthetic genomes that exercise the awkward cases (design §15.2):
//! several contigs (including ones shorter than k), N runs, soft-masked lowercase,
//! direct and inverted repeat families with point mutations, tandem arrays and
//! low-complexity runs.

use crate::fasta::Genome;
use crate::kmer::{BASE_CHAR, K};

/// SplitMix64: tiny, fast, deterministic.
#[derive(Clone, Debug)]
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Self {
        Self(seed)
    }

    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform in `lo..hi`.
    pub fn range(&mut self, lo: usize, hi: usize) -> usize {
        debug_assert!(lo < hi);
        lo + (self.next_u64() % (hi - lo) as u64) as usize
    }

    /// True with probability `p`.
    pub fn chance(&mut self, p: f64) -> bool {
        ((self.next_u64() >> 11) as f64) / ((1u64 << 53) as f64) < p
    }

    pub fn base(&mut self) -> u8 {
        BASE_CHAR[(self.next_u64() & 3) as usize]
    }

    /// A random 31-mer.
    pub fn kmer(&mut self) -> u64 {
        self.next_u64() & crate::kmer::KMER_MASK
    }
}

#[derive(Clone, Debug)]
pub struct SynthSpec {
    pub seed: u64,
    /// Lengths of the contigs; include some shorter than k.
    pub contig_lens: Vec<usize>,
    /// Repeat families, each inserted this many times (direct or inverted).
    pub families: usize,
    pub copies_per_family: usize,
    /// Per-base substitution rate when copying a family member.
    pub mutation_rate: f64,
    /// Number of N runs, low-complexity runs and tandem arrays to sprinkle in.
    pub n_runs: usize,
    pub low_complexity_runs: usize,
    pub tandem_arrays: usize,
    /// Fraction of the genome soft-masked (lowercased) in blocks.
    pub lowercase_fraction: f64,
}

impl Default for SynthSpec {
    fn default() -> Self {
        Self {
            seed: 1,
            contig_lens: vec![20_000, 7_500, 30, 12, 0, 3_000, 45_000],
            families: 6,
            copies_per_family: 8,
            mutation_rate: 0.03,
            n_runs: 10,
            low_complexity_runs: 10,
            tandem_arrays: 3,
            lowercase_fraction: 0.15,
        }
    }
}

/// A generated genome plus sequences worth querying.
pub struct Synth {
    pub genome: Genome,
    /// Family source sequences (unmutated); k-mers from these make good repeat queries.
    pub family_sources: Vec<Vec<u8>>,
}

pub fn revcomp_seq(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&c| match c {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'a' => b't',
            b'c' => b'g',
            b'g' => b'c',
            b't' => b'a',
            other => other,
        })
        .collect()
}

fn mutate(rng: &mut Rng, s: &[u8], rate: f64) -> Vec<u8> {
    s.iter()
        .map(|&c| {
            if rng.chance(rate) {
                loop {
                    let b = rng.base();
                    if b != c.to_ascii_uppercase() {
                        break b;
                    }
                }
            } else {
                c
            }
        })
        .collect()
}

/// Overwrite `dst[at..]` with `src`, clipped to `dst`.
fn paste(dst: &mut [u8], at: usize, src: &[u8]) {
    let end = (at + src.len()).min(dst.len());
    dst[at..end].copy_from_slice(&src[..end - at]);
}

pub fn generate(spec: &SynthSpec) -> Synth {
    let mut rng = Rng::new(spec.seed);
    let mut contigs: Vec<Vec<u8>> = spec.contig_lens.iter().map(|&n| (0..n).map(|_| rng.base()).collect()).collect();
    let big: Vec<usize> = (0..contigs.len()).filter(|&i| contigs[i].len() > 4 * K as usize).collect();
    let mut family_sources = Vec::new();

    if !big.is_empty() {
        let lens: Vec<usize> = contigs.iter().map(Vec::len).collect();
        let pick = |rng: &mut Rng, len: usize| -> (usize, usize) {
            let c = big[rng.range(0, big.len())];
            (c, rng.range(0, lens[c].saturating_sub(len).max(1)))
        };

        // Repeat families: direct and inverted copies with point mutations.
        for _ in 0..spec.families {
            let len = rng.range(60, 400);
            let src: Vec<u8> = (0..len).map(|_| rng.base()).collect();
            for _ in 0..spec.copies_per_family {
                let rate = spec.mutation_rate * rng.range(0, 3) as f64;
                let mut copy = mutate(&mut rng, &src, rate);
                if rng.chance(0.4) {
                    copy = revcomp_seq(&copy);
                }
                let (c, at) = pick(&mut rng, copy.len());
                paste(&mut contigs[c], at, &copy);
            }
            family_sources.push(src);
        }

        // Tandem arrays of a short unit, lightly mutated per copy.
        for _ in 0..spec.tandem_arrays {
            let unit: Vec<u8> = (0..rng.range(5, 40)).map(|_| rng.base()).collect();
            let copies = rng.range(10, 60);
            let arr: Vec<u8> = (0..copies).flat_map(|_| mutate(&mut rng, &unit, 0.02)).collect();
            let (c, at) = pick(&mut rng, arr.len());
            paste(&mut contigs[c], at, &arr);
            family_sources.push(arr);
        }

        // Low-complexity: homopolymers and dinucleotide runs.
        for _ in 0..spec.low_complexity_runs {
            let len = rng.range(20, 120);
            let a = rng.base();
            let b = if rng.chance(0.5) { a } else { rng.base() };
            let run: Vec<u8> = (0..len).map(|i| if i % 2 == 0 { a } else { b }).collect();
            let (c, at) = pick(&mut rng, len);
            paste(&mut contigs[c], at, &run);
        }

        // N runs (and the odd IUPAC code).
        for _ in 0..spec.n_runs {
            let len = rng.range(1, 60);
            let fill = if rng.chance(0.8) { b'N' } else { b"RYKMSWn"[rng.range(0, 7)] };
            let (c, at) = pick(&mut rng, len);
            paste(&mut contigs[c], at, &vec![fill; len]);
        }

        // Soft-masking in blocks.
        let total: usize = contigs.iter().map(Vec::len).sum();
        let mut masked = 0.0;
        while masked < spec.lowercase_fraction * total as f64 {
            let len = rng.range(50, 500);
            let (c, at) = pick(&mut rng, len);
            let end = (at + len).min(contigs[c].len());
            contigs[c][at..end].make_ascii_lowercase();
            masked += (end - at) as f64;
        }
    }

    let mut genome = Genome::default();
    for (i, seq) in contigs.iter().enumerate() {
        genome.push(format!("contig{i}"), format!("synthetic seed={}", spec.seed), seq).unwrap();
    }
    Synth { genome, family_sources }
}

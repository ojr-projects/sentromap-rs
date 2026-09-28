//! Hamming-ball sizes and the chance background (design §7.5, §9.7, §20.1).

use crate::kmer::K;

/// Binomial coefficient C(n, k), exact.
pub fn binom(n: u32, k: u32) -> u128 {
    if k > n {
        return 0;
    }
    let k = k.min(n - k);
    let mut r: u128 = 1;
    for i in 0..k {
        r = r * (n - i) as u128 / (i + 1) as u128;
    }
    r
}

/// Number of length-`d` sequences within `r` substitutions of a given one:
/// `Σ_{j ≤ r} C(d, j) · 3^j`.
pub fn ball(d: u32, r: u32) -> u128 {
    (0..=r.min(d)).map(|j| binom(d, j) * 3u128.pow(j)).sum()
}

/// Number of sequences at exactly distance `j` from a given one.
pub fn shell(d: u32, j: u32) -> u128 {
    binom(d, j) * 3u128.pow(j)
}

/// Probability that a random 31-mer lies within `n` of a given one.
pub fn p_within(n: u32) -> f64 {
    ball(K, n) as f64 / 4f64.powi(K as i32)
}

/// Expected chance hits for a random query against `distinct` random canonical k-mers,
/// counting both strands.
pub fn expected_chance_hits(distinct: u64, n: u32) -> f64 {
    2.0 * distinct as f64 * p_within(n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn small_values() {
        assert_eq!(binom(31, 0), 1);
        assert_eq!(binom(31, 1), 31);
        assert_eq!(binom(31, 15), 300_540_195);
        assert_eq!(ball(31, 0), 1);
        assert_eq!(ball(31, 1), 1 + 93);
        assert_eq!(ball(31, 31), 4u128.pow(31));
        assert_eq!(ball(15, 15), 4u128.pow(15));
    }

    #[test]
    fn matches_design_appendix() {
        // §20.1: P(dist ≤ 10) ≈ 6.7e-7, ≤ 15 ≈ 1.3e-3 (rounded in the table).
        let p10 = p_within(10);
        assert!((p10 / 6.7e-7 - 1.0).abs() < 0.05, "{p10}");
        // §9.7: D. mel N = 122,317,344 → ~160 chance hits at n = 10, ~318 k at n = 15.
        let e10 = expected_chance_hits(122_317_344, 10);
        let e15 = expected_chance_hits(122_317_344, 15);
        assert!((140.0..180.0).contains(&e10), "{e10}");
        assert!((290e3..340e3).contains(&e15), "{e15}");
    }
}

//! 2-bit k-mer encoding (design §5.1).
//!
//! A k-mer is a `u64` with the first base in bits 61–60 and the last base in bits 1–0.
//! Bases are A=0, C=1, G=2, T=3, so complement is `x ^ 3`.

use std::fmt;

/// k-mer length. Fixed: odd, so no k-mer is its own reverse complement.
pub const K: u32 = 31;
/// Mask of the 62 bits a k-mer occupies.
pub const KMER_MASK: u64 = (1 << (2 * K)) - 1;
/// Bases 15..30 (the back 16 bases) of a k-mer: its low 32 bits.
pub const BACK16: u64 = (1 << 32) - 1;
/// Number of leading bases in the pigeonhole front half.
pub const LEAD: u32 = 15;
/// Bases 0..14 (the lead 15 bases): bits 61..32.
pub const LEAD15: u64 = KMER_MASK & !BACK16;

const LOW_BITS: u64 = 0x5555_5555_5555_5555;

/// Sentinel returned by [`BASE_CODE`] for anything that is not A/C/G/T (either case).
pub const INVALID: u8 = 4;

/// ASCII to 2-bit code; lowercase is sequence (soft-masked), everything else is [`INVALID`].
pub static BASE_CODE: [u8; 256] = {
    let mut t = [INVALID; 256];
    t[b'A' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'T' as usize] = 3;
    t[b'a' as usize] = 0;
    t[b'c' as usize] = 1;
    t[b'g' as usize] = 2;
    t[b't' as usize] = 3;
    t
};

pub const BASE_CHAR: [u8; 4] = *b"ACGT";

/// Reverse complement of a k-mer, without loops.
#[inline]
pub fn revcomp(k: u64) -> u64 {
    let mut k = !k;
    k = ((k >> 2) & 0x3333_3333_3333_3333) | ((k & 0x3333_3333_3333_3333) << 2);
    k = ((k >> 4) & 0x0F0F_0F0F_0F0F_0F0F) | ((k & 0x0F0F_0F0F_0F0F_0F0F) << 4);
    k.swap_bytes() >> (64 - 2 * K)
}

/// The smaller of a k-mer and its reverse complement.
#[inline]
pub fn canonical(k: u64) -> u64 {
    k.min(revcomp(k))
}

/// One set bit (the low bit of the base's pair) per differing base.
/// Base `i` (0 = first) is bit `2 * (30 - i)`.
#[inline]
pub fn diff_mask(a: u64, b: u64) -> u64 {
    let x = a ^ b;
    (x | (x >> 1)) & LOW_BITS
}

/// Number of differing bases.
#[inline]
pub fn hamming(a: u64, b: u64) -> u32 {
    diff_mask(a, b).count_ones()
}

/// Rotate so the back 16 bases lead: bases 15..30 followed by bases 0..14 (design §6.3).
#[inline]
pub fn rotate(k: u64) -> u64 {
    ((k & BACK16) << 30) | (k >> 32)
}

/// Inverse of [`rotate`].
#[inline]
pub fn unrotate(r: u64) -> u64 {
    ((r & ((1 << 30) - 1)) << 32) | (r >> 30)
}

/// Base at index `i` (0 = first).
#[inline]
pub fn base_at(k: u64, i: u32) -> u8 {
    ((k >> (2 * (K - 1 - i))) & 3) as u8
}

/// Positions (0 = first base) at which `a` and `b` differ, ascending.
pub fn mismatch_positions(a: u64, b: u64) -> impl Iterator<Item = u32> {
    let mut m = diff_mask(a, b);
    std::iter::from_fn(move || {
        if m == 0 {
            return None;
        }
        let bit = 63 - m.leading_zeros();
        m &= !(1 << bit);
        Some(K - 1 - bit / 2)
    })
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum ParseKmerError {
    #[error("k-mer must be exactly {K} bases, got {0}")]
    Length(usize),
    #[error("invalid base {0:?} at position {1}; only A, C, G, T are allowed")]
    Base(char, usize),
}

/// Parse exactly 31 A/C/G/T bases (case-insensitive).
pub fn parse(s: &str) -> Result<u64, ParseKmerError> {
    let b = s.as_bytes();
    if b.len() != K as usize {
        return Err(ParseKmerError::Length(s.chars().count()));
    }
    let mut k = 0u64;
    for (i, &c) in b.iter().enumerate() {
        let code = BASE_CODE[c as usize];
        if code == INVALID {
            return Err(ParseKmerError::Base(c as char, i));
        }
        k = (k << 2) | code as u64;
    }
    Ok(k)
}

/// Uppercase string form.
pub fn to_string(k: u64) -> String {
    Kmer(k).to_string()
}

/// Display wrapper.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct Kmer(pub u64);

impl fmt::Display for Kmer {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut buf = [0u8; K as usize];
        for (i, c) in buf.iter_mut().enumerate() {
            *c = BASE_CHAR[base_at(self.0, i as u32) as usize];
        }
        f.write_str(std::str::from_utf8(&buf).unwrap())
    }
}

/// Every valid k-mer window in a sequence: `(offset of first base, forward k-mer)`.
/// Windows containing a non-ACGT byte are skipped.
pub fn windows(seq: &[u8]) -> impl Iterator<Item = (usize, u64)> + '_ {
    let mut k = 0u64;
    let mut valid = 0u32;
    seq.iter().enumerate().filter_map(move |(i, &c)| {
        let code = BASE_CODE[c as usize];
        if code == INVALID {
            valid = 0;
            return None;
        }
        k = ((k << 2) | code as u64) & KMER_MASK;
        valid += 1;
        if valid >= K { Some((i + 1 - K as usize, k)) } else { None }
    })
}

/// Every valid window with both orientations, rolled incrementally:
/// `(offset of first base, forward k-mer, reverse complement)`.
pub fn windows_both(seq: &[u8]) -> impl Iterator<Item = (usize, u64, u64)> + '_ {
    let mut fwd = 0u64;
    let mut rc = 0u64;
    let mut valid = 0u32;
    seq.iter().enumerate().filter_map(move |(i, &c)| {
        let code = BASE_CODE[c as usize];
        if code == INVALID {
            valid = 0;
            return None;
        }
        fwd = ((fwd << 2) | code as u64) & KMER_MASK;
        rc = (rc >> 2) | (((3 - code) as u64) << (2 * (K - 1)));
        valid += 1;
        if valid >= K { Some((i + 1 - K as usize, fwd, rc)) } else { None }
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_display_round_trip() {
        let s = "TCGTAAAATATGGAGACTTTTACTGGGTATC";
        let k = parse(s).unwrap();
        assert_eq!(to_string(k), s);
        assert_eq!(parse(&s.to_lowercase()).unwrap(), k);
        assert_eq!(parse("ACGT"), Err(ParseKmerError::Length(4)));
        assert!(matches!(parse(&"N".repeat(31)), Err(ParseKmerError::Base('N', 0))));
    }

    #[test]
    fn revcomp_known() {
        let k = parse("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAC").unwrap();
        assert_eq!(to_string(revcomp(k)), "GTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
        assert_eq!(canonical(k), k);
    }

    #[test]
    fn rotate_layout() {
        let k = parse("AAAAAAAAAAAAAAACCCCCCCCCCCCCCCG").unwrap();
        // Bases 15..30 then 0..14.
        assert_eq!(to_string(rotate(k)), "CCCCCCCCCCCCCCCGAAAAAAAAAAAAAAA");
        assert_eq!(unrotate(rotate(k)), k);
    }

    #[test]
    fn mismatch_positions_order() {
        let a = parse("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA").unwrap();
        let b = parse("CAAAAAAAAAGAAAAAAAAAAAAAAAAAAAT").unwrap();
        assert_eq!(mismatch_positions(a, b).collect::<Vec<_>>(), vec![0, 10, 30]);
        assert_eq!(hamming(a, b), 3);
    }

    #[test]
    fn windows_skip_invalid() {
        let mut seq = b"ACGT".repeat(10);
        seq[20] = b'N';
        let w: Vec<_> = windows(&seq).collect();
        assert!(w.is_empty(), "no 31-long run without N");
        let seq = b"acgtACGTacgtACGTacgtACGTacgtACGTa";
        let w: Vec<_> = windows(seq).collect();
        assert_eq!(w.len(), 3);
        assert_eq!(w[0].0, 0);
        assert_eq!(to_string(w[0].1), "ACGTACGTACGTACGTACGTACGTACGTACG");
    }

    #[test]
    fn windows_both_rolls_revcomp() {
        let seq = b"ACGGTTACGATTTAGGCATGCANNNACGTGGGGTACCATGGACTTTAGCAGCATCAGACTTAGAC";
        let a: Vec<_> = windows(seq).collect();
        let b: Vec<_> = windows_both(seq).collect();
        assert_eq!(a.len(), b.len());
        for ((i, f), (j, f2, r)) in a.into_iter().zip(b) {
            assert_eq!((i, f), (j, f2));
            assert_eq!(r, revcomp(f));
        }
    }
}

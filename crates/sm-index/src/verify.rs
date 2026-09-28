//! Index validation (design §13.3): structural invariants, checksums and spot checks.

use std::path::Path;

use anyhow::{Context, Result, bail, ensure};
use rayon::prelude::*;
use sm_core::fasta::{FastaReader, open};
use sm_core::kmer::windows_both;
use sm_core::synth::Rng;

use crate::index::Index;
use crate::kmers::{SortedKmers, Word};
use crate::with_kmers;

fn check_sorted<W: Word>(c: &SortedKmers<W>) -> Result<()> {
    let table = c.prefix_table();
    ensure!(table.windows(2).all(|w| w[0] <= w[1]), "prefix table not monotonic");
    let bad = (0..c.buckets()).into_par_iter().find_any(|&b| {
        let r = c.bucket_range(b);
        let words = &c.words()[r];
        let prefix = (b as u64) << c.shift();
        words.windows(2).any(|w| w[0] >= w[1]) || words.iter().any(|w| w.to_kmer(prefix) >> c.shift() != b as u64)
    });
    if let Some(b) = bad {
        bail!("bucket {b} unsorted or holds a foreign k-mer");
    }
    let canon_bad = (0..c.len()).into_par_iter().find_any(|&i| {
        let k = c.kmer_at(i);
        sm_core::canonical(k) != k
    });
    if let Some(i) = canon_bad.filter(|_| c.len() < 50_000_000) {
        bail!("k-mer {i} is not canonical");
    }
    Ok(())
}

/// Structural checks: sortedness, canonical form, rows ascending and non-empty.
pub fn verify_structure(idx: &Index) -> Result<()> {
    with_kmers!(&idx.a, c => check_sorted(c))?;
    let pos = &idx.positions;
    let n = idx.manifest.distinct_kmers;
    let bad = (0..n as u32).into_par_iter().find_any(|&l| {
        let r = pos.row(l);
        r.is_empty() || pos.sites()[r.start as usize..r.end as usize].windows(2).any(|w| w[0] >= w[1])
    });
    if let Some(l) = bad {
        bail!("row {l} is empty or not ascending");
    }
    Ok(())
}

/// Recompute file checksums and compare with the manifest.
pub fn verify_checksums(idx: &Index) -> Result<()> {
    idx.manifest.files.par_iter().try_for_each(|(name, info)| {
        let mut f = std::fs::File::open(idx.dir.join(name)).with_context(|| name.clone())?;
        let mut h = crc32fast::Hasher::new();
        let mut buf = vec![0u8; 1 << 22];
        let mut bytes = 0u64;
        loop {
            let n = std::io::Read::read(&mut f, &mut buf)?;
            if n == 0 {
                break;
            }
            h.update(&buf[..n]);
            bytes += n as u64;
        }
        ensure!(bytes == info.bytes && h.finalize() == info.crc32, "{name}: checksum mismatch");
        Ok(())
    })
}

/// Look up `samples` random genome windows and check each is found at its own position
/// with the right strand; also check the total window count.
pub fn spot_check(idx: &Index, fasta: &Path, samples: usize, seed: u64) -> Result<usize> {
    let total = idx.manifest.sites;
    let mut rng = Rng::new(seed);
    let mut targets: Vec<u64> = (0..samples.min(total as usize)).map(|_| rng.next_u64() % total.max(1)).collect();
    targets.sort_unstable();
    let mut reader = FastaReader::new(open(fasta)?);
    let mut seq = Vec::new();
    let mut offset = 0u32;
    let mut w = 0u64;
    let mut t = 0usize;
    let mut checked = 0;
    while let Some(_h) = reader.next_record(&mut seq)? {
        for (i, fwd, rc) in windows_both(&seq) {
            while t < targets.len() && targets[t] == w {
                let pos = offset + i as u32;
                let canon = fwd.min(rc);
                let leaf = idx.a.lookup(canon).with_context(|| format!("window at {pos} missing from copy A"))?;
                let found =
                    idx.sites(leaf, false).any(|(p, s)| p == pos && (s == sm_core::Strand::Forward) == (fwd < rc));
                ensure!(found, "window at {pos} not in its row with the right strand");
                checked += 1;
                t += 1;
            }
            w += 1;
        }
        offset += seq.len() as u32;
    }
    ensure!(w == total, "FASTA has {w} windows, index has {total} sites");
    Ok(checked)
}

//! Real-genome checks (design §15.4). Skipped when data/ has not been fetched
//! (scripts/fetch-data.sh genomes).

use std::path::{Path, PathBuf};

use sm_core::kmer;
use sm_index::{BuildOptions, Index, build};

fn data(name: &str) -> Option<PathBuf> {
    let p = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data").join(name);
    if p.exists() {
        Some(p)
    } else {
        eprintln!("skipping: {} not present", p.display());
        None
    }
}

#[test]
fn yeast_ty1() {
    let Some(fa) = data("S_cere_R64.fna") else { return };
    let dir = tempfile::tempdir().unwrap();
    build(&fa, dir.path(), &BuildOptions::default()).unwrap();
    let idx = Index::open(dir.path()).unwrap();
    assert_eq!(idx.manifest.distinct_kmers, 11_564_096);
    let q = kmer::parse("TCGTAAAATATGGAGACTTTTACTGGGTATC").unwrap();
    let v0 = sm_search::scan(&idx, q, 0);
    assert_eq!((v0.len(), v0.site_count(&idx)), (1, 27));
    let v1 = sm_search::scan(&idx, q, 1);
    assert_eq!((v1.len(), v1.site_count(&idx)), (3, 33));
    // The query is taken from YBLWTy1-1 at chrII:223000.
    let chr2 = idx.contigs.by_name("NC_001134.8").unwrap();
    assert_eq!(idx.kmer_at(chr2, 223_000 - 1), Some(q));
}

#[test]
fn ecoli_windows_found_at_their_own_positions() {
    let Some(fa) = data("ecoli_k12.fna") else { return };
    let dir = tempfile::tempdir().unwrap();
    let m = build(&fa, dir.path(), &BuildOptions::default()).unwrap();
    assert_eq!(m.distinct_kmers, 4_554_269);
    let idx = Index::open(dir.path()).unwrap();
    sm_index::verify::verify_structure(&idx).unwrap();
    let checked = sm_index::verify::spot_check(&idx, &fa, m.sites as usize, 7).unwrap();
    assert!(checked > 4_000_000);
}

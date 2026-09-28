use sm_agg::*;
use sm_core::kmer;
use sm_core::synth::{self, SynthSpec};
use sm_index::{BuildOptions, Index, build};

fn fixture() -> (tempfile::TempDir, Index, u64) {
    let dir = tempfile::tempdir().unwrap();
    let s = synth::generate(&SynthSpec { seed: 3, copies_per_family: 20, ..Default::default() });
    let fa = dir.path().join("g.fa");
    s.genome.write_fasta(&mut std::fs::File::create(&fa).unwrap()).unwrap();
    build(&fa, &dir.path().join("idx"), &BuildOptions::default()).unwrap();
    let idx = Index::open(&dir.path().join("idx")).unwrap();
    let q = kmer::parse(std::str::from_utf8(&s.family_sources[0][5..36]).unwrap()).unwrap();
    (dir, idx, kmer::canonical(q))
}

#[test]
fn result_set_and_views() {
    let (_d, idx, q) = fixture();
    let v = sm_search::scan(&idx, q, 8);
    let total_sites = v.site_count(&idx);
    let rs = ResultSet::new(&idx, v);
    assert!(rs.variants.mismatches.windows(2).all(|w| w[0] <= w[1]));
    assert_eq!(rs.sites.pos.len() as u64, total_sites);
    assert_eq!(*rs.sites_within.last().unwrap(), total_sites);
    assert!(rs.sites.pos.windows(2).all(|w| w[0] <= w[1]));
    let len = idx.manifest.genome_len as u32;
    for mode in [OrderMode::PerN, OrderMode::Frozen, OrderMode::Mst] {
        for n in [2u8, 5, 8] {
            let ov = rs.order(mode, n);
            assert_eq!(ov.rows(), rs.count(n));
            let mut o = ov.order.clone();
            o.sort_unstable();
            assert_eq!(o, (0..rs.count(n) as u32).collect::<Vec<_>>());
            // A raster covering everything draws every site within n exactly once.
            for width in [997, 5000, 60_000] {
                let view = View { n, start: 0, end: len, row0: 0, row1: ov.rows() as u32, width, height: 311 };
                let r = raster(&rs, &ov, &view);
                assert_eq!(r.sites_drawn, rs.sites_within[n as usize], "{mode:?} n={n} width={width}");
                let cells: u64 = r.count.iter().map(|&c| c as u64).sum();
                assert!(cells >= r.sites_drawn);
            }
            let h = histogram(&rs, n, 0, len, 50);
            assert_eq!(h.iter().map(|&x| x as u64).sum::<u64>(), rs.sites_within[n as usize]);
            // Clades nest within the row range.
            for c in &ov.clades {
                assert!(c.start < c.end && c.end as usize <= ov.rows());
            }
        }
    }
    // Frozen order at n is a subsequence of the max-n frozen order.
    let full = rs.order(OrderMode::Frozen, 8);
    let small = rs.order(OrderMode::Frozen, 3);
    let mut it = full.order.iter();
    assert!(small.order.iter().all(|x| it.any(|y| y == x)));
}

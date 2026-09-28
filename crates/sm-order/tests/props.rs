//! Ordering properties (design §15.3).

use proptest::prelude::*;
use sm_core::kmer::{KMER_MASK, hamming};
use sm_core::synth::Rng;
use sm_order::chartree::{characters, lcp};
use sm_order::*;

/// A synthetic family: random mutations of the query arranged in a real tree, plus noise.
fn family(seed: u64, size: usize, max_n: u32) -> (u64, Vec<u64>, Vec<u8>) {
    let mut rng = Rng::new(seed);
    let q = rng.kmer();
    let mut pool = vec![q];
    let mut out = std::collections::BTreeSet::new();
    while out.len() < size {
        let parent = pool[rng.range(0, pool.len())];
        let pos = rng.range(0, 31) as u32;
        let shift = 2 * (30 - pos);
        let child = (parent & !(3 << shift)) | ((rng.next_u64() & 3) << shift);
        let child = child & KMER_MASK;
        if child != q && hamming(q, child) <= max_n {
            pool.push(child);
            out.insert(child);
        }
    }
    let kmers: Vec<u64> = out.into_iter().collect();
    let mism = kmers.iter().map(|&k| hamming(q, k) as u8).collect();
    (q, kmers, mism)
}

fn is_permutation(order: &[u32], mism: &[u8], n: u8) -> bool {
    let mut o = order.to_vec();
    o.sort_unstable();
    let want: Vec<u32> = (0..mism.len() as u32).filter(|&i| mism[i as usize] <= n).collect();
    o == want
}

/// Tree length recounted from the tree: distinct non-empty key prefixes = trie nodes.
fn recount_tree_length(o: &Ordered, mism: &[u8]) -> u64 {
    let mut nodes = std::collections::HashSet::new();
    for (i, &k) in o.keys.iter().enumerate() {
        let len = mism[o.order[i] as usize] as u32;
        for d in 1..=len {
            let mask = if d * 7 >= 128 { !0u128 } else { !(!0u128 >> (d * 7)) };
            nodes.insert((d, k & mask));
        }
    }
    nodes.len() as u64
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(40))]
    #[test]
    fn chartree_invariants(seed in any::<u64>(), size in 1usize..400, max_n in 2u32..10) {
        let (q, kmers, mism) = family(seed, size, max_n);
        let top = max_n as u8;
        let (_r, full) = order_per_n(q, &kmers, &mism, top);
        prop_assert!(is_permutation(&full.order, &mism, top));
        prop_assert!(full.keys.windows(2).all(|w| w[0] < w[1]), "keys unique and sorted");
        prop_assert_eq!(full.tree_length, recount_tree_length(&full, &mism));
        prop_assert!(full.tree_length <= star_length(&mism, top));
        for n in 0..top {
            let frozen = filter_frozen(&full, &mism, n);
            prop_assert!(is_permutation(&frozen.order, &mism, n));
            // Frozen order at n is a subsequence of the order at n + 1 (and at max n).
            let next = filter_frozen(&full, &mism, n + 1);
            let mut it = next.order.iter();
            prop_assert!(frozen.order.iter().all(|x| it.any(|y| y == x)));
            let (_r, per_n) = order_per_n(q, &kmers, &mism, n);
            prop_assert!(is_permutation(&per_n.order, &mism, n));
            prop_assert_eq!(per_n.tree_length, recount_tree_length(&per_n, &mism));
        }
        // Clades: nested, contiguous, and members really share the clade's prefix.
        for c in full.clades() {
            prop_assert!(c.start < c.end || full.is_empty());
            for i in c.start + 1..c.end {
                prop_assert!(lcp(full.keys[c.start as usize], full.keys[i as usize]) >= c.depth);
            }
        }
    }

    #[test]
    fn mst_is_spanning_and_minimal(seed in any::<u64>(), size in 1usize..120) {
        let (q, kmers, mism) = family(seed, size, 8);
        let t = mst(q, &kmers, &mism, 8);
        prop_assert!(is_permutation(&t.order, &mism, 8));
        // Parents precede children in pre-order.
        let pos: std::collections::HashMap<u32, usize> = t.order.iter().enumerate().map(|(p, &i)| (i, p)).collect();
        for &i in &t.order {
            let p = t.parent[i as usize];
            if p != mst::ROOT {
                prop_assert!(pos[&p] < pos[&i]);
            }
        }
        // Kruskal on the full graph gives the same total weight.
        let m = kmers.len();
        let mut edges = Vec::new();
        for i in 0..m {
            edges.push((hamming(q, kmers[i]), m, i));
            for j in i + 1..m {
                edges.push((hamming(kmers[i], kmers[j]), i, j));
            }
        }
        edges.sort_unstable();
        let mut uf: Vec<usize> = (0..=m).collect();
        fn find(uf: &mut [usize], x: usize) -> usize { if uf[x] != x { let r = find(uf, uf[x]); uf[x] = r; } uf[x] }
        let mut total = 0u64;
        for (w, a, b) in edges {
            let (ra, rb) = (find(&mut uf, a), find(&mut uf, b));
            if ra != rb { uf[ra] = rb; total += w as u64; }
        }
        prop_assert_eq!(t.tree_length, total);
    }
}

#[test]
fn characters_decode() {
    let q = sm_core::kmer::parse("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA").unwrap();
    let v = sm_core::kmer::parse("CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAT").unwrap();
    let c: Vec<u8> = characters(q, v).collect();
    assert_eq!(c, vec![1, 30 * 4 + 3]);
}

#[test]
fn star_query_set() {
    // Variants sharing nothing: the tree is the star.
    let q = 0u64;
    let kmers: Vec<u64> = (0..31).map(|i| 1u64 << (2 * i)).collect();
    let mism = vec![1u8; 31];
    let (_r, o) = order_per_n(q, &kmers, &mism, 1);
    assert_eq!(o.tree_length, 31);
    assert_eq!(o.clades().len(), 1);
}

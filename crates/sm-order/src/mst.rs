//! Exact minimum spanning tree for small result sets (design §9.3): Prim's algorithm over
//! the variants plus the query, O(M²) time with no stored matrix, rooted at the query and
//! ordered by pre-order with children visited nearest-first.

use sm_core::hamming;

/// Above this many variants the MST is not attempted (1.2 s single-threaded in the prototype).
pub const MST_LIMIT: usize = 25_000;

pub const ROOT: u32 = u32::MAX;

#[derive(Clone, Debug, Default)]
pub struct Mst {
    /// Indices into the input variants, in pre-order.
    pub order: Vec<u32>,
    /// Parent of each input variant (by input index), `ROOT` for the query.
    pub parent: Vec<u32>,
    pub tree_length: u64,
}

/// MST over the variants with `mismatches ≤ n` (others get parent `ROOT` and are omitted).
pub fn mst(query: u64, kmers: &[u64], mismatches: &[u8], n: u8) -> Mst {
    let members: Vec<u32> = (0..kmers.len() as u32).filter(|&i| mismatches[i as usize] <= n).collect();
    let m = members.len();
    let vk: Vec<u64> = members.iter().map(|&i| kmers[i as usize]).collect();
    // Prim from the root (the query). One fused pass per step updates distances to the node
    // just added and finds the next nearest; sequential, since per-step work is tiny and
    // waking a thread pool 25k times costs more than the work.
    let mut dist: Vec<u32> = vk.iter().map(|&k| hamming(query, k)).collect();
    let mut par: Vec<u32> = vec![ROOT; m];
    let mut done = vec![false; m];
    let mut tree_length = 0u64;
    let mut next = (0..m).min_by_key(|&i| (dist[i], i));
    while let Some(j) = next {
        done[j] = true;
        tree_length += dist[j] as u64;
        let kj = vk[j];
        let mut best = (u32::MAX, usize::MAX);
        for i in 0..m {
            if done[i] {
                continue;
            }
            let h = hamming(kj, vk[i]);
            if h < dist[i] {
                dist[i] = h;
                par[i] = j as u32;
            }
            if (dist[i], i) < best {
                best = (dist[i], i);
            }
        }
        next = (best.1 != usize::MAX).then_some(best.1);
    }
    // Edge weights and children lists.
    let weight =
        |i: usize| -> u32 { if par[i] == ROOT { hamming(query, vk[i]) } else { hamming(vk[par[i] as usize], vk[i]) } };
    let mut children: Vec<Vec<u32>> = vec![Vec::new(); m + 1]; // index m = root
    for (i, &p) in par.iter().enumerate() {
        let p = if p == ROOT { m } else { p as usize };
        children[p].push(i as u32);
    }
    for c in &mut children {
        c.sort_by_key(|&i| (weight(i as usize), i));
    }
    let mut order = Vec::with_capacity(m);
    let mut stack: Vec<u32> = children[m].iter().rev().copied().collect();
    while let Some(i) = stack.pop() {
        order.push(members[i as usize]);
        stack.extend(children[i as usize].iter().rev());
    }
    let mut parent = vec![ROOT; kmers.len()];
    for (i, &p) in par.iter().enumerate() {
        if p != ROOT {
            parent[members[i] as usize] = members[p as usize];
        }
    }
    Mst { order, parent, tree_length }
}

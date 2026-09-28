//! Row × genome rasters at screen resolution, and per-bin mismatch histograms (design §10).
//!
//! The server renders the variant-by-position matrix directly at the client's pixel size, so
//! the browser only ever draws a screen-sized image, however many rows the result has.

use rayon::prelude::*;

use crate::order::OrderView;
use crate::result::{REVERSE, ResultSet};

/// A view: global `[start, end)` across `width` columns, rows `[row0, row1)` across `height`.
#[derive(Clone, Copy, Debug)]
pub struct View {
    pub n: u8,
    pub start: u32,
    pub end: u32,
    pub row0: u32,
    pub row1: u32,
    pub width: u32,
    pub height: u32,
}

/// One cell per pixel, row-major. Layout per cell: count (u16, saturating), the lowest mismatch
/// count among its sites (u8, 255 = empty), and strand bits (u8: 1 = '+', 2 = '-').
pub struct Raster {
    pub width: u32,
    pub height: u32,
    pub count: Vec<u16>,
    pub min_mm: Vec<u8>,
    pub strands: Vec<u8>,
    pub sites_drawn: u64,
}

impl Raster {
    /// Interleaved little-endian bytes: `[count_lo, count_hi, min_mm, strands]` per cell.
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(self.count.len() * 4);
        for i in 0..self.count.len() {
            out.extend_from_slice(&self.count[i].to_le_bytes());
            out.push(self.min_mm[i]);
            out.push(self.strands[i]);
        }
        out
    }
}

const STRIPES: u32 = 32;

/// Columns `[c0, c1)` rendered by one task.
struct Stripe {
    c0: u32,
    c1: u32,
    count: Vec<u16>,
    min_mm: Vec<u8>,
    strands: Vec<u8>,
    drawn: u64,
}

pub fn raster(rs: &ResultSet, ov: &OrderView, v: &View) -> Raster {
    let (w, h) = (v.width.max(1), v.height.max(1));
    let span = (v.end - v.start).max(1) as f64;
    let rows = (v.row1.saturating_sub(v.row0)).max(1) as f64;
    let m = rs.count(v.n) as u32;
    let k = sm_core::K;
    // Stripes of columns rendered independently (sites are position-sorted, so each stripe
    // reads a contiguous range), then stitched.
    let stripes = STRIPES.min(w);
    let parts: Vec<Stripe> = (0..stripes)
        .into_par_iter()
        .map(|si| {
            let c0 = si * w / stripes;
            let c1 = (si + 1) * w / stripes;
            let sw = (c1 - c0) as usize;
            let (mut cnt, mut mm, mut st) =
                (vec![0u16; sw * h as usize], vec![255u8; sw * h as usize], vec![0u8; sw * h as usize]);
            let g0 = v.start + ((c0 as f64) * span / w as f64) as u32;
            let g1 = v.start + ((c1 as f64) * span / w as f64).ceil() as u32;
            let mut drawn = 0u64;
            for i in rs.sites.range(g0, g1.min(v.end)) {
                let var = rs.sites.var[i];
                let vi = var & !REVERSE;
                if vi >= m {
                    continue;
                }
                let row = ov.row_of[vi as usize];
                if row < v.row0 || row >= v.row1 {
                    continue;
                }
                let pos = rs.sites.pos[i];
                let first = ((((pos.max(v.start) - v.start) as f64) * w as f64 / span) as u32).min(w - 1);
                let last = ((((pos + k).min(v.end) - v.start) as f64) * w as f64 / span).ceil() as u32;
                let (x0, x1) = (first.max(c0), last.max(first + 1).min(c1));
                if x0 >= x1 {
                    continue;
                }
                // Count each site once: in the stripe holding its first column.
                if (c0..c1).contains(&first) {
                    drawn += 1;
                }
                let y0 = ((row - v.row0) as f64 * h as f64 / rows) as u32;
                let y1 = ((((row - v.row0 + 1) as f64) * h as f64 / rows).ceil() as u32).clamp(y0 + 1, h);
                let mis = rs.variants.mismatches[vi as usize];
                let sbit = if var & REVERSE != 0 { 2 } else { 1 };
                for y in y0..y1 {
                    for x in x0..x1 {
                        let c = y as usize * sw + (x - c0) as usize;
                        cnt[c] = cnt[c].saturating_add(1);
                        mm[c] = mm[c].min(mis);
                        st[c] |= sbit;
                    }
                }
            }
            Stripe { c0, c1, count: cnt, min_mm: mm, strands: st, drawn }
        })
        .collect();
    let n = (w * h) as usize;
    let mut out =
        Raster { width: w, height: h, count: vec![0; n], min_mm: vec![255; n], strands: vec![0; n], sites_drawn: 0 };
    for p in parts {
        let sw = (p.c1 - p.c0) as usize;
        for y in 0..h as usize {
            let dst = y * w as usize + p.c0 as usize;
            out.count[dst..dst + sw].copy_from_slice(&p.count[y * sw..(y + 1) * sw]);
            out.min_mm[dst..dst + sw].copy_from_slice(&p.min_mm[y * sw..(y + 1) * sw]);
            out.strands[dst..dst + sw].copy_from_slice(&p.strands[y * sw..(y + 1) * sw]);
        }
        out.sites_drawn += p.drawn;
    }
    out
}

/// Sites per column by mismatch count, over all rows within n: `width × (n + 1)` counts,
/// column-major (`hist[col * (n + 1) + m]`).
pub fn histogram(rs: &ResultSet, n: u8, start: u32, end: u32, width: u32) -> Vec<u32> {
    let w = width.max(1) as usize;
    let stride = n as usize + 1;
    let span = (end - start).max(1) as f64;
    let m = rs.count(n) as u32;
    let mut out = vec![0u32; w * stride];
    for i in rs.sites.range(start, end) {
        let vi = rs.sites.var[i] & !REVERSE;
        if vi >= m {
            continue;
        }
        let pos = rs.sites.pos[i].max(start);
        let x = ((((pos - start) as f64) * w as f64 / span) as usize).min(w - 1);
        out[x * stride + rs.variants.mismatches[vi as usize] as usize] += 1;
    }
    out
}

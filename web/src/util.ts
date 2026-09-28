// Coordinates, formatting and colour palettes.

import type { Contig } from "./api";

export const K = 31;

/** Genome coordinate helpers over the global (concatenated) coordinate space. */
export class Coords {
  readonly contigs: Contig[];
  readonly length: number;
  constructor(contigs: Contig[]) {
    this.contigs = contigs;
    const last = contigs[contigs.length - 1];
    this.length = last ? last.offset + last.length : 0;
  }

  /** Index of the contig containing global position `pos`. */
  locate(pos: number): number {
    let lo = 0;
    let hi = this.contigs.length - 1;
    while (lo < hi) {
      const mid = (lo + hi + 1) >> 1;
      if (this.contigs[mid].offset <= pos) lo = mid;
      else hi = mid - 1;
    }
    return lo;
  }

  /** "2L:1,234,567" (1-based). */
  label(pos: number): string {
    const c = this.contigs[this.locate(pos)];
    return `${c.display}:${fmtInt(Math.floor(pos - c.offset) + 1)}`;
  }

  /** Parse "2L", "2L:100-200", "2L:1,000,000" or "all" into a global [start, end). */
  parseRegion(text: string): [number, number] | null {
    const t = text.trim().replace(/,/g, "");
    if (t === "" || t.toLowerCase() === "all") return [0, this.length];
    const m = /^([^:\s]+)(?::(\d+)(?:-(\d+))?)?$/.exec(t);
    if (!m) return null;
    const c = this.contigs.find((c) => c.display === m[1] || c.name === m[1]);
    if (!c) return null;
    if (!m[2]) return [c.offset, c.offset + c.length];
    const s = Math.max(1, Number(m[2]));
    const e = m[3] ? Number(m[3]) : s + 1000;
    return [c.offset + s - 1, c.offset + Math.min(c.length, Math.max(e, s + 1))];
  }

  /** "2L:1,000,000-1,100,000" for a global range (or "all"). */
  regionLabel(start: number, end: number): string {
    if (start <= 0 && end >= this.length) return "all";
    const ci = this.locate(start);
    const cj = this.locate(Math.max(start, end - 1));
    const c = this.contigs[ci];
    if (ci === cj) return `${c.display}:${fmtInt(start - c.offset + 1)}-${fmtInt(end - c.offset)}`;
    return `${this.label(start)} – ${this.label(end - 1)}`;
  }
}

export function fmtInt(n: number): string {
  return Math.round(n).toLocaleString("en-US");
}

export function fmtCount(n: number): string {
  if (n < 1e4) return fmtInt(n);
  if (n < 1e6) return (n / 1e3).toFixed(n < 1e5 ? 1 : 0) + "k";
  if (n < 1e9) return (n / 1e6).toFixed(n < 1e7 ? 2 : 1) + "M";
  return (n / 1e9).toFixed(2) + "G";
}

export function fmtExpected(x: number): string {
  if (x >= 100) return fmtCount(x);
  if (x >= 0.1) return x.toFixed(1);
  if (x === 0) return "0";
  return x.toExponential(0);
}

export function fmtBp(n: number): string {
  if (n >= 1e6) return (n / 1e6).toFixed(n >= 1e7 ? 1 : 2) + " Mb";
  if (n >= 1e3) return (n / 1e3).toFixed(n >= 1e4 ? 1 : 2) + " kb";
  return Math.round(n) + " bp";
}

export function fmtMs(ms: number): string {
  if (ms < 1) return ms.toFixed(2) + " ms";
  if (ms < 1000) return ms.toFixed(ms < 10 ? 1 : 0) + " ms";
  return (ms / 1000).toFixed(2) + " s";
}

/** A "nice" tick step for a span covering roughly `target` ticks. */
export function niceStep(span: number, target: number): number {
  const raw = span / Math.max(1, target);
  const p = Math.pow(10, Math.floor(Math.log10(raw)));
  for (const m of [1, 2, 5, 10]) if (m * p >= raw) return m * p;
  return 10 * p;
}

// --- colour -----------------------------------------------------------------------------

export type RGB = [number, number, number];

function hex(h: string): RGB {
  const v = parseInt(h.slice(1), 16);
  return [(v >> 16) & 255, (v >> 8) & 255, v & 255];
}

function ramp(stops: string[]): (t: number) => RGB {
  const cs = stops.map(hex);
  return (t: number) => {
    const x = Math.min(1, Math.max(0, t)) * (cs.length - 1);
    const i = Math.min(cs.length - 2, Math.floor(x));
    const f = x - i;
    const a = cs[i];
    const b = cs[i + 1];
    return [a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f, a[2] + (b[2] - a[2]) * f];
  };
}

/** Mismatch count → colour: exact matches bright and warm, distant variants cool. */
export const mismatchRamp = ramp(["#fff6c8", "#ffd23f", "#ff8c42", "#ef476f", "#b55cd6", "#6f7dff", "#3fa7c9"]);

export function mismatchColor(m: number, maxN: number): RGB {
  return mismatchRamp(maxN <= 0 ? 0 : m / maxN);
}

export function css([r, g, b]: RGB, a = 1): string {
  return `rgba(${r | 0},${g | 0},${b | 0},${a})`;
}

/** Signal tracks: viridis-like. */
export const signalRamp = ramp(["#1b1733", "#2d3f8e", "#1f8a8a", "#43bf71", "#c8e03a", "#fde725"]);

export const STRAND_PLUS: RGB = hex("#4cc3ff");
export const STRAND_MINUS: RGB = hex("#ff9a4c");
export const STRAND_BOTH: RGB = hex("#f2f2f2");

/** Distinct colours for feature kinds / repeat classes. */
const KIND_COLORS = ["#7aa2f7", "#e0af68", "#9ece6a", "#f7768e", "#bb9af7", "#2ac3de", "#ff9e64", "#73daca", "#c0caf5", "#db4b4b"];
export function kindColor(i: number): string {
  return KIND_COLORS[i % KIND_COLORS.length];
}

/** Stable colour per repeat class ("LTR/Gypsy" → class "LTR"). */
const CLASS_COLORS: Record<string, string> = {
  LTR: "#f7768e",
  LINE: "#e0af68",
  DNA: "#9ece6a",
  RC: "#73daca",
  Satellite: "#bb9af7",
  Simple_repeat: "#7aa2f7",
  Low_complexity: "#565f89",
  rRNA: "#2ac3de",
  Unknown: "#a9b1d6",
  ARTEFACT: "#db4b4b",
};
export function repeatClassColor(kind: string): string {
  const cls = kind.split("/")[0].replace("?", "");
  return CLASS_COLORS[cls] ?? "#a9b1d6";
}

export function clamp(x: number, lo: number, hi: number): number {
  return Math.min(hi, Math.max(lo, x));
}

/** Debounce that also cancels the previous pending call's AbortController. */
export class Latest {
  private ctrl: AbortController | null = null;
  next(): AbortSignal {
    this.ctrl?.abort();
    this.ctrl = new AbortController();
    return this.ctrl.signal;
  }
}

export function isAbort(e: unknown): boolean {
  return e instanceof DOMException && e.name === "AbortError";
}

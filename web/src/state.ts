// Shared application state and the canvas helpers every panel uses.

import type { Genome, OrderMode, QueryStatus } from "./api";
import { Coords, clamp } from "./util";

export type ColorMode = "mismatch" | "strand";

export interface Query {
  id: number;
  kmer: string;
  maxN: number;
  status: QueryStatus | null;
}

/** The visible window: genome span × row span, plus the n and ordering being shown. */
export interface View {
  start: number;
  end: number;
  row0: number;
  row1: number;
  n: number;
  order: OrderMode;
}

export class State {
  readonly genome: Genome;
  readonly coords: Coords;
  view: View;
  /** Rows available at the current n (0 until a query is ready). */
  rows = 0;
  query: Query | null = null;
  color: ColorMode = "mismatch";
  /** Log-scale intensity in the raster and histogram. */
  logScale = true;
  signalTracks: string[] = [];
  featureSets: string[] = [];
  private listeners: Array<(why: string) => void> = [];

  constructor(genome: Genome) {
    this.genome = genome;
    this.coords = new Coords(genome.contigs);
    this.view = { start: 0, end: this.coords.length, row0: 0, row1: 0, n: 0, order: "per_n" };
  }

  on(f: (why: string) => void): void {
    this.listeners.push(f);
  }

  emit(why: string): void {
    for (const f of this.listeners) f(why);
  }

  /** Set the genome window, clamped; keeps at least 40 bp visible. */
  setRegion(start: number, end: number, why = "region"): void {
    const len = this.coords.length;
    let span = clamp(end - start, 40, len);
    let s = clamp(start, 0, len - span);
    if (!Number.isFinite(s)) s = 0;
    span = Math.min(span, len - s);
    this.view.start = Math.round(s);
    this.view.end = Math.round(s + span);
    this.emit(why);
  }

  /** Set the row window, clamped to the rows available. */
  setRows(row0: number, row1: number, why = "rows"): void {
    const rows = this.rows;
    if (rows === 0) {
      this.view.row0 = 0;
      this.view.row1 = 0;
    } else {
      const span = clamp(row1 - row0, Math.min(rows, 4), rows);
      const r0 = clamp(row0, 0, rows - span);
      this.view.row0 = Math.round(r0);
      this.view.row1 = Math.round(r0 + span);
    }
    this.emit(why);
  }

  /** Zoom the genome window by `factor` (< 1 zooms in) around fraction `at` of the width. */
  zoomRegion(factor: number, at = 0.5): void {
    const { start, end } = this.view;
    const span = end - start;
    const pivot = start + span * at;
    const ns = clamp(span * factor, 40, this.coords.length);
    this.setRegion(pivot - ns * at, pivot - ns * at + ns, "zoom");
  }

  zoomRows(factor: number, at = 0.5): void {
    const { row0, row1 } = this.view;
    const span = row1 - row0;
    const pivot = row0 + span * at;
    const ns = clamp(span * factor, Math.min(4, this.rows), this.rows);
    this.setRows(pivot - ns * at, pivot - ns * at + ns, "zoom");
  }

  bpPerPx(widthCss: number): number {
    return (this.view.end - this.view.start) / Math.max(1, widthCss);
  }
}

/** A canvas that tracks its CSS size and device pixel ratio. */
export class Surface {
  readonly canvas: HTMLCanvasElement;
  readonly ctx: CanvasRenderingContext2D;
  /** CSS pixels. */
  width = 0;
  height = 0;
  dpr = 1;

  constructor(parent: HTMLElement, className: string, onResize: () => void) {
    this.canvas = document.createElement("canvas");
    this.canvas.className = className;
    parent.appendChild(this.canvas);
    this.ctx = this.canvas.getContext("2d")!;
    new ResizeObserver(() => {
      const r = this.canvas.getBoundingClientRect();
      const dpr = Math.min(window.devicePixelRatio || 1, 2);
      if (r.width === this.width && r.height === this.height && dpr === this.dpr) return;
      this.width = r.width;
      this.height = r.height;
      this.dpr = dpr;
      this.canvas.width = Math.max(1, Math.round(r.width * dpr));
      this.canvas.height = Math.max(1, Math.round(r.height * dpr));
      onResize();
    }).observe(this.canvas);
  }

  /** Device-pixel size. */
  get pw(): number {
    return this.canvas.width;
  }
  get ph(): number {
    return this.canvas.height;
  }

  /** Reset the transform to CSS pixels and clear. */
  begin(): CanvasRenderingContext2D {
    const c = this.ctx;
    c.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);
    c.clearRect(0, 0, this.width, this.height);
    return c;
  }
}

/**
 * Keeps the last rendered image with the view it was rendered for, so pans and zooms can be
 * previewed instantly by drawing it transformed while the fresh image is fetched.
 */
export class Snapshot {
  private img: HTMLCanvasElement | null = null;
  private v: View | null = null;

  store(source: HTMLCanvasElement, v: View): void {
    if (!this.img) this.img = document.createElement("canvas");
    this.img.width = source.width;
    this.img.height = source.height;
    this.img.getContext("2d")!.drawImage(source, 0, 0);
    this.v = { ...v };
  }

  clear(): void {
    this.v = null;
  }

  /** Draw the stored image mapped into view `v` (rows only if `rowsToo`). Returns false if none. */
  draw(target: Surface, v: View, rowsToo: boolean): boolean {
    if (!this.img || !this.v) return false;
    const o = this.v;
    const W = target.pw;
    const H = target.ph;
    const span = v.end - v.start;
    const dx = ((o.start - v.start) / span) * W;
    const dw = ((o.end - o.start) / span) * W;
    let dy = 0;
    let dh = H;
    if (rowsToo && o.row1 > o.row0 && v.row1 > v.row0) {
      const rs = v.row1 - v.row0;
      dy = ((o.row0 - v.row0) / rs) * H;
      dh = ((o.row1 - o.row0) / rs) * H;
    }
    const c = target.ctx;
    c.setTransform(1, 0, 0, 1, 0, 0);
    c.clearRect(0, 0, W, H);
    c.imageSmoothingEnabled = false;
    c.drawImage(this.img, dx, dy, dw, dh);
    return true;
  }
}

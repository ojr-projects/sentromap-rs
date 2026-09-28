// The main view: variants (rows, cladogram order) × genome (columns), rendered server-side at
// screen resolution; plus the cladogram gutter to its left.

import { api, type Clade, type Variant } from "./api";
import { Snapshot, State, Surface } from "./state";
import {
  Latest,
  STRAND_BOTH,
  STRAND_MINUS,
  STRAND_PLUS,
  clamp,
  fmtInt,
  isAbort,
  mismatchColor,
  type RGB,
} from "./util";

const LEVELS = 64;

/** RGBA lookup: [mismatches or strand class][intensity level] → packed little-endian u32. */
function buildLut(state: State): Uint32Array[] {
  const classes: RGB[] =
    state.color === "strand"
      ? [STRAND_PLUS, STRAND_MINUS, STRAND_BOTH]
      : Array.from({ length: state.view.n + 1 }, (_, m) => mismatchColor(m, Math.max(1, state.view.n)));
  return classes.map((rgb) => {
    const lut = new Uint32Array(LEVELS);
    for (let l = 0; l < LEVELS; l++) {
      const f = 0.3 + (0.7 * l) / (LEVELS - 1);
      const [r, g, b] = rgb.map((x) => Math.round(x * f));
      lut[l] = (255 << 24) | (b << 16) | (g << 8) | r;
    }
    return lut;
  });
}

export class RasterPanel {
  readonly el: HTMLElement;
  private surf: Surface;
  private overlay: Surface;
  private snap = new Snapshot();
  private latest = new Latest();
  private timer = 0;
  private hover: { x: number; y: number } | null = null;
  private rowCache = new Map<number, Variant>();
  private rowCacheKey = "";
  lastRenderMs = 0;
  lastSites = 0;
  onHover: (info: HoverInfo | null) => void = () => {};
  onSelectRow: (row: number) => void = () => {};

  constructor(
    parent: HTMLElement,
    private state: State,
  ) {
    this.el = document.createElement("div");
    this.el.className = "raster";
    parent.appendChild(this.el);
    this.surf = new Surface(this.el, "raster-img", () => this.update("resize"));
    this.overlay = new Surface(this.el, "raster-overlay", () => this.drawOverlay());
    this.bindMouse();
  }

  /** Map a CSS y within the panel to a row. */
  rowAt(y: number): number {
    const { row0, row1 } = this.state.view;
    return Math.floor(row0 + (y / Math.max(1, this.overlay.height)) * (row1 - row0));
  }

  posAt(x: number): number {
    const { start, end } = this.state.view;
    return Math.floor(start + (x / Math.max(1, this.overlay.width)) * (end - start));
  }

  update(why: string): void {
    const s = this.state;
    if (!s.query?.status || s.query.status.status !== "ready" || s.rows === 0) {
      this.snap.clear();
      const c = this.surf.begin();
      c.fillStyle = "#6b7089";
      c.font = "14px system-ui, sans-serif";
      c.textAlign = "center";
      const msg = !s.query
        ? "Enter a 31-mer (or double-click the genome) to search"
        : s.query.status?.status === "ready"
          ? `No variants within n = ${s.view.n}`
          : "Searching…";
      c.fillText(msg, this.surf.width / 2, this.surf.height / 2);
      this.drawOverlay();
      return;
    }
    // Instant preview for pans and zooms; fresh render shortly after.
    if (why === "zoom" || why === "pan" || why === "region" || why === "rows") {
      this.snap.draw(this.surf, s.view, true);
    }
    this.drawOverlay();
    clearTimeout(this.timer);
    this.timer = window.setTimeout(() => this.fetch(), why === "pan" ? 60 : 30);
  }

  private async fetch(): Promise<void> {
    const s = this.state;
    const q = s.query;
    if (!q) return;
    const v = { ...s.view };
    const signal = this.latest.next();
    try {
      const r = await api.raster(q.id, v, this.surf.pw, this.surf.ph, signal);
      if (r.width !== this.surf.pw || r.height !== this.surf.ph) return; // resized meanwhile
      this.paint(r.cells, r.width, r.height);
      this.snap.store(this.surf.canvas, v);
      this.lastRenderMs = r.renderMs;
      this.lastSites = r.sitesDrawn;
      s.emit("rendered");
    } catch (e) {
      if (!isAbort(e)) console.error(e);
    }
  }

  private paint(cells: Uint8Array, w: number, h: number): void {
    const s = this.state;
    const n = w * h;
    let max = 1;
    for (let i = 0; i < n; i++) {
      const c = cells[4 * i] | (cells[4 * i + 1] << 8);
      if (c > max) max = c;
    }
    const lut = buildLut(s);
    const img = new ImageData(w, h);
    const out = new Uint32Array(img.data.buffer);
    const logMax = Math.log1p(max);
    const strand = s.color === "strand";
    for (let i = 0; i < n; i++) {
      const c = cells[4 * i] | (cells[4 * i + 1] << 8);
      if (c === 0) continue;
      const t = s.logScale ? Math.log1p(c) / logMax : c / max;
      const level = Math.min(LEVELS - 1, Math.round(t * (LEVELS - 1)));
      const cls = strand ? [0, 0, 1, 2][cells[4 * i + 3] & 3] : Math.min(cells[4 * i + 2], lut.length - 1);
      out[i] = lut[cls][level];
    }
    const ctx = this.surf.ctx;
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.putImageData(img, 0, 0);
  }

  /** Contig boundaries and the hover crosshair. */
  drawOverlay(): void {
    const c = this.overlay.begin();
    const s = this.state;
    const { start, end } = s.view;
    const span = end - start;
    const W = this.overlay.width;
    const H = this.overlay.height;
    // Contig boundaries, when they are far enough apart to be meaningful.
    const contigs = s.coords.contigs;
    c.strokeStyle = "rgba(160,170,210,0.18)";
    c.lineWidth = 1;
    let lastX = -10;
    for (let i = s.coords.locate(start); i < contigs.length && contigs[i].offset < end; i++) {
      const x = ((contigs[i].offset - start) / span) * W;
      if (x - lastX < 6 || x <= 0) continue;
      c.beginPath();
      c.moveTo(Math.round(x) + 0.5, 0);
      c.lineTo(Math.round(x) + 0.5, H);
      c.stroke();
      lastX = x;
    }
    if (this.hover) {
      c.strokeStyle = "rgba(255,255,255,0.35)";
      c.setLineDash([3, 3]);
      c.beginPath();
      c.moveTo(Math.round(this.hover.x) + 0.5, 0);
      c.lineTo(Math.round(this.hover.x) + 0.5, H);
      const rows = s.view.row1 - s.view.row0;
      if (rows > 0) {
        const row = this.rowAt(this.hover.y);
        const y0 = ((row - s.view.row0) / rows) * H;
        const y1 = ((row + 1 - s.view.row0) / rows) * H;
        c.moveTo(0, Math.round((y0 + y1) / 2) + 0.5);
        c.lineTo(W, Math.round((y0 + y1) / 2) + 0.5);
      }
      c.stroke();
      c.setLineDash([]);
    }
  }

  /** The variant at a row (cached per query/n/order). */
  async variantAt(row: number): Promise<Variant | null> {
    const s = this.state;
    const q = s.query;
    if (!q || row < 0 || row >= s.rows) return null;
    const key = `${q.id}/${s.view.n}/${s.view.order}`;
    if (key !== this.rowCacheKey) {
      this.rowCache.clear();
      this.rowCacheKey = key;
    }
    const hit = this.rowCache.get(row);
    if (hit) return hit;
    const base = Math.max(0, row - 20);
    const r = await api.variants(q.id, s.view.n, s.view.order, base, 60);
    for (const v of r.variants) this.rowCache.set(v.row, v);
    if (this.rowCache.size > 20000) this.rowCache.clear();
    return this.rowCache.get(row) ?? null;
  }

  private bindMouse(): void {
    const el = this.overlay.canvas;
    const s = this.state;
    let drag: { x: number; y: number; start: number; end: number; row0: number; row1: number; moved: boolean } | null = null;

    el.addEventListener(
      "wheel",
      (e) => {
        e.preventDefault();
        const r = el.getBoundingClientRect();
        if (e.altKey) {
          window.dispatchEvent(new CustomEvent("sm-step-n", { detail: e.deltaY > 0 ? 1 : -1 }));
          return;
        }
        const f = Math.exp(clamp(e.deltaY, -200, 200) * 0.002);
        if (e.shiftKey) {
          s.zoomRows(f, (e.clientY - r.top) / r.height);
        } else if (Math.abs(e.deltaX) > Math.abs(e.deltaY)) {
          const span = s.view.end - s.view.start;
          const d = (e.deltaX / r.width) * span;
          s.setRegion(s.view.start + d, s.view.end + d, "pan");
        } else {
          s.zoomRegion(f, (e.clientX - r.left) / r.width);
        }
      },
      { passive: false },
    );

    el.addEventListener("pointerdown", (e) => {
      el.setPointerCapture(e.pointerId);
      drag = { x: e.clientX, y: e.clientY, ...s.view, moved: false };
    });
    el.addEventListener("pointermove", (e) => {
      const r = el.getBoundingClientRect();
      if (drag) {
        const dx = e.clientX - drag.x;
        const dy = e.clientY - drag.y;
        if (Math.abs(dx) + Math.abs(dy) > 3) drag.moved = true;
        if (drag.moved) {
          const span = drag.end - drag.start;
          const d = (-dx / r.width) * span;
          s.setRegion(drag.start + d, drag.end + d, "pan");
          const rs = drag.row1 - drag.row0;
          if (rs > 0) {
            const dr = (-dy / r.height) * rs;
            s.setRows(drag.row0 + dr, drag.row1 + dr, "pan");
          }
        }
      }
      this.hover = { x: e.clientX - r.left, y: e.clientY - r.top };
      this.drawOverlay();
      this.emitHover();
    });
    el.addEventListener("pointerup", (e) => {
      if (drag && !drag.moved) {
        const r = el.getBoundingClientRect();
        this.onSelectRow(this.rowAt(e.clientY - r.top));
      }
      drag = null;
    });
    el.addEventListener("pointerleave", () => {
      this.hover = null;
      this.drawOverlay();
      this.onHover(null);
    });
    el.addEventListener("dblclick", (e) => {
      const r = el.getBoundingClientRect();
      s.zoomRegion(0.25, (e.clientX - r.left) / r.width);
    });
  }

  private hoverSeq = 0;
  private async emitHover(): Promise<void> {
    if (!this.hover) return;
    const seq = ++this.hoverSeq;
    const pos = this.posAt(this.hover.x);
    const s = this.state;
    const row = s.rows > 0 ? this.rowAt(this.hover.y) : -1;
    this.onHover({ pos, row, variant: null });
    if (row >= 0) {
      try {
        const v = await this.variantAt(row);
        if (seq === this.hoverSeq) this.onHover({ pos, row, variant: v });
      } catch {
        /* transient */
      }
    }
  }
}

export interface HoverInfo {
  pos: number;
  row: number;
  variant: Variant | null;
}

/** Cladogram brackets for the rows in view, left of the raster. */
export class TreeGutter {
  readonly el: HTMLElement;
  private surf: Surface;
  private latest = new Latest();
  private timer = 0;
  private clades: Clade[] = [];
  private hovered: Clade | null = null;
  treeLength = 0;
  starLength = 0;
  onHoverClade: (c: Clade | null) => void = () => {};

  constructor(
    parent: HTMLElement,
    private state: State,
  ) {
    this.el = document.createElement("div");
    this.el.className = "tree";
    parent.appendChild(this.el);
    this.surf = new Surface(this.el, "tree-canvas", () => this.draw());
    this.bind();
  }

  update(why: string): void {
    if (why === "region" || why === "rendered") return; // independent of the genome window
    this.draw();
    clearTimeout(this.timer);
    this.timer = window.setTimeout(() => this.fetch(), 40);
  }

  private async fetch(): Promise<void> {
    const s = this.state;
    const q = s.query;
    if (!q?.status || q.status.status !== "ready" || s.rows === 0) {
      this.clades = [];
      this.draw();
      return;
    }
    const v = { ...s.view };
    const rowsPerPx = (v.row1 - v.row0) / Math.max(1, this.surf.height);
    try {
      const t = await api.tree(q.id, v, Math.max(2, Math.ceil(rowsPerPx * 4)), this.latest.next());
      this.clades = t.clades;
      this.treeLength = t.tree_length;
      this.starLength = t.star_length;
      this.draw();
      s.emit("tree");
    } catch (e) {
      if (!isAbort(e)) console.error(e);
    }
  }

  private y(row: number): number {
    const { row0, row1 } = this.state.view;
    return ((row - row0) / Math.max(1, row1 - row0)) * this.surf.height;
  }

  private x(depth: number): number {
    return Math.min(this.surf.width - 3, 4 + (depth - 1) * 6);
  }

  draw(): void {
    const c = this.surf.begin();
    const W = this.surf.width;
    const H = this.surf.height;
    c.fillStyle = "#12141d";
    c.fillRect(0, 0, W, H);
    const s = this.state;
    if (s.rows === 0) return;
    const maxD = Math.max(1, ...this.clades.map((k) => k.depth));
    for (const k of this.clades) {
      const y0 = Math.max(0, this.y(k.start));
      const y1 = Math.min(H, this.y(k.end));
      if (y1 - y0 < 1) continue;
      const x = this.x(k.depth);
      const hot = k === this.hovered;
      const t = (k.depth - 1) / maxD;
      c.strokeStyle = hot ? "#ffffff" : `hsla(${200 + 120 * t}, 55%, ${62 - 18 * t}%, ${hot ? 1 : 0.85})`;
      c.lineWidth = hot ? 2 : 1;
      c.beginPath();
      c.moveTo(x + 3, y0 + 0.5);
      c.lineTo(x + 0.5, y0 + 0.5);
      c.lineTo(x + 0.5, y1 - 0.5);
      c.lineTo(x + 3, y1 - 0.5);
      c.stroke();
    }
    // Row ticks when rows are tall enough to read.
    const rowPx = H / Math.max(1, s.view.row1 - s.view.row0);
    if (rowPx >= 12) {
      c.fillStyle = "#7c8199";
      c.font = "10px ui-monospace, monospace";
      c.textAlign = "right";
      for (let r = s.view.row0; r < s.view.row1; r++) c.fillText(fmtInt(r), W - 4, this.y(r) + rowPx / 2 + 3);
    }
  }

  private bind(): void {
    const el = this.surf.canvas;
    const s = this.state;
    const find = (e: PointerEvent | MouseEvent): Clade | null => {
      const r = el.getBoundingClientRect();
      const x = e.clientX - r.left;
      const y = e.clientY - r.top;
      const { row0, row1 } = s.view;
      const row = row0 + (y / r.height) * (row1 - row0);
      // Deepest clade containing the row whose bracket is left of the cursor.
      let best: Clade | null = null;
      for (const k of this.clades) {
        if (row >= k.start && row < k.end && this.x(k.depth) <= x + 3 && (!best || k.depth > best.depth)) best = k;
      }
      return best;
    };
    el.addEventListener("pointermove", (e) => {
      const k = find(e);
      if (k !== this.hovered) {
        this.hovered = k;
        this.draw();
        this.onHoverClade(k);
      }
    });
    el.addEventListener("pointerleave", () => {
      this.hovered = null;
      this.draw();
      this.onHoverClade(null);
    });
    el.addEventListener("click", (e) => {
      const k = find(e);
      if (k) {
        const pad = Math.max(1, Math.round((k.end - k.start) * 0.05));
        s.setRows(k.start - pad, k.end + pad);
      }
    });
    el.addEventListener(
      "wheel",
      (e) => {
        e.preventDefault();
        const r = el.getBoundingClientRect();
        const f = Math.exp(clamp(e.deltaY, -200, 200) * 0.002);
        s.zoomRows(f, (e.clientY - r.top) / r.height);
      },
      { passive: false },
    );
  }
}


// Genome-axis panels above the raster: ruler, annotation tracks and the mismatch histogram.

import { api, type FeatureSetInfo, type FeaturesResponse, type TrackInfo } from "./api";
import { Snapshot, State, Surface } from "./state";
import { Latest, clamp, css, fmtBp, fmtInt, isAbort, kindColor, mismatchColor, niceStep, repeatClassColor, signalRamp } from "./util";

/** Double-click anywhere on the genome axis picks the 31-mer starting there. */
export type PickHandler = (pos: number) => void;

function bindGenomeAxis(el: HTMLElement, state: State, onPick: PickHandler): void {
  let drag: { x: number; start: number; end: number } | null = null;
  el.addEventListener(
    "wheel",
    (e) => {
      e.preventDefault();
      const r = el.getBoundingClientRect();
      if (Math.abs(e.deltaX) > Math.abs(e.deltaY)) {
        const d = (e.deltaX / r.width) * (state.view.end - state.view.start);
        state.setRegion(state.view.start + d, state.view.end + d, "pan");
      } else {
        state.zoomRegion(Math.exp(clamp(e.deltaY, -200, 200) * 0.002), (e.clientX - r.left) / r.width);
      }
    },
    { passive: false },
  );
  el.addEventListener("pointerdown", (e) => {
    el.setPointerCapture(e.pointerId);
    drag = { x: e.clientX, start: state.view.start, end: state.view.end };
  });
  el.addEventListener("pointermove", (e) => {
    if (!drag) return;
    const r = el.getBoundingClientRect();
    const d = (-(e.clientX - drag.x) / r.width) * (drag.end - drag.start);
    state.setRegion(drag.start + d, drag.end + d, "pan");
  });
  el.addEventListener("pointerup", () => (drag = null));
  el.addEventListener("dblclick", (e) => {
    const r = el.getBoundingClientRect();
    const { start, end } = state.view;
    onPick(Math.floor(start + ((e.clientX - r.left) / r.width) * (end - start)));
  });
}

export class Ruler {
  private surf: Surface;
  constructor(
    parent: HTMLElement,
    private state: State,
    onPick: PickHandler,
  ) {
    this.surf = new Surface(parent, "ruler", () => this.draw());
    bindGenomeAxis(this.surf.canvas, state, onPick);
  }

  update(): void {
    this.draw();
  }

  draw(): void {
    const c = this.surf.begin();
    const s = this.state;
    const W = this.surf.width;
    const { start, end } = s.view;
    const span = end - start;
    const X = (p: number) => ((p - start) / span) * W;
    const contigs = s.coords.contigs;
    // Contig band.
    c.font = "11px system-ui, sans-serif";
    c.textBaseline = "middle";
    let alt = 0;
    for (let i = s.coords.locate(start); i < contigs.length && contigs[i].offset < end; i++) {
      const k = contigs[i];
      const x0 = Math.max(0, X(k.offset));
      const x1 = Math.min(W, X(k.offset + k.length));
      if (x1 - x0 < 0.5) continue;
      c.fillStyle = alt++ % 2 ? "#2a2f45" : "#343a55";
      c.fillRect(x0, 0, Math.max(1, x1 - x0), 14);
      if (x1 - x0 > 30) {
        c.fillStyle = "#c8cde6";
        c.save();
        c.beginPath();
        c.rect(x0, 0, x1 - x0, 14);
        c.clip();
        c.fillText(k.display, x0 + 4, 7.5);
        c.restore();
      }
    }
    // Ticks in contig-local coordinates when one contig fills the view, else global.
    const ci = s.coords.locate(start);
    const cj = s.coords.locate(Math.max(start, end - 1));
    const base = ci === cj ? contigs[ci].offset : 0;
    const step = niceStep(span, W / 110);
    c.strokeStyle = "#5b6283";
    c.fillStyle = "#9aa1c0";
    c.textBaseline = "alphabetic";
    const first = Math.ceil((start - base) / step) * step;
    for (let t = first; base + t < end; t += step) {
      const x = Math.round(X(base + t)) + 0.5;
      c.beginPath();
      c.moveTo(x, 18);
      c.lineTo(x, 24);
      c.stroke();
      c.fillText(ci === cj ? fmtTick(t) : fmtBp(t), x + 3, 31);
    }
  }
}

function fmtTick(bp: number): string {
  if (bp === 0) return "0";
  if (bp % 1e6 === 0) return bp / 1e6 + " Mb";
  if (bp % 1e3 === 0 && bp < 1e7) return fmtInt(bp / 1e3) + " kb";
  return fmtInt(bp);
}

interface Row {
  kind: "sequence" | "signal" | "features";
  name: string;
  height: number;
  label: string;
  title: string;
}

/** Stacked annotation rows on one canvas; labels are HTML in the gutter. */
export class TrackPanel {
  readonly el: HTMLElement;
  private surf: Surface;
  private snap = new Snapshot();
  private latest = new Latest();
  private timer = 0;
  private rows: Row[] = [];
  private tracks = new Map<string, TrackInfo>();
  private sets = new Map<string, FeatureSetInfo>();
  labels: HTMLElement;
  onLayout: () => void = () => {};

  constructor(
    parent: HTMLElement,
    labelParent: HTMLElement,
    private state: State,
    onPick: PickHandler,
  ) {
    this.el = document.createElement("div");
    this.el.className = "tracks";
    parent.appendChild(this.el);
    this.labels = document.createElement("div");
    this.labels.className = "track-labels";
    labelParent.appendChild(this.labels);
    this.surf = new Surface(this.el, "tracks-canvas", () => this.update("resize"));
    bindGenomeAxis(this.surf.canvas, state, onPick);
  }

  setCatalogue(tracks: TrackInfo[], sets: FeatureSetInfo[]): void {
    for (const t of tracks) this.tracks.set(t.name, t);
    for (const f of sets) this.sets.set(f.name, f);
  }

  private layout(): void {
    const s = this.state;
    const rows: Row[] = [];
    if (s.bpPerPx(this.surf.width || 1000) < 2) rows.push({ kind: "sequence", name: "sequence", height: 14, label: "sequence", title: "Reference sequence" });
    for (const name of s.featureSets) {
      const f = this.sets.get(name);
      if (f) rows.push({ kind: "features", name, height: 30, label: name, title: `${f.description} (${fmtInt(f.count)} features)` });
    }
    for (const name of s.signalTracks) {
      const t = this.tracks.get(name);
      if (t) rows.push({ kind: "signal", name, height: 20, label: name, title: `${t.description} — ${t.source}` });
    }
    const changed = rows.map((r) => r.name).join() !== this.rows.map((r) => r.name).join();
    this.rows = rows;
    if (changed) {
      this.labels.innerHTML = "";
      for (const r of rows) {
        const d = document.createElement("div");
        d.className = `track-label ${r.kind}`;
        d.style.height = r.height + "px";
        d.textContent = r.label;
        d.title = r.title;
        this.labels.appendChild(d);
      }
      const total = rows.reduce((a, r) => a + r.height, 0);
      this.el.style.height = total + "px";
      this.labels.style.height = total + "px";
      this.snap.clear();
      this.onLayout();
    }
  }

  update(why: string): void {
    if (why === "rows" || why === "tree" || why === "rendered") return;
    this.layout();
    if (why === "zoom" || why === "pan" || why === "region") this.snap.draw(this.surf, this.state.view, false);
    clearTimeout(this.timer);
    this.timer = window.setTimeout(() => this.fetch(), why === "pan" ? 60 : 30);
  }

  private async fetch(): Promise<void> {
    const s = this.state;
    const v = { ...s.view };
    const W = this.surf.pw;
    const signal = this.latest.next();
    try {
      const data = await Promise.all(
        this.rows.map((r) => {
          if (r.kind === "signal") return api.trackSummary(r.name, v.start, v.end, W, signal);
          if (r.kind === "features") {
            const kinds = r.name === "genes" ? (v.end - v.start < 200_000 ? "gene,exon,pseudogene,ncRNA_gene" : "gene,pseudogene,ncRNA_gene") : undefined;
            return api.features(r.name, v.start, v.end, Math.ceil(W / 2), 1500, kinds, signal);
          }
          return v.end - v.start <= 20_000 ? api.sequence(v.start, v.end, signal) : Promise.resolve("");
        }),
      );
      this.draw(v, data);
      this.snap.store(this.surf.canvas, v);
    } catch (e) {
      if (!isAbort(e)) console.error(e);
    }
  }

  private draw(v: { start: number; end: number }, data: Array<Float32Array | FeaturesResponse | string>): void {
    const c = this.surf.begin();
    const W = this.surf.width;
    const span = v.end - v.start;
    const X = (p: number) => ((p - v.start) / span) * W;
    let y = 0;
    this.rows.forEach((r, i) => {
      const d = data[i];
      c.save();
      c.beginPath();
      c.rect(0, y, W, r.height);
      c.clip();
      c.fillStyle = i % 2 ? "#161927" : "#191c2b";
      c.fillRect(0, y, W, r.height);
      if (r.kind === "signal") this.drawSignal(c, d as Float32Array, this.tracks.get(r.name)!, y, r.height, W);
      else if (r.kind === "features") this.drawFeatures(c, d as FeaturesResponse, r.name, y, r.height, W, X);
      else this.drawSequence(c, d as string, y, r.height, X, v.start);
      c.restore();
      y += r.height;
    });
  }

  private drawSignal(c: CanvasRenderingContext2D, vals: Float32Array, t: TrackInfo, y: number, h: number, W: number): void {
    const cols = vals.length / 2;
    const lo = Math.min(0, t.p01);
    const hi = t.p99 > lo ? t.p99 : t.max;
    const colW = W / cols;
    for (let i = 0; i < cols; i++) {
      const m = vals[2 * i];
      if (Number.isNaN(m)) continue;
      const f = clamp((m - lo) / (hi - lo), 0, 1);
      const bh = Math.max(1, f * (h - 2));
      c.fillStyle = css(signalRamp(0.15 + 0.85 * f));
      c.fillRect(i * colW, y + h - bh, Math.max(colW, 0.6), bh);
    }
  }

  private drawFeatures(
    c: CanvasRenderingContext2D,
    d: FeaturesResponse,
    name: string,
    y: number,
    h: number,
    W: number,
    X: (p: number) => number,
  ): void {
    if (d.mode === "density") {
      const n = d.density.length;
      const max = Math.max(1, ...d.density);
      const colW = W / n;
      c.fillStyle = name === "RepeatMasker" ? "#bb9af7" : name === "TRF" ? "#7aa2f7" : "#9ece6a";
      for (let i = 0; i < n; i++) {
        if (!d.density[i]) continue;
        const f = Math.log1p(d.density[i]) / Math.log1p(max);
        const bh = Math.max(1, f * (h - 4));
        c.globalAlpha = 0.4 + 0.6 * f;
        c.fillRect(i * colW, y + h - 2 - bh, Math.max(colW, 0.6), bh);
      }
      c.globalAlpha = 1;
      return;
    }
    // Glyphs packed greedily into lanes.
    const laneH = 9;
    const lanes = Math.floor((h - 2) / laneH);
    const laneEnd: number[] = new Array(lanes).fill(-Infinity);
    c.font = "9px system-ui, sans-serif";
    c.textBaseline = "middle";
    const sorted = [...d.features].sort((a, b) => (a.kind === b.kind ? a.start - b.start : a.start - b.start));
    for (const f of sorted) {
      const kind = d.kinds[f.kind] ?? "";
      const x0 = X(f.start);
      const x1 = Math.max(x0 + 1, X(f.end));
      let labelW = x1 - x0 > 40 ? 0 : Math.min(80, c.measureText(f.name).width + 4);
      let lane = laneEnd.findIndex((e) => e < x0 - 1);
      if (lane < 0) {
        // No free lane: overlap in the least-full lane, without a label.
        lane = laneEnd.indexOf(Math.min(...laneEnd));
        labelW = 0;
      }
      laneEnd[lane] = Math.max(laneEnd[lane], x1 + labelW);
      const ly = y + 1 + lane * laneH;
      const color = name === "RepeatMasker" ? repeatClassColor(kind) : name === "genes" ? (kind === "exon" ? "#e0af68" : "#9ece6a") : kindColor(f.kind);
      c.fillStyle = color;
      const thin = name === "genes" && kind !== "exon" && d.features.some((g) => d.kinds[g.kind] === "exon");
      c.fillRect(x0, ly + (thin ? 3 : 1), x1 - x0, thin ? 2 : laneH - 2);
      if (x1 - x0 > 40) {
        c.fillStyle = "#0f111a";
        c.save();
        c.beginPath();
        c.rect(x0, ly, x1 - x0, laneH);
        c.clip();
        if (!thin) c.fillText(`${f.strand === "-" ? "◂ " : ""}${f.name}${f.strand === "+" ? " ▸" : ""}`, x0 + 2, ly + laneH / 2);
        c.restore();
        if (thin) {
          c.fillStyle = "#c8cde6";
          c.fillText(f.name, Math.max(0, x0) + 2, ly + laneH / 2);
        }
      } else if (labelW) {
        c.fillStyle = "#9aa1c0";
        c.fillText(f.name, x1 + 2, ly + laneH / 2);
      }
    }
  }

  private drawSequence(c: CanvasRenderingContext2D, seq: string, y: number, h: number, X: (p: number) => number, start: number): void {
    if (!seq) return;
    const colors: Record<string, string> = { A: "#5fb86a", C: "#4c8fe0", G: "#e0b64c", T: "#e05c5c" };
    const pxPerBase = X(start + 1) - X(start);
    c.font = "bold 11px ui-monospace, monospace";
    c.textAlign = "center";
    c.textBaseline = "middle";
    for (let i = 0; i < seq.length; i++) {
      const b = seq[i];
      const col = colors[b.toUpperCase()] ?? "#555a70";
      const x = X(start + i);
      if (pxPerBase >= 9) {
        c.fillStyle = col;
        c.globalAlpha = b === b.toLowerCase() ? 0.55 : 1;
        c.fillText(b, x + pxPerBase / 2, y + h / 2);
      } else {
        c.fillStyle = col;
        c.globalAlpha = b === b.toLowerCase() ? 0.5 : 0.9;
        c.fillRect(x, y + 3, Math.max(0.8, pxPerBase), h - 6);
      }
    }
    c.globalAlpha = 1;
    c.textAlign = "start";
  }
}

/** Sites per column, stacked by mismatch count. */
export class Histogram {
  private surf: Surface;
  private latest = new Latest();
  private timer = 0;
  private snap = new Snapshot();

  constructor(
    parent: HTMLElement,
    private state: State,
    onPick: PickHandler,
  ) {
    this.surf = new Surface(parent, "histogram", () => this.update("resize"));
    bindGenomeAxis(this.surf.canvas, state, onPick);
  }

  update(why: string): void {
    if (why === "rows" || why === "tree" || why === "rendered") return;
    const s = this.state;
    if (!s.query?.status || s.query.status.status !== "ready") {
      this.snap.clear();
      this.surf.begin();
      return;
    }
    if (why === "zoom" || why === "pan" || why === "region") this.snap.draw(this.surf, s.view, false);
    clearTimeout(this.timer);
    this.timer = window.setTimeout(() => this.fetch(), 40);
  }

  private async fetch(): Promise<void> {
    const s = this.state;
    const q = s.query;
    if (!q) return;
    const v = { ...s.view };
    const width = Math.max(1, Math.floor(this.surf.width));
    try {
      const { stride, counts } = await api.histogram(q.id, v, width, this.latest.next());
      const c = this.surf.begin();
      const H = this.surf.height;
      let max = 1;
      const totals = new Float64Array(width);
      for (let i = 0; i < width; i++) {
        let t = 0;
        for (let m = 0; m < stride; m++) t += counts[i * stride + m];
        totals[i] = t;
        if (t > max) max = t;
      }
      const scale = (x: number) => (s.logScale ? Math.log1p(x) / Math.log1p(max) : x / max) * (H - 2);
      const colors = Array.from({ length: stride }, (_, m) => css(mismatchColor(m, Math.max(1, stride - 1))));
      for (let i = 0; i < width; i++) {
        if (!totals[i]) continue;
        const h = scale(totals[i]);
        let y = H;
        // Split the (possibly log-scaled) bar proportionally by mismatch count.
        for (let m = 0; m < stride; m++) {
          const n = counts[i * stride + m];
          if (!n) continue;
          const bh = (n / totals[i]) * h;
          c.fillStyle = colors[m];
          c.fillRect(i, y - bh, 1, bh);
          y -= bh;
        }
      }
      this.snap.store(this.surf.canvas, v);
    } catch (e) {
      if (!isAbort(e)) console.error(e);
    }
  }
}

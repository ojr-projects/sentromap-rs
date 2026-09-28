// sentromap web front end: wiring, query lifecycle, n-bar, side panel and URL state.

import { api, ApiError, type Clade, type FeatureSetInfo, type OrderMode, type QueryStatus, type TrackInfo, type Variant } from "./api";
import { RasterPanel, TreeGutter, type HoverInfo } from "./raster";
import { State, type ColorMode } from "./state";
import { Histogram, Ruler, TrackPanel } from "./tracks";
import { K, css, fmtBp, fmtCount, fmtExpected, fmtInt, fmtMs, mismatchColor } from "./util";

const $ = <T extends HTMLElement = HTMLElement>(sel: string) => document.querySelector(sel) as T;

function h<K extends keyof HTMLElementTagNameMap>(tag: K, attrs: Record<string, string> = {}, ...kids: (Node | string)[]): HTMLElementTagNameMap[K] {
  const e = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs)) {
    if (k === "class") e.className = v;
    else e.setAttribute(k, v);
  }
  for (const k of kids) e.append(k);
  return e;
}

/** A variant with its mismatches highlighted, coloured by position in the k-mer. */
function variantHtml(kmer: string, positions: number[]): HTMLElement {
  const set = new Set(positions);
  const span = h("span", { class: "seq" });
  for (let i = 0; i < kmer.length; i++) {
    const b = h("span", set.has(i + 1) ? { class: "mm" } : {}, kmer[i]);
    span.append(b);
  }
  return span;
}

const DEFAULT_SIGNAL = ["GC", "phastCons124way", "H3K9me3", "HP1a", "H3K27me3", "H3K4me3", "ATAC"];
const DEFAULT_FEATURES = ["genes", "RepeatMasker", "TRF"];

async function main(): Promise<void> {
  const genome = await api.genome();
  const catalogue = await api.tracks();
  const state = new State(genome);
  state.signalTracks = DEFAULT_SIGNAL.filter((n) => catalogue.tracks.some((t) => t.name === n));
  state.featureSets = DEFAULT_FEATURES.filter((n) => catalogue.feature_sets.some((f) => f.name === n));
  document.title = `sentromap · ${genome.name}`;
  $("#genome-name").textContent = genome.name;
  $("#genome-stats").textContent = `${fmtBp(genome.genome_len)} · ${fmtInt(genome.contigs.length)} contigs · ${fmtCount(genome.distinct_kmers)} distinct 31-mers`;

  // --- layout ---
  const gutter = $("#gutter");
  const main = $("#main");
  const ruler = new Ruler(main, state, pick);
  const tracks = new TrackPanel(main, gutter, state, pick);
  tracks.setCatalogue(catalogue.tracks, catalogue.feature_sets);
  const histLabel = h("div", { class: "track-label hist-label", title: "Sites per column, stacked by mismatch count" }, "sites by n");
  gutter.append(histLabel);
  const hist = new Histogram(main, state, pick);
  const body = h("div", { class: "body" });
  main.append(body);
  const tree = new TreeGutter(gutter, state);
  const raster = new RasterPanel(body, state);

  // --- query controls ---
  const input = $<HTMLInputElement>("#query");
  const maxNSel = $<HTMLSelectElement>("#max-n");
  const orderSel = $<HTMLSelectElement>("#order");
  const colorSel = $<HTMLSelectElement>("#color");
  const logBox = $<HTMLInputElement>("#log");
  for (let n = 0; n <= genome.max_n_limit; n++) maxNSel.append(h("option", { value: String(n) }, String(n)));
  maxNSel.value = "12";

  let pollTimer = 0;
  let clade: Clade | null = null;
  let hover: HoverInfo | null = null;
  let selectedRow: number | null = null;

  async function run(body: { kmer?: string; pos?: number; contig?: string; start?: number }): Promise<void> {
    const maxN = Number(maxNSel.value);
    const order = orderSel.value as OrderMode;
    setStatus("submitting…");
    try {
      if (state.query && state.query.status?.status !== "ready") api.cancel(state.query.id).catch(() => {});
      const c = await api.createQuery({ ...body, max_n: maxN, order });
      input.value = c.kmer;
      state.query = { id: c.id, kmer: c.kmer, maxN: c.max_n, status: null };
      state.view.order = order;
      state.view.n = Math.min(state.view.n || c.max_n, c.max_n);
      if (!body.kmer) state.view.n = c.max_n;
      selectedRow = null;
      setStatus(c.cached ? "cached" : `${c.engine}, estimated ${fmtMs(c.estimate_ms)} · chance hits at n=${c.max_n}: ${fmtExpected(c.expected_chance_hits)}`);
      writeHash();
      poll();
    } catch (e) {
      setStatus(e instanceof ApiError ? e.message : String(e), true);
    }
  }

  async function poll(): Promise<void> {
    clearTimeout(pollTimer);
    const q = state.query;
    if (!q) return;
    try {
      const st = await api.status(q.id);
      if (state.query !== q) return;
      q.status = st;
      if (st.status === "ready") {
        onReady(st);
      } else if (st.status === "failed" || st.status === "cancelled") {
        setStatus(`${st.status}${st.error ? ": " + st.error : ""}`, true);
      } else {
        setStatus(`${st.status}… ${fmtMs(st.elapsed_ms)} (estimate ${fmtMs(st.estimate_ms)})`);
        pollTimer = window.setTimeout(poll, 150);
      }
    } catch (e) {
      setStatus(String(e), true);
    }
  }

  function onReady(st: QueryStatus): void {
    const t = st.timings_ms;
    const parts = st.cached ? ["cached"] : [`${st.engine} ${fmtMs(t.search ?? 0)}`, `assemble ${fmtMs(t.assemble ?? 0)}`, `order ${fmtMs(t.order ?? 0)}`];
    setStatus(`ready · ${parts.join(" · ")}`);
    setN(Math.min(state.view.n, st.max_n), true);
    renderNBar();
    renderQueryInfo();
  }

  function setN(n: number, force = false): void {
    const st = state.query?.status;
    if (!st?.counts) return;
    n = Math.max(0, Math.min(st.max_n, n));
    if (n === state.view.n && !force) return;
    state.view.n = n;
    state.rows = st.counts[n].variants;
    state.view.row0 = 0;
    state.view.row1 = state.rows;
    selectedRow = null;
    renderNBar();
    writeHash();
    state.emit("n");
  }

  window.addEventListener("sm-step-n", (e) => setN(state.view.n + (e as CustomEvent<number>).detail));

  async function pick(pos: number): Promise<void> {
    const k = await api.kmerAt(pos);
    if (!k.valid) {
      setStatus(`no valid 31-mer at ${k.display}:${fmtInt(k.start)}`, true);
      return;
    }
    run({ pos });
  }

  $("#query-form").addEventListener("submit", (e) => {
    e.preventDefault();
    const v = input.value.trim();
    const m = /^([^:\s]+):([\d,]+)$/.exec(v);
    if (m) run({ contig: m[1], start: Number(m[2].replace(/,/g, "")) });
    else run({ kmer: v.replace(/\s+/g, "") });
  });
  orderSel.addEventListener("change", () => {
    state.view.order = orderSel.value as OrderMode;
    writeHash();
    state.emit("order");
  });
  colorSel.addEventListener("change", () => {
    state.color = colorSel.value as ColorMode;
    state.emit("color");
  });
  logBox.addEventListener("change", () => {
    state.logScale = logBox.checked;
    state.emit("color");
  });

  // --- region box ---
  const regionInput = $<HTMLInputElement>("#region");
  $("#region-form").addEventListener("submit", (e) => {
    e.preventDefault();
    const r = state.coords.parseRegion(regionInput.value);
    if (r) state.setRegion(r[0], r[1]);
    else regionInput.classList.add("bad");
  });
  regionInput.addEventListener("input", () => regionInput.classList.remove("bad"));
  const contigSel = $<HTMLSelectElement>("#contig");
  contigSel.append(h("option", { value: "all" }, "whole genome"));
  for (const c of [...genome.contigs].filter((c) => c.length >= 1_000_000).sort((a, b) => a.display.localeCompare(b.display, "en", { numeric: true })))
    contigSel.append(h("option", { value: c.display }, `${c.display} (${fmtBp(c.length)})`));
  contigSel.addEventListener("change", () => {
    const r = state.coords.parseRegion(contigSel.value);
    if (r) state.setRegion(r[0], r[1]);
  });

  // --- tracks picker ---
  buildTrackPicker(catalogue.tracks, catalogue.feature_sets, state);

  // --- n-bar ---
  const nbar = $("#nbar");
  function renderNBar(): void {
    nbar.innerHTML = "";
    const st = state.query?.status;
    if (!st?.counts) {
      nbar.append(h("div", { class: "nbar-empty" }, "n — substitutions allowed. Scroll here, use ↑/↓, or alt+scroll anywhere to change n."));
      return;
    }
    const maxV = Math.max(1, ...st.counts.map((c) => c.variants));
    for (const c of st.counts) {
      const frac = Math.log1p(c.variants) / Math.log1p(maxV);
      const chance = Math.min(1, Math.log1p(c.expected_chance_hits) / Math.log1p(maxV));
      const noisy = c.expected_chance_hits > 0.5 * c.variants;
      const cell = h(
        "div",
        {
          class: `ncell${c.n === state.view.n ? " sel" : ""}${noisy ? " noisy" : ""}`,
          title: `n = ${c.n}: ${fmtInt(c.variants)} variants at ${fmtInt(c.sites)} sites\nexpected by chance (random genome): ${fmtExpected(c.expected_chance_hits)}`,
        },
        h("div", { class: "nnum" }, String(c.n)),
        h(
          "div",
          { class: "nbar-track" },
          h("div", { class: "nbar-fill", style: `height:${(frac * 100).toFixed(1)}%;background:${css(mismatchColor(c.n, Math.max(1, st.max_n)))}` }),
          h("div", { class: "nbar-chance", style: `bottom:${(chance * 100).toFixed(1)}%` }),
        ),
        h("div", { class: "ncount" }, fmtCount(c.variants)),
      );
      cell.addEventListener("click", () => setN(c.n));
      nbar.append(cell);
    }
  }
  nbar.addEventListener(
    "wheel",
    (e) => {
      e.preventDefault();
      if (Math.abs(e.deltaY) > 2) setN(state.view.n + (e.deltaY > 0 ? 1 : -1));
    },
    { passive: false },
  );

  // --- side panel ---
  const side = $("#side");
  const qInfo = h("section", { class: "panel" });
  const hoverInfo = h("section", { class: "panel hover" });
  const rowInfo = h("section", { class: "panel" });
  side.append(qInfo, hoverInfo, rowInfo);

  function renderQueryInfo(): void {
    qInfo.innerHTML = "";
    const q = state.query;
    const st = q?.status;
    if (!q || !st?.counts) {
      qInfo.append(
        h("h3", {}, "Getting started"),
        h(
          "p",
          { class: "hint" },
          "Type or paste a 31-mer (or CONTIG:POS), or double-click anywhere on the ruler or tracks to search the 31-mer starting there. Then scroll through n to watch its relatives appear.",
        ),
      );
      return;
    }
    const c = st.counts[state.view.n];
    qInfo.append(
      h("h3", {}, "Query"),
      h("div", { class: "kmer" }, variantHtml(q.kmer, [])),
      h(
        "table",
        { class: "kv" },
        row("n", `${state.view.n} of ${st.max_n}`),
        row("variants", fmtInt(c.variants)),
        row("sites", fmtInt(c.sites)),
        row("chance", fmtExpected(c.expected_chance_hits)),
        row("tree", tree.starLength ? `${(tree.treeLength / Math.max(1, c.variants)).toFixed(2)} / variant (star ${(tree.starLength / Math.max(1, c.variants)).toFixed(2)})` : "–"),
      ),
      h(
        "div",
        { class: "exports" },
        h("a", { href: api.exportUrl(q.id, state.view.n, state.view.order, "tsv") }, "export TSV"),
        h("a", { href: api.exportUrl(q.id, state.view.n, state.view.order, "bed") }, "export BED"),
      ),
    );
  }

  function row(k: string, v: string | Node): HTMLTableRowElement {
    return h("tr", {}, h("th", {}, k), h("td", {}, v));
  }

  function renderHover(): void {
    hoverInfo.innerHTML = "";
    if (!hover && !clade) return;
    const t = h("table", { class: "kv" });
    if (hover) {
      t.append(row("position", state.coords.label(hover.pos)));
      if (hover.row >= 0 && hover.row < state.rows) t.append(row("row", fmtInt(hover.row)));
    }
    hoverInfo.append(t);
    const v = hover?.variant;
    if (v) hoverInfo.append(variantBlock(v));
    if (clade) hoverInfo.append(h("div", { class: "clade" }, h("b", {}, `clade · ${fmtInt(clade.end - clade.start)} variants`), h("div", {}, `shared: ${clade.label}`)));
  }

  function variantBlock(v: Variant): HTMLElement {
    return h(
      "div",
      { class: "variant" },
      variantHtml(v.kmer, v.positions),
      h("div", { class: "meta" }, `${v.mismatches} mismatch${v.mismatches === 1 ? "" : "es"} · ${fmtInt(v.sites)} site${v.sites === 1 ? "" : "s"}`),
    );
  }

  async function renderRow(): Promise<void> {
    rowInfo.innerHTML = "";
    const q = state.query;
    if (selectedRow === null || !q) return;
    const r = selectedRow;
    const v = await raster.variantAt(r);
    if (!v || selectedRow !== r) return;
    const sites = await api.sites(q.id, { ...state.view, start: 0, end: genome.genome_len, row0: r, row1: r + 1 }, 200);
    rowInfo.innerHTML = "";
    const search = h("button", { class: "small" }, "search this variant");
    search.addEventListener("click", () => run({ kmer: v.kmer }));
    rowInfo.append(h("h3", {}, `Row ${fmtInt(r)}`), variantBlock(v), search);
    const list = h("ul", { class: "sites" });
    for (const s of sites.sites) {
      const li = h("li", {}, `${s.display}:${fmtInt(s.start)} ${s.strand}`);
      li.addEventListener("click", () => state.setRegion(s.pos - 2000, s.pos + K + 2000));
      list.append(li);
    }
    rowInfo.append(list);
    if (sites.total > sites.sites.length) rowInfo.append(h("div", { class: "hint" }, `… ${fmtInt(sites.total - sites.sites.length)} more`));
  }

  raster.onHover = (hi) => {
    hover = hi;
    renderHover();
  };
  raster.onSelectRow = (r) => {
    selectedRow = r >= 0 && r < state.rows ? r : null;
    renderRow();
  };
  tree.onHoverClade = (c) => {
    clade = c;
    renderHover();
  };

  // --- footer ---
  const foot = $("#footer");
  function renderFooter(): void {
    const v = state.view;
    regionInput.value = state.coords.regionLabel(v.start, v.end);
    const rows = state.rows ? `rows ${fmtInt(v.row0)}–${fmtInt(v.row1)} of ${fmtInt(state.rows)}` : "";
    const px = state.bpPerPx(raster.el.clientWidth);
    foot.textContent = [
      `${fmtBp(v.end - v.start)} in view (${px < 1 ? (1 / px).toFixed(1) + " px/bp" : fmtBp(px) + "/px"})`,
      rows,
      state.rows ? `${fmtInt(raster.lastSites)} sites drawn in ${fmtMs(raster.lastRenderMs)}` : "",
      "wheel: zoom · shift+wheel: zoom rows · drag: pan · alt+wheel or ↑/↓: n · double-click axis: search there",
    ]
      .filter(Boolean)
      .join("   ·   ");
  }

  // --- keyboard ---
  window.addEventListener("keydown", (e) => {
    if ((e.target as HTMLElement).tagName === "INPUT") return;
    if (e.key === "ArrowUp") setN(state.view.n - 1);
    else if (e.key === "ArrowDown") setN(state.view.n + 1);
    else if (e.key === "+" || e.key === "=") state.zoomRegion(0.5);
    else if (e.key === "-") state.zoomRegion(2);
    else if (e.key === "ArrowLeft") state.setRegion(state.view.start - (state.view.end - state.view.start) * 0.2, state.view.end - (state.view.end - state.view.start) * 0.2, "pan");
    else if (e.key === "ArrowRight") state.setRegion(state.view.start + (state.view.end - state.view.start) * 0.2, state.view.end + (state.view.end - state.view.start) * 0.2, "pan");
    else if (e.key === "0") {
      state.setRegion(0, state.coords.length);
      state.setRows(0, state.rows);
    } else return;
    e.preventDefault();
  });

  // --- URL state ---
  let hashTimer = 0;
  function writeHash(): void {
    clearTimeout(hashTimer);
    hashTimer = window.setTimeout(() => {
      const p = new URLSearchParams();
      if (state.query) {
        p.set("q", state.query.kmer);
        p.set("max", String(state.query.maxN));
        p.set("n", String(state.view.n));
      }
      p.set("order", state.view.order);
      p.set("r", state.coords.regionLabel(state.view.start, state.view.end));
      history.replaceState(null, "", "#" + p.toString());
    }, 300);
  }

  // --- wiring ---
  state.on((why) => {
    ruler.update();
    tracks.update(why);
    hist.update(why);
    tree.update(why);
    if (why !== "rendered" && why !== "tree") raster.update(why);
    if (why === "tree" || why === "n" || why === "order") renderQueryInfo();
    if (why === "order" || why === "n") renderRow();
    renderFooter();
    if (why === "region" || why === "zoom" || why === "pan") writeHash();
  });
  tracks.onLayout = () => state.emit("layout");

  const params = new URLSearchParams(location.hash.slice(1));
  if (params.get("order")) orderSel.value = params.get("order")!;
  state.view.order = orderSel.value as OrderMode;
  if (params.get("max")) maxNSel.value = params.get("max")!;
  const r = params.get("r") && state.coords.parseRegion(params.get("r")!);
  if (r) state.setRegion(r[0], r[1]);
  else state.emit("init");
  renderNBar();
  renderQueryInfo();
  if (params.get("q")) {
    state.view.n = Number(params.get("n") ?? maxNSel.value);
    run({ kmer: params.get("q")! });
  }
}

function setStatus(msg: string, error = false): void {
  const el = $("#status");
  el.textContent = msg;
  el.classList.toggle("error", error);
}

function buildTrackPicker(tracks: TrackInfo[], sets: FeatureSetInfo[], state: State): void {
  const btn = $("#tracks-btn");
  const pop = $("#tracks-pop");
  const section = (title: string, items: { name: string; desc: string }[], selected: string[], set: (v: string[]) => void) => {
    const box = h("div", { class: "pick-section" }, h("h4", {}, title));
    for (const it of items) {
      const cb = h("input", { type: "checkbox" }) as HTMLInputElement;
      cb.checked = selected.includes(it.name);
      cb.addEventListener("change", () => {
        const cur = new Set(title === "Features" ? state.featureSets : state.signalTracks);
        if (cb.checked) cur.add(it.name);
        else cur.delete(it.name);
        set(items.map((i) => i.name).filter((n) => cur.has(n)));
        state.emit("tracks");
      });
      box.append(h("label", { title: it.desc }, cb, " ", it.name, h("span", { class: "desc" }, it.desc)));
    }
    return box;
  };
  pop.append(
    section(
      "Features",
      sets.map((s) => ({ name: s.name, desc: s.description })),
      state.featureSets,
      (v) => (state.featureSets = v),
    ),
    section(
      "Signal",
      tracks.map((t) => ({ name: t.name, desc: t.description })),
      state.signalTracks,
      (v) => (state.signalTracks = v),
    ),
  );
  btn.addEventListener("click", (e) => {
    e.stopPropagation();
    pop.classList.toggle("open");
  });
  document.addEventListener("click", (e) => {
    if (!pop.contains(e.target as Node)) pop.classList.remove("open");
  });
}

main().catch((e) => {
  document.body.innerHTML = `<pre class="fatal">Failed to start: ${String(e)}</pre>`;
});

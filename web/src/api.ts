// Typed bindings for the sentromap HTTP API (/v1). Binary endpoints return typed arrays.

export type OrderMode = "per_n" | "frozen" | "mst";
export type Status = "queued" | "searching" | "ordering" | "ready" | "failed" | "cancelled";

export interface Contig {
  name: string;
  display: string;
  offset: number;
  length: number;
}

export interface Genome {
  name: string;
  k: number;
  genome_len: number;
  distinct_kmers: number;
  sites: number;
  prefix_len: number;
  has_copy_b: boolean;
  max_n_limit: number;
  contigs: Contig[];
}

export interface TrackInfo {
  name: string;
  kind: string;
  description: string;
  source: string;
  p01: number;
  p99: number;
  min: number;
  max: number;
  coverage: number;
}

export interface FeatureSetInfo {
  name: string;
  description: string;
  source: string;
  count: number;
  kinds: string[];
  kind_counts: number[];
}

export interface Feature {
  start: number;
  end: number;
  kind: number;
  strand: string;
  name: string;
}

export type FeaturesResponse =
  | { mode: "features"; kinds: string[]; features: Feature[] }
  | { mode: "density"; kinds: string[]; total: number; density: number[] };

export interface CreatedQuery {
  id: number;
  kmer: string;
  max_n: number;
  order: OrderMode;
  engine: string;
  estimate_ms: number;
  expected_chance_hits: number;
  cached: boolean;
}

export interface CountAtN {
  n: number;
  variants: number;
  sites: number;
  expected_chance_hits: number;
}

export interface QueryStatus {
  id: number;
  status: Status;
  kmer: string;
  max_n: number;
  order: OrderMode;
  engine: string;
  estimate_ms: number;
  cached: boolean;
  elapsed_ms: number;
  timings_ms: Record<string, number>;
  error: string | null;
  counts: CountAtN[] | null;
}

export interface Clade {
  depth: number;
  start: number;
  end: number;
  label: string;
}

export interface TreeResponse {
  n: number;
  order: OrderMode;
  rows: number;
  tree_length: number;
  star_length: number;
  clades: Clade[];
}

export interface Variant {
  row: number;
  kmer: string;
  mismatches: number;
  positions: number[];
  sites: number;
}

export interface VariantsResponse {
  n: number;
  rows: number;
  query: string;
  variants: Variant[];
}

export interface Site {
  pos: number;
  contig: string;
  display: string;
  start: number;
  strand: "+" | "-";
  row: number;
  mismatches: number;
  kmer: string;
}

export interface Raster {
  width: number;
  height: number;
  rows: number;
  row0: number;
  row1: number;
  sitesDrawn: number;
  renderMs: number;
  /** 4 bytes per cell: count (u16 LE), min mismatches (255 = empty), strand bits (1 '+', 2 '-'). */
  cells: Uint8Array;
}

export class ApiError extends Error {
  constructor(
    public status: number,
    message: string,
  ) {
    super(message);
  }
}

const BASE = "/v1";

async function request(path: string, init?: RequestInit): Promise<Response> {
  const r = await fetch(BASE + path, init);
  if (!r.ok) {
    let msg = r.statusText;
    try {
      msg = (await r.json()).error ?? msg;
    } catch {
      /* not JSON */
    }
    throw new ApiError(r.status, msg);
  }
  return r;
}

async function json<T>(path: string, init?: RequestInit): Promise<T> {
  return (await request(path, init)).json() as Promise<T>;
}

function qs(params: Record<string, string | number | undefined>): string {
  const p = new URLSearchParams();
  for (const [k, v] of Object.entries(params)) if (v !== undefined) p.set(k, String(v));
  return p.toString();
}

export interface ViewParams {
  n: number;
  order: OrderMode;
  start: number;
  end: number;
  row0?: number;
  row1?: number;
}

export const api = {
  genome: () => json<Genome>("/genome"),
  tracks: () => json<{ tracks: TrackInfo[]; feature_sets: FeatureSetInfo[] }>("/tracks"),

  async trackSummary(name: string, start: number, end: number, width: number, signal?: AbortSignal): Promise<Float32Array> {
    const r = await request(`/tracks/${encodeURIComponent(name)}?${qs({ start, end, width })}`, { signal });
    return new Float32Array(await r.arrayBuffer());
  },

  features: (name: string, start: number, end: number, width: number, limit: number, kinds?: string, signal?: AbortSignal) =>
    json<FeaturesResponse>(`/features/${encodeURIComponent(name)}?${qs({ start, end, width, limit, kinds })}`, { signal }),

  async sequence(start: number, end: number, signal?: AbortSignal): Promise<string> {
    return (await request(`/sequence?${qs({ start, end })}`, { signal })).text();
  },

  kmerAt: (pos: number) =>
    json<{ pos: number; contig: string; display: string; start: number; kmer: string | null; valid: boolean }>(`/kmer?pos=${pos}`),

  createQuery: (body: { kmer?: string; pos?: number; contig?: string; start?: number; max_n: number; order: OrderMode }) =>
    json<CreatedQuery>("/queries", { method: "POST", headers: { "content-type": "application/json" }, body: JSON.stringify(body) }),

  status: (id: number) => json<QueryStatus>(`/queries/${id}`),
  cancel: (id: number) => json<unknown>(`/queries/${id}`, { method: "DELETE" }),

  async raster(id: number, v: ViewParams, width: number, height: number, signal?: AbortSignal): Promise<Raster> {
    const r = await request(`/queries/${id}/raster?${qs({ ...v, width, height })}`, { signal });
    const h = (k: string) => Number(r.headers.get(k) ?? 0);
    return {
      width: h("x-width"),
      height: h("x-height"),
      rows: h("x-rows"),
      row0: h("x-row0"),
      row1: h("x-row1"),
      sitesDrawn: h("x-sites-drawn"),
      renderMs: h("x-render-ms"),
      cells: new Uint8Array(await r.arrayBuffer()),
    };
  },

  async histogram(id: number, v: ViewParams, width: number, signal?: AbortSignal): Promise<{ stride: number; counts: Uint32Array }> {
    const r = await request(`/queries/${id}/histogram?${qs({ n: v.n, order: v.order, start: v.start, end: v.end, width })}`, { signal });
    return { stride: Number(r.headers.get("x-stride")), counts: new Uint32Array(await r.arrayBuffer()) };
  },

  tree: (id: number, v: ViewParams, minRows: number, signal?: AbortSignal) =>
    json<TreeResponse>(`/queries/${id}/tree?${qs({ n: v.n, order: v.order, row0: v.row0, row1: v.row1, min_rows: minRows, limit: 5000 })}`, {
      signal,
    }),

  variants: (id: number, n: number, order: OrderMode, offset: number, limit: number, signal?: AbortSignal) =>
    json<VariantsResponse>(`/queries/${id}/variants?${qs({ n, order, offset, limit })}`, { signal }),

  sites: (id: number, v: ViewParams, limit: number, signal?: AbortSignal) =>
    json<{ n: number; total: number; sites: Site[] }>(`/queries/${id}/sites?${qs({ ...v, limit })}`, { signal }),

  exportUrl: (id: number, n: number, order: OrderMode, format: "tsv" | "bed") => `${BASE}/queries/${id}/export?${qs({ n, order, format })}`,
};

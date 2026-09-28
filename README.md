# sentromap-rs

A genome browser for the neighbourhood of a sequence. Pick any 31-base k-mer in a reference
genome and see **every place within n substitutions of it, on either strand**, arranged as a
cladogram of variants and drawn against the genome's annotations. Scroll through n and watch a
unique sequence become a repeat family, then distant relatives, then chance look-alikes.

The design, the measurements behind it and the milestones are in [docs/design.md](docs/design.md).

## Quick start

Requires Rust (the toolchain is pinned in `rust-toolchain.toml`).

```sh
scripts/fetch-data.sh                 # genomes + dm6 annotation tracks into data/ (~5 GB)
cargo build --release
B=target/release/sentromap

$B build data/d_mel.fna -o indexes/dmel --name "D. melanogaster r6"     # ~11 s, 3.2 GB
$B annotate indexes/dmel --gff data/d_mel.gff \
    --tracks data/tracks/dm6/tracks.tsv --alias data/tracks/dm6/dm6.chromAlias.txt
$B serve indexes/dmel                  # http://127.0.0.1:8080
```

In the browser: type a 31-mer (or `CONTIG:POS`), or double-click anywhere on the ruler or tracks
to search the 31-mer starting there. Then scroll the n-bar (or press ↑/↓, or alt+scroll).

- wheel: zoom the genome · shift+wheel: zoom rows · drag: pan · double-click raster: zoom in
- click a row: its variant and sites · click a clade bracket: zoom rows to the clade
- the URL hash holds the query, n, ordering and region, so views can be shared

Command-line search works without the server:

```sh
$B search indexes/dmel AGTACGGGACCAGTACGGGACCAGTACGGGA -n 10              # summary per n
$B search indexes/dmel 3R:500000 -n 4 --output sites                      # k-mer from the genome
$B bench search indexes/dmel        # engine comparison over n (checks the engines agree)
$B bench order indexes/dmel         # ordering cost and tree quality
```

## How it works

| Stage | What | Where |
|---|---|---|
| Build | Partitioned external-memory builder: canonical k-mers → sorted copy A (+ prefix table), rotated copy B, positions tether with select, 2-bit sequence | `crates/sm-index` |
| Search | Full scan (memory-bandwidth bound) or the pigeonhole two-copy walk, chosen per query by a calibrated cost model; both strands | `crates/sm-search` |
| Order | Shared-derived-character tree (per-n or frozen ranking), LCP clade intervals, exact MST for small sets | `crates/sm-order` |
| Aggregate | Result sets with sites sorted by position; rasters rendered at screen resolution; per-column mismatch histograms | `crates/sm-agg` |
| Annotate | GFF3 / RepeatMasker / BED features with an interval index; bigWig and GC% tracks as mean/max pyramids | `crates/sm-annot` |
| Serve | axum API polled by the front end; job queue with admission control and cancellation; LRU result cache | `crates/sm-server` |
| Front end | TypeScript + canvas; the server renders rows × genome at pixel resolution, so millions of rows cost the browser one screen-sized image | `web/` |

## Performance (D. melanogaster r6, Ryzen 9 7900X, 24 threads, WSL2)

| | |
|---|---|
| Index build (122 M distinct 31-mers) | 11 s, 3.9 GB peak RSS, 3.2 GB on disk |
| Exact lookup (n = 0–1) | 7 µs |
| Search, n = 8 / 12 / 15 | 2 ms / 16 ms / 27 ms |
| Query ready (search + assemble + order), repeat family, n = 15 (178 k variants) | 57 ms |
| Query ready, unique k-mer, n = 17 (2.6 M variants) | 0.4 s |
| Raster, 2400 × 1400 px over 2.6 M rows | 26 ms (2.8 MB zstd) |

Full tables: [docs/bench](docs/bench).

## Development

```sh
cargo test --release                          # unit, property, oracle and HTTP tests
cargo test --release -- --ignored             # adds the P = 15 (u32 words) layout; writes 8 GB
cargo clippy --all-targets -- -D warnings

cd web && npm install && npm run build        # type-check and bundle into web/dist (committed)
npm run watch                                 # rebuild on change while `sentromap serve` runs
node shot.mjs http://127.0.0.1:8080/ out.png  # headless screenshot (npx playwright install chromium)
```

Every engine is checked against a brute-force oracle on synthetic genomes with contigs, N runs,
soft-masking, direct and inverted repeats, and low-complexity runs; the yeast Ty1 and E. coli
checks from the design run when `data/` is present.

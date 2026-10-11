# Doom roadmap: checks that try to kill the companion-track idea

**The idea under test.** For an annotation track T, align T on the *crease*: the 31-mer at a position p plus
every genome site within n substitutions of it, oriented by strand. Then ask whether the stack shows structure
better than chance. Scanned genome-wide, the scores form **companion tracks** marking sequence-interactive
regions next to T. These come in raw and null-adjusted versions, indexed by query and by site.

**Rules.**
- Each check states a kill or pivot threshold *before* it is run.
- Checks run on samples, and are capped at about 5% of a phase's compute.
- A phase starts only when the previous phase's checks pass, or a pivot is written down here.
- Each run gets one row in the ledger at the bottom: date, commit, result, verdict.

The backbone null is the **circular shift**. Every site of a neighbourhood moves by the same random offset within
its contig, which is equivalent to rotating the track. This keeps the track's autocorrelation, its regional
texture and the spacing between sites, but destroys any true alignment to sequence. Structure is scored as a
**shape** statistic (crease vs flanks) and a **level** statistic (whole window). Regional effects such as
heterochromatin move the level statistic, not the shape one.

## Kill test zero (before any new infrastructure)

| What | Plan |
|---|---|
| Queries | about 200: unique 31-mers, TE copies, k-mers containing a CTCF motif, random k-mers absent from the genome (negative control) |
| Tracks | CTCF, H3K9me3, phastCons124way, Input |
| n | 0, 4, 8, 10. The query's own site, and anything within the window of it, is left out |
| Null | 200 circular shifts per (query, n) |
| Kill if | no class except TE beats the shifted null on the **shape** statistic, at more than 3× the absent-k-mer control's hit rate, for any track and n |

## Phase checks

| # | Phase | Why this might be doomed | Test | Kill / pivot if |
|---|---|---|---|---|
| 5.3 | up front | Someone has already done it | literature search (see ledger) | close precedent: position against it, or stop |
| 0.1 | 0 | The GPU sweep is too slow | microbenchmark: comparisons/s and signature-accumulation bandwidth | < 10¹¹ comparisons/s: sampled sweeps only |
| 0.2 | 0 | Tracks have no crease-scale information | native interval length and autocorrelation length per bigWig | ≥ 500 bp effective for most tracks: kb-scale claims only |
| 0.3 | 0 | Big neighbourhoods only occur where mappability is poor | fraction of non-unique-31-mer bases in the windows of large-neighbourhood sites | > 70%: ChIP tracks adjusted-only, raw mode labelled as possible artefact |
| 0.4 | 0 | Almost all neighbourhoods are chance-level | census of 10 k random positions against the chance expectation | > 95% at chance at every n: claims narrowed to repeats and motif-like k-mers |
| 0.5 | 0 | Input shows the same structure as the treatment | kill test zero on Input | Input structure > ½ of the treatment's: score every track against Input |
| 1.1 | 1 | Self-leakage: the companion track copies the source | companion vs smoothed source, with and without leaving self out | r > 0.7 after leave-self-out |
| 1.2 | 1 | Miscalibrated statistic | QQ plot of 1 k random and shuffled queries | λ > 1.2 |
| 1.3 | 1 | It only measures base composition | regress the score on the query's GC, CpG and dinucleotide content | R² > 0.5: lead with adjusted mode |
| 1.4 | 1 | Hamming distance adds nothing over plain motifs | the same queries scored with their exact 8–12-mers and a HOMER-style motif | baselines ≥ ours in > 80% of hits: no paper |
| 1.5 | 1 | Tandem arrays inflate significance (pseudoreplication) | effective sample size of significant results | > ½ of hits have effective size < 5: down-weight, report as arrays |
| 1.6 | 1 | It can't recover known biology | planted profiles in synthetic genomes; CTCF-motif k-mers vs CTCF ChIP | planted structure missed, or no CTCF signal at any n: stop and debug |
| 1.7 | 1 | Underpowered | sites needed for the CTCF effect vs typical neighbourhood size | typical < ¼ of needed: restrict to high n and repeats |
| 1.8 | 1 | The GPU shortcut (Haar signatures) distorts scores | shortcut vs exact score on 1 k queries | ρ < 0.9: change basis or exact-score candidates |
| 2.1 | 2 | Visual structure isn't real | blind real-vs-shifted stacks shown to a person | ≤ 60% correct on significant stacks: drop the visual claim |
| 2.2 | 2 | The rendering creates the structure | randomised row order, linear intensity | structure vanishes |
| 3.1 | 3 | Discoveries are just false discoveries (FDR) | full pipeline on 2L with shifted tracks | discoveries > 2× the FDR target |
| 3.2 | 3 | It measures noise | companion tracks from rep1 and rep2 of the same mark on 2L | locus Jaccard < 0.3 |
| 3.3 | 3 | Redundant with tracks that already exist | regress on source track, copy number, GC, RepeatMasker coverage | R² > 0.8 |
| 3.4 | 3 | A few TE families explain everything | family composition of the adjusted loci | > 90% from < 10 families: pivot to a TE-subfamily atlas |
| 3.5 | 3 | It doesn't generalise across the genome | structured queries found on chr2 replicate on chr3 sites? | < 50% replicate |
| 3.6 | 3 | Useless for prediction | predicting masked 10% blocks of the source track | no gain over smoothing the source |
| 4.1 | 4 | Synteny is just RepeatMasker | overlap of synteny loci with annotated repeats | > 95% |
| 4.2 | 4 | Sketches don't measure homology | sketch score vs anchor chaining on 500 site pairs | ρ < 0.7 |
| 4.3 | 4 | Structured backgrounds dominate | raw minus adjusted synteny | adjusted mode leaves nothing: report raw only, as a repeat landscape |
| 5.1 | 5 | Only works for a few tracks | fraction of sampled TF tracks with adjusted hits beyond the motif baseline | < 20% |
| 5.2 | 5 | Sequence-to-function models already capture it | gapped-k-mer SVM (gkm-SVM) on one track explains the companion? | most of it: reframe as interpretability and homology visualisation |

## Ledger

### 5.3 Prior art (2026-10-11): PASS, with positioning required

Web literature search covering: Hamming/mismatch k-mer neighbourhoods tied to chromatin signal, repeat-aware
and TE-consensus metaplots, genome-wide k-mer–epigenome scans, sequence-anchored aggregate profiles, and k-mer
motif discovery from ChIP signal. Nothing found that defines the anchor set as the **Hamming neighbourhood of a
long k-mer**, scans every genome position, and outputs a **companion track**. The closest precedents, and how we
differ:

| Precedent | What it does | Difference |
|---|---|---|
| **NoPeak** (Menzel et al., *Bioinformatics* 2021) | profile of ChIP read 5′ ends in a ±250 bp window around every occurrence of each *exact* short k-mer (8-mers; up to 12 tested); scored by profile height; output is a ranked motif list | exact short k-mers, no Hamming neighbourhood, reads only, motif output with no genome track. **Warning: they found shorter k-mers gave more distinctive profiles. That is direct support for doom check 1.4.** |
| **MOCCS2 / POSMO** | k-mer appearance relative to ChIP peak summits; per-k-mer specificity scores | the opposite anchoring (peaks → k-mers); motif output |
| deepTools / EnrichedHeatmap stacks | heatmaps anchored on user-given features | the anchors are features, not sequence neighbourhoods; no scan |
| TE consensus metaplots (e.g. T2T repeat epigenome papers; PAtChER) | chromatin and methylation profiled per repeat copy in consensus coordinates | anchors from RepeatMasker; restricted to TEs. Our TE-class results will partly re-derive these (check 3.4) |
| kmerRRR (preprint 2026) | local vs genome-wide k-mer counts to rescue mappability in repeats | a mappability tool; useful as a covariate for 0.3 |
| EpiAlignment; WashU comparative browser | pairwise sequence + epigenome alignment; cross-species synteny | pairwise or cross-genome, no neighbourhood aggregation |
| gkm-SVM, Enformer-class models | predict tracks from sequence | predictive models; we produce an aligned, inspectable view. Check 5.2 |

Caveat: a web search is not a systematic review. Re-run 5.3 against Google Scholar/bioRxiv before writing.

### KT0 Kill test zero (2026-10-11): KILL as pre-registered; the test itself failed its own controls

`sentromap doom kt0 indexes/dmel` (154 queries, 3.3 s). Per-query results: `docs/doom/kt0.tsv`.

- **Verdict as run: KILL.** No non-TE class beat the shifted null on shape at > 3× the absent-k-mer rate. The
  best was unique/phastCons at n = 10: 62% vs 42% for absent k-mers, only 1.5×.
- **The negative control "passes" for phastCons (42%) and Input (12–16%).** GC explains it: across queries,
  corr(query GC, shape z) = +0.84 for phastCons and +0.75 to +0.81 for Input, for both the absent and unique
  classes. Random k-mers are 49% GC against the genome's 42%. Their chance neighbours are GC-biased, so they sit
  in exon-like, conserved sequence that also maps well. **This confirms check 1.3 already: raw-mode shape
  structure is largely composition. The circular-shift null preserves the track but not composition, so it
  cannot be the only null.**
- **The positive control was broken.** The 13-mer CTCF motif occurs only 4 times in the genome, which is chance
  level, so the class had 4 queries. Two queries containing the full 15-mer `GCCAGCAGGGGGCGC` scored CTCF shape
  z = 8.4 and 6.2 at n = 10, with their own site left out. Two with a truncated motif scored about 0.
  Suggestive, not evidence.
- **TE class.** There is no crease-shape structure for H3K9me3 at this scale; the mark is broad, which the
  multiscale plan addresses. There is phastCons structure (42–55%), but uncontrolled for composition.

### KT0b Kill test zero, revised (pre-registered 2026-10-11, before running)

KT0's verdict stands on record. KT0b fixes the two defects the test exposed in itself, and nothing else.
Forking-paths warning: this is the second attempt. If KT0b fails, there is no KT0c of the same design.

- Same tracks, n values, window, shape/level statistics, 100 shifts, and α = 1/101.
- **Classes:** unique, TE, CTCF motif (13-mer with ≤ 1 mismatch, either strand), **CTCF summit** (top local
  maxima of the CTCF track on major contigs, with the query 31-mer centred on the summit; its own site is left
  out, so only neighbours count as evidence), and absent.
- **Composition-matched controls:** every genomic query gets a base-shuffled twin that is absent from the
  genome, with identical base composition. Each class's baseline is the significance rate of its own twins.
- **Kill if:** no non-TE genomic class beats its own shuffled twins on the shape statistic, for CTCF, H3K9me3
  or phastCons, at > 3× the twins' rate (floored at α) with ≥ 5 significant queries. Input is a diagnostic only:
  a survivor whose Input rate is as high as its target-track rate is flagged as a mappability or composition
  artefact.

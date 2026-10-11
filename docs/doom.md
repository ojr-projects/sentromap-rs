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

### KT0b result (2026-10-11): CONTINUE, narrowly, on one mechanism

`sentromap doom kt0 indexes/dmel --b` (450 queries, 3.7 s; `docs/doom/kt0b.tsv`). Pre-registration committed
first, in 026c832.

**Verdict as registered: CONTINUE.** Unflagged survivors: ctcf_summit/CTCF at n = 8 (5/43, against 0/26 for
the twins), and ctcf_summit/phastCons at n = 4 and n = 8. Each survivor class needed a closer look:

- **Low-n survivors are floor-inflated.** Shuffled twins are absent k-mers with no neighbours at n ≤ 4, so the
  baseline is the α floor. At n = 0–4 the summit queries' "neighbours" are exact repeat copies, and Input is
  high there too. Not credible.
- **phastCons survivors are power-mismatched.** At n = 8 the twins have about a third of the effective sites
  (5 vs 15). The twins match composition but not neighbourhood size. **Design fix for every later check:
  controls must match neighbourhood size as well as composition.**
- **CTCF is the real survivor.** Per query, the class-level Input flag was misleading:
  - n = 10: 6 of 6 significant summit queries show a positive CTCF peak at the crease (z = 3.5 to 7.0).
    Input is flat in 4 of them and only modest in 2 (z = 1.8 and 3.2, against CTCF z = 3.9 and 7.0).
    Composition-matched twins: 0/50, with a comparable number of sites (median 171 vs 292 effective).
    Composition and power both fail to explain it.
  - n = 8: 3 clean peaks (z = 7.3, 15.0, 6.1, with Input z near 0); 1 repeat/mappability artefact (CTCF z = 31.7
    with Input z = 29.7); 1 dip.
  - The effect is a minority one: about 12% of strong CTCF summits. One survivor contains a CTCF-like motif
    (`CTGGCCGCTGGGGGGCGC`).
- **Unchanged from KT0:** phastCons and Input "structure" follows GC. ctcf_motif/phastCons at n = 10 is 88%
  against 78% for its twins: composition, not motif.

**What this means for the plan:**
1. Adjusted mode has to control composition *and* neighbourhood size. The circular shift alone is
   insufficient (check 1.3, confirmed).
2. Input-concordance has to be scored per query, not per class.
3. **The top threat is now check 1.4.** The surviving signal is the kind a short exact motif could carry. If
   the exact 8–12-mer sites shared by these queries reproduce the CTCF peaks equally well, the Hamming
   neighbourhood adds nothing. Run 1.4 (pre-registered, cheap) before investing in the GPU sweep.

### C1.4 Does the Hamming neighbourhood beat plain exact short motifs? (pre-registered 2026-10-11, before running)

The same queries as KT0b (same seed) for the ctcf_summit, ctcf_motif and unique classes, on the **CTCF** track.
Window, statistics, deduplication, the own-site exclusion and 100 circular shifts are all as in KT0b. For each
query, CTCF shape z is computed for these site sets:

| Set | Sites |
|---|---|
| **H8, H10** (ours) | Hamming-neighbourhood sites at n = 8 and n = 10 |
| **E8c, E10c, E12c** (baseline) | exact genome-wide occurrences, either strand, of the query's *central* k-bp window (k = 8, 10, 12), each aligned so the window sits where it does in the query |
| **Ebest** (stronger baseline, report only) | the best z over *every* k-bp window of the query, for each k. Its upward bias from picking the maximum is not corrected |
| **H10-residual** (report only) | H10 sites whose variant does *not* contain the query's exact central 12-mer: is there structure beyond that motif? |

**Kill if:** among ctcf_summit queries where ours is significant (max of H8 and H10 has p ≤ α), the central
baseline (max of E8c, E10c, E12c) reaches a z at least as high in **more than 80%** of them. Each side gets
its best of its own options: 2 for ours, 3 for the baseline.
**Pivot (narrow the claim to "mostly motif") if:** 50–80%.
**Continue if:** ≤ 50%.
**Inconclusive if:** fewer than 5 queries qualify. No re-run of the same design.

HOMER is not installed, and the exact-k-mer baselines stand in for it. A HOMER comparison stays on the list for
Phase 5.

### C1.4 result (2026-10-11): CONTINUE as registered, with a serious caveat

`sentromap doom c14 indexes/dmel` (150 queries, 12.8 s; `docs/doom/c14.tsv`).

- **Verdict as registered: CONTINUE.** The central baseline matched or beat ours in 3 of 8 qualifying summit
  queries (38%). The count is conservative: one qualifier is a dip (z = −3.9), and one is the
  repeat-mappability artefact flagged in KT0b (z = 22.7). Among the 6 genuine positive peaks the baseline wins
  2, both by hairlines (2.54 vs 2.53, 5.77 vs 5.65).
- **The structure is not carried by the central 12-mer.** H10-residual, with all sites containing the query's
  exact central 12-mer removed, keeps essentially the full z in the clean queries (5.89 → 5.89, 7.00 → 7.00,
  4.89 → 4.89, 4.22 → 3.76). Only the artefact collapses (22.7 → 0.7).
- **Caveat: Ebest dominates in raw z** (45–150 against our 4–7). The best exact 8-mer window, scanned
  genome-wide, gives thousands of sites, and z grows with site count. **As a motif detector the Hamming
  neighbourhood loses badly, and a simpler competitor exists: a companion track built from exact short
  k-mers (NoPeak-style).** The registered comparison used the central window, so this does not change the
  verdict, but it makes "structure beyond the best motif" the claim that has to survive.

### C1.4b Does the structure survive removal of the best short motif? (pre-registered 2026-10-11, before running)

The same 150 queries and CTCF track. For each ctcf_summit query where ours is significant (as in C1.4):

- **Exclude** queries whose Input shape p ≤ α at the same n (per-query concordance, the lesson of KT0b).
- **Best motif:** the query's exact k-bp window (k ∈ {8, 10, 12}, any offset) whose genome-wide occurrences give
  the highest CTCF z. This is the strongest motif proxy available, chosen with access to the CTCF track, so it
  is deliberately harsh on us.
- **Best-motif residual:** the H sites (at whichever n gave our max z) whose variant 31-mer does not contain
  that k-mer at any offset, in either orientation. Score its CTCF shape z.
- Also report the per-site effect size (mean shape per effective site) for H and for the best motif. Raw z is
  not comparable across very different site counts.

**Kill if:** in more than 80% of qualifying queries the best-motif residual has z < 2 (the structure is
entirely the best short motif). **Pivot (to "motif plus context") if:** 50–80%. **Continue if:** ≤ 50%.
**Inconclusive if:** fewer than 5 qualify.

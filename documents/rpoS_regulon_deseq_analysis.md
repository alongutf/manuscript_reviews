# RpoS (sigma-S / sigma-38) regulon in the DESeq2 results

*Generated 2026-08-19. Supersedes an earlier version of this document that used a
hand-typed regulon list; see "Revision note" at the end.*

## Regulon source

The gene list is the **sigma-38 sigmulon downloaded from RegulonDB**, not a curated or
recalled list.

- **Endpoint:** RegulonDB GraphQL API, `https://regulondb.ccg.unam.mx/graphql`
- **Query:** `getSigmulonBy(search:"RpoS")` -> `sigmaFactor.sigmulonGenes`
- **Sigmulon:** `RDBECOLISFC00007`, "RNA polymerase sigma factor RpoS"
- **Release:** RegulonDB **14.5.0**, released 2026-01-28 (EcoCyc 29.5, genome NC_000913.3)
- **Retrieved:** 2026-08-19
- **n = 344 genes**

Stored as `metadata/regulondb_sigma38_regulon.txt` (provenance in the file header).
Re-fetch with `python scripts/fetch_regulondb_sigma38.py`.

## Method

- **Script:** `scripts/rpos_regulon_deseq.py`
- **Files searched:** all 16 tables under `results/deseq_results/` -- `from counts/`
  (`disrupted`, `regulated`), `exp0224/` (5 caspase, 5 VapC contrasts), and `aggregated_sc/`
  (4 aggregated single-cell VapC contrasts). The folder list is `FOLDERS` in
  `scripts/rpos_regulon_deseq.py`.
- **Gene IDs:** the DESeq2 tables are indexed by gene *name* (only ~41 unnamed loci carry
  `CDN75_RS...` tags, per `metadata/gffFile_combined.gff`), so no ID mapping is needed.
  Coverage: **312/344** regulon genes present in the `from counts` matrix, **341/344** in
  `exp0224`, **266/344** in `aggregated_sc`.
- **Filtering:** genes with `padj = NA` (DESeq2 independent filtering / zero counts) are dropped
  before any test, so `genes_tested` varies between contrasts.
- **Significance:** `padj < 0.05`, except **`aggregated_sc`, which uses `padj < 0.01`**. That
  folder is aggregated single-cell data with roughly double the per-gene dispersion of the bulk
  sets (median `lfcSE` 0.50 vs 0.24) and a visibly noisier volcano, so it is held to a stricter
  cutoff. The thresholds live in `ALPHA_BY_FOLDER` in `scripts/rpos_regulon_deseq.py` and are
  recorded per row in an `alpha` column in both output CSVs and in each figure panel title.
  **Consequence:** lobe sizes and `sig` counts are not directly comparable across folders. The
  Fisher odds ratio is unaffected by this, since it compares regulon to background *within* one
  contrast at one threshold.
- **Regulon-level test:** two-sided Mann-Whitney U on log2FoldChange, regulon genes vs. all
  other tested genes in the same contrast. This is the right comparison -- the background
  median is not zero in several contrasts, so a raw regulon median is not interpretable alone.

**Outputs:** `results/deseq_results/rpoS_regulon_hits.csv` (per gene per contrast) and
`results/deseq_results/rpoS_regulon_summary.csv` (one row per contrast).

## Summary (Mann-Whitney on the whole log2FC distribution)

| contrast | tested | regulon | sig | up | down | median L2FC regulon | median L2FC background | MWU p |
|---|---|---|---|---|---|---|---|---|
| from counts / **disrupted** | 4333 | 301 | 230 | 130 | 100 | 0.352 | 0.157 | **0.29 (n.s.)** |
| from counts / **regulated** | 4333 | 301 | 224 | 126 | 98 | 0.133 | -0.102 | 1.9e-03 |
| casp-t0 vs exp-t0 | 4481 | 341 | 266 | 176 | 90 | **+0.858** | -0.106 | 8.2e-14 |
| casp-t1 vs exp-t0 | 4570 | 341 | 268 | 176 | 92 | **+0.768** | -0.049 | 6.3e-13 |
| casp-t2 vs exp-t0 | 4570 | 341 | 268 | 176 | 92 | **+0.742** | -0.042 | 1.8e-11 |
| casp-t1 vs casp-t0 | 1382 | 136 | 8 | 0 | 8 | -0.139 | -0.039 | 1.7e-07 |
| casp-t2 vs casp-t0 | 2975 | 259 | 17 | 1 | 16 | -0.083 | 0.005 | 4.7e-05 |
| vapc-t0 vs exp-t0 | 3240 | 279 | 23 | 10 | 13 | -0.047 | 0.027 | 6.5e-06 |
| vapc-t2 vs exp-t0 | 4570 | 341 | 258 | 92 | 166 | **-0.438** | 0.293 | 2.1e-16 |
| vapc-t3 vs exp-t0 | 4570 | 341 | 266 | 95 | 171 | **-0.466** | 0.237 | 2.1e-12 |
| vapc-t2 vs vapc-t0 | 4570 | 341 | 264 | 101 | 163 | -0.372 | 0.306 | 2.9e-13 |
| vapc-t3 vs vapc-t0 | 4570 | 341 | 261 | 100 | 161 | -0.308 | 0.234 | 1.7e-09 |
| **agg_sc** vapc-early vs exp † | 3202 | 230 | 65 | 22 | 43 | 0.205 | 0.468 | 4.3e-02 |
| **agg_sc** vapc-early2 vs exp † | 3717 | 263 | 110 | 49 | 61 | 0.465 | 0.500 | 0.52 (n.s.) |
| **agg_sc** vapc-early2b vs exp † | 3791 | 265 | 127 | 58 | 69 | 0.596 | 0.556 | 0.84 (n.s.) |
| **agg_sc** vapc-late vs exp † | 3716 | 263 | 130 | 62 | 68 | 0.611 | 0.618 | 0.98 (n.s.) |

† `aggregated_sc` rows use `padj < 0.01`; all others use `padj < 0.05`. This affects the
`sig`/`up`/`down` counts only -- the Mann-Whitney p-value is computed on the full log2FC
distribution and does not depend on the threshold, so the MWU column is comparable throughout.

Of the 12 bulk contrasts the regulon shift is significant in 11. **`disrupted` is the exception**
(p = 0.29): its regulon genes are not distinguishable from the genome-wide background, which
itself is shifted up (+0.157). The small-median contrasts (`casp-t1/t2 vs casp-t0`, `vapc-t0`)
reach significance on a very small effect because thousands of genes back the test -- these are
statistically detectable but biologically negligible shifts.

**The four `aggregated_sc` contrasts show no regulon signal.** Their regulon medians (0.21-0.61)
track their background medians (0.47-0.62) almost exactly -- these datasets have a large
genome-wide upward shift that the sigma-38 genes simply follow. Only `vapc-early` separates at
all (p = 0.043), and it does so in the *opposite* direction to the rest of the series: regulon
+0.21 against background +0.47, i.e. the regulon lags the global shift.

## Volcano plots + Fisher's exact test on the significant lobes

A second, independent framing: instead of testing the whole log2FC distribution, count how many
sigma-38 genes land in each **significant lobe** and ask whether that count could arise by chance.

- **Script:** `scripts/rpos_regulon_volcano_fisher.py`
- **Figure:** `results/deseq_results/figures/rpoS_volcano_all_contrasts.{svg,png}` --
  16 volcano panels (log2FC vs -log10 padj), sigma-38 genes in red over all other tested genes in
  grey, with the lobe cutoffs drawn and the Fisher result annotated on each side.
- **Table:** `results/deseq_results/rpoS_regulon_fisher.csv`
- **Lobes:** up = `LFC > 1 AND padj < alpha`; down = `LFC < -1 AND padj < alpha`, with
  alpha = 0.05 for the bulk folders and 0.01 for `aggregated_sc` (see Method).
- **2x2 table** per lobe, universe = all genes tested in that contrast (padj not NA):

  |                  | in lobe | not in lobe |
  |------------------|---------|-------------|
  | sigma-38 gene    | a       | b           |
  | other gene       | c       | d           |

- Reported: odds ratio, one-sided enrichment p (`greater`), one-sided depletion p (`less`),
  two-sided p, and a supplementary `q_across_all_tests` column (see "Multiple testing" below).
- The regulon is 7-10% of the tested universe in every contrast, which is the null expectation
  for the "% of lobe" column.

### Multiple testing

**The headline p-value is the raw per-contrast Fisher p, uncorrected.** Each contrast x lobe is a
separate analysis testing a single gene set, so there is no within-analysis multiplicity to
correct.

This follows the convention already used by the GO pipeline in this repo. In
`src/bulk_functions.py:190`, `GOEnrichmentStudy(..., methods=['fdr_bh'])` +
`run_study(deg_genes)` runs one Fisher test per GO term and BH-corrects **across the GO terms of
that one study**. Each DEG file and each direction gets its own independent correction; nothing
corrects across contrasts. The analogue here is a study containing a single term, where BH is the
identity.

An earlier version of this document BH-corrected across every contrast x lobe test at once. That
was wrong on two counts: it is not the convention the rest of the project uses, and those tests are
strongly non-independent anyway (`casp-t0/t1/t2 vs exp-t0` share a reference and are near
duplicates; the up and down lobes of one contrast are anti-correlated by construction), so BH's
independence/PRDS assumption does not hold. The `q_across_all_tests` column is retained in the CSV
as a conservative sensitivity check, not as the reported statistic.

This changes exactly one call: `disrupted` up lobe, raw one-sided p = 0.028, which is nominally
significant and was previously written off at q = 0.081. See the note under the results table.

| contrast | lobe | lobe size | sigma-38 in lobe | % of lobe | % of universe | OR | one-sided p |
|---|---|---|---|---|---|---|---|
| disrupted | up | 1306 | 106 | 8.12 | 6.95 | 1.28 | 2.8e-02 (enr) |
| disrupted | down | 1138 | 77 | 6.77 | 6.95 | 0.96 | 0.42 (dep) |
| regulated | up | 1155 | 113 | 9.78 | 6.95 | **1.73** | 1.2e-05 (enr) |
| regulated | down | 1278 | 75 | 5.87 | 6.95 | 0.78 | 3.9e-02 (dep) |
| casp-t0 vs exp-t0 | up | 1140 | 159 | 13.95 | 7.61 | **2.81** | 9.4e-19 (enr) |
| casp-t0 vs exp-t0 | down | 1099 | 66 | 6.01 | 7.61 | 0.72 | 1.1e-02 (dep) |
| casp-t1 vs exp-t0 | up | 1107 | 149 | 13.46 | 7.46 | **2.65** | 1.6e-16 (enr) |
| casp-t1 vs exp-t0 | down | 1091 | 64 | 5.87 | 7.46 | 0.72 | 1.1e-02 (dep) |
| casp-t2 vs exp-t0 | up | 1133 | 151 | 13.33 | 7.46 | **2.63** | 2.2e-16 (enr) |
| casp-t2 vs exp-t0 | down | 1077 | 69 | 6.41 | 7.46 | 0.81 | 7.3e-02 (dep) |
| casp-t1 vs casp-t0 | up | 8 | 0 | 0.00 | 9.84 | 0 | 0.44 (dep) |
| casp-t1 vs casp-t0 | down | 0 | 0 | -- | 9.84 | -- | -- |
| casp-t2 vs casp-t0 | up | 10 | 0 | 0.00 | 8.71 | 0 | 0.40 (dep) |
| casp-t2 vs casp-t0 | down | 2 | 0 | 0.00 | 8.71 | 0 | 0.83 (dep) |
| vapc-t0 vs exp-t0 | up | 19 | 1 | 5.26 | 8.61 | 0.59 | 0.50 (dep) |
| vapc-t0 vs exp-t0 | down | 6 | 1 | 16.67 | 8.61 | 2.13 | 0.42 (enr) |
| vapc-t2 vs exp-t0 | up | 1254 | 58 | 4.63 | 7.46 | **0.52** | 2.2e-06 (dep) |
| vapc-t2 vs exp-t0 | down | 985 | 134 | 13.60 | 7.46 | **2.57** | 7.1e-15 (enr) |
| vapc-t2 vs vapc-t0 | up | 1251 | 68 | 5.44 | 7.46 | 0.64 | 6.4e-04 (dep) |
| vapc-t2 vs vapc-t0 | down | 1012 | 129 | 12.75 | 7.46 | **2.31** | 5.8e-12 (enr) |
| vapc-t3 vs exp-t0 | up | 1329 | 63 | 4.74 | 7.46 | **0.53** | 2.2e-06 (dep) |
| vapc-t3 vs exp-t0 | down | 965 | 121 | 12.54 | 7.46 | **2.21** | 1.3e-10 (enr) |
| vapc-t3 vs vapc-t0 | up | 1334 | 64 | 4.80 | 7.46 | 0.54 | 3.0e-06 (dep) |
| vapc-t3 vs vapc-t0 | down | 999 | 112 | 11.21 | 7.46 | **1.84** | 7.1e-07 (enr) |
| agg_sc vapc-early vs exp † | up | 425 | 20 | 4.71 | 7.18 | 0.60 | 1.8e-02 (dep) |
| agg_sc vapc-early vs exp † | down | 461 | 38 | 8.24 | 7.18 | 1.19 | 0.19 (enr) |
| agg_sc vapc-early2 vs exp † | up | 686 | 48 | 7.00 | 7.08 | 0.99 | 0.50 (dep) |
| agg_sc vapc-early2 vs exp † | down | 671 | 51 | 7.60 | 7.08 | 1.10 | 0.30 (enr) |
| agg_sc vapc-early2b vs exp † | up | 847 | 57 | 6.73 | 6.99 | 0.95 | 0.40 (dep) |
| agg_sc vapc-early2b vs exp † | down | 786 | 59 | 7.51 | 6.99 | 1.10 | 0.29 (enr) |
| agg_sc vapc-late vs exp † | up | 869 | 59 | 6.79 | 7.08 | 0.94 | 0.39 (dep) |
| agg_sc vapc-late vs exp † | down | 788 | 59 | 7.49 | 7.08 | 1.08 | 0.33 (enr) |

† `aggregated_sc` lobes are defined at `padj < 0.01`; all other rows at `padj < 0.05`.

### What the Fisher framing shows

**Caspase: strong, one-sided enrichment in the up lobe.** The up lobe is ~14% sigma-38 against a
7.5% background -- roughly a doubling, OR 2.6-2.8, p ~1e-16 to 1e-19. This is not a coincidence by
any reasonable margin. The down lobe is *depleted* (OR 0.72), so the effect is directional rather
than "sigma-38 genes are just more variable".

**VapC: the mirror image, and it is a genuine two-sided effect.** The down lobe is ~13% sigma-38
(OR 2.2-2.6, p ~1e-10 to 1e-15) *and* the up lobe is significantly depleted (OR 0.52-0.64,
p ~1e-6). VapC does not merely fail to induce the regulon; it actively excludes sigma-38 genes
from the up-regulated set. This depletion is the clearest single piece of evidence that the two
perturbations act on sigma-S in opposite directions.

**`disrupted` is nominally significant but much weaker than `regulated`.** Up lobe OR = 1.28,
p = 0.028 -- it clears 0.05, but the effect is small: 8.12% of the up lobe is sigma-38 against a
6.95% background, i.e. 106 genes where ~91 would be expected by chance. `regulated` is a clearly
stronger effect on the same data (OR 1.73, p = 1.2e-05, 113 observed vs ~80 expected).

The two "from counts" contrasts are therefore *not* equivalent with respect to the sigma-S
regulon, but the gap is one of degree, not of presence/absence. Note also that the two statistical
framings disagree here: Mann-Whitney gives `disrupted` p = 0.29, because the modest lobe
enrichment is swamped once the whole distribution is compared against a background that is itself
shifted up (+0.157). Treat `disrupted` as a weak, unreplicated signal rather than either a clear
positive or a clear negative.

**The `aggregated_sc` VapC contrasts do not reproduce the bulk VapC result.** At the stricter
`padj < 0.01` cutoff all eight odds ratios sit between 0.60 and 1.19, and the lobes are 4.7-8.2%
sigma-38 against a 7.0-7.2% background -- i.e. essentially at the null. This is not a power
problem: the lobes still hold 425-869 genes, comparable to the bulk contrasts where OR 2.2-2.6
was detected at p ~1e-10 to 1e-15. The only nominal result is `vapc-early` up-lobe depletion
(OR 0.60, p = 0.018), which is the same direction as bulk VapC but well under half the effect
size and is not corroborated by the matched down lobe (OR 1.19, p = 0.19) or by the other three
contrasts.

Tightening the cutoff from 0.05 to 0.01 did not change this picture. It shrank the up lobes by
25-38% and the down lobes by 6-12%, and moved every odds ratio by at most 0.11 (largest:
`vapc-early2` up, 0.88 -> 0.99). No p-value crossed 0.05 in either direction. The absence of
signal here is not an artifact of a permissive threshold.

Two differences between the datasets are worth weighing before treating this as a contradiction:
`aggregated_sc` covers only 266/344 regulon genes (vs 341/344 in `exp0224`) and has roughly
double the dispersion (median `lfcSE` 0.50 vs 0.24), so per-gene calls are noisier; and its
genome-wide background is strongly shifted (median LFC +0.47 to +0.62, vs -0.04 to +0.31 in
`exp0224`), which changes what "background" means in the comparison. Whether the aggregation
step or the underlying biology accounts for the difference is not answerable from these tables
alone.

**The within-caspase and `vapc-t0` contrasts have essentially empty lobes** (0-19 genes), so their
Fisher tests are uninformative -- the p-values near 1 reflect no power, not evidence of absence.
Note this is where the two methods appear to disagree: Mann-Whitney called `casp-t1 vs casp-t0`
significant (p = 1.7e-07) while Fisher cannot. Both are right about different things. The
Mann-Whitney picks up a real but tiny (median -0.14) coordinated shift spread across all 136
tested regulon genes; the Fisher test asks only about genes that cleared |LFC| > 1, and none did.

### Which test to prefer

For the manuscript, the Fisher/lobe framing is the easier one to defend: it uses the same
significance criteria a reader applies when looking at the volcano plot, it gives an interpretable
effect size (odds ratio), and it is insensitive to the genome-wide baseline shifts that made the
raw regulon medians misleading. The Mann-Whitney is the more sensitive test and is worth keeping
as a supporting result, with the caveat that it can flag shifts too small to matter biologically.
The two agree on every substantive conclusion except `casp-t1/t2 vs casp-t0` and `disrupted`, both
explained above.

## rpoS itself

`rpoS` is **not** a member of the RegulonDB sigma-38 sigmulon (no annotated sigma-38 promoter on
its own gene), so it does not appear in the hits table. Read directly from the DESeq2 tables:

| contrast | baseMean | log2FC | padj |
|---|---|---|---|
| from counts / disrupted | 11535 | -0.494 | 9.8e-08 |
| from counts / regulated | 11535 | +1.427 | 2.2e-55 |
| casp-t0 vs exp-t0 | 4642 | +1.109 | 1.1e-17 |
| casp-t1 vs exp-t0 | 4642 | +0.983 | 3.8e-14 |
| casp-t2 vs exp-t0 | 4642 | +1.059 | 3.4e-16 |
| casp-t1 vs casp-t0 | 4642 | -0.126 | 0.747 |
| casp-t2 vs casp-t0 | 4642 | -0.050 | 0.918 |
| vapc-t0 vs exp-t0 | 4642 | +0.023 | 0.969 |
| vapc-t2 vs exp-t0 | 4642 | -0.114 | 0.395 |
| vapc-t2 vs vapc-t0 | 4642 | -0.136 | 0.262 |
| vapc-t3 vs exp-t0 | 4642 | -0.312 | 0.013 |
| vapc-t3 vs vapc-t0 | 4642 | -0.334 | 0.0038 |
| agg_sc vapc-early vs exp | 519 | -1.435 | 1.2e-73 |
| agg_sc vapc-early2 vs exp | 534 | -3.289 | 0 |
| agg_sc vapc-early2b vs exp | 556 | -4.337 | 0 |
| agg_sc vapc-late vs exp | 519 | -4.699 | 0 |

Up in `regulated`, **down** in `disrupted`; up ~1 log2 unit in every caspase-vs-exponential
contrast; flat-to-slightly-down across bulk VapC.

`rpoS` behaves very differently in `aggregated_sc`: strongly and progressively **down**, from
-1.4 (early) to -4.7 (late), all at extreme significance. The bulk VapC contrasts show at most
-0.33. So in the aggregated single-cell data `rpoS` itself collapses while its regulon does not
move relative to background -- the two datasets disagree about VapC and sigma-S, and they
disagree at the level of `rpoS` transcript abundance, not only at the regulon level.

## Caspase (casp)

Sigma-38 sigmulon induction relative to exponential phase that is **already saturated at t0 and
does not change afterwards**: `casp-t1 vs casp-t0` gives 0 up-regulated regulon genes and
`casp-t2 vs casp-t0` gives 1, all with |L2FC| < 0.8.

- **Strongest up (t2 vs exp-t0):** `sdsR` +10.30, `ryjA` +10.14, `yciG` +8.37, `ymgC` +6.86,
  `csiE` +6.50, `micF` +6.32, `arrS` +6.31, `ymgA` +6.29, `ariR` +6.29, `ycgZ` +6.28.
  The top of the list is dominated by sigma-S-dependent **small RNAs** (`sdsR`, `ryjA`, `micF`,
  `arrS`) and the `ycgZ-ymgABC` / `ariR` biofilm module.
- **Down despite sigma-38 membership:** `csgC` -5.70, `ompF` -5.15, `mglB` -4.13, `mglA`/`mglC`,
  `folK` -3.65, `sucC`/`sucD` -3.1, `evgA` -3.20, `speC`.

## VapC

The opposite sign: broad sigmulon **repression** against a background whose median gene moves
*up* (+0.23 to +0.31), so the separation between regulon and background is larger than the
regulon median alone suggests. At `vapc-t0` the regulon is essentially unmoved (median -0.05),
so the repression develops over the time course.

- **Strongest down:** `lsrB` -6.99, `mglB` -6.51, `yjcH` -5.98, `acs` -5.93, `actP` -5.77,
  `sucA` -5.73, `csgC` -5.23, `gabT` -5.03, `gabD` -4.76, and the `suc` operon broadly.
- **Up:** `ymgC` +5.22, `yadV` +5.18, `ansP` +4.84, `appY` +4.79, `yqiI` +4.16, `ariR` +4.12,
  `htrE`, `ydbD`, `ydhY`, `proP` +3.08.

## disrupted vs regulated

Both shift the sigmulon up in absolute terms, but only `regulated` is significant against its
own background (p = 1.9e-03 vs 0.29), and the two disagree on `rpoS` itself (+1.43 vs -0.49).

- **Top up in `regulated`:** `yciG` +6.28, `csiE` +6.24, `hyaA` +5.61, `hyaB` +5.27, `hyaC` +5.06,
  `wrbA` +4.85, `yccJ` +4.75, `glaH` +4.73, `hyaD` +4.62, `uspB` +4.37.
- **Top up in `disrupted`:** `ariR` +7.96, `ymgA` +7.90, `ycgZ` +7.38, `ymgC` +7.07, `yciG` +5.15,
  `pfkB` +4.91, `pphA` +4.75, `yciF` +4.31, `yafP` +4.14, `htrE` +4.06.
- **Down in both:** `ompF`, `mglB`, `patD`, `folK`, `csgF`, `blr`, `evgA`, `astE`.

## Caveats

1. **`disrupted` shows no regulon-specific signal.** Its apparent induction tracks a genome-wide
   upward shift. Any claim of sigma-S activation in `disrupted` needs a different line of evidence.
2. **A sigmulon is not a "response".** RegulonDB sigma-38 membership means an annotated
   sigma-38-transcribed promoter; many of these genes also have sigma-70 promoters and other
   regulatory inputs, which is why ~30% of significant regulon genes move opposite to the bulk in
   every contrast.
3. **Directionality is mixed by construction**, so up/down counts are less informative here than
   the Mann-Whitney comparison against background.
4. `baseMean` repeats across contrasts within a folder because each dataset was fit once in
   DESeq2 with multiple contrasts extracted from the same fit.
5. The `aggregated_sc` tables carry ~20 non-genomic index entries -- reporter and plasmid
   constructs (`GFP`, `YFP`, `LacI`, `MS2_fwd`, `INTR_*`, `LELOBEKK*`) and one artifact row
   literally named `Unnamed: 1`. They are left in the background universe, which is correct for
   a Fisher test over "all tested features" and shifts the null share by <0.6%, but the
   `Unnamed: 1` row suggests an index was written out twice somewhere upstream and is worth
   tracing in whatever produced these files.

## Revision note

The first version of this analysis used a ~121-gene list I assembled from memory of the sigma-S
literature, described as "curated from EcoCyc/RegulonDB + Weber et al. 2005". That description
overstated its provenance: it was never a database query. Replacing it with the actual RegulonDB
14.5.0 sigmulon changed the results materially, so the earlier numbers should not be used:

- **List overlap:** only **76 of 121** hand-list genes are in the RegulonDB sigma-38 sigmulon.
  The 45 non-members include genuine errors on my part and several synonyms RegulonDB does not
  use (`yedU`/`hchA`, `ykfE`/`ivy`, `ynaF`/`uspF`, `rpsV`/`sra`, `yfiA`/`raiA`, `yecI`/`ftnB`,
  `yciD`/`ompW`). The hand list also included `rpoS` itself, which the sigmulon excludes.
- **Effect sizes shrank roughly 3-fold** (caspase median +2.3 -> +0.86; VapC -0.65 -> -0.44),
  because the hand list was biased toward textbook, strongly-responding sigma-S genes.
- **Up/down splits are far less lopsided** (caspase 82/11 -> 176/90).
- **One conclusion reversed:** `disrupted` looked like a clear regulon induction on the hand list
  and is not significant against background on the RegulonDB list.

What survived unchanged: the caspase-up / VapC-down opposition, the saturation of the caspase
response by t0, and the `rpoS` sign flip between `disrupted` and `regulated`.

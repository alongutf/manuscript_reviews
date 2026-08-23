# Figure S11 — Implementation Details

Supplementary figure showing **volcano plots for three DESeq2 contrasts with the RpoS
(sigma-38) regulon highlighted**, plus a Fisher's exact test of regulon over/under-representation
in the significant lobes annotated on each panel.

**Layout (1×3):** A Dis-Arrest · B Reg-Arrest · C Early VapC.

Panel titles are the short experiment names only — the `padj` cutoff is deliberately **not** printed
in the titles. It is `padj < 0.05` for all three panels (see "Significance threshold" below).

The statistics are identical to those in `scripts/rpos_regulon_volcano_fisher.py` /
`documents/rpoS_regulon_deseq_analysis.md` — this figure is a three-panel, publication-styled
extract of that 16-panel diagnostic figure, not a new analysis.

## How to run

The script lives in `scripts/supplementary_figures/`.

```bash
cd scripts/supplementary_figures
python figure_s11.py       # writes figure_s11.pdf + figure_s11_preview.png next to the script
```

Output paths are absolute (built from `_HERE`), so the script can also be run from the repo root.

Import bootstrap (path-independent), same as the other S-figures plus one extra entry so the
shared regulon/contrast helpers can be imported:

```python
_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))                # repo root
sys.path.insert(0, _REPO)                                      # -> import src.*
sys.path.insert(0, os.path.join(_REPO, 'scripts', 'figures'))  # -> figure_functions
sys.path.insert(0, os.path.join(_REPO, 'scripts'))             # -> rpos_regulon_deseq
```

---

## Inputs

| Input | Path | Used for |
|---|---|---|
| Regulon gene list | `metadata/regulondb_sigma38_regulon.txt` | Purple points in all panels. RegulonDB 14.5.0 sigma-38 sigmulon `RDBECOLISFC00007`, n = 344 genes. Loaded via `load_regulon()`; the file header also supplies the release string printed in the legend. Re-fetch with `scripts/fetch_regulondb_sigma38.py`. |
| DESeq2 contrast (A) | `results/deseq_results/from counts/deseq2_results_disrupted.csv` | Panel A |
| DESeq2 contrast (B) | `results/deseq_results/from counts/deseq2_results_regulated.csv` | Panel B |
| DESeq2 contrast (C) | `results/deseq_results/aggregated_sc/deseq2_results_vapc-early_vs_exp.csv` | Panel C |

Reused from `scripts/rpos_regulon_deseq.py`: `DESEQ` (results dir), `REGULON` (list path),
`load_regulon()`. The regulon list and its provenance are therefore never redefined locally. The
significance threshold **is** set locally (`ALPHA`), deliberately overriding that module's
`alpha_for()` / `ALPHA_BY_FOLDER` — see below.

---

## Global settings

| Variable | Value | Meaning |
|---|---|---|
| `fsize` | `10` | base font size (panel letters 12pt bold via `PanelFigure`; titles 10pt; axis labels 8pt; ticks 7pt; annotation boxes 6pt; legend 8pt) |
| `figsize` | `(10.5, 4.0)` in | single-row, three-column landscape layout |
| `ALPHA` | `0.05` | significance threshold, **the same for all three panels** |
| `LFC_CUT` | `1.0` | lobe boundary, `|log2FoldChange| > 1` |
| `SIG_COLOR` | `#5BC8DC` (light cyan) | significant genes |
| `NS_COLOR` | `#c8c8c8` (grey) | non-significant genes |
| `REG_COLOR` | `#7B3FA0` (purple) | **significant** rpoS (sigma-38) regulon genes |

Panel rects (figure-normalized `[l, b, w, h]`):
`A [0.070, 0.165, 0.245, 0.585]`, `B [0.395, 0.165, 0.245, 0.585]`, `C [0.720, 0.165, 0.245, 0.585]`,
with `label_offset=(-0.045, 0.015)`. Legend is a figure-level legend at `(0.5, 0.005)`, 3 columns,
no frame.

---

## Point classification (the three colors)

Genes are drawn in three passes so the ordering is deterministic. All three are semi-transparent so
the dense core of the volcano does not read as a solid block:

1. **grey** — `~significant`, `s=4`, `alpha=0.45`, rasterized
2. **light cyan** — `significant & ~in_regulon`, `s=4`, `alpha=0.55`, rasterized
3. **purple** — `significant & in_regulon`, `s=9`, `alpha=0.75`, `zorder=3`

**Significance takes precedence over regulon membership**: a sigma-38 gene that does not clear both
cutoffs is grey, exactly like any other non-significant gene. Purple therefore marks *significant
regulon genes only* — which is also what the Fisher test counts, so the purple points in each lobe
are precisely the `a` cell of that lobe's 2×2 table.

**"Significant" means lobe membership**, `|log2FC| > LFC_CUT` **and** `padj < alpha` — i.e. the same
criterion the Fisher test uses, so the cyan cloud is exactly the union of the two lobes being
tested. Genes that pass `padj` but not the fold-change cut are grey.

Only genes with `padj` and `log2FoldChange` both non-NA are plotted or tested (DESeq2 independent
filtering drops the rest), which is why `genes_tested` differs between panels.

### Significance threshold

`ALPHA = 0.05` for **all three panels**, so lobe sizes, regulon shares and odds ratios are directly
comparable across them and one dashed `-log10(alpha)` line means the same thing in each.

This deliberately overrides `ALPHA_BY_FOLDER` in `scripts/rpos_regulon_deseq.py`, which holds the
`aggregated_sc` folder (panel C, Early VapC) to `padj < 0.01` because that data has roughly double
the per-gene dispersion of the bulk sets (median `lfcSE` 0.50 vs 0.24). That stricter cutoff still
governs the 16-panel diagnostic figure from `rpos_regulon_volcano_fisher.py` and the tables in
`documents/rpoS_regulon_deseq_analysis.md`, which are **not** changed by this figure — so the
Early VapC row there (OR 0.60) and panel C here (OR 0.67) are the same test at two thresholds, not
a discrepancy.

Effect of the change on panel C, the only panel affected:

| lobe | alpha | lobe size | regulon in lobe | % of lobe | OR | one-sided p |
|---|---|---|---|---|---|---|
| up | 0.01 | 425 | 20 | 4.71 | 0.60 | 0.018 (dep) |
| up | **0.05** | **680** | **36** | **5.29** | **0.67** | **0.017 (dep)** |
| down | 0.01 | 461 | 38 | 8.24 | 1.19 | 0.19 (enr) |
| down | **0.05** | **525** | **44** | **8.38** | **1.23** | **0.14 (enr)** |

The up lobe grows 60% and the down lobe 14%, the odds ratios move by <0.07, and the call is
unchanged in both direction and significance: the regulon is depleted from the up-regulated genes
(p ~0.017 either way) and not significantly enriched in the down-regulated ones. The panel's
conclusion does not depend on the threshold. The regulon's share of the tested universe (7.18%,
230/3202 genes) is unaffected — the universe is every gene with a non-NA `padj`, which no
threshold changes.

### y-axis flooring

A few genes have `padj` exactly 0. Rather than let `-log10` hit the float floor (~300), padj is
floored at the smallest **nonzero** padj in that same contrast. This keeps the y-axis honest per
panel; it is the same treatment used in `rpos_regulon_volcano_fisher.py`. Visible as a flat row of
points at the top of panels A and C.

---

## Statistics — Fisher's exact test (`_fisher`)

Per lobe, universe = all genes tested in that contrast:

|               | in lobe | not in lobe |
|---------------|---------|-------------|
| sigma-38 gene | a       | b           |
| other gene    | c       | d           |

`scipy.stats.fisher_exact` is called three times (two-sided for the odds ratio, `greater` for
enrichment, `less` for depletion).

**No multiple-testing correction is applied.** Each contrast × lobe is a separate analysis testing
a single gene set — the same convention the GO pipeline uses
(`src/bulk_functions.py:190`, where goatools BH-corrects across the GO terms of *one* study and
never across studies). Rationale and the sensitivity check are in
`documents/rpoS_regulon_deseq_analysis.md`.

---

## The annotation box (`_lobe_call` / `_best_lobe`)

**One box per panel**, on the lobe with the stronger result only. `_lobe_call` returns the
one-sided p-value and the text for a lobe; `_best_lobe` picks whichever of the two has the smaller
p and reports which lobe it was. The box is then placed over that side of the plot — up lobe at the
top right (`x=0.97`, `ha='right'`), down lobe at the top left (`x=0.03`, `ha='left'`) — at
`y=0.985`, `va='top'`, `multialignment='left'`. White face, **black** rounded border
(`boxstyle='round,pad=0.35'`, `linewidth=0.7`, `alpha=0.92`). `ax.set_ylim(top=...*1.42)` adds the
headroom that keeps the box clear of the points.

**In all three panels the up lobe wins**, so every box sits at the top right: A and B report
enrichment in the up-regulated genes, C reports depletion from them. The selection is computed, not
hard-coded — if the underlying tables change, the box follows the stronger lobe and moves side
accordingly.

Each box has three lines:

```
rpoS regulon <enriched|depleted>
OR = x.xx
p = …
```

The word and the p-value follow the direction the odds ratio points, so the two are never
inconsistent:

- `OR >= 1` → `enriched`, one-sided `p_greater`
- `OR < 1`  → `depleted`, one-sided `p_less`

This matters most in panel C, whose up lobe is *depleted* of regulon genes (OR 0.67) rather than
enriched — quoting an enrichment p there would misstate the test. A lobe with zero genes would
print `"rpoS regulon / not testable"` and score `p = inf` so the other lobe is chosen; that does
not arise for these three contrasts, where all lobes hold 525–1306 genes.

The unannotated lobe is still computed and printed to the console — it is simply not drawn, since
in each panel it is the weaker and (in A and C) non-significant direction. Its values are in the
table below.

`_fmt_p` formats to **two decimals** (`p = 0.42`) down to 0.01, and switches to one decimal with an
explicit decade below that (`p = 1.2×10⁻⁵`). OR is always two decimals.

---

## Values in the current figure

All panels at `alpha = 0.05`. A and B match `results/deseq_results/rpoS_regulon_fisher.csv` exactly;
C is recomputed here at 0.05 and so differs from that CSV's 0.01 rows (see "Significance threshold").

| Panel | contrast | alpha | tested | regulon tested | lobe | lobe size | regulon in lobe | OR | one-sided p |
|---|---|---|---|---|---|---|---|---|---|
| A | disrupted | 0.05 | 4333 | 301 | up | 1306 | 106 | 1.28 | 0.028 (enr) |
| A | disrupted | 0.05 | 4333 | 301 | down | 1138 | 77 | 0.96 | 0.42 (dep) |
| B | regulated | 0.05 | 4333 | 301 | up | 1155 | 113 | **1.72** | 1.2e-05 (enr) |
| B | regulated | 0.05 | 4333 | 301 | down | 1278 | 75 | 0.78 | 0.039 (dep) |
| C | agg_sc vapc-early | 0.05 | 3202 | 230 | up | 680 | 36 | **0.67** | 0.017 (dep) |
| C | agg_sc vapc-early | 0.05 | 3202 | 230 | down | 525 | 44 | 1.23 | 0.14 (enr) |

As printed in the box (up lobe in all three panels, top right):

| Panel | box text |
|---|---|
| A Dis-Arrest | rpoS regulon enriched / OR = 1.28 / p = 0.03 |
| B Reg-Arrest | rpoS regulon enriched / OR = 1.72 / p = 1.2×10⁻⁵ |
| C Early VapC | rpoS regulon depleted / OR = 0.67 / p = 0.02 |

Rounding to two decimals compresses the weaker calls (A 0.0283 → 0.03, C 0.0171 → 0.02); the
unrounded values are in the table above and in the console output.

Interpretation caveats — in particular that `disrupted` is a weak, unreplicated signal (Mann-Whitney
on the same contrast gives p = 0.29) and that the `aggregated_sc` VapC contrasts do not reproduce
the bulk VapC regulon repression — are in `documents/rpoS_regulon_deseq_analysis.md`.

---

## Console output

The script prints, per panel, the contrast path, alpha, `genes_tested`, `regulon_tested`, and for
each lobe the lobe size, regulon count, odds ratio, and both one-sided p-values — enough to check
the figure without reopening the PDF. Note the printed `alpha` is 0.05 on every line, including
Early VapC — that is the check that this figure is not using the pipeline's 0.01 default.

---

## Outputs

| File | Notes |
|---|---|
| `scripts/supplementary_figures/figure_s11.pdf` | 300 dpi, vector except the rasterized grey/cyan scatter |
| `scripts/supplementary_figures/figure_s11_preview.png` | 200 dpi preview |

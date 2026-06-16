# Figure 4 — Implementation Details

## How to run

Run from `scripts/figures/` so that `os.getcwd()` resolves correctly:

```bash
cd scripts/figures
python figure4.py          # writes figure4.svg and figure4_preview.png
```

`figure4.py` loads the GO ontology (`metadata/go-basic.obo`) and associations
(`metadata/ecocyc.gaf`) via goatools for panel A, so it prints goatools load
messages and takes a few seconds.

---

## Global settings

| Variable | Value | Where to change |
|---|---|---|
| `fsize` | `10` | assembly block — base font size |
| `figsize` | `(7, 6.5)` inches | `PanelFigure(figsize=...)` in assembly block |
| `label_offset` | `(0, 0.03)` | `PanelFigure(label_offset=...)` |
| `root_dir` | 2 dirs above `os.getcwd()` | resolves to repo root when run from `scripts/figures/` |

The figure is saved both as vector `figure4.svg` (dpi=300, transparent) and as a
raster preview `figure4_preview.png` (dpi=300).

---

## Layout map

Panel positions are `[left, bottom, width, height]` in **figure-normalized coordinates** (0–1).
The figure is 7 × 6.5 inches.

```
 0.0                                                          1.0
 ┌──────────────────────────────────────────────────────────────┐ ~0.93
 │ A (0.08,0.30) [0.10×0.63]   B (0.35,0.66) [0.60×0.27]        │
 │ chemotaxis heatmap          GO-term enrichment bars + inset   │
 │ (genes × condition)                                          │
 │                            ↳ long GO-term-name x labels       │
 │                              occupy the gap below panel B     │
 ├────────────────────────────────────────────────────────────── ~0.27
 │ C (0.30,0.07) [0.30×0.20]   D (0.70,0.07) [0.28×0.20]        │
 │ Relative enrichment vs time Lag-time survival (time in SHX)   │
 └──────────────────────────────────────────────────────────────┘ 0.07
```

Panel B sits high (bottom = 0.66) specifically so the 45°-rotated GO-term-name
x-axis labels have room to drop into the gap above panels C/D without colliding.

---

## Panel A — Chemotaxis gene heatmap

**Function:** `panel_A(ax)`

| Detail | Value |
|---|---|
| DE data | `results/deseq_results/from counts/deseq2_results_disrupted.csv` (`log2FoldChange_x`) and `..._regulated.csv` (`log2FoldChange_y`) |
| Genes shown | genes annotated to GO terms listed in `scripts/figures/figure4/GO_terms_heatmap.csv` (GO:0000027 ribosomal large-subunit assembly + GO:0006935 chemotaxis) |
| GO-ID → name map | `scripts/figures/figure4/GO_ID_name.csv` |
| Plot | `sns.heatmap`, coolwarm centered at 0, two columns: Dis-Arrest, Reg-Arrest |
| Side bracket | drawn per GO-term group with the term name (e.g. "chemotaxis") |

---

## Panel B — GO-term enrichment bars + inset

**Function:** `panel_B(ax)`

**Data files** (goatools over-representation output, `down`-regulated genes):

- `results/GO_results/from_counts/GOATOOLS_GO_enrichment_results_disrupted_down.csv`
- `results/GO_results/from_counts/GOATOOLS_GO_enrichment_results_regulated_down.csv`

| Detail | Value |
|---|---|
| Terms plotted | GO terms significant (in file) in **both** conditions (`common_go_terms`), sorted by Reg-Arrest FDR |
| Bar height | `-log10(FDR)`, labeled "Enrichment score" (this is an ORA −log10(FDR), **not** a GSEA enrichment score) |
| Bars | Dis-Arrest `#de2d26`, Reg-Arrest `#9ecae1`, `bar_width=0.35` |
| **x-axis labels** | **GO term name + gene count**, e.g. `chemotaxis (n=34)` — built from the `Term` column and the first element of `Ratio in Population`; `rotation=45`, `ha='right'`, fontsize `fsize-5` |
| x-axis title | `GO term` |
| Significance `*` | per-term `p_adj` from `scripts/figures/figure4/GO_diff_pvals.csv` (looked up by GO_ID), brackets drawn above the taller bar |
| Inset | boxplot of `-log10(FDR)` distributions (Dis vs Reg), Mann–Whitney U significance |

**Gene count source:** the `Ratio in Population` field is a string like `"(34, 3938)"`;
`get_set_size()` parses the first integer (number of background genes annotated to the
term). This is the gene-set size used both for the `(n=…)` label and for the
gene-set-size sensitivity analysis in `scripts/go_param_sensitivity.py`
(Reviewer #1, Comment 4).

> **Reviewer #1, Comment 4 note:** the x-axis previously showed numeric GO IDs
> (`val[3:]`). It now shows GO term names with gene counts. The robustness of the
> Dis- vs Reg-Arrest difference to gene-set-size choices, and the specificity of the
> panel A program, are quantified in `scripts/go_param_sensitivity.py` →
> `results/GO_results/param_sensitivity/`.

---

## Panel C — Relative enrichment over time

**Function:** `panel_C(ax)`
**Data:** `results/GO_results/time_series11/GOATOOLS_GO_enrichment_results_time_series{1..8}.csv`

Tracks `-log10(FDR)` of GO terms common to all 8 SHX time points (filtered to
`t1 < 1e-6`), plotted relative to t1, with a mean ± SE ribbon. Time points (min):
`[218, 318, 426, 529, 586, 1609, 1794, 1904]`.

---

## Panel D — Lag-time survival vs time in SHX

**Function:** `panel_D(ax)`
**Data:** `scanlag_data/bulk time in shx/*.csv` (one survival curve per time point)

Log–log survival functions colored by time in SHX (Reds colormap, `0–2200` min),
with a colorbar (`×10³` units).

---

## Data files used in figure

| File | Used in |
|---|---|
| `results/deseq_results/from counts/deseq2_results_disrupted.csv` | Panel A |
| `results/deseq_results/from counts/deseq2_results_regulated.csv` | Panel A |
| `scripts/figures/figure4/GO_terms_heatmap.csv` | Panel A |
| `scripts/figures/figure4/GO_ID_name.csv` | Panel A |
| `results/GO_results/from_counts/GOATOOLS_GO_enrichment_results_disrupted_down.csv` | Panel B |
| `results/GO_results/from_counts/GOATOOLS_GO_enrichment_results_regulated_down.csv` | Panel B |
| `scripts/figures/figure4/GO_diff_pvals.csv` | Panel B |
| `results/GO_results/time_series11/GOATOOLS_GO_enrichment_results_time_series{1..8}.csv` | Panel C |
| `scanlag_data/bulk time in shx/*.csv` | Panel D |

---

## Common tweaks quick-reference

| Goal | What to change |
|---|---|
| Resize the whole figure | `figsize=(W, H)` in `PanelFigure(...)` |
| Move a panel | edit `panel_pos[i]` — `[left, bottom, width, height]` |
| Shrink/grow x-label font in panel B | `fontsize=fsize-5` in `ax.set_xticklabels(...)` |
| Give panel B labels more vertical room | raise `panel_pos[1]` bottom and/or lower panels C/D |
| Switch panel B labels back to GO IDs | replace the `labels = [...]` line with `[v[3:] for v in df.index]` |
| Change bar colors | `#de2d26` (Dis) / `#9ecae1` (Reg) in `panel_B` |
| Save outputs | `pf.save("figure4.svg", ...)` and `pf.save("figure4_preview.png", ...)` at bottom |

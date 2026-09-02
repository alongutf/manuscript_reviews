# Figure S6 — Implementation Details

Supplementary figure characterising Dis-Arrest vs Reg-Arrest cells by microscopy:
cell-length distributions (panel A), constitutive promoter activity (panel B), and
representative phase images (panel C).

Script: `scripts/supplementary_figures/figure_s6.py`

## How to run

```bash
cd scripts/supplementary_figures
python figure_s6.py        # interactive window; also writes figure_s6.pdf + figure_s6_preview.png
```

The same import bootstrap as the other supplementary scripts makes it path-independent
(`_REPO` on `sys.path` for `src.*`, plus `scripts/figures` for `figure_functions`).

---

## Inputs

| Input | Path | Used for |
|---|---|---|
| Exponential cells | `microscopy/all_positions_exp.csv` | panel A (`Exponential` group) |
| SHX experiment | `microscopy/all_positions_shx.csv` | panel A (`SHX` → Dis-Arrest, `CASP` → Reg-Arrest) |
| VapC experiment | `microscopy/all_positions_vapc.csv` | panels A + B (`VapC` → Dis-Arrest, `Reg-Arrest` → Reg-Arrest) |
| Representative images | `microscopy/images/*.png` | panel C |

All three CSVs are filtered to `kept == True` before use.

Relevant columns: `length_px` (panel A), `mcherry_bgsub_median` (panel B), `condition`, `kept`.

---

## Panel A — cell length distributions

- **Quantity:** `length_px`, converted to µm with the microscope calibration
  `PX_PER_UM = 15.15` (15.15 px = 1 µm).
- **Groups (in order):** `Exponential`, `SHX⁺` (Dis-Arrest SHX), `SHX⁻` (Reg-Arrest SHX / CASP),
  `VapC⁺ 24h`, `VapC⁻` (Reg-Arrest VapC).
- **Subsampling:** every group is randomly subsampled **without replacement** to
  `N_SUBSAMPLE = min(group sizes)` so all violins carry equal weight. The draw uses
  `np.random.default_rng(RNG_SEED)` with `RNG_SEED = 0`, so the figure is reproducible.
  With the current data the smallest group is the exponential sample (**n = 233** kept cells).
- **Annotations:** `mean` (2 decimals, µm) and `CV` above each violin. The `n=` line is
  **not** shown here — after subsampling it is identical for every group.
- **Colours:** exponential grey (`EXP_COLOR = '0.6'`), Dis-Arrest warm (`DIS_COLOR = '#E07B54'`),
  Reg-Arrest cool (`REG_COLOR = 'steelblue'`).

## Panel B — promoter activity

Violin of `mcherry_bgsub_median` for `VapC⁺ 24h` vs `VapC⁻`, same colours/order as panel A.
Not subsampled, so the annotation keeps `n=`, `mean` (0 decimals) and `CV`.

## Panel C — representative images

Single row of four phase images (`shx2.png`, `shx_reg1.png`, `vapc.png`, `vapc_reg1.png`),
converted to grayscale, in the same order and with the same labels/colours as the arrest
conditions of panel A. Frames are coloured by condition.

---

## Shared helper

`_violin_with_stats(ax, data_dict, order, colors, ylabel, title, xlabels=None, fmt='.0f', show_n=True)`
draws the violins, an IQR bar with a median dot, and the per-group annotation block.
`fmt` controls the mean's format string (`'.2f'` for panel A µm values) and `show_n`
toggles the `n=` line (`False` for the subsampled panel A).

## Layout

`PanelFigure(figsize=(7, 5), label_offset=(-0.04, 0.04))`

| Panel | Rect `[left, bottom, width, height]` |
|---|---|
| A | `[0.08, 0.62, 0.52, 0.30]` |
| B | `[0.70, 0.62, 0.26, 0.30]` |
| C | `[0.10, 0.10, 0.80, 0.40]` (axis hidden; images drawn as inset axes) |

## Outputs

- `scripts/supplementary_figures/figure_s6.pdf` (300 dpi, transparent)
- `scripts/supplementary_figures/figure_s6_preview.png` (200 dpi)

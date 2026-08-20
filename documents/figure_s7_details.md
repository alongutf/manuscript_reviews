# Figure S7 — details

Script: `scripts/supplementary_figures/figure_s7.py`
Outputs: `figure_s7.pdf`, `figure_s7_preview.png` (next to the script)

Microfluidics analysis of the SHX experiment. 1 frame = 10 minutes; SHX is added
at frame 18 (= 180 min), marked by a red dashed line in every panel.

## Layout

`PanelFigure(figsize=(7, 7), label_offset=(-0.02, 0.04))`

| Panel | Rect (`[left, bottom, width, height]`) | Content |
|-------|----------------------------------------|---------|
| A | `[0.05, 0.60, 0.88, 0.36]` | Kymograph image `microscopy/images/lineage_trench24_timeaxis.png`, with an annotated arrow at the SHX-addition line (axes-fraction coords, x ≈ 0.485) |
| B | `[0.10, 0.37, 0.85, 0.18]` | Histogram of halt times (short and wide) |
| C | `[0.10, 0.06, 0.85, 0.18]` | Histogram of division times, x-axis matched to B |

## Data

Both histograms come from `microscopy/true_events.csv` — the manually validated
event table (`confirmed_by` column). Rows are split by `event_type`
(`halt` / `division`); rows with no `frame` value are dropped (17 halt rows).
Counts used: 155 halt events, 62 division events.

`halt_times_all.csv` is no longer used by this figure.

- Bin width: 2 frames = 20 min.
- Bin edges and `xlim` are shared between B and C, computed from the maximum
  event time across both sets (halt reaches frame 40 → 400 min), so the two
  histograms are directly comparable.
- Colors: halt `#E07B54`, division `#4C86A8`.

## Reading

Divisions occur throughout the pre-drug window and stop shortly after SHX
addition; halts are concentrated just after the SHX line, peaking around
220–240 min with a tail to 400 min.
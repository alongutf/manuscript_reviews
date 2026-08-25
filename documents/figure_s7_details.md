# Figure S7 — details

Script: `scripts/supplementary_figures/figure_s7.py`
Outputs: `figure_s7.pdf`, `figure_s7_preview.png` (next to the script)

Microfluidics analysis of the SHX experiment. 1 frame = 10 minutes; SHX is added
at frame 18 (= 180 min), marked by a red dashed line in every panel.

**Time is normalized to the drug-addition time in all three panels**: 180 min is
subtracted from every time, so t = 0 is SHX addition and the pre-drug window is
negative. All three x-axes are labelled `Time from SHX addition (min)`.

## Layout

`PanelFigure(figsize=(7, 7), label_offset=(-0.02, 0.04))`

| Panel | Rect (`[left, bottom, width, height]`) | Content |
|-------|----------------------------------------|---------|
| A | `[0.05, 0.64, 0.88, 0.30]` | Kymograph image `microscopy/images/lineage_trench24_timeaxis.png`, cropped to the kymograph box (rows 18–1024, cols 18–3723) so its baked-in time axis is discarded and re-drawn in SHX-relative time; annotated arrow at t = 0 |
| B | `[0.10, 0.37, 0.85, 0.18]` | Histogram of halt times (short and wide) |
| C | `[0.10, 0.06, 0.85, 0.18]` | Histogram of division times, x-axis matched to B |

## Data

Both histograms come from `microscopy/true_events.csv` — the manually validated
event table (`confirmed_by` column). Rows are split by `event_type`
(`halt` / `division`); rows with no `frame` value are dropped (17 halt rows).
Counts used: 155 halt events, 62 division events.

`halt_times_all.csv` is no longer used by this figure.

- Bin width: 2 frames = 20 min.
- Bin edges start at t = −180 (the first frame) and `xlim` is (−180, 320); both
  are shared between B and C, computed from the maximum event time across both
  sets, so the two histograms are directly comparable.
- The dashed reference line in B and C is at t = 0.
- Colors: halt `#E07B54`, division `#4C86A8`.

## Reading

Divisions occur throughout the pre-drug window and stop shortly after SHX
addition; halts are concentrated just after t = 0, peaking around
+40–60 min with a tail to +220 min.

Note: the panel A image has a red dashed line baked in at 175 min (the frame
boundary before the first post-SHX frame). `panel_A` erases it — it lies on the
black separator between frames, so a neighbouring separator column is copied over
it — and redraws the marker with `axvline(0, lw=2.5)`, centred on the first
post-SHX frame rather than offset by half a frame.
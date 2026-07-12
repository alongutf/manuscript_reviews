"""
1-CDF (fraction of non-growing bacteria) vs. lag time for the scanlag data.

Reads scanlag_data.xlsx where:
  row 0  -> column label (condition name)
  row 2  -> t0 (subtracted from every value in that column)
  row 4+ -> lag-time measurements (minutes)

Log-log axes. Reg-Arrest is drawn in steelblue; the remaining conditions in
shades of red with distinct line styles. Output: transparent SVG next to the data.
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

DATA = r'G:\Other computers\My MacBook Air\Alon\PhD\documents\GRC conference\figures\scanlag_data.xlsx'
OUT_DIR = os.path.dirname(DATA)

# ------------------------------------------------------------------
# Load
# ------------------------------------------------------------------
raw = pd.read_excel(DATA, header=None)
labels = raw.iloc[0].tolist()
t0 = raw.iloc[2].astype(float).tolist()

# ------------------------------------------------------------------
# Styling: Reg-Arrest = steelblue; others = shades of red w/ line styles
# ------------------------------------------------------------------
RED_SHADES = ['#8B0000', '#D62728', '#F08080']   # dark -> light red
RED_STYLES = ['-', '--', ':']                     # solid / dashes / dots

fig, ax = plt.subplots(figsize=(5, 4))
xmin = 10
red_i = 0
for j, name in enumerate(labels):
    vals = pd.to_numeric(raw.iloc[4:, j], errors='coerce').dropna().values - t0[j]
    vals = vals[vals > 0]              # drop non-positive (log axis)
    x = np.sort(vals)
    # complementary CDF: fraction of bacteria with lag time > x
    y = 1.0 - np.arange(1, len(x) + 1) / len(x)
    x = np.append(xmin, x)
    y = np.append(1,y)
    if str(name).strip() == 'Reg-Arrest':
        ax.step(x, y, where='post', color='steelblue', lw=2, label=name)
    else:
        ax.step(x, y, where='post', color=RED_SHADES[red_i], lw=2,
                linestyle=RED_STYLES[red_i], label=name)
        red_i += 1

#ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('time (minutes)', fontsize=14)
ax.set_ylabel('Fraction of non-growing bacteria', fontsize=14)
ax.legend(frameon=False, fontsize=12)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xlim(left=10)
fig.tight_layout()

out = os.path.join(OUT_DIR, 'scanlag_1cdf.svg')
fig.savefig(out, transparent=True)
fig.savefig(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), 'scanlag_1cdf_preview.png'),
    dpi=200)
print('wrote', out)

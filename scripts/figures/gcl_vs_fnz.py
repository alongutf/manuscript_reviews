"""GCL vs fraction of non-zero entries, coloured by category."""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
import os

fsize = 10
root_dir = os.path.dirname(os.path.dirname(os.getcwd()))
REG_COLOR = 'steelblue'
DIS_COLOR = '#E07B54'

data = pd.read_csv(os.path.join(root_dir, 'results', 'data_metrics', 'data_metricsGCL.csv'), index_col=0)
to_exclude = ['adam_matrix_filtered.csv', 'deb_Ec_CDS_untreated.csv', 'deb_KP_CDS_untreated.csv']
data = data[~data['file_name'].isin(to_exclude)].dropna(subset=['GCL'])
x, y = data['fraction_non-zero'], data['GCL']

fig, ax = plt.subplots(figsize=(3.2, 2.8))
for cat, color, label in [('r', REG_COLOR, 'Regulated'), ('d', DIS_COLOR, 'Dis-Arrest')]:
    sub = data[data['category'] == cat]
    ax.scatter(sub['fraction_non-zero'], sub['GCL'], color=color, s=18, alpha=0.8,
               edgecolors='none', label=label)
slope, intercept = np.polyfit(x, y, 1)
xx = np.linspace(x.min(), x.max(), 50)
ax.plot(xx, slope * xx + intercept, color='grey', lw=0.75, ls='--')
rho, p_s = stats.spearmanr(x, y)
r, p_p = stats.pearsonr(x, y)
print(f'n={len(data)}, Spearman rho={rho:.3f} p={p_s:.4f}; Pearson r={r:.3f} p={p_p:.4f}')
ax.text(0.97, 0.05, f'Spearman ρ = {rho:.2f}, p = {p_s:.2f}', transform=ax.transAxes,
        ha='right', va='bottom', fontsize=fsize - 3)
ax.set_xlabel('Fraction of non-zeros', fontsize=fsize - 2)
ax.set_ylabel('GCL', fontsize=fsize - 2)
ax.tick_params(axis='both', which='major', labelsize=fsize - 2)
ax.legend(fontsize=fsize - 3, frameon=False, loc='upper left')
fig.tight_layout()
fig.savefig('gcl_vs_fnz.svg')
fig.savefig('gcl_vs_fnz_preview.png', dpi=300)

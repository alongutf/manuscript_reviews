# revision plan for figure 3
figure 3 builds on the theoretical understanding from figure 2 to show and explain the differences in overall correlations between Reg-Arrest and Dis-Arrest conditions

**main points that need to be revised**
- only show 2 plots for experimental data: 1 Reg-Arrest, 1 Dis-Arrest
- plot eigenvalue distributions in the ccdf format
- Need to show simulation results for dysregulated (low $\rho$) and regulated (high $\rho$)
- Present the $\rho$ to GMP-Cor calibration curve. Show where Dis-Arrest/Reg-Arrest datasets fall
- Show the box-plot for GMP-Cor values for all analyzed datasets. add significance test

**New panels**
- Panel A-B replace panels A-C removing the exponential data plot (panel A) and changing all plots to ccdf format. Use data from sample_15b_filtered (regulated) and sample_15a_filtered (disrupted). The data is in `ev_data\`.
  - Panel C-D simulated data results with high and low $\rho$ respectively. The data is in `ev_data\`.
- Panel E: $\rho$ sweep calibration curve. the data from the sweep is saved in `C:\Users\owner\Documents\Projects\manuscript_reviews\results\simulation_results\raw\rho_sweep_raw.txt`
Use this code:
```python
# Calibration curve: GMP-Cor vs. correlation strength (rho)
import pandas as pd
import os
import matplotlib.pyplot as plt
RESULTS_DIR = os.path.join(os.path.dirname(os.getcwd()), 'results', 'simulation_results')
summary = pd.read_csv(os.path.join(RESULTS_DIR, 'raw', 'rho_sweep_summary.txt'),
                      sep=r'\s+', comment='#', index_col=0)

rho_vals = summary.index.values
medians  = summary['median'].values
stds     = summary['std'].values

fig, ax = plt.subplots(figsize=(7, 4))
ax.errorbar(rho_vals, medians, yerr=stds, fmt='o-', color='steelblue',
            capsize=4, linewidth=1.8, markersize=6, label='median ± SD')
ax.fill_between(rho_vals, medians - stds, medians + stds, alpha=0.15, color='steelblue')
ax.set_xlabel(r'Correlation strength ($\rho$)', fontsize=13)
ax.set_ylabel('GMP-Cor', fontsize=13)
ax.set_title('GMP-Cor calibration curve', fontsize=14)
ax.set_xticks(rho_vals)
ax.tick_params(axis='x', rotation=45)
ax.set_xlim(-0.02, 1.02)
ax.set_ylim(bottom=-1)
ax.grid(True, linestyle='--', alpha=0.2)
ax.legend(fontsize=10)
plt.tight_layout()
plt.savefig(os.path.join(RESULTS_DIR, 'figures', 'rho_sweep_calibration.svg'), dpi=150)
plt.show()
```
- Panel F: Box plot with GMP-cor values from all datasets. use the following code: 
```python 
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
def format_p(p):
    if p < 0.0001:
        return '****'
    elif p < 0.001:
        return '***'
    elif p < 0.01:
        return '**'
    elif p < 0.05:
        return '*'
    else:
        return 'NS'
fsize = 12
fig,ax = plt.subplots(figsize=(8, 6))
path = r"C:\Users\owner\Documents\Projects\manuscript_reviews\results\data_metrics\test8.csv"
data = pd.read_csv(path, index_col=0)
ranking_param = 'sum_denoised_ev'
# Sort values and compute rank
data['Rank'] = data[ranking_param].rank(method='min').astype(int)
# Split data into two groups based on labels
group1 = data[data['category'] == 'r'][ranking_param]
group0 = data[data['category'] == 'd'][ranking_param]
c = "k"
box = ax.boxplot([group1, group0], meanline=True, showmeans=True, patch_artist=True,
                 boxprops=dict(facecolor="None", color=c), whiskerprops=dict(color=c), capprops=dict(color=c),
                 flierprops=dict(markeredgecolor=c, markersize=2), medianprops=dict(color=c))
for element in ['boxes', 'whiskers', 'caps', 'medians']:
    for item in box[element]:
        item.set_linewidth(1)
for mean_line in box['means']:
    mean_line.set_linewidth(1)
    mean_line.set_color(c)
    mean_line.set_linestyle('solid')

ax.set_xticklabels(['Regulated', 'Dis-Arrest'], fontsize=fsize, rotation=0, ha='center')
#ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1])
#ax.set_ylim([0.2, 1])
#ax.set_yticklabels([0.2, 0.4, 0.6, 0.8, 1], fontsize=fsize - 2)
# add u-test results
u_stat, u_p = stats.mannwhitneyu(group1, group0)

# Define the level of significance
asterisks = format_p(u_p)
# Add significance annotation
x1, x2 = 1, 2  # x-coordinates of the box plots
y, h = max(max(group1), max(group0))+1, 1  # y-position and height of the annotation
ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=.75, color='black')
ax.text((x1 + x2) * 0.5, y + h, asterisks, ha='center', va='bottom', color='black', fontsize=fsize - 2)
ax.set_ylabel(r'Metric', fontsize=fsize)
ax.grid(False)
```

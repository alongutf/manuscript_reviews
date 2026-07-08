import numpy as np
import matplotlib.pyplot as plt
import os

project_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
output_dir = "/Users/along/Documents/Alon/PhD/documents/IPS 2026"

fsize = 12
bins = np.linspace(0, 1000, 40)

hist_data = np.loadtxt(os.path.join(project_dir, 'scanlag_data', 'kaplan_shx', 't0.txt'))
hist_data2 = np.loadtxt(os.path.join(project_dir, 'scanlag_data', 'kaplan_shx', 't1346.txt'))

x0 = hist_data[0]
hist_data = hist_data - x0
hist_data2 = hist_data2 - x0

def make_hist_fig(data, color, xlabel=False):
    fig, ax = plt.subplots(figsize=(3.5, 2.5))
    ax.hist(data, bins=bins, color=color, alpha=0.7)
    ax.set_ylabel('# of cells', fontsize=fsize, labelpad=0)
    ax.set_xlim(0, 1000)
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='both', which='major', labelsize=fsize - 2)
    if xlabel:
        ax.set_xlabel('Lag time (min)', fontsize=fsize, labelpad=0)
    else:
        ax.set_xticks([])
    fig.tight_layout()
    return fig

fig1 = make_hist_fig(hist_data, color='steelblue', xlabel=True)
fig1.savefig(os.path.join(output_dir, 'panel_B_t0.png'), dpi=300, bbox_inches='tight')
plt.close(fig1)

fig2 = make_hist_fig(hist_data2, color='#CD5C5C', xlabel=True)
fig2.savefig(os.path.join(output_dir, 'panel_B_t1346.png'), dpi=300, bbox_inches='tight')
plt.close(fig2)

print("Saved panel_B_t0.png and panel_B_t1346.png to", output_dir)
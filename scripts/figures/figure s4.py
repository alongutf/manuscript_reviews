import matplotlib.pyplot as plt
import matplotlib as mpl
from figure_functions import PanelFigure
import numpy as np
import pandas as pd
import os
from scipy import stats
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def panel_A(ax):
    # Load the DEG results
    file_name = 'deseq2_results_regulated.csv'
    path = os.path.join(root_dir, 'results', 'deseq_results', 'from counts', file_name)
    deg_results = pd.read_csv(path, index_col=0)
    # Add a new column for -log10(padj)
    deg_results['-log10(padj)'] = -np.log10(deg_results['padj'].replace(0, 1e-300))  # Avoid log(0)
    # plot the volcano plot
    # values for cutoffs
    p_cutoff = 0.01
    lfc_cutoff = 1
    # Separate the points above the cutoffs
    not_significant = np.logical_or(deg_results['padj'] > p_cutoff, deg_results['log2FoldChange'].abs() < lfc_cutoff)
    gray = deg_results[not_significant]
    highlight = deg_results[~not_significant]
    # Plot gray points (background)
    ax.scatter(gray['log2FoldChange'], gray['-log10(padj)'],
               color='gray', alpha=0.3, s=marker_size)

    # Plot orange (highlighted) points
    ax.scatter(highlight['log2FoldChange'], highlight['-log10(padj)'],
               color='orange', edgecolor='orange', alpha=0.5, s=marker_size, zorder=3)
    ax.set_title('Reg-Arrest', fontsize=fsize)
    # plot second GO term:

    # Add cutoff lines
    ax.axhline(y=-np.log10(p_cutoff), linestyle='--', color='black', label='P-value cutoff', lw=1)
    ax.axvline(x=lfc_cutoff, linestyle='--', color='black', label='Log2FC cutoff', lw=1)
    ax.axvline(x=-lfc_cutoff, linestyle='--', color='black', lw=1)
    ax.set_xticks([-5, 0, 5])
    ax.set_xticklabels([-5, 0, 5], fontsize=fsize-2)
    ax.set_yticks([0, 100, 200])
    ax.set_yticklabels([0, 1, 2], fontsize=fsize-2)
    ax.set_xlabel(r'Log$_2$ Fold Change', fontsize=fsize, labelpad=0)
    ax.set_xlim([-10, 10])
    ax.set_ylim([-10, 250])
    ax.set_ylabel(r'$-log_{10}{p_{adj}}$', fontsize=fsize, labelpad=0)
    ax.text(-10,245, r'$\times10^2$', fontsize=fsize-2, color='black', ha='left', va='top')

def panel_B(ax):
    # Load the DEG results
    file_name = 'deseq2_results_disrupted.csv'
    path = os.path.join(root_dir, 'results', 'deseq_results', 'from counts', file_name)
    deg_results = pd.read_csv(path, index_col=0)
    # Add a new column for -log10(padj)
    deg_results['-log10(padj)'] = -np.log10(deg_results['padj'].replace(0, 1e-300))  # Avoid log(0)
    # plot the volcano plot
    # values for cutoffs
    p_cutoff = 0.01
    lfc_cutoff = 1
    # Separate the points above the cutoffs
    not_significant = np.logical_or(deg_results['padj'] > p_cutoff, deg_results['log2FoldChange'].abs() < lfc_cutoff)
    gray = deg_results[not_significant]
    highlight = deg_results[~not_significant]
    # Plot gray points (background)
    ax.scatter(gray['log2FoldChange'], gray['-log10(padj)'],
               color='gray', alpha=0.3, s=marker_size)

    # Plot orange (highlighted) points
    ax.scatter(highlight['log2FoldChange'], highlight['-log10(padj)'],
               color='orange', edgecolor='orange', alpha=0.5, s=marker_size, zorder=3)
    ax.set_title('Dis-Arrest', fontsize=fsize)
    # plot second GO term:

    # Add cutoff lines
    ax.axhline(y=-np.log10(p_cutoff), linestyle='--', color='black', label='P-value cutoff', lw=1)
    ax.axvline(x=lfc_cutoff, linestyle='--', color='black', label='Log2FC cutoff', lw=1)
    ax.axvline(x=-lfc_cutoff, linestyle='--', color='black', lw=1)
    ax.set_xticks([-5, 0, 5])
    ax.set_xticklabels([-5, 0, 5], fontsize=fsize-2)
    ax.set_yticks([0, 100, 200])
    ax.set_yticklabels([0, 1, 2], fontsize=fsize-2)
    ax.set_xlabel(r'Log$_2$ Fold Change', fontsize=fsize, labelpad=0)
    ax.set_xlim([-10, 10])
    ax.set_ylim([-10, 250])
#    ax.set_ylabel(r'$-log_{10}{p_{adj}}$', fontsize=fsize, labelpad=0)
    ax.text(-10,245, r'$\times10^2$', fontsize=fsize-2, color='black', ha='left', va='top')

def panel_C(ax):
    # compare p-values of GO terms
    # Load the GO enrichment results
    file_name = 'GOATOOLS_GO_enrichment_results_disrupted_down.csv'
    path = os.path.join(root_dir, 'results', 'GO_results', 'from_counts', file_name)
    go_results = pd.read_csv(path)
    # load the second condition
    file_name = 'GOATOOLS_GO_enrichment_results_regulated_down.csv'
    path = os.path.join(root_dir, 'results', 'GO_results', 'from_counts', file_name)
    go_results2 = pd.read_csv(path)
    # filter go terms that appear in both conditions
    common_go_terms = set(go_results['GO_ID']).intersection(set(go_results2['GO_ID']))
    # get the p-values of the common go terms
    go_dict = {}
    for go_term in common_go_terms:
        p_value1 = go_results[go_results['GO_ID'] == go_term]['FDR'].values[0]
        p_value2 = go_results2[go_results2['GO_ID'] == go_term]['FDR'].values[0]
        go_dict[go_term] = [p_value1, p_value2]

    # create a dataframe
    df = pd.DataFrame(go_dict).T
    df.columns = ['SHX', 'Casp']
    df = df.sort_values('Casp', ascending=True)
    # plot bar plot of -log10(p-values) for each GO term
    # Define bar width
    bar_width = 0.35

    # Create bar positions for each group
    x = np.arange(len(df.index))  # the label locations
    shx_values = -np.log10(df['SHX'])
    casp_values = -np.log10(df['Casp'])

    # Plot each bar next to each other
    ax.bar(x - bar_width / 2, shx_values, width=bar_width, label='Dis-Arrest', color='#de2d26')
    ax.bar(x + bar_width / 2, casp_values, width=bar_width, label='Reg-Arrest', color='#9ecae1')

    # Set title and labels
    ax.set_ylabel('$-\log_{10}{p_{adj}}$', fontsize=fsize, labelpad=0)
    ax.set_xlabel('GO ID', fontsize=fsize, labelpad=0)
    ax.set_yticks([4, 8, 12, 16])
    ax.set_xticks(x)  # Set x-ticks to be at the center of each pair
    labels = [val[3:] for val in df.index]
    ax.set_xticklabels(labels, fontsize=fsize - 2, rotation=45, ha='right')  # Set x-tick labels to match the DataFrame index
    # Add legend
    ax.legend(fontsize=fsize-2, loc='upper center', bbox_to_anchor=(0.35, 1.015), ncol=1)
    ax.set_yticklabels([4, 8, 12, 16], fontsize=fsize-2)
    # add inset
    # get samples
    sample1 = np.flip(np.sort(-np.log10(go_results['FDR']).values))
    sample2 = np.flip(np.sort(-np.log10(go_results2['FDR']).values))
    n = np.minimum(len(sample1), len(sample2))
    sample1 = sample1[:n]
    sample2 = sample2[:n]
    # plot boxplots
    inset_ax = inset_axes(ax, width="30%", height="55%", loc='upper right')
    c = "k"
    box = inset_ax.boxplot([sample1, sample2], meanline=True, showmeans=True, patch_artist=True,
                     boxprops=dict(facecolor="None", color=c), whiskerprops=dict(color=c), capprops=dict(color=c),
                     flierprops=dict(markeredgecolor=c, markersize=1), medianprops=dict(color=c))

    for element in ['boxes', 'whiskers', 'caps', 'medians']:
        for item in box[element]:
            item.set_linewidth(0.5)

    for mean_line in box['means']:
        mean_line.set_linewidth(1)
        mean_line.set_color(c)
        mean_line.set_linestyle('solid')

    inset_ax.set_xticklabels(['Dis-Arrest', 'Reg-Arrest'], fontsize=fsize - 2, rotation=15, ha='center')
    inset_ax.set_ylabel('$-\log_{10}{p_{adj}}$', fontsize=fsize - 2, labelpad=0)
    inset_ax.set_yticks([4, 8, 12, 16], )
    inset_ax.set_yticklabels([4, 8, 12, 16], fontsize=fsize - 2)
    # add u-test results
    u_stat, u_p = stats.mannwhitneyu(sample1, sample2)

    # Define the level of significance
    def get_asterisks(p_value):
        if p_value < 0.001:
            return '***'
        elif p_value < 0.01:
            return '**'
        elif p_value < 0.05:
            return '*'
        else:
            return 'ns'  # Not significant

    asterisks = get_asterisks(u_p)
    # Add significance annotation
    x1, x2 = 1, 2  # x-coordinates of the box plots
    y, h = max(max(sample1), max(sample2)) , 0.2  # y-position and height of the annotation
    inset_ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=.5, color='black')
    inset_ax.text((x1 + x2) * 0.5, y + h / 2, asterisks, ha='center', va='bottom', color='black', fontsize=fsize-4)
    inset_ax.set_ylim(0, 20)


def panel_F(ax):
    # compare p-values of GO terms
    # Load the GO enrichment results
    file_name = 'GOATOOLS_GO_enrichment_results_disrupted_down.csv'
    path = os.path.join(root_dir, 'results', 'GO_results', 'ercc_norm', file_name)
    go_results = pd.read_csv(path)
    # load the second condition
    file_name = 'GOATOOLS_GO_enrichment_results_regulated_down.csv'
    path = os.path.join(root_dir, 'results', 'GO_results', 'ercc_norm', file_name)
    go_results2 = pd.read_csv(path)
    # filter go terms that appear in both conditions
    common_go_terms = set(go_results['GO_ID']).intersection(set(go_results2['GO_ID']))
    # get the p-values of the common go terms
    go_dict = {}
    for go_term in common_go_terms:
        p_value1 = go_results[go_results['GO_ID'] == go_term]['FDR'].values[0]
        p_value2 = go_results2[go_results2['GO_ID'] == go_term]['FDR'].values[0]
        go_dict[go_term] = [p_value1, p_value2]

    # create a dataframe
    df = pd.DataFrame(go_dict).T
    df.columns = ['SHX', 'Casp']
    df = df.sort_values('Casp', ascending=True)
    # plot bar plot of -log10(p-values) for each GO term
    # Define bar width
    bar_width = 0.35

    # Create bar positions for each group
    x = np.arange(len(df.index))  # the label locations
    shx_values = -np.log10(df['SHX'])
    casp_values = -np.log10(df['Casp'])

    # Plot each bar next to each other
    ax.bar(x - bar_width / 2, shx_values, width=bar_width, label='Dis-Arrest', color='#de2d26')
    ax.bar(x + bar_width / 2, casp_values, width=bar_width, label='Reg-Arrest', color='#9ecae1')

    # Set title and labels
    ax.set_ylabel('$-\log_{10}{p_{adj}}$', fontsize=fsize, labelpad=0)
    ax.set_xlabel('GO ID', fontsize=fsize, labelpad=0)
    ax.set_yticks([4, 8, 12, 16])
    ax.set_xticks(x)  # Set x-ticks to be at the center of each pair
    labels = [val[3:] for val in df.index]
    ax.set_xticklabels(labels, fontsize=fsize - 2, rotation=45, ha='right')  # Set x-tick labels to match the DataFrame index
    # Add legend
    ax.legend(fontsize=fsize-2, loc='upper center', bbox_to_anchor=(0.35, 1.015), ncol=1)
    ax.set_yticklabels([4, 8, 12, 16], fontsize=fsize-2)
    # add inset
    # get samples
    sample1 = np.flip(np.sort(-np.log10(go_results['FDR']).values))
    sample2 = np.flip(np.sort(-np.log10(go_results2['FDR']).values))
    n = np.minimum(len(sample1), len(sample2))
    sample1 = sample1[:n]
    sample2 = sample2[:n]
    # plot boxplots
    inset_ax = inset_axes(ax, width="30%", height="55%", loc='upper right')
    c = "k"
    box = inset_ax.boxplot([sample1, sample2], meanline=True, showmeans=True, patch_artist=True,
                     boxprops=dict(facecolor="None", color=c), whiskerprops=dict(color=c), capprops=dict(color=c),
                     flierprops=dict(markeredgecolor=c, markersize=1), medianprops=dict(color=c))

    for element in ['boxes', 'whiskers', 'caps', 'medians']:
        for item in box[element]:
            item.set_linewidth(0.5)

    for mean_line in box['means']:
        mean_line.set_linewidth(1)
        mean_line.set_color(c)
        mean_line.set_linestyle('solid')

    inset_ax.set_xticklabels(['Dis-Arrest', 'Reg-Arrest'], fontsize=fsize - 2, rotation=15, ha='center')
    inset_ax.set_ylabel('$-\log_{10}{p_{adj}}$', fontsize=fsize - 2, labelpad=0)
    inset_ax.set_yticks([4, 8, 12, 16], )
    inset_ax.set_yticklabels([4, 8, 12, 16], fontsize=fsize - 2)
    # add u-test results
    u_stat, u_p = stats.mannwhitneyu(sample1, sample2)

    # Define the level of significance
    def get_asterisks(p_value):
        if p_value < 0.001:
            return '***'
        elif p_value < 0.01:
            return '**'
        elif p_value < 0.05:
            return '*'
        else:
            return 'ns'  # Not significant

    asterisks = get_asterisks(u_p)
    # Add significance annotation
    x1, x2 = 1, 2  # x-coordinates of the box plots
    y, h = max(max(sample1), max(sample2)) , 0.2  # y-position and height of the annotation
    inset_ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=.5, color='black')
    inset_ax.text((x1 + x2) * 0.5, y + h / 2, asterisks, ha='center', va='bottom', color='black', fontsize=fsize-4)
    inset_ax.set_ylim(0, 20)

def panel_D(ax):
    # Load the DEG results
    file_name = 'deseq2_results_regulated.csv'
    path = os.path.join(root_dir, 'results', 'deseq_results', 'ercc_norm', file_name)
    deg_results = pd.read_csv(path, index_col=0)
    # Add a new column for -log10(padj)
    deg_results['-log10(padj)'] = -np.log10(deg_results['padj'].replace(0, 1e-300))  # Avoid log(0)
    # plot the volcano plot
    # values for cutoffs
    p_cutoff = 0.01
    lfc_cutoff = 1
    # Separate the points above the cutoffs
    not_significant = np.logical_or(deg_results['padj'] > p_cutoff, deg_results['log2FoldChange'].abs() < lfc_cutoff)
    gray = deg_results[not_significant]
    highlight = deg_results[~not_significant]
    # Plot gray points (background)
    ax.scatter(gray['log2FoldChange'], gray['-log10(padj)'],
               color='gray', alpha=0.3, s=marker_size)

    # Plot orange (highlighted) points
    ax.scatter(highlight['log2FoldChange'], highlight['-log10(padj)'],
               color='orange', edgecolor='orange', alpha=0.5, s=marker_size, zorder=3)
    ax.set_title('Reg-Arrest', fontsize=fsize)
    # plot second GO term:

    # Add cutoff lines
    ax.axhline(y=-np.log10(p_cutoff), linestyle='--', color='black', label='P-value cutoff', lw=1)
    ax.axvline(x=lfc_cutoff, linestyle='--', color='black', label='Log2FC cutoff', lw=1)
    ax.axvline(x=-lfc_cutoff, linestyle='--', color='black', lw=1)
    ax.set_xticks([-5, 0, 5])
    ax.set_xticklabels([-5, 0, 5], fontsize=fsize-2)
    ax.set_yticks([0, 100, 200])
    ax.set_yticklabels([0, 1, 2], fontsize=fsize-2)
    ax.set_xlabel(r'Log$_2$ Fold Change', fontsize=fsize, labelpad=0)
    ax.set_xlim([-10, 10])
    ax.set_ylim([-10, 250])
    ax.set_ylabel(r'$-log_{10}{p_{adj}}$', fontsize=fsize, labelpad=0)
    ax.text(-10,245, r'$\times10^2$', fontsize=fsize-2, color='black', ha='left', va='top')

def panel_E(ax):
    # Load the DEG results
    file_name = 'deseq2_results_disrupted.csv'
    path = os.path.join(root_dir, 'results', 'deseq_results', 'ercc_norm', file_name)
    deg_results = pd.read_csv(path, index_col=0)
    # Add a new column for -log10(padj)
    deg_results['-log10(padj)'] = -np.log10(deg_results['padj'].replace(0, 1e-300))  # Avoid log(0)
    # plot the volcano plot
    # values for cutoffs
    p_cutoff = 0.01
    lfc_cutoff = 1
    # Separate the points above the cutoffs
    not_significant = np.logical_or(deg_results['padj'] > p_cutoff, deg_results['log2FoldChange'].abs() < lfc_cutoff)
    gray = deg_results[not_significant]
    highlight = deg_results[~not_significant]
    # Plot gray points (background)
    ax.scatter(gray['log2FoldChange'], gray['-log10(padj)'],
               color='gray', alpha=0.3, s=marker_size)

    # Plot orange (highlighted) points
    ax.scatter(highlight['log2FoldChange'], highlight['-log10(padj)'],
               color='orange', edgecolor='orange', alpha=0.5, s=marker_size, zorder=3)
    ax.set_title('Dis-Arrest', fontsize=fsize)
    # plot second GO term:

    # Add cutoff lines
    ax.axhline(y=-np.log10(p_cutoff), linestyle='--', color='black', label='P-value cutoff', lw=1)
    ax.axvline(x=lfc_cutoff, linestyle='--', color='black', label='Log2FC cutoff', lw=1)
    ax.axvline(x=-lfc_cutoff, linestyle='--', color='black', lw=1)
    ax.set_xticks([-5, 0, 5])
    ax.set_xticklabels([-5, 0, 5], fontsize=fsize-2)
    ax.set_yticks([0, 100, 200])
    ax.set_yticklabels([0, 1, 2], fontsize=fsize-2)
    ax.set_xlabel(r'Log$_2$ Fold Change', fontsize=fsize, labelpad=0)
    ax.set_xlim([-10, 10])
    ax.set_ylim([-10, 250])
#   ax.set_ylabel(r'$-log_{10}{p_{adj}}$', fontsize=fsize, labelpad=0)
    ax.text(-10,245, r'$\times10^2$', fontsize=fsize-2, color='black', ha='left', va='top')
###
# Build figure s4:
fsize = 10
marker_size=2.5
plt.close("all")
root_dir = os.path.dirname(os.path.dirname(os.getcwd()))
pf = PanelFigure(figsize=(7, 5), label_offset=(0, 0.04))
panel_pos = [
    [0.05, 0.62, 0.17, 0.32],  # A
    [0.25, 0.62, 0.17, 0.32],  # B
    [0.5, 0.62, 0.48, 0.32],  # C
    [0.05, 0.13, 0.17, 0.32],  # D
    [0.25, 0.13, 0.17, 0.32],  # E
    [0.5, 0.13, 0.48, 0.32],  # F
]
# panel A:
pf.add_panel(panel_pos[0], draw_func=panel_A, label="A")
# panel B:
pf.add_panel(panel_pos[1], draw_func=panel_B, label="B")
# panel C:
pf.add_panel(panel_pos[2], draw_func=panel_C, label="C")
# panel D:
pf.add_panel(panel_pos[3], draw_func=panel_D, label="D")
# panel E:
pf.add_panel(panel_pos[4], draw_func=panel_E, label="E")
# panel F:
pf.add_panel(panel_pos[5], draw_func=panel_F, label="F")
pf.save("figure s4.pdf", dpi=300)
plt.show()
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from scipy.cluster.hierarchy import linkage, leaves_list
from figure_functions import PanelFigure
import numpy as np
import pandas as pd
import os
from scipy import stats
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from goatools import obo_parser
from goatools.associations import read_gaf

# ------------------------------------------------------------------
# BUILD FIGURE
# ------------------------------------------------------------------
def get_gene_dict():
    file_name = 'genomic.gtf'
    # get project path

    path = os.path.join(root_dir, 'metadata', file_name)
    df = pd.read_csv(path, sep='\t', comment='#', header=None)
    df.columns = ['seqid', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    df = df[df['feature'] == 'CDS']
    # create dictionary to store gene id and gene name
    gene_id_name = {}
    for index, row in df.iterrows():
        gene_id = row['attributes'].split(';')[2].split(' ')[2].replace('"', '')
        # remove 'GeneID:' from gene_id
        gene_id = gene_id.split(':')[-1]
        gene_name = row['attributes'].split(';')[4].split(' ')[2].replace('"', '')
        gene_id_name[gene_name.lower()] = gene_id
    return gene_id_name


def get_GO_gene_list(go_term):
    # path to metadata files
    dir = os.path.join(root_dir, 'metadata')
    # Define file paths
    GO_OBO = "go-basic.obo"
    GAF_FILE = "ecocyc.gaf"  # Replace with your GAF file
    # Load GO DAG (ontology structure)
    go_dag = obo_parser.GODag(os.path.join(dir, GO_OBO))

    # Load gene-GO associations (you can download GAF files from Gene Ontology site)
    assoc = read_gaf(os.path.join(dir, GAF_FILE))

    # Get genes annotated to the term or its children
    def get_genes_for_term(term, associations, go_dag):
        terms = {t for t in go_dag[term].get_all_children()}
        terms.add(term)
        genes = {gene for gene, go_terms in associations.items() if go_terms & terms}
        return genes

    genes = get_genes_for_term(go_term, assoc, go_dag)
    return list(genes)


def panel_A(ax):
    file_name = 'deseq2_results_disrupted.csv'
    path = os.path.join(root_dir, 'results', 'deseq_results', 'from counts', file_name)
    deg_results1 = pd.read_csv(path, index_col=0)
    file_name = 'deseq2_results_regulated.csv'
    path = os.path.join(root_dir, 'results', 'deseq_results', 'from counts', file_name)
    deg_results2 = pd.read_csv(path, index_col=0)

    merged_df = pd.merge(deg_results1, deg_results2, left_index=True, right_index=True, how='inner')
    GO_terms_to_plot = ['GO:0006935']
    gene_id_name = get_gene_dict()
    gene_names = merged_df.index.str.lower()
    gene_ID = []
    for gene in gene_names:
        try:
            gene_ID.append(gene_id_name[gene.lower()])
        except:
            gene_ID.append('not_found')
            pass
    merged_df['gene_id'] = gene_ID
    GO_term = np.full(len(merged_df), 'GO:0000000')
    # Add a new column for GO term
    for GO in GO_terms_to_plot:
        gene_list = get_GO_gene_list(GO)
        for i, gene in enumerate(gene_ID):
            if gene in gene_list:
                GO_term[i] = GO
    merged_df['GO_term'] = GO_term
    merged_df_GO = merged_df[merged_df['GO_term'] != 'GO:0000000']
    # sort by GO term
    merged_df_GO = merged_df_GO.sort_values('GO_term')
    # leave only log2foldchange columns
    merged_df_GO = merged_df_GO[['log2FoldChange_x', 'log2FoldChange_y', 'GO_term']]
    # plot heatmap of log2foldchange
    # Load and sort by GO term
    go_id_name = pd.read_csv(os.path.join(root_dir, 'scripts', 'figures', 'figure4', 'GO_ID_name.csv'), index_col=0,
                             header=0, usecols=[0, 1])
    go_id_name = go_id_name.transpose()

    df = merged_df_GO
    df_sorted = df.sort_values(by='GO_term')
    go_terms = df_sorted['GO_term'].values

    heatmap_data = df_sorted.drop(columns=['GO_term'])
    # Prepare for new ordering
    new_order = []
    # Get the unique GO terms in order
    unique_go_terms = pd.unique(go_terms)
    # Cluster genes within each GO group
    for term in unique_go_terms:
        # Extract the genes in this GO group
        group = heatmap_data.loc[df_sorted['GO_term'] == term, ['log2FoldChange_y']]
        # If group has only 1 gene, just take it
        if group.shape[0] == 1:
            new_order.append(group.index[0])
        else:
            # Perform hierarchical clustering within the group
            linkage_matrix = linkage(group.values, method='average', metric='euclidean')
            ordered_idx = leaves_list(linkage_matrix)
            new_order.extend(group.index[ordered_idx])

    # Reorder the heatmap data
    heatmap_data_reordered = heatmap_data.loc[new_order]
    go_terms_reordered = df_sorted.loc[new_order, 'GO_term'].values

    # Now find new boundaries and labels
    boundaries = []
    group_labels = []
    group_positions = []
    prev_term = go_terms_reordered[0]
    start_idx = 0

    for i, term in enumerate(go_terms_reordered):
        if term != prev_term:
            boundaries.append(i)
            group_labels.append(prev_term)
            group_positions.append((start_idx + i - 1) / 2)
            start_idx = i
            prev_term = term

    group_labels.append(prev_term)
    group_positions.append((start_idx + len(go_terms_reordered) - 1) / 2)

    # Plot
    heatmap = sns.heatmap(
        heatmap_data_reordered,
        cmap='coolwarm',
        center=0,
        ax=ax,
        cbar_kws={'label': 'log$_2$ Fold Change', 'shrink': 0.5},
        linewidths=0.5,
        yticklabels=True,
        vmax=2,
        vmin=-4
    )
    # Access the colorbar
    cbar = heatmap.collections[0].colorbar

    # Now modify the colorbar
    cbar.ax.tick_params(labelsize=fsize-2)  # smaller font size for tick labels
    cbar.set_label('log2 Fold Change', fontsize=fsize-2)  # smaller label font size

    # (Optional) You can also shrink the colorbar
    #cbar.ax.set_aspect(10)  # Make it thinner vertically
    # Format x-axis
    ax.set_xticks(np.arange(heatmap_data.shape[1]) + 0.5)
    ax.set_xticklabels(['Dis-Arrest', 'Reg-Arrest'], rotation=45, ha='right', fontsize=fsize)
    ax.tick_params(axis='y', labelsize=fsize-3)
    for b in boundaries:
        ax.axhline(b, color='black', linewidth=lw)
    # Draw side brackets and labels
    for i, (start, end) in enumerate(zip([0] + boundaries, boundaries + [heatmap_data.shape[0]])):
        x = -1.5
        ax.plot([x, x], [start + .5, end - .5], color='black', linewidth=lw, clip_on=False)
        tick_size = 0.1
        ax.plot([x, x + tick_size], [start + .5, start + .5], color='black', linewidth=lw, clip_on=False)
        ax.plot([x, x + tick_size], [end - .5, end - .5], color='black', linewidth=lw, clip_on=False)

        mid = (start + end) / 2
        ax.plot([x, x - tick_size], [mid, mid], color='black', linewidth=lw, clip_on=False)
        ax.text(x - tick_size, mid, go_id_name[group_labels[i]].values[0],
                va='center', ha='right', fontsize=fsize-2, clip_on=False, rotation=90)

def panel_B(ax):
    def get_asterisks(p_value):
        if p_value < 0.001:
            return '***'
        elif p_value < 0.01:
            return '**'
        elif p_value < 0.05:
            return '*'
        else:
            return 'ns'  # Not significant
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
    h=0.1
    # Create bar positions for each group
    x = np.arange(len(df.index))  # the label locations
    shx_values = -np.log10(df['SHX'])
    casp_values = -np.log10(df['Casp'])
    labels = [val[3:] for val in df.index]
    # Plot each bar next to each other
    ax.bar(x - bar_width / 2, shx_values, width=bar_width, label='Dis-Arrest', color='#de2d26')
    ax.bar(x + bar_width / 2, casp_values, width=bar_width, label='Reg-Arrest', color='#9ecae1')
    # significance annotations:
    GO_pvals = pd.read_csv(os.path.join(root_dir,'scripts','figures','figure4','GO_diff_pvals.csv'),index_col=0,header=0)
    for pos,ID,y_shx,y_casp in zip(x,labels,shx_values,casp_values):
        y = max(y_shx,y_casp)
        p_adj = GO_pvals['p_adj'][GO_pvals.index == 'GO:'+ID].values[0]
        if p_adj<0.05:
            ax.plot([pos-bar_width/2,pos-bar_width/2, pos+bar_width/2, pos+bar_width/2], [y+h,y+2*h,y+2*h,y+h],color='k', linewidth=0.5)
            ax.text(pos,y+2*h,get_asterisks(p_adj), ha='center', fontsize = fsize-3)

    # Set title and labels
    ax.set_ylabel('Enrichment score', fontsize=fsize, labelpad=0)
    ax.set_xlabel('GO ID', fontsize=fsize, labelpad=0)
    ax.set_ylim([0,19])
    ax.set_yticks([4, 8, 12, 16])
    ax.set_xticks(x)  # Set x-ticks to be at the center of each pair
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
    inset_ax.set_ylabel('Enrichment score', fontsize=fsize - 2, labelpad=0)
    inset_ax.set_yticks([4, 8, 12, 16], )
    inset_ax.set_yticklabels([4, 8, 12, 16], fontsize=fsize - 2)
    # add u-test results
    u_stat, u_p = stats.mannwhitneyu(sample1, sample2)

    # Define the level of significance
    asterisks = get_asterisks(u_p)
    # Add significance annotation
    x1, x2 = 1, 2  # x-coordinates of the box plots
    y, h = max(max(sample1), max(sample2)) , 0.2  # y-position and height of the annotation
    inset_ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=.5, color='black')
    inset_ax.text((x1 + x2) * 0.5, y + h / 2, asterisks, ha='center', va='bottom', color='black', fontsize=fsize-3)
    inset_ax.set_ylim(0, 20)



def panel_C(ax):
    # compare GO term FDR's over time
    # Load the GO enrichment results
    file_name = 'GOATOOLS_GO_enrichment_results_time_series'
    time = [218, 318, 426, 529, 586, 1609, 1794, 1904]
    path = os.path.join(root_dir, 'results', 'GO_results','time_series11', file_name)
    go_results = pd.read_csv(path + '1.csv')
    common_go_terms = set(go_results['GO_ID'])
    for i in range(2, 9):
        go_results = pd.read_csv(path + f'{i}.csv')
        common_go_terms = set(go_results['GO_ID']).intersection(common_go_terms)
    common_go_df = pd.DataFrame(index=list(common_go_terms))
    for i in range(1, 9):
        go_results = pd.read_csv(path + f'{i}.csv')
        go_results = go_results[go_results['GO_ID'].isin(common_go_terms)]
        go_results = go_results.set_index('GO_ID')
        go_results = go_results['FDR']
        go_results.name = f't{i}'
        common_go_df = common_go_df.join(go_results)

    common_go_df = common_go_df.loc[common_go_df['t1'] < 1e-6]
    # plot average FDR of GO terms over time with error ribbon
    data = common_go_df
    data = -np.log10(data)
    data = data - np.tile(np.array(data.iloc[:, 0]), (data.shape[1], 1)).T
    # plot the average FDR of GO terms over time with error ribbon
    # add the individual GO terms in light grey
    for _, row in data.iterrows():
        plt.plot(time, row, marker='.', markersize=6, color='grey', alpha=0.3)
    y = np.mean(data, axis=0)
    n_terms = data.shape[0]
    t_crit = stats.t.ppf(0.84, df=n_terms - 1)
    err = t_crit * np.std(data, axis=0) / np.sqrt(n_terms)
    ax.fill_between(time, y - err, y + err, color='#de2d26', alpha=0.3)
    ax.plot(time, y, marker='o', markersize=4, color='#de2d26', alpha=0.7)
    ax.set_xlabel('Time (min)', fontsize=fsize-2)
    ax.set_ylabel('Relative enrichment', fontsize=fsize-2, labelpad=0)
    ax.set_yticks([-6,-4,-2, 0])
    ax.set_xticks([500,1000,1500])
    ax.set_xticklabels([500, 1000, 1500], fontsize=fsize)
    # set tick size
    ax.tick_params(axis='both', which='major', labelsize=fsize)

def panel_D(ax):
    data_dir = os.path.join(root_dir, 'scanlag_data', 'bulk time in shx')
    x_min = 900
    x_max = 3000
    cmap = plt.get_cmap('Reds')
    v_min = -400
    v_max = 2200

    def get_normalized_value(data, v_min, v_max):
        data = data - v_min
        return data / (v_max - v_min)

    for file in os.listdir(data_dir):
        data = pd.read_csv(os.path.join(data_dir, file), header=0)
        time_point = int(file[2:file.find('Min')])
        ax.plot(data['X'], data['Y'], color=cmap(get_normalized_value(time_point, v_min, v_max)), label=time_point,
                linewidth=1)

    norm = mpl.colors.Normalize(vmin=v_min, vmax=v_max)
    mappable = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    mappable.set_array([])  # required, even if empty
    # add the colorbar
    cbar = plt.colorbar(mappable,
                        ax=ax,
                        orientation='vertical',
                        pad=0.01,
                        aspect=10,
                        shrink=0.6)
    cbar.set_label('Time in SHX(min)', fontsize=fsize-2, labelpad=5)
    cbar.set_ticks([0, 1000, 2000])
    cbar.set_ticklabels(['0', '1', '2'])
    cbar.ax.tick_params(labelsize=fsize - 2)
    cbar.ax.set_title(r'$\times10^3$', pad=1, fontsize=fsize-3, loc='left')

    ax.set_xlabel('Lag time(min)', fontsize=fsize - 2)
    ax.set_ylabel('Survival Function', fontsize=fsize, labelpad=0)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(x_min, x_max)
    ax.set_xticks([1000, 2000, 3000])
    ax.set_xticklabels([1000, 2000, ''])
    # set tick fontsize
    ax.tick_params(axis='both', which='major', labelsize=fsize - 2)
    ax.set_ylim(0.0005, 2)
###
# Build figure 1:
fsize = 10
lw=.5
marker_size=2.5
plt.close("all")
root_dir = os.path.dirname(os.path.dirname(os.getcwd()))
pf = PanelFigure(figsize=(7, 5), label_offset=(0, 0.04))
panel_pos = [
    [0.08, 0.18, 0.1, 0.75],  # A
    [0.35, 0.58, 0.55, 0.35],  # B
    [0.3, 0.09, 0.3, 0.32],  # C
    [0.7, 0.09, 0.28, 0.32],  # D
]
# panel A:
pf.add_panel(panel_pos[0], draw_func=panel_A, label="A")
# panel B:
pf.add_panel(panel_pos[1], draw_func=panel_B, label="B")
# panel C:
pf.add_panel(panel_pos[2], draw_func=panel_C, label="C")
# panel D:
pf.add_panel(panel_pos[3], draw_func=panel_D, label="D")
pf.save("figure4.svg", dpi=300, transparent=True)
plt.show()
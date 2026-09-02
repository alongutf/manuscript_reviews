import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from figure_functions import PanelFigure
import permutation_pvalues as pv
import numpy as np
import pandas as pd
import os
import base64
import io
from PIL import Image
from scipy.stats import t
import src.analysis_functions as af
import src.data_functions as df

# ------------------------------------------------------------------
# BUILD FIGURE
# ------------------------------------------------------------------
def _load_svg_image(svg_path):
    """Extract the embedded PNG from a BioRender SVG (base64 data URI)."""
    import xml.etree.ElementTree as ET
    tree = ET.parse(svg_path)
    root = tree.getroot()
    for elem in root.iter():
        href = (elem.get('{http://www.w3.org/1999/xlink}href') or
                elem.get('href', ''))
        if href.startswith('data:image/png;base64,'):
            img_bytes = base64.b64decode(href[len('data:image/png;base64,'):])
            return Image.open(io.BytesIO(img_bytes))
    return None


def _read_norm_block(sheet, label, n_reps=3, n_times=5):
    """Return (times, replicates) for a *normalised* MPN block, located by its row label.

    Some labels occur twice in the sheet — once for the raw MPN counts and once for the
    ratios normalised to each replicate's own t=0. The normalised block is the one whose
    t=0 row is exactly 1 for every replicate.
    """
    rows = np.where(sheet[0].astype(str).str.strip() == label)[0]
    for r in rows:
        block = (sheet.iloc[r + 1:r + 1 + n_times, 0:1 + n_reps]
                 .apply(pd.to_numeric, errors='coerce').dropna())
        if len(block) == 0:
            continue
        reps = block.iloc[:, 1:].to_numpy(dtype=float)
        if np.allclose(reps[0], 1.0):
            return block.iloc[:, 0].to_numpy(dtype=float), reps
    raise ValueError(f'no normalised block labelled {label!r} in the MPN sheet')


def _plot_ccdf(ax, npy_file, title, signal_color='skyblue'):
    arr = np.load(os.path.join(ev_data_dir, npy_file))
    data1 = arr[0, :];  data1 = data1[data1 > 0]
    data2 = arr[1, :];  data2 = data2[data2 > 0]
    x2 = float(np.max(data2))
    d1s = np.sort(data1)
    d2s = np.sort(data2)
    p1 = len(d1s)
    ccdf1 = 1 - np.arange(1, p1 + 1) / p1 + 1 / p1
    p2 = len(d2s)
    ccdf2 = 1 - np.arange(1, p2 + 1) / p2 + 1 / p2
    noise = d1s < x2
    ax.loglog(d1s[noise], ccdf1[noise], '.', linestyle='-',
              color='darkgray', alpha=0.7, label='spurious', markersize=3)
    ax.loglog(d1s[~noise], ccdf1[~noise], '.', linestyle='-',
              color=signal_color, label='signal', markersize=3)
    ax.loglog(d2s, ccdf2, '.', linestyle='-',
              color='black', alpha=0.5, label='scrambled', markersize=3)
    ax.set_xlim([0.1, 30])
    ax.axvline(x2, color='k', linestyle='--', alpha=0.6)
    ax.text(x2 * 1.1, 0.8, r'$\lambda_\mathrm{max}^\mathrm{scr}$',
            fontsize=fsize - 2, va='center', ha='left', color='k', alpha=0.7,
            transform=ax.get_xaxis_transform())
    ax.set_xlabel(r'$\lambda$', fontsize=fsize - 2, labelpad=0)
    ax.set_ylabel('CCDF', fontsize=fsize - 2, labelpad=0)
    ax.set_title(title, fontsize=fsize)
    # p-value rides in the legend box, under the three series keys
    pv.legend_with_p(ax, npy_file, fontsize=fsize - 2)
    ax.tick_params(labelsize=fsize - 2)


def get_data_for_plot(path, norm=True, log=False, norm_method='sum', norm_sum=1):
    # get annotated matrix from file
    amat = df.read_from_csv(path)
    # calculate the eigenvalues and plot:
    m = amat.get_filtered_matrix().m
    pcs, pcs1 = af.get_eig_dist(m, norm=norm, log=log, norm_method=norm_method, norm_sum=norm_sum)
    return pcs, pcs1, m.shape[0]


def plot_eigvals(ax, pcs, pcs1, N, x_max, y_max, n_bins, x_label=True, y_label=True):
    # plot the eigenvalue distribution of the normalized filtered matrix
    # define limits and bin number
    P = len(pcs)
    scale = 1  # scale factor for the Marchenko-Pastur distribution
    edges = np.linspace(-0.1, x_max, num=n_bins)

    # remove zeros in pcs and pcs1
    # if alpha>1 adjust the scale factor to match theoretical results
    if P / N > 1:
        scale = N / P
        pcs = pcs[pcs != 0]
        pcs1 = pcs1[pcs1 != 0]

    # first plot
    counts, bins = np.histogram(pcs, bins=edges, density=True)
    ax.plot(bins[1:], scale * counts, color='#3182bd', linewidth=.75, label='original data')
    ax.fill_between(bins[1:], scale * counts, 0, color='#9ecae1', alpha=.4)
    # second plot
    counts, bins = np.histogram(pcs1, bins=edges, density=True)
    ax.plot(bins[1:], scale * counts, color='#de2d26', linewidth=.5, label='scrambled data')
    ax.fill_between(bins[1:], scale * counts, 0, color='#fc9272', alpha=.4)
    # plot analytical Marchenko-Pastur distribution
    x = np.linspace(-0.1, x_max, 100)
    y = [af.mp_distribution(val, P / N) for val in x]
    ax.plot(x, y, color='#756bb1', linestyle='dashed', label='MP')
    # labels and limits
    if x_label:
        ax.set_xlabel("$\lambda$", fontsize=fsize, labelpad=0)
    if y_label:
        ax.set_ylabel(r"$\rho(\lambda)$", fontsize=fsize, labelpad=0)
    ax.set_ylim(0, y_max)
    ax.set_xlim(0, x_max)
    # set x_ticks with difference of 2
    ax.set_xticks(np.arange(0, (x_max // 2) * 2 + 2, 2))
    # set y_ticks with difference of 0.1
    ax.set_yticks([0.1,0.2])
    ax.legend(facecolor='white', framealpha=1, fontsize=fsize-2, loc='upper right')
    # set the font size of the ticks
    ax.tick_params(axis='both', which='major', labelsize=fsize)


def panel_A(ax):
    svg_path = os.path.join(root_dir, 'scripts', 'figures', 'figure5', 'biorender2.svg')
    img = _load_svg_image(svg_path)
    ax.imshow(img)
    ax.axis('off')


def panel_B(ax):
    """Kill curve: VapC cells vs. control, in the same format as figure 1C.

    Survival fractions are ratios spanning several decades, so the panel plots the
    geometric mean with multiplicative geometric-SD whiskers (GM/GSD to GM*GSD). An
    arithmetic mean +/- SD is unusable here: at 48 h the VapC point is 0.0069 +/- 0.0070,
    putting the lower bound below zero, which cannot be drawn on a log axis.

    GM and GSD are computed directly from the three biological replicates (sample SD of
    the logs, ddof=1). `kill curves/VapC.xlsx` holds the same curves but only as summary
    mean/SD, so it is not used here.
    """
    path = os.path.join(root_dir, 'kill curves', '20260719_VIGA24h_TolwoATC_Prep.xlsx')
    sheet = pd.read_excel(path, sheet_name='MPN', header=None)

    floor = 1e-5  # detection limit; points below it get a label only, no marker or line
    series = [
        ('CASP dilAMP (20260705)', 'Reg-Arrest', '#2166ac'),
        ('Norm. vapC dilAMP', 'VapC 24h', '#b2182b'),
    ]
    for block_label, label, color in series:
        x, reps = _read_norm_block(sheet, block_label)
        # a replicate of exactly 0 is an undetected plating, not a real zero: it has no
        # logarithm, so the time point is treated as below the detection limit
        measurable = (reps > 0).all(axis=1)
        gm = np.full(len(x), np.nan)
        gsd = np.full(len(x), np.nan)
        log_reps = np.log(reps[measurable])
        gm[measurable] = np.exp(log_reps.mean(axis=1))
        gsd[measurable] = np.exp(log_reps.std(axis=1, ddof=1))

        detected = measurable & (gm >= floor)
        # multiplicative whiskers: GM / GSD  to  GM * GSD
        yerr = np.vstack([gm[detected] * (1 - 1 / gsd[detected]),
                          gm[detected] * (gsd[detected] - 1)])
        ax.errorbar(x[detected], gm[detected], yerr=yerr, fmt='o-', markersize=4,
                    capsize=2, linewidth=1, color=color, label=label)

        censored = ~detected
        if censored.any():
            # below the detection limit: label only, so the curve ends at the last
            # measurable point instead of being drawn down onto the axis cut-off
            ax.text(x[censored].min(), floor, r'$<10^{-5}$', fontsize=fsize - 2,
                    color=color, ha='center', va='bottom')

    ax.set_yscale('log')
    #ax.set_ylim(floor * 0.6, 4)
    #ax.set_yticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1])
    ax.set_xlim(-1, 50)
    ax.set_xticks([0, 24, 48])
    ax.set_xlabel('Time in Ampicillin (h)', fontsize=fsize - 2, labelpad=0)
    ax.set_ylabel('Survival fraction', fontsize=fsize - 2, labelpad=0)
    ax.tick_params(axis='both', labelsize=fsize - 2)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(fontsize=fsize - 3, frameon=False, loc='upper right', handlelength=1.2,
              borderpad=0.1, labelspacing=0.2, handletextpad=0.4)


def panel_C(ax):
    data = pd.read_csv(os.path.join(root_dir, 'scanpy', 'umap_coordinates_vapc.csv'), index_col=0,
                       header=0)

    # scatter plot
    exp_data = data[data['batch'] == 'exp']
    t2_data = data[data['batch'] == 'T2']
    t5a_data = data[np.logical_or(data['batch'] == 'T5A', data['batch'] == 'T5B')]
    dis_data = data[data['batch'] == 'TON']
    colors = ['#4393c3', '#f4a582', '#d6604d', '#b2182b']
    ax.scatter(exp_data.UMAP_1, exp_data.UMAP_2, color=colors[0], alpha=.8, s=.5, label='Exponential')
    ax.scatter(t2_data.UMAP_1, t2_data.UMAP_2, color=colors[1], alpha=.8, s=.5, label='VapC: 2h')
    ax.scatter(t5a_data.UMAP_1, t5a_data.UMAP_2, color=colors[2], alpha=.8, s=.5, label='VapC: 5h')
    ax.scatter(dis_data.UMAP_1, dis_data.UMAP_2, color=colors[3], alpha=.8, s=.5, label='VapC: 24h')
    ax.legend(fontsize=fsize-2, loc='lower right', bbox_to_anchor=(1.2,-0.2), markerscale=4, frameon=False)
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])

def panel_D(ax):
    # UMAP panel:
    # Load data
    # project directory
    data = pd.read_csv(os.path.join(root_dir, 'scanpy', 'umap_coordinates_vapc.csv'), index_col=0,
                       header=0)
    # scatter plot
    colors = ['#8073ac', '#b2182b', '#4393c3', '#d6604d', '#92c5de']
    # color by cluster
    for i in range(max(data['cluster']) + 1):
        ax.scatter(data[data['cluster'] == i].UMAP_1, data[data['cluster'] == i].UMAP_2, color=colors[i], alpha=.6, s=.5,
                   label=f"Cluster {i + 1}")
        ax.text(data[data['cluster'] == i].UMAP_1.mean(), data[data['cluster'] == i].UMAP_2.mean(), str(i),
                fontsize=fsize, color='k', ha='center', va='center')
    #remove axis spines
    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])

def panel_E(ax):
    _plot_ccdf(ax, 'VapC_biorep_t2A_filtered.npy', 'Early VapC (2h)', signal_color=REG_COLOR)


def panel_F(ax):
    _plot_ccdf(ax, 'VapC_biorep_tONA_filtered.npy', 'Late VapC (24h)', signal_color=DIS_COLOR)

def panel_G(ax):
    # data_metrics.csv is the current metrics table (18 datasets, and the source of
    # the GMP-Cor CI); test8.csv is an older 15-dataset scramble realisation
    data = pd.read_csv(os.path.join(root_dir, 'results', 'data_metrics',
                                    'data_metrics.csv'), index_col=0)
    sample_map = [
        ('Expira_biorep_t0A_filtered.csv', 'Exponential', '#4393c3'),
        ('VapC_biorep_t2A_filtered.csv',   'VapC\n2h',   '#f4a582'),
        ('VapC_biorep_t5A_filtered.csv',   'VapC\n5h',   '#d6604d'),
        ('VapC_biorep_tONA_filtered.csv',  'VapC\n24h',  '#b2182b'),
    ]
    labels, values, colors, errs = [], [], [], []
    for fname, lbl, col in sample_map:
        row = data[data['file_name'] == fname]
        if len(row) > 0:
            labels.append(lbl)
            values.append(row['sum_denoised_ev'].iloc[0])
            colors.append(col)
            # GMP-Cor uncertainty sqrt(N)*sigma, propagated from the noise in the
            # scrambled threshold (see scripts/add_permutation_metrics.py)
            errs.append(pv.gmp_cor_ci(fname) or 0.0)

    bar_width = 0.25
    gap_between_bars = 0.4
    positions = [i * (bar_width + gap_between_bars) for i in range(len(values))]
    ax.bar(positions, values, color=colors, edgecolor='black', alpha=0.7, width=bar_width,
           yerr=errs, capsize=3, error_kw=dict(ecolor='black', elinewidth=1, capthick=1))
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, rotation=0, fontsize=fsize-2, ha='center')
    ax.set_ylabel('GMP-Cor', fontsize=fsize, labelpad=0)
    ax.tick_params(axis='both', which='major', labelsize=fsize-2)
    ax.set_ylim([0, 45])
    ax.set_xlim([positions[0] - gap_between_bars/2, positions[-1] + gap_between_bars/2])


def panel_H(ax):
    path = os.path.join(root_dir,'scripts','figures','figure5','normalizedOD_at_20h.csv')
    data = pd.read_csv(path, header=0,
                       index_col=0)

    data_dict = {}  # will hold (group_mean, group_SE) for plotting
    bio_means = {}  # stores list of biological-replicate means per group (for tests)

    for col in data.columns:
        X = data[col].dropna().to_numpy()
        means = []
        errors = []
        for i in range(int(len(X) / 2)):
            biorep = X[2 * i:2 * i + 2]
            means.append(np.mean(biorep))
            errors.append(np.std(biorep, ddof=1) / np.sqrt(2))
        bio_means[col] = means[:]  # keep biological-replicate means

        # your original rule to choose SE (biological vs technical)
        if np.std(means, ddof=1) < np.mean(errors):
            data_dict[col] = (np.mean(means), np.mean(errors))
        else:
            data_dict[col] = (np.mean(means), np.std(means, ddof=1) / np.sqrt(len(means)))

    # --- plotting ---
    labels = list(data_dict.keys())
    plot_colors = ['#2166ac', '#92c5de', '#f4a582', '#d6604d', '#b2182b']

    bar_means = [data_dict[col][0] for col in labels]
    bar_errs = [data_dict[col][1] for col in labels]
    bar_labels = ['Reg\n2h', 'Reg\n24h', 'VapC\n2h', 'VapC\n5h', 'VapC\n24h']
    bars = ax.bar(bar_labels, bar_means, yerr=bar_errs, capsize=3,
                  color=plot_colors[:len(labels)], alpha=0.8, edgecolor='black', width=0.4)

    ax.set_ylabel('Normalized OD', fontsize=fsize, labelpad=0)

    ax.set_yscale('log')
    ax.set_yticks([0.1,1])
    ax.set_yticklabels(['0.1','1'])
    ax.tick_params(axis='both', which='major', labelsize=fsize - 2)
    ax.set_title('SDS added to culture', fontsize=fsize-2, pad=0)
    ax.margins(y=0.2)  # add headroom for significance caps

    # --- significance helpers ---
    def get_pval(means, errors, n1, n2):
        tstat = (means[0] - means[1]) / np.sqrt(errors[0] ** 2 + errors[1] ** 2)
        df = n1 + n2 - 2
        return t.sf(tstat, df)

    def p_to_stars(p):
        if np.isnan(p):   return 'n.s.'
        if p < 1e-4:      return '****'
        if p < 1e-3:      return '***'
        if p < 0.01:      return '**'
        if p < 0.05:      return '*'
        return 'n.s.'

    def add_sig_between(ax, bar1, bar2, text, y_mult=1.2, cap_mult=1.1, text_mult=1.01):
        x1 = bar1.get_x() + bar1.get_width() * 0.6
        x2 = bar2.get_x() + bar2.get_width() * 0.4
        y_base = max(bar1.get_height(), bar2.get_height())
        y = y_base * y_mult
        y_cap = y * cap_mult
        ax.plot([x1, x1, x2, x2], [y, y_cap, y_cap, y], lw=1, c='black')
        ax.text((x1 + x2) / 2, y_cap * text_mult, text, ha='center', va='bottom', fontsize=fsize-2)

    # --- compute Welch t-tests between neighbor bars and annotate ---
    for i in range(len(labels) - 1):
        g1, g2 = labels[i], labels[i + 1]
        means = [data_dict[g1][0], data_dict[g2][0]]
        errors = [data_dict[g1][1], data_dict[g2][1]]

        # Welch t-test on biological-replicate means (needs n>=2 ideally)
        p = get_pval(means, errors, len(bio_means[g1]), len(bio_means[g2]))
        print(p)
        add_sig_between(ax, bars[i], bars[i + 1], p_to_stars(p), y_mult=1.3, cap_mult=1.2)



def panel_I(ax):
    # vapc lag time distribution -- vertical violins, one per condition
    conditions = ['CTRLt0', 'CTRLt1400', 'VAPCt240', 'VAPCt1400']
    labels = ['Exp', 'Reg-\nArrest', 'Early\nVapC', 'Late\nVapC']
    colors = ['#2166ac', '#9ecae1', '#fb6a4a', '#a50f15']
    plt.style.use('default')

    data = {}
    for condition in conditions:
        path = os.path.join(root_dir, 'scripts', 'figures', 'figure5', condition + '.csv')
        data[condition] = pd.read_csv(path, index_col=False, header=None).to_numpy()
        if condition == 'CTRLt0':
            t0 = np.median(data[condition])
    print(t0)
    # left to right in the listed order
    series = [data[c].flatten() - t0 for c in conditions]
    positions = list(range(1, len(conditions) + 1))

    parts = ax.violinplot(series, positions=positions, vert=True,
                          widths=0.75, showextrema=False, showmedians=False)
    for body, color in zip(parts['bodies'], colors):
        body.set_facecolor(color)
        body.set_edgecolor('k')
        body.set_linewidth(0.5)
        body.set_alpha(1)

    # median + interquartile range on top of each violin
    for x, pos in zip(series, positions):
        q1, med, q3 = np.percentile(x, [25, 50, 75])
        ax.vlines(pos, q1, q3, color='k', linewidth=1.5, zorder=3)
        ax.plot(pos, med, 'o', color='w', markersize=3,
                markeredgecolor='k', markeredgewidth=0.5, zorder=4)

    ax.set_ylabel('Lag time (min)', fontsize=fsize, labelpad=2)
    ax.set_ylim([-100, 750])
    ax.set_yticks([0, 200, 400, 600])
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=fsize - 2)
    ax.set_xlim([0.4, len(conditions) + 0.6])
    ax.tick_params(axis='y', which='major', labelsize=fsize)
    ax.tick_params(axis='x', which='major', length=0, pad=1)
    ax.spines[['top', 'right']].set_visible(False)


def panel_I_3d(ax):
    # vapc lag time distribution -- 3D waterfall histogram, one condition per depth slice
    # Load the data
    conditions = ['CTRLt0', 'CTRLt1400', 'VAPCt240', 'VAPCt1400']
    labels = ['Exp', 'Reg-Arrest', 'Early VapC', 'Late VapC']
    data = {}
    colors = ['#2166ac','#9ecae1', '#fb6a4a', '#a50f15']
    plt.style.use('default')

    for condition in conditions:
        path = os.path.join(root_dir, 'scripts', 'figures', 'figure5', condition + '.csv')
        data[condition] = pd.read_csv(path, index_col=False, header=None).to_numpy()
        if condition=='CTRLt0':
            t0 = np.min(data[condition])

    # replace the 2D axes with a 3D one at the same position
    pos = ax.get_position()
    fig = ax.figure
    ax.remove()
    ax = fig.add_axes(pos, projection='3d')

    edges = np.linspace(0, 700, 51)
    centers = 0.5 * (edges[:-1] + edges[1:])
    width = edges[1] - edges[0]
    # draw back-to-front so the near slices overlay the far ones
    for i in range(len(conditions) - 1, -1, -1):
        condition = conditions[i]
        x = data[condition].flatten() - t0
        h, _ = np.histogram(x, bins=edges, density=True)
        ax.bar(centers, h, zs=i, zdir='y', width=width,
               color=colors[i], edgecolor='k', linewidth=0.2, alpha=1,
               label=labels[i])

    ax.set_xlabel('Lag time (min)', fontsize=fsize - 2, labelpad=-4)
    ax.set_zlabel(r'Frequency ($\times10^{-2}$)', fontsize=fsize - 2, labelpad=-6)
    ax.set_xlim([0, 750])
    ax.set_ylim([-0.5, len(conditions) - 0.5])
    ax.set_xticks([200, 400, 600])
    ax.set_yticks(range(len(conditions)))
    ax.set_yticklabels(labels)
    ax.set_zticks([0, 0.01, 0.02])
    ax.set_zticklabels([0, 1, 2])
    ax.tick_params(axis='both', which='major', labelsize=fsize - 3, pad=-2)
    ax.tick_params(axis='y', pad=-3)
    ax.view_init(elev=22, azim=-58)
    ax.set_box_aspect((1.5, 1.1, 0.85), zoom=1.0)
    ax.xaxis.pane.set_alpha(0.0)
    ax.yaxis.pane.set_alpha(0.0)
    ax.zaxis.pane.set_alpha(0.0)
    ax.grid(False)

###
# Build figure 5:
fsize = 10
plt.close("all")
root_dir = os.path.dirname(os.path.dirname(os.getcwd()))
ev_data_dir = os.path.join(root_dir, 'ev_data')
REG_COLOR = 'steelblue'
DIS_COLOR = '#E07B54'
pf = PanelFigure(figsize=(7, 6), label_offset=(0, 0.03))
panel_pos = [
    [0.01, 0.7, 0.34, 0.24],  # A
    [0.42, 0.74, 0.21, 0.2],  # B
    [0.075, 0.45, 0.24, 0.2],  # C
    [0.4, 0.45, 0.24, 0.2],  # D
    [0.075, 0.08, 0.24, 0.28],  # E
    [0.4, 0.08, 0.24, 0.28],  # F
    [0.7, 0.72, 0.275, 0.22],  # G
    [0.7, 0.41, 0.275, 0.22],  # H
    [0.745, 0.08, 0.23, 0.24],  # I
]
# panel A:
pf.add_panel(panel_pos[0], draw_func=panel_A)
# panel B:
pf.add_panel(panel_pos[1], draw_func=panel_B)
# panel C:
pf.add_panel(panel_pos[2], draw_func=panel_C)
# panel D:
pf.add_panel(panel_pos[3], draw_func=panel_D)
# panel E:
pf.add_panel(panel_pos[4], draw_func=panel_E)
# panel F:
pf.add_panel(panel_pos[5], draw_func=panel_F)
# panel G:
pf.add_panel(panel_pos[6], draw_func=panel_G)
# panel H:
pf.add_panel(panel_pos[7], draw_func=panel_H)
# panel I:
pf.add_panel(panel_pos[8], draw_func=panel_I)
pf.save("figure5.pdf", dpi=300)
pf.save("figure5_preview.png", dpi=200)
plt.show()
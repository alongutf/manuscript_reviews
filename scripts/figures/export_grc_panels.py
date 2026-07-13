"""
Export selected figure panels as standalone transparent SVGs for the GRC poster.

Reuses the exact panel code from figureN.py by exec'ing each script up to (but
not including) its `pf = PanelFigure(...)` assembly block. That runs all the
panel-function definitions and their module-level globals (fsize, root_dir,
colours, precomputed medians, ...) without building/saving the full figure.
We then raise fsize by 2 and render only the requested panels.

Requested panels:
    figure2 C, D | figure3 A, B, E | figure4 B | figure5 E

All CCDF panels (figure3 A/B, figure5 E) share one figure size => same aspect.
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, REPO)     # for `import src.*`
sys.path.insert(0, HERE)     # for `from figure_functions import ...`
os.chdir(HERE)               # module root_dir is derived from os.getcwd()

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT_DIR = r'G:\Other computers\My MacBook Air\Alon\PhD\documents\GRC conference\figures'
FSIZE_BUMP = 2

# Shared aspect ratio for all CCDF panels
CCDF_SIZE = (3.2, 2.8)


def load_prefix_namespace(fig_filename, bump=True):
    """Exec a figureN.py up to its `pf = PanelFigure` line; return the namespace."""
    path = os.path.join(HERE, fig_filename)
    with open(path, 'r', encoding='utf-8') as fh:
        lines = fh.readlines()
    cut = next(i for i, ln in enumerate(lines) if ln.lstrip().startswith('pf = PanelFigure'))
    prefix = ''.join(lines[:cut])
    ns = {'__name__': '_grc_export', '__file__': path}
    exec(compile(prefix, path, 'exec'), ns)
    if bump:
        ns['fsize'] = ns['fsize'] + FSIZE_BUMP   # increase all font sizes by 2 pt
    return ns


def save(fig, name, tight_bbox=True):
    bbox = 'tight' if tight_bbox else None
    svg = os.path.join(OUT_DIR, name + '.svg')
    fig.savefig(svg, transparent=True, bbox_inches=bbox)
    fig.savefig(os.path.join(HERE, '_grc_preview_' + name + '.png'),
                dpi=200, bbox_inches=bbox)
    plt.close(fig)
    print('wrote', svg)


def render_single(ns, panel, name, figsize, tight_bbox=True, adjust=None, post=None):
    fig, ax = plt.subplots(figsize=figsize)
    ns[panel](ax)
    if post is not None:
        post(ax)
    if adjust is not None:
        fig.subplots_adjust(**adjust)
    elif not tight_bbox:
        fig.tight_layout()
    save(fig, name, tight_bbox=tight_bbox)


def render_grid(ns, panel, name, figsize):
    fig, axes = plt.subplots(2, 1, squeeze=False, figsize=figsize)
    ns[panel](axes)
    fig.tight_layout()
    save(fig, name)


def _wrap_two_rows(text):
    """Split a label across two rows at the space that best balances line length."""
    words = text.split(' ')
    if len(words) < 2:
        return text
    best_i, best_diff = 1, None
    for i in range(1, len(words)):
        diff = abs(len(' '.join(words[:i])) - len(' '.join(words[i:])))
        if best_diff is None or diff < best_diff:
            best_diff, best_i = diff, i
    return ' '.join(words[:best_i]) + '\n' + ' '.join(words[best_i:])


def _tidy_go_label(text):
    """Shorten a GO-term x-tick label for the poster panel 4B."""
    # drop the trailing gene count, e.g. " (n=34)"
    text = re.sub(r'\s*\(n=\d+\)', '', text)
    # rule 1: remove 'bacterial-type'
    text = re.sub(r'\bbacterial-type\s*', '', text)
    # rules 2 & 3: abbreviations
    text = re.sub(r'type III secretion system', 'T3SS', text, flags=re.IGNORECASE)
    text = re.sub(r'proton motive force', 'PMF', text, flags=re.IGNORECASE)
    text = re.sub(r'\s+', ' ', text).strip()
    # rule 4: if still > 5 words (a dash counts as 2 words), break into two rows
    if len(re.split(r'[\s-]+', text)) > 5:
        text = _wrap_two_rows(text)
    return text


def strip_gene_counts(ax):
    """Apply the panel-4B label rules to each x-tick label, keeping font/rotation."""
    old = ax.get_xticklabels()
    if not old:
        return
    fs = old[0].get_fontsize()
    labels = [_tidy_go_label(t.get_text()) for t in old]
    ax.set_xticklabels(labels, fontsize=fs, rotation=45, ha='right')


# ------------------------------------------------------------------
ns2 = load_prefix_namespace('figure2.py')
render_grid(ns2, 'panel_C', 'figure2_panelC', (1.9, 3.6))
# Panel D: keep the original font size (no bump) and the original panel size so
# the inset CCDF legend fits as it does in figure2.py.
ns2_orig = load_prefix_namespace('figure2.py', bump=False)
render_grid(ns2_orig, 'panel_D', 'figure2_panelD', (3.08, 3.38))

ns3 = load_prefix_namespace('figure3.py')
render_single(ns3, 'panel_A', 'figure3_panelA', CCDF_SIZE)
render_single(ns3, 'panel_B', 'figure3_panelB', CCDF_SIZE)
render_single(ns3, 'panel_C', 'figure3_panelC', CCDF_SIZE)
render_single(ns3, 'panel_D', 'figure3_panelD', CCDF_SIZE)
render_single(ns3, 'panel_E', 'figure3_panelE', (2.4, 3.2))

ns4 = load_prefix_namespace('figure4.py')
# inset_axes (mpl_toolkits) is incompatible with bbox_inches='tight'; the panel
# is a wide/short bar chart with long rotated x labels -> explicit margins.
render_single(ns4, 'panel_B', 'figure4_panelB', (7.6, 5.2), tight_bbox=False,
              adjust=dict(left=0.16, right=0.97, top=0.95, bottom=0.5),
              post=strip_gene_counts)

ns5 = load_prefix_namespace('figure5.py')
render_single(ns5, 'panel_D', 'figure5_panelD', CCDF_SIZE)
render_single(ns5, 'panel_E', 'figure5_panelE', CCDF_SIZE)

print('done')

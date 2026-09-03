"""Cluster-annotated marker-gene heatmaps for the scRNA-seq experiments.

Used by ``scripts/heatmaps.ipynb``.  The figure is a genes x cells z-score
heatmap with three layers of annotation:

* a coloured bar above the heatmap naming each Leiden cluster by its biological
  identity, not the bare Leiden number,
* a colour strip on the y axis marking, for every gene row, the cluster the
  gene was picked as a marker for (same colours as the cluster bar),
* curated gene-group labels on the left, each bracketed to its block of rows.

Gene rows are sectioned by the cluster they mark, so the strip is a run of solid
blocks.  A curated group whose genes mark more than one cluster is split into one
labelled block per section; a block holding a single gene shows only the gene
name, with no group title.
"""

import textwrap

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from scipy.cluster.hierarchy import linkage, leaves_list

# Okabe-Ito blue / vermillion / bluish-green.  Validated against the categorical
# colour checks: lightness band, chroma floor, all-pairs CVD separation
# (min protan/deutan dE >= 11), normal-vision separation, contrast vs. white.
CLUSTER_COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00"]


def select_markers(adata, cats, topN=100, p_cut=0.05, groupby="leiden",
                   layer="counts", use_raw=False):
    """Up to ``topN`` upregulated markers per cluster, deduplicated.

    Only genes that are significantly *up* in a cluster are kept, so a cluster
    contributes fewer than ``topN`` rows whenever it has fewer than ``topN``
    significant upregulated genes.  A gene that makes more than one cluster's
    list is kept once, assigned to the cluster where its score is strongest, so
    the total is at most ``topN * len(cats)`` and usually less.
    """
    # reuse a rank_genes_groups run already stashed on adata.uns instead of
    # recomputing the (slow) Wilcoxon test on every call
    if "rank_genes_groups" not in adata.uns:
        sc.tl.rank_genes_groups(adata, groupby=groupby, method="wilcoxon",
                                layer=layer, use_raw=use_raw)

    df = sc.get.rank_genes_groups_df(adata, None).copy()
    # a gene with zero counts in one group can blow up the Wilcoxon score to +/-inf;
    # such rows can't be ranked meaningfully, so they are dropped rather than kept
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["pvals_adj", "scores"])

    picked = []
    for grp in cats:
        # "scores" > 0 keeps only genes that are UP in this cluster, not just
        # differentially expressed in either direction
        sub = df[(df["group"] == grp) & (df["pvals_adj"] < p_cut) & (df["scores"] > 0)]
        pick = sub.nlargest(topN, "scores").copy()
        pick["chosen_group"] = grp
        picked.append(pick)

    sel = pd.concat(picked, ignore_index=False)
    # resolve genes shared between clusters' top-N lists: keep the single row where
    # the gene's score magnitude is largest, so each gene is assigned once, to
    # whichever cluster it marks most strongly
    sel["abs_score"] = sel["scores"].abs()
    sel = sel.loc[sel.groupby("names")["abs_score"].idxmax()].sort_index()
    return sel


def build_matrix(adata, sel, cats, gene_groups, groupby="leiden", layer="counts",
                 use_raw=False, metric="correlation", method="average",
                 order_cells_within_cluster="umap"):
    """Return the genes x cells z-scored matrix plus its ordering metadata."""
    A = adata.raw.to_adata() if use_raw else adata

    Xfull = A.layers[layer] if layer is not None else A.X
    if sp.issparse(Xfull):
        Xfull = Xfull.toarray()

    # restrict to the marker genes selected upstream, in adata's own var order
    genes_present = [g for g in sel["names"].tolist() if g in A.var_names]
    Xsel = Xfull[:, A.var_names.get_indexer(genes_present)]

    # per-gene z-score (cells x genes here; transposed to genes x cells below).
    # a gene with zero variance across cells would divide by zero, so its sd is
    # pinned to 1 -- the z-scored column is then all zeros instead of NaN/inf
    mu, sd = np.nanmean(Xsel, axis=0), np.nanstd(Xsel, axis=0)
    sd[sd == 0] = 1.0
    Xz = (Xsel - mu) / sd

    # ---- order cells, grouped by cluster -------------------------------------
    obs_grp = A.obs[groupby].astype(str).values
    cells_order, cluster_boundaries, start = [], [], 0
    for grp in cats:
        idx = np.where(obs_grp == grp)[0]
        if len(idx) == 0:
            continue
        if order_cells_within_cluster == "umap" and "X_umap" in A.obsm:
            # cheap proxy for within-cluster structure: sort along UMAP-1 instead
            # of clustering (which would be expensive at cell count)
            sub = idx[np.argsort(A.obsm["X_umap"][idx, 0])]
        elif order_cells_within_cluster == "cluster":
            sub = idx[leaves_list(linkage(Xz[idx, :], method=method, metric=metric))]
        else:
            sub = np.sort(idx)
        cells_order.extend(sub.tolist())
        cluster_boundaries.append((grp, start, start + len(sub)))
        start += len(sub)

    Xz_cells = Xz[cells_order, :]

    # ---- order genes: sectioned by the cluster they mark ---------------------
    # Within a section genes are clustered hierarchically and the curated groups
    # are pulled together into blocks.  A curated group whose genes mark more
    # than one cluster therefore appears once per section, as separate blocks.
    gene2grp = {g: lab for lab, genes in gene_groups.items() for g in genes}
    gene2cluster = dict(zip(sel["names"], sel["chosen_group"].astype(str)))

    def order_section(section_genes):
        """Hierarchical order of one section, with curated groups blockified."""
        cols = [genes_present.index(g) for g in section_genes]
        # linkage on 2 or fewer rows is degenerate/uninformative, so just keep them
        # in their existing order instead of clustering
        if len(cols) > 2:
            Z = linkage(Xz_cells[:, cols].T, method=method, metric=metric)
            base = [section_genes[i] for i in leaves_list(Z)]
        else:
            base = list(section_genes)

        in_section = set(section_genes)
        order, ranges, placed, placed_groups = [], [], set(), set()
        for g in base:
            if g in placed:
                continue
            grp = gene2grp.get(g)
            if grp is not None and grp not in placed_groups:
                members = [x for x in gene_groups[grp]
                           if x in in_section and x not in placed]
                if members:
                    idx = [genes_present.index(x) for x in members]
                    # same >2 guard as above: clustering needs at least 3 rows
                    if len(idx) > 2:
                        Zg = linkage(Xz_cells[:, idx].T, method=method, metric=metric)
                        idx = np.array(idx)[leaves_list(Zg)].tolist()
                    genes_local = [genes_present[j] for j in idx]
                    s = len(order)
                    order.extend(genes_local)
                    ranges.append((grp, s, len(order)))
                    placed.update(genes_local)
                    placed_groups.add(grp)
                continue
            order.append(g)
            placed.add(g)
        return order, ranges

    final_order, group_ranges, gene_sections = [], [], []
    for grp in cats:
        section = [g for g in genes_present if gene2cluster.get(g) == grp]
        if not section:
            continue
        order, ranges = order_section(section)
        offset = len(final_order)
        final_order.extend(order)
        group_ranges.extend((lab, s + offset, e + offset) for lab, s, e in ranges)
        gene_sections.append((grp, offset, len(final_order)))

    # genes whose cluster is unknown (should not happen) go last
    leftover = [g for g in genes_present if g not in set(final_order)]
    if leftover:
        final_order.extend(leftover)

    H = Xz_cells[:, [genes_present.index(g) for g in final_order]].T
    return H, final_order, group_ranges, cluster_boundaries, gene_sections


def _wrap_genes(genes, width=30):
    """Comma-join gene names and word-wrap to ``width`` chars per line for a label."""
    return textwrap.wrap(", ".join(genes), width=width) or [""]


def _place_labels(blocks, n_genes, row_h_pt, line_h_pt, title_h_pt):
    """Greedy top-to-bottom placement so stacked labels never overlap.

    ``blocks`` is a list of (label, start_row, end_row, n_lines, has_title).
    Returns (label, start, end, text_centre_row, half_height_rows) per block.
    """
    placed, cursor = [], -0.5
    for label, s, e, n_lines, has_title in blocks:
        h_pt = (title_h_pt if has_title else 0.0) + n_lines * line_h_pt
        half = (h_pt / 2.0) / row_h_pt
        mid = max((s + e - 1) / 2.0, cursor + half)
        # 1.2x line height as a fixed gap between consecutive stacked labels
        cursor = mid + half + (1.2 * line_h_pt / row_h_pt)
        placed.append([label, s, e, mid, half])

    # if the stack ran past the bottom, slide the whole stack up rigidly
    if placed:
        overflow = (placed[-1][3] + placed[-1][4]) - (n_genes - 0.5)
        if overflow > 0:
            for blk in placed:
                blk[3] -= overflow
    return [tuple(b) for b in placed]


def plot_marker_heatmap(adata, cats, cluster_labels, gene_groups,
                        outfile, topN=100, p_cut=0.05,
                        groupby="leiden", layer="counts", use_raw=False,
                        metric="correlation", method="average",
                        order_cells_within_cluster="umap",
                        vmin=-2, vmax=2, cmap="viridis", fontsize=10,
                        label_width=30):
    """Draw and save the annotated marker heatmap.

    Selects up to ``topN`` marker genes per cluster in ``cats`` (see
    ``select_markers``), builds the genes x cells z-score matrix (see
    ``build_matrix``), and renders it with the cluster-identity bar, "marker of"
    colour strip and curated group labels described in the module docstring.
    ``vmin``/``vmax`` clip the z-score colour scale (default +/-2 std). The figure
    is written to ``outfile``. Returns (fig, info dict), where info carries the
    row/column ordering and boundary metadata used to draw the panel, for anyone
    who wants to reuse or annotate it further.
    """
    sel = select_markers(adata, cats, topN=topN, p_cut=p_cut,
                         groupby=groupby, layer=layer, use_raw=use_raw)
    H, gene_order, group_ranges, cluster_boundaries, gene_sections = build_matrix(
        adata, sel, cats, gene_groups, groupby=groupby, layer=layer,
        use_raw=use_raw, metric=metric, method=method,
        order_cells_within_cluster=order_cells_within_cluster)

    n_genes, n_cells = H.shape
    gene2cluster = dict(zip(sel["names"], sel["chosen_group"].astype(str)))
    color_of = {grp: CLUSTER_COLORS[i % len(CLUSTER_COLORS)]
                for i, grp in enumerate(cats)}

    # ---- canvas -------------------------------------------------------------
    # figure size scales with the row/column counts (with a floor) so labels and
    # cells stay legible regardless of how many genes/cells this call ends up with
    fig_h = max(10, 0.012 * n_genes + 4)
    fig_w = max(9, 0.004 * n_cells + 4)
    fig = plt.figure(figsize=(fig_w, fig_h))
    # 2x4 grid: (label col, colour-strip col, heatmap col, colourbar col) x
    # (top identity-bar row, main row); ratios were tuned by eye for this layout
    gs = fig.add_gridspec(
        2, 4, width_ratios=[2.5, 0.16, 7.5, 0.32], height_ratios=[0.6, 10],
        wspace=0.03, hspace=0.02, left=0.02, right=0.93, top=0.95, bottom=0.05)

    ax = fig.add_subplot(gs[1, 2])                   # heatmap
    ax_top = fig.add_subplot(gs[0, 2], sharex=ax)    # cluster identity bar
    ax_strip = fig.add_subplot(gs[1, 1], sharey=ax)  # "marker of" colour strip
    ax_lab = fig.add_subplot(gs[1, 0], sharey=ax)    # gene-group labels
    cax = fig.add_subplot(gs[1, 3])

    im = ax.imshow(H, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax,
                   interpolation="nearest", origin="upper",
                   extent=(-0.5, n_cells - 0.5, n_genes - 0.5, -0.5))
    # rasterize the (potentially huge) pixel grid so the saved SVG stays a
    # manageable size; vector elements (labels, brackets) are unaffected
    im.set_rasterized(True)
    ax.set_xticks([])
    ax.set_yticks([])
    for side in ax.spines.values():
        side.set_visible(False)

    # ---- cluster identity bar above the heatmap -----------------------------
    ax_top.set_ylim(0, 1)
    for grp, s, e in cluster_boundaries:
        ax_top.add_patch(plt.Rectangle((s - 0.5, 0.0), e - s, 0.42,
                                       facecolor=color_of[grp], edgecolor="none"))
        ax_top.text((s + e - 1) / 2.0, 0.52,
                    f"{grp} - {cluster_labels.get(grp, grp)}",
                    ha="center", va="bottom", color="black",
                    fontsize=fontsize + 3, fontweight="bold")
        if e < n_cells:
            ax.vlines([e - 0.5], -0.5, n_genes - 0.5, linewidth=1.2, color="white")
    ax_top.set_xlim(-0.5, n_cells - 0.5)
    ax_top.axis("off")

    # ---- "marker of" colour strip on the y axis -----------------------------
    # one categorical colour per gene, encoded as its index into cats; genes with
    # no recognised cluster (should not happen) fall back to grey via set_bad
    strip = np.array([[cats.index(gene2cluster[g]) if gene2cluster.get(g) in cats
                       else np.nan] for g in gene_order], dtype=float)
    cmap_strip = ListedColormap([color_of[c] for c in cats])
    cmap_strip.set_bad("#dddddd")
    # integer-valued bins centred on each cluster index, so imshow maps each
    # discrete category to its own solid colour rather than interpolating
    norm = BoundaryNorm(np.arange(-0.5, len(cats) + 0.5), len(cats))
    ax_strip.imshow(np.ma.masked_invalid(strip), aspect="auto", cmap=cmap_strip,
                    norm=norm, interpolation="nearest", origin="upper",
                    extent=(0, 1, n_genes - 0.5, -0.5))
    ax_strip.set_xticks([])
    ax_strip.set_yticks([])
    for side in ax_strip.spines.values():
        side.set_visible(False)
    ax_strip.set_title("marker\nof", fontsize=fontsize - 1, pad=3,
                       linespacing=1.1)

    # ---- curated gene-group labels ------------------------------------------
    ax_lab.set_xlim(0, 1)
    ax_lab.set_ylim(n_genes - 0.5, -0.5)
    ax_lab.axis("off")

    # label geometry is worked out in points (72 pt = 1 inch) so label heights can
    # be compared directly against the heatmap's per-row height in points
    line_h = (fontsize - 1) * 1.35
    title_h = (fontsize + 1) * 1.45
    row_h = (fig_h * gs[1, 2].get_position(fig).height) * 72 / n_genes

    # a single-gene block carries no group title -- the gene name says it all
    blocks = [(label, s, e, _wrap_genes(gene_order[s:e], label_width), e - s > 1)
              for label, s, e in sorted(group_ranges, key=lambda t: t[1])]
    placed = _place_labels([(l, s, e, len(w), t) for l, s, e, w, t in blocks],
                           n_genes, row_h, line_h, title_h)

    bracket_x = 0.985
    for (label, s, e, mid, half), (_, _, _, wrapped, has_title) in zip(placed, blocks):
        ax.hlines([s - 0.5, e - 0.5], -0.5, n_cells - 0.5,
                  linewidth=0.5, color="black")
        # bracket spanning the block of rows + one leader line to the text
        ax_lab.plot([bracket_x, bracket_x], [s - 0.5, e - 0.5],
                    color="black", linewidth=1.4, clip_on=False,
                    solid_capstyle="butt")
        ax_lab.plot([bracket_x - 0.035, bracket_x], [mid, (s + e - 1) / 2.0],
                    color="black", linewidth=0.8, clip_on=False)
        if has_title:
            ax_lab.text(bracket_x - 0.06, mid - half, label, ha="right", va="top",
                        fontsize=fontsize + 1, fontweight="bold", color="black")
        ax_lab.text(bracket_x - 0.06,
                    mid - half + (title_h / row_h if has_title else 0.0),
                    "\n".join(wrapped), ha="right", va="top",
                    fontsize=fontsize - 1, color="#333333", style="italic",
                    linespacing=1.35)

    # ---- separators between the gene sections -------------------------------
    for _, _, e in gene_sections[:-1]:
        ax.hlines([e - 0.5], -0.5, n_cells - 0.5, linewidth=1.2, color="white")

    # ---- colourbar and axis label -------------------------------------------
    cbar = fig.colorbar(im, cax=cax)
    cbar.ax.tick_params(labelsize=fontsize + 2)
    cbar.ax.set_ylabel("Expression z-score", rotation=90, labelpad=14,
                       fontsize=fontsize + 4)
    cbar.solids.set_rasterized(True)

    ax.set_xlabel(f"Single cells, grouped by {groupby} cluster",
                  fontsize=fontsize + 5, labelpad=12)

    fig.savefig(outfile, bbox_inches="tight", dpi=300)
    info = dict(n_genes=n_genes, n_cells=n_cells, gene_order=gene_order,
                group_ranges=group_ranges, cluster_boundaries=cluster_boundaries,
                gene_sections=gene_sections, markers=sel)
    return fig, info

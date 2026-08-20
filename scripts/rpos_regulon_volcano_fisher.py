"""Volcano plots + Fisher's exact enrichment for the RegulonDB sigma-38 (RpoS) sigmulon.

For every DESeq2 contrast:
  * a volcano plot (log2FoldChange vs -log10 padj) with sigma-38 genes highlighted;
  * Fisher's exact test asking whether sigma-38 genes are over-represented among the
    up-regulated lobe (LFC > 1, padj < alpha), and separately among the down-regulated
    lobe (LFC < -1, padj < alpha), relative to all other tested genes. alpha is 0.05
    except for aggregated_sc, which is held to 0.01 -- see ALPHA_BY_FOLDER in
    rpos_regulon_deseq.py.

The 2x2 table for each lobe is

                      in lobe    not in lobe
    sigma-38 gene        a            b
    other gene           c            d

with the universe restricted to genes DESeq2 actually tested in that contrast
(padj not NA). One-sided p-values are reported in both directions (enrichment and
depletion) alongside the two-sided value.

Each contrast/lobe is a separate analysis testing a single gene set, so the reported
p-value needs no within-analysis multiple-testing correction -- the same convention
goatools uses in src/bulk_functions.py, where BH runs across the GO terms of one
study and never across studies. A BH column across all contrast x lobe tests is
included as a conservative sensitivity check only.

Regulon list: metadata/regulondb_sigma38_regulon.txt, see fetch_regulondb_sigma38.py.

Outputs:
    results/deseq_results/figures/rpoS_volcano_all_contrasts.{svg,png}
    results/deseq_results/rpoS_regulon_fisher.csv

Run from the repo root or from scripts/:
    python scripts/rpos_regulon_volcano_fisher.py
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy import stats

from rpos_regulon_deseq import (DESEQ, REGULON, alpha_for, contrast_files,
                                load_regulon)

LFC_CUT = 1.0
FIGDIR = os.path.join(DESEQ, "figures")

SIG_COLOR = "#d62728"     # sigma-38 genes
BG_COLOR = "#c8c8c8"      # everything else


def bh(p):
    """Benjamini-Hochberg adjusted p-values (scipy in this env is too old to ship it)."""
    p = np.asarray(p, float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def short_name(path):
    tag = os.path.basename(path)
    tag = tag.replace("deseq2_results_", "").replace(".csv", "")
    folder = os.path.basename(os.path.dirname(path))
    label = {"from counts": "counts", "aggregated_sc": "agg_sc"}.get(folder, folder)
    return "%s: %s" % (label, tag)


def fisher_for_lobe(in_regulon, lobe):
    """2x2 Fisher for one lobe. Returns dict of counts and p-values."""
    a = int((in_regulon & lobe).sum())
    b = int((in_regulon & ~lobe).sum())
    c = int((~in_regulon & lobe).sum())
    dd = int((~in_regulon & ~lobe).sum())
    table = [[a, b], [c, dd]]
    odds, p_two = stats.fisher_exact(table, alternative="two-sided")
    _, p_greater = stats.fisher_exact(table, alternative="greater")
    _, p_less = stats.fisher_exact(table, alternative="less")
    frac_reg = a / (a + b) if (a + b) else np.nan
    frac_bg = c / (c + dd) if (c + dd) else np.nan
    # fraction of the lobe that is sigma-38, vs the regulon's share of the universe
    share_lobe = a / (a + c) if (a + c) else np.nan
    share_universe = (a + b) / (a + b + c + dd)
    return dict(
        regulon_in_lobe=a, regulon_out=b, other_in_lobe=c, other_out=dd,
        lobe_size=a + c,
        pct_regulon_in_lobe=round(100 * frac_reg, 2) if np.isfinite(frac_reg) else np.nan,
        pct_other_in_lobe=round(100 * frac_bg, 2) if np.isfinite(frac_bg) else np.nan,
        pct_of_lobe_that_is_regulon=round(100 * share_lobe, 2) if np.isfinite(share_lobe) else np.nan,
        pct_of_universe_that_is_regulon=round(100 * share_universe, 2),
        odds_ratio=round(odds, 3) if np.isfinite(odds) else np.inf,
        p_two_sided=p_two, p_greater=p_greater, p_less=p_less,
    )


def main():
    regulon, header = load_regulon(REGULON)
    print("RegulonDB release %s, %d sigma-38 genes"
          % (header.get("RegulonDB release", "?"), len(regulon)))
    os.makedirs(FIGDIR, exist_ok=True)

    files = contrast_files()
    ncol = 4
    nrow = int(np.ceil(len(files) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.0 * ncol, 3.4 * nrow))
    axes = axes.ravel()

    rows = []
    for ax, path in zip(axes, files):
        tag = os.path.relpath(path, DESEQ).replace(os.sep, "/")
        d = pd.read_csv(path, index_col=0)
        d = d[d.padj.notna() & d.log2FoldChange.notna()]

        alpha = alpha_for(path)
        in_regulon = d.index.isin(regulon)
        up = (d.log2FoldChange > LFC_CUT) & (d.padj < alpha)
        down = (d.log2FoldChange < -LFC_CUT) & (d.padj < alpha)

        for lobe_name, lobe in (("up", up), ("down", down)):
            r = fisher_for_lobe(pd.Series(in_regulon, index=d.index), lobe)
            r.update(contrast=tag, lobe=lobe_name, alpha=alpha,
                     genes_tested=len(d), regulon_tested=int(in_regulon.sum()))
            rows.append(r)

        # ---- volcano ----
        # padj is exactly 0 for a few genes; floor it at the smallest nonzero
        # padj in this contrast so the y axis stays honest instead of hitting 300.
        pad = d.padj.values.copy()
        nz = pad[pad > 0]
        floor = nz.min() if nz.size else 1e-300
        y = -np.log10(np.where(pad > 0, pad, floor))

        ax.scatter(d.log2FoldChange[~in_regulon], y[~in_regulon], s=4,
                   c=BG_COLOR, alpha=0.5, linewidths=0, rasterized=True)
        ax.scatter(d.log2FoldChange[in_regulon], y[in_regulon], s=9,
                   c=SIG_COLOR, alpha=0.85, linewidths=0)
        ax.axhline(-np.log10(alpha), color="k", lw=0.6, ls="--", alpha=0.5)
        ax.axvline(LFC_CUT, color="k", lw=0.6, ls="--", alpha=0.5)
        ax.axvline(-LFC_CUT, color="k", lw=0.6, ls="--", alpha=0.5)

        def annot(res):
            """Report the one-sided p in the direction the odds ratio actually points."""
            if not np.isfinite(res["odds_ratio"]) or res["lobe_size"] == 0:
                return "n/a"
            if res["odds_ratio"] >= 1:
                return "OR=%.1f\nenrich p=%.0e" % (res["odds_ratio"], res["p_greater"])
            return "OR=%.1f\ndeplete p=%.0e" % (res["odds_ratio"], res["p_less"])

        up_row, down_row = rows[-2], rows[-1]
        ax.set_title("%s   (padj < %g)" % (short_name(path), alpha), fontsize=9)
        ax.text(0.97, 0.96, "up lobe\n" + annot(up_row), transform=ax.transAxes,
                ha="right", va="top", fontsize=7)
        ax.text(0.03, 0.96, "down lobe\n" + annot(down_row), transform=ax.transAxes,
                ha="left", va="top", fontsize=7)
        ax.set_xlabel("log2 fold change", fontsize=8)
        ax.set_ylabel("-log10 padj", fontsize=8)
        ax.tick_params(labelsize=7)

    for ax in axes[len(files):]:
        ax.axis("off")

    handles = [Line2D([], [], marker="o", ls="", color=SIG_COLOR,
                      label="sigma-38 sigmulon (RegulonDB)"),
               Line2D([], [], marker="o", ls="", color=BG_COLOR, label="other genes")]
    fig.legend(handles=handles, loc="lower center", fontsize=9, ncol=2,
               bbox_to_anchor=(0.5, 0.005), frameon=False)
    fig.suptitle("RpoS (sigma-38) sigmulon on DESeq2 volcano plots  "
                 "[lobes: |LFC| > %g and padj below the per-panel cutoff; "
                 "Fisher exact vs all other tested genes]" % LFC_CUT, fontsize=10)
    fig.tight_layout(rect=[0, 0.035, 1, 0.965])
    for ext in ("svg", "png"):
        out = os.path.join(FIGDIR, "rpoS_volcano_all_contrasts." + ext)
        fig.savefig(out, dpi=200)
        print("wrote", out)
    plt.close(fig)

    F = pd.DataFrame(rows)
    # Each contrast/lobe is its own analysis testing a single gene set, so the
    # reported p-values need no within-analysis correction -- this mirrors the GO
    # convention in src/bulk_functions.py, where goatools BH-corrects across the GO
    # terms of one study and never across studies. The column below applies BH
    # across all contrast x lobe tests as a conservative sensitivity check only;
    # it is NOT the headline statistic. See the doc for why.
    F["q_across_all_tests"] = bh(F.p_two_sided.values)
    F = F[[
        "contrast", "lobe", "alpha", "genes_tested", "regulon_tested", "lobe_size",
        "regulon_in_lobe", "regulon_out", "other_in_lobe", "other_out",
        "pct_regulon_in_lobe", "pct_other_in_lobe",
        "pct_of_lobe_that_is_regulon", "pct_of_universe_that_is_regulon",
        "odds_ratio", "p_greater", "p_less", "p_two_sided", "q_across_all_tests"]]
    out = os.path.join(DESEQ, "rpoS_regulon_fisher.csv")
    F.to_csv(out, index=False)
    print("wrote", out)

    show = F[["contrast", "lobe", "alpha", "lobe_size", "regulon_in_lobe",
              "pct_of_lobe_that_is_regulon", "pct_of_universe_that_is_regulon",
              "odds_ratio", "p_greater", "p_less", "q_across_all_tests"]]
    with pd.option_context("display.width", 250):
        print(show.to_string(index=False))


if __name__ == "__main__":
    main()

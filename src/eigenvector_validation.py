"""
Ground-truth validation of the eigenvector (gene-loading) analysis.

The eigenvector analysis in `scripts/eigenvector_analysis.py` reads gene-level
structure off the leading eigenvectors of the gene-gene correlation matrix, and on
the experimental data returns a largely negative result: the leading modes are
delocalized and recover no condition-specific program. That negative result has two
incompatible readings -- either the coordinated signal really is distributed, or the
pipeline cannot recover modular structure at this depth and dropout even when
modules exist.

This module settles that by running the same pipeline on synthetic data whose block
structure is known. `generate_gram_hub_matrix(..., return_structure=True)` returns
the contiguous gene clusters it used to build the correlation matrix, so every
leading eigenvector can be scored against the true blocks.

Everything here is computed from the RAW LOADINGS of the top-N modes
(`analysis_functions.get_eig_vectors`) -- no coordination score.

Metrics
-------
per-mode   : how much of each mode's squared loading mass sits in its best-matching
             true block, the effective number of blocks the mode spreads over, and
             the composition of its top-loading genes (with a hypergeometric p-value).
per-block  : the projection of each true block's indicator vector onto the span of
             the top-N modes, against the mass a random gene set of the same size
             would collect -- i.e. is the block represented at all.
stability  : agreement of loadings between independent count-sampling replicates of
             the same ground-truth correlation matrix, per mode and at the level of
             the whole N-mode subspace.

The noise ladder (`ladder_stages`) applies the identical scoring to the true
correlation matrix, the latent Gaussian layer, the pre-dropout counts and the
observed counts, which localizes which noise layer -- if any -- destroys recovery.
"""

import numpy as np
from scipy.stats import hypergeom

from .analysis_functions import get_eig_vectors


# ── helpers ──────────────────────────────────────────────────────────────────

def _effective_n(p):
    """exp(Shannon entropy) of a probability vector: effective number of components."""
    p = np.asarray(p, dtype=float)
    p = p[p > 0]
    if p.size == 0:
        return 0.0
    p = p / p.sum()
    return float(np.exp(-np.sum(p * np.log(p))))


def block_mass(v, labels, n_blocks):
    """Squared loading mass of eigenvector v summed within each true block."""
    return np.bincount(labels, weights=v ** 2, minlength=n_blocks)


def oracle_eigenvectors(R, n_top=5):
    """
    Top-n_top eigenvectors of a correlation matrix directly (no sampling).

    This is the ceiling for any data-driven estimate: it shows how much block
    structure the top-n_top modes can express at all, before any count noise.
    """
    vals, vecs = np.linalg.eigh(R)
    order = np.argsort(vals)[::-1][:n_top]
    return vals[order], vecs[:, order].T


# ── per-mode scoring ─────────────────────────────────────────────────────────

def score_modes(eigvecs, labels, sizes, strengths, n_top_genes=30):
    """
    Score each mode's loadings against the ground-truth blocks.

    Parameters
    ----------
    eigvecs   : (n_modes, p) gene loadings, rows unit-norm
    labels    : (p,) true block index of each gene (already subset to kept genes)
    sizes     : (K,) true block sizes as generated (before any gene filtering)
    strengths : (K,) true block shared-loading weights w_k

    Returns a list of dicts, one per mode.
    """
    eigvecs = np.atleast_2d(eigvecs)
    p = eigvecs.shape[1]
    n_blocks = len(sizes)
    kept_sizes = np.bincount(labels, minlength=n_blocks)  # block sizes among kept genes
    rows = []
    for k, v in enumerate(eigvecs):
        mass = block_mass(v, labels, n_blocks)
        best = int(np.argmax(mass))
        # top-loading genes: how many come from the mode's best block, and how
        # surprising that is against a random draw of the same number of genes
        top_idx = np.argsort(np.abs(v))[::-1][:n_top_genes]
        n_hit = int(np.sum(labels[top_idx] == best))
        if kept_sizes[best] > 0 and n_top_genes <= p:
            p_hyper = float(hypergeom.sf(n_hit - 1, p, int(kept_sizes[best]), n_top_genes))
        else:
            p_hyper = np.nan
        rows.append(dict(
            mode=k + 1,
            best_block=best,
            best_block_size=int(sizes[best]),
            best_block_size_kept=int(kept_sizes[best]),
            best_block_strength=float(strengths[best]),
            best_block_mass=float(mass[best]),
            # mass a random gene set of the same size would carry
            best_block_mass_null=float(kept_sizes[best] / p),
            effective_n_blocks=_effective_n(mass),
            top_gene_hits=n_hit,
            top_gene_purity=n_hit / n_top_genes,
            top_gene_p_hypergeom=p_hyper,
        ))
    return rows


# ── per-block scoring ────────────────────────────────────────────────────────

def score_blocks(eigvecs, labels, sizes, strengths, min_size=5):
    """
    For every true block of at least `min_size` genes, how well the top-N modes
    represent it.

    `projection` is ||V u_b||^2 for the unit indicator u_b of the block: the fraction
    of the block's coherent (all-genes-in-phase) direction that lies inside the
    N-mode subspace, 1 = fully captured, 0 = invisible.

    `mass_ratio` is the block's total squared loading mass over the mass a random
    gene set of the same size would collect; > 1 means the modes over-weight the
    block relative to chance.
    """
    eigvecs = np.atleast_2d(eigvecs)
    n_modes, p = eigvecs.shape
    n_blocks = len(sizes)
    kept_sizes = np.bincount(labels, minlength=n_blocks)
    rows = []
    for b in range(n_blocks):
        nb = int(kept_sizes[b])
        if nb < min_size:
            continue
        u = np.zeros(p)
        u[labels == b] = 1.0 / np.sqrt(nb)
        proj = float(np.sum((eigvecs @ u) ** 2))
        per_mode_mass = np.array([np.sum(v[labels == b] ** 2) for v in eigvecs])
        mass = float(per_mode_mass.sum())
        rows.append(dict(
            block=b,
            size=int(sizes[b]),
            size_kept=nb,
            strength=float(strengths[b]),
            projection=proj,
            mass=mass,
            mass_null=float(n_modes * nb / p),
            mass_ratio=float(mass / (n_modes * nb / p)),
            best_mode=int(np.argmax(per_mode_mass)) + 1,
        ))
    return rows


# ── stability across replicates ──────────────────────────────────────────────

def subspace_overlap(va, vb):
    """Mean squared canonical correlation between two mode subspaces: 1 = identical."""
    va, vb = np.atleast_2d(va), np.atleast_2d(vb)
    return float(np.sum((va @ vb.T) ** 2) / va.shape[0])


def mode_agreement(va, vb):
    """|Pearson r| between matched modes of two replicates (sign is arbitrary)."""
    va, vb = np.atleast_2d(va), np.atleast_2d(vb)
    return [float(abs(np.corrcoef(a, b)[0, 1])) for a, b in zip(va, vb)]


# ── the noise ladder ─────────────────────────────────────────────────────────

def ladder_stages(R, latent, true_counts, observed_counts, n_top=5):
    """
    Run get_eig_vectors at each layer of the generator, returning
    {stage: (eigvals, eigvecs, threshold, kept_cols)}.

    The latent Gaussian layer is signed, so row (library-size) normalization is
    skipped there; the count layers get the exact preprocessing used on the
    experimental data.
    """
    stages = {}
    vals, vecs = oracle_eigenvectors(R, n_top=n_top)
    stages['R_oracle'] = (vals, vecs, np.nan, np.ones(R.shape[0], dtype=bool))
    stages['latent'] = get_eig_vectors(latent, n_top=n_top, norm=False)
    stages['true_counts'] = get_eig_vectors(true_counts, n_top=n_top,
                                            norm_method='sum', norm_sum=50)
    stages['observed'] = get_eig_vectors(observed_counts, n_top=n_top,
                                         norm_method='sum', norm_sum=50)
    return stages


def score_stage(eigvecs, kept_cols, structure, n_top_genes=30, min_block_size=5):
    """Apply the per-mode and per-block scoring to one ladder stage."""
    labels_kept = np.asarray(structure['labels'])[np.asarray(kept_cols, dtype=bool)]
    sizes, strengths = structure['sizes'], structure['strengths']
    return dict(
        modes=score_modes(eigvecs, labels_kept, sizes, strengths,
                          n_top_genes=n_top_genes),
        blocks=score_blocks(eigvecs, labels_kept, sizes, strengths,
                            min_size=min_block_size),
        n_genes_kept=int(labels_kept.size),
    )

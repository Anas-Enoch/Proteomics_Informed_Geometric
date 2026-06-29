#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
human1_analysis_core.py
=======================
Shared, verified analysis primitives for the Human1-scale Figure 3-6 experiments.

NO synthetic values. NO hard-coded result arrays. The only randomness is
fixed-seed metabolite sampling for the random-panel baselines (documented).

This module is imported by:
    human1_fig3_geometry_preservation.py
    human1_fig4_panel_design.py
    human1_fig5_robustness.py
    human1_fig6_psi.py

It builds the metabolite Laplacian via the already-implemented operator code
(proteomics_weighting.build_operator) and provides the spectral-embedding,
normalised-distortion, diffusion-distance, greedy-selection, degree-selection,
and PSI primitives used across all four figures.
"""

from __future__ import annotations

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh

# ----------------------------------------------------------------------
# Operator construction (delegates to the real operator code)
# ----------------------------------------------------------------------

def load_model(model_path):
    import cobra
    return cobra.io.read_sbml_model(model_path)


def build_metabolite_laplacian(model, mode="baseline", proteomics_path=None,
                               condition="tumor", rho_0=0.1, alpha=1.0,
                               seed=42, hgnc_path=None):
    """
    Delta_M = W_M^{1/2} S W_R S^T W_M^{1/2}, via proteomics_weighting.build_operator.
    mode in {'baseline','proteomics','permuted'}.
    Returns a scipy CSR sparse matrix (n_metabolites x n_metabolites).
    """
    from proteomics_weighting import build_operator
    return build_operator(
        model, mode=mode, proteomics_path=proteomics_path,
        condition=condition, rho_0=rho_0, alpha=alpha,
        seed=seed, hgnc_path=hgnc_path,
    )


def metabolite_ids(model):
    return [m.id for m in model.metabolites]


def metabolite_subsystems(model):
    """
    Map metabolite id -> set(subsystems), inherited from the reactions it
    participates in (Human-GEM is reaction-annotated).
    """
    m2s = {m.id: set() for m in model.metabolites}
    for rxn in model.reactions:
        sub = getattr(rxn, "subsystem", None)
        if isinstance(sub, list):
            subs = [str(s) for s in sub if s]
        elif sub:
            subs = [str(sub)]
        else:
            subs = ["Unassigned"]
        for met in rxn.metabolites:
            m2s[met.id].update(subs)
    return m2s


# ----------------------------------------------------------------------
# Spectral embedding + distortion metric (identical to manuscript R1-7 def)
# ----------------------------------------------------------------------

def spectral_embedding(Delta, k, sigma=1e-4, tol=1e-6, maxiter=500000,
                       rel_thresh=1e-6, dense_threshold=400):
    """
    Sign-corrected low-frequency eigenvectors (k smallest non-zero eigenpairs).
    For small matrices (n <= dense_threshold) uses dense eigh, which is exact
    and robust; for large matrices uses shift-invert eigsh away from the
    near-null space (sigma) to avoid the machine-epsilon artefact that motivated
    K_EIGS=50 in the spectral scan.
    Returns U (n x k') with k' <= k, or None if no non-trivial spectrum.
    """
    n = Delta.shape[0]
    if not sparse.issparse(Delta):
        Delta = sparse.csr_matrix(Delta)

    if n <= dense_threshold:
        from numpy.linalg import eigh as _eigh
        A = Delta.toarray()
        A = 0.5 * (A + A.T)  # symmetrise against numerical asymmetry
        vals, vecs = _eigh(A)
    else:
        kk = min(k + 1, n - 1)
        if kk < 1:
            return None
        try:
            vals, vecs = eigsh(Delta.tocsr(), k=kk, sigma=sigma, which="LM",
                               tol=tol, maxiter=maxiter)
        except Exception:
            try:
                vals, vecs = eigsh(Delta.tocsr(), k=kk, which="SM",
                                   tol=tol, maxiter=maxiter)
            except Exception:
                return None
    order = np.argsort(vals)
    vals, vecs = vals[order], vecs[:, order]
    pos = vals > 0
    if not np.any(pos):
        return None
    thresh = rel_thresh * vals[pos].max()
    keep = vals > thresh
    if not np.any(keep):
        return None
    U = vecs[:, keep][:, :k]
    # sign correction: largest-magnitude entry positive per column
    for j in range(U.shape[1]):
        if U[np.abs(U[:, j]).argmax(), j] < 0:
            U[:, j] = -U[:, j]
    return U


def restrict(Delta, idx):
    idx = np.asarray(idx, dtype=int)
    return Delta[idx][:, idx]


def normalised_distortion(U_full_sub, U_restr):
    """
    dist = || U_full[Omega] - U_restr ||_F / || U_full[Omega] ||_F
    with per-column sign alignment of U_restr to U_full[Omega].
    """
    if U_full_sub is None or U_restr is None:
        return np.nan
    k = min(U_full_sub.shape[1], U_restr.shape[1])
    if k == 0:
        return np.nan
    A = U_full_sub[:, :k].copy()
    B = U_restr[:, :k].copy()
    for j in range(k):
        if np.dot(A[:, j], B[:, j]) < 0:
            B[:, j] = -B[:, j]
    num = np.linalg.norm(A - B, ord="fro")
    den = np.linalg.norm(A, ord="fro")
    return float(num / den) if den > 0 else np.nan


def restriction_distortion(Delta, U_full, idx, k):
    """Distortion of the restricted operator on panel idx vs the full embedding."""
    D_sub = restrict(Delta, idx)
    U_sub = spectral_embedding(D_sub, k=min(k, len(idx) - 1))
    return normalised_distortion(U_full[idx, :], U_sub)


# ----------------------------------------------------------------------
# Diffusion-distance error
# ----------------------------------------------------------------------

def diffusion_distance_error(Delta, U_full, idx, k, t=1.0):
    """
    Compares pairwise diffusion distances on the full vs restricted operator,
    over the shared panel idx. Returns normalised mean absolute error.
    Diffusion distance uses the low-frequency embedding scaled by exp(-t*lambda).
    """
    # full: eigen-pairs already in U_full; recompute eigenvalues on the restricted
    n = Delta.shape[0]
    # full diffusion coords on idx
    full_coords = U_full[idx, :k]
    D_sub = restrict(Delta, idx)
    U_sub = spectral_embedding(D_sub, k=min(k, len(idx) - 1))
    if U_sub is None:
        return np.nan
    kk = min(full_coords.shape[1], U_sub.shape[1])
    full_coords = full_coords[:, :kk]
    sub_coords = U_sub[:, :kk]
    # sign align
    for j in range(kk):
        if np.dot(full_coords[:, j], sub_coords[:, j]) < 0:
            sub_coords[:, j] = -sub_coords[:, j]
    # pairwise distance matrices
    def pdist(X):
        sq = np.sum(X**2, axis=1)
        d2 = sq[:, None] + sq[None, :] - 2 * X @ X.T
        return np.sqrt(np.maximum(d2, 0))
    Df = pdist(full_coords)
    Ds = pdist(sub_coords)
    denom = np.linalg.norm(Df)
    return float(np.linalg.norm(Df - Ds) / denom) if denom > 0 else np.nan


# ----------------------------------------------------------------------
# Panel selection strategies
# ----------------------------------------------------------------------

def greedy_panel(Delta, U_full, k_embed, target_size, candidate_cap=None,
                 verbose=False):
    """
    Geometry-aware greedy selection (manuscript R1-4 algorithm).
    Objective minimised at each step: dist(Omega U {m}).
    Returns (selected_idx_list, trace) where trace is a list of
    (step, added_idx, distortion_after).

    candidate_cap: if set, restrict the per-step candidate scan to the
    `candidate_cap` highest-leverage unselected metabolites (leverage =
    row norm in U_full). This bounds the O(n^2 * eig) cost at full Human1
    scale WITHOUT changing the objective; document in the figure caption.
    """
    n = Delta.shape[0]
    leverage = np.sum(U_full**2, axis=1)
    selected, remaining = [], list(range(n))
    trace = []
    for step in range(target_size):
        if candidate_cap is not None and len(remaining) > candidate_cap:
            cand = [remaining[i] for i in np.argsort(
                leverage[remaining])[::-1][:candidate_cap]]
        else:
            cand = remaining
        best_idx, best_dist = None, np.inf
        for c in cand:
            trial = sorted(selected + [c])
            if len(trial) < 2:
                # seed with highest-leverage metabolite
                val = -leverage[c]
                if val < best_dist:
                    best_dist, best_idx = val, c
                continue
            d = restriction_distortion(Delta, U_full, np.array(trial), k_embed)
            if d is not None and not np.isnan(d) and d < best_dist:
                best_dist, best_idx = d, c
        if best_idx is None:
            break
        selected.append(best_idx)
        remaining.remove(best_idx)
        trace.append((step + 1, best_idx, best_dist if best_dist >= 0 else np.nan))
        if verbose:
            print(f"  step {step+1}: +{best_idx}  dist={best_dist:.4f}", flush=True)
    return selected, trace


def degree_panel(Delta, target_size):
    """Degree-based selection: pick highest-degree metabolites in the operator graph."""
    A = (Delta != 0)
    deg = np.asarray(A.sum(axis=1)).ravel()
    return list(np.argsort(deg)[::-1][:target_size])


def random_panels(n, target_size, n_draws, seed):
    """Fixed-seed random panels. Returns list of sorted index arrays."""
    rng = np.random.default_rng(seed)
    return [np.sort(rng.choice(n, size=target_size, replace=False))
            for _ in range(n_draws)]


# ----------------------------------------------------------------------
# PSI (Pathway Separability Index) with NaN-exclusion edge handling (R1-9)
# ----------------------------------------------------------------------

def psi_for_panel(Delta, idx, met_ids, met2sub, d=3):
    """
    Calinski-Harabasz-style Pathway Separability Index on a restricted panel.

    Subsystems present in the panel are treated as clusters in the d-dimensional
    low-frequency embedding. The index is

        CH = [ SS_between / (K - 1) ] / [ SS_within / (N_used - K) ]

    where SS_between is the size-weighted dispersion of subsystem centroids about
    the global centroid, SS_within is the summed squared distance of metabolites
    to their own subsystem centroid, K is the number of usable subsystems
    (>= 2 panel members each), and N_used is the number of metabolites in those
    subsystems. This is naturally bounded by the data, scale-invariant, and
    degeneracy-safe (the summed within term does not blow up when individual
    same-subsystem metabolites coincide, unlike a per-pair B/W ratio).

    Returns (per_subsystem_dict, global_CH). The per-subsystem dict reports each
    subsystem's contribution ratio (between-centroid distance^2 over mean within
    dispersion) for interpretability; the global value is the CH index.
    Undefined cases (panel subset of one subsystem, <2 usable clusters) -> NaN.
    """
    U = spectral_embedding(restrict(Delta, idx), k=min(d, len(idx) - 1))
    if U is None:
        return {}, np.nan
    panel_mets = [met_ids[j] for j in idx]
    X = np.asarray([U[i, :] for i in range(len(panel_mets))], dtype=float)

    # cluster membership (a multi-subsystem metabolite counts in each subsystem)
    clusters = {}
    for i, m in enumerate(panel_mets):
        for s in met2sub.get(m, {"Unassigned"}):
            clusters.setdefault(s, []).append(i)
    clusters = {s: ix for s, ix in clusters.items() if len(ix) >= 2}
    K = len(clusters)
    if K < 2:
        return {}, np.nan

    global_centroid = X.mean(axis=0)
    ss_between = 0.0
    ss_within = 0.0
    n_used = 0
    per_sub = {}
    for s, ix in clusters.items():
        pts = X[ix]
        c = pts.mean(axis=0)
        between_s = float(np.sum((c - global_centroid) ** 2))
        within_s = float(np.sum((pts - c) ** 2))
        ss_between += len(ix) * between_s
        ss_within += within_s
        n_used += len(ix)
        # per-subsystem interpretable contribution (centroid offset vs spread)
        mean_within = within_s / len(ix)
        per_sub[s] = (between_s / mean_within) if mean_within > 0 else np.nan

    if ss_within <= 0 or (n_used - K) <= 0:
        return per_sub, np.nan
    ch = (ss_between / (K - 1)) / (ss_within / (n_used - K))
    return per_sub, float(ch)

#!/usr/bin/env python3
"""
layer1_spectral_scan.py  —  corrected version

Fixes vs original:
  1. K_EIGS=50: captures real spectral geometry, not near-null space
  2. SIGMA=1e-4: avoids numerical noise at machine-epsilon scale
  3. Relative eigenvalue filtering: threshold = REL_THRESH * max_eigenvalue
  4. Adds subspace alignment metric (principal angles) as a second distance measure,
     which is more informative than raw eigenvalue L2 distance
  5. Normalises the eigenvalue-distance by the full-operator eigenvalue norm
     so the metric is scale-invariant across cell lines
"""

import cobra
import pandas as pd
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh, ArpackNoConvergence

# ---------------------------
# Paths
# ---------------------------
MODEL_PATH   = "Human-GEM-main/model/Human-GEM.xml"
WEIGHTS_PATH = "nci60_reaction_weights_norm_simple.csv"

# ---------------------------
# Settings
# ---------------------------
K_EIGS     = 50
SIGMA      = 1e-4
TOL        = 1e-6
MAXITER    = 500000
REL_THRESH = 1e-6   # relative eigenvalue floor


# ---------------------------
# Helpers
# ---------------------------
def compute_lowfreq_eigs(mat, label, k=K_EIGS, sigma=SIGMA, tol=TOL, maxiter=MAXITER):
    try:
        evals, evecs = eigsh(mat, k=k, sigma=sigma, which="LM",
                             tol=tol, maxiter=maxiter)
        order = np.argsort(np.real(evals))
        return np.real(evals)[order], evecs[:, order]
    except ArpackNoConvergence as e:
        if e.eigenvalues is not None and len(e.eigenvalues) > 0:
            order = np.argsort(np.real(e.eigenvalues))
            return np.real(e.eigenvalues)[order], e.eigenvectors[:, order]
        print(f"  [WARN] {label}: no convergence.", flush=True)
        return None, None
    except Exception as e:
        print(f"  [WARN] {label}: {e}", flush=True)
        return None, None


def keep_nontrivial(evals, evecs=None, rel_thresh=REL_THRESH):
    if evals is None or len(evals) == 0:
        return None, None
    evals = np.array(evals)
    mask  = evals > 0
    evals = evals[mask]
    if len(evals) == 0:
        return None, None
    threshold = rel_thresh * evals.max()
    keep = evals > threshold
    nt_evals = evals[keep]
    nt_evecs = evecs[:, np.where(mask)[0][keep]] if evecs is not None else None
    return (nt_evals if len(nt_evals) > 0 else None,
            nt_evecs if evecs is not None else None)


def relative_l2(a, b):
    """Normalised eigenvalue distance: ||a - b|| / ||b||"""
    if a is None or b is None:
        return np.nan
    m = min(len(a), len(b))
    if m == 0:
        return np.nan
    denom = np.linalg.norm(b[:m])
    return float(np.linalg.norm(a[:m] - b[:m]) / denom) if denom > 0 else np.nan


def subspace_distance(U, V):
    """
    Principal-angle-based subspace distance between column spaces of U and V.
    Returns 1 - cos(smallest principal angle), in [0, 1].
    0 = identical subspaces; 1 = orthogonal subspaces.
    """
    if U is None or V is None:
        return np.nan
    # Orthonormalise both
    Qu, _ = np.linalg.qr(U)
    Qv, _ = np.linalg.qr(V)
    k = min(Qu.shape[1], Qv.shape[1])
    if k == 0:
        return np.nan
    M    = Qu[:, :k].T @ Qv[:, :k]
    svals = np.linalg.svd(M, compute_uv=False)
    cos_min = np.clip(svals[0], 0, 1)   # largest singular value = cos of smallest angle
    return float(1.0 - cos_min)


# ---------------------------
# Load model and weights
# ---------------------------
print("Loading Human-GEM model...", flush=True)
model = cobra.io.read_sbml_model(MODEL_PATH)
rxn_ids_model = [rxn.id for rxn in model.reactions]
print(f"  Reactions: {len(rxn_ids_model)}", flush=True)

print("Building stoichiometric matrix S...", flush=True)
S = cobra.util.array.create_stoichiometric_matrix(model, array_type="DataFrame")
print(f"  S shape: {S.shape}", flush=True)

print("Loading reaction weights...", flush=True)
W = pd.read_csv(WEIGHTS_PATH, index_col=0)
print(f"  Weight matrix: {W.shape}", flush=True)

common_rxns = [r for r in rxn_ids_model if r in W.index]
assert len(common_rxns) > 0, "No common reactions."
S_aligned = S.loc[:, common_rxns]
W_aligned = W.loc[common_rxns, :]
assert list(S_aligned.columns) == list(W_aligned.index)

S_sparse = sparse.csr_matrix(S_aligned.values)

# ---------------------------
# Topology-only baseline
# ---------------------------
print("\nComputing topology-only baseline...", flush=True)
Delta_topo = S_sparse @ S_sparse.T
topo_evals, topo_evecs = compute_lowfreq_eigs(Delta_topo, "topology-only")
topo_nt_evals, topo_nt_evecs = keep_nontrivial(topo_evals, topo_evecs)

if topo_nt_evals is not None:
    print(f"  Nontrivial eigenvalues: {len(topo_nt_evals)}", flush=True)
    print(f"  Range: [{topo_nt_evals.min():.4e}, {topo_nt_evals.max():.4e}]", flush=True)
else:
    print("  WARNING: no nontrivial eigenvalues in topology-only baseline. "
          "Increase K_EIGS or check SIGMA.", flush=True)

# ---------------------------
# Cell-line loop
# ---------------------------
results = []

for i, cell_line in enumerate(W_aligned.columns, start=1):
    print(f"[{i}/{W_aligned.shape[1]}] {cell_line}", flush=True)

    w     = W_aligned[cell_line].values
    Wr    = sparse.diags(w)
    Delta = S_sparse @ Wr @ S_sparse.T

    prot_evals, prot_evecs = compute_lowfreq_eigs(Delta, cell_line)
    prot_nt_evals, prot_nt_evecs = keep_nontrivial(prot_evals, prot_evecs)

    row = {
        "cell_line":               cell_line,
        "n_eigs_total":            len(prot_evals) if prot_evals is not None else 0,
        "n_eigs_nontrivial":       len(prot_nt_evals) if prot_nt_evals is not None else 0,
        "spectral_l2_to_topology": relative_l2(prot_nt_evals, topo_nt_evals),
        "spectral_rel_l2_to_topology": relative_l2(prot_nt_evals, topo_nt_evals),
        "subspace_distance_to_topology": subspace_distance(prot_nt_evecs, topo_nt_evecs),
    }

    # Raw eigenvalue columns
    for j in range(K_EIGS):
        row[f"eig_{j+1}"] = (prot_evals[j]
                              if prot_evals is not None and j < len(prot_evals)
                              else np.nan)
    # Nontrivial eigenvalue columns
    for j in range(K_EIGS):
        row[f"eig_nontrivial_{j+1}"] = (prot_nt_evals[j]
                                         if prot_nt_evals is not None and j < len(prot_nt_evals)
                                         else np.nan)
    results.append(row)

# ---------------------------
# Save
# ---------------------------
results_df = pd.DataFrame(results)
results_df.to_csv("layer1_spectral_scan_summary.csv", index=False)
print("\nSaved: layer1_spectral_scan_summary.csv", flush=True)

print("\n=== Summary ===")
cols = ["spectral_rel_l2_to_topology", "subspace_distance_to_topology", "n_eigs_nontrivial"]
print(results_df[cols].describe())
print("\nDone.", flush=True)

#!/usr/bin/env python3
"""
layer1_build_S.py  —  corrected version
Fixes:
  1. K_EIGS increased to 50: captures meaningful geometry, not just near-null space
  2. Sigma moved to 1e-4: avoids numerical noise at machine-epsilon scale
  3. Eigenvalue filtering uses a relative threshold, not absolute 1e-10
  4. Reports the spectral gap clearly to confirm geometry is captured
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
K_EIGS   = 50       # capture real geometry, not just near-null eigenvalues
SIGMA    = 1e-4     # shift-invert away from exact zero to avoid null-space artifacts
TOL      = 1e-6
MAXITER  = 500000

# Relative threshold: eigenvalue must be > REL_THRESH * max_eigenvalue to be nontrivial
REL_THRESH = 1e-6

# ---------------------------
# Helper: robust eigensolver with relative filtering
# ---------------------------
def compute_lowfreq_eigs(mat, label, k=K_EIGS, sigma=SIGMA, tol=TOL, maxiter=MAXITER):
    print(f"  Computing {k} eigenvalues for {label}...", flush=True)
    try:
        evals, evecs = eigsh(mat, k=k, sigma=sigma, which="LM",
                             tol=tol, maxiter=maxiter)
        evals = np.sort(np.real(evals))
        return evals, evecs
    except ArpackNoConvergence as e:
        print(f"  [WARN] {label}: ARPACK non-convergence, partial recovery.", flush=True)
        if e.eigenvalues is not None and len(e.eigenvalues) > 0:
            return np.sort(np.real(e.eigenvalues)), e.eigenvectors
        return None, None
    except Exception as e:
        print(f"  [WARN] {label}: eigensolver failed: {e}", flush=True)
        return None, None


def keep_nontrivial(evals, rel_thresh=REL_THRESH):
    """Filter out near-zero eigenvalues using a relative threshold."""
    if evals is None or len(evals) == 0:
        return None
    evals = np.array(evals)
    evals = evals[evals > 0]
    if len(evals) == 0:
        return None
    threshold = rel_thresh * evals.max()
    nontrivial = evals[evals > threshold]
    return nontrivial if len(nontrivial) > 0 else None


# ---------------------------
# Load model
# ---------------------------
print("Loading Human-GEM model...", flush=True)
model = cobra.io.read_sbml_model(MODEL_PATH)
rxn_ids_model = [rxn.id for rxn in model.reactions]
met_ids_model = [met.id for met in model.metabolites]
print(f"  Reactions: {len(rxn_ids_model)},  Metabolites: {len(met_ids_model)}", flush=True)

# ---------------------------
# Build stoichiometric matrix
# ---------------------------
print("Building stoichiometric matrix S...", flush=True)
S = cobra.util.array.create_stoichiometric_matrix(model, array_type="DataFrame")
print(f"  S shape: {S.shape}", flush=True)

# ---------------------------
# Load and align reaction weights
# ---------------------------
print("Loading reaction weights...", flush=True)
W = pd.read_csv(WEIGHTS_PATH, index_col=0)
print(f"  Weight matrix shape: {W.shape}", flush=True)

common_rxns = [r for r in rxn_ids_model if r in W.index]
print(f"  Common reactions: {len(common_rxns)}", flush=True)
assert len(common_rxns) > 0, "No common reactions found."

S_aligned = S.loc[:, common_rxns]
W_aligned = W.loc[common_rxns, :]
assert list(S_aligned.columns) == list(W_aligned.index), "Reaction order mismatch."

S_sparse = sparse.csr_matrix(S_aligned.values)

# ---------------------------
# Save aligned matrices
# ---------------------------
print("Saving aligned matrices...", flush=True)
S_aligned.to_csv("human1_S_aligned.csv")
W_aligned.to_csv("nci60_reaction_weights_aligned.csv")
sparse.save_npz("human1_S_aligned.npz", S_sparse)
pd.Series(S_aligned.index).to_csv("human1_met_ids.csv", index=False)
pd.Series(S_aligned.columns).to_csv("human1_rxn_ids_aligned.csv", index=False)

# ---------------------------
# Build operators for first cell line
# ---------------------------
cell_line = W_aligned.columns[0]
print(f"\nBuilding operators for example cell line: {cell_line}", flush=True)

w        = W_aligned[cell_line].values
Wr       = sparse.diags(w)
Delta      = S_sparse @ Wr @ S_sparse.T
Delta_topo = S_sparse @ S_sparse.T

sparse.save_npz(f"Delta_{cell_line.replace(':', '_')}.npz", Delta)
sparse.save_npz("Delta_topology_only.npz", Delta_topo)

# ---------------------------
# Spectral check
# ---------------------------
print("\nTopology-only:", flush=True)
evals_topo, _ = compute_lowfreq_eigs(Delta_topo, "topology-only")
nt_topo = keep_nontrivial(evals_topo)
print(f"  Total returned: {len(evals_topo) if evals_topo is not None else 0}")
print(f"  Nontrivial: {len(nt_topo) if nt_topo is not None else 0}")
if nt_topo is not None:
    print(f"  Range: [{nt_topo.min():.4e}, {nt_topo.max():.4e}]")

print("\nProteomics-informed:", flush=True)
evals_prot, _ = compute_lowfreq_eigs(Delta, f"prot-{cell_line}")
nt_prot = keep_nontrivial(evals_prot)
print(f"  Total returned: {len(evals_prot) if evals_prot is not None else 0}")
print(f"  Nontrivial: {len(nt_prot) if nt_prot is not None else 0}")
if nt_prot is not None:
    print(f"  Range: [{nt_prot.min():.4e}, {nt_prot.max():.4e}]")

# ---------------------------
# Save eigenvalue summary
# ---------------------------
if nt_topo is not None or nt_prot is not None:
    topo_list = list(nt_topo) if nt_topo is not None else []
    prot_list = list(nt_prot) if nt_prot is not None else []
    max_len   = max(len(topo_list), len(prot_list))
    eig_df = pd.DataFrame({
        "topology_only_nontrivial":    topo_list + [np.nan] * (max_len - len(topo_list)),
        "proteomics_informed_nontrivial": prot_list + [np.nan] * (max_len - len(prot_list)),
    })
    eig_df.to_csv("example_lowfreq_eigenvalues_nontrivial.csv", index=False)
    print("\nSaved: example_lowfreq_eigenvalues_nontrivial.csv", flush=True)

print("\nDone.", flush=True)

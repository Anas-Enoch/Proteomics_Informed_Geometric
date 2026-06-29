#!/usr/bin/env python3

import cobra
import pandas as pd
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh, ArpackNoConvergence

# ---------------------------
# Paths
# ---------------------------
MODEL_PATH = "Human-GEM-main/model/Human-GEM.xml"
WEIGHTS_PATH = "nci60_reaction_weights_norm_simple.csv"

# ---------------------------
# Settings
# ---------------------------
K_EIGS = 50
SIGMA = 1e-4
TOL = 1e-5
MAXITER = 200000
REL_THRESH = 1e-6       # keep eigenvalues > REL_THRESH * max(eigs)
TARGET_RANK = 10        # rank used for panel/subspace comparisons

PANEL_SIZES = [50, 100, 200, 500]
N_RANDOM_TRIALS = 10
RANDOM_SEED = 42

# ---------------------------
# Helpers
# ---------------------------
def compute_eigs_and_vecs(mat, k=K_EIGS, sigma=SIGMA, tol=TOL, maxiter=MAXITER):
    try:
        evals, evecs = eigsh(
            mat,
            k=min(k, mat.shape[0] - 1),
            sigma=sigma,
            which="LM",
            tol=tol,
            maxiter=maxiter,
        )
        order = np.argsort(np.real(evals))
        evals = np.real(evals[order])
        evecs = np.real(evecs[:, order])
        return evals, evecs
    except ArpackNoConvergence as e:
        if e.eigenvalues is not None and len(e.eigenvalues) > 0:
            order = np.argsort(np.real(e.eigenvalues))
            evals = np.real(e.eigenvalues[order])
            evecs = np.real(e.eigenvectors[:, order]) if e.eigenvectors is not None else None
            return evals, evecs
        return None, None
    except Exception:
        return None, None


def keep_nontrivial(evals, evecs, rel_thresh=REL_THRESH):
    if evals is None or evecs is None or len(evals) == 0:
        return None, None
    max_eval = np.max(np.abs(evals))
    if max_eval <= 0:
        return None, None
    mask = evals > (rel_thresh * max_eval)
    if np.sum(mask) == 0:
        return None, None
    return evals[mask], evecs[:, mask]


def restrict_operator(Delta, idx):
    return Delta[idx, :][:, idx]


def subspace_distortion(U_full_rows, U_restr, target_rank=TARGET_RANK):
    """
    Compare the subspace spanned by:
      - restricted rows of the full operator basis
      - basis of the restricted operator

    Uses projection-matrix Frobenius distance, which is sign- and rotation-invariant.
    """
    if U_full_rows is None or U_restr is None:
        return np.nan

    r = min(target_rank, U_full_rows.shape[1], U_restr.shape[1], U_full_rows.shape[0], U_restr.shape[0])
    if r < 1:
        return np.nan

    Uf = U_full_rows[:, :r]
    Ur = U_restr[:, :r]

    # Orthonormalize both on the restricted panel
    try:
        Qf, _ = np.linalg.qr(Uf)
        Qr, _ = np.linalg.qr(Ur)
    except Exception:
        return np.nan

    # Projection matrices
    Pf = Qf @ Qf.T
    Pr = Qr @ Qr.T

    denom = np.linalg.norm(Pf, ord="fro")
    if denom == 0:
        return np.nan

    return float(np.linalg.norm(Pf - Pr, ord="fro") / denom)


def leverage_scores(U_full, target_rank=TARGET_RANK):
    if U_full is None:
        return None
    r = min(target_rank, U_full.shape[1])
    if r < 1:
        return None
    return np.sum(U_full[:, :r] ** 2, axis=1)


def geometry_panel_from_leverage(U_full, panel_size):
    """
    Deterministic geometry-aware panel:
    select metabolites with largest low-frequency leverage scores.
    """
    lev = leverage_scores(U_full, TARGET_RANK)
    if lev is None:
        return np.array([], dtype=int)
    idx = np.argsort(-lev)[:panel_size]
    return np.sort(idx)


def random_panel(n_mets, panel_size, rng):
    return np.sort(rng.choice(n_mets, size=panel_size, replace=False))


def summarize_and_save(df, group_col, metric_cols, out_csv, label):
    if df.empty:
        print(f"{label}: no valid rows to summarize.", flush=True)
        return None
    summary = (
        df.groupby(group_col)[metric_cols]
        .agg(["mean", "median", "std", "count"])
    )
    summary.to_csv(out_csv)
    print(f"Saved: {out_csv}", flush=True)
    print(f"\n=== {label} summary ===", flush=True)
    print(summary, flush=True)
    return summary


# ---------------------------
# Load model and weights
# ---------------------------
print("Loading Human-GEM model...", flush=True)
model = cobra.io.read_sbml_model(MODEL_PATH)

rxn_ids_model = [rxn.id for rxn in model.reactions]
met_ids_model = [met.id for met in model.metabolites]

print(f"Model reactions: {len(rxn_ids_model)}", flush=True)
print(f"Model metabolites: {len(met_ids_model)}", flush=True)

print("Building stoichiometric matrix S...", flush=True)
S = cobra.util.array.create_stoichiometric_matrix(model, array_type="DataFrame")

print("Loading reaction weights...", flush=True)
W = pd.read_csv(WEIGHTS_PATH, index_col=0)

common_rxns = [r for r in rxn_ids_model if r in W.index]
if len(common_rxns) == 0:
    raise ValueError("No common reactions between Human-GEM and weight matrix.")

S_aligned = S.loc[:, common_rxns]
W_aligned = W.loc[common_rxns, :]

assert list(S_aligned.columns) == list(W_aligned.index), "Reaction order mismatch."

S_sparse = sparse.csr_matrix(S_aligned.values)
n_mets = S_sparse.shape[0]

# ---------------------------
# Topology-only baseline
# ---------------------------
print("Building topology-only operator...", flush=True)
Delta_topo = S_sparse @ S_sparse.T

topo_evals, topo_evecs = compute_eigs_and_vecs(Delta_topo)
topo_evals, topo_evecs = keep_nontrivial(topo_evals, topo_evecs)

if topo_evecs is None:
    raise RuntimeError("Topology-only operator did not yield usable nontrivial spectrum.")

print(f"Topology-only nontrivial eigenvalue count: {len(topo_evals)}", flush=True)

# ---------------------------
# Experiment A
# Correct claim test:
# geometry-aware panel vs random panels on topology-only operator
# ---------------------------
print("\n=== Experiment A: geometry-aware vs random on topology-only ===", flush=True)
rng = np.random.default_rng(RANDOM_SEED)

results_A = []

for panel_size in PANEL_SIZES:
    print(f"Panel size {panel_size}...", flush=True)

    # deterministic geometry-aware panel from full topology basis
    idx_geom = geometry_panel_from_leverage(topo_evecs, panel_size)

    if len(idx_geom) != panel_size:
        print(f"[WARN] Geometry-aware panel failed at size {panel_size}.", flush=True)
        continue

    Delta_geom = restrict_operator(Delta_topo, idx_geom)
    geom_evals, geom_evecs = compute_eigs_and_vecs(Delta_geom)
    geom_evals, geom_evecs = keep_nontrivial(geom_evals, geom_evecs)

    geom_dist = subspace_distortion(topo_evecs[idx_geom, :], geom_evecs, TARGET_RANK)

    rand_dists = []
    for trial in range(1, N_RANDOM_TRIALS + 1):
        idx_rand = random_panel(n_mets, panel_size, rng)
        Delta_rand = restrict_operator(Delta_topo, idx_rand)
        rand_evals, rand_evecs = compute_eigs_and_vecs(Delta_rand)
        rand_evals, rand_evecs = keep_nontrivial(rand_evals, rand_evecs)

        rand_dist = subspace_distortion(topo_evecs[idx_rand, :], rand_evecs, TARGET_RANK)
        if np.isfinite(rand_dist):
            rand_dists.append(rand_dist)

        results_A.append({
            "panel_size": panel_size,
            "trial": trial,
            "panel_type": "random",
            "distortion": rand_dist,
        })

    rand_mean = float(np.mean(rand_dists)) if len(rand_dists) > 0 else np.nan
    advantage = rand_mean - geom_dist if np.isfinite(rand_mean) and np.isfinite(geom_dist) else np.nan

    results_A.append({
        "panel_size": panel_size,
        "trial": 0,
        "panel_type": "geometry_aware",
        "distortion": geom_dist,
    })

    results_A.append({
        "panel_size": panel_size,
        "trial": -1,
        "panel_type": "summary",
        "distortion": np.nan,
        "geometry_aware_distortion": geom_dist,
        "random_mean_distortion": rand_mean,
        "advantage": advantage,
        "n_valid_random_trials": len(rand_dists),
    })

df_A = pd.DataFrame(results_A)
df_A.to_csv("layer1_experiment_A_raw.csv", index=False)
print("Saved: layer1_experiment_A_raw.csv", flush=True)

if df_A.empty:
    print("Experiment A produced no rows.", flush=True)
else:
    df_A_summary = df_A[df_A["panel_type"] == "summary"].copy()
    df_A_summary.to_csv("layer1_experiment_A_summary_rows.csv", index=False)
    print("Saved: layer1_experiment_A_summary_rows.csv", flush=True)

    if not df_A_summary.empty:
        summarize_and_save(
            df_A_summary,
            "panel_size",
            ["geometry_aware_distortion", "random_mean_distortion", "advantage"],
            "layer1_experiment_A_summary.csv",
            "Experiment A"
        )
    else:
        print("Experiment A: no valid summary rows.", flush=True)

# ---------------------------
# Experiment B
# Secondary operator-type comparison under random masking
# ---------------------------
print("\n=== Experiment B: proteomics-informed vs topology-only under random masking ===", flush=True)

results_B = []

for i, cell_line in enumerate(W_aligned.columns, start=1):
    print(f"[{i}/{W_aligned.shape[1]}] {cell_line}", flush=True)

    w = W_aligned[cell_line].values
    Wr = sparse.diags(w)
    Delta_prot = S_sparse @ Wr @ S_sparse.T

    prot_evals, prot_evecs = compute_eigs_and_vecs(Delta_prot)
    prot_evals, prot_evecs = keep_nontrivial(prot_evals, prot_evecs)

    if prot_evecs is None:
        print(f"[WARN] {cell_line}: no usable proteomics-informed spectrum; skipping.", flush=True)
        continue

    for panel_size in PANEL_SIZES:
        for trial in range(1, N_RANDOM_TRIALS + 1):
            idx = random_panel(n_mets, panel_size, rng)

            Delta_topo_mask = restrict_operator(Delta_topo, idx)
            Delta_prot_mask = restrict_operator(Delta_prot, idx)

            topo_mask_evals, topo_mask_evecs = compute_eigs_and_vecs(Delta_topo_mask)
            topo_mask_evals, topo_mask_evecs = keep_nontrivial(topo_mask_evals, topo_mask_evecs)

            prot_mask_evals, prot_mask_evecs = compute_eigs_and_vecs(Delta_prot_mask)
            prot_mask_evals, prot_mask_evecs = keep_nontrivial(prot_mask_evals, prot_mask_evecs)

            topo_dist = subspace_distortion(topo_evecs[idx, :], topo_mask_evecs, TARGET_RANK)
            prot_dist = subspace_distortion(prot_evecs[idx, :], prot_mask_evecs, TARGET_RANK)

            results_B.append({
                "cell_line": cell_line,
                "panel_size": panel_size,
                "trial": trial,
                "topology_distortion": topo_dist,
                "proteomics_distortion": prot_dist,
                "delta_distortion": prot_dist - topo_dist if np.isfinite(prot_dist) and np.isfinite(topo_dist) else np.nan,
            })

df_B = pd.DataFrame(results_B)
df_B.to_csv("layer1_experiment_B_raw.csv", index=False)
print("Saved: layer1_experiment_B_raw.csv", flush=True)

if df_B.empty:
    print("Experiment B produced no rows.", flush=True)
else:
    summarize_and_save(
        df_B,
        "panel_size",
        ["topology_distortion", "proteomics_distortion", "delta_distortion"],
        "layer1_experiment_B_summary.csv",
        "Experiment B"
    )

print("\nDone.", flush=True)

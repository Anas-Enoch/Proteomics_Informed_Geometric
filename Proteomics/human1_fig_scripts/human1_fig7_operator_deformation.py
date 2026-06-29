#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
human1_fig7_operator_deformation.py
===================================
Real Figure 7: quantify how CPTAC proteomics-derived reaction weights deform the
metabolite operator relative to topology, at full Human1 scale.

NO synthetic values, NO illustrative weights, NO hard-coded arrays. The reaction
weights W_R^(c) are computed by the SAME functions used for Figures 3-6
(proteomics_weighting._compute_reaction_weights etc.); this script only EXPOSES
them (build_operator returns Delta only, so we call the underlying helpers to
also obtain W_R and S) and then computes deformation diagnostics.

Operator:  Delta_M = W_M^{1/2} S W_R^(c) S^T W_M^{1/2}   (W_M = I here, as in the
cohort/NCI-60 work), compared against the topology baseline Delta_topo = S S^T.

Outputs (all under --outdir):
  human1_fig7_reaction_weights.csv   (every reaction's W_R)
  human1_fig7_spectral_shift.csv     (topology vs proteomics eigenvalues)
  human1_fig7_subspace_rotation.csv  (principal angles, Grassmann, projection dist)
  human1_fig7_leverage_changes.csv   (per-metabolite Δ leverage + subsystem)
  human1_fig7_summary.csv            (coverage + W_R stats)

Usage:
  python human1_fig7_operator_deformation.py \
    --model Human-GEM-main/model/Human-GEM.xml \
    --proteomics data/proteomics/cptac_breast_tumor_only.tsv \
    --hgnc data/hgnc_complete_set.txt --outdir results --k 50
"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.linalg import eigsh

import human1_analysis_core as core

# Import the EXISTING weight pipeline (no reimplementation)
import proteomics_weighting as pw


def compute_WR_and_S(model, proteomics_path, condition, hgnc_path,
                     rho_0=0.1, alpha=1.0):
    """
    Reproduce build_operator's proteomics branch using ONLY the existing helper
    functions, returning (W_R, S, gene_expr, coverage_dict). This is the accessor
    that build_operator lacks (it returns Delta only).
    """
    S = pw._get_stoichiometric_matrix(model)
    gene_expr = pw._load_and_normalize_proteomics(proteomics_path, condition)
    gene_expr = pw._remap_symbols_to_ensembl(model, gene_expr, hgnc_path=hgnc_path)
    pw._report_gene_coverage(model, gene_expr)
    W_R = pw._compute_reaction_weights(model, gene_expr, rho_0, alpha)

    # coverage numbers
    all_model_genes = {g.id for rxn in model.reactions for g in rxn.genes}
    matched = all_model_genes & set(gene_expr.keys())
    observed = set(gene_expr.keys())
    n_rxn_covered = sum(
        1 for rxn in model.reactions
        if rxn.genes and (len({g.id for g in rxn.genes} & observed) > 0)
    )
    coverage = dict(
        n_proteomics_genes=len(gene_expr),
        n_model_gpr_genes=len(all_model_genes),
        n_genes_matched=len(matched),
        pct_genes_matched=100.0 * len(matched) / max(len(all_model_genes), 1),
        n_reactions=len(model.reactions),
        n_reactions_covered=n_rxn_covered,
        pct_reactions_covered=100.0 * n_rxn_covered / len(model.reactions),
    )
    return W_R, S, gene_expr, coverage


def low_freq_eigs(Delta, k):
    """k smallest non-trivial eigenvalues + vectors (dense-safe via core helper
    for the vectors; eigenvalues via the same shift-invert path)."""
    U = core.spectral_embedding(Delta, k=k)
    # eigenvalues aligned to U: Rayleigh quotients u^T Delta u
    if U is None:
        return None, None
    D = Delta.tocsr() if sparse.issparse(Delta) else sparse.csr_matrix(Delta)
    vals = np.array([float(U[:, j].T @ (D @ U[:, j])) for j in range(U.shape[1])])
    order = np.argsort(vals)
    return vals[order], U[:, order]


def principal_angles(U1, U2):
    """Principal angles between two orthonormal-column subspaces via SVD of U1^T U2."""
    # orthonormalise columns (QR) to be safe
    Q1, _ = np.linalg.qr(U1)
    Q2, _ = np.linalg.qr(U2)
    M = Q1.T @ Q2
    s = np.linalg.svd(M, compute_uv=False)
    s = np.clip(s, -1.0, 1.0)
    return np.arccos(s)  # radians, ascending cos -> descending angle


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True)
    ap.add_argument("--proteomics", required=True)
    ap.add_argument("--hgnc", default=None)
    ap.add_argument("--condition", default="tumor")
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--k", type=int, default=50, help="low-frequency subspace dim")
    ap.add_argument("--top", type=int, default=20)
    args = ap.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    print("Loading model...", flush=True)
    model = core.load_model(args.model)
    met_ids = core.metabolite_ids(model)
    met2sub = core.metabolite_subsystems(model)

    print("Computing W_R^(c) from CPTAC via existing pipeline...", flush=True)
    W_R, S, gene_expr, coverage = compute_WR_and_S(
        model, args.proteomics, args.condition, args.hgnc)

    # ---- VALIDATION PRINTS ------------------------------------------------
    print("\n=== COVERAGE ===")
    print(f"  proteomics genes loaded : {coverage['n_proteomics_genes']}")
    print(f"  model GPR genes         : {coverage['n_model_gpr_genes']}")
    print(f"  genes matched           : {coverage['n_genes_matched']} "
          f"({coverage['pct_genes_matched']:.1f}%)")
    print(f"  reactions               : {coverage['n_reactions']}")
    print(f"  reactions covered       : {coverage['n_reactions_covered']} "
          f"({coverage['pct_reactions_covered']:.1f}%)")
    print("\n=== W_R STATISTICS ===")
    print(f"  min={W_R.min():.4f}  max={W_R.max():.4f}  "
          f"mean={W_R.mean():.4f}  std={W_R.std():.4f}")

    # ---- 1. reaction-weight distribution ----------------------------------
    rxn_rows = []
    for i, rxn in enumerate(model.reactions):
        rxn_rows.append(dict(
            reaction_id=rxn.id,
            subsystem=str(getattr(rxn, "subsystem", "") or "Unassigned"),
            W_R=float(W_R[i]),
            n_genes=len(rxn.genes),
        ))
    pd.DataFrame(rxn_rows).to_csv(outdir / "human1_fig7_reaction_weights.csv", index=False)

    # ---- build operators --------------------------------------------------
    print("\nBuilding operators (topology vs proteomics)...", flush=True)
    Delta_topo = (S @ S.T).tocsr()
    Delta_prot = (S @ sparse.diags(W_R, format="csr") @ S.T).tocsr()

    print(f"Computing {args.k} low-frequency eigenpairs for each...", flush=True)
    eig_topo, U_topo = low_freq_eigs(Delta_topo, args.k)
    eig_prot, U_prot = low_freq_eigs(Delta_prot, args.k)
    if eig_topo is None or eig_prot is None:
        raise RuntimeError("Eigendecomposition failed; cannot produce Figure 7.")

    # ---- 2. spectral deformation ------------------------------------------
    kk = min(len(eig_topo), len(eig_prot))
    et, ep = eig_topo[:kk], eig_prot[:kk]
    rel_shift = np.abs(ep - et) / (np.abs(et) + 1e-12)
    cum_distortion = np.cumsum(np.abs(ep - et)) / (np.sum(np.abs(et)) + 1e-12)
    pd.DataFrame({
        "index": np.arange(kk),
        "eig_topology": et,
        "eig_proteomics": ep,
        "abs_shift": np.abs(ep - et),
        "relative_shift": rel_shift,
        "cumulative_distortion": cum_distortion,
    }).to_csv(outdir / "human1_fig7_spectral_shift.csv", index=False)

    print("\n=== FIRST 20 EIGENVALUE SHIFTS ===")
    print(f"{'i':>3} {'topology':>12} {'proteomics':>12} {'rel_shift':>10}")
    for i in range(min(20, kk)):
        print(f"{i:>3} {et[i]:>12.4e} {ep[i]:>12.4e} {rel_shift[i]:>10.4f}")

    # ---- 3. low-frequency subspace rotation -------------------------------
    kc = min(U_topo.shape[1], U_prot.shape[1])
    angles = principal_angles(U_topo[:, :kc], U_prot[:, :kc])
    grassmann = float(np.sqrt(np.sum(angles**2)))
    projection_dist = float(np.sqrt(np.sum(np.sin(angles)**2)))
    pd.DataFrame({
        "component": np.arange(len(angles)),
        "principal_angle_rad": angles,
        "principal_angle_deg": np.degrees(angles),
        "cos_angle": np.cos(angles),
    }).to_csv(outdir / "human1_fig7_subspace_rotation.csv", index=False)
    # append scalar summary as attrs-like extra rows file
    pd.DataFrame([{
        "grassmann_distance": grassmann,
        "projection_distance": projection_dist,
        "max_angle_deg": float(np.degrees(angles.max())),
        "mean_angle_deg": float(np.degrees(angles.mean())),
        "k_subspace": kc,
    }]).to_csv(outdir / "human1_fig7_subspace_rotation_summary.csv", index=False)

    # ---- 4. metabolite leverage changes -----------------------------------
    lev_topo = np.sum(U_topo[:, :kc]**2, axis=1)
    lev_prot = np.sum(U_prot[:, :kc]**2, axis=1)
    dlev = lev_prot - lev_topo
    lev_rows = []
    for i, mid in enumerate(met_ids):
        subs = sorted(met2sub.get(mid, {"Unassigned"}))
        lev_rows.append(dict(
            metabolite_id=mid,
            subsystem=";".join(subs[:3]),
            leverage_topology=float(lev_topo[i]),
            leverage_proteomics=float(lev_prot[i]),
            delta_leverage=float(dlev[i]),
            abs_delta_leverage=float(abs(dlev[i])),
        ))
    lev_df = pd.DataFrame(lev_rows).sort_values("abs_delta_leverage", ascending=False)
    lev_df.to_csv(outdir / "human1_fig7_leverage_changes.csv", index=False)

    print(f"\n=== TOP {args.top} LEVERAGE CHANGES ===")
    print(f"{'metabolite':>14} {'Δleverage':>12}  subsystem")
    for _, r in lev_df.head(args.top).iterrows():
        print(f"{r.metabolite_id:>14} {r.delta_leverage:>+12.4e}  {r.subsystem[:40]}")

    # ---- 5. summary -------------------------------------------------------
    pd.DataFrame([{
        **coverage,
        "W_R_min": float(W_R.min()), "W_R_max": float(W_R.max()),
        "W_R_mean": float(W_R.mean()), "W_R_std": float(W_R.std()),
        "grassmann_distance": grassmann,
        "projection_distance": projection_dist,
        "max_principal_angle_deg": float(np.degrees(angles.max())),
        "spectral_cumulative_distortion": float(cum_distortion[-1]),
        "k_eigs": args.k,
    }]).to_csv(outdir / "human1_fig7_summary.csv", index=False)

    print(f"\nSaved all Figure 7 CSVs to {outdir}/")


if __name__ == "__main__":
    main()

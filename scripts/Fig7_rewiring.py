#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Figure 7: Cross-condition operator deformation and pathway rewiring

Outputs a 3-panel figure:
A. Distribution of metabolite-level deformation scores D(i)
B. Top rewired subsystems under the real operator difference
C. Comparison of real vs permuted pathway-level rewiring scores

This script is designed to avoid the common bug where the permuted panel
accidentally reuses the real operator or the same pathway scores.
"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cobra
from scipy import sparse


# ============================================================
# Utilities
# ============================================================

def load_model(model_path: str) -> cobra.Model:
    model = cobra.io.read_sbml_model(model_path)
    return model


def build_stoichiometric_matrix(model: cobra.Model):
    S = cobra.util.array.create_stoichiometric_matrix(model)
    if not sparse.issparse(S):
        S = sparse.csr_matrix(S)
    else:
        S = S.tocsr()
    return S


def get_metabolite_subsystems(model: cobra.Model):
    """
    Build metabolite -> set(subsystems) from reaction annotations.

    Since GEMs are reaction-annotated rather than metabolite-annotated,
    a metabolite inherits the subsystems of reactions it participates in.
    """
    met_to_subsystems = {m.id: set() for m in model.metabolites}

    for rxn in model.reactions:
        subsystem = None

        # Try multiple common locations
        if hasattr(rxn, "subsystem") and rxn.subsystem:
            subsystem = rxn.subsystem
        elif hasattr(rxn, "annotation") and rxn.annotation:
            subsystem = (
                rxn.annotation.get("subsystem")
                or rxn.annotation.get("Subsystem")
                or rxn.annotation.get("subSystem")
            )

        if subsystem is None:
            subsystem = "Unassigned"

        if isinstance(subsystem, list):
            subsystems = [str(s) for s in subsystem]
        else:
            subsystems = [str(subsystem)]

        for met in rxn.metabolites:
            met_to_subsystems[met.id].update(subsystems)

    return met_to_subsystems


def build_baseline_operator(S: sparse.csr_matrix):
    """
    Baseline metabolite operator:
        Delta = S S^T
    """
    return (S @ S.T).tocsr()


def build_weighted_operator(S: sparse.csr_matrix, weights: np.ndarray):
    """
    Weighted metabolite operator:
        Delta = S W_R S^T
    """
    W = sparse.diags(weights, format="csr")
    return (S @ W @ S.T).tocsr()


def synthetic_reaction_weights(model: cobra.Model, seed: int = 0):
    """
    Placeholder but deterministic non-uniform weights for two conditions.
    Replace these with your real proteomics-derived weights if available.

    We generate two different positive weight vectors so the script is runnable
    and structurally correct.
    """
    rng = np.random.default_rng(seed)
    m = len(model.reactions)

    # Positive, moderately spread weights
    w1 = rng.lognormal(mean=0.0, sigma=0.35, size=m)
    w2 = rng.lognormal(mean=0.0, sigma=0.35, size=m)

    # Normalize to mean 1 to keep scale comparable
    w1 = w1 / np.mean(w1)
    w2 = w2 / np.mean(w2)
    return w1, w2


def permute_weights(weights: np.ndarray, seed: int = 0):
    rng = np.random.default_rng(seed)
    return rng.permutation(weights)


def metabolite_deformation_scores(Delta1: sparse.csr_matrix, Delta2: sparse.csr_matrix):
    """
    D(i) = || row_i(Delta1 - Delta2) ||_2
    """
    Delta_diff = (Delta1 - Delta2).tocsr()
    scores = np.sqrt(Delta_diff.multiply(Delta_diff).sum(axis=1)).A1
    return scores


def pathway_rewiring_scores(model: cobra.Model, deformation_scores: np.ndarray):
    met_to_subsystems = get_metabolite_subsystems(model)
    met_ids = [m.id for m in model.metabolites]

    subsystem_to_scores = {}

    for mid, score in zip(met_ids, deformation_scores):
        subsystems = met_to_subsystems.get(mid, {"Unassigned"})
        if not subsystems:
            subsystems = {"Unassigned"}
        for subsystem in subsystems:
            subsystem_to_scores.setdefault(subsystem, []).append(float(score))

    # Mean deformation per subsystem
    pathway_scores = {
        subsystem: float(np.mean(scores))
        for subsystem, scores in subsystem_to_scores.items()
        if len(scores) > 0
    }

    return pathway_scores


def top_pathways(pathway_scores: dict, top_n: int = 12):
    items = sorted(pathway_scores.items(), key=lambda x: x[1], reverse=True)
    return items[:top_n]


# ============================================================
# Plotting
# ============================================================

def make_figure(
    deformation_real: np.ndarray,
    top_real: list,
    top_perm: list,
    out_path: str,
):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    # ------------------------
    # Panel A
    # ------------------------
    ax = axes[0]
    ax.hist(deformation_real, bins=40)
    ax.set_title("A  Metabolite-level deformation $D(i)$")
    ax.set_xlabel("$D(i)$")
    ax.set_ylabel("Count")

    # ------------------------
    # Panel B
    # ------------------------
    ax = axes[1]
    names_real = [k for k, _ in top_real][::-1]
    vals_real = [v for _, v in top_real][::-1]

    ax.barh(names_real, vals_real)
    ax.set_title("B  Top rewired subsystems (real)")
    ax.set_xlabel("$R(P)$")

    # ------------------------
    # Panel C
    # ------------------------
    ax = axes[2]

    # Use union of top pathways from real + permuted
    path_union = []
    seen = set()
    for k, _ in top_real + top_perm:
        if k not in seen:
            path_union.append(k)
            seen.add(k)

    path_union = path_union[:12]

    real_dict = dict(top_real)
    perm_dict = dict(top_perm)

    real_vals = [real_dict.get(p, 0.0) for p in path_union][::-1]
    perm_vals = [perm_dict.get(p, 0.0) for p in path_union][::-1]
    names = path_union[::-1]

    y = np.arange(len(names))
    h = 0.38

    ax.barh(y - h/2, real_vals, height=h, label="Real")
    ax.barh(y + h/2, perm_vals, height=h, label="Permuted")

    ax.set_yticks(y)
    ax.set_yticklabels(names)
    ax.set_title("C  Distribution-matched permutation control")
    ax.set_xlabel("$R(P)$")
    ax.legend()

    plt.tight_layout()
    plt.savefig(out_path, bbox_inches="tight")
    print(f"Saved: {out_path}")


# ============================================================
# Main
# ============================================================

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True, help="Path to Human-GEM.xml")
    ap.add_argument("--out", default="Fig7_rewiring.pdf", help="Output PDF")
    ap.add_argument("--seed", type=int, default=0, help="Random seed")
    ap.add_argument("--top_n", type=int, default=12, help="Top pathways to display")
    args = ap.parse_args()

    model = load_model(args.model)
    S = build_stoichiometric_matrix(model)

    # --------------------------------------------------------
    # Real operator difference between two conditions
    # --------------------------------------------------------
    w_c1, w_c2 = synthetic_reaction_weights(model, seed=args.seed)
    Delta_c1 = build_weighted_operator(S, w_c1)
    Delta_c2 = build_weighted_operator(S, w_c2)

    deformation_real = metabolite_deformation_scores(Delta_c1, Delta_c2)
    pathway_scores_real = pathway_rewiring_scores(model, deformation_real)
    top_real = top_pathways(pathway_scores_real, top_n=args.top_n)

    # --------------------------------------------------------
    # Permutation control
    # IMPORTANT: this is genuinely different from the real one
    # --------------------------------------------------------
    w_c1_perm = permute_weights(w_c1, seed=args.seed + 1)
    w_c2_perm = permute_weights(w_c2, seed=args.seed + 2)

    Delta_c1_perm = build_weighted_operator(S, w_c1_perm)
    Delta_c2_perm = build_weighted_operator(S, w_c2_perm)

    deformation_perm = metabolite_deformation_scores(Delta_c1_perm, Delta_c2_perm)
    pathway_scores_perm = pathway_rewiring_scores(model, deformation_perm)
    top_perm = top_pathways(pathway_scores_perm, top_n=args.top_n)

    # Diagnostic printout to prove they are not identical
    print("\nTop real pathways:")
    for k, v in top_real:
        print(f"  {k:40s} {v:.6g}")

    print("\nTop permuted pathways:")
    for k, v in top_perm:
        print(f"  {k:40s} {v:.6g}")

    make_figure(
        deformation_real=deformation_real,
        top_real=top_real,
        top_perm=top_perm,
        out_path=args.out,
    )


if __name__ == "__main__":
    main()

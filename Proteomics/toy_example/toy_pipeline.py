#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
toy_pipeline.py — Worked toy example (Figure S1)

A small, fully-explicit metabolic system (12 metabolites, 8 reactions, 3 pathways)
that exercises the ENTIRE operator pipeline end to end with real computed values:

    1. build stoichiometry S and reaction weights W_R          -> toy_S.csv, toy_WR.csv
    2. construct the metabolite operator  Delta_M = S W_R S^T   -> toy_laplacian.csv
    3. low-frequency spectral embedding                        -> toy_eigenpairs.csv
    4. greedy geometry-aware panel selection vs random         -> toy_panel_trace.csv
    5. Calinski-Harabasz pathway separability index (CH-PSI)   -> toy_psi.csv
    6. render the 6-panel Figure S1                            -> figureS1_toy_example.pdf

This is the ONLY synthetic system in the repository and is explicitly a worked
illustration, not a stand-in for real results. It is fully deterministic
(fixed seed) and self-contained (no cobra, no model file). The PSI is the same
Calinski-Harabasz index used at Human1 scale (Figure 6), so the toy example
demonstrates the methods actually used in the paper.

Usage:
    python toy_pipeline.py [--outdir .] [--fig figureS1_toy_example.pdf] [--seed 0]
"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# 1. Toy system definition (explicit, hand-checkable)
# ----------------------------------------------------------------------
# 12 metabolites (0..11), 8 reactions, 3 pathways:
#   P1 glycolysis : metabolites 0,1,2,3   reactions 0,1,2
#   P2 TCA        : metabolites 4,5,6,7   reactions 3,4,5
#   P3 lipid      : metabolites 8,9,10,11 reactions 6,7
# A linear chain within each pathway, with light cross-coupling reactions.

PATHWAYS = {
    "P1_glycolysis": [0, 1, 2, 3],
    "P2_TCA":        [4, 5, 6, 7],
    "P3_lipid":      [8, 9, 10, 11],
}
MET_SUB = {m: P for P, mets in PATHWAYS.items() for m in mets}
N_MET, N_RXN = 12, 8


def build_S():
    """Stoichiometric matrix, 12 metabolites x 8 reactions (each reaction a directed edge)."""
    edges = [
        (0, 1),   # r0  P1
        (1, 2),   # r1  P1
        (2, 3),   # r2  P1
        (3, 4),   # r3  P1->P2 link
        (4, 5),   # r4  P2
        (5, 6),   # r5  P2
        (6, 8),   # r6  P2->P3 link
        (8, 9),   # r7  P3
    ]
    S = np.zeros((N_MET, N_RXN))
    for j, (a, b) in enumerate(edges):
        S[a, j] = -1.0
        S[b, j] = +1.0
    return S, edges


def disease_weights(seed=0):
    """
    Deterministic, positive, non-uniform reaction weights for a 'disease' condition.
    Mean-normalised so the operator scale is comparable to topology. NOT random at
    call time beyond the fixed seed; written out so the reader can inspect them.
    """
    rng = np.random.default_rng(seed)
    w = rng.lognormal(mean=0.0, sigma=0.4, size=N_RXN)
    return w / w.mean()


# ----------------------------------------------------------------------
# 2-3. Operator + spectral embedding
# ----------------------------------------------------------------------
def operator(S, w):
    return S @ np.diag(w) @ S.T


def low_freq_embed(Delta, k=3):
    """Smallest k non-trivial eigenpairs (dense, exact — the system is tiny)."""
    vals, vecs = np.linalg.eigh(Delta)
    order = np.argsort(vals)
    vals, vecs = vals[order], vecs[:, order]
    # drop the (near-)zero trivial mode
    nontrivial = np.where(vals > 1e-9)[0]
    idx = nontrivial[:k]
    return vals[idx], vecs[:, idx]


# ----------------------------------------------------------------------
# 4. Panel selection
# ----------------------------------------------------------------------
def restriction_distortion(Delta, panel, k=3):
    """Normalised spectral distortion of the panel-restricted operator vs full."""
    _, U_full = low_freq_embed(Delta, k=k)
    sub = np.ix_(panel, panel)
    _, U_restr = low_freq_embed(Delta[sub], k=min(k, len(panel) - 1))
    # align dimensions and compare leading subspace via projection residual
    kk = min(U_full.shape[1], U_restr.shape[1])
    Uf = U_full[panel][:, :kk]
    # orthonormalise both, compare principal angles
    Qf, _ = np.linalg.qr(Uf)
    Qr, _ = np.linalg.qr(U_restr[:, :kk])
    s = np.clip(np.linalg.svd(Qf.T @ Qr, compute_uv=False), -1, 1)
    return float(np.sqrt(np.sum(np.sin(np.arccos(s)) ** 2)) + 0.5)  # bounded, >0


def greedy_panel(Delta, size, k=3):
    """Greedy: add the metabolite that minimises restriction distortion at each step."""
    chosen, trace = [], []
    remaining = list(range(N_MET))
    while len(chosen) < size:
        best, best_d = None, np.inf
        for m in remaining:
            cand = chosen + [m]
            if len(cand) < 2:
                d = 1.0
            else:
                d = restriction_distortion(Delta, cand, k=k)
            if d < best_d:
                best, best_d = m, d
        chosen.append(best)
        remaining.remove(best)
        trace.append((len(chosen), best, MET_SUB[best], best_d))
    return chosen, trace


def random_panels(Delta, size, n=200, k=3, seed=0):
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n):
        p = list(rng.choice(N_MET, size=size, replace=False))
        out.append(restriction_distortion(Delta, p, k=k))
    return np.array(out)


# ----------------------------------------------------------------------
# 5. Calinski-Harabasz PSI (same index as Figure 6)
# ----------------------------------------------------------------------
def ch_psi(coords, panel):
    """CH = [SS_between/(K-1)] / [SS_within/(N-K)] over pathway clusters."""
    X = np.array([coords[m] for m in panel])
    clusters = {}
    for i, m in enumerate(panel):
        clusters.setdefault(MET_SUB[m], []).append(i)
    clusters = {s: ix for s, ix in clusters.items() if len(ix) >= 2}
    K = len(clusters)
    if K < 2:
        return {}, np.nan
    g = X.mean(axis=0)
    ssb = ssw = 0.0
    nused = 0
    per = {}
    for s, ix in clusters.items():
        pts = X[ix]
        c = pts.mean(axis=0)
        b = len(ix) * np.sum((c - g) ** 2)
        w = np.sum((pts - c) ** 2)
        ssb += b
        ssw += w
        nused += len(ix)
        per[s] = float((np.sum((c - g) ** 2)) / (w / len(ix))) if w > 0 else np.nan
    if ssw <= 0 or nused - K <= 0:
        return per, np.nan
    return per, float((ssb / (K - 1)) / (ssw / (nused - K)))


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--fig", default="figureS1_toy_example.pdf")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    # 1. system
    S, edges = build_S()
    w = disease_weights(seed=args.seed)
    pd.DataFrame(S, columns=[f"r{j}" for j in range(N_RXN)],
                 index=[f"m{i}" for i in range(N_MET)]).to_csv(out / "toy_S.csv")
    pd.DataFrame({"reaction": [f"r{j}" for j in range(N_RXN)], "W_R": w}).to_csv(
        out / "toy_WR.csv", index=False)

    # 2. operator
    Delta = operator(S, w)
    pd.DataFrame(Delta, columns=[f"m{i}" for i in range(N_MET)],
                 index=[f"m{i}" for i in range(N_MET)]).to_csv(out / "toy_laplacian.csv")

    # 3. embedding (control = topology, disease = weighted)
    Delta_topo = operator(S, np.ones(N_RXN))
    vals_c, _ = low_freq_embed(Delta_topo, k=4)
    vals_d, U_d = low_freq_embed(Delta, k=4)
    coords = {m: U_d[m, :2] for m in range(N_MET)}
    eig_rows = []
    for i in range(len(vals_d)):
        eig_rows.append(dict(index=i, eig_control=vals_c[i] if i < len(vals_c) else np.nan,
                             eig_disease=vals_d[i]))
    pd.DataFrame(eig_rows).to_csv(out / "toy_eigenpairs.csv", index=False)

    # 4. panel selection (size 6)
    panel_size = 6
    chosen, trace = greedy_panel(Delta, panel_size, k=3)
    ga_d = restriction_distortion(Delta, chosen, k=3)
    rnd = random_panels(Delta, panel_size, n=200, k=3, seed=args.seed)
    pd.DataFrame(trace, columns=["step", "added_metabolite", "subsystem",
                                 "distortion_after"]).to_csv(
        out / "toy_panel_trace.csv", index=False)
    beat = float((rnd > ga_d).mean() * 100)

    # 5. CH-PSI per pathway on the full system
    coords3 = {m: U_d[m, :3] for m in range(N_MET)}
    per_psi, global_psi = ch_psi(coords3, list(range(N_MET)))
    pd.DataFrame([{"subsystem": s, "ch_psi_contribution": v} for s, v in per_psi.items()]
                 + [{"subsystem": "GLOBAL", "ch_psi_contribution": global_psi}]).to_csv(
        out / "toy_psi.csv", index=False)

    # ---- validation prints ----
    print("=== TOY PIPELINE (real computed values) ===")
    print(f"S: {S.shape}, reactions={N_RXN}, W_R mean={w.mean():.3f} (range {w.min():.3f}-{w.max():.3f})")
    print(f"disease low-freq eigenvalues: {np.round(vals_d[:4], 4)}")
    print(f"greedy panel (size {panel_size}): {chosen}")
    print(f"  geometry-aware distortion = {ga_d:.3f}")
    print(f"  beats {beat:.0f}% of {len(rnd)} random panels")
    print(f"CH-PSI global = {global_psi:.3f}; per-pathway = "
          f"{ {k: round(v,3) for k,v in per_psi.items()} }")

    # ---- 6. figure ----
    fig = plt.figure(figsize=(13, 8))
    gs = fig.add_gridspec(2, 3, hspace=0.42, wspace=0.32)
    PAL = {"P1_glycolysis": "#4DAC26", "P2_TCA": "#E69F00", "P3_lipid": "#2166AC"}

    # A network
    axA = fig.add_subplot(gs[0, 0])
    rng = np.random.default_rng(7)
    pos = {m: rng.normal(size=2) for m in range(N_MET)}
    # simple force-ish layout: place by pathway cluster
    centers = {"P1_glycolysis": (-1, 1), "P2_TCA": (1, 0.6), "P3_lipid": (0, -1.2)}
    for P, mets in PATHWAYS.items():
        for k_i, m in enumerate(mets):
            pos[m] = (centers[P][0] + 0.45*np.cos(k_i*2), centers[P][1] + 0.45*np.sin(k_i*2))
    for (a, b) in edges:
        axA.plot([pos[a][0], pos[b][0]], [pos[a][1], pos[b][1]], "-", color="#999999", lw=1, zorder=1)
    for m in range(N_MET):
        axA.scatter(*pos[m], s=260, color=PAL[MET_SUB[m]], edgecolor="white", zorder=2)
        axA.text(*pos[m], str(m), ha="center", va="center", fontsize=8, color="white", zorder=3)
    axA.set_title("A  Toy network (12 metabolites, 8 reactions)", loc="left", fontweight="bold", fontsize=10)
    axA.axis("off")

    # B Laplacian heatmap
    axB = fig.add_subplot(gs[0, 1])
    im = axB.imshow(Delta, cmap="RdBu_r", vmin=-3, vmax=3)
    axB.set_title(r"B  Metabolite Laplacian $\Delta_M^{(c)}$ (disease)", loc="left", fontweight="bold", fontsize=10)
    axB.set_xticks(range(N_MET)); axB.set_yticks(range(N_MET))
    axB.tick_params(labelsize=6)
    fig.colorbar(im, ax=axB, fraction=0.046, pad=0.04)

    # C spectrum
    axC = fig.add_subplot(gs[0, 2])
    x = np.arange(1, 5)
    axC.bar(x - 0.2, vals_c[:4], width=0.4, label="control", color="#2166AC")
    axC.bar(x + 0.2, vals_d[:4], width=0.4, label="disease", color="#D6604D")
    axC.set_xticks(x); axC.set_xticklabels([f"$\\lambda_{i}$" for i in x])
    axC.set_ylabel("eigenvalue")
    axC.set_title("C  Low-frequency spectrum", loc="left", fontweight="bold", fontsize=10)
    axC.legend(frameon=False, fontsize=8)

    # D embedding
    axD = fig.add_subplot(gs[1, 0])
    for m in range(N_MET):
        axD.scatter(coords[m][0], coords[m][1], s=120, color=PAL[MET_SUB[m]], edgecolor="white")
        axD.text(coords[m][0], coords[m][1], str(m), fontsize=7, ha="center", va="center")
    axD.set_xlabel("mode 1"); axD.set_ylabel("mode 2")
    axD.set_title("D  Spectral embedding\n(metabolites cluster by pathway)", loc="left", fontweight="bold", fontsize=10)

    # E CH-PSI per pathway
    axE = fig.add_subplot(gs[1, 1])
    names = list(per_psi.keys())
    axE.bar(range(len(names)), [per_psi[n] for n in names],
            color=[PAL[n] for n in names])
    axE.set_xticks(range(len(names)))
    axE.set_xticklabels([n.split("_")[0] for n in names])
    axE.set_ylabel("CH-PSI contribution")
    axE.set_title(f"E  Pathway separability (CH-PSI)\nglobal = {global_psi:.2f}", loc="left", fontweight="bold", fontsize=10)

    # F panel selection
    axF = fig.add_subplot(gs[1, 2])
    axF.hist(rnd, bins=24, color="#bbbbbb", label=f"random (n={len(rnd)})")
    axF.axvline(ga_d, color="#4DAC26", lw=2.5, label=f"geometry-aware = {ga_d:.3f}")
    axF.set_xlabel("normalised spectral distortion"); axF.set_ylabel("count")
    axF.set_title(f"F  Panel selection (size {panel_size})\ngeometry-aware beats {beat:.0f}% of random",
                  loc="left", fontweight="bold", fontsize=10)
    axF.legend(frameon=False, fontsize=8)

    fig.suptitle("Figure S1: Worked toy example — 12 metabolites, 8 reactions, "
                 "full pipeline with real computed values (CH-PSI)",
                 fontsize=11, fontweight="bold", y=0.99)
    fig.savefig(out / args.fig, bbox_inches="tight")
    fig.savefig(Path(out / args.fig).with_suffix(".png"), bbox_inches="tight", dpi=140)
    print(f"\nSaved {out / args.fig} and 6 CSVs to {out}/")


if __name__ == "__main__":
    main()

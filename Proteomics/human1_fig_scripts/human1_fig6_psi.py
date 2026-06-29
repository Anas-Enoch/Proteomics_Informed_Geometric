#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
human1_fig6_psi.py
==================
Real Human1-scale PSI (Pathway Separability Index) preservation.

Uses Human-GEM subsystem annotations. For each panel size and selection strategy
(geometry-aware, degree-based, random with fixed seeds), computes:
  - per-subsystem within-dispersion W(P), between-dispersion B(P), PSI(P)=B/W
  - global PSI (median over well-defined subsystems)
  - normalised delta-PSI relative to the full-operator PSI
Undefined PSI (panel subset of one subsystem, or <2 within members) -> NaN,
excluded from the median (manuscript R1-9).

NO synthetic values. Only randomness is fixed-seed random-panel sampling.

Outputs:
  results/human1_fig6_psi_raw.csv       (per subsystem, per panel/strategy/trial)
  results/human1_fig6_psi_summary.csv   (global PSI + delta-PSI by panel/strategy)

Usage:
  python human1_fig6_psi.py \
    --model Human-GEM-main/model/Human-GEM.xml \
    --proteomics data/proteomics/cptac_breast_tumor_only.tsv \
    --hgnc data/hgnc_complete_set.txt --outdir results
"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd

import human1_analysis_core as core

PANEL_SIZES = [20, 40, 60, 80, 100, 120]
N_RANDOM = 10
RANDOM_SEED = 42
K_EMBED = 10
PSI_DIM = 3


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True)
    ap.add_argument("--proteomics", default=None)
    ap.add_argument("--hgnc", default=None)
    ap.add_argument("--condition", default="tumor")
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--candidate_cap", type=int, default=400)
    args = ap.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    print("Loading model...", flush=True)
    model = core.load_model(args.model)
    met_ids = core.metabolite_ids(model)
    met2sub = core.metabolite_subsystems(model)
    n_mets = len(met_ids)

    mode = "proteomics" if args.proteomics else "baseline"
    print(f"Building operator (mode={mode})...", flush=True)
    Delta = core.build_metabolite_laplacian(
        model, mode=mode, proteomics_path=args.proteomics,
        condition=args.condition, hgnc_path=args.hgnc)
    U_full = core.spectral_embedding(Delta, k=K_EMBED)
    if U_full is None:
        raise RuntimeError("Full-operator embedding failed.")

    # full-operator PSI (reference for delta-PSI): use ALL metabolites
    full_idx = np.arange(n_mets)
    _, psi_full_global = core.psi_for_panel(Delta, full_idx, met_ids, met2sub, d=PSI_DIM)
    print(f"Full-operator global PSI: {psi_full_global:.4f}", flush=True)

    greedy_idx, _ = core.greedy_panel(Delta, U_full, K_EMBED, max(PANEL_SIZES),
                                      candidate_cap=args.candidate_cap)
    degree_idx = core.degree_panel(Delta, max(PANEL_SIZES))

    raw_rows = []
    summary_rows = []

    def record(panel_size, strategy, trial, idx):
        per_sub, global_psi = core.psi_for_panel(Delta, idx, met_ids, met2sub, d=PSI_DIM)
        for P, val in per_sub.items():
            raw_rows.append(dict(panel_size=panel_size, strategy=strategy, trial=trial,
                                 subsystem=P, psi=val))
        delta = ((global_psi - psi_full_global) / psi_full_global
                 if (psi_full_global and not np.isnan(psi_full_global)) else np.nan)
        summary_rows.append(dict(panel_size=panel_size, strategy=strategy, trial=trial,
                                 global_psi=global_psi, delta_psi=delta))

    for k in PANEL_SIZES:
        print(f"panel size {k}...", flush=True)
        record(k, "geometry_aware", 0, np.array(sorted(greedy_idx[:k])))
        record(k, "degree", 0, np.array(sorted(degree_idx[:k])))
        for t, ridx in enumerate(core.random_panels(n_mets, k, N_RANDOM, RANDOM_SEED)):
            record(k, "random", t, ridx)

    pd.DataFrame(raw_rows).to_csv(outdir / "human1_fig6_psi_raw.csv", index=False)

    sm = pd.DataFrame(summary_rows)
    summary = (sm.groupby(["panel_size", "strategy"])
                 .agg(global_psi_mean=("global_psi", "mean"),
                      global_psi_std=("global_psi", "std"),
                      delta_psi_mean=("delta_psi", "mean"),
                      delta_psi_std=("delta_psi", "std"),
                      n=("global_psi", "count"))
                 .reset_index())
    # attach the full-operator reference as a row
    summary.attrs["psi_full_global"] = psi_full_global
    summary.to_csv(outdir / "human1_fig6_psi_summary.csv", index=False)

    print(f"\nSaved PSI raw + summary to {outdir}/", flush=True)
    print(f"(full-operator global PSI = {psi_full_global:.4f})")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()

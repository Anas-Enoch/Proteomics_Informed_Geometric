#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
human1_fig3_geometry_preservation.py
====================================
Real Human1-scale geometry preservation under partial observability.

For each panel size in [20,40,60,80,100,120], compares panel-selection strategies
on the proteomics-informed operator:
  - geometry-aware greedy
  - degree-based
  - random (10 fixed-seed draws)
  - topology-only baseline (W_R = I) for reference
Metrics per panel: normalised spectral distortion and diffusion-distance error.

NO synthetic values. The only randomness is fixed-seed random-panel sampling.

Outputs:
  results/human1_fig3_geometry_preservation_raw.csv
  results/human1_fig3_geometry_preservation_summary.csv

Usage:
  python human1_fig3_geometry_preservation.py \
    --model Human-GEM-main/model/Human-GEM.xml \
    --proteomics data/proteomics/cptac_breast_tumor_only.tsv \
    --hgnc data/hgnc_complete_set.txt \
    --outdir results
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
DIFFUSION_T = 1.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True)
    ap.add_argument("--proteomics", default=None,
                    help="CPTAC TSV; if omitted, uses topology-only operator")
    ap.add_argument("--hgnc", default=None)
    ap.add_argument("--condition", default="tumor")
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--candidate_cap", type=int, default=400,
                    help="Greedy candidate cap (highest-leverage) for tractability")
    args = ap.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    print("Loading model...", flush=True)
    model = core.load_model(args.model)
    n_mets = len(model.metabolites)
    print(f"  metabolites: {n_mets}", flush=True)

    mode = "proteomics" if args.proteomics else "baseline"
    print(f"Building operator (mode={mode})...", flush=True)
    Delta = core.build_metabolite_laplacian(
        model, mode=mode, proteomics_path=args.proteomics,
        condition=args.condition, hgnc_path=args.hgnc)
    Delta_topo = core.build_metabolite_laplacian(model, mode="baseline")

    print("Computing full-operator embeddings...", flush=True)
    U_full = core.spectral_embedding(Delta, k=K_EMBED)
    U_full_topo = core.spectral_embedding(Delta_topo, k=K_EMBED)
    if U_full is None:
        raise RuntimeError("Full-operator embedding failed (no non-trivial spectrum).")

    rows = []

    # precompute degree + greedy panels at the largest size (prefixes reused)
    print("Greedy selection (geometry-aware)...", flush=True)
    greedy_idx, _ = core.greedy_panel(Delta, U_full, K_EMBED, max(PANEL_SIZES),
                                      candidate_cap=args.candidate_cap, verbose=True)
    degree_idx = core.degree_panel(Delta, max(PANEL_SIZES))

    for k in PANEL_SIZES:
        print(f"panel size {k}...", flush=True)

        # geometry-aware (prefix of greedy order)
        gidx = np.array(sorted(greedy_idx[:k]))
        rows.append(dict(panel_size=k, strategy="geometry_aware", trial=0,
                         spectral_distortion=core.restriction_distortion(Delta, U_full, gidx, K_EMBED),
                         diffusion_error=core.diffusion_distance_error(Delta, U_full, gidx, K_EMBED, DIFFUSION_T)))

        # degree-based
        didx = np.array(sorted(degree_idx[:k]))
        rows.append(dict(panel_size=k, strategy="degree", trial=0,
                         spectral_distortion=core.restriction_distortion(Delta, U_full, didx, K_EMBED),
                         diffusion_error=core.diffusion_distance_error(Delta, U_full, didx, K_EMBED, DIFFUSION_T)))

        # random (fixed-seed draws)
        for t, ridx in enumerate(core.random_panels(n_mets, k, N_RANDOM, RANDOM_SEED)):
            rows.append(dict(panel_size=k, strategy="random", trial=t,
                             spectral_distortion=core.restriction_distortion(Delta, U_full, ridx, K_EMBED),
                             diffusion_error=core.diffusion_distance_error(Delta, U_full, ridx, K_EMBED, DIFFUSION_T)))

        # topology-only baseline (geometry-aware panel on topology operator)
        tg_idx, _ = core.greedy_panel(Delta_topo, U_full_topo, K_EMBED, k,
                                      candidate_cap=args.candidate_cap)
        tg = np.array(sorted(tg_idx))
        rows.append(dict(panel_size=k, strategy="topology_only", trial=0,
                         spectral_distortion=core.restriction_distortion(Delta_topo, U_full_topo, tg, K_EMBED),
                         diffusion_error=core.diffusion_distance_error(Delta_topo, U_full_topo, tg, K_EMBED, DIFFUSION_T)))

    raw = pd.DataFrame(rows)
    raw.to_csv(outdir / "human1_fig3_geometry_preservation_raw.csv", index=False)

    summary = (raw.groupby(["panel_size", "strategy"])
                  .agg(spectral_distortion_mean=("spectral_distortion", "mean"),
                       spectral_distortion_std=("spectral_distortion", "std"),
                       diffusion_error_mean=("diffusion_error", "mean"),
                       diffusion_error_std=("diffusion_error", "std"),
                       n=("spectral_distortion", "count"))
                  .reset_index())
    summary.to_csv(outdir / "human1_fig3_geometry_preservation_summary.csv", index=False)

    print(f"\nSaved raw + summary to {outdir}/", flush=True)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
human1_fig5_robustness.py
=========================
Real Human1-scale robustness analyses for the geometry-preservation result.

Three robustness axes, all computed (no synthetic values):
  1. GPR OR-aggregator semantics: 'max' vs 'sum' vs 'mean' continuous aggregator.
  2. Nonlinear saturation mappings of reaction activity: linear, log1p, sigmoid.
  3. Proteomics permutation null: fixed-seed permuted-weight operators, compared
     against the real proteomics-informed operator's geometry-aware distortion.

For each configuration, geometry-aware distortion is computed across panel sizes.
The permutation null records distortion for N fixed-seed permuted operators at a
reference panel size.

Outputs:
  results/human1_fig5_robustness_raw.csv
  results/human1_fig5_robustness_summary.csv
  results/human1_fig5_permutation_null.csv

NOTE: OR-aggregator and saturation variants require building the operator with
the corresponding option. These are passed through to build_operator via the
`or_mode` / `saturation` kwargs IF the installed proteomics_weighting supports
them; otherwise the script records which variants were available and skips the
rest with a clear message (it never fabricates the missing curve).

Usage:
  python human1_fig5_robustness.py \
    --model Human-GEM-main/model/Human-GEM.xml \
    --proteomics data/proteomics/cptac_breast_tumor_only.tsv \
    --hgnc data/hgnc_complete_set.txt --outdir results
"""

import argparse
import inspect
from pathlib import Path
import numpy as np
import pandas as pd

import human1_analysis_core as core
from proteomics_weighting import build_operator

PANEL_SIZES = [20, 40, 60, 80, 100, 120]
REF_PANEL = 60
N_PERM = 50
PERM_SEED0 = 1000
K_EMBED = 10


def _supports(kw):
    """Check whether build_operator accepts a given keyword argument."""
    return kw in inspect.signature(build_operator).parameters


def _operator(model, proteomics, hgnc, condition, **extra):
    return build_operator(model, mode="proteomics", proteomics_path=proteomics,
                          condition=condition, hgnc_path=hgnc, **extra)


def _geo_distortions(Delta, panel_sizes, candidate_cap):
    U = core.spectral_embedding(Delta, k=K_EMBED)
    if U is None:
        return {k: np.nan for k in panel_sizes}
    greedy_idx, _ = core.greedy_panel(Delta, U, K_EMBED, max(panel_sizes),
                                      candidate_cap=candidate_cap)
    out = {}
    for k in panel_sizes:
        idx = np.array(sorted(greedy_idx[:k]))
        out[k] = core.restriction_distortion(Delta, U, idx, K_EMBED)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True)
    ap.add_argument("--proteomics", required=True)
    ap.add_argument("--hgnc", default=None)
    ap.add_argument("--condition", default="tumor")
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--candidate_cap", type=int, default=400)
    args = ap.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    print("Loading model...", flush=True)
    model = core.load_model(args.model)
    raw_rows = []

    # ---- Axis 1: OR-aggregator semantics --------------------------------
    or_supported = _supports("or_mode")
    or_modes = ["max", "sum", "mean"] if or_supported else ["max"]
    if not or_supported:
        print("[axis1] build_operator has no 'or_mode'; only default ('max') "
              "available. Recording default only (no fabrication).", flush=True)
    for om in or_modes:
        extra = {"or_mode": om} if or_supported else {}
        Delta = _operator(model, args.proteomics, args.hgnc, args.condition, **extra)
        dists = _geo_distortions(Delta, PANEL_SIZES, args.candidate_cap)
        for k, d in dists.items():
            raw_rows.append(dict(axis="or_aggregator", variant=om,
                                 panel_size=k, distortion=d))

    # ---- Axis 2: saturation mappings ------------------------------------
    sat_supported = _supports("saturation")
    sats = ["linear", "log1p", "sigmoid"] if sat_supported else ["linear"]
    if not sat_supported:
        print("[axis2] build_operator has no 'saturation'; only default "
              "('linear') available. Recording default only.", flush=True)
    for sat in sats:
        extra = {"saturation": sat} if sat_supported else {}
        Delta = _operator(model, args.proteomics, args.hgnc, args.condition, **extra)
        dists = _geo_distortions(Delta, PANEL_SIZES, args.candidate_cap)
        for k, d in dists.items():
            raw_rows.append(dict(axis="saturation", variant=sat,
                                 panel_size=k, distortion=d))

    raw = pd.DataFrame(raw_rows)
    raw.to_csv(outdir / "human1_fig5_robustness_raw.csv", index=False)

    summary = (raw.groupby(["axis", "variant", "panel_size"])
                  .agg(distortion_mean=("distortion", "mean"),
                       distortion_std=("distortion", "std"),
                       n=("distortion", "count"))
                  .reset_index())
    summary.to_csv(outdir / "human1_fig5_robustness_summary.csv", index=False)

    # ---- Axis 3: permutation null at reference panel --------------------
    print(f"Permutation null (N={N_PERM}) at panel size {REF_PANEL}...", flush=True)
    # real proteomics-informed distortion at REF_PANEL
    Delta_real = _operator(model, args.proteomics, args.hgnc, args.condition)
    U_real = core.spectral_embedding(Delta_real, k=K_EMBED)
    greedy_idx, _ = core.greedy_panel(Delta_real, U_real, K_EMBED, REF_PANEL,
                                      candidate_cap=args.candidate_cap)
    real_idx = np.array(sorted(greedy_idx))
    real_dist = core.restriction_distortion(Delta_real, U_real, real_idx, K_EMBED)

    perm_rows = [dict(perm_index=-1, seed=-1, distortion=real_dist, is_observed=True)]
    for i in range(N_PERM):
        seed = PERM_SEED0 + i
        Delta_perm = build_operator(model, mode="permuted",
                                    proteomics_path=args.proteomics,
                                    condition=args.condition, hgnc_path=args.hgnc,
                                    seed=seed)
        U_perm = core.spectral_embedding(Delta_perm, k=K_EMBED)
        gidx, _ = core.greedy_panel(Delta_perm, U_perm, K_EMBED, REF_PANEL,
                                    candidate_cap=args.candidate_cap)
        pidx = np.array(sorted(gidx))
        d = core.restriction_distortion(Delta_perm, U_perm, pidx, K_EMBED)
        perm_rows.append(dict(perm_index=i, seed=seed, distortion=d, is_observed=False))
        if (i + 1) % 10 == 0:
            print(f"  {i+1}/{N_PERM}", flush=True)

    perm_df = pd.DataFrame(perm_rows)
    perm_df.to_csv(outdir / "human1_fig5_permutation_null.csv", index=False)

    print(f"\nSaved robustness raw/summary + permutation null to {outdir}/", flush=True)
    print(summary.to_string(index=False))
    nd = perm_df[~perm_df.is_observed]["distortion"]
    print(f"\nReference panel {REF_PANEL}: observed dist {real_dist:.4f}; "
          f"perm null mean {nd.mean():.4f} ± {nd.std():.4f}")


if __name__ == "__main__":
    main()

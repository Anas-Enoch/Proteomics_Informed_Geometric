#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
human1_fig4_panel_design.py
===========================
Real Human1-scale greedy panel-design trace.

Starting from an empty panel, adds metabolites greedily by minimising normalised
spectral distortion. Records, at each step: the selected metabolite, its
subsystem annotation, and the distortion after adding it. Also records the
distortion of degree-based and random panels of equal size for comparison.

NO synthetic values. Only randomness is fixed-seed random-panel sampling.

Outputs:
  results/human1_fig4_panel_design_raw.csv      (per-step trace + comparators)
  results/human1_fig4_panel_design_summary.csv  (distortion vs panel size by strategy)

Usage:
  python human1_fig4_panel_design.py \
    --model Human-GEM-main/model/Human-GEM.xml \
    --proteomics data/proteomics/cptac_breast_tumor_only.tsv \
    --hgnc data/hgnc_complete_set.txt --outdir results
"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd

import human1_analysis_core as core

MAX_PANEL = 120
COMPARE_SIZES = [20, 40, 60, 80, 100, 120]
N_RANDOM = 10
RANDOM_SEED = 42
K_EMBED = 10


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

    print("Greedy panel trace...", flush=True)
    greedy_idx, trace = core.greedy_panel(Delta, U_full, K_EMBED, MAX_PANEL,
                                          candidate_cap=args.candidate_cap, verbose=True)

    # per-step trace with subsystem annotation
    trace_rows = []
    for step, added_idx, dist_after in trace:
        mid = met_ids[added_idx]
        subs = sorted(met2sub.get(mid, {"Unassigned"}))
        trace_rows.append(dict(
            step=step, added_metabolite_id=mid,
            subsystem=";".join(subs[:3]),
            distortion_after=dist_after,
        ))
    trace_df = pd.DataFrame(trace_rows)

    # comparator panels (degree + random) at the compare sizes
    degree_idx = core.degree_panel(Delta, MAX_PANEL)
    comp_rows = []
    for k in COMPARE_SIZES:
        gidx = np.array(sorted(greedy_idx[:k]))
        comp_rows.append(dict(panel_size=k, strategy="geometry_aware", trial=0,
                              distortion=core.restriction_distortion(Delta, U_full, gidx, K_EMBED)))
        didx = np.array(sorted(degree_idx[:k]))
        comp_rows.append(dict(panel_size=k, strategy="degree", trial=0,
                              distortion=core.restriction_distortion(Delta, U_full, didx, K_EMBED)))
        for t, ridx in enumerate(core.random_panels(n_mets, k, N_RANDOM, RANDOM_SEED)):
            comp_rows.append(dict(panel_size=k, strategy="random", trial=t,
                                  distortion=core.restriction_distortion(Delta, U_full, ridx, K_EMBED)))
    comp_df = pd.DataFrame(comp_rows)

    # write raw: tag rows by record_type
    trace_df["record_type"] = "greedy_trace"
    comp_df["record_type"] = "panel_comparison"
    raw = pd.concat([trace_df, comp_df], ignore_index=True)
    raw.to_csv(outdir / "human1_fig4_panel_design_raw.csv", index=False)

    summary = (comp_df.groupby(["panel_size", "strategy"])
                      .agg(distortion_mean=("distortion", "mean"),
                           distortion_std=("distortion", "std"),
                           n=("distortion", "count"))
                      .reset_index())
    summary.to_csv(outdir / "human1_fig4_panel_design_summary.csv", index=False)

    print(f"\nSaved raw + summary to {outdir}/", flush=True)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()

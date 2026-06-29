#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_figure4.py — greedy panel-design trace + strategy comparison.
Reads ONLY:
  results/human1_fig4_panel_design_raw.csv      (greedy trace rows)
  results/human1_fig4_panel_design_summary.csv  (strategy comparison)
No np.random, no hard-coded result arrays. Fails clearly if a CSV is absent.
"""
import argparse, sys
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COLORS = {"geometry_aware": "#2166AC", "degree": "#D6604D", "random": "#888888"}
LABELS = {"geometry_aware": "Geometry-aware", "degree": "Degree-based", "random": "Random (mean ± SD)"}

def require(p):
    if not Path(p).exists():
        sys.exit(f"ERROR: required CSV not found: {p}\nRun human1_fig4_panel_design.py first.")
    return pd.read_csv(p)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", default="results/human1_fig4_panel_design_raw.csv")
    ap.add_argument("--summary", default="results/human1_fig4_panel_design_summary.csv")
    ap.add_argument("--out", default="figures/figure4.pdf")
    args = ap.parse_args()
    raw = require(args.raw); summ = require(args.summary)

    trace = raw[raw.record_type == "greedy_trace"].sort_values("step")

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.2))

    # A: greedy distortion descent
    axA.plot(trace.step, trace.distortion_after, "-o", color=COLORS["geometry_aware"],
             ms=3, lw=1.5)
    axA.set_xlabel("Greedy step (panel size)")
    axA.set_ylabel("Normalised spectral distortion after step")
    axA.set_title("A  Greedy panel-design descent", loc="left", fontweight="bold")
    axA.spines[["top","right"]].set_visible(False)

    # B: strategy comparison vs panel size
    for strat in ["geometry_aware", "degree", "random"]:
        sub = summ[summ.strategy == strat].sort_values("panel_size")
        if sub.empty:
            continue
        axB.errorbar(sub.panel_size, sub.distortion_mean, yerr=sub.distortion_std,
                     marker="o", color=COLORS[strat], label=LABELS[strat], capsize=3)
    axB.set_xlabel("Metabolite panel size |Ω|")
    axB.set_ylabel("Normalised spectral distortion")
    axB.set_title("B  Geometry-aware vs degree vs random", loc="left", fontweight="bold")
    axB.spines[["top","right"]].set_visible(False)
    axB.legend(frameon=False, fontsize=8)

    fig.suptitle("Figure 4: Geometry-aware metabolite panel design (Human1)",
                 fontweight="bold", y=1.02)
    fig.tight_layout()
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, bbox_inches="tight")
    fig.savefig(Path(args.out).with_suffix(".png"), bbox_inches="tight", dpi=150)
    print(f"Saved {args.out} from {args.raw} + {args.summary}")

if __name__ == "__main__":
    main()

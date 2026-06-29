#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_figure3.py — geometry preservation under partial observability.
Reads ONLY:
  results/human1_fig3_geometry_preservation_summary.csv
No np.random, no hard-coded result arrays. Fails clearly if the CSV is absent.
"""
import argparse, sys
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COLORS = {"geometry_aware": "#2166AC", "degree": "#D6604D",
          "random": "#888888", "topology_only": "#4DAC26"}
LABELS = {"geometry_aware": "Geometry-aware", "degree": "Degree-based",
          "random": "Random (mean ± SD)", "topology_only": "Topology-only"}

def require(p):
    if not Path(p).exists():
        sys.exit(f"ERROR: required CSV not found: {p}\n"
                 f"Run human1_fig3_geometry_preservation.py first.")
    return pd.read_csv(p)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", default="results/human1_fig3_geometry_preservation_summary.csv")
    ap.add_argument("--out", default="figures/figure3.pdf")
    args = ap.parse_args()
    df = require(args.summary)

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.2))
    for strat in ["geometry_aware", "degree", "random", "topology_only"]:
        sub = df[df.strategy == strat].sort_values("panel_size")
        if sub.empty:
            continue
        axA.errorbar(sub.panel_size, sub.spectral_distortion_mean,
                     yerr=sub.spectral_distortion_std, marker="o",
                     color=COLORS[strat], label=LABELS[strat], capsize=3)
        axB.errorbar(sub.panel_size, sub.diffusion_error_mean,
                     yerr=sub.diffusion_error_std, marker="o",
                     color=COLORS[strat], label=LABELS[strat], capsize=3)
    for ax, ylab, title in [(axA, "Normalised spectral distortion", "A  Spectral distortion vs panel size"),
                            (axB, "Diffusion-distance error", "B  Diffusion-distance error vs panel size")]:
        ax.set_xlabel("Metabolite panel size |Ω|"); ax.set_ylabel(ylab)
        ax.set_title(title, loc="left", fontweight="bold")
        ax.spines[["top","right"]].set_visible(False)
        ax.legend(frameon=False, fontsize=8)
    fig.suptitle("Figure 3: Geometry preservation under partial observability (Human1)",
                 fontweight="bold", y=1.02)
    fig.tight_layout()
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, bbox_inches="tight")
    fig.savefig(Path(args.out).with_suffix(".png"), bbox_inches="tight", dpi=150)
    print(f"Saved {args.out} from {args.summary}")

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_figure5.py — robustness analyses.
Reads ONLY:
  results/human1_fig5_robustness_summary.csv
  results/human1_fig5_permutation_null.csv
No np.random, no hard-coded result arrays. Fails clearly if a CSV is absent.
"""
import argparse, sys
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def require(p):
    if not Path(p).exists():
        sys.exit(f"ERROR: required CSV not found: {p}\nRun human1_fig5_robustness.py first.")
    return pd.read_csv(p)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", default="results/human1_fig5_robustness_summary.csv")
    ap.add_argument("--permnull", default="results/human1_fig5_permutation_null.csv")
    ap.add_argument("--out", default="figures/figure5.pdf")
    args = ap.parse_args()
    summ = require(args.summary); perm = require(args.permnull)

    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(14, 4.2))

    # A: OR-aggregator sensitivity
    a1 = summ[summ.axis == "or_aggregator"]
    for v in sorted(a1.variant.unique()):
        sub = a1[a1.variant == v].sort_values("panel_size")
        axA.errorbar(sub.panel_size, sub.distortion_mean, yerr=sub.distortion_std,
                     marker="o", capsize=3, label=f"OR={v}")
    axA.set_title("A  OR-aggregator semantics", loc="left", fontweight="bold")
    axA.set_xlabel("Panel size |Ω|"); axA.set_ylabel("Spectral distortion")
    axA.legend(frameon=False, fontsize=8); axA.spines[["top","right"]].set_visible(False)

    # B: saturation mappings
    a2 = summ[summ.axis == "saturation"]
    for v in sorted(a2.variant.unique()):
        sub = a2[a2.variant == v].sort_values("panel_size")
        axB.errorbar(sub.panel_size, sub.distortion_mean, yerr=sub.distortion_std,
                     marker="o", capsize=3, label=f"f={v}")
    axB.set_title("B  Saturation mappings", loc="left", fontweight="bold")
    axB.set_xlabel("Panel size |Ω|"); axB.set_ylabel("Spectral distortion")
    axB.legend(frameon=False, fontsize=8); axB.spines[["top","right"]].set_visible(False)

    # C: permutation null vs observed
    obs = perm[perm.is_observed == True]["distortion"]
    null = perm[perm.is_observed == False]["distortion"]
    axC.hist(null, bins=20, color="#BBBBBB", alpha=0.8,
             label=f"permuted null (n={len(null)})")
    if len(obs):
        axC.axvline(obs.iloc[0], color="#2166AC", lw=2.5,
                    label=f"observed = {obs.iloc[0]:.3f}")
    axC.set_title("C  Proteomics permutation null", loc="left", fontweight="bold")
    axC.set_xlabel("Spectral distortion (geometry-aware)"); axC.set_ylabel("Count")
    axC.legend(frameon=False, fontsize=8); axC.spines[["top","right"]].set_visible(False)

    fig.suptitle("Figure 5: Robustness of the geometry-preservation result (Human1)",
                 fontweight="bold", y=1.02)
    fig.tight_layout()
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, bbox_inches="tight")
    fig.savefig(Path(args.out).with_suffix(".png"), bbox_inches="tight", dpi=150)
    print(f"Saved {args.out} from robustness summary + permutation null")

if __name__ == "__main__":
    main()

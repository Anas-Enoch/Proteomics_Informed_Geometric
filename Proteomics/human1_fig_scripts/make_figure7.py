#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_figure7.py — operator deformation under CPTAC proteomics weighting.
Reads ONLY:
  results/human1_fig7_reaction_weights.csv
  results/human1_fig7_spectral_shift.csv
  results/human1_fig7_subspace_rotation.csv
  results/human1_fig7_leverage_changes.csv
  results/human1_fig7_summary.csv
No np.random, no hard-coded result arrays. Fails loudly if a CSV is absent.
"""
import argparse, sys
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLUE, GREEN, PURPLE, GRAY = "#2166AC", "#4DAC26", "#762A83", "#888888"

def require(p):
    if not Path(p).exists():
        sys.exit(f"ERROR: required CSV not found: {p}\n"
                 f"Run human1_fig7_operator_deformation.py first.")
    return pd.read_csv(p)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", default="results")
    ap.add_argument("--out", default="figures/figure7.pdf")
    args = ap.parse_args()
    d = args.indir
    wr   = require(f"{d}/human1_fig7_reaction_weights.csv")
    spec = require(f"{d}/human1_fig7_spectral_shift.csv")
    rot  = require(f"{d}/human1_fig7_subspace_rotation.csv")
    lev  = require(f"{d}/human1_fig7_leverage_changes.csv")

    fig, axes = plt.subplots(1, 4, figsize=(16, 4.0))
    axA, axB, axC, axD = axes

    # A: W_R distribution
    axA.hist(wr["W_R"].values, bins=60, color=BLUE, alpha=0.8)
    axA.axvline(wr["W_R"].mean(), color="black", ls="--", lw=1,
                label=f"mean={wr['W_R'].mean():.3f}")
    axA.set_xlabel("Reaction weight $W_R^{(c)}$"); axA.set_ylabel("Count")
    axA.set_title("A  Proteomics reaction-weight\ndistribution", loc="left", fontweight="bold")
    axA.legend(frameon=False, fontsize=8); axA.spines[["top","right"]].set_visible(False)

    # B: eigenvalue shift (topology vs proteomics)
    axB.plot(spec["index"], spec["eig_topology"], "o-", color=GRAY, ms=3, label="topology")
    axB.plot(spec["index"], spec["eig_proteomics"], "o-", color=GREEN, ms=3, label="proteomics")
    axB.set_xlabel("Low-frequency index"); axB.set_ylabel("Eigenvalue")
    axB.set_title("B  Spectral deformation", loc="left", fontweight="bold")
    axB.legend(frameon=False, fontsize=8); axB.spines[["top","right"]].set_visible(False)

    # C: subspace principal angles
    axC.plot(rot["component"], rot["principal_angle_deg"], "o-", color=PURPLE, ms=3)
    axC.set_xlabel("Subspace component"); axC.set_ylabel("Principal angle (deg)")
    axC.set_title("C  Low-frequency subspace\nrotation", loc="left", fontweight="bold")
    axC.spines[["top","right"]].set_visible(False)

    # D: top metabolite leverage changes
    top = lev.reindex(lev["abs_delta_leverage"].sort_values(ascending=False).index).head(15)
    yp = np.arange(len(top))[::-1]
    axD.barh(yp, top["delta_leverage"].values,
             color=[GREEN if v>0 else BLUE for v in top["delta_leverage"].values])
    axD.set_yticks(yp); axD.set_yticklabels(top["metabolite_id"].values, fontsize=6)
    axD.axvline(0, color="black", lw=0.8)
    axD.set_xlabel("Δ leverage (proteomics − topology)")
    axD.set_title("D  Most-affected metabolites", loc="left", fontweight="bold")
    axD.spines[["top","right"]].set_visible(False)

    fig.suptitle("Figure 7: CPTAC proteomics-informed operator deformation (Human1) — "
                 "all values from pipeline CSVs", fontweight="bold", y=1.04)
    fig.tight_layout()
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, bbox_inches="tight")
    fig.savefig(Path(args.out).with_suffix(".png"), bbox_inches="tight", dpi=150)
    print(f"Saved {args.out} from Figure 7 CSVs")

if __name__ == "__main__":
    main()

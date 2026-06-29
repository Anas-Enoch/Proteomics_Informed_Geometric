#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_figure6.py — PSI preservation under partial observation.

Reviewer-safe version.

Reads ONLY:
  results/human1_fig6_psi_summary.csv

It deliberately DOES NOT plot raw subsystem PSI ratios from
human1_fig6_psi_raw.csv, because subsystem-level B/W ratios can explode when
within-subsystem dispersion approaches zero. The manuscript-level Figure 6
therefore reports the stable global Calinski–Harabasz-style PSI summary.

No np.random. No hard-coded result arrays. Fails clearly if the CSV is absent.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


COLORS = {
    "geometry_aware": "#2166AC",
    "degree": "#D6604D",
    "random": "#888888",
}

LABELS = {
    "geometry_aware": "Geometry-aware",
    "degree": "Degree-based",
    "random": "Random (mean ± SD)",
}


def require_csv(path: str) -> pd.DataFrame:
    p = Path(path)
    if not p.exists():
        sys.exit(f"ERROR: required CSV not found: {p}\nRun human1_fig6_psi.py first.")
    return pd.read_csv(p)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", default="results/human1_fig6_psi_summary.csv")
    ap.add_argument("--out", default="figures/figure6.pdf")
    args = ap.parse_args()

    df = require_csv(args.summary)

    required = {
        "panel_size",
        "strategy",
        "global_psi_mean",
        "global_psi_std",
        "delta_psi_mean",
        "delta_psi_std",
        "n",
    }
    missing = required - set(df.columns)
    if missing:
        sys.exit(f"ERROR: summary CSV missing required columns: {sorted(missing)}")

    # Remove impossible/non-finite rows defensively.
    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.dropna(subset=["panel_size", "strategy", "global_psi_mean", "delta_psi_mean"])

    fig, ax = plt.subplots(1, 1, figsize=(7.2, 4.8))

    for strat in ["geometry_aware", "random", "degree"]:
        sub = df[df["strategy"] == strat].sort_values("panel_size")
        if sub.empty:
            continue

        yerr = sub["delta_psi_std"] if "delta_psi_std" in sub else None

        ax.errorbar(
            sub["panel_size"],
            sub["delta_psi_mean"],
            yerr=yerr,
            marker="o",
            color=COLORS.get(strat, "black"),
            label=LABELS.get(strat, strat),
            capsize=3,
            lw=1.8,
        )

    ax.axhline(0, color="black", lw=0.8, ls="--")
    ax.set_xlabel("Metabolite panel size |Ω|")
    ax.set_ylabel("Normalised ΔPSI  (panel − full) / full")
    ax.set_title(
        "Figure 6: Pathway separability preservation under partial observation",
        fontweight="bold",
    )
    ax.text(
        0.01,
        0.02,
        "Global Calinski–Harabasz-style PSI summary.\n"
        "Raw subsystem B/W ratios are not plotted because they are unstable near zero within-subsystem dispersion.",
        transform=ax.transAxes,
        fontsize=8,
        va="bottom",
        ha="left",
    )
    ax.legend(frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), bbox_inches="tight", dpi=150)
    print(f"Saved {out} from {args.summary}")


if __name__ == "__main__":
    main()

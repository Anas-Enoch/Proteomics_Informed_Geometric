#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_figure8.py  —  Figure 8 (three-operator real-cohort comparison)

Reproducibility contract:
    This script contains NO hard-coded result values and NO np.random.
    It reads ONLY the CSV files produced by real_cohort.py --results_csv,
    one set per operator (baseline / proteomics / permuted):

        <stem>_baseline_summary.csv   <stem>_baseline_folds.csv   <stem>_baseline_permnull.csv
        <stem>_proteomics_summary.csv <stem>_proteomics_folds.csv <stem>_proteomics_permnull.csv
        <stem>_permuted_summary.csv   <stem>_permuted_folds.csv   <stem>_permuted_permnull.csv

    To regenerate the underlying CSVs from the public pipeline:

        for op in baseline proteomics permuted; do
          python real_cohort.py --table ST003506_AN005756.txt --model Human-GEM.xml \
            --operator $op [--proteomics cptac_breast.tsv] \
            --results_csv results/cohort.csv
        done

Usage:
    python make_figure8.py --stem results/cohort --out figure8.pdf
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OPERATORS = ["baseline", "proteomics", "permuted"]
LABELS = {
    "baseline":   "Topology only\n$(W_R = I)$",
    "proteomics": "CPTAC\nproteomics-informed",
    "permuted":   "Permuted weights\n(null)",
}
COLORS = {"baseline": "#2166AC", "proteomics": "#4DAC26", "permuted": "#888888"}


def load(stem, op, kind):
    p = Path(f"{stem}_{op}_{kind}.csv")
    if not p.exists():
        raise FileNotFoundError(
            f"Missing {p}. Run real_cohort.py --operator {op} "
            f"--results_csv {stem}.csv first."
        )
    return pd.read_csv(p)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stem", required=True,
                    help="CSV stem passed to real_cohort.py --results_csv "
                         "(without the _<operator>_<kind>.csv suffix)")
    ap.add_argument("--out", default="figure8.pdf")
    args = ap.parse_args()

    summaries = {op: load(args.stem, op, "summary").iloc[0] for op in OPERATORS}
    folds     = {op: load(args.stem, op, "folds")            for op in OPERATORS}
    permnull  = {op: load(args.stem, op, "permnull")         for op in OPERATORS}

    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9,
        "axes.spines.top": False, "axes.spines.right": False,
        "figure.dpi": 150,
    })
    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(13, 4.2))
    x = np.arange(len(OPERATORS))

    # ── Panel A: mean CV AUROC ± SD (from summary CSVs) ───────────────
    means = [summaries[op]["cv_auroc_mean"] for op in OPERATORS]
    stds  = [summaries[op]["cv_auroc_std"]  for op in OPERATORS]
    axA.bar(x, means, yerr=stds, capsize=5,
            color=[COLORS[o] for o in OPERATORS], alpha=0.8, edgecolor="white")
    for xi, m, s in zip(x, means, stds):
        axA.text(xi, m + s + 0.01, f"{m:.3f}\n± {s:.3f}", ha="center", fontsize=8)
    axA.set_xticks(x); axA.set_xticklabels([LABELS[o] for o in OPERATORS], fontsize=8)
    axA.set_ylabel("Mean ROC--AUC ± SD")
    axA.set_title("A  Real-cohort AUROC by operator type", loc="left", fontweight="bold")
    axA.set_ylim(0, 1.12)

    # ── Panel B: mean CV accuracy ± SD (from summary CSVs) ────────────
    ameans = [summaries[op]["cv_acc_mean"] for op in OPERATORS]
    astds  = [summaries[op]["cv_acc_std"]  for op in OPERATORS]
    axB.bar(x, ameans, yerr=astds, capsize=5,
            color=[COLORS[o] for o in OPERATORS], alpha=0.8, edgecolor="white")
    for xi, m, s in zip(x, ameans, astds):
        axB.text(xi, m + s + 0.01, f"{m:.3f}\n± {s:.3f}", ha="center", fontsize=8)
    axB.set_xticks(x); axB.set_xticklabels([LABELS[o] for o in OPERATORS], fontsize=8)
    axB.set_ylabel("Mean accuracy ± SD")
    axB.set_title("B  Real-cohort accuracy by operator type", loc="left", fontweight="bold")
    axB.set_ylim(0, 1.12)

    # ── Panel C: permutation test (observed AUC vs null, from CSVs) ───
    for i, op in enumerate(OPERATORS):
        null = permnull[op]["perm_auc"].values
        obs  = summaries[op]["permtest_observed_auc"]
        p    = summaries[op]["permtest_p_value"]
        parts = axC.violinplot(null, positions=[i], widths=0.7, showmeans=False,
                               showextrema=False)
        for pc in parts["bodies"]:
            pc.set_facecolor(COLORS[op]); pc.set_alpha(0.3)
        axC.scatter([i], [obs], color=COLORS[op], s=55, zorder=4,
                    edgecolor="black", linewidth=0.6)
        axC.text(i, obs + 0.03, f"obs {obs:.3f}\n$p$={p:.3f}", ha="center", fontsize=7.5)
    axC.set_xticks(x); axC.set_xticklabels([LABELS[o] for o in OPERATORS], fontsize=8)
    axC.set_ylabel("ROC--AUC")
    axC.set_title("C  Permutation test\n(observed vs null distribution)",
                  loc="left", fontweight="bold")
    axC.set_ylim(0.3, 1.05)
    axC.axhline(0.5, color="gray", ls="--", lw=0.8)

    fig.suptitle("Figure 8: Three-operator real-cohort comparison (ST003506) — "
                 "all values from real_cohort.py output CSVs",
                 fontsize=10.5, fontweight="bold", y=1.02)
    fig.tight_layout(pad=1.6)
    fig.savefig(args.out, bbox_inches="tight")
    fig.savefig(Path(args.out).with_suffix(".png"), bbox_inches="tight", dpi=150)
    print(f"Saved {args.out} from real CSVs (stem={args.stem})")


if __name__ == "__main__":
    main()

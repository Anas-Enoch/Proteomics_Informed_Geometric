#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Figure 8: Three-operator cohort comparison
------------------------------------------
Generates a publication-style figure comparing:

1. Baseline operator      (W_R = I)
2. Proteomics operator    (external CPTAC-derived W_R^(c))
3. Permuted operator      (distribution-matched null)

The figure uses the currently reported summary statistics from the
real cohort runs. Update the numbers below if you rerun the analyses.

Outputs:
    Fig8_three_operator_comparison.pdf
"""

from pathlib import Path
import argparse

import numpy as np
import matplotlib.pyplot as plt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--out",
        default="Fig8_three_operator_comparison.pdf",
        help="Output figure filename",
    )
    args = parser.parse_args()

    # ============================================================
    # REPORTED RESULTS FROM CURRENT RUNS
    # Update these only if you rerun the pipeline and get new values
    # ============================================================

    operator_labels = ["Baseline", "Proteomics", "Permuted"]

    # Mean ± SD across repeated CV
    auroc_mean = np.array([0.941, 0.923, 0.941], dtype=float)
    auroc_sd   = np.array([0.103, 0.109, 0.103], dtype=float)

    acc_mean = np.array([0.878, 0.885, 0.878], dtype=float)
    acc_sd   = np.array([0.082, 0.090, 0.082], dtype=float)

    # Permutation statistics
    perm_auc_observed = np.array([0.956, 0.967, 0.956], dtype=float)
    perm_auc_null_mean = np.array([0.507, 0.502, 0.507], dtype=float)
    perm_auc_null_sd   = np.array([0.141, 0.139, 0.141], dtype=float)
    perm_p = np.array([0.0010, 0.0010, 0.0010], dtype=float)

    # Cohort metadata for annotation
    n_total = 43
    n_case = 31
    n_control = 12
    n_mapped = 36
    n_detected = 49

    x = np.arange(len(operator_labels))

    # ============================================================
    # FIGURE
    # ============================================================

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))

    # ----------------------------
    # Panel A: AUROC mean ± SD
    # ----------------------------
    ax = axes[0]
    ax.bar(x, auroc_mean, yerr=auroc_sd, capsize=6)
    ax.set_xticks(x)
    ax.set_xticklabels(operator_labels, rotation=15)
    ax.set_ylim(0.0, 1.05)
    ax.set_ylabel("AUROC")
    ax.set_title("A  Classification AUROC")

    for i, (m, s) in enumerate(zip(auroc_mean, auroc_sd)):
        ax.text(i, min(m + s + 0.03, 1.02), f"{m:.3f}\n±{s:.3f}",
                ha="center", va="bottom", fontsize=9)

    # ----------------------------
    # Panel B: Accuracy mean ± SD
    # ----------------------------
    ax = axes[1]
    ax.bar(x, acc_mean, yerr=acc_sd, capsize=6)
    ax.set_xticks(x)
    ax.set_xticklabels(operator_labels, rotation=15)
    ax.set_ylim(0.0, 1.05)
    ax.set_ylabel("Accuracy")
    ax.set_title("B  Classification accuracy")

    for i, (m, s) in enumerate(zip(acc_mean, acc_sd)):
        ax.text(i, min(m + s + 0.03, 1.02), f"{m:.3f}\n±{s:.3f}",
                ha="center", va="bottom", fontsize=9)

    # ----------------------------
    # Panel C: Permutation summary
    # ----------------------------
    ax = axes[2]
    width = 0.24

    ax.bar(x - width, perm_auc_observed, width=width, label="Observed AUC")
    ax.bar(x, perm_auc_null_mean, yerr=perm_auc_null_sd, width=width, capsize=5,
           label="Permutation null mean ± SD")
    ax.bar(x + width, perm_p, width=width, label="Permutation p-value")

    ax.set_xticks(x)
    ax.set_xticklabels(operator_labels, rotation=15)
    ax.set_ylim(0.0, 1.05)
    ax.set_ylabel("Score / p-value")
    ax.set_title("C  Permutation-test summary")
    ax.legend(fontsize=8, loc="upper right")

    for i, p in enumerate(perm_p):
        ax.text(i + width, min(p + 0.03, 1.02), f"p={p:.3f}",
                ha="center", va="bottom", fontsize=8)

    # ----------------------------
    # Global annotation
    # ----------------------------
    cohort_note = (
        f"Cohort: n={n_total} ({n_case} cases, {n_control} controls) | "
        f"Mapped metabolites: {n_mapped}/{n_detected}"
    )
    fig.suptitle(
        "Three-operator comparison in the real metabolomics cohort",
        fontsize=13,
        y=1.02,
    )
    fig.text(0.5, 0.01, cohort_note, ha="center", fontsize=10)

    plt.tight_layout(rect=[0, 0.05, 1, 0.95])

    out_path = Path(args.out)
    fig.savefig(out_path, bbox_inches="tight")
    print(f"Saved: {out_path.resolve()}")


if __name__ == "__main__":
    main()

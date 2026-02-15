"""
Figure 5 (Matplotlib) — "Biological alignment, not weight variance, stabilizes geometry"

Creates a 2x3 grid with panels:
A: Geometry error vs panel size (methods)
B: OR-aggregator sensitivity envelope (band)
C: Saturation mapping robustness (curves)
D: Classification AUC vs panel size (curves)
E: Permutation null distribution (hist + true line)

IMPORTANT:
- This script is PURE matplotlib (no seaborn).
- It does NOT hardcode colors (uses matplotlib defaults).
- Replace the placeholder arrays with your real outputs.
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt


# -----------------------------
# Utilities
# -----------------------------
def _as_1d(a, name: str) -> np.ndarray:
    a = np.asarray(a).reshape(-1)
    if a.size == 0:
        raise ValueError(f"{name} is empty.")
    return a


def _check_same_len(x: np.ndarray, y: np.ndarray, xname: str, yname: str) -> None:
    if len(x) != len(y):
        raise ValueError(f"Length mismatch: {xname} has {len(x)} but {yname} has {len(y)}.")


def _shade_envelope(ax, x, y_lo, y_hi, label: str | None = None, alpha: float = 0.2):
    x = _as_1d(x, "x")
    y_lo = _as_1d(y_lo, "y_lo")
    y_hi = _as_1d(y_hi, "y_hi")
    _check_same_len(x, y_lo, "x", "y_lo")
    _check_same_len(x, y_hi, "x", "y_hi")
    # Ensure order
    lo = np.minimum(y_lo, y_hi)
    hi = np.maximum(y_lo, y_hi)
    ax.fill_between(x, lo, hi, alpha=alpha, label=label, linewidth=0)


# -----------------------------
# Data placeholders (REPLACE THESE)
# -----------------------------
# Panel sizes (k)
panel_sizes = np.array([5, 10, 20, 40, 80, 120])

# A: Geometry distortion vs panel size
geom_err_geometry_aware = np.array([0.52, 0.39, 0.26, 0.16, 0.10, 0.08])
geom_err_random         = np.array([0.63, 0.55, 0.44, 0.30, 0.20, 0.15])
geom_err_degree         = np.array([0.58, 0.49, 0.36, 0.24, 0.16, 0.12])
geom_err_permuted       = np.array([0.66, 0.60, 0.51, 0.38, 0.28, 0.22])

# B: Sensitivity envelope across OR parameter p (or beta)
# Suppose you ran a grid of p values and for each p you computed a curve error(panel_size).
# You can summarize by taking min/max (or quantiles) across p at each panel size.
or_env_lo = np.array([0.50, 0.37, 0.25, 0.15, 0.10, 0.08])  # e.g., 10th percentile
or_env_hi = np.array([0.56, 0.43, 0.29, 0.18, 0.12, 0.10])  # e.g., 90th percentile
# Optional: a representative curve (median)
or_env_med = np.array([0.53, 0.40, 0.27, 0.16, 0.11, 0.09])

# C: Saturation mapping robustness (linear vs log vs sigmoid)
geom_err_linear  = np.array([0.52, 0.39, 0.26, 0.16, 0.10, 0.08])
geom_err_log1p   = np.array([0.51, 0.38, 0.25, 0.16, 0.10, 0.08])
geom_err_sigmoid = np.array([0.53, 0.40, 0.27, 0.17, 0.11, 0.09])

# D: Classification AUC vs panel size (same strategies as A)
auc_geometry_aware = np.array([0.62, 0.70, 0.78, 0.84, 0.88, 0.90])
auc_random         = np.array([0.57, 0.62, 0.68, 0.73, 0.77, 0.79])
auc_degree         = np.array([0.59, 0.65, 0.71, 0.77, 0.81, 0.83])
auc_permuted       = np.array([0.55, 0.59, 0.64, 0.68, 0.72, 0.74])

# Optional: confidence intervals for AUC (replace with your CI)
auc_ci_geometry_aware = 0.02 * np.ones_like(auc_geometry_aware)
auc_ci_random         = 0.02 * np.ones_like(auc_random)
auc_ci_degree         = 0.02 * np.ones_like(auc_degree)
auc_ci_permuted       = 0.02 * np.ones_like(auc_permuted)

# E: Permutation null distribution for a single panel size (pick a k of interest)
# Example: distortion at k=40 across N permutations
np.random.seed(0)
perm_null_dist = np.random.normal(loc=0.38, scale=0.03, size=200)  # REPLACE with your permuted runs
true_dist_at_k = 0.16  # REPLACE with your true proteomics geometry-aware distortion at the same k
k_for_null = 40


# -----------------------------
# Figure construction
# -----------------------------
def make_figure5(savepath: str = "Fig5_permutation_robustness.pdf") -> None:
    # Layout: 2 rows x 3 cols, last cell empty
    fig, axes = plt.subplots(2, 3, figsize=(12, 7))
    axA, axB, axC = axes[0]
    axD, axE, axEmpty = axes[1]

    # --- Panel A: Geometry preservation curves
    x = panel_sizes
    axA.plot(x, geom_err_geometry_aware, label="Geometry-aware", marker="o")
    axA.plot(x, geom_err_random,         label="Random",        marker="o", linestyle="--")
    axA.plot(x, geom_err_degree,         label="Degree-based",  marker="o", linestyle="-.")
    axA.plot(x, geom_err_permuted,       label="Proteomics permuted", marker="o", linestyle=":")
    axA.set_title("A  Geometry distortion vs panel size")
    axA.set_xlabel("Panel size |Ω|")
    axA.set_ylabel("Geometry distortion (lower is better)")
    axA.grid(True, alpha=0.3)
    axA.legend(frameon=False)

    # --- Panel B: OR-aggregator sensitivity envelope
    _shade_envelope(axB, x, or_env_lo, or_env_hi, label="OR semantics sensitivity band", alpha=0.25)
    axB.plot(x, or_env_med, label="Representative (median)", marker="o")
    axB.set_title("B  OR-aggregator sensitivity envelope")
    axB.set_xlabel("Panel size |Ω|")
    axB.set_ylabel("Geometry distortion")
    axB.grid(True, alpha=0.3)
    axB.legend(frameon=False)

    # --- Panel C: Saturation mapping robustness
    axC.plot(x, geom_err_linear,  label="Linear f(x)=x", marker="o")
    axC.plot(x, geom_err_log1p,   label="Log f(x)=log(1+x)", marker="o", linestyle="--")
    axC.plot(x, geom_err_sigmoid, label="Sigmoid f(x)=σ(β(x-τ))", marker="o", linestyle="-.")
    axC.set_title("C  Robustness across saturation mappings")
    axC.set_xlabel("Panel size |Ω|")
    axC.set_ylabel("Geometry distortion")
    axC.grid(True, alpha=0.3)
    axC.legend(frameon=False)

    # --- Panel D: AUC vs panel size (with optional CI)
    axD.errorbar(x, auc_geometry_aware, yerr=auc_ci_geometry_aware, label="Geometry-aware", marker="o", capsize=3)
    axD.errorbar(x, auc_random,         yerr=auc_ci_random,         label="Random",        marker="o", linestyle="--", capsize=3)
    axD.errorbar(x, auc_degree,         yerr=auc_ci_degree,         label="Degree-based",  marker="o", linestyle="-.", capsize=3)
    axD.errorbar(x, auc_permuted,       yerr=auc_ci_permuted,       label="Proteomics permuted", marker="o", linestyle=":", capsize=3)
    axD.set_title("D  Cohort classification AUC vs panel size")
    axD.set_xlabel("Panel size |Ω|")
    axD.set_ylabel("AUC (higher is better)")
    axD.set_ylim(0.5, 1.0)
    axD.grid(True, alpha=0.3)
    axD.legend(frameon=False)

    # --- Panel E: Permutation null distribution
    axE.hist(perm_null_dist, bins=20, density=True, alpha=0.6, label="Permutation null")
    axE.axvline(true_dist_at_k, linewidth=2, linestyle="-", label=f"True (k={k_for_null})")
    axE.set_title("E  Permutation null at fixed panel size")
    axE.set_xlabel("Geometry distortion")
    axE.set_ylabel("Density")
    axE.grid(True, alpha=0.3)
    axE.legend(frameon=False)

    # --- Empty panel (turn off)
    axEmpty.axis("off")

    # Global title (optional)
    fig.suptitle("Figure 5. Robustness and permutation controls", y=0.98)

    # Tight layout while keeping suptitle
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    # Save
    fig.savefig(savepath, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {savepath}")


if __name__ == "__main__":
    make_figure5("Fig5_permutation_robustness.pdf")

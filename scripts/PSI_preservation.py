"""
Figure 6 / Fig 5F: PSI preservation
- Panel A: ΔPSI vs panel size (curves)
- Panel B: Heatmap of ΔPSI(P) across pathways/subsystems and methods

Replace the placeholder arrays with your computed outputs.
No seaborn. No hardcoded colors.
"""

from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt


def _as_1d(a, name: str) -> np.ndarray:
    a = np.asarray(a).reshape(-1)
    if a.size == 0:
        raise ValueError(f"{name} is empty.")
    return a


def plot_psi_figure(
    savepath: str = "Fig6_PSI_preservation.pdf",
) -> None:
    # -----------------------------
    # REPLACE THESE WITH REAL OUTPUTS
    # -----------------------------
    panel_sizes = np.array([5, 10, 20, 40, 80, 120])

    # Panel A: global ΔPSI vs panel size
    # ΔPSI = (PSI_Omega - PSI_full) / (PSI_full + eps)
    dpsi_geometry_aware = np.array([-0.22, -0.15, -0.08, -0.03, -0.01, -0.005])
    dpsi_degree         = np.array([-0.30, -0.22, -0.14, -0.08, -0.04, -0.02])
    dpsi_random         = np.array([-0.35, -0.28, -0.20, -0.12, -0.07, -0.04])
    dpsi_permuted       = np.array([-0.40, -0.33, -0.27, -0.20, -0.14, -0.10])

    # Optional: error bars (e.g., across conditions or bootstrap resamples)
    dpsi_err_geometry_aware = 0.01 * np.ones_like(dpsi_geometry_aware)
    dpsi_err_degree         = 0.01 * np.ones_like(dpsi_degree)
    dpsi_err_random         = 0.01 * np.ones_like(dpsi_random)
    dpsi_err_permuted       = 0.01 * np.ones_like(dpsi_permuted)

    # Panel B: subsystem-level ΔPSI(P) heatmap for a fixed panel size
    # Choose one panel size to display (e.g., 40)
    k_fixed = 40

    pathways = [
        "Glycolysis",
        "TCA cycle",
        "Pentose phosphate",
        "Urea cycle",
        "Fatty acid ox.",
        "Cholesterol",
        "Amino acid metab.",
        "Nucleotide metab.",
        "Glutathione",
        "OxPhos / ETC",
        "Bile acid",
        "One-carbon",
    ]

    methods = ["Geometry-aware", "Degree-based", "Random", "Permuted"]

    # Shape: (n_pathways, n_methods)
    # Values are ΔPSI(P) at fixed k. (Less negative is better preservation.)
    heat = np.array([
        [-0.05, -0.10, -0.15, -0.22],
        [-0.03, -0.09, -0.14, -0.25],
        [-0.06, -0.12, -0.18, -0.28],
        [-0.04, -0.11, -0.17, -0.24],
        [-0.07, -0.13, -0.20, -0.30],
        [-0.08, -0.14, -0.22, -0.32],
        [-0.05, -0.12, -0.19, -0.27],
        [-0.06, -0.13, -0.21, -0.29],
        [-0.04, -0.10, -0.16, -0.23],
        [-0.03, -0.08, -0.13, -0.24],
        [-0.07, -0.15, -0.23, -0.33],
        [-0.05, -0.11, -0.18, -0.26],
    ], dtype=float)

    # -----------------------------
    # Figure layout
    # -----------------------------
    fig = plt.figure(figsize=(12, 5))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.05, 1.15])

    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])

    # --- Panel A: ΔPSI curves (errorbars optional)
    x = panel_sizes
    axA.errorbar(x, dpsi_geometry_aware, yerr=dpsi_err_geometry_aware, marker="o", capsize=3, label="Geometry-aware")
    axA.errorbar(x, dpsi_degree,         yerr=dpsi_err_degree,         marker="o", capsize=3, linestyle="--", label="Degree-based")
    axA.errorbar(x, dpsi_random,         yerr=dpsi_err_random,         marker="o", capsize=3, linestyle="-.", label="Random")
    axA.errorbar(x, dpsi_permuted,       yerr=dpsi_err_permuted,       marker="o", capsize=3, linestyle=":", label="Permuted")

    axA.axhline(0.0, linewidth=1, linestyle="-")
    axA.set_title("A  Pathway separability preservation")
    axA.set_xlabel("Panel size |Ω|")
    axA.set_ylabel(r"Normalized change $\Delta \mathrm{PSI}$")
    axA.grid(True, alpha=0.3)
    axA.legend(frameon=False)

    # --- Panel B: heatmap of ΔPSI(P) at fixed panel size
    im = axB.imshow(heat, aspect="auto")  # default colormap

    axB.set_title(f"B  Subsystem-level $\Delta$PSI at |Ω|={k_fixed}")
    axB.set_xticks(np.arange(len(methods)))
    axB.set_xticklabels(methods, rotation=30, ha="right")
    axB.set_yticks(np.arange(len(pathways)))
    axB.set_yticklabels(pathways)

    # annotate heatmap values (optional but helpful)
    for i in range(heat.shape[0]):
        for j in range(heat.shape[1]):
            axB.text(j, i, f"{heat[i, j]:.2f}", ha="center", va="center", fontsize=8)

    # colorbar
    cbar = fig.colorbar(im, ax=axB, fraction=0.046, pad=0.04)
    cbar.set_label(r"$\Delta \mathrm{PSI}(P)$")

    fig.suptitle("PSI-based pathway modularity preservation under partial observation", y=1.02)
    fig.tight_layout()

    fig.savefig(savepath, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {savepath}")


if __name__ == "__main__":
    plot_psi_figure("Fig6_PSI_preservation.pdf")

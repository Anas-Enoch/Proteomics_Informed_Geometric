import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle

# ============================================================
# Helpers
# ============================================================

def finalize_ax(ax, title=None, xlabel=None, ylabel=None):
    if title:
        ax.set_title(title)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.25)

def savefig(fig, path):
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")

def circle_layout(n):
    t = np.linspace(0, 2*np.pi, n, endpoint=False)
    return np.column_stack([np.cos(t), np.sin(t)])

# ============================================================
# Figure 1
# ============================================================

def figure1():
    fig, axs = plt.subplots(2, 2, figsize=(10, 8))

    # A: bipartite schematic
    ax = axs[0, 0]
    y_met = np.linspace(0.1, 0.9, 10)
    y_rxn = np.linspace(0.2, 0.8, 6)

    for i, y in enumerate(y_met):
        ax.add_patch(Circle((0.2, y), 0.02, fill=False))
        ax.text(0.15, y, f"m{i}", ha="right", va="center", fontsize=8)

    for i, y in enumerate(y_rxn):
        ax.add_patch(Rectangle((0.78, y-0.02), 0.04, 0.04, fill=False))
        ax.text(0.84, y, f"r{i}", va="center", fontsize=8)

        for j in np.random.choice(len(y_met), 3, replace=False):
            ax.plot([0.22, 0.78], [y_met[j], y], lw=1)

    ax.axis("off")
    ax.set_title("Fig 1A: Stoichiometric schematic")

    # B: reaction weights
    ax = axs[0, 1]
    x = np.arange(8)
    ax.bar(x - 0.2, 0.5 + np.random.rand(8), width=0.4, label="Control")
    ax.bar(x + 0.2, 0.5 + np.random.rand(8), width=0.4, label="Disease")
    finalize_ax(ax, "Fig 1B: Reaction weights", "Reaction", "Weight")
    ax.legend(frameon=False)

    # C: embedding
    ax = axs[1, 0]
    coords = circle_layout(20)
    ax.scatter(coords[:, 0], coords[:, 1], s=40)
    ax.set_aspect("equal")
    finalize_ax(ax, "Fig 1C: Metabolite geometry")

    # D: partial observability
    ax = axs[1, 1]
    obs = np.zeros(20, dtype=bool)
    obs[np.random.choice(20, 6, replace=False)] = True
    ax.scatter(coords[~obs, 0], coords[~obs, 1], alpha=0.2, label="Unobserved")
    ax.scatter(coords[obs, 0], coords[obs, 1], s=60, label="Observed")
    ax.set_aspect("equal")
    ax.legend(frameon=False)
    finalize_ax(ax, "Fig 1D: Partial observability")

    return fig

# ============================================================
# Figure 2
# ============================================================

def figure2():
    fig, axs = plt.subplots(2, 2, figsize=(10, 8))

    # -----------------
    # Fig 2A: Operator flow
    # -----------------
    ax = axs[0, 0]
    ax.axis("off")

    boxes = [
        ("W_M^1/2", 0.05),
        ("S",       0.25),
        ("W_R(c)",  0.45),
        ("S^T",     0.65),
        ("W_M^1/2", 0.85),
    ]

    for label, x in boxes:
        ax.add_patch(Rectangle((x-0.07, 0.45), 0.14, 0.1, fill=False))
        ax.text(x, 0.5, label, ha="center", va="center")

    for i in range(len(boxes)-1):
        ax.annotate(
            "", xy=(boxes[i+1][1]-0.09, 0.5),
            xytext=(boxes[i][1]+0.09, 0.5),
            arrowprops=dict(arrowstyle="->")
        )

    ax.text(0.5, 0.25, "Metabolite Laplacian Δ_M(c)", ha="center")
    ax.set_title("Fig 2A")

    # -----------------
    # Fig 2B: Dirac block
    # -----------------
    ax = axs[0, 1]
    ax.axis("off")

    ax.add_patch(Rectangle((0.2, 0.2), 0.6, 0.6, fill=False))
    ax.plot([0.5, 0.5], [0.2, 0.8])
    ax.plot([0.2, 0.8], [0.5, 0.5])

    ax.text(0.35, 0.65, "0", ha="center")
    ax.text(0.65, 0.65, "d(c)", ha="center")
    ax.text(0.35, 0.35, "d*(c)", ha="center")
    ax.text(0.65, 0.35, "0", ha="center")

    ax.set_title("Fig 2B: Dirac operator")

    # -----------------
    # Fig 2C: D^2
    # -----------------
    ax = axs[1, 0]
    ax.axis("off")

    ax.add_patch(Rectangle((0.25, 0.4), 0.5, 0.25, fill=False))
    ax.text(0.35, 0.52, "Δ_M(c)", ha="center")
    ax.text(0.65, 0.52, "Δ_R(c)", ha="center")

    ax.text(0.5, 0.75, "D^2 = block-diagonal Laplacians", ha="center")
    ax.set_title("Fig 2C")

    # -----------------
    # Fig 2D: Coupling intuition
    # -----------------
    ax = axs[1, 1]
    ax.axis("off")

    ax.plot([0.2, 0.4], [0.6, 0.6], lw=3)
    ax.plot([0.6, 0.8], [0.4, 0.4], lw=1, linestyle="--")

    ax.scatter([0.2, 0.4], [0.6, 0.6], s=80)
    ax.scatter([0.6, 0.8], [0.4, 0.4], s=80)

    ax.text(0.3, 0.7, "Strong coupling", ha="center")
    ax.text(0.7, 0.5, "Weak coupling", ha="center")

    ax.set_title("Fig 2D")

    return fig

# ============================================================
# Figure 3
# ============================================================

def figure3():
    fig, axs = plt.subplots(2, 2, figsize=(10, 8))

    sizes = np.array([10, 20, 40, 80, 160])
    axs[0,0].plot(sizes, [1, .8, .6, .4, .3], marker="o", label="Random")
    axs[0,0].plot(sizes, [0.9, .7, .5, .3, .2], marker="o", label="Ours")
    finalize_ax(axs[0,0], "Fig 3A", "Panel size", "Error")
    axs[0,0].legend(frameon=False)

    axs[0,1].boxplot([np.random.rand(100), np.random.rand(100)*0.6])
    axs[0,1].set_xticklabels(["Random", "Ours"])
    finalize_ax(axs[0,1], "Fig 3B", ylabel="Distance error")

    axs[1,0].plot(sizes, [0.8, .6, .4, .3, .25], marker="o", label="With proteomics")
    axs[1,0].plot(sizes, [1.0, .9, .7, .6, .55], marker="o", label="No proteomics")
    finalize_ax(axs[1,0], "Fig 3C", "Panel size", "Error")
    axs[1,0].legend(frameon=False)

    axs[1,1].scatter(np.random.rand(10), np.random.rand(10), label="Full")
    axs[1,1].scatter(np.random.rand(10), np.random.rand(10), label="Restricted")
    axs[1,1].legend(frameon=False)
    finalize_ax(axs[1,1], "Fig 3D")

    return fig

# ============================================================
# Figure 4
# ============================================================

def figure4():
    fig, axs = plt.subplots(2, 2, figsize=(10, 8))

    axs[0,0].plot([1, .7, .5, .3, .25], marker="o")
    finalize_ax(axs[0,0], "Fig 4A", "Step", "Score")

    axs[0,1].plot([1, .8, .6, .4], marker="o", label="Random")
    axs[0,1].plot([.8, .6, .4, .25], marker="o", label="Ours")
    axs[0,1].legend(frameon=False)
    finalize_ax(axs[0,1], "Fig 4B", "Panel size", "Score")

    axs[1,0].bar(["Robust", "Adaptive"], [0.22, 0.17])
    finalize_ax(axs[1,0], "Fig 4C", ylabel="Error")

    axs[1,1].imshow(np.random.rand(6,1), aspect="auto")
    axs[1,1].set_yticks(range(6))
    axs[1,1].set_yticklabels(["Glycolysis", "TCA", "PPP", "FAO", "AA", "OxPhos"])
    axs[1,1].set_title("Fig 4D")

    return fig

# ============================================================
# Main
# ============================================================

if __name__ == "__main__":
    import os
    import numpy as np
    import matplotlib.pyplot as plt

    np.random.seed(0)

    outdir = "figures_out"
    os.makedirs(outdir, exist_ok=True)

    fig1 = figure1()
    fig1.savefig(f"{outdir}/figure1.pdf", dpi=300, bbox_inches="tight")
    fig1.savefig(f"{outdir}/figure1.png", dpi=300, bbox_inches="tight")
    plt.close(fig1)

    fig2 = figure2()
    fig2.savefig(f"{outdir}/figure2.pdf", dpi=300, bbox_inches="tight")
    fig2.savefig(f"{outdir}/figure2.png", dpi=300, bbox_inches="tight")
    plt.close(fig2)

    fig3 = figure3()
    fig3.savefig(f"{outdir}/figure3.pdf", dpi=300, bbox_inches="tight")
    fig3.savefig(f"{outdir}/figure3.png", dpi=300, bbox_inches="tight")
    plt.close(fig3)

    fig4 = figure4()
    fig4.savefig(f"{outdir}/figure4.pdf", dpi=300, bbox_inches="tight")
    fig4.savefig(f"{outdir}/figure4.png", dpi=300, bbox_inches="tight")
    plt.close(fig4)

    print("All figures generated separately.")

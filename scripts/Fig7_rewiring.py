import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import scipy.sparse as sp

import cobra
from cobra.util.array import create_stoichiometric_matrix


def load_model(model_path: str):
    p = Path(model_path)
    if not p.exists():
        raise FileNotFoundError(f"Model file not found: {p.resolve()}")

    suf = p.suffix.lower()
    if suf in [".xml", ".sbml"]:
        return cobra.io.read_sbml_model(str(p))
    if suf in [".json"]:
        return cobra.io.load_json_model(str(p))
    if suf in [".mat"]:
        # requires scipy (you installed it). But you must actually have the .mat file.
        return cobra.io.load_matlab_model(str(p))
    raise ValueError("Unsupported model format. Use .xml/.sbml, .json, or .mat")


def read_reaction_weights_csv(csv_path: str, reaction_ids: list[str]) -> np.ndarray:
    """
    CSV format: reaction_id,rho
    """
    df = pd.read_csv(csv_path)
    if "reaction_id" not in df.columns or "rho" not in df.columns:
        raise ValueError("CSV must contain columns: reaction_id, rho")

    m = {rid: float(rho) for rid, rho in zip(df["reaction_id"], df["rho"])}
    rho = np.ones(len(reaction_ids), dtype=float)
    missing = 0
    for j, rid in enumerate(reaction_ids):
        if rid in m:
            rho[j] = m[rid]
        else:
            missing += 1
    if missing > 0:
        print(f"[warn] {missing}/{len(reaction_ids)} reactions missing in {csv_path}; using rho=1 for them.")
    return rho


def get_reaction_subsystem(rxn) -> str | None:
    # COBRApy: many models expose rxn.subsystem
    sub = getattr(rxn, "subsystem", None)
    if sub:
        return str(sub)

    # Some models store subsystem in annotations/notes
    ann = getattr(rxn, "annotation", {}) or {}
    for key in ["subsystem", "subSystem", "subsystems", "Subsystem"]:
        if key in ann and ann[key]:
            val = ann[key]
            if isinstance(val, list):
                return str(val[0])
            return str(val)

    return None


def build_subsystem_maps(model):
    """
    Returns:
      met_index_map: met_id -> index
      sub_to_met_indices: subsystem -> list[int]
    """
    mets = list(model.metabolites)
    met_index_map = {m.id: i for i, m in enumerate(mets)}

    sub_to_mets = {}
    for rxn in model.reactions:
        sub = get_reaction_subsystem(rxn)
        if not sub:
            continue
        # add all metabolites participating in this reaction
        for met in rxn.metabolites.keys():
            sub_to_mets.setdefault(sub, set()).add(met.id)

    # Convert to index lists, filter tiny subsystems
    sub_to_met_indices = {}
    for sub, met_ids in sub_to_mets.items():
        idx = [met_index_map[mid] for mid in met_ids if mid in met_index_map]
        if len(idx) >= 10:  # avoid noisy tiny subsystems
            sub_to_met_indices[sub] = idx

    return met_index_map, sub_to_met_indices


def compute_deformation_scores(S, rho_c1, rho_c2, W_M_diag):
    """
    Compute D(i) = || row_i( Delta_diff ) ||_2
    Delta_diff = W_M^{1/2} S diag(rho_c1 - rho_c2) S^T W_M^{1/2}
    done in sparse form without densifying.
    """
    d_rho = rho_c1 - rho_c2  # (m,)
    WMs = sp.diags(np.sqrt(W_M_diag), format="csr")
    dWR = sp.diags(d_rho, format="csr")

    Delta_diff = WMs @ S @ dWR @ S.T @ WMs  # sparse (n x n)
    # row-wise L2 norm: sqrt(sum_j Delta_ij^2)
    D = np.sqrt(Delta_diff.multiply(Delta_diff).sum(axis=1)).A1
    return D


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True, help="Path to Human1 model (.xml/.sbml/.json/.mat)")
    ap.add_argument("--rho_c1", default=None, help="CSV with reaction_id,rho for condition 1")
    ap.add_argument("--rho_c2", default=None, help="CSV with reaction_id,rho for condition 2")
    ap.add_argument("--out", default="Fig7_rewiring.pdf", help="Output PDF filename")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    np.random.seed(args.seed)

    model = load_model(args.model)
    reaction_ids = [r.id for r in model.reactions]

    # Stoichiometric matrix (sparse CSR)
    S = create_stoichiometric_matrix(model, array_type="dok").tocsr()

    # W_M: if you don’t have metabolite variances yet, use identity
    n = S.shape[0]
    W_M_diag = np.ones(n, dtype=float)

    # Reaction weights: either from CSV, or default to ones (so it runs now)
    if args.rho_c1:
        rho_c1 = read_reaction_weights_csv(args.rho_c1, reaction_ids)
    else:
        rho_c1 = np.ones(len(reaction_ids), dtype=float)

    if args.rho_c2:
        rho_c2 = read_reaction_weights_csv(args.rho_c2, reaction_ids)
    else:
        # if not provided, create a slightly perturbed condition to test pipeline
        rho_c2 = rho_c1.copy()
        rho_c2 *= (1.0 + 0.05 * np.random.randn(len(rho_c2)))

    # Subsystem maps
    met_index_map, sub_to_met_indices = build_subsystem_maps(model)

    # Real deformation
    D = compute_deformation_scores(S, rho_c1, rho_c2, W_M_diag)

    # Pathway rewiring index R(P)
    R = {sub: float(np.mean(D[idx])) for sub, idx in sub_to_met_indices.items()}
    R_sorted = sorted(R.items(), key=lambda x: x[1], reverse=True)

    # Permutation null: permute reaction weights across reactions (distribution matched)
    rho_c1_perm = rho_c1.copy()
    np.random.shuffle(rho_c1_perm)
    D_perm = compute_deformation_scores(S, rho_c1_perm, rho_c2, W_M_diag)

    R_perm = {sub: float(np.mean(D_perm[idx])) for sub, idx in sub_to_met_indices.items()}

    # ---------------- Plot ----------------
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # A: distribution of D(i)
    sns.histplot(D, bins=60, ax=axes[0], kde=True)
    axes[0].set_title("A  Metabolite-level deformation $D(i)$")
    axes[0].set_xlabel("$D(i)$")
    axes[0].set_ylabel("Count")

    # B: top 12 pathways
    top = R_sorted[:12]
    names = [x[0] for x in top]
    vals = [x[1] for x in top]
    sns.barplot(x=vals, y=names, ax=axes[1])
    axes[1].set_title("B  Top rewired subsystems (real)")
    axes[1].set_xlabel("$R(P)$")

    # C: real vs permuted for those top subsystems
    real_vals = vals
    perm_vals = [R_perm.get(k, 0.0) for k in names]
    y = np.arange(len(names))
    h = 0.35
    axes[2].barh(y - h/2, real_vals, height=h, label="Real")
    axes[2].barh(y + h/2, perm_vals, height=h, label="Permuted")
    axes[2].set_yticks(y)
    axes[2].set_yticklabels(names)
    axes[2].invert_yaxis()
    axes[2].set_title("C  Distribution-matched permutation control")
    axes[2].set_xlabel("$R(P)$")
    axes[2].legend()

    plt.tight_layout()
    plt.savefig(args.out)
    print(f"Saved: {args.out}")
    plt.show()


if __name__ == "__main__":
    main()

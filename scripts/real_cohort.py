#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import cobra
from scipy import sparse
from scipy.sparse.linalg import eigsh
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, accuracy_score
from proteomics_weighting import build_operator

# -------------------------
# Text parsing (Workbench)
# -------------------------

def parse_workbench_txt(path: str):
    """
    Robust parser for Metabolomics Workbench TXT exports containing:
    #NMR_METABOLITE_DATA
    ...
    NMR_METABOLITE_DATA_START
    Samples <tab> sample_ids...
    Factors <tab> Factor strings...
    metabolite rows...

    Returns:
      df:        DataFrame [metabolites x samples] floats
      groups:    dict sample_id -> group label (uses 'Group:...' if present)
      met_ann:   DataFrame indexed by metabolite raw names, with ID fields if present
    """
    path = Path(path)
    lines = path.read_text(errors="ignore").splitlines()

    # locate NMR section
    start = None
    for i, ln in enumerate(lines):
        if ln.strip().startswith("#NMR_METABOLITE_DATA"):
            start = i
            break
    if start is None:
        raise ValueError("Could not find #NMR_METABOLITE_DATA section.")

    # locate 'Samples' header within the next ~400 lines
    header_idx = None
    for j in range(start + 1, min(start + 500, len(lines))):
        s = lines[j].strip()
        if not s:
            continue
        if s.startswith("#"):
            continue
        if s.lower().startswith("samples\t") or s.lower() == "samples":
            header_idx = j
            break
    if header_idx is None:
        raise ValueError("Could not locate 'Samples' row after #NMR_METABOLITE_DATA.")

    # next meaningful row should be Factors / Factor
    def next_meaningful(idx):
        for k in range(idx + 1, len(lines)):
            t = lines[k].strip()
            if not t:
                continue
            if t.startswith("#"):
                continue
            return k, lines[k]
        return None, None

    samples_line = lines[header_idx].strip()
    factors_idx, factors_line = next_meaningful(header_idx)
    if factors_line is None:
        raise ValueError("Found 'Samples' row but no subsequent 'Factors' row.")

    if not (factors_line.strip().lower().startswith("factors") or factors_line.strip().lower().startswith("factor")):
        raise ValueError(f"Expected 'Factors' row after 'Samples', got: {factors_line[:120]}")

    samples = samples_line.split("\t")[1:]
    factors = factors_line.split("\t")[1:]

    if len(samples) == 0:
        raise ValueError("Parsed zero samples from 'Samples' row.")
    if len(factors) != len(samples):
        # keep only aligned portion
        m = min(len(samples), len(factors))
        samples = samples[:m]
        factors = factors[:m]

    # parse group labels from factors
    groups = {}
    for sid, fac in zip(samples, factors):
        m = re.search(r"Group:([^|]+)", fac)
        grp = m.group(1).strip() if m else fac.strip()
        groups[sid] = grp

    # metabolite table starts after factors row until end marker
    mets = []
    values = []
    raw = []

    for ln in lines[factors_idx + 1:]:
        s = ln.strip()
        if not s:
            continue
        if s.startswith("#"):
            continue
        if s.strip().startswith("NMR_METABOLITE_DATA_END"):
            break

        parts = ln.split("\t")
        if len(parts) < 2:
            continue
        met = parts[0].strip()
        nums = parts[1:]

        if len(nums) != len(samples):
            continue

        try:
            row = [float(x) if x != "" else np.nan for x in nums]
        except ValueError:
            continue

        mets.append(met)
        values.append(row)
        raw.append(met)

    df = pd.DataFrame(values, index=mets, columns=samples, dtype=float)

    # crude ID extraction if present in raw label
    def extract_id(s: str, key: str):
        m = re.search(rf"{re.escape(key)}\s*:\s*([A-Za-z0-9_\-]+)", s)
        return m.group(1) if m else None

    ann = []
    for met in raw:
        ann.append(
            dict(
                raw=met,
                hmdb=extract_id(met, "HMDB ID"),
                kegg=extract_id(met, "KEGG ID"),
                chebi=extract_id(met, "ChEBI ID"),
                pubchem=extract_id(met, "PubChem ID"),
            )
        )
    met_ann = pd.DataFrame(ann, index=mets)

    return df, groups, met_ann


# -------------------------
# Name normalization + mapping
# -------------------------

def normalize_met_name(s: str) -> str:
    if s is None:
        return ""
    s = s.strip().lower()
    # remove bracket punctuation, keep alphanum/spaces
    s = re.sub(r"[\(\)\[\]\{\}]", " ", s)
    s = re.sub(r"[^a-z0-9\s\-]", " ", s)
    s = s.replace("-", " ")
    s = re.sub(r"\s+", " ", s).strip()

    # normalize common chemical naming endings
    # (IMPORTANT: do NOT blindly drop 'acid' because model might use 'citric acid' vs 'citrate'
    # we handle that via aliases)
    return s


ALIASES = {
    # acids / salts / common swaps
    "citric acid": "citrate",
    "lactic acid": "lactate",
    "vitamin c": "ascorbate",
    "3 hydroxybutyric acid": "beta hydroxybutyrate",
    "3 hydroxybutyrate": "beta hydroxybutyrate",
    "ketoisovaleric acid": "3 methyl 2 oxobutanoate",
    # your dataset odd label
    "car 2 0": "acetate",
    "acetic acid": "acetate",
}

def build_name_to_ids(model: cobra.Model):
    """
    Build normalized name -> candidate metabolite IDs.
    Human-GEM IDs end with compartment letters (c,e,m,n,...) e.g. MAM01621e
    """
    name_to_ids = {}
    for met in model.metabolites:
        if met.name:
            k = normalize_met_name(met.name)
            name_to_ids.setdefault(k, []).append(met.id)
    # sort candidates: prefer cytosol 'c' then extracellular 'e'
    def rank(mid: str):
        if mid.endswith("c"):
            return 0
        if mid.endswith("e"):
            return 1
        return 9
    name_to_ids = {k: sorted(set(v), key=rank) for k, v in name_to_ids.items()}
    return name_to_ids


def map_workbench_to_model(met_name: str, met_ann_row: pd.Series, name_to_ids: dict):
    """
    For ST003506: IDs are missing, so mapping is primarily name-based with aliases + light fuzzy fallback.
    Returns chosen model metabolite id or None.
    """
    raw = met_name
    key = normalize_met_name(raw)

    # apply alias (on normalized keys)
    key2 = ALIASES.get(key, key)

    # exact name match
    cand = name_to_ids.get(key2, [])
    if cand:
        return cand[0]

    # very light fuzzy: try contains match ONLY if unique
    # (prevents garbage mappings)
    hits = []
    for k in name_to_ids.keys():
        if key2 and (key2 in k or k in key2):
            hits.append(k)
    hits = sorted(set(hits), key=len)
    if len(hits) == 1:
        return name_to_ids[hits[0]][0]

    return None


# -------------------------
# Operator construction
# -------------------------
# build_operator() from proteomics_weighting.py is used directly.
# It supports three modes: baseline (W_R=I), proteomics (CPTAC-informed),
# and permuted (distribution-matched null). See proteomics_weighting.py.

def submatrix_sparse(A: sparse.spmatrix, idx: np.ndarray):
    """
    Extract A[idx, idx] as CSR.
    """
    idx = np.asarray(idx, dtype=int)
    return A[idx][:, idx].tocsr()


def spectral_embedding(L: sparse.spmatrix, k: int, seed: int = 0):
    """
    Returns k smallest nonzero eigenvectors of L (as embedding coordinates).
    """
    n = L.shape[0]
    if k >= n:
        k = max(1, n - 1)

    # add tiny ridge to avoid numerical issues
    L = L + 1e-12 * sparse.eye(n, format="csr")

    # compute k+1 smallest eigenpairs, drop the first (near-zero) component
    vals, vecs = eigsh(L, k=min(k + 1, n - 1), which="SM", tol=1e-6, maxiter=5000)
    order = np.argsort(vals)
    vals = vals[order]
    vecs = vecs[:, order]

    # drop first eigenvector if it's the trivial one
    if vecs.shape[1] > 1:
        vecs = vecs[:, 1:k + 1]
        vals = vals[1:k + 1]
    else:
        vecs = vecs[:, :k]
        vals = vals[:k]

    return vals, vecs

def permutation_test_auc(X_raw, y, U, n_perm=1000, n_splits=5, random_state=42):
    rng = np.random.default_rng(random_state)

    def cv_auc(labels):
        skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)
        fold_aucs = []

        for tr, te in skf.split(X_raw, labels):
            imputer = SimpleImputer(strategy="mean")
            X_tr_raw = imputer.fit_transform(X_raw[tr])
            X_te_raw = imputer.transform(X_raw[te])

            X_tr = X_tr_raw @ U
            X_te = X_te_raw @ U

            pipe = Pipeline([
                ("scaler", StandardScaler()),
                ("lr", LogisticRegression(max_iter=5000, solver="lbfgs"))
            ])

            pipe.fit(X_tr, labels[tr])
            prob = pipe.predict_proba(X_te)[:, 1]
            fold_aucs.append(roc_auc_score(labels[te], prob))

        return np.mean(fold_aucs)

    observed_auc = cv_auc(y)

    perm_aucs = np.array([
        cv_auc(rng.permutation(y))
        for _ in range(n_perm)
    ])

    p_value = (np.sum(perm_aucs >= observed_auc) + 1) / (n_perm + 1)
    return observed_auc, perm_aucs, p_value
# -------------------------
# Main experiment
# -------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", required=True, help="Metabolomics Workbench TXT table (e.g., ST003506_AN005756.txt)")
    ap.add_argument("--model", required=True, help="Human-GEM SBML XML path")
    ap.add_argument("--out", default="Fig6_real_cohort_classification.pdf")
    ap.add_argument("--k_eigs", type=int, default=10)
    ap.add_argument("--repeats", type=int, default=50)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--condition", default="tumor", help="Proteomics condition to use")
    ap.add_argument(
        "--operator",
        choices=["baseline", "proteomics", "permuted"],
        default="baseline",
        help=(
            "baseline   = S S^T  (W_R = I, topology only)\n"
            "proteomics = S W_R^(c) S^T  (CPTAC-informed)\n"
            "permuted   = S W_R_perm S^T  (distribution-matched null)"
        ),
    )
    ap.add_argument(
        "--proteomics",
        default=None,
        help="Path to CPTAC TSV file (required when --operator != baseline)",
    )
    ap.add_argument(
        "--hgnc",
        default=None,
        help=(
            "Path to HGNC complete set TSV for symbol->Ensembl remapping. "
            "Required when --operator=proteomics and model gene.name fields are empty. "
            "Download: curl -o data/hgnc_complete_set.txt "
            "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
        ),
    )
    args = ap.parse_args()

    df, sample_groups, met_ann = parse_workbench_txt(args.table)

    # keep only two groups (Control + Lymphedema for this dataset)
    groups_series = pd.Series(sample_groups)
    keep = groups_series[groups_series.isin(["Control", "Lymphedema"])].index.tolist()
    df = df[keep]
    y = groups_series.loc[keep].map({"Control": 0, "Lymphedema": 1}).astype(int).values

    print(f"Samples kept: {len(keep)} (Control={(y==0).sum()}, Lymphedema={(y==1).sum()})")
    print(f"Measured metabolites: {df.shape[0]}")

    # load model
    model = cobra.io.read_sbml_model(args.model)

    # build operator (baseline)
    Delta_full = build_operator(
        model,
        mode=args.operator,
        proteomics_path=args.proteomics,
        condition=args.condition,
        rho_0=0.1,
        alpha=1.0,
        seed=args.seed,
        hgnc_path=args.hgnc,
    )


    # build mapping: Workbench metabolite -> model metabolite id
    name_to_ids = build_name_to_ids(model)
    model_met_index = {m.id: i for i, m in enumerate(model.metabolites)}

    mapped = {}
    unmapped = []
    for met_name in df.index:
        chosen = map_workbench_to_model(met_name, met_ann.loc[met_name], name_to_ids)
        if chosen is None:
            unmapped.append(met_name)
            continue
        if chosen not in model_met_index:
            unmapped.append(met_name)
            continue
        mapped[met_name] = chosen

    print(f"Mapped to model: {len(mapped)}")
    print(f"Unmapped: {len(unmapped)}")
    if len(mapped) < 15:
        # print a few unmapped examples to debug
        print("Example unmapped metabolites:", unmapped[:20])
        raise RuntimeError("Too few metabolites mapped (<15). Mapping failed; expand aliases or relax fuzzy matching.")

    # build restricted Laplacian on mapped metabolites (in model index space)
    wb_mets = list(mapped.keys())
    model_ids = [mapped[m] for m in wb_mets]
    idx = np.array([model_met_index[mid] for mid in model_ids], dtype=int)

    Delta_obs = submatrix_sparse(Delta_full, idx)

    # compute spectral embedding basis (metabolite-space)
    _, U = spectral_embedding(Delta_obs, k=args.k_eigs, seed=args.seed)  # shape: (#mets, k)

    # Build sample features by projecting metabolite concentrations onto eigenvectors
    X_raw = df.loc[wb_mets].T.values  # shape: (#samples, #mets)

    # repeated stratified K-fold cross-validation
    rskf = RepeatedStratifiedKFold(
        n_splits=5,
        n_repeats=10,
        random_state=args.seed
    )

    aucs = []
    accs = []

    for tr, te in rskf.split(X_raw, y):
        # fold-specific imputation
        imputer = SimpleImputer(strategy="mean")
        X_tr_raw = imputer.fit_transform(X_raw[tr])
        X_te_raw = imputer.transform(X_raw[te])

        # project into the fixed operator-derived spectral basis
        X_tr = X_tr_raw @ U
        X_te = X_te_raw @ U

        pipe = Pipeline([
            ("scaler", StandardScaler()),
            ("lr", LogisticRegression(max_iter=5000, solver="lbfgs"))
        ])

        pipe.fit(X_tr, y[tr])
        prob = pipe.predict_proba(X_te)[:, 1]
        pred = (prob >= 0.5).astype(int)

        aucs.append(roc_auc_score(y[te], prob))
        accs.append(accuracy_score(y[te], pred))

    aucs = np.array(aucs)
    accs = np.array(accs)

    print(f"AUROC: mean={aucs.mean():.3f} ± {aucs.std():.3f}")
    print(f"ACC:   mean={accs.mean():.3f} ± {accs.std():.3f}")

    observed_auc, perm_aucs, perm_p = permutation_test_auc(
        X_raw, y, U, n_perm=1000, n_splits=5, random_state=args.seed
    )

    print(f"Permutation test observed AUC: {observed_auc:.3f}")
    print(f"Permutation p-value: {perm_p:.4f}")
    print(f"Permutation null mean ± sd: {perm_aucs.mean():.3f} ± {perm_aucs.std():.3f}")    # plot figure
    plt.figure(figsize=(7.2, 3.6))
    plt.boxplot([aucs, accs], tick_labels=["AUROC", "Accuracy"])
    plt.ylim(0.0, 1.0)
    plt.title("Real-cohort classification using operator spectral embeddings (ST003506)")
    plt.ylabel("Score")
    plt.tight_layout()
    plt.savefig(args.out)
    print(f"Saved: {args.out}")


if __name__ == "__main__":
    main()

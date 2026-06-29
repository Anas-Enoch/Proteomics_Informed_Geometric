#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
proteomics_weighting.py
=======================
Drop-in proteomics weighting module for real_cohort.py.

Implements Sections 4.2-4.5 and 4.7 of the paper:

    baseline   : Delta_M = S S^T              (W_R = I, topology only)
    proteomics : Delta_M = S W_R^(c) S^T      (CPTAC-informed)
    permuted   : Delta_M = S W_R_perm S^T     (distribution-matched null)

NOTE: This module sets W_M = I (no metabolite reliability weighting).
For the cohort classification experiment, metabolite reliability weights
were set to identity. Missing proteomics values do not zero out reactions;
instead they revert toward baseline coupling via fallback=1.0 and
confidence blending (Eq. 8).

Usage
-----
Replace build_metabolite_laplacian(model) in real_cohort.py with:

    from proteomics_weighting import build_operator
    Delta_full = build_operator(model, mode=args.operator,
                                proteomics_path=args.proteomics,
                                seed=args.seed)
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import cobra
from scipy import sparse


# ==============================================================================
# Public API
# ==============================================================================

def build_operator(
    model: cobra.Model,
    mode: str = "baseline",
    proteomics_path: str | None = None,
    condition: str = "tumor",
    rho_0: float = 0.1,
    alpha: float = 1.0,
    seed: int = 42,
    hgnc_path: str | None = None,
) -> sparse.csr_matrix:
    """
    Build the metabolite Laplacian Delta_M for the requested operator mode.

    Parameters
    ----------
    model           : COBRApy model (Human-GEM)
    mode            : 'baseline' | 'proteomics' | 'permuted'
    proteomics_path : path to CPTAC TSV (required for proteomics / permuted)
    condition       : which condition to extract from proteomics file
    rho_0           : baseline coupling, keeps all weights > 0 (Eq. 9)
    alpha           : proteomics dynamic range (Eq. 9)
    seed            : RNG seed for permutation control
    hgnc_path       : optional path to HGNC complete set TSV for symbol->Ensembl
                      mapping when model gene.name fields are empty.
                      Download from: https://www.genenames.org/download/statistics-and-files/
                      (file: hgnc_complete_set.txt)

    Returns
    -------
    Delta_M : scipy CSR sparse matrix, shape (n_metabolites, n_metabolites)
    """
    mode = mode.lower().strip()
    if mode not in ("baseline", "proteomics", "permuted"):
        raise ValueError(
            "mode must be 'baseline', 'proteomics', or 'permuted', got '{}'".format(mode)
        )

    S = _get_stoichiometric_matrix(model)
    n_rxn = S.shape[1]

    if mode == "baseline":
        W_R = np.ones(n_rxn)
        print("[Operator] mode=baseline  (W_R = I)")

    else:
        if proteomics_path is None:
            raise ValueError(
                "proteomics_path is required when mode='proteomics' or 'permuted'"
            )
        gene_expr = _load_and_normalize_proteomics(proteomics_path, condition)
        gene_expr = _remap_symbols_to_ensembl(model, gene_expr, hgnc_path=hgnc_path)
        _report_gene_coverage(model, gene_expr)
        W_R = _compute_reaction_weights(model, gene_expr, rho_0, alpha)

        if mode == "permuted":
            rng = np.random.default_rng(seed)
            W_R = rng.permutation(W_R)
            print("[Operator] mode=permuted  (distribution-matched null, seed={})".format(seed))
        else:
            print("[Operator] mode=proteomics  (CPTAC W_R^(c))")

    # Delta_M = S diag(W_R) S^T
    W_diag = sparse.diags(W_R, format="csr")
    Delta = (S @ W_diag @ S.T).tocsr()

    _report_operator(Delta, W_R)
    return Delta


# ==============================================================================
# Proteomics loading  (Sections 4.2.1 - 4.2.2)
# ==============================================================================

def _load_and_normalize_proteomics(filepath: str, condition: str) -> dict:
    """
    Load CPTAC TSV, aggregate to gene level via median (Eq. 3),
    then apply within-condition median normalization (Eq. 4).

    Supported file formats
    ----------------------
    Wide (pre-averaged per condition):
        gene_symbol | tumor_mean | normal_mean | ...

    Long (per-sample rows):
        gene_symbol | sample_id | condition | log2_ratio | ...

    Returns
    -------
    dict {gene_symbol: normalized_abundance}
    """
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError("Proteomics file not found: {}".format(filepath))

    df = pd.read_csv(str(path), sep="\t", low_memory=False)
    df.columns = [c.strip() for c in df.columns]
    col_lower = {c.lower(): c for c in df.columns}

    # find gene column
    gene_col = (
        col_lower.get("gene_symbol")
        or col_lower.get("gene")
        or col_lower.get("genename")
        or col_lower.get("gene.symbol")
    )
    if gene_col is None:
        raise ValueError(
            "Cannot find gene column. Available columns:\n{}".format(list(df.columns))
        )

    # try wide format first
    expr_col = None
    for candidate in [
        condition,
        condition + "_mean",
        condition.capitalize() + "_mean",
        "mean_" + condition,
    ]:
        if candidate.lower() in col_lower:
            expr_col = col_lower[candidate.lower()]
            break

    if expr_col is not None:
        df_work = df[[gene_col, expr_col]].copy()
        df_work.columns = ["gene", "abundance"]

    else:
        # long format: filter rows by condition column
        cond_col = (
            col_lower.get("condition")
            or col_lower.get("group")
            or col_lower.get("type")
        )
        val_col = (
            col_lower.get("log2_ratio")
            or col_lower.get("abundance")
            or col_lower.get("intensity")
            or col_lower.get("expression")
        )
        if cond_col is None or val_col is None:
            raise ValueError(
                "Could not identify expression column for condition='{}'.\n"
                "Available columns: {}\n"
                "Rename columns to include 'tumor_mean'/'normal_mean', "
                "or ensure 'condition' and 'abundance' columns exist.".format(
                    condition, list(df.columns)
                )
            )
        mask = df[cond_col].astype(str).str.lower().str.contains(
            condition.lower(), na=False
        )
        if mask.sum() == 0:
            raise ValueError(
                "No rows matched condition='{}' in column '{}'.".format(
                    condition, cond_col
                )
            )
        df_work = df.loc[mask, [gene_col, val_col]].copy()
        df_work.columns = ["gene", "abundance"]

    # clean
    df_work = df_work.dropna(subset=["gene", "abundance"])
    df_work["abundance"] = pd.to_numeric(df_work["abundance"], errors="coerce")
    df_work = df_work[df_work["abundance"] > 0].copy()

    if len(df_work) == 0:
        raise ValueError(
            "No valid (gene, abundance) pairs remain after cleaning. "
            "Check that abundance values are positive and numeric."
        )

    # aggregate to gene level via median (Eq. 3)
    gene_expr = df_work.groupby("gene")["abundance"].median().to_dict()

    # within-condition median normalization (Eq. 4)
    eps = 1e-9
    median_val = float(np.median(list(gene_expr.values())))
    gene_expr_norm = {
        g: float(v) / (median_val + eps)
        for g, v in gene_expr.items()
    }

    print(
        "[Proteomics] Loaded {} genes  (condition='{}', file={})".format(
            len(gene_expr_norm), condition, path.name
        )
    )
    return gene_expr_norm


# ==============================================================================
# Gene identifier remapping  (symbol -> Ensembl)
# ==============================================================================

def _remap_symbols_to_ensembl(
    model: cobra.Model,
    gene_expr: dict,
    hgnc_path: str | None = None,
) -> dict:
    """
    Remap gene_expr keys from gene symbols (e.g. 'TP53') to Ensembl IDs
    (e.g. 'ENSG00000141510') so they match Human-GEM GPR rule identifiers.

    Strategy (tried in order):
    1. Model gene.name fields  -- works if the XML is annotated
    2. HGNC complete set TSV   -- works when gene.name is empty (your case)
    3. Identity passthrough    -- used if proteomics already has Ensembl IDs

    Obtaining the HGNC file (one-time download, ~10 MB):
        https://www.genenames.org/download/statistics-and-files/
        -> "Complete HGNC dataset" -> Download TXT
        Save as e.g.  data/hgnc_complete_set.txt

    Parameters
    ----------
    model      : COBRApy model
    gene_expr  : dict {gene_symbol_or_id: abundance}
    hgnc_path  : path to HGNC complete set TSV (optional but recommended
                 when model gene.name fields are empty)

    Returns
    -------
    dict {ensembl_id: abundance}
    """
    # ── build symbol -> [ensembl_id, ...] map ──────────────────────────────

    symbol_to_ensembl: dict[str, list[str]] = {}

    # Strategy 1: model gene.name (works when XML is annotated)
    n_named = 0
    for g in model.genes:
        sym = (g.name or "").strip()
        if sym:
            symbol_to_ensembl.setdefault(sym, []).append(g.id)
            n_named += 1
        # always allow Ensembl passthrough
        symbol_to_ensembl.setdefault(g.id, []).append(g.id)

    if n_named == 0:
        print("[Remap] model gene.name fields are empty -- trying alternative sources")

        # Strategy 2: HGNC complete set file
        if hgnc_path is not None:
            hgnc_file = Path(hgnc_path)
            if not hgnc_file.exists():
                raise FileNotFoundError(
                    "HGNC file not found: {}".format(hgnc_path)
                )
            hgnc_df = pd.read_csv(str(hgnc_file), sep="\t", low_memory=False,
                                  usecols=lambda c: c in ("symbol", "ensembl_gene_id"))
            if "symbol" not in hgnc_df.columns or "ensembl_gene_id" not in hgnc_df.columns:
                raise ValueError(
                    "HGNC file must contain columns 'symbol' and 'ensembl_gene_id'. "
                    "Found: {}".format(list(hgnc_df.columns))
                )
            hgnc_df = hgnc_df.dropna(subset=["symbol", "ensembl_gene_id"])

            # restrict to Ensembl IDs that actually appear in the model
            model_ensembl = {g.id for g in model.genes}
            hgnc_df = hgnc_df[hgnc_df["ensembl_gene_id"].isin(model_ensembl)]

            for _, row in hgnc_df.iterrows():
                sym = str(row["symbol"]).strip()
                eid = str(row["ensembl_gene_id"]).strip()
                if sym and eid:
                    symbol_to_ensembl.setdefault(sym, []).append(eid)

            print("[Remap] HGNC table loaded: {} symbol->Ensembl mappings "
                  "restricted to model genes".format(
                      sum(len(v) for v in symbol_to_ensembl.values())
                  ))

        else:
            # Strategy 3: identity passthrough only
            # Useful if proteomics file already uses Ensembl IDs
            print(
                "[Remap] No HGNC file provided. "
                "Attempting identity passthrough (works only if proteomics "
                "keys are already Ensembl IDs). "
                "To fix properly, re-run with --hgnc data/hgnc_complete_set.txt"
            )

    # ── remap ──────────────────────────────────────────────────────────────

    remapped: dict[str, float] = {}
    n_mapped = 0
    n_already_ensembl = 0

    for sym, abund in gene_expr.items():
        targets = symbol_to_ensembl.get(sym, [])
        if not targets:
            continue
        for eid in targets:
            if eid not in remapped or abund > remapped[eid]:
                remapped[eid] = abund
        n_mapped += 1
        if sym.startswith("ENSG"):
            n_already_ensembl += 1

    print(
        "[Remap] {}/{} proteomics genes mapped to model Ensembl IDs "
        "({} already had Ensembl format)".format(
            n_mapped, len(gene_expr), n_already_ensembl
        )
    )

    if len(remapped) == 0:
        raise RuntimeError(
            "Symbol-to-Ensembl remapping produced zero matches.\n\n"
            "Your proteomics keys  : {}\n"
            "Your model gene IDs   : {}\n\n"
            "Fix: download the HGNC complete set file and pass it via --hgnc:\n"
            "  curl -o data/hgnc_complete_set.txt "
            "https://storage.googleapis.com/public-download-files/"
            "hgnc/tsv/tsv/hgnc_complete_set.txt\n"
            "Then re-run with:\n"
            "  --hgnc data/hgnc_complete_set.txt".format(
                list(gene_expr.keys())[:5],
                [g.id for g in list(model.genes)[:5]],
            )
        )

    return remapped


# ==============================================================================
# Gene coverage diagnostic
# ==============================================================================

def _report_gene_coverage(model: cobra.Model, gene_expr: dict) -> None:
    """
    Print what fraction of model GPR genes are covered by proteomics.

    If coverage is very low, most reactions fall back to baseline and
    the proteomics-informed operator will be numerically near-identical
    to the topology-only operator despite the mode being set to 'proteomics'.
    This is almost always a gene identifier mismatch (e.g. symbols vs Ensembl).
    """
    all_model_genes = {
        g.id
        for rxn in model.reactions
        for g in rxn.genes
    }
    matched = all_model_genes & set(gene_expr.keys())
    pct = 100.0 * len(matched) / max(len(all_model_genes), 1)

    print(
        "[GPR] Gene coverage: {}/{} model GPR genes matched to proteomics "
        "({:.1f}%)".format(len(matched), len(all_model_genes), pct)
    )

    if len(matched) == 0:
        raise RuntimeError(
            "Zero model genes matched to proteomics. "
            "This is a gene identifier mismatch. "
            "Check: list(gene_expr.keys())[:5] vs "
            "[g.id for g in list(model.reactions)[0].genes]"
        )

    if pct < 10.0:
        print(
            "[GPR] WARNING: fewer than 10% of model genes matched. "
            "The proteomics operator will be near-identical to baseline. "
            "Check whether your proteomics file uses gene symbols, "
            "Ensembl IDs, or another identifier format."
        )


# ==============================================================================
# GPR evaluation  (Section 4.3)
# ==============================================================================

def _evaluate_gpr(rule_str: str, gene_expr: dict, fallback: float = 1.0) -> float:
    """
    Recursively evaluate a GPR rule string using gene expression values.

        AND  ->  min   (limiting subunit principle, Eq. 5)
        OR   ->  max   (most expressed isoenzyme, Eq. 6)

    Parameters
    ----------
    rule_str  : GPR rule as plain string, e.g. "(G1 and G2) or G3"
    gene_expr : normalized gene abundance dict
    fallback  : value for genes absent from proteomics

    Returns
    -------
    float  reaction activity proxy kappa_r^(c)
    """
    rule_str = rule_str.strip()
    if not rule_str:
        return fallback

    rule_str = _strip_outer_parens(rule_str)

    # split on top-level OR (lower precedence than AND)
    tokens = _split_top_level(rule_str, " or ")
    if len(tokens) > 1:
        return max(
            _evaluate_gpr(t.strip(), gene_expr, fallback)
            for t in tokens
        )

    # split on top-level AND
    tokens = _split_top_level(rule_str, " and ")
    if len(tokens) > 1:
        return min(
            _evaluate_gpr(t.strip(), gene_expr, fallback)
            for t in tokens
        )

    # base case: single gene ID
    gene_id = rule_str.strip("() ")
    return float(gene_expr.get(gene_id, fallback))


def _strip_outer_parens(s: str) -> str:
    """Remove a matching outer parenthesis pair only if it wraps the whole string."""
    s = s.strip()
    if not (s.startswith("(") and s.endswith(")")):
        return s
    depth = 0
    for i, ch in enumerate(s):
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        # depth returns to 0 before the last char: parens do not wrap fully
        if depth == 0 and i < len(s) - 1:
            return s
    return s[1:-1].strip()


def _split_top_level(expr: str, delimiter: str) -> list:
    """Split expr on delimiter only at parenthesis depth 0."""
    parts = []
    current = []
    depth = 0
    i = 0
    dlen = len(delimiter)
    expr_lower = expr.lower()
    delim_lower = delimiter.lower()

    while i < len(expr):
        ch = expr[i]
        if ch == "(":
            depth += 1
            current.append(ch)
            i += 1
        elif ch == ")":
            depth -= 1
            current.append(ch)
            i += 1
        elif depth == 0 and expr_lower[i:i + dlen] == delim_lower:
            parts.append("".join(current).strip())
            current = []
            i += dlen
        else:
            current.append(ch)
            i += 1

    if current:
        parts.append("".join(current).strip())

    return parts if len(parts) > 1 else [expr]


# ==============================================================================
# Reaction weight construction  (Sections 4.4 - 4.5)
# ==============================================================================

def _compute_reaction_weights(
    model: cobra.Model,
    gene_expr: dict,
    rho_0: float,
    alpha: float,
) -> np.ndarray:
    """
    Build reaction weight vector rho_r^(c)  (Eq. 7-9).

    For each reaction:
      1. Compute coverage q_r = |observed genes| / |all GPR genes|  (Eq. 7)
      2. Evaluate GPR rule on observed gene abundances
      3. Blend with prior via q_r  (Eq. 8)
      4. Apply rho_0 + alpha * kappa  (Eq. 9)
    """
    observed = set(gene_expr.keys())
    n_rxn = len(model.reactions)
    W_R = np.empty(n_rxn)
    zero_coverage_count = 0

    for i, rxn in enumerate(model.reactions):
        rule = (rxn.gene_reaction_rule or "").strip()
        genes = rxn.genes

        if not genes or not rule:
            W_R[i] = rho_0 + alpha * 1.0
            continue

        gene_ids = {g.id for g in genes}
        q_r = len(gene_ids & observed) / len(gene_ids)     # Eq. 7

        kappa = _evaluate_gpr(rule, gene_expr, fallback=1.0)
        kappa_blended = q_r * kappa + (1.0 - q_r) * 1.0   # Eq. 8
        W_R[i] = rho_0 + alpha * kappa_blended             # Eq. 9

        if q_r == 0.0:
            zero_coverage_count += 1

    n_covered = n_rxn - zero_coverage_count
    print(
        "[GPR] Reaction coverage: {}/{} reactions have >= 1 matched gene "
        "({:.1f}%)".format(n_covered, n_rxn, 100.0 * n_covered / n_rxn)
    )
    print(
        "[GPR] W_R  min={:.4f}  max={:.4f}  mean={:.4f}  std={:.4f}".format(
            float(W_R.min()), float(W_R.max()),
            float(W_R.mean()), float(W_R.std())
        )
    )

    if W_R.std() < 1e-6:
        print(
            "[GPR] WARNING: W_R has near-zero variance. "
            "The proteomics operator is numerically identical to baseline. "
            "Verify gene identifier matching."
        )

    return W_R


# ==============================================================================
# Stoichiometric matrix helper
# ==============================================================================

def _get_stoichiometric_matrix(model: cobra.Model) -> sparse.csr_matrix:
    """Return sparse stoichiometric matrix, compatible across COBRApy versions."""
    S = cobra.util.array.create_stoichiometric_matrix(model)
    if not sparse.issparse(S):
        S = sparse.csr_matrix(S)
    return S.tocsr()


# ==============================================================================
# Operator diagnostics
# ==============================================================================

def _report_operator(Delta: sparse.spmatrix, W_R: np.ndarray) -> None:
    diag = Delta.diagonal()
    print(
        "[Operator] shape={}  nnz={}  "
        "diag=[{:.3f}, {:.3f}]  W_R=[{:.4f}, {:.4f}]".format(
            Delta.shape, Delta.nnz,
            float(diag.min()), float(diag.max()),
            float(W_R.min()), float(W_R.max()),
        )
    )

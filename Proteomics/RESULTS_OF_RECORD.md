# Results of Record

This document defines the **canonical numerical outputs** of the repository — the
"results of record." Every quantitative figure in the manuscript is generated
from these CSV files. They are the authoritative computational outputs; if a
figure and a results-of-record CSV ever disagree, the CSV is correct and the
figure must be regenerated.

All filenames below were verified against the repository tree. Paths are relative
to the repository root.

---

## What "result of record" means

Each experiment produces:
- a **raw CSV** — one row per trial / cell line / fold (the granular data),
- a **summary CSV** — aggregated statistics (means, SDs) used directly by the figure.

A figure may be regenerated from the summary CSV alone, but the raw CSV is
retained so the aggregation itself is reproducible and auditable. For Figure 9,
the raw → summary aggregation has been independently verified (per-trial rows
reproduce the published advantages to machine precision).

---

## Experiment A — topology-only panel design (single representative cell line)

| Item | File |
|------|------|
| Analysis script | `nci60/nci60_script/layer1_masking_benchmark.py` |
| Raw CSV | `nci60/nci60_csv/Raw experiment CSVs/layer1_experiment_A_raw.csv` |
| Summary CSV | `nci60/nci60_csv/Core manuscript-result CSVs/layer1_experiment_A_summary.csv` |
| Summary rows | `nci60/nci60_csv/Raw experiment CSVs/layer1_experiment_A_summary_rows.csv` |
| Plotting script | `make_figure9.py` (panel context) |
| Manuscript figure | Figure 9 (supporting; establishes the panel-design principle on the topology-only operator) |

Establishes that geometry-aware panels beat random even without proteomics.

---

## Experiment A′ — proteomics-informed panel design (57 cell lines)

| Item | File |
|------|------|
| Analysis script | `nci60/nci60_script/layer1_experiment_Aprime.py` |
| Raw CSV | `nci60/nci60_csv/Raw experiment CSVs/layer1_experiment_Aprime_raw.csv` |
| Summary CSV | `nci60/nci60_csv/Core manuscript-result CSVs/layer1_experiment_Aprime_summary.csv` |
| Summary rows | `nci60/nci60_csv/Raw experiment CSVs/layer1_experiment_Aprime_summary_rows.csv` |
| Plotting script | `make_figure9.py` (panels A, B) |
| Manuscript figure | Figure 9A, 9B |

The primary multi-cell-line benchmark. Advantage +0.152 → +0.499 across panel
sizes 50 → 500. Verified: raw rows aggregate to these summary advantages to 1e-16.

---

## Experiment B — operator-type comparison under random masking (57 cell lines)

| Item | File |
|------|------|
| Analysis script | `nci60/nci60_script/layer1_masking_benchmark.py` |
| Raw CSV | `nci60/nci60_csv/Raw experiment CSVs/layer1_experiment_B_raw.csv` |
| Summary CSV | `nci60/nci60_csv/Core manuscript-result CSVs/layer1_experiment_B_summary.csv` |
| Plotting script | `make_figure9.py` (panel D) |
| Manuscript figure | Figure 9D |

Near-null result (|mean Δ| < 0.006): panel selection dominates operator type.

---

## Spectral Scan — proteomics-induced operator deformation (57 cell lines)

| Item | File |
|------|------|
| Analysis script | `nci60/nci60_script/layer1_spectral_scan.py` |
| Baseline eigenvalues | `nci60/nci60_csv/Core manuscript-result CSVs/layer1_topology_baseline_eigs.csv` |
| Summary CSV | `nci60/nci60_csv/Core manuscript-result CSVs/layer1_spectral_scan_summary.csv` |
| Plotting script | `make_figure9.py` (panel C) |
| Manuscript figure | Figure 9C |

Bimodal finding: eigenvalue-magnitude deviation median 7.32; eigenvector subspace
distance median 0.000 / mean 0.140; 10 of 57 cell lines rotate. The
`layer1_topology_baseline_eigs.csv` documents the near-null baseline eigenvalues
(~1e-10) that motivated the corrected K_EIGS=50 / subspace-distance metric.

---

## Coverage bridge — SWATH → Human1 (reaction level)

| Item | File |
|------|------|
| Analysis | GPR evaluation against Human-GEM (see `scripts/proteomics_weighting.py`) |
| Coverage CSV | `nci60/nci60_csv/Core manuscript-result CSVs/human1_reaction_swath_coverage.csv` |
| Manuscript figure | Methods / Section 2.9 text (43.9% coverage) |

12,971 reactions; 8,043 GPR-linked; 3,533 SWATH-covered → 43.9%. Verified by
recomputation from the `has_gpr` / `has_swath_overlap` columns.

---

## Real Cohort — three-operator classification (ST003506)

| Item | File |
|------|------|
| Analysis script | `scripts/real_cohort.py` (patched version adds `--results_csv`) |
| Inputs | `Table/ST003506_AN005756.txt`, `Human-GEM-main/model/Human-GEM.xml`, `data/proteomics/cptac_breast_tumor_only.tsv`, `data/hgnc_complete_set.txt` |
| Summary CSV (per operator) | `cohort_{baseline,proteomics,permuted}_summary.csv` ‡ |
| Raw CSV (per operator) | `cohort_{baseline,proteomics,permuted}_folds.csv` ‡ |
| Permutation null (per operator) | `cohort_{baseline,proteomics,permuted}_permnull.csv` ‡ |
| Plotting script | `make_figure8.py` (reads the CSVs above only) |
| Manuscript figure | Figure 8 |

‡ These CSVs are produced by running the **patched** `real_cohort.py` with
`--results_csv` (the patch is provided; once merged, the script writes them).
The current repository `real_cohort.py` computes the values (AUROC 0.941 / 0.923
/ 0.941, perm p = 0.001) but prints them; the patch persists them so Figure 8
becomes reproducible from CSVs rather than hard-coded arrays.

---

## NCI-60 Layer 1 — full benchmark (umbrella)

The NCI-60 Layer 1 benchmark is the umbrella over Experiments A, A′, B, the
Spectral Scan, and the Coverage bridge above. Its complete pipeline and the
operator/weight-bridge inputs are documented in `README_NCI60_LAYER1.md`.

| Item | File(s) |
|------|---------|
| Pipeline scripts | `nci60/nci60_script/layer1_build_S.py`, `layer1_spectral_scan.py`, `layer1_masking_benchmark.py`, `layer1_experiment_Aprime.py` |
| Operator/weight inputs | `nci60/nci60_csv/Operator : weight bridge CSVs/nci60_reaction_weights_norm_simple.csv`, `nci60_reaction_weights_aligned.csv`, `nci60_reaction_observed_mask.csv` |
| Clean processed inputs | `nci60/nci60_csv/Clean processed input CSVs/*.csv` |
| Result CSVs | all six in `nci60/nci60_csv/Core manuscript-result CSVs/` |
| Manuscript figure | Figure 9 |

---

## Quick index: which CSV backs which figure

| Figure | Results-of-record CSV(s) |
|--------|--------------------------|
| Fig 8 | `cohort_{op}_summary.csv`, `_folds.csv`, `_permnull.csv` ‡ |
| Fig 9A/9B | `layer1_experiment_Aprime_summary.csv` |
| Fig 9C | `layer1_spectral_scan_summary.csv` (+ `layer1_topology_baseline_eigs.csv`) |
| Fig 9D | `layer1_experiment_B_summary.csv` |
| Fig S1 | `toy_model_outputs/01–09_*.csv` (generated by `toy_model_worked_example.py`) |
| Coverage (text) | `human1_reaction_swath_coverage.csv` |

Figures 0, 1, 2 are conceptual schematics with no results of record.
Figures 3, 4, 5, 6, 7 require Human1 analysis CSVs not yet committed (see
FIGURE_PROVENANCE.md).

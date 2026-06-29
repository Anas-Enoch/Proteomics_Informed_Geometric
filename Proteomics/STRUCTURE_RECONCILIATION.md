# Structure Reconciliation Notes

The READMEs are written to match the **proposed** directory structure. Four naming
details differ from the current working artifacts. None is a computational problem — each
is a rename or a small flag change so the committed filenames match the structure the
READMEs (and reviewers) expect. Resolve these when you lay down the repo.

## 1. real_cohort/ — Figure 8 CSV names
**Proposed:** `figure8_summary.csv`, `figure8_cv_results.csv`, `figure8_permutation_null.csv`
**Actual (patched `real_cohort.py --results_csv`):** writes per operator —
`cohort_<op>_summary.csv`, `cohort_<op>_folds.csv`, `cohort_<op>_permnull.csv`
(op ∈ baseline, proteomics, permuted), and `make_figure8.py --stem` reads those.

The proposed structure flattens 3 operators × 3 kinds into 3 files. Two options:
- **(a)** Keep the real per-operator files and update the figures/ + real_cohort/ README
  to list `cohort_{baseline,proteomics,permuted}_{summary,folds,permnull}.csv` (9 files).
- **(b)** Add a tiny post-step that concatenates the three operators into the three
  proposed files (`figure8_summary.csv` with an `operator` column, etc.) and point
  `make_figure8.py` at those.
Recommended: **(a)** — it matches the verified plotting contract with no code change.

## 2. nci60/scripts/ — experiment drivers
**Proposed:** `layer1_experiment_A.py`, `layer1_experiment_Aprime.py`, `layer1_experiment_B.py`
**Actual:** `layer1_masking_benchmark.py` writes Experiment A **and** B;
`layer1_experiment_Aprime.py` writes A'. Standalone `layer1_experiment_A.py` /
`_B.py` do not currently exist.
Resolve by either splitting `layer1_masking_benchmark.py` into A/B drivers, or keeping
`layer1_masking_benchmark.py` and updating the nci60 README to name it (the README as
written already lists `layer1_masking_benchmark.py` as the A/B producer).

## 3. toy_example/ — CSV filenames
**Proposed:** `toy_S.csv`, `toy_WR.csv`, `toy_laplacian.csv`, `toy_eigenpairs.csv`,
`toy_panel_trace.csv`, `toy_psi.csv`
**Actual:** the toy script currently writes 9 CSVs under `toy_model_outputs/` with
ordinal prefixes (`01_..09_`). Rename the toy script's outputs to the proposed names (and
drop the subfolder), or update the toy_example README to the actual names.

## 4. figureS1 filename
**Proposed:** `figureS1_toy_example.pdf`
**Actual:** `figureS1_toy_model.pdf`
Rename the output (one-line change in the toy plotting call) or update the figures README.

---
Once these four are reconciled, every path named in every README resolves, and the
consistency audit returns 0 broken references.

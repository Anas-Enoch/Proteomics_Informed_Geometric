# Reproduction Protocol

This document gives the exact commands to reproduce every quantitative figure
from a fresh clone. Each step is labelled with a **verification status**:

- **VERIFIED** — executed during this audit; confirmed to produce output.
- **RUNNABLE (not executed here)** — all inputs are committed to the repo and the
  command is correct, but the step was not run in the audit environment (e.g.
  it needs `cobra` + the 43 MB Human-GEM model and minutes-to-hours of compute).
  It is *not* marked PASS because that would misrepresent an unrun step.

Honesty note: a genuine "fresh-clone PASS for every figure" can only be claimed
by someone who runs the heavy Human1 steps end-to-end. This protocol is written
so that you (or a reviewer) can do exactly that; the audit verified the
lightweight steps and the data-availability of the heavy ones.

---

## 0. Environment

```bash
git clone https://github.com/Anas-Enoch/Proteomics_Informed_Geometric.git
cd Proteomics_Informed_Geometric
python3 -m venv venv && source venv/bin/activate
pip install -r requirements.txt   # cobra numpy scipy pandas scikit-learn matplotlib networkx
```
**Status: RUNNABLE** — `requirements.txt` is present and complete; core imports
(numpy/scipy/pandas/sklearn/matplotlib) **VERIFIED** importable. `cobra` install
**RUNNABLE (not executed here)**.

Committed inputs confirmed present (no external download needed for the model):
- `Human-GEM-main/model/Human-GEM.xml` (43 MB) ✅
- `Table/ST003506_AN005756.txt` (cohort) ✅
- `data/proteomics/cptac_breast_tumor_only.tsv` ✅
- `data/hgnc_complete_set.txt` ✅

---

## 1. NCI-60 Layer 1 — Figure 9

The result CSVs are **already committed**, so the figure reproduces with no heavy
compute:
```bash
python make_figure9.py \
  --csvdir "nci60/nci60_csv/Core manuscript-result CSVs" \
  --out figures/figure9.pdf
```
**Status: VERIFIED** — regenerated `figure9.pdf` from committed CSVs during this
audit. Advantages +0.152/+0.335/+0.388/+0.499 and the bimodal spectral panel
reproduce exactly.

To regenerate the CSVs themselves from scratch (optional, heavier):
```bash
python nci60/nci60_script/layer1_build_S.py
python nci60/nci60_script/layer1_spectral_scan.py
python nci60/nci60_script/layer1_masking_benchmark.py
python nci60/nci60_script/layer1_experiment_Aprime.py
```
**Status: RUNNABLE (not executed here)** — needs `cobra` + Human-GEM.

---

## 2. Real cohort — Figure 8

```bash
mkdir -p results
for op in baseline proteomics permuted; do
  python scripts/real_cohort.py \
    --table Table/ST003506_AN005756.txt \
    --model Human-GEM-main/model/Human-GEM.xml \
    --operator $op \
    $( [ "$op" != baseline ] && echo "--proteomics data/proteomics/cptac_breast_tumor_only.tsv --hgnc data/hgnc_complete_set.txt" ) \
    --results_csv results/cohort.csv
done
python scripts/make_figure8.py --stem results/cohort --out figures/figure8.pdf
```
**Status: RUNNABLE (not executed here)** — needs `cobra` + Human-GEM. The
plotting contract (`make_figure8.py` reading the three CSV kinds) was **VERIFIED**
against schema-conformant CSVs during the audit. Requires the patched
`real_cohort.py` (adds `--results_csv`).

---

## 3. Toy worked example — Figure S1

```bash
python scripts/toy_model_worked_example.py    # writes toy_model_outputs/*.csv
python scripts/make_figureS1.py               # reads them
```
**Status: VERIFIED** — runs end-to-end with no external data; wrote 10 CSVs and
the figure. Deterministic (seed = 0).

---

## 4. Human1 masking / panel / PSI / rewiring — Figures 3, 4, 6, 7

These require the corrected analysis scripts to be run on Human-GEM to emit CSVs,
after which the (CSV-only) plotting scripts produce the figures.
```bash
# analysis (emit CSVs) — to be wired per FIGURE_PROVENANCE.md schemas
python scripts/<analysis>.py --model Human-GEM-main/model/Human-GEM.xml --out_csv <figureN>_summary.csv
# plot from CSV
python scripts/plotting/make_figure<N>.py --csv <figureN>_summary.csv --out figures/figure<N>.pdf
```
**Status: BLOCKED** — these figures do **not** yet have committed CSVs, and the
legacy plotting scripts (`make_figures.py`, `PSI_preservation.py`,
`Biological_alignment.py`) use synthetic values. See REPOSITORY_INTEGRITY_AUDIT.md.
Until the analysis is run and CSVs committed, these are **FAIL** for
reproducibility.

---

## Reproduction status summary

| Figure | Command status |
|--------|----------------|
| 0,1,2 (schematics) | VERIFIED (rebuilt, no data) |
| 9 | VERIFIED (from committed CSVs) |
| S1 | VERIFIED (self-contained) |
| 8 | RUNNABLE (data committed; plotting contract verified) |
| 3,4,6,7 | BLOCKED (CSVs not committed / legacy plotting) |

---

## What a reviewer would still hit

1. `cobra` model load is memory-heavy (43 MB SBML); expect ~minutes to parse.
2. The permutation test in `real_cohort.py` runs 1000 permutations × 5-fold CV per
   operator — budget compute time accordingly.
3. Figures 3, 4, 6, 7 cannot currently be reproduced from the public repo because
   their CSVs are not committed; this is the last blocking gap.

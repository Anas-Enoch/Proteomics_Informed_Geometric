# Proteomics-Informed Geometric Framework for Identifiability and Panel Design in Genome-Scale Metabolic Networks

Author: **Anas Enoch, MD**  
Mohammed VI University of Health Sciences (UM6SS), Casablanca, Morocco

---

## Overview

This repository contains the code accompanying the manuscript:

**"A Proteomics-Informed Geometric Framework for Identifiability and Panel Design in Genome-Scale Metabolic Networks."**

The framework constructs a proteomics-informed metabolite operator derived from the stoichiometric structure of genome-scale metabolic models.

The metabolite Laplacian:

    Δ_M = W_M^(1/2) S W_R S^T W_M^(1/2)

defines a condition-specific metabolic geometry that enables:

- mechanistic identifiability analysis
- geometry-preserving metabolite panel design
- robustness analysis under partial metabolite observability
- downstream disease classification evaluation

---

## Repository Structure

```
scripts/          figure generation and analysis scripts
figures/          generated figures used in the manuscript
Human-GEM-main/   Human genome-scale metabolic model
Table/            metabolomics cohort data
data/proteomics/  external proteomics data (see DATA_SOURCES.md)
data/             auxiliary data files (HGNC mapping table)
references.bib    bibliography
FIGURE_CAPTIONS.md
DATA_SOURCES.md
requirements.txt
```

---

## Core Components

- stoichiometric network loading via COBRApy
- proteomics-informed reaction weighting via GPR rule evaluation
- metabolite Laplacian construction (baseline, proteomics-informed, permuted)
- spectral geometry analysis
- greedy metabolite panel selection
- robustness experiments (permutation tests, OR aggregation sensitivity)
- proof-of-concept disease classification

---

## Dependencies

Python >= 3.10

Install dependencies:

```bash
pip install -r requirements.txt
```

---

## Example Execution

### Baseline cohort experiment

```bash
./venv/bin/python scripts/real_cohort.py \
  --table Table/ST003506_AN005756.txt \
  --model Human-GEM-main/model/Human-GEM.xml \
  --operator baseline \
  --out results_baseline.pdf
```

### Proteomics-weighted cohort experiment

```bash
./venv/bin/python scripts/real_cohort.py \
  --table Table/ST003506_AN005756.txt \
  --model Human-GEM-main/model/Human-GEM.xml \
  --operator proteomics \
  --proteomics data/proteomics/cptac_breast_tumor_only.tsv \
  --hgnc data/hgnc_complete_set.txt \
  --out results_proteomics.pdf
```

### Permuted-weight cohort experiment

```bash
./venv/bin/python scripts/real_cohort.py \
  --table Table/ST003506_AN005756.txt \
  --model Human-GEM-main/model/Human-GEM.xml \
  --operator permuted \
  --proteomics data/proteomics/cptac_breast_tumor_only.tsv \
  --hgnc data/hgnc_complete_set.txt \
  --out results_permuted.pdf
```

### Operator rewiring analysis

```bash
python scripts/Fig7_rewiring.py \
  --model Human-GEM-main/model/Human-GEM.xml
```

---

## Reproducibility and Methodological Clarification

The real-cohort experiment provided in `scripts/real_cohort.py` reproduces the
proof-of-concept classification analysis reported in the manuscript using the
publicly available Metabolomics Workbench dataset `ST003506` together with the
Human-GEM genome-scale metabolic model.

The cohort contains **43 samples** (**31 cases**, **12 controls**) and
**49 measured metabolites**, of which **36** were harmonized to Human-GEM
metabolite identifiers after preprocessing and identifier mapping.

Three operator modes are implemented and reproducible from the repository:

| Mode | Description |
|------|-------------|
| `baseline` | Stoichiometric operator with W_R = I (Delta_M = S S^T) |
| `proteomics` | CPTAC-weighted operator using external breast-tumor proteomics mapped to Human-GEM genes via GPR rules and HGNC identifier remapping |
| `permuted` | Distribution-matched null obtained by randomly permuting reaction weights while preserving their empirical distribution |

### Observed results

| Operator | Mean AUROC | Std | Permutation p |
|----------|-----------|-----|---------------|
| Baseline (W_R = I) | 0.941 | 0.103 | 0.001 |
| CPTAC proteomics-informed | 0.923 | 0.109 | 0.001 |
| Permuted weights (null) | 0.941 | 0.103 | 0.001 |

All three conditions substantially exceeded the permutation null
(mean null AUROC = 0.50 ± 0.14), confirming that stoichiometric operator
geometry captures genuine disease-relevant structure in the spectral embedding.

The proteomics-weighted operator did not outperform the baseline in this cohort.
The permuted operator matched baseline performance, indicating that the
classification signal is dominated by stoichiometric coupling and metabolite
abundances rather than by the external proteomics prior at this panel size.

Because matched proteomics measurements are not available for ST003506, this
cohort experiment should be interpreted as a proof-of-concept evaluation of
operator-derived geometry under partial observability rather than as a definitive
matched multi-omics validation. The effect of proteomics-informed operator
deformation on geometry preservation is investigated separately in the masking,
PSI, and rewiring analyses described in the manuscript.

---

## Data Availability

All datasets used in this study are publicly available.  
See `DATA_SOURCES.md` for full details and download instructions.

---

## Citation

If you use this framework, please cite the associated manuscript.

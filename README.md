# Proteomics-Informed Geometric Framework for Identifiability and Panel Design in Genome-Scale Metabolic Networks

Author: **Anas Enoch, MD**  
Mohammed VI University of Health Sciences (UM6SS), Casablanca, Morocco

---

## Overview

This repository contains the code accompanying the manuscript:

**“A Proteomics-Informed Geometric Framework for Identifiability and Panel Design in Genome-Scale Metabolic Networks.”**

The framework constructs a proteomics-informed metabolite operator derived from the stoichiometric structure of genome-scale metabolic models.

The metabolite Laplacian:

Δ_M = W_M^(1/2) S W_R Sᵀ W_M^(1/2)

defines a condition-specific metabolic geometry that enables:

• mechanistic identifiability analysis  
• geometry-preserving metabolite panel design  
• robustness analysis under partial metabolite observability  
• downstream disease classification evaluation

---

## Repository Structure
scripts/
figure generation and analysis scripts

figures/
generated figures used in the manuscript

Human-GEM-main/
Human genome-scale metabolic model

ST003506_AN005756.txt
breast cancer metabolomics cohort

references.bib
bibliography

FIGURE_CAPTIONS.md
figure captions

DATA_SOURCES.md
dataset references

---

## Core Components

The repository implements:

• stoichiometric network loading via COBRApy  
• proteomics-informed reaction weighting  
• metabolite Laplacian construction  
• spectral geometry analysis  
• greedy metabolite panel selection  
• robustness experiments (permutation tests, OR aggregation sensitivity)  
• proof-of-concept disease classification

---

## Dependencies

Python ≥ 3.10

Required packages:

cobra  
numpy  
scipy  
pandas  
scikit-learn  
matplotlib

Install dependencies:

pip install cobra numpy scipy pandas scikit-learn matplotlib

---

## Example Execution

Example operator analysis:

python scripts/Fig7_rewiring.py --model Human-GEM-main/model/Human-GEM.xml

Example cohort classification:

python scripts/classify_ST003506_operator.py \
--table Table/ST003506_AN005756.txt \
--model Human-GEM-main/model/Human-GEM.xml



---

## Data Availability

All datasets used in this study are publicly available.

See:

DATA_SOURCES.md

---

## Citation

If you use this framework please cite the associated manuscript.

### 🔴 <h3 style="color:red;">Reproducibility and methodological clarification</h3>

The real-cohort experiment provided in scripts/real_cohort.py reproduces the proof-of-concept classification analysis reported in Section 2.8 of the manuscript using the publicly available Metabolomics Workbench dataset ST003506 together with the Human-GEM genome-scale metabolic model. The cohort contains 43 samples (31 breast-cancer-related lymphedema cases and 12 controls) with 49 measured metabolites, of which 36 could be harmonized to Human-GEM metabolite identifiers after preprocessing and identifier mapping. The metabolite operator used in this experiment is the unweighted stoichiometric Laplacian
\Delta_M = S S^\top
constructed from the Human-GEM stoichiometric matrix. Condition-specific proteomics weights W_R^{(c)} are not applied here because matched proteomics measurements are not available for ST003506; the classification experiment therefore evaluates the intrinsic stoichiometric geometry under partial metabolite observability, while the effects of proteomics-informed deformation are investigated separately in the masking and separability analyses described in earlier sections of the manuscript. Spectral features are obtained by projecting metabolite abundance vectors onto the leading non-trivial eigenvectors of the restricted operator. Classification performance is evaluated using repeated stratified cross-validation with fold-specific imputation and standardization to avoid data leakage, and statistical significance of the observed ROC–AUC is assessed using a permutation test with 1000 label permutations. The current implementation yields a mean AUROC ≈ 0.94 with a permutation-based significance of p ≈ 0.001, indicating that the discriminative signal captured by the operator-based embedding is unlikely to arise from class imbalance or feature dimensionality alone. All scripts, data inputs, and model files included in this repository allow the analysis and figures to be reproduced directly.

<h3 style="color:red;">One-command reproduction</h3>

The real-cohort classification experiment can be reproduced directly from this repository using a single command.

First install the required Python environment:

pip install -r requirements.txt

Then run the experiment:

./venv/bin/python scripts/real_cohort.py \
--table Table/ST003506_AN005756.txt \
--model Human-GEM-main/model/Human-GEM.xml

This script will:

• load the ST003506 metabolomics cohort  
• harmonize metabolite identifiers to the Human-GEM namespace  
• construct the metabolite Laplacian operator  
• compute the spectral embedding of the observed metabolite subspace  
• perform repeated stratified cross-validation classification  
• run a 1000-permutation significance test  
• generate the figure `Fig6_real_cohort_classification.pdf`

Typical output (reproducible with the provided dataset):

Samples kept: 43 (Control=12, Lymphedema=31)  
Measured metabolites: 49  
Mapped to model: 36  
Unmapped: 13  

AUROC: mean ≈ 0.94  
Permutation test p-value ≈ 0.001  

The resulting figure summarizes AUROC and accuracy distributions across cross-validation splits and corresponds to the real-cohort experiment reported in the manuscript.
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
--table ST003506_AN005756.txt \
--model Human-GEM-main/model/Human-GEM.xml

---

## Data Availability

All datasets used in this study are publicly available.

See:

DATA_SOURCES.md

---

## Citation

If you use this framework please cite the associated manuscript.
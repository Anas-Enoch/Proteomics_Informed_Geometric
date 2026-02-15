# Proteomics-Informed Geometric Framework for Identifiability and Panel Design in Genome-Scale Metabolic Networks

## Overview

This repository contains code and analysis pipelines for constructing a proteomics-informed, metabolite-centric geometric representation of genome-scale metabolic models (GEMs).

The framework introduces a condition-specific metabolite Laplacian derived from:

- Stoichiometric structure
- Gene–protein–reaction (GPR) rules
- Gene-level proteomics
- Metabolite reliability weighting

Partial metabolomic observability is modeled as an operator restriction problem, and identifiability is defined as stability of spectral geometry under metabolite masking.

This repository reproduces all figures and ablation analyses presented in the manuscript.

---

## Core Concept

For each biological condition \( c \), we construct:

\[
\Delta_M^{(c)} = W_M^{1/2} S W_R^{(c)} S^\top W_M^{1/2}
\]

where:

- \( S \) = stoichiometric matrix
- \( W_R^{(c)} \) = reaction weight matrix derived from proteomics
- \( W_M \) = metabolite reliability weighting

The induced spectral geometry is analyzed using:

- Low-frequency eigenvalues
- Heat trace
- Diffusion distances

---

## Repository Structure
├── data/                 # Proteomics and model data
├── scripts/              # Figure generation and analysis scripts
├── results/              # Output files
├── figures/              # Final manuscript figures
├── README.md
├── FIGURE_CAPTIONS.md
---

## Key Analyses

- Operator construction from stoichiometry + proteomics
- OR-aggregator parameter sensitivity
- Nonlinear saturation mapping robustness
- Permutation ablation control
- Geometry-preserving metabolite panel selection
- Cohort-level classification demonstration

---

## Reproducibility

All scripts are written in Python and rely on:

- numpy
- scipy
- matplotlib
- COBRApy

To reproduce Figure 5:

```bash
python scripts/figure5_permutation_robustness.py
python scripts/make_figures.py
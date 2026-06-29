# Proteomics-Informed Geometric Framework for Identifiability and Panel Design in Genome-Scale Metabolic Networks

> **All quantitative manuscript figures are generated from committed CSV outputs produced by public analysis scripts. Conceptual workflow schematics are explicitly labelled as non-data figures.**

---

## What is this repository?

Computational companion to the manuscript *"A Proteomics-Informed Geometric Framework for
Identifiability and Panel Design in Genome-Scale Metabolic Networks."* It contains the
analysis code, public input data, committed numerical outputs ("results of record"), and
figure-generation scripts needed to reproduce every quantitative figure.

The framework builds a metabolite-centric operator
$\Delta_M = W_M^{1/2}\, S\, W_R^{(c)}\, S^\top\, W_M^{1/2}$
from genome-scale stoichiometry $S$ and proteomics-derived reaction weights $W_R^{(c)}$,
and studies its low-frequency spectral geometry for metabolite identifiability and panel
design. Evidence is organised in two layers: model-scale (Human-GEM / Human1 and the
NCI-60 cell-line benchmark) and cohort (ST003506 breast-cancer metabolomics).

## How do I install?

```bash
git clone https://github.com/Anas-Enoch/Proteomics_Informed_Geometric.git
cd Proteomics_Informed_Geometric
python3 -m venv venv && source venv/bin/activate
pip install -r requirements.txt
```
Dependencies: `cobra`, `numpy`, `scipy`, `pandas`, `scikit-learn`, `matplotlib`, `networkx`.
The Human-GEM SBML model and all public input data are included, so the pipeline runs
from a fresh clone without external downloads.

## How do I reproduce every figure?

Each quantitative figure is produced in two stages: an **analysis script** that writes
CSV results of record, and a **plotting script** that reads only those CSVs.

```bash
# Model-scale figures (3-7)
cd human1_fig_scripts
bash run_human1_figures.sh          # -> ../human1_fig_results/*.csv and ../figures/figure3-7.pdf

# Cohort figure (8)
cd ../real_cohort
for op in baseline proteomics permuted; do
  python real_cohort.py --table ../data/ST003506_AN005756.txt \
    --model ../human_gem/model/Human-GEM.xml --operator $op \
    $( [ "$op" != baseline ] && echo "--proteomics ../data/cptac_breast_tumor_only.tsv --hgnc ../data/hgnc_complete_set.txt" ) \
    --results_csv figure8.csv
done
python make_figure8.py --stem figure8 --out ../figures/figure8.pdf

# NCI-60 validation figure (9)
cd ../nci60
python ../human1_fig_scripts/make_figure9.py --csvdir csv --out ../figures/figure9.pdf

# Toy worked example (S1)
cd ../toy_example
python toy_pipeline.py
```
Per-figure commands and verification status: [`REPRODUCTION_PROTOCOL.md`](REPRODUCTION_PROTOCOL.md).

## Where are the results-of-record CSVs?

| Figure(s) | Location |
|-----------|----------|
| 3-7 (model scale) | `human1_fig_results/` |
| 8 (cohort) | `real_cohort/figure8_*.csv` |
| 9 (NCI-60) | `nci60/csv/` |
| S1 (toy) | `toy_example/toy_*.csv` |

[`RESULTS_OF_RECORD.md`](RESULTS_OF_RECORD.md) maps each experiment to its scripts and
CSVs; [`FIGURE_PROVENANCE.md`](FIGURE_PROVENANCE.md) gives the full provenance chain.

## Where are the public data sources?

In `data/` and `human_gem/`, documented in [`DATA_SOURCES.md`](DATA_SOURCES.md):
Human-GEM (`human_gem/`), CPTAC breast proteomics (`data/cptac_breast_tumor_only.tsv`),
ST003506 cohort (`data/ST003506_AN005756.txt`), HGNC mapping
(`data/hgnc_complete_set.txt`), NCI-60/CellMiner SWATH proteomics (processed,
`nci60/csv/`).

## Which figures are conceptual vs computed?

**Conceptual schematics (no numerical data; non-data figures):**
Figure 1 (operator-centric workflow), Figure 2 (Dirac-to-Laplacian construction).

**Computed (every value from a committed CSV):**
Figure 3 (geometry preservation), Figure 4 (panel design), Figure 5 (robustness),
Figure 6 (CH-PSI pathway separability), Figure 7 (CPTAC operator deformation),
Figure 8 (ST003506 cohort), Figure 9 (NCI-60 validation), Figure S1 (toy example).

## Repository documentation

| Document | Purpose |
|----------|---------|
| [`RESULTS_OF_RECORD.md`](RESULTS_OF_RECORD.md) | Canonical CSV outputs; which CSV backs which figure |
| [`FIGURE_PROVENANCE.md`](FIGURE_PROVENANCE.md) | Per-figure provenance; conceptual vs computed |
| [`REPRODUCTION_PROTOCOL.md`](REPRODUCTION_PROTOCOL.md) | Commands to reproduce every figure from a clone |
| [`DATA_SOURCES.md`](DATA_SOURCES.md) | Provenance and licensing of external datasets |
| [`REPOSITORY_INTEGRITY_AUDIT.md`](REPOSITORY_INTEGRITY_AUDIT.md) | Per-script integrity audit |

Each subdirectory has its own `README.md`.

# Data Sources

This repository does not host large external datasets.
All datasets used in this study are publicly available from their original repositories.

---

## Genome-scale metabolic model

**Human-GEM (Human1)**

Repository: https://github.com/SysBioChalmers/Human-GEM

Reference:
Robinson JL et al.
An atlas of human metabolism.
*Science Signaling* (2020)

File used in this repository:

    Human-GEM-main/model/Human-GEM.xml

This file provides the stoichiometric matrix and gene--protein--reaction (GPR)
associations used to construct the metabolite operator.

---

## Metabolomics cohort dataset

**Metabolomics Workbench — Study ST003506**

Repository: https://www.metabolomicsworkbench.org  
Study ID: ST003506  
Analysis ID: AN005756

Description:
Serum NMR metabolomics dataset from breast cancer-related lymphedema patients
and healthy controls (43 samples: 31 cases, 12 controls; 49 measured metabolites).

File included in this repository:

    Table/ST003506_AN005756.txt

This dataset is used for the proof-of-concept classification experiment
described in Section 2.8 of the manuscript.

---

## Proteomics dataset

**CPTAC Breast Cancer Proteomics (PDC Study PDC000120)**

Repository: https://pdc.cancer.gov  
Study: CPTAC2 Breast Cancer — Proteome  
PDC Study ID: PDC000120

Description:
Gene-level tumor protein abundances from the CPTAC breast cancer cohort.
Used to construct condition-specific reaction weights W_R^(c) via GPR rule
evaluation for the proteomics-informed operator mode.
These measurements originate from a different patient population than ST003506
and represent a population-level condition prior, not sample-matched proteomics.

File used in this repository:

    data/proteomics/cptac_breast_tumor_only.tsv

To reproduce the proteomics operator mode, download the protein abundance
matrix from the PDC portal and pre-process to a two-column TSV with columns
`gene_symbol` and `tumor_mean`.

---

## Gene identifier mapping table

**HGNC Complete Gene Set**

Source: https://www.genenames.org/download/statistics-and-files/  
File: Complete HGNC dataset (TXT format)

Description:
Used to remap CPTAC gene symbols to Human-GEM Ensembl gene identifiers,
which are required for GPR rule evaluation. Human-GEM gene IDs are Ensembl
identifiers (e.g. ENSG00000141510) while CPTAC proteomics uses HGNC gene
symbols (e.g. TP53).

File used in this repository:

    data/hgnc_complete_set.txt

One-time download:

    curl -o data/hgnc_complete_set.txt \
      "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"

This file is not included in the repository due to size (~16 MB) and is
downloaded separately as described above.

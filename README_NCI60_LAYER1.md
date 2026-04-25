# NCI-60 Layer 1 External Validation Pipeline

## Overview

This repository contains a complete **Layer 1 external validation pipeline** for the manuscript:

**A Proteomics-Informed Geometric Framework for Identifiability and Panel Design in Genome-Scale Metabolic Networks**

The goal of this pipeline is to test, in a real multi-omics benchmark, whether:

1. **proteomics-informed operators are spectrally distinct** from the topology-only baseline \(SS^\top\),
2. **geometry-aware panel design** outperforms random panel selection,
3. and **proteomics-informed geometry-aware selection** outperforms random selection on the proteomics-informed operators themselves.

The benchmark uses:

- **Human-GEM / Human1** as the metabolic scaffold
- **NCI-60 / CellMiner SWATH proteomics** as the external real-data layer
- **57 aligned cell lines**
- **43.9% GPR-linked reaction coverage** by SWATH-informed reaction weighting

---

## Project structure

### Scripts

The main Layer 1 scripts are stored in:

```bash
nci60/nci60_script/
```

They are:

```bash
nci60/nci60_script/layer1_build_S.py
nci60/nci60_script/layer1_spectral_scan.py
nci60/nci60_script/layer1_masking_benchmark.py
nci60/nci60_script/layer1_experiment_Aprime.py
```

### CSV outputs

The generated CSVs are organized in:

```bash
nci60/nci60_csv/
```

The categories visible in the repository are:

```bash
nci60/nci60_csv/Core manuscript-result CSVs
nci60/nci60_csv/Raw experiment CSVs
nci60/nci60_csv/Operator : weight bridge CSVs
nci60/nci60_csv/Clean processed input CSVs
```

### Important root-level files

Key output files in the project root include:

```bash
human1_reaction_swath_coverage.csv
layer1_spectral_scan_summary.csv
layer1_topology_baseline_eigs.csv
layer1_experiment_A_summary.csv
layer1_experiment_Aprime_summary.csv
layer1_experiment_B_summary.csv
layer1_experiment_A_raw.csv
layer1_experiment_Aprime_raw.csv
layer1_experiment_B_raw.csv
nci60_reaction_weights_norm_simple.csv
nci60_reaction_weights_aligned.csv
nci60_reaction_observed_mask.csv
nci60_clean_sample_table.csv
nci60_swa_gene_level_median.csv
nci60_prot_ensembl_normalized.csv
nci60_drug_panel_strict.csv
nci60_drug_panel_55.csv
```

---

## Biological and computational rationale

This Layer 1 benchmark does **not** replace the main Human1 masking experiments in the manuscript. It extends them.

The logic is:

- the manuscript already shows, inside Human1, that the operator-centric framework supports identifiability analysis and panel design;
- the NCI-60 benchmark asks whether this survives **real external multi-omics heterogeneity**.

The benchmark is therefore a **real-data externalization** of the operator framework.

---

## Environment setup

From the project root:

```bash
python3 -m venv venv
source venv/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install cobra pandas numpy scipy
```

If the scripts use additional packages already listed in `requirements.txt`, you can instead do:

```bash
source venv/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install -r requirements.txt
```

---

## Full Layer 1 pipeline

Run the scripts in this order.

### 1. Build aligned stoichiometric/operator inputs

```bash
source venv/bin/activate
python nci60/nci60_script/layer1_build_S.py
```

### What this step does

- loads `Human-GEM-main/model/Human-GEM.xml`
- builds the stoichiometric matrix \(S\)
- loads `nci60_reaction_weights_norm_simple.csv`
- aligns reaction weights to Human-GEM reaction order
- constructs an example proteomics-informed operator and the topology-only baseline
- saves aligned matrices and baseline eigensummaries

### Key outputs

```bash
human1_S_aligned.csv
human1_S_aligned.npz
human1_met_ids.csv
human1_rxn_ids_aligned.csv
nci60_reaction_weights_aligned.csv
example_lowfreq_eigenvalues.csv
example_lowfreq_eigenvalues_nontrivial.csv
layer1_topology_baseline_eigs.csv
```

---

### 2. Spectral scan across the 57 NCI-60 cell lines

```bash
source venv/bin/activate
python nci60/nci60_script/layer1_spectral_scan.py
```

### What this step does

For each of the 57 aligned NCI-60 cell lines:

- constructs the proteomics-informed operator  
  \[
  \Delta^{(c)} = S W_R^{(c)} S^\top
  \]
- compares its low-frequency spectrum to the topology-only operator  
  \[
  \Delta^{\mathrm{topo}} = SS^\top
  \]
- computes spectral deformation metrics across the full panel

### Key output

```bash
layer1_spectral_scan_summary.csv
```

### Main result

This is a **real positive result**.

Across **57/57** cell lines:

- **median relative spectral deviation from topology-only**: **0.564**
- **mean relative spectral deviation**: **0.720**
- range: **0.038–4.57**

### Interpretation

This means:

> Proteomics-informed reaction weighting induces substantial low-frequency deformation of the metabolic operator relative to the topology-only baseline across the NCI-60 panel.

This is important because it demonstrates that the proteomics-informed operator is **not collapsing back to \(SS^\top\)** in real data — the proteomics signal is actively deforming the mechanistic geometry in a cell-line-specific way.

---

### 3. Masking benchmark (Experiment A and Experiment B)

```bash
source venv/bin/activate
python nci60/nci60_script/layer1_masking_benchmark.py
```

### What this step does

This script runs **two different experiments**.

#### Experiment A — topology-only geometry-aware panel benchmark

- builds geometry-aware panels from the low-frequency topology-only operator
- compares them to random panels of equal size
- evaluates normalised spectral distortion under restriction

This tests the **panel-design principle itself**, independently of proteomics.

> **Note:** the `layer1_experiment_A_summary.csv` reflects a single representative cell line run.
> The full 57-cell-line benchmark is in Experiment A′ (Step 4).
> Experiment A establishes that the topology-only geometry supports principled panel design;
> Experiment A′ is the primary multi-cell-line benchmark.

#### Experiment B — operator-type comparison under random masking

- compares proteomics-informed vs topology-only operators under the **same random masks**
- this is a secondary diagnostic, not the core claim

### Key outputs

```bash
layer1_experiment_A_raw.csv
layer1_experiment_A_summary.csv
layer1_experiment_A_summary_rows.csv

layer1_experiment_B_raw.csv
layer1_experiment_B_summary.csv
layer1_experiment_B_summary_rows.csv
```

### Main result — Experiment A

Geometry-aware panels outperform random panels at all tested panel sizes:

- panel size **50** → advantage **+0.168**
- panel size **100** → advantage **+0.369**
- panel size **200** → advantage **+0.434**
- panel size **500** → advantage **+0.475**

where:

\[
\text{advantage} = \text{random mean distortion} - \text{geometry-aware distortion}
\]

So **positive** means geometry-aware beats random.

### Interpretation — Experiment A

This validates the claim that:

> Low-frequency metabolic operator geometry supports principled, non-random panel construction — even on the topology-only operator, without any proteomics information.

This result is important because it cannot be dismissed as circular: the topology-only operator has no proteomics component, yet geometry-aware panels derived from its spectral structure substantially outperform random panels.

### Interpretation — Experiment B

Experiment B is near-neutral.

The operator-type comparison under random masking does **not** show a strong systematic advantage of proteomics-informed over topology-only under blind random restriction (|mean Δ| < 0.006 at all panel sizes).

This is not a failure of the framework. It means:

> The main benefit of proteomics lies in enabling **structured panel design**, not in a generic guarantee that proteomics-informed operators are always more robust to arbitrary random subsampling.

This result directly motivates the geometry-aware panel design criterion: when a richer, more structured operator is used, the choice of which metabolites to measure becomes more consequential, not less.

---

### 4. Proteomics-informed Experiment A′

```bash
source venv/bin/activate
python nci60/nci60_script/layer1_experiment_Aprime.py
```

### What this step does

This is the primary multi-cell-line benchmark.

For each of the 57 cell lines:

- builds the **proteomics-informed operator**
- derives geometry-aware panels directly from its low-frequency eigenbasis
- compares those panels against random panels of equal size under normalised spectral distortion

### Key outputs

```bash
layer1_experiment_Aprime_raw.csv
layer1_experiment_Aprime_summary.csv
layer1_experiment_Aprime_summary_rows.csv
```

### Main result — Experiment A′

Proteomics-informed geometry-aware selection beats random at every tested panel size across all 57 cell lines:

- panel size **50** → mean advantage **+0.152** (median +0.175)
- panel size **100** → mean advantage **+0.335** (median +0.375)
- panel size **200** → mean advantage **+0.388** (median +0.426)
- panel size **500** → mean advantage **+0.499** (median +0.562)

Geometry-aware distortion improves monotonically with panel size:

- **1.128 → 1.025 → 1.002 → 0.903**

while random mean distortion remains elevated and approximately constant:

- **1.280 → 1.360 → 1.390 → 1.402**

### Interpretation — Experiment A′

This is the strongest Layer 1 result.

It shows that:

> In real multi-omics data across 57 cancer cell lines, proteomics-informed geometry-aware panel selection preserves the proteomics-informed operator geometry substantially better than random selection under equal panel budgets — and the advantage grows monotonically as panel budget increases.

The divergence between geometry-aware and random distortion as panel size grows is the key signal: geometry-aware panels get better with budget, random panels do not.

---

## Coverage bridge from SWATH proteomics to Human1

### Key file

```bash
human1_reaction_swath_coverage.csv
```

### Main numbers

- SWATH genes mapped to Human1 Ensembl IDs: **722**
- Human1 GPR-linked reactions: **8043**
- reactions with SWATH support: **3533**
- reaction coverage: **43.9%**

### Interpretation

The NCI-60 benchmark is therefore a **partial but real** proteomics-informed modulation of Human1. Unsupported reactions retain baseline coupling via confidence blending (see manuscript Section 4.4).

This is sufficient for external validation, but should be described honestly as:

> partial proteomic modulation of the GEM — 43.9% reaction coverage — not full proteome saturation of Human1.

---

## Minimal command-only pipeline

If everything is already prepared, the Layer 1 pipeline is:

```bash
source venv/bin/activate
python nci60/nci60_script/layer1_build_S.py
python nci60/nci60_script/layer1_spectral_scan.py
python nci60/nci60_script/layer1_masking_benchmark.py
python nci60/nci60_script/layer1_experiment_Aprime.py
```

Or chained:

```bash
source venv/bin/activate && \
python nci60/nci60_script/layer1_build_S.py && \
python nci60/nci60_script/layer1_spectral_scan.py && \
python nci60/nci60_script/layer1_masking_benchmark.py && \
python nci60/nci60_script/layer1_experiment_Aprime.py
```

---

## Which CSVs must be kept for reproducibility

### Core manuscript-result CSVs

These are mandatory:

```bash
human1_reaction_swath_coverage.csv
layer1_spectral_scan_summary.csv
layer1_topology_baseline_eigs.csv
layer1_experiment_A_summary.csv
layer1_experiment_Aprime_summary.csv
layer1_experiment_B_summary.csv
```

### Raw experiment CSVs

These should also be kept:

```bash
layer1_experiment_A_raw.csv
layer1_experiment_A_summary_rows.csv
layer1_experiment_Aprime_raw.csv
layer1_experiment_Aprime_summary_rows.csv
layer1_experiment_B_raw.csv
layer1_experiment_B_summary_rows.csv
```

### Operator / weight bridge CSVs

These are important for reproducibility of the Human1–NCI-60 bridge:

```bash
nci60_reaction_weights_norm_simple.csv
nci60_reaction_weights_aligned.csv
nci60_reaction_observed_mask.csv
```

### Clean processed inputs

These should be kept if you want the benchmark to be easy to rerun:

```bash
nci60_clean_sample_table.csv
nci60_swa_gene_level_median.csv
nci60_prot_ensembl_normalized.csv
nci60_drug_panel_strict.csv
nci60_drug_panel_55.csv
```

---

## Final interpretation

The Layer 1 benchmark now supports the following hierarchy of claims:

### Claim 1 — spectral externalization
Proteomics-informed operators are materially and consistently different from topology-only operators across real NCI-60 biological variation (57/57 cell lines, median relative spectral deviation 0.564).

### Claim 2 — topology-only panel design validity
Geometry-aware panel design outperforms random selection even on the topology-only operator, demonstrating that the panel design principle holds independently of proteomics availability.

### Claim 3 — proteomics-informed panel design validity
Geometry-aware selection derived from proteomics-informed operators outperforms random selection for preserving those operators under restriction, with the advantage strengthening as panel budget increases.

### Claim 4 — panel selection dominates operator type
Under random metabolite masking, proteomics-informed and topology-only operators show near-identical distortion (|Δ| < 0.006), confirming that the primary driver of geometry preservation is panel selection strategy, not operator type — which directly motivates the geometry-aware design criterion.

Together these claims mean:

> Real multi-omics proteomic heterogeneity deforms metabolic operator geometry in a structured way, and that structure can be exploited to design better metabolite panels than naive random selection — both with and without proteomics information.

---

## Recommended citation of the Layer 1 result in the manuscript

A compact paragraph for use in the manuscript:

> In the NCI-60 / Human1 benchmark, SWATH-derived proteomics supported partial modulation
> of 43.9% of GPR-linked reactions and yielded cell-line-specific operators across 57 aligned
> cell lines. These operators deviated substantially from the topology-only baseline (median
> relative spectral deviation 0.564, range 0.038–4.57). Geometry-aware panels derived from
> the topology-only operator (Experiment A) and from the proteomics-informed operators
> (Experiment A′) consistently outperformed random panels of equal size in preserving
> restricted operator geometry. For the proteomics-informed case, the mean distortion
> advantage increased monotonically from +0.152 at panel size 50 to +0.499 at panel size 500
> across all 57 cell lines. Under random masking, proteomics-informed and topology-only
> operators showed near-identical distortion (|mean Δ| < 0.006), confirming that panel
> selection strategy rather than operator type is the primary determinant of geometry
> preservation at fixed panel budgets.

---

# Remediation Plan — Figures 3–7 Integrity Issue

## Problem (confirmed by audit)

Manuscript Figures 3, 4, 5, 6, and 7 are generated from hard-coded arrays,
`np.random` placeholder values, or synthetic reaction weights — not from
computed results. No committed CSV backs any of them. This is a research-
integrity issue, not a presentation issue.

| Fig | Label | Script | Failure |
|-----|-------|--------|---------|
| 3 | `fig:identifiability` | `make_figures.py::figure3()` | `np.random.rand` + hand-typed curves |
| 4 | `fig:panel` | `make_figures.py::figure4()` | `np.random.rand(6,1)` + hand-typed curves |
| 5 | `fig:Biological_alignment` | `Biological_alignment.py` | 15 hand-typed arrays + `np.random.normal` null; "# REPLACE" |
| 6 | `fig:psi` | `PSI_preservation.py` | hand-typed `dpsi_*` + heatmap arrays |
| 7 | `fig:rewiring` | `Fig7_rewiring.py` | `synthetic_reaction_weights()` — not proteomics-derived |

## Actions already taken

1. **Quarantine script** (`quarantine_figures_3to6.sh`): moves the five scripts to
   `scripts/_quarantine_synthetic/` and the placeholder PDFs to
   `figures/deprecated_placeholders/`, each with a DO-NOT-USE README.
2. **Manuscript marked**: each dependent claim now carries
   "[UNSUPPORTED — pending Human1 run]" and an integrity banner heads the Results
   section. (Both are internal markers to remove once figures are real.)

No fabricated values were introduced and no synthetic figure was generated.

---

## Option A — Regenerate Figures 3–7 from real Human1-scale analyses (preferred)

Run the genuine model-scale experiments and emit raw + summary CSVs, then plot
from those CSVs only. This restores the figures as originally intended.

**Required runs (on Human-GEM, using the real operator code in
`proteomics_weighting.py`):**

| Fig | Experiment | Produces (results of record) | Plot from |
|-----|-----------|------------------------------|-----------|
| 3 | Masking / partial-observation sweep across panel sizes; geometry-aware vs random vs degree vs permuted | `human1_masking_raw.csv`, `human1_masking_summary.csv` | new `make_figure3.py` (CSV-only) |
| 4 | Greedy panel design; distortion vs panel size; selected-metabolite trace | `human1_panel_raw.csv`, `human1_panel_summary.csv` | new `make_figure4.py` (CSV-only) |
| 5 | Permutation robustness (OR-aggregator, saturation mappings, permutation null) | `human1_robustness_raw.csv`, `human1_robustness_summary.csv` | new `make_figure5.py` (CSV-only) |
| 6 | PSI preservation across selection strategies and subsystems | `human1_psi_raw.csv`, `human1_psi_summary.csv` | new `make_figure6.py` (CSV-only) |
| 7 | Cross-condition rewiring using REAL proteomics-derived weights (CPTAC), not synthetic | `human1_rewiring_real.csv`, `human1_rewiring_permuted.csv` | revised `Fig7_rewiring.py` (reads CSVs; real weights) |

**Contract for every replacement plotting script:**
- reads only committed CSVs,
- contains no `np.random`,
- contains no hard-coded result arrays,
- fails loudly if its CSV is absent.

**Effort:** the model-scale runs are the real cost (cobra + 43 MB Human-GEM,
minutes-to-hours each). The plotting layer is straightforward once CSVs exist;
the Figure 8/9 scripts are the template.

**Critical note for Figure 7:** the current script's `synthetic_reaction_weights()`
must be replaced with the real CPTAC-derived `W_R^(c)` (the same operator used
for the cohort and NCI-60 work). Synthetic weights cannot stand in.

---

## Option B — Replace Figures 3–7 with currently reproducible figures

If the Human1 runs cannot be completed in time, drop the unsupported figures and
rescope the manuscript around what is genuinely reproducible today:

- **Figure 9 (NCI-60 Layer 1)** — already reproducible from committed CSVs; covers
  panel design (A′), the random-vs-geometry-aware claim, and operator-type
  comparison (B). This subsumes much of what Figures 3–4 claimed, at model scale,
  on 57 real cell lines.
- **Figure 8 (ST003506 cohort)** — reproducible once the `--results_csv` patch and
  `make_figure8.py` are committed; covers the classification claim.
- **Figure S1 (toy model)** — fully reproducible; illustrates masking, PSI, and
  panel selection on a worked example with real computed values.

**Rescoping consequence:** the manuscript would no longer claim full-Human1
masking/PSI/rewiring results (Figs 3,4,5,6,7) and would instead present the
NCI-60 external validation + cohort + toy model as the evidence base. The
geometry-vs-flux limitation already conceded in the Discussion makes this a
coherent, honest narrowing. The two-layer claim survives because Layer 1
(NCI-60) and Layer 2 (cohort) are both real.

---

## Recommendation

**Option A if the Human1 runs are feasible before submission** — it preserves the
manuscript's scope. **Option B if not** — it is the honest fallback and still
yields a publishable, fully-reproducible paper built only on real figures.

Either way: do not submit with Figures 3–7 in their current placeholder form.

---

## Verification gate (before un-marking any figure)

A figure may have its "[UNSUPPORTED]" marker removed only when:
1. its results-of-record CSV is committed,
2. its plotting script reads only that CSV (no random, no hard-coded arrays),
3. the figure regenerates from a fresh clone, and
4. RESULTS_OF_RECORD.md and FIGURE_PROVENANCE.md list the CSV and script.

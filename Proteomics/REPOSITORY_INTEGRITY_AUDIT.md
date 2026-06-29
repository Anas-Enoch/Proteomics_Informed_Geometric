# Repository Integrity Audit

**Scope:** every Python script in the repository, classified by whether any
synthetic or hard-coded value reaches a *quantitative manuscript figure*.
Conceptual workflow diagrams are permitted to use illustrative values.
Quantitative manuscript figures must never depend on synthetic values.

**Method:** static trace (`grep`) of `np.random.*`, hard-coded result arrays,
and figure output paths, followed by manual classification of each occurrence
as (a) legitimate analysis randomness, (b) conceptual illustration, or
(c) fabrication reaching a quantitative figure.

**Verdict key:** PASS = clean / WARNING = non-fatal issue / FAIL = synthetic
values reach a quantitative manuscript figure.

---

## Summary table

| Script | random | hard-coded | Reaches quantitative figure? | Verdict |
|--------|:------:|:----------:|------------------------------|:------:|
| `nci60/nci60_script/layer1_build_S.py` | 0 | 0 | builds operators (no figure) | **PASS** |
| `nci60/nci60_script/layer1_spectral_scan.py` | 0 | 0 | writes spectral CSV (Fig 9) | **PASS** |
| `nci60/nci60_script/layer1_experiment_Aprime.py` | seeded | 0 | writes A′ CSV (Fig 9) | **PASS** |
| `nci60/nci60_script/layer1_masking_benchmark.py` | seeded | 0 | writes A/B CSVs (Fig 9) | **PASS** |
| `scripts/proteomics_weighting.py` | seeded | 0 | permuted-weight null operator | **PASS** |
| `scripts/real_cohort.py` | seeded | 0 | permutation test (Fig 8 data) | **PASS** |
| `scripts/Fig7_rewiring.py` | seeded | 0 | reads CSVs; seeded sampling in analysis | **WARNING** |
| `scripts/PSI_preservation.py` | 0 | 1 | Fig 6 — no CSV input | **WARNING** |
| `scripts/setup/HumanGEMInstaller.py` | 0 | 0 | setup utility | **PASS** |
| `scripts/patch_real_cohort.py` | 0 | 0 | patch helper | **PASS** |
| `scripts/Biological_alignment.py` | 1 | 15 | placeholder perm-null (line 90) | **FAIL** |
| `scripts/Fig8_three_operator_comparison.py` | 0 | 8 | **hard-coded AUROC arrays → Fig 8** | **FAIL** |
| `scripts/make_figures.py` | 9 | 4 | **np.random → figures 1–4** | **FAIL** |

---

## Legitimate randomness (NOT fabrication)

These uses of RNG are correct and must be kept. They are seeded analysis
operations, not invented figure values:

- **`real_cohort.py:286,314`** — `permutation_test_auc`: the permutation null
  is *supposed* to permute labels; this is the statistical test itself.
- **`proteomics_weighting.py:99-100`** — `rng.permutation(W_R)` builds the
  distribution-matched permuted-weight null operator. This is the control that
  breaks the PSI circularity argument; randomness is the point.
- **`layer1_experiment_Aprime.py`, `layer1_masking_benchmark.py`** — seeded
  random-panel baselines (the "random selection" arm of the real experiment).
- **`Fig7_rewiring.py:104,118`** — `np.random.default_rng(seed)` for seeded
  sampling inside the analysis. WARNING only because it computes inside the
  plotting script rather than reading a pre-computed CSV (see refactor note).

---

## FAIL details (must fix before submission)

### `scripts/make_figures.py` — FAIL (most serious)
Generates `figures/figure1.pdf`–`figure4.pdf` using `np.random.rand`,
`np.random.choice`, and `np.random.seed(0)` (lines 47, 56-57, 71, 177, 186-187,
211). These are the figures referenced in the manuscript. Every numerical mark
in figures 1–4 as committed is synthetic.
**Resolution applied in this work:** figures 1 and 2 were rebuilt as genuine
schematics (no fabricated data — they are conceptual by design and contain no
numerical claims). Figures 3 and 4 are quantitative and must be regenerated
from real Human1 masking / panel-design CSVs (computation scripts provided;
CSVs must be produced by running on Human-GEM). Until then they are FAIL.

### `scripts/Fig8_three_operator_comparison.py` — FAIL
Lines 44–54 hard-code the AUROC/accuracy/permutation arrays
(`auroc_mean = np.array([0.941, 0.923, 0.941])`, etc.). The values are *correct*
(they match `real_cohort.py` output) but they are typed into the plotting
script, not read from it.
**Resolution applied:** replaced by `make_figure8.py`, which reads only the
`cohort_{op}_summary/folds/permnull.csv` files emitted by the patched
`real_cohort.py`. The original script should be deleted or quarantined.

### `scripts/Biological_alignment.py` — FAIL / REPLACE
Line 90: `perm_null_dist = np.random.normal(loc=0.38, scale=0.03, size=200)`
with an inline comment `# REPLACE with your permuted runs`. This is an
acknowledged placeholder. If this script feeds any manuscript figure, the
placeholder null must be replaced with the real permuted-run distribution;
if it does not feed a manuscript figure, the script should be moved out of the
figure pipeline or deleted.

---

## Orphan / obsolete files

- **`nci60/layer1_masking_expA_panel_strategy.csv`** — present but **no script
  reads or writes it** (confirmed: `grep -rln masking_expA --include=*.py`
  returns nothing). It is superseded by `layer1_experiment_A_raw.csv` +
  `layer1_experiment_A_summary_rows.csv`. **Action: delete.** (See empty-CSV
  investigation note.)
- **`scripts/__pycache__/`** — build artifact; add to `.gitignore`.

---

## Bottom line

No FAIL remains in the *analysis* layer: every CSV that underpins Figure 9 and
the Figure 8 data is produced by seeded, real computation. All three FAILs are
in the *legacy plotting* layer (`make_figures.py`, `Fig8_three_operator_comparison.py`,
`Biological_alignment.py`). The fix is mechanical: delete/quarantine the legacy
plotting scripts and regenerate quantitative figures from committed CSVs. The
provided `make_figure8.py` and `make_figure9.py` already do this correctly.

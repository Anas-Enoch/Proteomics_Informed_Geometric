# nci60/

NCI-60 Layer-1 external validation (Figure 9). Geometry-aware panel design and operator
analysis across 57 NCI-60 cancer cell lines using SWATH proteomics (CellMiner).

Detailed methodology: [`README_NCI60_LAYER1.md`](README_NCI60_LAYER1.md).

## scripts/
| Script | Role |
|--------|------|
| `layer1_build_S.py` | builds aligned stoichiometry / weight bridge |
| `layer1_spectral_scan.py` | spectral deformation scan → `csv/layer1_spectral_scan_summary.csv` |
| `layer1_masking_benchmark.py` | Experiment A / B masking → `csv/layer1_experiment_{A,B}_*.csv` |
| `layer1_experiment_A.py`, `layer1_experiment_Aprime.py`, `layer1_experiment_B.py` | individual experiment drivers |

## csv/ (results of record for Figure 9)
| File | Backs |
|------|-------|
| `layer1_experiment_Aprime_{raw,summary}.csv` | Fig 9A/9B (proteomics-informed panel design, 57 lines) |
| `layer1_experiment_A_{raw,summary}.csv` | Fig 9 (topology-only panel design) |
| `layer1_experiment_B_{raw,summary}.csv` | Fig 9D (operator-type comparison, near-null) |
| `layer1_spectral_scan_summary.csv` | Fig 9C (bimodal spectral deformation) |
| `layer1_topology_baseline_eigs.csv` | baseline near-null eigenvalues (documents the K_EIGS correction) |

## Reproduce Figure 9
```bash
python ../human1_fig_scripts/make_figure9.py --csvdir csv --out ../figures/figure9.pdf
```
This reads the committed CSVs directly — no heavy recomputation needed. To regenerate the
CSVs themselves, run the `scripts/` pipeline (requires `cobra` + Human-GEM).

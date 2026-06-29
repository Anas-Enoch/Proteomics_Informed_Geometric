# human1_fig_scripts/

Model-scale (Human-GEM / Human1) analysis and plotting scripts for Figures 3-7.

## Two-stage contract
Each figure has an **analysis** script (writes CSVs to `../human1_fig_results/`) and a
**plotting** script (reads only those CSVs). Plotting scripts contain no `np.random`,
no hard-coded result arrays, and fail loudly if a required CSV is missing.

| Analysis script | Plotting script | Figure | Outputs (in `../human1_fig_results/`) |
|-----------------|-----------------|--------|----------------------------------------|
| `human1_fig3_geometry_preservation.py` | `make_figure3.py` | 3 | `..._fig3_geometry_preservation_{raw,summary}.csv` |
| `human1_fig4_panel_design.py` | `make_figure4.py` | 4 | `..._fig4_panel_design_{raw,summary}.csv` |
| `human1_fig5_robustness.py` | `make_figure5.py` | 5 | `..._fig5_robustness_{raw,summary}.csv`, `..._fig5_permutation_null.csv` |
| `human1_fig6_psi.py` | `make_figure6.py` | 6 | `..._fig6_psi_{raw,summary}.csv` |
| `human1_fig7_operator_deformation.py` | `make_figure7.py` | 7 | `..._fig7_{reaction_weights,spectral_shift,subspace_rotation,leverage_changes,summary}.csv` |

`make_figure9.py` (Figure 9) also lives here and reads `../nci60/csv/`.

## Shared core
`human1_analysis_core.py`: operator construction (via `proteomics_weighting.build_operator`),
spectral embedding, normalised distortion, diffusion-distance error,
greedy/degree/random panel selection, and the Calinski–Harabasz pathway separability
index (PSI; bounded, scale-invariant — it replaced a numerically degenerate ratio metric).

## Run everything
```bash
bash run_human1_figures.sh
```
Runs all five analyses and renders Figures 3-7. Requires `cobra` and the Human-GEM model.
The only randomness is fixed-seed metabolite sampling for random-panel baselines (seed 42)
and fixed-seed permuted-weight operators for the Figure 5 null (seeds 1000+) — documented
experimental randomness, not figure fabrication.

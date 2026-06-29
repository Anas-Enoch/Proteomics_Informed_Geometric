# human1_fig_results/

Results of record for the model-scale figures (3-7). These CSVs are the canonical
numerical outputs from which Figures 3-7 are generated; the plotting scripts in
`../human1_fig_scripts/` read only these files.

## Figure 3 — geometry preservation
- `human1_fig3_geometry_preservation_raw.csv` — per panel size / strategy / trial: spectral distortion, diffusion error
- `human1_fig3_geometry_preservation_summary.csv` — mean ± SD per panel size / strategy

## Figure 4 — panel design
- `human1_fig4_panel_design_raw.csv` — greedy trace (step, added metabolite, subsystem, distortion) + comparator panels
- `human1_fig4_panel_design_summary.csv` — distortion vs panel size by strategy

## Figure 5 — robustness
- `human1_fig5_robustness_raw.csv` — OR-aggregator and saturation variants (deterministic, n=1)
- `human1_fig5_robustness_summary.csv` — summary of the above
- `human1_fig5_permutation_null.csv` — 50 fixed-seed permuted-weight distortions + observed

## Figure 6 — pathway separability (Calinski–Harabasz PSI)
- `human1_fig6_psi_raw.csv` — per subsystem / panel size / strategy / trial
- `human1_fig6_psi_summary.csv` — global CH-PSI + normalised ΔPSI by panel size / strategy

## Figure 7 — CPTAC operator deformation
- `human1_fig7_reaction_weights.csv` — every reaction's W_R(c) + subsystem
- `human1_fig7_spectral_shift.csv` — topology vs proteomics eigenvalues, relative shift, cumulative distortion
- `human1_fig7_subspace_rotation.csv` — principal angles between low-frequency subspaces
- `human1_fig7_leverage_changes.csv` — per-metabolite Δ leverage + subsystem, ranked
- `human1_fig7_summary.csv` — coverage + W_R statistics + deformation scalars

To regenerate: run the matching analysis scripts in `../human1_fig_scripts/`
(see that folder's README).

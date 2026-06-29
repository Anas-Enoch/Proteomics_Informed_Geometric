#!/usr/bin/env bash
# ==============================================================================
# run_human1_figures.sh — generate real Figures 3-6 from Human1
#
# Runs the four analysis scripts (which write results/*.csv) then the four
# CSV-only plotting scripts. Requires cobra + the Human-GEM model + proteomics.
# ==============================================================================
set -euo pipefail

MODEL="Human-GEM-main/model/Human-GEM.xml"
PROT="data/proteomics/cptac_breast_tumor_only.tsv"
HGNC="data/hgnc_complete_set.txt"
OUT="results"
mkdir -p "$OUT" figures

echo "=== Figure 3: geometry preservation ==="
python3 human1_fig3_geometry_preservation.py --model "$MODEL" --proteomics "$PROT" --hgnc "$HGNC" --outdir "$OUT"
echo "=== Figure 4: panel design ==="
python3 human1_fig4_panel_design.py        --model "$MODEL" --proteomics "$PROT" --hgnc "$HGNC" --outdir "$OUT"
echo "=== Figure 5: robustness ==="
python3 human1_fig5_robustness.py          --model "$MODEL" --proteomics "$PROT" --hgnc "$HGNC" --outdir "$OUT"
echo "=== Figure 6: PSI preservation ==="
python3 human1_fig6_psi.py                 --model "$MODEL" --proteomics "$PROT" --hgnc "$HGNC" --outdir "$OUT"

echo "=== Plotting (CSV-only) ==="
python3 scripts/plotting/make_figure3.py --out figures/figure3.pdf
python3 scripts/plotting/make_figure4.py --out figures/figure4.pdf
python3 scripts/plotting/make_figure5.py --out figures/figure5.pdf
python3 scripts/plotting/make_figure6.py --out figures/figure6.pdf

echo "Done. Real Figures 3-6 generated from Human1 results of record."

# ══════════════════════════════════════════════════════════════════════════════
# PATCH for real_cohort.py
# Apply these three changes.  Nothing else in your script needs to change.
# ══════════════════════════════════════════════════════════════════════════════


# ─────────────────────────────────────────────
# CHANGE 1 — add import at top of real_cohort.py
# (after your existing imports)
# ─────────────────────────────────────────────

from proteomics_weighting import build_operator


# ─────────────────────────────────────────────
# CHANGE 2 — add two arguments to your argparse block
# (after your existing ap.add_argument calls)
# ─────────────────────────────────────────────

# OLD (nothing — these are new lines)
# NEW
ap.add_argument(
    "--operator",
    choices=["baseline", "proteomics", "permuted"],
    default="baseline",
    help=(
        "baseline  = S S^T  (W_R = I, topology only)\n"
        "proteomics = S W_R^(c) S^T  (CPTAC-informed)\n"
        "permuted   = S W_R_perm S^T  (distribution-matched null)"
    ),
)
ap.add_argument(
    "--proteomics",
    default=None,
    help="Path to CPTAC TSV file (required when --operator != baseline)",
)


# ─────────────────────────────────────────────
# CHANGE 3 — replace the single operator construction call
# ─────────────────────────────────────────────

# OLD
Delta_full = build_metabolite_laplacian(model)

# NEW  (drop-in replacement — same variable name, same type)
Delta_full = build_operator(
    model,
    mode=args.operator,
    proteomics_path=args.proteomics,
    condition="tumor",          # which CPTAC condition to use
    rho_0=0.1,                  # baseline coupling (Eq. 9)
    alpha=1.0,                  # proteomics dynamic range (Eq. 9)
    seed=args.seed,
)


# ══════════════════════════════════════════════════════════════════════════════
# That's it.  The rest of real_cohort.py is untouched.
#
# Run the three conditions:
#
#   python real_cohort.py \
#       --table  ST003506_AN005756.txt \
#       --model  Human-GEM.xml \
#       --operator baseline \
#       --out    results_baseline.pdf
#
#   python real_cohort.py \
#       --table  ST003506_AN005756.txt \
#       --model  Human-GEM.xml \
#       --operator proteomics \
#       --proteomics cptac_breast.tsv \
#       --out    results_proteomics.pdf
#
#   python real_cohort.py \
#       --table  ST003506_AN005756.txt \
#       --model  Human-GEM.xml \
#       --operator permuted \
#       --proteomics cptac_breast.tsv \
#       --out    results_permuted.pdf
#
# Collect the three mean AUC values and report them in a single table.
# ══════════════════════════════════════════════════════════════════════════════

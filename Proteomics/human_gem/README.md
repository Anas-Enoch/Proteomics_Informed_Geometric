# human_gem/

The Human-GEM genome-scale metabolic model and accessory tables.

| File | Description |
|------|-------------|
| `model/Human-GEM.xml` | SBML model (~43 MB) loaded by `cobra`. Canonical input for all model-scale operator construction. |
| `Human-GEM.txt`, `Human-GEM.xlsx` | Human-readable reaction/metabolite tables (reference only). |
| `genes.tsv` | Gene annotation table (reference only). |

The pipeline reads `model/Human-GEM.xml`; stoichiometry $S$, GPR rules, and subsystem
annotations are extracted via `cobra`. Version and licensing: [`DATA_SOURCES.md`](../DATA_SOURCES.md).

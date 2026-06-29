# real_cohort/

ST003506 breast-cancer cohort classification (Figure 8). Three-operator comparison
(topology-only / CPTAC proteomics-informed / permuted-weight null) using spectral
embeddings of the restricted metabolite operator.

| File | Description |
|------|-------------|
| `real_cohort.py` | analysis: parses ST003506, builds the operator, runs repeated stratified CV + a permutation test, writes results CSVs via `--results_csv` |
| `make_figure8.py` | plotting: reads only the CSVs below, no hard-coded values |
| `figure8_summary.csv` | per-operator mean ± SD AUROC/accuracy + permutation-test stats |
| `figure8_cv_results.csv` | per-fold AUROC/accuracy (the distribution behind the plot) |
| `figure8_permutation_null.csv` | full permutation null distribution per operator |

## Reproduce Figure 8
```bash
for op in baseline proteomics permuted; do
  python real_cohort.py --table ../data/ST003506_AN005756.txt \
    --model ../human_gem/model/Human-GEM.xml --operator $op \
    $( [ "$op" != baseline ] && echo "--proteomics ../data/cptac_breast_tumor_only.tsv --hgnc ../data/hgnc_complete_set.txt" ) \
    --results_csv figure8.csv
done
python make_figure8.py --stem figure8 --out ../figures/figure8.pdf
```
The reported AUROC values are computed by `real_cohort.py` (real cross-validation), not
hard-coded; `make_figure8.py` only renders the persisted CSVs.

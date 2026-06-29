# toy_example/

A small, fully-specified worked example (Figure S1) that exercises the entire pipeline on
a synthetic-but-explicit toy metabolic system. Unlike a black-box test, every matrix is
written out, so the operator construction, spectral embedding, panel selection, and PSI
can be checked by hand.

| File | Description |
|------|-------------|
| `toy_pipeline.py` | runs the toy system end-to-end and writes the CSVs below + `../figures/figureS1_toy_example.pdf` |
| `toy_S.csv` | toy stoichiometric matrix |
| `toy_WR.csv` | toy reaction weights |
| `toy_laplacian.csv` | resulting metabolite operator $\Delta_M$ |
| `toy_eigenpairs.csv` | low-frequency eigenvalues / eigenvectors |
| `toy_panel_trace.csv` | greedy panel-selection trace |
| `toy_psi.csv` | per-subsystem PSI on the toy system |

This is the only place a synthetic system appears, and it is **explicitly** a worked
example: the toy model is part of the illustration, not a stand-in for real results. It is
deterministic (fixed seed). Run:
```bash
python toy_pipeline.py
```

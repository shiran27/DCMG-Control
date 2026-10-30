# Regenerated MATLAB results
Scope: Robust CT integral-data local design; known-B global composition
Run `report = run_all_paper_results` for the complete comparison and checks; `run_paper_results` rebuilds this core record alone.

| Method | Directed links | Certified gamma | Settling (ms) | Max real pole |
|---|---:|---:|---:|---:|
| MB stabilizing | 30 | NaN | 5.2 | -408.93 |
| MB co-design | 0 | 5.0872 | 3.3 | -315.55 |
| DD local + known-B global | 0 | 5.0975 | 3.2 | -331.97 |
| Conventional droop (5%) | 0 | NaN | 6.3 | -213.63 |

The general unknown-B global theorem is not validated by this implementation. Integral-data local synthesis is a sufficient continuous-time counterpart. Model-based and data-driven co-designs share the same plant, scaling, tuning grid, budget, and certificate gates.

See `paper_results.mat` for all raw data, unrounded gains and certificates; `metrics.csv` and `droop_sensitivity.csv` for metrics; `design_attempts.mat` for accepted and rejected attempts; `numerical_results.tex` for generated manuscript text; and `.fig`/`.pdf`/`.png` for figures.

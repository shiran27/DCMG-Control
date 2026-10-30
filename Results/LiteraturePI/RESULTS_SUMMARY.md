# Literature PI baseline comparison
Source: Nahata et al., Automatica 113 (2020), 108770. https://doi.org/10.1016/j.automatica.2019.108770
Run `run_all_paper_results` or `main.mlx` to rebuild all controllers, results, and regression checks. To reuse the verified core gains while regenerating the PI design and every comparison response, run `run_literature_pi_comparison`.

| Method | Settling ms | Pulse voltage RMS V | Pulse sharing RMS pu | Pulse command RMS V | Load-step tail voltage RMS V | Load-step tail sharing RMS pu |
|---|---:|---:|---:|---:|---:|---:|
| MB stabilizing | 5.20 | 0.07497 | 0.01763 | 0.07127 | 0.201046 | 0.01914 |
| MB co-design | 3.30 | 0.15998 | 0.00339 | 0.26108 | 1.13791 | 0.02445 |
| DD local + known-B global | 3.20 | 0.15204 | 0.00336 | 0.25367 | 1.08696 | 0.02458 |
| Conventional droop (5%) | 6.30 | 0.12301 | 0.01663 | 0.04163 | 0.0982666 | 0.01184 |
| Passivity PI (Nahata et al.) | 7.40 | 0.04838 | 0.00820 | 0.14630 | 6.80253e-06 | 0.01994 |

Selected kp=2, kd=0.5, integral fraction=0.25; 60 of 125 candidates pass. Nominal peak-command budget: 4.08865 V.
Independent ODE residual: 2.27e-11. Absolute-law residual: 2.84e-14. Local/network passivity checks, constant-load equilibrium check, matched input check, and peak budget check: passed.
The constant-load test holds the pre-change references fixed. The PI controller regulates each voltage to that reference; proportional sharing under a changed load would need a separate secondary layer. No L2 bound is claimed for PI.

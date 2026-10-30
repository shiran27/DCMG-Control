# MATLAB literature comparison completion — 2026-09-11

The duplicate task edits have been reconciled into one reproducible experiment.
Both literature baselines are retained: conventional decentralized V–I droop
and the passivity-based multivariable PI controller of Nahata et al. The three
existing methods are freshly synthesized MB stabilizing, MB co-design, and
integral-data local / known-B global co-design. All five use the existing
six-DG physical model and the same verified operating point.

## One complete reproduction command

```matlab
report = run_all_paper_results; % main.mlx calls this same function
sync_paper_results('/path/to/Data_Driven_Design_ACC2025');
```

The first command rebuilds every controller, generates the data and all
simulation responses, exports figures/tables/discussion, and runs the complete
regression suite. It uses the existing MATLAB classes and simulation setup.
The synchronization command verifies that the core and PI records share the
same source generation and trajectories before copying both sets of publication
assets. Repeated synchronization preserves exactly one PI manuscript include.

For focused work, `run_paper_results` rebuilds the four core methods, and
`run_literature_pi_comparison` subsequently tunes the PI controller and
recomputes all five methods' comparison responses from the verified core gains.
Run both before synchronizing if the core experiment changes.

## Results and interpretation

| Controller | Nominal settling (ms) | Disturbance voltage RMS (V) | Disturbance sharing RMS (p.u.) | Disturbance command RMS (V) |
|---|---:|---:|---:|---:|
| MB stabilizing | 5.2 | 0.07497 | 0.01763 | 0.07127 |
| MB co-design | 3.3 | 0.15998 | 0.00339 | 0.26108 |
| DD local / known-B global | 3.2 | 0.15204 | 0.00336 | 0.25367 |
| Conventional droop (5%) | 6.3 | 0.12301 | 0.01663 | 0.04163 |
| Passivity PI | 7.4 | 0.04838 | 0.00820 | 0.14630 |

The co-designs have faster nominal settling and smaller transient sharing error
than either literature baseline. PI gives the smallest disturbance voltage
error. Droop uses the least disturbance command effort. In the separate
permanent-load test, PI restores the fixed voltage references, with a
100–120 ms voltage RMS of about 6.8e-6 V; this does not guarantee proportional
current sharing after the load changes. The manuscript reports these tradeoffs
without claiming uniform superiority.

Droop uses a declared 5% rated-current slope, with a saved 1/3/5/10% sensitivity
sweep; it implements the conventional primary law cited in Guerrero et al.
(2011) and Lu et al. (2014), without claiming their secondary restoration loops.
PI follows Nahata et al., *Automatica* 113 (2020), 108770,
DOI 10.1016/j.automatica.2019.108770. Its analytical gain region and storage
are checked on the static-resistor ZI specialization. The nominal-only grid
selects kp=2, kd=0.5, and integral fraction=0.25; 60 of 125 candidates pass
the declared stability, gain, convergence, and command-effort gates.
See `LITERATURE_PI_BASELINE.md` for the exact source equations and mapping.

## Executed validation

- The complete `main.mlx` entry point executed successfully in MATLAB R2025a
  with YALMIP and SeDuMi 1.3.8.
- Original MB/DD certificates, saved data and gains, rollback after a rejected
  design, graph construction, and droop sensitivity checks passed.
- PI local/network passivity and absolute-law checks passed. Its independent
  ODE residual is 2.27e-11; the absolute/equilibrium-centered law residual is
  2.84e-14.
- All shared trajectories and commands agree between the core record and the
  PI comparison. Exported CSV values agree with the saved MAT record. A
  deliberately stale comparison is rejected.
- The retained three-DG and archived six-DG regressions completed without
  failures. The expected historical 15-link versus current 0-link topology
  difference remains a labeled warning.
- Synchronization is idempotent; both manuscript inputs and all three baseline
  bibliography entries appear once. All 19 synchronized publication files are
  byte-identical to the MATLAB outputs.
- The updated paper compiles to 21 pages without undefined or multiply defined
  references. The numerical-results and comparison pages were rendered and
  visually checked for complete tables, legible figures, and unclipped curves.

`Results/Regenerated` contains the core record, four-controller figures,
droop sensitivity, and generated numerical section. `Results/LiteraturePI`
contains all five methods' matched metrics and trajectories, PI tuning,
publication comparison, `suite_checks.mat`, and `full_suite.log`.

The paper repository includes the generated sections under these same two
result-directory names. Its canonical compiled delivery is
`Revision/Paper_with_Regenerated_Results.pdf`. Uploaded files under
`output/source` and older root-level `Results` figures remain historical inputs,
clearly labeled as such.

## Remaining paper scope

This completes the literature-comparison implementation and results task.
The DD experiment remains explicitly local integral-data synthesis followed by
known-input-matrix global composition. The previously documented obstruction
in the general unknown-B global theorem is separate from this comparison and
remains unresolved; these results do not claim to validate that construction.

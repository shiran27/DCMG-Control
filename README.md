# DCMG control and reproducible paper results

For author review, open [the paper checklist](Review/01_PAPER_REVIEW.md) and
[the implementation checklist](Review/02_IMPLEMENTATION_REVIEW.md) side by side.
They share seven decision IDs, source locations, and alternatives.

Run in MATLAB from this repository:

```matlab
report = run_all_paper_results;
```

`main.mlx` calls that same complete entry point. It regenerates the proposed
controllers, conventional droop, the published PI baseline, all figures and
tables, and the regression checks. The implementation extends the existing
`DG`, `TransmissionLine`, `DCMicrogrid`, six-DG factory, and simulation setup.
It does not load old controller gains or old trajectories to produce new results.

See [the comparison completion record](COMPARISON_COMPLETION_20260911.md) for
the reconciled baseline additions, measured tradeoffs, and executed checks.

## Validated scope

- Continuous-time six-DG, seven-resistor network with ZI loads.
- Conventional decentralized V--I droop literature baseline with a fixed 5%
  rated-current slope and a declared 1--10% sensitivity sweep.
- Published decentralized passivity-based PI comparison (Nahata et al.,
  Automatica 113, 2020), with its explicit gain region and storage function
  checked on the same physical benchmark.
- Fresh dense model-based stabilizing controller.
- Hard local and global model-based dissipativity co-design.
- Robust local synthesis from exact continuous-time integral data, followed by
  the **known-input-matrix global** composition LMI.
- Original local, supply-dominance, and aggregate storage inequalities are
  independently checked. Failed solves do not install gains. Graph thresholding
  proposes structural zeros; a second hard solve must certify the graph.

The general unknown-B global lift in the original paper is not a validated
implementation: its final principal block forces a zero transformed gain under
strict QMI interior. The new DD experiment is a sufficient continuous-time
integral-data counterpart, not a claim to have repaired that theorem. See the
paper repository's `Revision/IMPLEMENTATION_COMPLETION_20260910.md`.

## Reproduction and dependencies

MATLAB R2025a and YALMIP were used. `support/dcmgSolverOptions` probes MOSEK,
SeDuMi and SDPT3. The installed MOSEK license was expired; this run used SeDuMi
1.3.8. See `vendor/README.md` to reproduce the solver installation. No saved
MATLAB path or license configuration is changed.

The core controllers and droop artifacts are in `Results/Regenerated/`:

- `paper_results.mat`: full verified run, parameters, seeds, data, source hashes,
  controller gains, decision matrices, certificates, trajectories and metrics.
- `design_attempts.mat`: every nominal tuning attempt, including failures.
- `metrics.csv`, `droop_sensitivity.csv`, `RESULTS_SUMMARY.md`,
  `numerical_results.tex`: values and prose generated directly from the current
  result record.
- `paper_*.pdf`, `.png`, `.fig`: publication and editable MATLAB figures.
- `run.log`: complete run log; `regression_checks.mat`: serialization regression.

The workflow includes the archived large initial state, a moderate nominal
perturbation, matched primitive disturbances/state jumps, bounded-noise responses,
and an empirical physical-parameter variation study with offline MB/DD redesign.
The conventional droop law is `V_ti = V_0i - R_di I_ti`, with
`R_di = 0.05 V_rated,i/I_rated,i`; its no-load intercept is shifted to preserve
the same verified equilibrium used by every comparison controller. It is labeled
as an empirical, non-certified primary-control baseline.

`run_literature_pi_comparison` writes a separate verified comparison to
`Results/LiteraturePI/`. Its 125-candidate gain grid is selected using nominal
settling under the proposed controllers' command-effort budget; pulse, noise,
and permanent-load scenarios are withheld from tuning. The saved checks cover
the published local storage, the coupled-network dissipation matrix, the
absolute/equilibrium-centered control-law equivalence, closed-loop poles, and
an independent ODE trajectory.
See `LITERATURE_PI_BASELINE.md` for the source equations and coordinate mapping.
Clipping is disabled for certified linear runs. Exceedances of historical limits
are still logged. Set `DG.EnableLimiters=true` only for separately labeled
nonlinear experiments, which are outside these certificates.

To replot without optimization:

```matlab
addpath('support')
load('Results/Regenerated/paper_results.mat','report')
plotPaperResults(report,'Results/Regenerated')
```

To copy publication assets to the existing manuscript:

```matlab
sync_paper_results('/path/to/Data_Driven_Design_ACC2025')
```

This verifies that both result records come from the same core run, copies both
asset folders, and ensures exactly one PI comparison input in the manuscript.
It regenerates publication figures, text, and CSVs from those records before
copying, including after `MakePlots=false` runs, and rejects stale combinations.
Raw MAT records remain in the MATLAB repository.

## Data contract

`collectIntegralData` simulates the existing continuously coupled plant with held
converter commands and logs interval-average state/input measurements. The exact
regression is `Delta=A*X+B*U+W`, even while physical line currents vary within an
interval. Measurement noise is zero. A declared input-matrix norm bound and a
bounded primitive disturbance define the energy bound before data generation.
True model matrices are used only by the generator and independent verifier;
the DD solver consumes only the measured regression/QMI. A test deliberately
changes the withheld A/B matrices and checks that the DD gain does not change.

Legacy forward-Euler/fitted-residual data ingestion and the old seven-stage DD
redesign path are rejected with explicit contract errors. Use the new reproducible
entry point for certified results. The old Live Script is preserved in
`archive/legacy/main_before_results_20260910.mlx` for provenance.

## Verification

```matlab
addpath('tests')
checks = test_regenerated_results;
droopChecks = test_conventional_droop_baseline;
% Also exercise the retained three-DG and archived six-DG benchmarks:
regressions = run_all_regressions;
```

The main run also cross-checks exact propagation against `DCMicrogrid.dynamics`
with `ode45`, checks command-only ZOH against the existing stage simulator,
verifies the equilibrium/rank/bound/storage conditions, and retains every run's
settings. The older three-DG simulator and archived six-DG configuration tests
remain available; a topology difference from archived gains is expected.

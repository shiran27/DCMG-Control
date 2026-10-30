# Results Directory

Current publication results are in `Regenerated/` and `LiteraturePI/`.
Run `run_all_paper_results` or `main.mlx` to recreate both folders and run all
regression checks. The root-level files listed below are historical artifacts.

- `TrajResults.mat` stores per-method trajectories used by the archived Live Script.
- `.fig` files preserve editable MATLAB figures.
- `.png` files are rendered figure exports.
- `ModelBasedSimulatorBaseline.mat` stores the structured report from the independent deterministic three-DG model-based test.
- `ModelBasedSimulatorBaseline.fig` and `.png` compare open-loop, model-based stabilizing, LQR, and model-based dissipative responses from that test.
- `ModelBasedSimulatorNetworkStates.fig` and `.png` show the example DCMG at its initial state and the terminal state reached by each model-based controller, including physical-line currents and communication links.
- `PaperConfigurationModelBasedComparison.mat` stores the recovered six-DG benchmark report, controller gains, poles, trajectories, and topology checks.
- `PaperConfigurationModelBasedComparison.fig` and `.png` compare current model-based stabilizing, archived paper dissipative, and current dissipative responses.
- `PaperConfigurationNetworkStates.fig` and `.png` show the recovered initial network and each controller's terminal state/communication topology.

These are generated artifacts, but they remain relevant to reproducing and displaying the paper results.

## Regenerated results (2026-09-10)

`Regenerated/` contains the newly synthesized, independently verified experiment.
It supersedes historical paper plots for the explicitly documented integral-DD
local / known-B global scope and includes the conventional 5% V--I droop
literature baseline plus its declared coefficient-sensitivity table. Run
`run_paper_results` to recreate this core record alone.
The older artifacts above remain as provenance and are not inputs to the new run.

`LiteraturePI/` contains the separately reproduced Nahata et al. decentralized
passivity-based PI comparison. Its MAT record includes all 125 tuning candidates,
the selected gains and storage checks, matched trajectories, persistent-load
diagnostics, source hashes, and independent ODE/control-law residuals.

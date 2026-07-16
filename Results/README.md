# Results Directory

This directory contains generated paper figures and the accumulated trajectory dataset.

- `TrajResults.mat` stores per-method trajectories used by the final comparison section of `main.mlx`.
- `.fig` files preserve editable MATLAB figures.
- `.png` files are rendered figure exports.
- `ModelBasedSimulatorBaseline.mat` stores the structured report from the independent deterministic three-DG model-based test.
- `ModelBasedSimulatorBaseline.fig` and `.png` compare open-loop, model-based stabilizing, LQR, and model-based dissipative responses from that test.
- `ModelBasedSimulatorNetworkStates.fig` and `.png` show the example DCMG at its initial state and the terminal state reached by each model-based controller, including physical-line currents and communication links.
- `PaperConfigurationModelBasedComparison.mat` stores the recovered six-DG benchmark report, controller gains, poles, trajectories, and topology checks.
- `PaperConfigurationModelBasedComparison.fig` and `.png` compare current model-based stabilizing, archived paper dissipative, and current dissipative responses.
- `PaperConfigurationNetworkStates.fig` and `.png` show the recovered initial network and each controller's terminal state/communication topology.

These are generated artifacts, but they remain relevant to reproducing and displaying the paper results.

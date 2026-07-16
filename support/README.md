# Main Workflow Support Functions

`main.mlx` adds this directory to the MATLAB path when it starts. These
functions were moved out of the Live Script so its experiment stages remain
readable and the utilities can be tested independently.

| File | Role | Inputs | Output |
|---|---|---|---|
| `getARandomDCMicrogrid.m` | Builds the evolving random DCMG used by `main.mlx`. | Number of DGs. | Configured `DCMicrogrid`. |
| `createPaperDCMicrogrid.m` | Reconstructs the fixed six-DG, seven-line model-based paper benchmark. | Optional final time, noise step/factor, and initial-state flag. | `DCMicrogrid` and benchmark metadata/reference gains. |
| `drawStateTrajectories.m` | Converts stacked trajectories to per-unit values and draws state/input plots. | Time, states, inputs, operating points, DG count, save flag, title. | Per-unit trajectory structure. |
| `saveFigureHighRes.m` | Saves the current or selected MATLAB figure as PNG and FIG under `Results/`. | Base name and formatting options. | None. |
| `adjustAxesPadding.m` | Expands axes limits using caller-provided padding. | Axes and four padding values. | None. |
| `loadOrInitResults.m` | Loads a stored result structure or initializes an empty one. | MAT filename and reset flag. | Result structure array. |
| `appendResult.m` | Adds or replaces one named simulation result. | Results, method name, time, per-unit data. | Updated result structure array. |
| `plotComparisonAcrossMethods.m` | Plots average per-unit states and inputs across stored methods. | Results and four-row axis-limit matrix. | None. |

The paper benchmark constants come from Git revision `6a0091d` and were
cross-checked against `archive/data/matlab.mat`, the accepted six-DG
dissipative-run workspace.

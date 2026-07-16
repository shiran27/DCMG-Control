# DCMG Control Simulation

This folder contains the MATLAB implementation used to generate the controller comparisons and figures for the paper. `main.mlx` is the canonical entry point.

## Key runtime files

| File | Role |
| --- | --- |
| `main.mlx` | Sole canonical workflow. Orchestrates network generation/loading, controller design, staged simulation, result storage, and comparison plots. Reusable functions are in `support/`. |
| `support/` | Runtime helpers used by `main.mlx`, plus the fixed six-DG paper-benchmark factory. `main.mlx` adds this directory to the MATLAB path. |
| `DCMicrogrid.m` | Network-level model, steady-state calculation, staged simulation, data assembly, controller co-design, and topology drawing. |
| `DG.m` | Per-DG physical model, runtime dynamics, and local model-based/data-driven dissipativity designs. |
| `TransmissionLine.m` | Resistive line model and line-current convention used to assemble the network Laplacian. |
| `tests/test_model_based_simulator.m` | Deterministic three-DG baseline test. Rebuilds the plant for open-loop, model-based stabilizing, LQR, and model-based dissipative cases; checks equilibrium, poles, convergence, and limiter activation. |
| `tests/test_paper_model_based_configuration.m` | Regression test for the recovered six-DG, seven-line paper setup. Compares current stabilizing, archived paper dissipative, and current dissipative controllers without overwriting paper figures. |
| `netFile.mat` | Generated network snapshot loaded by the independent sections of `main.mlx`. Regenerate it by running the first live-script section. |
| `Results/TrajResults.mat` | Accumulated trajectories loaded by the final cross-method comparison section. |
| `Results/*.fig`, `Results/*.png` | Figure artifacts used to inspect or include simulation results in the paper. |
| `output/source/P6_Data_Driven_Co_Design.pdf` | Uploaded 18-page paper used as the primary theory reference. |
| `output/source/P6_Data_Driven_Co_Design_with_notes.pdf` | Distinct older 16-page annotated paper retained as background. |
| `output/source/P6.tex` | Uploaded paper TeX source. |
| `output/review/IMPLEMENTATION_REVIEW_CHECKLIST.md` | Editable implementation-review checklist, prioritized by criticality and revision effort. |
| `output/review/APPLICATION_SPECIFIC_THEORY_CODE_CROSS_REFERENCE.md` | DCMG-specific equations, dimensions, controller transformations, data contracts, and function acceptance checks. |
| `output/review/*.pdf` | Print-ready PDF counterparts colocated with both review Markdown documents. |
| `output/` | Review documents under `review/` and uploaded paper/background inputs under `source/`; no MATLAB runtime files are duplicated there. |

## Function reference

The class tables omit the implicit `obj` input. All three classes are handle
classes, so methods can update the existing object even when they have no
explicit output. `X` uses the stacked order `[V1; It1; ...; VN; ItN]`; `U`
uses `[Inet1; Vt1; ...; InetN; VtN]`.

### `main.mlx`

| Function/section | Role | Inputs | Outputs |
| --- | --- | --- | --- |
| Top-level Live Script sections | Generate/load a DCMG, design controllers, simulate, save results, and plot comparisons. | Workspace settings and generated MAT files. | `netFile.mat`, `Results/TrajResults.mat`, figures, and workspace variables. |

### `support/`

| Function | Role | Inputs | Outputs |
| --- | --- | --- | --- |
| `getARandomDCMicrogrid` | Build the evolving random connected DCMG used by `main.mlx`. | `N`: number of DGs. | Configured `DCMicrogrid`. |
| `createPaperDCMicrogrid` | Reconstruct the fixed six-DG, seven-line paper benchmark and archived reference design. | Name-value horizon, noise, and initial-state options. | `net`, benchmark metadata/reference gains. |
| `drawStateTrajectories` | Plot states/inputs and convert trajectories to per unit. | `t`, `X`, `U`, `Xs`, `Us`, `N`, save flag, title. | Per-unit trajectory struct; optional figures. |
| `saveFigureHighRes` | Save a selected figure as PNG and FIG. | Base name and formatting options. | Files under `Results/`. |
| `adjustAxesPadding` | Expand axes limits around plotted objects. | Axes and four padding values. | Updated axes. |
| `loadOrInitResults` | Load an existing result array or create an empty one. | Results filename and reset flag. | Result struct array. |
| `appendResult` | Add/replace one method and compute mean/std traces. | Results, method, time, per-unit trajectories. | Updated result struct array. |
| `plotComparisonAcrossMethods` | Compare mean per-unit states and inputs across methods. | Results and four-row axis limits. | Comparison figures. |

### `DCMicrogrid.m`

| Method | Role | Inputs | Outputs |
| --- | --- | --- | --- |
| `DCMicrogrid` | Construct the network and its conductance model. | DG array, transmission-line array. | `DCMicrogrid` handle. |
| `buildConductance` | Assemble `Y` and Laplacian `YBar` from line conductances. | None. | Updates network/DG properties. |
| `getStateVector` | Stack all local DG states. | None. | `X`: `2N x 1` state vector. |
| `setStateVector` | Distribute a stacked state to the DG objects. | `X`: `2N x 1`. | Updates DG states. |
| `simulate` | Run the seven-stage event simulation and optional DD redesigns. | `tspan`, `x0`, `useData`, gain-pruning threshold. | `t`, `X`, `U`, `Xs`, `Us`. |
| `simulateStageZOH` | Simulate one interval with held local/global feedback. | Start/end time, `x0`, output step `dt`. | Segment `t`, `X`, `U`, `Xs`, `Us`. |
| `loadDataMatrices` | Sample trajectories and build local/global DD QMI data. | `t`, `X`, `Xs`, `Utilde`, `Wtilde`. | Data stored on DG/network objects; declared `out` is currently unassigned. |
| `compute_Qw_from_wtilde` | Fit a quadratic disturbance description by SDP. | `W`: disturbance/residual samples. | `Qw_val`: disturbance matrix. |
| `computeUWTrajectories` | Reconstruct total, steady, deviation, and disturbance inputs. | `t`, `X`, `useData`. | `U`, component array `Uc`, `Utilde`, `Wtilde`, `Us`, `Xs`. |
| `dynamics` | Evaluate coupled continuous-time DCMG dynamics. | Time `t`, stacked state `X`, held-data flag. | `dX`: state derivative. |
| `setupNoise` | Generate per-DG sampled disturbance sequences. | `tspan`, noise step, `2 x 2` noise factor. | Updated object/noise properties. |
| `draw` | Draw physical topology, state labels, optional line currents, and communication links. | Axes plus optional title/display settings. | Graphics in the supplied axes. |
| `drawComm` | Draw directed communication links implied by nonzero gain blocks. | Axes and optional arguments. | Graphics in the supplied axes. |
| `buildSystemMatrices` | Assemble aggregate `A`, `B`, `BBar`, `E`, `D`, and `DBar`. | None. | Updates network matrices and `wBar`. |
| `solveSteadyState` | Solve equilibrium voltage/current sharing and commands. | Optional `opts` argument (currently unused). | `ss`; also updates `x_s`, `u_s`, and DG operating points. |
| `design_MB_GSC` | Design the model-based global stabilizing controller. | Pruning threshold and data-collection flag. | `AdjMat`, pruned `KMat`, solver `out`; updates `obj.K`. |
| `buildCommAdjFromK` | Threshold `1 x 2` gain blocks into a communication graph. | Relative threshold `thr`. | `Adj`, pruned `K`; updates network gain/adjacency. |
| `codesign_MB_DRC` | Co-design local and global model-based dissipative controllers. | Pruning threshold. | `AdjMat`, `KMat`, solver `out`; updates local/global gains. |
| `codesign_DD_DRC` | Co-design local and global data-driven dissipative controllers. | Pruning threshold. | `AdjMat`, `KMat`, solver `out`; updates local/global gains. |
| `criticalLeadingMinor` | Compute determinant-based leading-minor diagnostics. | Matrix `M`. | First critical index pair, minimum minor, all minors. |

### `DG.m`

| Method | Role | Inputs | Outputs |
| --- | --- | --- | --- |
| `DG` | Construct one converter/load subsystem. | DG `id`, parameter struct. | `DG` handle. |
| `updateModel` | Rebuild local `A`, `B`, `E`, and `BBar` from physical parameters. | None. | Updates model matrices. |
| `setState` | Set local `[V; It]`. | Two-state vector. | Updates `x`. |
| `getState` | Read local `[V; It]`. | None. | Two-state vector `x`. |
| `dynamics` | Evaluate local converter/load dynamics with limiting and noise. | `t`, network current, global command, held-data flag. | `dx`: two-state derivative. |
| `draw` | Draw a DG marker and compact electrical-state label. | Axes plus optional current/label/style settings. | Graphics in the supplied axes. |
| `designLocalXiDissipative` | Solve the model-based local Xi-dissipativity LMIs. | Stored physical model. | Solver/result struct `out`; updates `K`, `nu`, and `rho`. |
| `designLocalXiDissipative_DataDriven` | Solve the data-driven local Xi-dissipativity QMIs/LMIs. | Stored sampled data and disturbance matrices. | Solver/result struct `out`; updates `K`, `nu`, and `rho`. |
| `criticalLeadingMinor` | Compute local leading-minor diagnostics. | Matrix `M`. | Critical index pair, minimum minor, all minors. |

### `TransmissionLine.m`

| Method | Role | Inputs | Outputs |
| --- | --- | --- | --- |
| `TransmissionLine` | Construct one resistive physical link. | Line `id`, endpoint indices `i/j`, resistance `R`. | `TransmissionLine` handle. |
| `current` | Evaluate oriented line current from endpoint `i` to `j`. | Endpoint voltages `vCi`, `vCj`. | `i_ij = (vCi-vCj)/R`. |
| `draw` | Draw the line and optional resistance/current label. | Axes, endpoint positions, optional display settings. | Graphics in the supplied axes. |

### `tests/test_model_based_simulator.m`

Only `test_model_based_simulator` is the public entry point; the remaining
functions are local helpers in the same file.

| Function | Role | Inputs | Outputs |
| --- | --- | --- | --- |
| `test_model_based_simulator` | Run open-loop, MB stabilizing, LQR, and MB dissipative baselines. | Name-value options for plots, artifacts, failure handling, horizon, and step. | Structured `report` and optional result artifacts. |
| `makeExampleDCMG` | Build the deterministic heterogeneous three-DG ring. | Final time and output step. | `net`, initial state `x0`, state normalization scale. |
| `configureController` | Install the selected model-based controller. | Network, method key, state scale. | Design status struct and total command gain. |
| `localCommandGain` | Embed local DG gains in the aggregate gain matrix. | Network. | Block-local gain matrix. |
| `evaluateTrajectory` | Compute convergence, voltage/current, command, and limiter metrics. | Network, `t`, `X`, command gain, state scale. | Metrics struct. |
| `classifyResult` | Assign PASS/WARN/FAIL from residual, poles, and metrics. | Equilibrium residual, poles, metrics. | Status string. |
| `buildSummary` | Convert case results into the displayed table. | Case struct array. | Summary table. |
| `plotComparison` | Plot response errors, commands, and pole comparison. | Case struct array. | Response-comparison figure. |
| `plotNetworkStates` | Plot the initial DCMG and each controller's terminal network. | Cases and retained network handles. | Network-state figure. |
| `makeTimeGrid` | Build a grid that includes the exact final time. | Final time and output step. | Time vector. |
| `emptyCaseResult` | Create a fixed-schema placeholder for one case. | None. | Empty case struct. |

### `tests/test_paper_model_based_configuration.m`

Only `test_paper_model_based_configuration` is public. Its local helpers
validate the frozen benchmark, configure each controller, simulate the linear
paper model, calculate topology/gain diagnostics, and create the comparison
and terminal-network figures.

| Function | Role | Inputs | Outputs |
| --- | --- | --- | --- |
| `test_paper_model_based_configuration` | Run the recovered paper benchmark and save an isolated report. | Plot/save/failure/horizon/step/threshold options. | Structured `report`. |
| `validateBenchmark` | Check exact topology, parameters, initial state, and equilibrium. | Network and benchmark metadata. | Assertions. |
| `configureController` | Install current stabilizing, archived dissipative, or current dissipative gains. | Network, benchmark, method, threshold. | Design status and aggregate command gain. |
| `simulateLinearPaperCase` | Integrate the continuous-time linear error dynamics used by the old paper script. | Network, benchmark, gain, horizon, step. | Time, state, and metrics. |
| `referenceGainDistance` / `referenceTopologyMatch` | Compare a design with the accepted archived dissipative result. | Network, benchmark, method. | Relative distance / logical match. |
| `buildSummary` / `classifyResult` | Assign status and form the report table. | Case metrics. | Status / summary table. |
| `plotPaperTrajectories` / `plotPaperNetworkStates` | Plot per-unit responses and initial/terminal network states. | Completed cases and networks. | MATLAB figures. |

## Canonical main file

`main.mlx` is the only main file in the project root. A stale historical text export is retained as `archive/legacy/main_legacy_export.m` for provenance, but it is not a runtime entry point and should not be synchronized with the Live Script.

## Archive

Historical experiments, debug snapshots, and legacy artifacts are under
`archive/`. They are not normal runtime entry points. The paper regression
test uses `archive/data/matlab.mat` only as provenance for the frozen constants
already encoded by `createPaperDCMicrogrid`. See `archive/README.md`.

## Generated metadata

`.DS_Store` and `Results/.DS_Store` are macOS metadata only. They are safe to remove and should not be used by the simulation.

## Model-based baseline test

Run the independent simulator check from the project root with:

```matlab
addpath('tests')
report = test_model_based_simulator;
```

The test does not modify `netFile.mat` or `Results/TrajResults.mat`. Its isolated report and response comparison are written as `Results/ModelBasedSimulatorBaseline.*`; the initial/terminal DCMG visualization is written as `Results/ModelBasedSimulatorNetworkStates.*`.

## Recovered paper-configuration test

Run the six-DG model-based regression independently with:

```matlab
addpath('tests')
report = test_paper_model_based_configuration;
```

The report and figures use the `PaperConfiguration*` prefix and do not replace
the existing paper figures or `TrajResults.mat`. The archived dissipative
controller is the reproduction reference. The test warns when a current
design is stable but its communication topology differs from that reference.

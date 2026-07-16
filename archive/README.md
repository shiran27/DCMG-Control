# Archive Manifest

Files in this folder are retained for provenance and are not normal runtime
entry points. `data/matlab.mat` is important evidence for the recovered paper
benchmark, but the regression test does not load it at runtime.

| File | Description |
| --- | --- |
| `experiments/basicTest.m` | Standalone two-state identification/LQR example used to test a basic data-driven workflow; it does not instantiate the DC microgrid classes. |
| `experiments/temp.m` | Standalone block-diagonal stabilizability and obstruction experiment using YALMIP/MOSEK. |
| `experiments/debug.mlx` | Interactive debugging notebook for inspecting saved DG data matrices, disturbance estimates, and local/global data-driven designs. |
| `data/temp.mat` | High-resolution debugging workspace snapshot formerly consumed by `debug.mlx`. |
| `data/temp2.mat` | Lower-resolution debugging workspace snapshot formerly consumed by `debug.mlx`. |
| `data/matlab.mat` | Accepted six-DG, seven-line dissipative paper-run workspace (`N=6`, `M=7`, `t=0.05`), including the 15-link communication topology and controller gains used as the frozen regression reference. |
| `legacy/main_legacy_export.m` | Plain-text export of the later simplified four-DG `main.mlx`, including live-output annotations and embedded image data. Retained for comparison/history. |

None of these files should be added to the MATLAB path. The exact paper
constants are exposed through `support/createPaperDCMicrogrid.m`.

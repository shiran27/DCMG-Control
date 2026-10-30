# MATLAB support functions

`run_paper_results` adds this directory to the MATLAB path.

| Function | Role |
|---|---|
| `createPaperDCMicrogrid` | Existing fixed six-DG/seven-line benchmark; legacy equilibrium is optional. |
| `designSteadyState` | Exact bounded scalar convex equilibrium solve with proportional sharing. |
| `designConventionalDroop` | Fixed conventional V--I droop literature baseline and equilibrium-preserving bias. |
| `designLiteraturePI`, `tuneLiteraturePI` | Nahata et al. passivity-based PI controller and declared nominal-only tuning grid. |
| `dcmgSolverOptions` | Probe working MATLAB SDP solvers; project-local SeDuMi fallback. |
| `collectIntegralData` | Exact integral experiment, declared bounds and independent truth validation. |
| `designLocalIntegralDissipativity` | Hard MB or robust integral-DD local LMI. |
| `designGlobalDissipativity` | Hard known-B global composition with physical current rows fixed. |
| `verifyLocalDissipativity` | Original local storage inequality, using withheld truth for verification. |
| `verifyGlobalStorage` | Original aggregate bounded-real storage inequality. |
| `certificateDiagnostics` | Minimum eigenvalue, scale, symmetry and acceptance checks. |
| `captureControllerState`, `restoreControllerState` | Failed co-designs preserve the installed gains. |
| `simulatePaperCase` | Exact propagation of the existing linear plant and matched disturbances. |
| `verifyPaperImplementation` | Existing ODE/ZOH cross-checks, data-rank guard, model-withholding test. |
| `plotPaperResults` | Publication figures from the verified serialized result struct. |
| `writePaperResultsText` | Numerical table, manuscript prose and summary generated from results. |
| `verifyLiteraturePI`, `plotLiteraturePIComparison`, `writeLiteraturePIComparison` | Independent checks and publication artifacts for the literature PI baseline. |
| `verifyCompleteComparison` | Verify common provenance and matching trajectories before synchronization. |
| `syncLiteraturePIComparison` | Copy verified PI publication assets and keep a single manuscript include. |

The earlier random-network and trajectory-plot helpers remain available for
exploration. The paper workflow uses the fixed factory and a fresh output folder;
older figures and `TrajResults.mat` are historical artifacts.

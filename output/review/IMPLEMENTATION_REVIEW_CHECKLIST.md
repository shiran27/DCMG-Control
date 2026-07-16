# DCMG Controller Implementation Review Checklist

Prepared from the theory-to-code review of `P6_Data_Driven_Co_Design.pdf` and the current MATLAB implementation. This document combines the controller-design, DG/line-model, simulation, data-generation, and main-workflow findings.

Date: 2026-07-15

## How to use this document

- `[ ]` means a decision or change remains to be made.
- `[x]` means the implementation was checked and matches the stated paper model at the reviewed level.
- **P1** should be resolved before claiming the theoretical certificate is implemented.
- **P2** can materially change trajectories or the model/data-driven comparison.
- **P3** concerns numerical robustness, reproducibility, diagnostics, or maintenance.
- Search the MATLAB files for `% REVIEW PROPOSAL <ID>` to find disabled code suggestions.
- `ADD` proposals can normally be uncommented directly.
- For `REPLACE` proposals, disable the identified active line/block before uncommenting the replacement.
- `DESIGN` and `CHECK` proposals require a modeling choice rather than a blind code edit.

No proposal has been enabled. The baseline simulation behavior is unchanged.

## 1. Prioritized revision queue

Use this section as the implementation order. The later sections retain the
same revisions grouped by subsystem and provide the detailed rationale.

Effort estimates:

- **E1 - Quick:** localized guard, assertion, option, or diagnostic; normally under one hour including a focused test.
- **E2 - Moderate:** coordinated edits across a few methods or a new regression experiment; normally several hours.
- **E3 - Deep:** requires a theoretical/modeling decision, experiment redesign, retuning, and broad revalidation.

The groups are ordered by a combined criticality-and-ease score; ties favor
criticality. Items inside each group are also in recommended execution order.
Some E1 changes may expose infeasibility. That is a useful result, not a reason
to restore an invalid certificate.

### A. P1 + E1 - critical quick wins

- [ ] Abort MB/DD global design immediately when any local design fails; reject nonfinite or incorrectly signed `nu_i`. See `MB-G01` and `DD-G01`.
- [ ] Move every `value(...)`, minor check, and gain recovery after `sol.problem==0`.
- [ ] Assert recovered `Y22<0`, `X_i11>0`, `X_i22<0`, and finite/nonzero local dissipativity parameters.
- [ ] Assert the recovered model-based local closed loop is Hurwitz with a numerical margin. See `MB-L03`.
- [ ] Restore the hard, unsoftened MB and DD global LMIs for certificate-producing runs. See `MB-G03` and `DD-G05`.
- [ ] Report and require the minimum eigenvalue of each unsoftened certificate; enforce `gammaSq>0`, `lambda>=0`, and `Q>0`.
- [ ] Remove the redundant `1e6` multiplier from the DD physical equality. See `DD-G06`.
- [ ] Decouple DD `epsilon` from `gammaSq`; the current bounds force both to `1e-3`. See `DD-G04`.
- [ ] Stop perturbing the stabilizing GSC gain by up to 200% during data generation. See `GSC-01`.
- [ ] Disable independent network-current clipping in the certified linear run and retain it only as a separately labeled stress mode. See `DG-03` and `DATA-05`.
- [ ] Disable or explicitly label VSC command saturation as outside the linear certificate. See `DG-04` and `DATA-06`.
- [ ] Assert equilibrium dynamics and line-current residuals after every successful steady-state solve. See `SS-03`.

### B. P1 + E2 - critical moderate revisions

- [ ] Check `rank([U_i;X_i])` and its smallest singular value before creating DD YALMIP variables; include the uncontrollable line-current row. See `DATA-01` and `DD-L01`.
- [ ] Use the exact per-DG `QBar_w` normalization again when assembling Proposition 8. See `DD-L02` and `DD-G02`.
- [ ] Use one documented `Ts` and hold every input component over each identification interval, including Stages 1 and 5.
- [ ] Resolve whether `sigma_noise` is a covariance or a factor and make generation, logging, and documentation consistent. See `DATA-07`, `DATA-08`, and `MAIN-01`.
- [ ] Add Eq. (45) bounds on equilibrium line current and VSC command, and replace `+epsI` with a stated sharing target. See `SS-01` and `SS-02`.
- [ ] Verify the `QBar_w` sign convention and all required positivity conditions after scaling.
- [ ] Confirm that output data meet the paper's noise-free-output assumption or explicitly adopt the augmented measurement-noise model.

### C. P2 + E1 - important quick wins

- [ ] Handle an all-zero global gain and recheck stability after communication-gain thresholding. See `GSC-02` and `GSC-03`.
- [ ] Use exact-ZOH matrices for DD diagnostic poles and require the recovered sampled loop to be Schur stable. See `DD-L05` and `DD-L06`.
- [ ] Re-enable the physical consistency residual `D'*K-P*YBar*D'` and fail when it exceeds tolerance.
- [ ] Replace determinant/leading-minor diagnostics with minimum eigenvalues. See `NUM-01` and `NUM-02`.
- [ ] Use `norm(diag(epsilonSlack),1)` when an L1 sum of diagonal slacks is intended.
- [ ] Correct event times for nonzero `tspan(1)`. See `SIM-01`.
- [ ] Make the state-disturbance draws symmetric in both states. See `SIM-02` and `SIM-05`.
- [ ] Decide and encode whether heightened noise is a short pulse or a full-stage interval. See `SIM-03` and `SIM-06`.
- [ ] Initialize `comps`, deduplicate proximity/MST edges, and handle parallel lines consistently. See `MAIN-04`, `MAIN-05`, and `SYS-02`.
- [ ] Guard line-current per-unit division by zero or a near-zero equilibrium value. See `MAIN-07`.
- [ ] Create `Results/` before saving and assert compatible horizons before comparison. See `MAIN-08` and `MAIN-10`.
- [ ] Apply axis limits to explicit axes and make `adjustAxesPadding()` match its documented fractional units. See `MAIN-09` and `MAIN-11`.
- [ ] Add missing-noise guards, physical-parameter validation, endpoint validation, and coincident-position protection. See `DG-01`, `DG-02`, `LINE-01`, and `DRAW-01`.
- [ ] Remove the inactive DD BMI warning and guard dimension-dependent debug slices.
- [ ] Correct the `epsI` numeric format and resolve the documented `Irated=2*Pr/Vr` factor. See `MAIN-02`.

### D. P1 + E3 - critical design/theory work

- [ ] Define one theorem-consistent sampled experiment satisfying `Xbar=A_d X+B_d U+W` exactly.
- [ ] Specify an a-priori `Q_w` independently of the realized data and unknown `A/B`, with a justified reserve margin.
- [ ] Remove true-model `A/BBar` from DD synthesis inputs; use logged disturbance information with a consistent exact-discrete interpretation.
- [ ] Verify matrix S-lemma regularity/strict feasibility for every local QMI and the aggregate QMI.
- [ ] Decide whether Proposition 5 remains relaxed or Remark 10 is enforced; if enforced, iterate the BMI surrogate until the equality residual is acceptable. See `MB-L01`.

### E. P2 + E2 - important moderate revisions

- [ ] Add a bounded persistently exciting voltage-command signal instead of relying on initial transients/noise.
- [ ] Nondimensionalize states and inputs, then define dimensionless communication budgets and gain bounds shared by MB and DD.
- [ ] Correct local/global objective directions for `epsilon`, `trace(P)`, `gammaSq`, and `lambda`; add required upper bounds before maximizing margins.
- [ ] Use the same communication-cost budget for MB and DD before comparing performance. See `DD-G03`.
- [ ] Perturb physical `L`, `C`, `Rf`, and `RL` and rebuild the model instead of perturbing arbitrary entries of `A`. See `SIM-04`.
- [ ] Derive stage times and solver output spacing from inputs rather than hard-coding `dt=1e-5`.
- [ ] Add usable steady-state weights/bounds through `opts` and explicit converter command/current ratings.
- [ ] Prevent public line `R` and cached `g` from diverging. See `LINE-02`.
- [ ] Decide the active load ripple, initial-state range, bounded `alphaZ`, and converter-rating interpretation used in reported experiments. See `DG-05`, `MAIN-03`, and `MAIN-06`.
- [ ] Integrate logged process disturbance over each sample interval consistently rather than blindly multiplying by `Ts`.
- [ ] Store local/global gains, minimum eigenvalues, solver status, physical residuals, rank diagnostics, and spectral radii in returned `out` structs.

### F. P3 + E1 - robustness and documentation quick wins

- [ ] Record MATLAB/YALMIP/MOSEK versions, seed, noise interpretation, disturbance bound, and full DG/line parameters with every run.
- [ ] Correct or remove unused integral-action dimensions and fields (`z`, `K_I`).
- [ ] Remove clearly unused variables only after the final theory form is selected.
- [ ] Store unpruned and pruned gains separately with topology, cost, and recertification results.

### G. P2 + E3 - important experiment/model redesign

- [ ] If limiters are part of the claimed model, replace clipping with a physical limiter model and extend the analysis; otherwise keep limiter tests separate.
- [ ] Define structured component uncertainty, regenerate DD data after the event, and keep uncertainty experiments separate from nominal validation.
- [ ] Run a fair MB/DD comparison using the same plant, normalization, budgets, disturbance schedule, and communication constraints.

### H. P3 + E2/E3 - final reproducibility and paper pass

- [ ] Complete the controlled experiment protocol in Section 12 after Groups A-G are resolved.
- [ ] Store the minimum result record in Section 13 for every accepted run.
- [ ] Generate paper plots only from runs whose configuration and certificate diagnostics are stored with the trajectories.

## 2. Confirmed theory-to-code matches

- [x] `DG.updateModel()` implements paper Eq. (39) with state `x_i=[V_i;I_ti]`:
  - `A_i=[-1/(R_L C), 1/C; -1/L, -R_f/L]`.
  - `B_i=B_wi=[-1/C,0;0,1/L]` through `BBar=[E,B]`.
  - `theta_i=E*Ibar` has the correct constant-current-load sign.
- [x] `TransmissionLine.current()` implements `(V_i-V_j)/R_ij`.
- [x] `DCMicrogrid.dynamics()` computes exported network current as `YBar*V`, consistent with Eq. (41).
- [x] The local controller controls only the VSC voltage-command row; the uncontrollable line-current row is fixed to zero in local-gain decision variables.
- [x] `designLocalXiDissipative()` contains the continuous-time local LMI corresponding to Eq. (24) with `C_i=I`.
- [x] Its second local LMI follows the relaxed necessary condition in Eq. (21).
- [x] `designLocalXiDissipative_DataDriven()` follows the main block structure of Proposition 7, Eqs. (33)-(34).
- [x] Data-driven local gain recovery `K_i=Ktilde_i*(X_i*Khat_i)^(-1)` is implemented as right division.
- [x] `codesign_MB_DRC()` structurally enforces the physical line-current interconnection in the full pre-scaled gain.
- [x] `codesign_DD_DRC()` follows the major `Q`, `S`, and `R` block construction of Proposition 8, Eqs. (35)-(36).
- [x] `solveSteadyState()` enforces the equilibrium physics and proportional internal-current sharing in Eqs. (44), (46), and (47).

## 3. Shared theoretical assumptions

### Assumptions 1-3

- [ ] **P1:** Explicitly assert `Y22<0` and `X_i11>0`, `X_i22<0` after recovery, even when principal LMI blocks imply them numerically.
- [ ] **P1:** Confirm that the scalar-inner-block parameterization required by Assumption 3 is the intended restriction for every experiment.
- [ ] **P1:** Abort the global design if any local design fails. Current code warns and continues, which can lead to division by `nu_i=0`. See `MB-G01` and `DD-G01`.
- [ ] **P1:** Check every recovered `nu_i` is finite, nonzero, and has the sign required by `X_i11=-nu_i I>0`.

### Assumptions 4-5

- [ ] **P1:** Supply an a-priori quadratic disturbance bound satisfying Eq. (26). Do not infer the only bound from the same realized disturbance trajectory used for design.
- [ ] **P1:** Check the matrix S-lemma regularity/strict-feasibility condition for every local QMI and the aggregate QMI.
- [ ] **P1:** Verify `rank([U_i;X_i])=n_ui+n_xi` and report the smallest singular value for each DG. See `DATA-01` and `DD-L01`.
- [ ] **P1:** Confirm output data are noise-free as assumed by the paper. Process noise is allowed; measurement noise needs the augmented uncertainty model discussed in Remark 13.
- [ ] **P2:** Add a persistently exciting, bounded additive voltage-command signal during data collection instead of relying on initial transients and process noise.
- [ ] **P2:** Confirm the uncontrollable line-current component still provides sufficient independent data for the four-row `[U_i;X_i]` rank test.

## 4. Model-based local design

Function: `DG.designLocalXiDissipative()`

- [ ] **P1:** Decide whether to use the relaxed Proposition-5 condition or enforce the optional Remark-10 equality.
- [ ] **P1:** If enforcing Remark 10, iterate `xBar11ValGuess` until `xBar11*KHat=x21*KTilde` is satisfied. A warning alone does not restore the equality. See `MB-L01`.
- [ ] **P1:** Move all `value(...)`, minor, and gain recovery operations after checking `sol.problem==0`; failed solves currently produce meaningless post-processing.
- [ ] **P1:** Verify the recovered continuous closed loop is Hurwitz with a numerical margin. See `MB-L03`.
- [ ] **P2:** Revisit the objective. Minimizing `+epsilon` drives the strictness margin to its lower bound. See `MB-L02`.
- [ ] **P2:** Nondimensionalize voltage/current states and both input channels before choosing objective weights.
- [ ] **P3:** Remove unused variables (`YBar`, unused output selector variables) after the theory form is finalized.
- [ ] **P3:** Store `minEig(LMI1)`, `minEig(LMI2)`, BMI residual, gain norm, and solver residuals in `out`.

## 5. Data-driven local design

Function: `DG.designLocalXiDissipative_DataDriven()`

- [ ] **P1:** Make data satisfy one discrete equation `Xbar=A_d X+B_d U+W` exactly for the selected interpretation.
- [ ] **P1:** Hold all input components over `[kTs,(k+1)Ts)` during the identification/design trajectory.
- [ ] **P1:** Do not compute the required disturbance description from true `A` and `BBar` in a direct data-driven claim. See `DATA-02` and `DATA-03`.
- [ ] **P1:** Use the exact per-DG `QBar_w` scaling again in Proposition 8. See `DD-L02` and `DD-G02`.
- [ ] **P1:** Check Assumption 5 before creating YALMIP variables. See `DD-L01`.
- [ ] **P1:** Confirm the paper's `QBar_w` sign convention and strict S-lemma regularity numerically after scaling.
- [ ] **P2:** Replace the `1e-9` feasibility margin with a tolerance meaningful at normalized scale. See `DD-L03`.
- [ ] **P2:** Rework the objective so it does not minimize the margin. See `DD-L04`.
- [ ] **P2:** Use exact ZOH matrices for diagnostic eigenvalues, not forward Euler, when validating the ode45/ZOH implementation. See `DD-L05`.
- [ ] **P2:** Require the recovered sampled closed loop to be Schur stable. See `DD-L06`.
- [ ] **P2:** Remove the unused `xBar11ValGuess` warning while its BMI constraint is disabled; it currently suggests a consistency requirement that is not imposed.
- [ ] **P3:** Remove debugging slices such as `Q_w(1:5,1:5)` or guard them by dimensions.
- [ ] **P3:** Store `rank`, singular values, `lambda1`, `lambda2`, LMI minimum eigenvalues, and spectral radius in `out`.

## 6. Model-based global co-design

Function: `DCMicrogrid.codesign_MB_DRC()`

- [ ] **P1:** Abort immediately if any local Proposition-5 problem fails. See `MB-G01`.
- [ ] **P1:** Remove the negative soft diagonal slack from certificate-producing runs. `W-epsilonSlack>=epsilon I` is easier than `W>=0` when `epsilonSlack<0`. See `MB-G03`.
- [ ] **P1:** Only label a result Proposition-4 feasible if the unsoftened `W` has nonnegative minimum eigenvalue within solver tolerance.
- [ ] **P1:** Enforce/verify `gammaSq>0` when it represents `Y11`.
- [ ] **P2:** Optimize `gammaSq` if figures or text claim an optimized L2 gain; its objective coefficient is currently zero.
- [ ] **P2:** Reconsider `+trace(P)`. Remark 6 proposes `-phi(p)` to encourage larger `p_i`; minimizing trace has the opposite direction. Add `P<=pMax I` before using a negative trace term. See `MB-G02` and `MB-G04`.
- [ ] **P2:** Use `norm(diag(epsilonSlack),1)` if an L1 sum is intended; `norm(epsilonSlack,1)` is an induced matrix norm.
- [ ] **P2:** Confirm the communication-cost budget and gain bounds are dimensionless after state/input scaling.
- [ ] **P3:** Re-enable and assert the physical consistency residual `D'*K-P*YBar*D'`.
- [ ] **P3:** Store the unpruned gain, pruned gain, certificate eigenvalue, physical residual, and communication cost separately.

## 7. Data-driven global co-design

Function: `DCMicrogrid.codesign_DD_DRC()`

- [ ] **P1:** Abort if any local Proposition-7 design fails. See `DD-G01`.
- [ ] **P1:** Assemble aggregate `QBar_w` from the same normalized blocks used locally. See `DD-G02`.
- [ ] **P1:** Restore the hard Proposition-8 LMI for certificate-producing runs. See `DD-G05`.
- [ ] **P1:** Remove the `1e6` multiplier on the physical equality. It is algebraically redundant and can damage solver scaling. See `DD-G06`.
- [ ] **P1:** Decouple `epsilon` from `gammaSq`. Current bounds force both to exactly `1e-3`. See `DD-G04`.
- [ ] **P1:** Verify `Q>0`, `lambda>=0`, and the minimum eigenvalue of the unsoftened final LMI.
- [ ] **P2:** Use the same dimensionless communication budget as the model-based design before comparing performance. Current values are `1` and `1e-3`. See `DD-G03`.
- [ ] **P2:** Align the objective with Proposition 8/Remark 6; current `+trace(P)` and `+epsilon` directions are questionable. See `DD-G07`.
- [ ] **P2:** Either penalize `lambda` as documented or remove the claim that the objective makes it small.
- [ ] **P2:** Replace determinant-based diagnostics with minimum eigenvalues. See `NUM-01`.
- [ ] **P3:** Report how QMI normalization changes the recovered `lambda`; do not compare raw lambda values across differently scaled runs.
- [ ] **P3:** Validate the recovered gain against the known simulated model as a diagnostic only, clearly separated from the data-driven synthesis inputs.

## 8. Model-based stabilizing baseline

Function: `DCMicrogrid.design_MB_GSC()`

- [ ] **P1:** Do not perturb a stabilizing gain by up to 200% to generate data without rechecking stability. Prefer bounded additive PE input. See `GSC-01`.
- [ ] **P2:** Confirm the dense full gain satisfies the physical line-current equality before extracting the VSC command rows.
- [ ] **P2:** Rework the objective if epsilon is intended to measure robustness; `epsilon+trace(P)` minimizes epsilon.
- [ ] **P2:** Handle the all-zero gain case in communication thresholding. See `GSC-02`.
- [ ] **P2:** Recheck stability after thresholding small communication blocks. See `GSC-03`.

## 9. Equilibrium design

Function: `DCMicrogrid.solveSteadyState()`

- [ ] **P1:** Add Eq. (45) bounds on equilibrium exported current and VSC command. See `SS-01`.
- [ ] **P1:** Replace `+epsI`, which drives current utilization toward zero, with a stated target such as `(epsI-1)^2` or a justified alternative. See `SS-02`.
- [ ] **P1:** Assert equilibrium and line-current residuals after a successful solve. See `SS-03`.
- [ ] **P2:** Add optional weights/bounds through the currently unused `opts` input.
- [ ] **P2:** Add explicit converter command/current ratings rather than relying on later simulation clipping.
- [ ] **P3:** Correct the `z`/`K_I` documented dimensions or remove unused integral-action fields.
- [ ] **P3:** Use a numeric format such as `%.4f` when printing `epsI`, rather than `%s`.

## 10. DG and transmission-line runtime model

- [ ] **P1:** Remove independent line-current clipping for the linear/theory-consistent run. It can violate KCL because each DG clips a current already determined by `YBar*V`. See `DG-03` and `DATA-05`.
- [ ] **P1:** Treat VSC command saturation as an uncertified nonlinear stress test, or include a limiter model in the analysis. See `DG-04` and `DATA-06`.
- [ ] **P2:** Guard missing noise data for deterministic unit tests. See `DG-02`.
- [ ] **P2:** Decide whether `iload_fun` is active; it is configured but unused. See `DG-05`.
- [ ] **P2:** Validate positive `L`, `C`, `RL`, line `R`, and valid endpoints. See `DG-01` and `LINE-01`.
- [ ] **P2:** Prevent public `R` and cached `g` from diverging. See `LINE-02`.
- [ ] **P2:** Accumulate or prohibit parallel lines consistently. See `SYS-02` and `MAIN-05`.
- [x] `DG.draw()` and `TransmissionLine.draw()` now attach all text to the caller-provided axes and retain their original call forms.
- [ ] **P3:** Guard coincident DG positions in communication-topology drawing. See `DRAW-01`.

## 11. Simulation and main workflow

### Staged simulation

- [ ] **P1:** Use ZOH data collection in Stages 1 and 5 when the design interprets samples as discrete-time inputs.
- [ ] **P1:** Use a single documented `Ts` shared by data collection, design, ZOH simulation, and validation.
- [ ] **P2:** Correct event times for nonzero `tspan(1)`. See `SIM-01`.
- [ ] **P2:** Correct the state-disturbance expression; it currently only decreases voltage and never perturbs current. See `SIM-02` and `SIM-05`.
- [ ] **P2:** Decide whether heightened noise is a 2% pulse or a full-stage interval. See `SIM-03` and `SIM-06`.
- [ ] **P2:** Perturb physical `L`, `C`, `Rf`, and `RL`, then call `updateModel()`, instead of perturbing every `A` entry. See `SIM-04`.
- [ ] **P2:** Derive stage times and solver output spacing from inputs rather than hard-coding `dt=1e-5` inside `simulate()`.
- [ ] **P3:** Add defaults/validation for `linkFiltThresh`, `useData`, `Ts`, and `tspan`.

### Noise and disturbance data

- [ ] **P1:** Resolve whether `sigma_noise` is a covariance or a factor. Current code multiplies `randn` directly by it. See `DATA-07`, `DATA-08`, and `MAIN-01`.
- [ ] **P1:** Do not fit `Q_w` tightly to one observed residual without an a-priori margin. See `DATA-04`.
- [ ] **P2:** If using logged process disturbance in simulation, integrate it over the sample interval consistently rather than blindly multiplying by `Ts`.
- [ ] **P3:** Record the random seed, covariance/factor, realized disturbance norm, and chosen bound in every stored result.

### Main file and plots

- [x] `main.mlx` is the sole canonical workflow and contains the disabled `MAIN-*` review proposals directly.
- [x] The redundant root `main.m` was removed; the stale output-heavy historical export remains archived as `archive/legacy/main_legacy_export.m`.
- [ ] **P2:** Initialize `comps` before graph-growth attempts. See `MAIN-04`.
- [ ] **P2:** Deduplicate proximity and MST edges. See `MAIN-05`.
- [ ] **P2:** Decide whether `Irated=2*Pr/Vr` is intentional and correct the displayed rating comments. See `MAIN-02`.
- [ ] **P2:** Prefer bounded heterogeneous `alphaZ` if random load composition is restored. See `MAIN-03`.
- [ ] **P2:** Decide whether the 50%-100% initial-state range is too severe for the intended comparison. See `MAIN-06`.
- [ ] **P2:** Do not divide line current by a zero/near-zero equilibrium current. See `MAIN-07`.
- [ ] **P3:** Create `Results/` before saving on a clean checkout. See `MAIN-08`.
- [ ] **P3:** Make `adjustAxesPadding()` implement its documented fractional units. See `MAIN-09`.
- [ ] **P3:** Verify all compared result horizons and event schedules agree. See `MAIN-10`.
- [ ] **P3:** Apply axis limits to explicit axes handles. See `MAIN-11`.

## 12. Controlled experiment protocol

Use this order so each result isolates one source of behavior.

1. [ ] Save the current baseline outputs and solver logs unchanged.
2. [ ] Run deterministic model checks with noise, saturation, pruning, and parameter events disabled.
3. [ ] Verify `A`, `BBar`, `YBar`, equilibrium residuals, and current-sign conventions numerically.
4. [ ] Run MB-GSC without gain corruption and verify the exact closed-loop eigenvalues before/after pruning.
5. [ ] Run the model-based local/global design with hard LMIs and record all minimum eigenvalues.
6. [ ] Generate a ZOH PE dataset with fixed `Ts`, known seed, and independently specified `Q_w`.
7. [ ] Check Assumption 5 singular values and S-lemma regularity for every DG.
8. [ ] Run data-driven local designs and verify gain recovery, LMI eigenvalues, and sampled-model spectral radii.
9. [ ] Run data-driven global design with hard LMI, consistent QMI scaling, and the same communication budget as MB-DRC.
10. [ ] Compare MB and DD in the nominal linear case before adding parameter changes.
11. [ ] Add structured component uncertainty and regenerate DD data after the event.
12. [ ] Add saturation/noise stress tests separately and label them as outside the linear certificate where applicable.
13. [ ] Generate paper plots only from runs whose configuration and certificate diagnostics were stored with the trajectories.

## 13. Minimum result record for the paper

For every method/run, store:

- [ ] MATLAB, YALMIP, and MOSEK versions.
- [ ] Random seed and complete DG/line parameters.
- [ ] `Ts`, ODE tolerances, noise interpretation, and disturbance bound.
- [ ] Data ranks and singular values for every DG.
- [ ] Local `nu`, `rho`, gain, LMI minimum eigenvalues, and solver status.
- [ ] Global unpruned/pruned gains, topology, cost, `gammaSq`, lambdas, and physical residual.
- [ ] Unsoftened certificate minimum eigenvalue.
- [ ] Equilibrium states/inputs and residuals.
- [ ] Saturation activation counts.
- [ ] Event definitions and structured parameter changes.
- [ ] Trajectories using a documented, nonzero-safe per-unit base.

## 14. File organization completed

- [x] Key runtime files are documented in `README.md`.
- [x] Temporary experiments and data snapshots were moved under `archive/` and described in `archive/README.md`.
- [x] Paper figures and `TrajResults.mat` remain under `Results/` with a directory description.
- [x] macOS `.DS_Store` files are ignored for future additions.
- [x] Class headers describe the runtime/design role of each key file.

## 15. Final acceptance questions

- [ ] Can every theoretical assumption used in the paper be pointed to in code and to a stored numerical check?
- [ ] Does the data-driven design avoid using the true model except for explicitly labeled post-design diagnostics?
- [ ] Does the implemented sampled experiment satisfy the same data equation used by the theorem?
- [ ] Are reported certificates evaluated without negative feasibility slack and before any uncertified gain pruning?
- [ ] Are model-based and data-driven methods compared with the same plant, budgets, normalization, disturbances, and event schedule?
- [ ] Are nonlinear saturation and current-limiting effects either modeled or clearly separated from certified linear results?

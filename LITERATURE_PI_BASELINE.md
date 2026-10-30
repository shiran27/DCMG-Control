# Published PI baseline: implementation and comparison

The baseline implements the multivariable PI controller and explicit gain
region in Nahata, Soloperto, Tucci, Martinelli, and Ferrari-Trecate,
*Automatica* 113 (2020), 108770,
[DOI 10.1016/j.automatica.2019.108770](https://doi.org/10.1016/j.automatica.2019.108770).
The paper already has this entry under `Nahata2020` in its bibliography.
The architecture originates in Tucci et al. (2018); this implementation uses
Nahata et al.'s analytical passivity conditions, not Tucci et al.'s LMI design.

## Reproduce

```matlab
% Complete five-controller experiment, figures, tables, and regression checks:
report = run_all_paper_results; % main.mlx calls the same function
% Reuse the verified proposed gains, but regenerate every comparison response:
report = run_literature_pi_comparison;
% Publication assets and generated discussion:
sync_paper_results('/path/to/Data_Driven_Design_ACC2025');
```

The complete entry point builds the core methods and conventional droop first,
then the PI comparison. Every core method is included automatically. The source
record's date and hashes are retained, and shared trajectories are checked
before publication synchronization. Outputs go to `Results/LiteraturePI`,
including a complete MAT record, all 125 tuning attempts, CSV metrics, editable
MATLAB figures, PDF/PNG figures, and generated LaTeX/Markdown discussion.

## Mapping to the existing plant

For each DG, the published law (18)-(19) is

```
dz_abs/dt = Vref - V
Vt = k1*V + k2*It + k3*z_abs
k1 < 1, k2 < Rf, 0 < k3 < (1-k1)*(Rf-k2)/L.
```

We use `Vref=V_ss` and initialize the integral state at
`z_ss=(u_ss-k1*V_ss-k2*It_ss)/k3`. With `h=(z_abs-z_ss)/t0`,
`t0=1e-3 s`, the exact equivalent realization is

```
dh/dt = -(V-V_ss)/t0
Vt = u_ss + k1*(V-V_ss) + k2*(It-It_ss) + t0*k3*h.
```

The code augments the existing plant with six integrators. The physical
resistor network, ZI loads, converter command input, and disturbance channels
remain the existing ones. It uses continuous feedback, no saturation, and no
remote measurements. The original absolute-state controller is checked
numerically against this realization. Integral deviations are zero initially,
consistent with the operating-point initialization used by the other methods.

For `alpha=(k1-1)/L`, `beta=(k2-Rf)/L`, `gamma=k3/L`, and
`omega=gamma-alpha*beta`, the source's storage (38) is

```
Pphysical = [C, 0, 0;
             0, beta/omega, gamma/omega;
             0, gamma/omega, alpha*gamma/omega].
```

We transform it by `diag(1,1,t0)`, check its positivity and the local passive
port identity, and independently evaluate the full network dissipation matrix
and closed-loop poles. For the present static resistor network, physical
coupling contributes `-Verror'*YBar*Verror` to the storage derivative. This
checks the static-line specialization directly; no dynamic-line simulation
claim is made. The ZI load case has no constant-power nonlinearity.

## Fair tuning and interpretation

The declared grid is `kp,kd in {0.25,0.5,1,2,4}` and
`fraction in {0.1,0.25,0.5,0.75,0.9}`. It sets `k1=-kp`, `k2=-kd`, and
`k3=fraction*(1-k1)*(Rf-k2)/L`. Each augmented normalized gain norm is bounded
by 20. Candidates must be stable, converge over the nominal horizon, and obey
the larger nominal peak command deviation of the MB and DD co-designs.
Selection minimizes nominal 2% settling time, then state energy, then command
peak. The pulse, noise, and persistent-load test results are not used to select
gains. The final choice is benchmark-specific, not claimed to be the authors'
published numerical tuning.

Every method receives identical initial errors, pulses, state jumps, and noise
realizations. A further diagnostic adds a permanent 0.35 A current-load
increment at each bus at 25 ms. References stay fixed; tail metrics use
100--120 ms. This step is a regulation diagnostic, not an L2-gain test.
Restoring each fixed voltage after a changed load need not preserve proportional
sharing. The PI baseline has a voltage-stability/passivity certificate; a
global disturbance-to-state L2 bound is not synthesized for it.

Validation in `verifyLiteraturePI` compares exact propagation against the
existing `DCMicrogrid.dynamics` plus the published integral law, checks the
constant-load equilibrium analytically, confirms identical disturbance arrays,
and enforces the declared effort budget. All raw diagnostics are saved.

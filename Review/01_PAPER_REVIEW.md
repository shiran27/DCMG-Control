# Paper review — seven decisions

Reviewed 12 September 2026. Read R1–R7 in order beside the companion report; matching IDs refer to the same decision. **Approval means accepting the recommendation, not certifying the unchanged draft.** Automated checks support the evidence; scientific approval remains yours.

Open [the implementation checklist](</Users/shiranthamacbook/Library/CloudStorage/OneDrive-stevens.edu/Research/Papers/P6_DataDriven/MatlabCodesDataDriven/DCMG-Control/Review/02_IMPLEMENTATION_REVIEW.md>) alongside this file. Use [the current paper PDF](</Users/shiranthamacbook/Documents/GitHub/Data_Driven_Design_ACC2025/Revision/Paper_with_Regenerated_Results.pdf>) for figures; source links open the exact revision locations. Older uploaded papers and audit checkboxes describe historical versions.

## R1 — Resolve the global data-driven claim first · publication blocker

**Inspect:** [abstract / contributions](</Users/shiranthamacbook/Documents/GitHub/Data_Driven_Design_ACC2025/main.tex:113>), [Proposition 8 and its block definitions](</Users/shiranthamacbook/Documents/GitHub/Data_Driven_Design_ACC2025/main.tex:1173>) and [S-lemma proof step](</Users/shiranthamacbook/Documents/GitHub/Data_Driven_Design_ACC2025/main.tex:1352>); PDF pp. 1–3, 12–13. Also check [the matrix S-lemma hypotheses](</Users/shiranthamacbook/Documents/GitHub/Data_Driven_Design_ACC2025/main.tex:206>).

**Verify / recommendation:** Reject the current unrestricted global-DD claim until resolved. The last diagonal block is `−λ E Q̄w Eᵀ`; a positive QMI direction forces `λ=0`, then the zero block forces the gain-carrying cross block to vanish. The reported implementation avoids this construction. Reconcile the opening, theorem body, and conclusion around the supported local-DD / known-B global scope.

**Alternative:** Repair the lifted uncertainty argument, including strict/non-strict and repeated-parameter hypotheses, and then implement and validate the repaired theorem before retaining the broader claim.

**Decision:** ☐ Approve recommendation ☐ Reject ☐ Alternative / notes: ______

## R2 — Approve the actual continuous-time certificate

**Inspect:** [application contract](</Users/shiranthamacbook/Documents/GitHub/Data_Driven_Design_ACC2025/main.tex:1399>), [integral regression and certificate, Eqs. (52)–(54)](</Users/shiranthamacbook/Documents/GitHub/Data_Driven_Design_ACC2025/Results/Regenerated/numerical_results.tex:7>), and [short derivation](</Users/shiranthamacbook/Documents/GitHub/Data_Driven_Design_ACC2025/Revision/INTEGRAL_DATA_CERTIFICATE.md:53>); PDF pp. 15–16.

**Verify / recommendation:** Approve the exact integral-data S-procedure as a **sufficient continuous-time counterpart**. Check substitution, Schur complement, and gain recovery. The experiments use the basic MB local LMI, not the revisited compatibility/BMI construction. Decide whether to promote this implemented result into the theory section or a formal appendix.

**Alternative:** Implement the original discrete-time propositions with a justified sampling/input contract and regenerate their evidence.

**Decision:** ☐ Approve recommendation ☐ Reject ☐ Alternative / notes: ______

## R3 — Approve what “data-driven” assumes here

**Inspect:** [data acquisition and declared bounds](</Users/shiranthamacbook/Documents/GitHub/Data_Driven_Design_ACC2025/Results/Regenerated/numerical_results.tex:8>) and [data/uncertainty contract](</Users/shiranthamacbook/Documents/GitHub/Data_Driven_Design_ACC2025/Revision/INTEGRAL_DATA_CERTIFICATE.md:7>); PDF p. 16.

**Verify / recommendation:** Accept a simulation study with exact interval averages, noise-free measurements, 400 transitions, rank four, and a declared input-matrix norm bound of 1.25. Process noise is bounded before acquisition; measurement/quadrature error is not covered. Global composition and equilibrium design retain model knowledge.

**Alternative:** Require measured-data realism; add measurement/integration uncertainty to the certificate and rerun, rather than merely adding sensor noise to plots.

**Decision:** ☐ Approve recommendation ☐ Reject ☐ Alternative / notes: ______

## R4 — Approve the objective and zero-link co-design result

**Inspect:** [normalization and tuning](</Users/shiranthamacbook/Documents/GitHub/Data_Driven_Design_ACC2025/Results/Regenerated/numerical_results.tex:5>), [graph interpretation and Table I](</Users/shiranthamacbook/Documents/GitHub/Data_Driven_Design_ACC2025/Results/Regenerated/numerical_results.tex:23>); PDF pp. 16–17, Fig. 6.

**Verify / recommendation:** Retain the measured zero cyber links for both co-designs. Gamma ≈5.09 bounds normalized additive-state disturbance to state error; it is not a primitive-noise or sharing-error bound. The distance-weighted transformed-gain penalty is a surrogate for communication cost, not link-count optimization.

**Alternative:** Add a predeclared tradeoff sweep or a benchmark requiring remote feedback if topology design needs stronger empirical emphasis; do not insert artificial links.

**Decision:** ☐ Approve recommendation ☐ Reject ☐ Alternative / notes: ______

## R5 — Approve the literature comparison and its fairness

**Inspect:** [droop law and sensitivity](</Users/shiranthamacbook/Documents/GitHub/Data_Driven_Design_ACC2025/Results/Regenerated/numerical_results.tex:19>), [PI source, tuning, and Table II](</Users/shiranthamacbook/Documents/GitHub/Data_Driven_Design_ACC2025/Results/LiteraturePI/comparison.tex:3>); PDF pp. 16–19, Tables I–II and Fig. 9.

**Verify / recommendation:** Retain both baselines and the tradeoff discussion. Check cited equations, equilibrium-preserving droop bias, PI integral initialization, and nominal-only tuning under the 4.089 V peak-command budget. DD settles faster (3.2 versus 7.4 ms for PI) and has lower pulse sharing error; PI has lower voltage error. Droop uses less command effort.

**Alternative:** Use a common actuator/effort objective for every controller, or add a secondary voltage-restoration baseline. Either requires new tuning and comparisons.

**Decision:** ☐ Approve recommendation ☐ Reject ☐ Alternative / notes: ______

## R6 — Approve the operating envelope and experiment claims

**Inspect:** [ZI plant and equilibrium specification](</Users/shiranthamacbook/Documents/GitHub/Data_Driven_Design_ACC2025/main.tex:1375>), [archived initial state / pulse experiment](</Users/shiranthamacbook/Documents/GitHub/Data_Driven_Design_ACC2025/Results/Regenerated/numerical_results.tex:20>), [permanent-load interpretation](</Users/shiranthamacbook/Documents/GitHub/Data_Driven_Design_ACC2025/Results/LiteraturePI/comparison.tex:12>), and [variation study](</Users/shiranthamacbook/Documents/GitHub/Data_Driven_Design_ACC2025/Results/Regenerated/numerical_results.tex:35>); PDF pp. 14–19.

**Verify / recommendation:** Retain explicit linear-model limits: no saturation, static resistive lines, no CPLs, and an archived initial state containing negative voltages. Permanent-load PI voltage restoration does not imply proportional sharing; fixed-reference MB/DD offsets are about 1.14/1.09 V. Parameter variation is one empirical offline-redesign study, not online adaptation or a robustness certificate.

**Alternative:** Move the archived large-error figure to an appendix; require operating-envelope, multiple-seed, or nonlinear experiments before making broader claims.

**Decision:** ☐ Approve recommendation ☐ Reject ☐ Alternative / notes: ______

## R7 — Approve the publication package

**Inspect:** [generated section inputs and conclusion](</Users/shiranthamacbook/Documents/GitHub/Data_Driven_Design_ACC2025/main.tex:1575>); PDF pp. 15–21. Cross-check [the five-controller metrics](</Users/shiranthamacbook/Library/CloudStorage/OneDrive-stevens.edu/Research/Papers/P6_DataDriven/MatlabCodesDataDriven/DCMG-Control/Results/LiteraturePI/metrics.csv>) against both tables.

**Verify / recommendation:** After R1–R6, approve the generated figures/text and choose the final venue length. The current working paper is 21 pages. Edit numerical prose through its MATLAB writer so regeneration preserves revisions.

**Alternative:** Keep the full reproducibility material as supplementary content and shorten the main paper; defer submission while R1 remains unresolved.

**Decision:** ☐ Approve recommendation ☐ Reject ☐ Alternative / notes: ______

## Minor fixes already handled

The Schur lemma and corollary now state the positive-definite pivot required by their ordinary inverse. No global-DD theorem rewrite was made. Publication exports now regenerate from verified MAT records before synchronization; older audit documents point to this current review pair.

**Verification, 12 September:** The full MATLAB run and regression gates passed. Both numerical sections are included once; all 19 synchronized assets match. The corrected 21-page PDF compiled without undefined references and its changed/result pages were visually checked. R1–R7 remain for your approval.

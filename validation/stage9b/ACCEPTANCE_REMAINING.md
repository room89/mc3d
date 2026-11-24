# Stage 9b: remaining acceptance work

Physical acceptance is **not complete**. Software tests, synthetic coverage
screens and completed solver executions are separate evidence. Original gates
and the full PHYSICS_PLAN.md scope remain in force.

## Current evidence

- The primary Alsmeyer PDF is available. Its hash and both committed primary
  digitizations were reverified during pilot preparation. Source acquisition is
  no longer a blocker. Digitization uncertainty, scan deformation and matched-
  condition comparison remain unqualified; pointwise interval overlap is not
  experimental acceptance. See [primary evidence](EXPERIMENT_SOURCE_COMPARISON.md).
- #70–77 are merged. The fixed-coordinate arithmetic estimator and conservative
  band passed prescribed synthetic screens. Those results do not qualify
  physical noise, stationarity, intrinsic shock shape or gradient thickness.
- Historical fixed-coordinate profiles show appreciable run variation and a
  mean-front offset. Reservoir input checks and isolated flux checks found no
  gross normalization error at their tested precision. Three-seed historical
  C/E differences do not establish time-step convergence.
- The [new independent C/E pilot](TIMESTEP_PILOT.md) is specified to measure
  run variation and cost, not to pass acceptance. Inspect its actual process
  and report state before claiming progress or completion.
- The matched 40-experimental-mean-free-path matrix has completed short runs
  and a verified 1024-window SPARTA C run. The long MC3D C160 comparison covers
  896 matched windows in the [matrix report](ALSMEYER_MATRIX_RESULTS.md).
  The original window plateau criterion remains open. Gradient thickness and
  numerical convergence also remain open.

## Required next decisions and evidence

1. Complete the frozen C/E pilot; verify all runs, source hashes and full-window
   coverage. Report per-run means and half changes, treating runs rather than
   windows as independent units. Retain negative and inconclusive results.
2. Design a prospective comparison with a justified uncertainty procedure,
   sample size and stopping rule. Do not tune on confirmation results. Resolve
   the relationship between fixed-frame profiles and the intrinsic thickness
   required by experimental validation; the density band alone is insufficient.
3. Establish stationarity and grid, time, population and domain convergence,
   followed by inter-code and primary-experiment comparisons in matched conditions.
   Preserve numerical, statistical and experimental uncertainties separately.
   Verify the complete long C160 history and match it to SPARTA, then compare
   against the long C320 run. The uniform downstream controls will test
   whether the outlet excursions also occur in the absence of the shock.
4. After stage 9, complete the isothermal wedge and dynamic-wall checks in
   stages 10–11 and the final limitations/reproducibility report. They remain
   required for the user's full physical acceptance objective.

Apply the current user-approved load limits to prepared launches. Inspect
the solver and output manifests before resuming interrupted work.

## Evidence and provenance

Negative estimator results and historical investigation remain in the
[original PR snapshot](https://github.com/room89/mc3d/tree/2c58c95701e45c7e3b196e7e86c441fc8496b594/validation/stage9b).
See [archive and file disposition](PR59_ARCHIVE.md). Historical estimates are not
fresh confirmation samples. No live progress counters are recorded here.

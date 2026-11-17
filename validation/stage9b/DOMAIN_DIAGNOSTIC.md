# Longitudinal boundary diagnostic

Delayed sampling at [20,28] original-domain transit times did not establish
stationarity. A descriptive translation reduces the late inter-code density
maximum from 0.080742 to 0.005402, but shifts are fitted from the same noisy
profiles. This is not acceptance or evidence that boundary effects are absent.

Freeze a six-run E diagnostic: both pinned binaries from WARMUP_DIAGNOSTIC.md,
seeds 210001, 210002, 210003, warmup 20, domain multiplier 2. Double longitudinal
length and cell count (800 to 1600), preserving cell size, model-particle weight,
reservoirs, collision/transport settings, transverse dimensions and physical
sampling interval. The reference transit time remains that of the original
40-collision-mean-free-path domain; it is not recomputed from the extended box.
The new box has 80 collision mean free paths. This is an engineering diagnostic,
not the experiment-matched matrix.

Command: `shock_validation.py --variants E --seeds 210001 210002 210003
--warmup-transits 20 --domain-multiplier 2 --load-mode full --no-time-budget`,
with explicit pinned binary paths and a new output directory. Retain all runs,
64 windows and unsuccessful results; no favorable-window selection, automatic
retry or early stopping. Changing domain on resume must fail before rewriting
existing evidence.

Compare only common original-coordinate nodes in [-14,14] collision mean free
paths. Report within-run half changes, whole means, descriptive front/width,
and differences from the completed multiplier-1 late campaign. No interpolation
is needed for common nodes. Keep translation diagnostics separate from the
fixed-grid comparison. Shared seed numbers do not imply identical initial
particle arrangements or paired random trajectories after the domain changes.
Three seeds supply exploratory evidence, not qualified coverage or confirmation.

A reduced discrepancy would motivate a prospectively specified domain-convergence
campaign; it would not establish convergence from two boxes. No reduction would
redirect investigation toward front motion, statistical variation and reservoir
implementation. Neither outcome alone qualifies intrinsic thickness or closes
stages 9–11. No physical tolerance is relaxed.

## Completed diagnostic

All six runs in `output/shock-domain2-warmup20-v1` completed. All recorded
run hashes and the 64 expected SPARTA sampling files per run were checked.
The fixed-grid analysis rechecked 1621 evidence files. Its report is
`output/shock-domain2-analysis-v1.json`, SHA256
`9df652991179d36f5ef1f681591d944f2904548832754586e99c32279013e265`.

| Descriptive density metric | MC3D | SPARTA |
| --- | ---: | ---: |
| Maximum change of ensemble mean, large minus small domain | 0.048354 | 0.029211 |
| Maximum ensemble half-run change, large domain | 0.011360 | 0.025962 |
| Maximum ensemble half-run change, small domain | 0.014218 | 0.016342 |

The maximum fixed-grid inter-code difference increased from 0.080742 to
0.099707. Doubling the domain did not remove the discrepancy or establish
stationarity. Three seeds and two domains do not establish a systematic
boundary bias or domain convergence.

A separate exploratory translation diagnostic
(`output/domain2-translation-diagnostic-v1.json`) finds an ensemble front
offset of 0.404690 collision mean free paths. Translating the SPARTA mean by
this fitted offset reduces the maximum difference on the common interior
[-11.975,11.975] to 0.005713. Ensemble 10–90 widths are 4.673090 (MC3D) and
4.650044 (SPARTA), in collision mean free paths. These are descriptive
crossings of isotonic mean profiles, not qualified intrinsic widths.
Half-run front changes range from -0.110746 to 0.019731 in MC3D and from
-0.013431 to 0.235715 in SPARTA. Translated half-run residual maxima remain
0.0112–0.0130. Same-data translation supplies neither an acceptance correction
nor uncertainty coverage.

The next investigation must distinguish front motion, between-run variation
and reservoir effects before another domain enlargement. Preserve the raw
fixed-coordinate comparison and existing tolerances. Stages 9–11 remain open;
this engineering diagnostic is not an experiment-matched acceptance result.

### Boundary bookkeeping follow-up

`output/domain2-boundary-exchange-v1.json` audits all three MC3D runs using
whole steps over the full run and eight fixed time sections. Particle-count
balances close exactly; momentum and energy residuals are below 1e-10 of
the summed absolute incoming/outgoing exchange. Integrated incoming counts
differ from the analytic reservoir expectation by -0.755741 model particles
in every run, consistent with a sub-particle injection remainder.
Full-run incoming momentum deviations in analytic standard-error units are
-2.5263, -0.2704 and 0.2104; energy deviations are -1.0677, -1.7264 and 0.1473.
The largest absolute sectional momentum deviation is 3.8163. These are
uncalibrated descriptive diagnostics under adaptive stepping, not simultaneous
significance tests. They do not identify a bookkeeping defect or exclude a
boundary-model effect. Equivalent SPARTA per-step exchange evidence is absent.

### Full-domain density decomposition

The follow-up `output/domain2-full-regions-v1.json` (SHA256
`566a6794bb3ea9204b175e973df229f209ee092f894492cba6597e5110566554`)
re-extracts all 64 windows and 1600 cells per run and verifies the 1621
campaign evidence hashes before and after analysis. The table integrates
second-half minus first-half normalized density in collision-mean-free-path
coordinates. Central means [-14,14]; the outer regions complete the box.

| Code / seed | Outer upstream | Central | Outer downstream | Total |
| --- | ---: | ---: | ---: | ---: |
| mc3d / 210001 | -0.005859 | 0.125281 | -0.046933 | 0.072489 |
| mc3d / 210002 | 0.002940 | -0.007488 | 0.028557 | 0.024010 |
| mc3d / 210003 | 0.001590 | 0.013746 | 0.010566 | 0.025902 |
| sparta / 210001 | -0.004129 | -0.084020 | -0.003601 | -0.091750 |
| sparta / 210002 | -0.004721 | -0.179137 | 0.018604 | -0.165253 |
| sparta / 210003 | 0.001036 | -0.016047 | -0.018822 | -0.033833 |

Outer downstream changes can compensate central changes or reverse the
sign of their total (MC3D seed 210002). Total density integrals increase in
these three MC3D runs and decrease in these three SPARTA runs. This small
exploratory sample does not establish a systematic drift or a reservoir bug.
A central mass estimate cannot be treated as an independent front position.
The next control must distinguish population fluctuations and boundary
exchange from front translation, retaining full-domain evidence. No physical
gate is closed by this decomposition.

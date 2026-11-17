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

# Delayed sampling diagnostic

The completed C/E pilot samples 64 equal windows over transit times [4,12].
Its three-seed SPARTA E mean changes by up to 0.02545 in normalized density
between halves. This does not establish a stationary fixed-coordinate profile.

Run E again with the same seeds 210001, 210002, 210003 and both pinned pilot
binaries, but sample [20,28]. Keep the original grid, reservoir states, particle
weight, timestep settings and eight-transit sampling duration. This is a
paired diagnostic using previously inspected seeds, not independent confirmation.
The MC3D binary hash is
`330d9ed03dea9c36e39c281edbd4150aeb1c4b8e186003316ffa72f432f489b1`;
SPARTA is `0d6dedda6d90b1f703d5cd0f745dbeb976cac7f74bf57e9a13016b26e8b10afb`.

Use `shock_validation.py --variants E --seeds 210001 210002 210003
--warmup-transits 20 --load-mode full --no-time-budget` with the pinned binary
paths and a new output directory. The new option defaults to four transit times;
production always has 64 windows spanning eight transit times. Resume rejects
changed warmup. SPARTA metadata records the warmup window count; both readers
use it, with 32 as the historical fallback.

Retain all six runs and all windows. Compare original-coordinate whole means,
first/second-half changes, front locations and descriptive widths with the
completed pilot. Report both absolute translation and shape evolution; do not
silently align profiles or select favorable windows. Completion and legacy
analyzer exit codes do not establish physical acceptance. No automatic retries
or stopping based on a favorable profile.

A smaller late drift would support a transient explanation, not prove
stationarity or boundary independence. A subsequent domain comparison remains
necessary. The experiment-matched matrix and wedge stages remain outstanding.

## Completed diagnostic

All six solver runs completed in `output/shock-warmup20-v1`; supervisor returned
2 (physical acceptance incomplete). The frozen comparison completed successfully
and verified 1621 campaign files before and after reading. Its report is
`output/shock-warmup20-analysis-v1.json`, SHA-256 `464924e067f74fe046e77d03d4ec4e0750f09c9747ede8547160cb98a58bdf2a`.
Separate input audits confirmed that all three MC3D SI inputs change only
observation/end times and all three SPARTA decks change only the warmup run.
Species and collision files are unchanged.

| Descriptive fixed-grid quantity | MC3D | SPARTA |
| --- | --- | --- |
| Max absolute ensemble half change, early [4,12] | 0.010701 | 0.025446 |
| Max absolute ensemble half change, late [20,28] | 0.014218 | 0.016342 |
| Max absolute late minus early ensemble profile | 0.009100 | 0.056767 |

The late inter-code maximum difference is 0.080742, compared with 0.028512
in the early E pilot. These are normalized-density differences on the same
560 original nodes; no alignment, interpolation or confidence claim is applied.
Later sampling alone does not establish stationarity or solve the fixed-frame
comparison. The SPARTA late run-mean fronts are 0.35769, 0.66316, 0.52447
collision mean free paths; MC3D fronts are 0.15519, 0.11642, 0.30274.
Isotonic crossings are descriptive and cannot qualify intrinsic shock width.

The same seeds are reused, not fresh confirmation. Wider run variation and
front displacement require distinguishing stochastic front motion, residual
transients and boundary effects. Do not increase warmup repeatedly until a
favorable comparison appears. Next investigate translation versus shape and
boundary sensitivity under an explicit protocol; retain all current evidence.
The experiment-matched convergence matrix and stages 10–11 remain outstanding.

The host elapsed clock showed a large jump during the second SPARTA run;
no cause was established. Do not interpret total elapsed time as CPU cost.

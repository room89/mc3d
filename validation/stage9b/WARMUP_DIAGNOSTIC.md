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

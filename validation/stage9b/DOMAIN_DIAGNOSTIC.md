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

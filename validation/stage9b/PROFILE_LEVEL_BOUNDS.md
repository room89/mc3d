# Width bounds from a simultaneous nodal profile band

This provides a conditional uncertainty transformation, not a new physical
acceptance gate. It avoids estimating a translation from noisy data and avoids
assuming a linear profile between grid nodes. The target is the width of the
same continuous, nondecreasing mean profile covered by the input band. It is
not the average instantaneous width or an intrinsic width after random front
motion has broadened the mean.

## Construction and coverage argument

Given ordered nodes x_i and simultaneous bounds L_i <= f(x_i) <= U_i, tighten
L by the prefix maximum and U by the suffix minimum. Every nondecreasing
admissible profile satisfies these tightened bounds. Reject an empty feasible
set; do not isotonic-fit incompatible limits into a feasible band.

For crossing q(p)=inf{x: f(x)>=p}, the last node with U_i<p is a lower outer
bound, and the first node with L_i>=p is an upper bound. Both must exist;
otherwise report insufficient bracketing rather than extrapolate. Continuity
and monotonicity constrain the crossing between those nodes, with no chosen
subgrid interpolant. Closed endpoints deliberately enclose limiting cases.

If crossing intervals for p=.1 and .9 are [a,b] and [c,d], the width lies in
[max(0,c-b),d-a]. This can be conservative because the crossings share a profile.
Both are transformed from the same simultaneous coverage event; no extra
multiple-testing correction or independence between crossings is needed.
Consequently coverage is at least the input coverage **conditional on the
assumed target monotonicity and input coverage**. The tool neither establishes
those assumptions nor assigns a new confidence level.

Constant translation of all nodes translates crossing bounds and leaves width
bounds unchanged. Random shifts inside an ensemble mean can change its width;
this invariance must not be mistaken for recovery of a latent intrinsic shape.
A nodal band cannot bound arbitrarily large subgrid slopes without an additional
regularity model, so gradient thickness remains explicitly outside this tool.

## Resolution/precision example

For f(x)=(1+tanh(x/2))/2, exact width is 2 log(9)=4.394449 λ₁.
Using x from -8 to 8, the outer intervals are:

| Node spacing | Exact nodal values | Nodal band ±0.005 |
| --- | --- | --- |
| 0.25 | [4.0, 4.5] | [4.0, 5.0] |
| 0.05 | [4.3, 4.4] | [4.2, 4.6] |
| 0.025 | [4.35, 4.40] | [4.25, 4.55] |

These are deterministic illustration bands containing the truth, not empirical
coverage estimates. They show why a ±0.005 density band alone need not meet the
existing width precision requirement. Increasing runs cannot remove subgrid
ambiguity; interpolation or smoothness assumptions must be justified separately.
No tolerance is relaxed and no physical profile is accepted by this example.

## Validation and next use

All 228 Python tests passed. New checks cover exact logistic width, translation
invariance, subgrid ambiguity, monotonicity tightening and impossible bands,
unbracketed levels, and an enumerated family of admissible piecewise-linear
profiles. `tools/profile_level_bounds.py` does not read physical data.

Next evaluate this transformation under independently specified synthetic
motion and time-evolution scenarios, preserving the distinction between mean
width and intrinsic shape. Apply it to physical data only with justified band
coverage, domain and monotonicity assumptions. A small translated residual from
the pilot is not such justification.

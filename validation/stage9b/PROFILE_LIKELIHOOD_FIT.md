# Unknown-profile likelihood prototype

Physical acceptance remains incomplete. The existing alignment-order estimator
and acceptance gates are unchanged. This prototype fits a profile with **known**
noise covariance and shift probabilities; it does not estimate either from
physical data and does not produce confidence bands.

## Model and algorithm

Independent runs contain correlated time windows. The Gaussian block covariance
is supplied as `T ⊗ S`. The fitter enumerates the supplied joint shift
configurations and marginalizes each whole block; it never treats dependent
windows as independent observations. Prediction matrices interpolate the latent
profile, not noisy observations.

The profile has 33 knots on [-8, 8], fixed endpoint values 0 and 1, monotone
interior values and constant plateau continuation. EM alternates posterior shift
weights with a constrained quadratic profile update. The update uses projected
gradient descent and bounded isotonic projection; it adds no ridge term.

The M-step limit is 10,000 iterations, with projected-gradient residual divided
by the Hessian spectral norm <= 1e-8. The EM limit is 200 iterations; convergence
requires three consecutive changes below 1e-7 in the mean negative log
likelihood (constant omitted) and 1e-6 in the profile. Objective increases,
iteration limits and covariance/support failures remain explicit. Three starts
are retained: linear, projected observation mean, and its three-point smoothed
version. Converged starts within 1e-7 objective but over 1e-4 profile difference
are flagged ambiguous. These numerical checks do not prove global optimality.

`fit(observations, x, temporal, spatial, configurations, weights)` returns the
profile and all start diagnostics. Observations have dimensions run × time ×
space. Common weights or per-run weights are accepted; the latter support a
known-shift oracle. No true-profile or family argument is accepted. The hard
assignment control optimizes its classification objective, separately from the
marginal objective.

## Frozen synthetic comparison

The runner freezes source hashes, NumPy version, seed and specifications in
`plan.json` before computation. Resume requires an identical plan and validates
checkpoint identities, profile constraints and recomputed error metrics. Results
are exclusive-write files; existing evidence is never overwritten.

Main matrix: 192 trials = 2 shapes × 3 motion amplitudes × 2 noise levels × 16
repetitions. Each trial has 32 independent runs, two windows and 41 observation
nodes on [-10, 10]. Shapes are sampled at profile knots and given explicit
plateau endpoints, making the main truth representable by the model. A separate
16-trial analytic-curve control tests model approximation; it cannot establish
matched-model superiority.

Each window shift is uniform on {-a, 0, a}, with a = 0, .05 or .5; duplicate zero
configurations are removed. Noise has temporal covariance
`run_sigma² + window_sigma² * .8^|t-u|` and spatial covariance
`scale[i] * scale[j] * .5^|i-j|`, where scale increases linearly from 1 to 3
across the observation grid. Low noise uses (.01, .02) run/window amplitudes;
high noise doubles both. SeedSequence uses [202611120000, trial_index] with
separate motion and noise child streams. These are specified synthetic laws,
not fitted physical fluctuation models.

Marginal, hard-assignment and known-shift methods see identical observations.
Truth is used only by the generator, error evaluator and known-shift oracle.
All starts, signed error curves, maximum error, integrated squared error (ISE)
and failures are retained. Comparison uses paired ISE differences, stratified
Monte Carlo standard errors and descriptive normal intervals. Advancement
requires lower mean ISE for both shapes at nonzero motion, no additional
failures, complete successful pairs, and intervals wholly below zero. Failure
or inconclusive results do not trigger radius inflation or threshold changes.

## Reproduction and load

From the repository root with NumPy 2 installed:

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 \
python tools/run_profile_likelihood_screen.py \
  --directory output/profile-likelihood-screen --prepare-only
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 \
python tools/run_profile_likelihood_screen.py \
  --directory output/profile-likelihood-screen --resume
```

The runner uses one process and lowers its scheduling priority with `nice(19)`.
`--max-new-trials 1` limits a pilot without changing the frozen matrix. The
full screen comprises 208 trials. Raw results stay in local `output`; only a
compact outcome belongs in this document. Old expensive series are not resumed.

## First frozen screen outcome

The local `output/profile-likelihood-screen-v2` completed all 192 matched trials
and 16 analytic controls. Every method converged in all 208 trials. All checkpoint
hashes, error metrics and the aggregate report were read back and verified.

For the 64 nonzero-motion matched trials per shape:

| Shape | Marginal mean ISE | Hard mean ISE | Known-shift mean ISE | Paired difference ± 1.96 MCSE |
| --- | ---: | ---: | ---: | --- |
| Logistic | 0.00121451 | 0.00134306 | 0.00109147 | -0.00012855 ± 0.00009524 |
| Gompertz | 0.00121179 | 0.00119385 | 0.00110483 | +0.00001795 ± 0.00005662 |

Decision: **no_demonstrated_improvement** across both families. Logistic improves
in this specified comparison; Gompertz has no demonstrated benefit. Its interval
also does not establish degradation. The prototype is retained for research,
but this screen does not justify advancing it to noise-model fitting or replacing
the current estimator. No criteria were changed after observing results.

Plan SHA-256:
`3408280936f4f44163b16fe23a2784077a8c06a8bb3480b20988fcfe02d4a86c`.
Complete report SHA-256:
`f4a094c4b46172430797d6507457723ddea438065ade03558908e84198ef0be0`.
The report contains individual hashes for all 208 trial files. The preliminary
one-trial `v1` run only checked execution and zero-motion agreement; it is not
additional independent evidence. The v2 protocol adds checkpoint validation.

## Post-hoc diagnosis: Gompertz

The read-only diagnostic verifies every source checkpoint hash and recomputes the
screen summary before comparing scenarios. It introduces no new evidence and
leaves the original decision unchanged. All intervals below are descriptive,
unadjusted for multiple scenario comparisons, and based on 16 paired trials.

| Motion | Noise | Marginal minus hard mean ISE | Difference ± 1.96 MCSE |
| --- | --- | ---: | ---: |
| .05 | low | +0.00001574 | ±0.00002367 |
| .05 | high | -0.00005396 | ±0.00008816 |
| .5 | low | -0.00000432 | ±0.00000579 |
| .5 | high | +0.00011432 | ±0.00020718 |

Every Gompertz scenario interval spans zero. The positive aggregate difference
is concentrated in the high-motion, high-noise case, where 9/16 paired differences
are positive and 7/16 negative. Trial 188 contributes +0.00141130 to the sum of
paired differences; it remains included, along with every other trial. Selecting
or deleting trials after observing their error is not a correction.

In that case the known-shift oracle has mean ISE 0.00136903, versus 0.00173362
for marginal fitting and 0.00161929 for hard assignment. Marginal-minus-oracle is
+0.00036459 ± 0.00022777. The oracle uses the same profile representation and
observations: the gap is consistent with loss from unknown shifts in this
specified model, not evidence of an error in the physical solver.

The three marginal starts differ by at most 0.00001527 across all matched
Gompertz cases; maximum spread is under 0.00000027 at motion .05. This gives no
visible evidence of distinct-start minima explaining the Gompertz outcome;
three starts cannot certify global optimality. The hard-assignment control is
more initialization-sensitive, reaching a spread of 0.03873822. Separately,
logistic at high motion/noise does show distinct marginal solutions in trials
82, 90 and 95 (maximum spread 0.01494741), so global optimization is not established
for the overall method.

For high-motion/high-noise Gompertz, the empirical error decomposition gives a
marginal-minus-hard change of -0.00003037 in squared mean-error energy and
+0.00014469 in scatter energy. This is an exact decomposition of the observed
sample ISE, not an unbiased squared-bias estimate. It does not support describing
the positive ISE difference as demonstrated extra systematic bias.

Conclusion: retain draft status and the failed advancement decision. There is
no identified Gompertz implementation defect to repair from these observations.
A future confirmation would require a separately frozen independent comparison,
with sample size based on Monte Carlo precision and all four nonzero-motion
Gompertz scenarios retained. Do not tune this estimator on the current held-out
trials or infer qualification of density coverage from mean ISE.

Reproduce without rerunning the screen:

```sh
python tools/diagnose_profile_likelihood_screen.py \
  --directory output/profile-likelihood-screen-v2 \
  --output output/profile-likelihood-diagnosis.json
```

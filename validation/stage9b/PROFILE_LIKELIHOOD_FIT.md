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

# Independent confirmation of the profile likelihood prototype

Physical acceptance is incomplete. This is a fixed-size synthetic comparison of
mean integrated squared error (ISE), not qualification of density confidence
bands. The estimator, generator and their settings are unchanged from PR #67.

## Frozen design

The pilot `profile-likelihood-screen-v2` is used only for variance/precision
planning. Its 208 observations are never pooled with the confirmation sample.
SeedSequence inputs use the new seed family 202611140000 and trial index, with
independent motion/noise child streams. All methods in each trial share the same
observations. The frozen plan records source hashes, NumPy version and pilot
plan/report hashes before the first confirmation trial.

Both logistic and Gompertz retain all four combinations of motion (.05, .5) and
noise (low, high). Each trial has 32 independent runs of two correlated windows;
the known covariance, monotone 33-knot model and supplied shift probabilities
are exactly those of the pilot. Marginal fitting, hard assignment and the
known-shift oracle are compared. Zero motion and analytic-curve controls were
already checked in the pilot; the confirmation targets nonzero-motion matched
models only.

The precision target is a confidence-interval halfwidth equal to 2.5% of the
pilot hard-assignment mean ISE for each family. This is a Monte Carlo design
choice, **not a new physical acceptance tolerance**. For four equally weighted
strata with pilot paired-difference variances s_j², the required repeats per
stratum are `ceil(z² * sum(s_j²) / (16 * target²))`. Here
`z = NormalDist().inv_cdf(.9875) = 2.2414027276` provides a two-family Bonferroni
adjustment for approximate two-sided normal intervals. Counts are rounded up to
multiples of 16 with a minimum of 32.

| Family | Target halfwidth | Formula count | Frozen repeats/stratum | Trials |
| --- | ---: | ---: | ---: | ---: |
| Logistic | 0.0000335765 | 169 | 176 | 704 |
| Gompertz | 0.0000298462 | 76 | 80 | 320 |

Total: **1024 new trials**. The count is fixed; there is no optional success
stopping, post-hoc stratum exclusion or automatic sample-size extension.

## Decision and integrity

The family mean difference is the equally weighted mean of four paired ISE
stratum means. Its Monte Carlo SE is `sqrt(sum(s_j²/n))/4`. A research advancement
requires all trials complete, every method converged, both observed interval
halfwidths within their targets, and both adjusted upper endpoints below zero.
Otherwise the result explicitly records fit failures, inconclusive precision,
no advantage for all families, or inconclusive advantage. Failure rows are never
silently excluded to obtain a pass. Normal intervals remain approximate;
passing this comparison is not physical acceptance or coverage qualification.

A resume validates checkpoint identities and recomputed profile error metrics.
The completed report contains every checkpoint hash; `--audit-only` reads all
files back and verifies hashes, metrics and the complete decision. Source or
NumPy changes prevent continuation. Existing reports and trials cannot be
overwritten. Interim progress is operational information, not an early decision.

## Reproduction and load

With NumPy 2 installed, set `OPENBLAS_NUM_THREADS=1`,
`VECLIB_MAXIMUM_THREADS=1` and `OMP_NUM_THREADS=1`, then run:

```sh
python tools/confirm_profile_likelihood.py \
  --pilot output/profile-likelihood-screen-v2 \
  --directory output/profile-likelihood-confirmation-v1 --prepare-only
python tools/confirm_profile_likelihood.py \
  --directory output/profile-likelihood-confirmation-v1 --resume
python tools/confirm_profile_likelihood.py \
  --directory output/profile-likelihood-confirmation-v1 --audit-only
```

The runner uses one process at nice 19. `--max-new-trials` permits an operational
pause without changing the design; incomplete output cannot pass the decision.
The older physical-calibration series are not resumed by this tool. Detailed
JSON results stay in local output; this document carries the compact outcome.

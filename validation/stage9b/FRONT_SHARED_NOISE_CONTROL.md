# Synthetic control for shared-field front associations

This follows [the training front/balance diagnostic](PHYSICAL_PLATEAU_FRONT_BALANCE.md).
It reads no physical data and fits no correction. Physical acceptance remains incomplete.

## Fixed experiment

32 independent seeded trials each contain three runs of 64 windows on the
560-node bank grid. The true profile is `(1+tanh((x-shift)/2))/2`, with stationary
Gaussian AR(1) shifts (SD 0.08, phi 0.7). Node noise has SD 0.003 and temporal
phi 0.7. An additional independent AR(1) scalar with SD 0.003 is multiplied by
one of three fixed spatial loadings: zero, one, or `tanh((x-4)/2)`.
The last case has opposite signs near the front and on the downstream plateau.
These illustrative laws were not fitted to observed physical covariances.

For each trial, generate two independent noisy measurements of the same true
shifted profiles. Estimate both fronts with the existing isotonic 50% crossing.
Average downstream noise over x in [8,12], subtracting the exact shifted profile
(an oracle operation that isolates measurement noise). Correlate the three-run
mean noise series with the same-field front, independently measured front,
and known true shift. Also record a fixed unshifted-profile residual sensitivity
check and correlation with front estimation error. No per-run correction is fit.

All three cases share random draws within a trial. Report paired same-field
minus independent-field correlations and Monte Carlo standard errors across
independent trials, not a standard error treating windows as independent.

## Result

| Spatial noise case | Same-field front | Independent-field front | True shift | Paired difference (MCSE) |
| --- | ---: | ---: | ---: | ---: |
| Independent nodes | +0.0372 | +0.0358 | +0.0402 | +0.0014 (0.0068) |
| Same-sign common | -0.1231 | +0.0184 | +0.0278 | -0.1415 (0.0085) |
| Opposite-sign common | +0.1713 | +0.0300 | +0.0278 | +0.1412 (0.0093) |

Shared density noise alone can induce a positive or negative association,
depending on spatial covariance. The positive physical group correlations
therefore do not by themselves justify a correction based on front motion.
This does not establish that the physical association is an artifact: the
illustrative covariance is not a validated physical noise law. Finite-trial
true-shift correlations are not identically zero (MCSE about 0.03).

The independent measurement control cannot be implemented by simply splitting
physical windows: autocorrelation and changing true shifts would violate its
construction. Nor is disjoint spatial support automatically independent under
common spatial noise. Do not subtract the synthetic correlation from physical
measurements or select the noise case that best matches the physical result.

## Next decision

Keep the front-dependent correction disabled. Establish whether the stored
physical outputs support a genuinely independent estimator, or whether new
particle-level independent measurement replicates are necessary. Specify that
measurement design and its shared-state/independence checks before another
physical campaign. The synthetic result closes this mechanism check only;
stationarity, confidence coverage and physical acceptance remain unqualified.

## Reproduction and checks

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 tools/control_front_shared_noise.py \
  --output output/front-shared-noise-new.json
python3 -m unittest discover -s tests/python
```

The output must not already exist; the protocol and per-trial statistics are
stored with source hashes and NumPy version. Local evidence:
`output/front-shared-noise-v1.json`, SHA-256
`3618ff004501684af861d8e45ce86ff90004e462798e84779c4511f07b0e0fc1`.
Full recomputation matched exactly. All 215 Python tests passed, including exact
noise-free shifts, opposite crossing responses to spatial common noise,
stationary AR recursion and Monte Carlo standard-error calculations.

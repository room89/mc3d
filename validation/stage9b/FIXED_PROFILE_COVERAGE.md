# Fixed-coordinate ensemble mean: first coverage screen

This implements the baseline decision in [the measurement design](FRONT_MEASUREMENT_DESIGN.md).
The [protocol](FIXED_PROFILE_PROTOCOL.json) was written before execution and is
hashed in the result. No physical data, front alignment or front-dependent
correction are used. Physical acceptance remains incomplete.

## Estimator and experiment

The estimand is the laboratory-frame ensemble mean density, including Gaussian
front motion. Its truth is the Gaussian-shift expectation of the tanh profile,
not the unshifted shock. Gauss-Hermite orders 64 and 128 agree to 5.56e-16.

Each of 256 independent trials has eight independent groups of three runs,
64 temporally correlated windows and 113 fixed spatial nodes. Persistent
run-level and group-level offsets are retained. All three illustrative spatial
noise cases share random draws within a trial. They are not calibrated physical
noise laws. The selected evaluation nodes were fixed at x=0 and x=10.

Average all windows within each run, then average runs within each group.
Estimate the profile by the mean of the eight group profiles. At each node use
`mean +/- 2.3646242515927844 * sd(group means)/sqrt(8)`, the nominal 95% Student
interval with seven degrees of freedom. This is an approximate pointwise
procedure for the nonlinear synthetic profile distribution, not an exact
normality assertion. No fitted quantities or training sample are needed for
this frozen arithmetic estimator; evaluation trials do not tune it.

The report preserves per-trial errors, coverage indicators and half-widths.
Bias MCSE and coverage MCSE use independent trials. The mean-node coverage MCSE
uses each trial's fraction of covered nodes, retaining spatial dependence.

## Results

| Noise case | Mean node coverage (MCSE) | Coverage at x=0 | Coverage at x=10 | Whole-grid coverage |
| --- | ---: | ---: | ---: | ---: |
| Independent nodes | 95.50% (0.91 pp) | 97.27% | 95.70% | 74.61% |
| Same-sign common | 94.58% (1.02 pp) | 96.88% | 94.92% | 75.39% |
| Opposite-sign common | 95.36% (0.87 pp) | 97.27% | 94.92% | 71.48% |

Pointwise coverage is near nominal in this limited screen. Precision at the
selected nodes is about 1.0–1.4 percentage points; whole-grid coverage MCSE is
about 2.7–2.8 points. The pointwise intervals are not a simultaneous 95% band.
The whole-grid result is descriptive, not failure of a promised band guarantee.

Measured bias at x=0 ranges from -0.000197 to -0.000179 (MCSE about 0.000084–
0.000085); at x=10 it ranges from -0.000148 to -0.000139 (MCSE about 0.000054).
Report these finite-trial deviations rather than claiming demonstrated zero
bias. The arithmetic estimator is unbiased under the stipulated zero-mean
noise and stationary shift laws by linearity of expectation; this theoretical
property neither validates the physical laws nor removes finite-screen error.
The cases share draws, so their similar bias estimates are not independent
replications of a systematic effect. No tuning or rerun with a more favorable
seed was performed.

This screen does not qualify physical covariance, intrinsic shock shape,
gradient thickness, simultaneous bands or experimental agreement.

## Next decision

Retain this unaligned estimator as a reference. If acceptance needs a whole
profile uncertainty band, predeclare a group-level simultaneous-band method
and test it independently; do not reuse pointwise intervals as a band. Freeze
the band construction, precision target and independent confirmation seed
before evaluating it. Additional synthetic trials may be needed to resolve
coverage near 95%, but no physical campaign is justified by this screen alone.

## Reproduction and validation

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 tools/screen_fixed_profile_coverage.py \
  --output output/fixed-profile-coverage-new.json
python3 -m unittest discover -s tests/python
```

The output must not exist. Local evidence `output/fixed-profile-coverage-v1.json`,
SHA-256 `4774ecb8c896a4f95f7a177d022b9f056e1cd68e34568b22e6e97ba88576a8ec`.
Full recomputation matched exactly. All 219 Python tests passed. New tests check
shift-averaged truth, persistent group offsets, invariance to duplicated windows,
data translation and invalid grouping. Large generated trial arrays stay local.

## Research PR chain

At the 2026-11-14 check, all three CI jobs passed for each of #69–75; all remain
draft. #69 and #70 target master independently; #71 targets #70, then #72–75
each target their predecessor. Green CI does not establish scientific acceptance.
Review the chain in dependency order and refresh child bases/checks after merges;
no merge or draft-state changes were made as part of this experiment.

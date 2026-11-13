# Training-only plateau mean diagnostic

This research step follows [the physical plateau audit](PHYSICAL_PLATEAUS.md).
Physical acceptance remains incomplete. No solver settings or acceptance gates
change, and no new simulations are required.

## Frozen model and checks

Restore original normalized density from each bank's mean and fluctuations.
Use the training-selected upstream/downstream regions from the preceding audit.
Fit `beta0 + beta1*s + beta2*s^2 + beta3*tau`, where `s` spans each region
and `tau` spans the 64 ordered windows, both from -1 to 1. Window index is
not a physical time law. There are no per-run intercepts, held-out fitting,
profile alignment, clipping, or adaptive choice of model complexity.

Each solver has 24 training runs in eight equal-size groups and 12 held-out
runs in four disjoint groups. Excluding one entire training group at a time
provides residuals that preserve excluded run offsets and a jackknife
coefficient covariance. Regions remain selected from the full training sample:
this is not fully independent region-selection cross-validation.

Held-out mean, spatial-half difference and temporal-half difference are
reported with both conditional group standard errors and descriptive errors
combining held-out group variance with training prediction jackknife variance.
Normal +/-1.96 intervals are diagnostics with only eight/four groups, not
calibrated confidence intervals or simultaneous bands. Residual second moments
and lag-one products retain run offsets and model error; they are not estimates
of a qualified separable noise covariance.

## Historical-bank result

| Solver / region | Held-out node MSE / constant-baseline MSE | Held-out minus training mean |
| --- | ---: | ---: |
| MC3D upstream | 0.946248 | -0.000000864 |
| MC3D downstream | 1.000191 | -0.000047436 |
| SPARTA upstream | 0.951719 | -0.000011954 |
| SPARTA downstream | 0.999651 | -0.000664354 |

The common fitted surface reduces upstream residual MSE by about 5%,
descriptively, while downstream MSE barely changes. Subtracting that common
surface preserves the original held-out/training mean difference exactly.
MC3D downstream temporal half-change remains about +0.000796; SPARTA downstream
spatial half-change remains about +0.000625. The latter's descriptive interval,
including training uncertainty, is [0.000064, 0.001186]. SPARTA downstream mean
transfer has descriptive interval [-0.001133, -0.000196]. These observations
motivate further investigation; they do not establish a solver defect.

The banks are historical engineering runs, previously inspected, and are not
fresh qualification data or a matched production reference campaign. Do not
feed these residuals into the synthetic temporal-covariance estimator as if
its known-mean/stationarity assumptions had been validated. The next research
step is to inspect downstream mean evolution against actual sampling times
and run groups before freezing a richer model and collecting fresh validation.

## Reproduction

From the repository root, with the historical banks available locally:

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 tools/fit_physical_plateau_mean.py \
  --manifest output/shock-physical-fluctuation-banks/manifest.json \
  --output output/physical-plateau-mean-new.json
python3 -m unittest discover -s tests/python
```

The output must not already exist. The report records bank, manifest and
source hashes and verifies evidence remains unchanged during analysis.
Local evidence: `output/physical-plateau-mean-v3.json`, SHA-256
`d5079d15e2c3e94378f9a8b31650ed3c841a15b6ec118b67fdb33250a6929149`.
Recomputation matched the complete saved JSON. All 204 Python tests passed,
including exact model recovery, held-out isolation, excluded-group offsets,
training uncertainty and invalid-input checks. Large banks and generated
prediction arrays are intentionally not committed.

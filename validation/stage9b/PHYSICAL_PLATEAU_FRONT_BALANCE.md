# Training-only front and particle-balance diagnostic

This follows the [time/group audit](PHYSICAL_PLATEAU_TIME.md). Only the 24
training runs per solver are read. No held-out arrays are read, no predictor is
fitted to the diagnostic associations, and physical acceptance remains incomplete.

## Method and scope

Reconstruct density on the bank's fixed spatial grid. Use the unchanged
training-selected downstream region and quadratic-space/linear-time training
mean surface. For each window, compare its downstream residual average with
the isotonic 50% front crossing on the bank grid. Report contemporaneous Pearson
correlations separately for each original three-run group, each run, and the
ensemble curve. Pearson centering is only part of the diagnostic statistic;
the saved density residuals and front positions retain offsets.

For MC3D, verify the historical output manifest and conservation-file hashes,
then check exact integer particle balance and population continuity on every
step. Compute duration-weighted pre-step population in each configured window.
Assign each whole-step net particle flux to a window by its step midpoint.
This avoids inventing fractional particle transfers but is not exact window
boundary flux accounting. The largest step/window ratio is 0.003762. Counts
need no SI conversion; configured windows share the conservation log's time unit.

This step audits particle number only, not momentum or energy. It does not
provide an equivalent SPARTA flux audit. Historical evidence and bank hashes
are checked again after analysis. The full training mean fit is in-sample.

## Observations

| Association | Ensemble correlation | Range across eight group curves | Positive groups |
| --- | ---: | ---: | ---: |
| MC3D residual / front | +0.368 | +0.054 to +0.369 | 8 / 8 |
| SPARTA residual / front | +0.059 | +0.070 to +0.256 | 8 / 8 |
| MC3D residual / population | +0.403 | -0.193 to +0.458 | 6 / 8 |
| MC3D residual / net particle flux | -0.033 | -0.175 to +0.189 | 4 / 8 |
| MC3D front / population | -0.312 | -0.475 to +0.235 | 2 / 8 |

All 630,556 MC3D steps satisfy exact particle-count accounting and continuity.
No particle-count conservation error explains the observed plateau structure.
This is not proof that the boundary fluxes are physically accurate.

The residual/front association has a consistent positive sign across groups,
but is modest. Correlation of ensemble curves is not the average of group
correlations. Front and plateau are estimated from the same density field,
and the series are autocorrelated; no independence, significance, causal
explanation, or deterministic transient is established. Net particle flux
shows no consistent contemporaneous association. Lagged or spatially resolved
transport is not tested here.

## Next decision

Do not add a front-dependent correction yet. First use controlled synthetic
profiles with known shifts and noise to quantify the association induced by
extracting both observables from one field. Compare those with a known-shift
observable; freeze that diagnostic before examining additional physical data.
If the physical association survives that control, evaluate a training-only
front-coupled candidate with whole-group exclusion, then freeze its independent
validation protocol. Fresh physical data will still be required for acceptance.

## Reproduction

From this worktree, supplying the primary checkout that stores the evidence:

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 tools/diagnose_plateau_front_balance.py \
  --root ../.. --output ../../output/physical-plateau-front-new.json
python3 -m unittest discover -s tests/python
```

The output must not already exist. Generated series remain local; the report
records the source, manifest, training-bank and 57 historical input hashes.
Tests cover exact particle accounting, broken continuity/coverage, fractional
step overlap, discrete flux assignment and constant/invalid correlation series.

Local report: `output/physical-plateau-front-v2.json`, SHA-256
`542fc39239b791111153a01774e2987aeb3c17787409c642a3150d3c4f644a7f`.
Complete report recomputation matched exactly; all 211 Python tests passed.

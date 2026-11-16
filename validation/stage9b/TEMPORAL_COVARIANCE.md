# Two-window temporal covariance from independent calibration runs

Physical acceptance remains incomplete. This first noise-estimation step learns
only the total stationary 2×2 temporal covariance. Spatial covariance and plateau
means are supplied. It does not fit separate persistent/AR noise components:
`a=.0001, b=.0004, phi=.8` and `a=.0003, b=.0002, phi=.6` both give lag-zero
covariance .0005 and lag-one covariance .00042. Two windows cannot distinguish
these parameter sets.

`shock_temporal_covariance.estimate` accepts run × two windows × plateau-node
observations, a known plateau mean vector and a positive definite spatial
covariance. It spatially whitens residuals, averages their temporal outer products
across independent runs and whitened nodes, and averages the two diagonal entries
to impose stationarity. It retains persistent run offsets: it never subtracts a
separate mean from each run. Degenerate estimates are rejected without a ridge.
The output includes the unconstrained sample matrix and the stationary matrix.

## Synthetic diagnostic

The unchanged piecewise-linear generator from the profile screen supplies 32
independent calibration runs and 32 separate evaluation runs per trial. Seed
families 202611160000 (calibration) and 202611150000 (evaluation), each combined
with trial index, keep them independent. All four motion/noise combinations for
logistic and Gompertz are retained, with 16 repetitions each: 128 trials total.

Calibration uses nodes `abs(x) >= 8 + motion_amplitude`: these are exact constant
plateaus for every allowed shift in this **specified synthetic model**. Their
means are zero and one. The spatial submatrix used for whitening is known.
Estimated covariance is applied to separate evaluation runs; the same marginal
profile fitter with true covariance is the paired control. No confidence bands,
physical-data calibration or advancement criterion are asserted by this screen.

All 128 trials completed with both profile fits converged. Mean relative
Frobenius covariance error was 0.0717852. Mean profile ISE was 0.00113403 with
estimated covariance versus 0.00113722 with known covariance. Paired difference
was -0.00000319, MCSE 0.00000595; this does not establish superiority or formal
noninferiority. This diagnostic supplies feasibility evidence only.

The audit checks all checkpoint hashes, regenerates calibration inputs, verifies
the estimated matrices and profile error arithmetic, and recomputes the summary.
It does not rerun the expensive profile optimization. Plan SHA-256:
`e4eb3c3b2a179153db16ae69a4aaf829865a8b821a4a98da1e0caa0d996b52e4`.
Report SHA-256:
`77cc25e41349b43ca01ad44487af090b3543e996a9dc02e7906b99017a4e8ec0`.

## Reproduction and limits

Set `OPENBLAS_NUM_THREADS=1`, `VECLIB_MAXIMUM_THREADS=1`, `OMP_NUM_THREADS=1`:

```sh
python tools/screen_temporal_covariance.py --directory output/temporal-covariance-screen-v1
python tools/audit_temporal_covariance.py --directory output/temporal-covariance-screen-v1
```

Use `--resume` for an interrupted screen with the same frozen sources. The screen
runs sequentially at nice 19 and stops on invalid covariance inputs. Existing
reports must match recomputed summaries; completed files are never overwritten.
Raw evidence remains under local output. No existing acceptance runner is changed.

Before physical use, plateau selection and means, spatial covariance, temporal
stationarity, and separation of calibration/evaluation runs require validation
on actual observations. Exact plateau assumptions here must not silently be
transferred to the tails of analytic curves or physical shock profiles. Estimated
covariance uncertainty must eventually enter profile confidence bands.

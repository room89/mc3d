# Stage 9b: remaining acceptance work

Physical acceptance is **not complete**. Passing software tests does not qualify
an estimator or validate the physical model. Existing acceptance gates remain unchanged.

- The current order-envelope candidate fails density coverage in conditional
  calibration cases. Increasing the sample count alone is not an established fix.
- Historical MC3D/SPARTA comparisons used a different domain from the matched
  experimental reference; their apparent agreement is not a final acceptance result.
- The matched 40-experimental-mean-free-path matrix is prepared, but its production
  run is blocked on estimator qualification.
- Global shift search and latent-shift likelihood are separate research candidates.
  The likelihood prototype assumes known covariance; it does not yet fit and qualify
  the uncertainty model.

## Next steps

1. Finish and audit the ongoing coverage diagnostics with the agreed reduced load.
2. Resolve estimator bias/coverage failures and qualify on independent scenarios.
3. Run the matched reference matrix with the qualified estimator.
4. Complete convergence, stationarity, cross-code and experimental uncertainty
   checks before declaring acceptance. Do not relax gates to obtain a pass.

## Evidence and provenance

The full investigation, including negative results and intermediate reports, is
preserved in the [original PR snapshot](https://github.com/room89/mc3d/tree/2c58c95701e45c7e3b196e7e86c441fc8496b594/validation/stage9b).
See [archive and file disposition](PR59_ARCHIVE.md) for checksums and reproduction.
This compact status deliberately records no live progress counters: local ongoing
runs must be inspected before reporting their completion.

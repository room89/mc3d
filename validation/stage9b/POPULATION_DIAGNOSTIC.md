# Population sensitivity control

The completed domain diagnostic did not establish stationary fixed-coordinate
profiles. Test particle-number sensitivity independently using variant E,
longitudinal domain multiplier 2, warmup 20 original-domain transit times,
64 equal sampling windows over [20,28], and population multiplier 4.
This changes 160 to 640 model particles per upstream cell and divides model
particle weight by four. Grid, reservoirs, physical domain, gas model and
physical time schedule remain unchanged. Use the same pinned MC3D and SPARTA
binaries as DOMAIN_DIAGNOSTIC.md, seeds 210001, 210002 and 210003, with a new
output directory. Same seed labels do not imply paired random trajectories.

Compare with the completed multiplier-1 doubled-domain campaign. Retain all
runs and windows, with no automatic retry or favorable-window selection.
Report original-coordinate mean density differences, run-level front and
width descriptions, half-run changes and full-domain density inventories.
Do not fit alignment as an acceptance correction. A fourfold particle count
suggests a factor-two noise reduction only under an independent-particle
approximation; it is not a promised variance reduction or a pass criterion.
Three seeds provide exploratory evidence, not calibrated uncertainty.

The runner option is `--population-multiplier 4`. Its default is 1, and changing
it on resume is rejected before existing evidence is rewritten. This control
does not itself close numerical convergence, experimental comparison or
stages 9–11.

## Completed exploratory result

All six runs completed in `output/shock-population4-domain2-v1`.
The analysis `output/shock-population4-analysis-v1.json` has SHA-256
`02d5c090da09c3195991584dcb15a904340b74339e2073b1be5c9c63ae434b13`.
It records 1,621 input hashes and retains all three seeds and 64 windows.

Maximum absolute changes below use density normalized by the reservoir
jump, in original coordinates, without profile alignment:

| Quantity | Population 1 | Population 4 |
| --- | ---: | ---: |
| MC3D ensemble second-half minus first-half | 0.01136045 | 0.00388104 |
| SPARTA ensemble second-half minus first-half | 0.02596195 | 0.00562813 |
| MC3D versus SPARTA ensemble mean | 0.09970725 | 0.04910187 |

The maximum population-4 minus population-1 mean change is 0.03592533
for MC3D and 0.01654682 for SPARTA. Reduced half-run changes do not establish
stationarity or population convergence from two populations and three seeds.

An exploratory translation diagnostic locates ensemble fronts at 0.14493149
and 0.34281224 collision mean free paths for MC3D and SPARTA. Their 10–90
widths are 4.64562758 and 4.63499274. Fitting a shift of 0.19788075 to these
same profiles reduces the maximum difference on [-11.975,11.975] to
0.00327738. This is a description of the discrepancy, not an acceptance
correction or an independently validated shape comparison.

Full-domain integrated normalized density decreases between the two halves
in all six runs. The changes are [-0.01240802,-0.02963447,-0.02015064] for
MC3D and [-0.03222984,-0.02381037,-0.03283226] for SPARTA, in seed order.
Changes outside the central region also contribute; a rigid front translation
alone is not established. No significance claim follows from these signs.

MC3D particle-count balances close exactly over each complete run and eight
fixed time sections; momentum and energy residuals are below 1e-10 times
the sum of absolute incoming and outgoing fluxes. This verifies bookkeeping,
not the physical correctness of outflow. Initial particle-count discrepancies
correspond to less than 0.00039 mean free paths under a rigid-step mass model;
that simple initial count error cannot explain the observed front positions.
Initial momentum and energy fluctuations are not thereby excluded.

Supporting local artifacts are `output/population4-translation-diagnostic-v1.json`,
`output/population4-full-regions-v1.json`,
`output/population4-boundary-exchange-v1.json`, and
`output/population4-initial-inventory-v1.json`.
Physical acceptance remains incomplete. The remaining fixed-coordinate
cross-code discrepancy requires a prospective uncertainty and convergence
check; fitting shifts or selecting favorable windows must not turn it into a pass.

## Next isolated control: inflow count rounding

Pinned MC3D source `398078186333600d9a704ce0c4da447468eac8a6` carries the
fractional expected inflow across steps. Pinned SPARTA source
`95b9abaa8bd548991cc3c3f1c58b34722f7ade74` instead uses `floor(expected + U)`
with a new uniform variate per step. Their collision/motion order also differs.
Neither difference alone identifies a defect. The next control will change only
MC3D hyperfree inflow count rounding; collision/motion order stays fixed.

Implement an explicit experimental configuration mode with accumulated rounding
as the unchanged default. Stochastic rounding must preserve the analytic expected
incoming count and use the seeded cell random stream. Record the mode in effective
configuration and prevent resume from silently mixing modes. Do not apply this
mode to the legacy free-boundary batch temperature normalization.

Before production, check counts below one and above one per step, varying time
steps, reproducibility, invalid modes, serialization, conservation, and default
behavior. Verify the new binary against the pinned control with identical
physical inputs. Freeze both the new source and binary hashes.

Run the E, domain-2, population-4, warmup-20 control for three seeds 210001–210003
in a new directory, retaining the complete [20,28] interval and 64 windows.
Use the completed accumulated-rounding MC3D and SPARTA runs as historical
comparators. Reused seed labels do not make trajectories paired. A changed random
stream is part of this exploratory intervention, so three runs cannot certify
causality or acceptance. No automatic retries, selected windows or fitted-shift
acceptance correction are permitted.

Report original-coordinate profile means, front and width descriptions, half-run
changes, full-domain inventory and boundary balances regardless of outcome.
A smaller discrepancy alone must not promote the experimental mode to default.
If discrepancy persists, retain the negative result and isolate time-step/order
sensitivity next. Any promising result still requires fresh-seed confirmation,
qualified uncertainty and the remaining numerical and experimental gates.

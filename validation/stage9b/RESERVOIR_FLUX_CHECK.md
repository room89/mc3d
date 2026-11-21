# MC3D engineering-shock reservoir flux check

The historical fixed-coordinate audit found a mean-front offset between codes.
Before attributing it to boundary sampling, check the actual MC3D hyperfree
boundary against analytical flux moments at both M=3.8 engineering reservoirs.

`HyperFreeSamplingTest.ShockReservoirMomentumAndEnergyFluxes` uses the saved
normalized upstream and downstream temperatures, drift speeds and density ratios.
It calls `HyperFreeBoundary::GenerateInflow` for 1,000 alternating unequal steps,
with 199,750 expected particles per reservoir. There are no collisions or walls;
large test step sizes exercise sampling and count accumulation, not transport
accuracy or production stability. No solver implementation is changed.

Let M_k be the integral of w^k times the normal Gaussian density over w>0.
The test obtains M_0 and M_1 analytically and uses Gaussian integration by parts
for higher moments. Incoming mean normal momentum is M_2/M_1; mean kinetic
energy is (M_3/M_1 + 2T)/2. Analytical sampling variances use M_3 and M_5,
including both independent tangential Gaussian components. These are flux-weighted
moments, not the unweighted reservoir velocity moments.

Checks require cumulative count discrepancy at most one particle and momentum
and energy means within six analytical standard errors, with a fixed seed chosen
before execution. Both reservoirs passed, along with the existing hyperfree,
open-boundary and statistical-moment tests (22 tests total).

This establishes a local sampling regression check for MC3D only. It does not
validate the SPARTA boundary, full shock transport, stationarity, physical
confidence bands, or experimental agreement. The corresponding SPARTA measurement is recorded below; time-discretization
and full-shock sensitivity remain open.

Reproduction after building `mc3d_tests`:

```sh
build/mc3d_tests --gtest_filter=HyperFreeSamplingTest.ShockReservoirMomentumAndEnergyFluxes
```

## Pinned SPARTA measurement

Executed the historical SPARTA binary (SHA-256
`0d6dedda6d90b1f703d5cd0f745dbeb976cac7f74bf57e9a13016b26e8b10afb`)
for one step per reservoir, starting empty with no collisions or walls.
The box is 2×1×1 m, with periodic transverse faces and injection only from
the selected x face. Reservoir states and time step come from the historical
C reference. Particle weight is chosen for 200,000.25 expected arrivals; this
is an isolated measurement configuration, not the historical shock population.
All particles survive this short step, so the particle dump measures arrivals.

| Reservoir | Actual count | Momentum deviation / analytic SE | Energy deviation / analytic SE |
| --- | ---: | ---: | ---: |
| Upstream | 200000 | -0.831 | -0.714 |
| Downstream | 200000 | -0.234 | -0.347 |

Both cases meet the preselected one-particle count and six-standard-error
moment checks. No gross boundary sign, normalization or moment discrepancy
appears at this precision. The known small normal-tail truncation is not
excluded by these finite samples; passing is not a claim of exact equality
to the untruncated Maxwellian law or equivalent count fluctuations.

Local script: `output/check_sparta_reservoir_flux.py`; run from the primary
checkout using its Python environment. It refuses to overwrite the output
folder. Decks, particle dumps and logs are in `output/sparta-reservoir-flux-v1`.
The report records binary, reference, species, script, input, dump and log hashes.
Report SHA-256: `b03fccff63cdcae789f73426d007a5cac5c31a3c10b9759ab2cebe87f23c65a8`.

With both isolated flux checks passing, do not change reservoir parameters to
remove the historical front offset. Next compare historical C and reduced-time-step
E profiles in fixed coordinates to quantify sensitivity before a longer run.

## Historical C/E time-step sensitivity

Re-extracted all 64 windows for seeds 123, 456 and 789 in C and E for both
codes (12 runs). E halves MC3D's transport/collision limits and SPARTA's fixed
time step. Recover the original source coordinates by adding each window's
stored scalar front offset back to `x_aligned`; the density arrays themselves
are not shifted or interpolated. Compare the 560 original nodes in [-14,14] λ₁.
All 834 inspected file hashes were unchanged after reading; MC3D window files
also match their saved output manifests, and the existing parser verifies the
per-step conservation balances.

| Comparison | Maximum absolute mean density difference | Mean-front difference, λ₁ |
| --- | ---: | ---: |
| MC3D C − SPARTA C | 0.004931 | -0.016654 |
| MC3D E − SPARTA E | 0.010100 | -0.039574 |
| MC3D C − MC3D E | 0.015969 | +0.057075 |
| SPARTA C − SPARTA E | 0.010172 | +0.034155 |

Reducing the step does not monotonically improve the observed cross-code mean
agreement in this three-seed subset. The C result here differs from the 24-run
training audit because it uses only the original three seeds; it cannot replace
that larger-sample result. There is no simultaneous uncertainty estimate here,
so these differences neither establish time-discretization bias nor prove
convergence. Same seed numbers across cases do not provide common random paths.
No zero-step extrapolation or acceptance pass is inferred.

Local read-only script: `output/audit_fixed_timestep.py`; report:
`output/fixed-timestep-audit-v1.json`, SHA-256 `4909cd75e488bfc8a9834cbf6bece1d612c8ddf604afc7661dcf4c5ee65e6ef1`.
Run from the primary checkout using its Python environment; output must not exist.
The report retains all mean profiles, per-seed front locations, half-profile
changes, input hashes and script hash. The next design must estimate independent
run uncertainty for C/E comparisons and temporal changes; three-seed point
estimates do not select a production time step reliably.

## SPARTA initial-state check

A separate run-zero measurement executes the exact initialization prefix of the
saved `sparta-C-123/in.shock` deck with the pinned binary, before any boundary
injection, transport or collision. The particle dump contains 64,000 upstream
particles (64,000 expected) and 211,963 downstream particles (211,963.303
expected). Regional mean velocity components differ from their prescribed
values by at most 2.801 analytical standard errors; thermal energy deviations
are -0.351 and +0.146 standard errors. Both satisfy the preselected six-SE check.

This checks the realized regional population and first/second velocity moments
for one seed. It does not establish equivalent cell-population covariance or
prove that initial transients have decayed. MC3D's initial-step implementation
uses cell-wise stochastic count rounding and unnormalized Gaussian velocities;
SPARTA distributes a regional particle total across cells. These need not yield
identical finite-particle fluctuations.

Local script `output/check_sparta_initial_state.py` creates
`output/sparta-initial-state-v1` from the primary checkout and refuses overwrite.
The report records source/reference/deck/binary/dump hashes and observed values.
Report SHA-256: `4cd9aa658e9ee1f022a9415f358ee0c81e48c4fedf4b417433337c096f948a7f`. No production configuration was modified.

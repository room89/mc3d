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
confidence bands, or experimental agreement. The next investigation remains a
corresponding isolated SPARTA flux measurement and time-discretization check.

Reproduction after building `mc3d_tests`:

```sh
build/mc3d_tests --gtest_filter=HyperFreeSamplingTest.ShockReservoirMomentumAndEnergyFluxes
```

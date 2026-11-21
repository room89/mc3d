# Historical fixed-coordinate precision audit

After merging #70–77, apply the unaligned group-mean estimator to the 24
historical training runs per solver. Restore raw profiles from the verified
bank mean plus fluctuations. These banks retain original fixed physical nodes;
there is no moving interpolation, per-run centering or fitted correction.
Held-out banks are not read. This is planning evidence, not physical acceptance.

Keep all 560 shared source nodes rather than interpolate onto the synthetic
113-node grid. Average windows within runs and three runs within each of eight
original groups. Allocate nominal family error 0.025 to each solver over all
560 nodes using t7 Bonferroni intervals. Sum the two half-widths for the
MC3D-minus-SPARTA difference. This construction needs no cross-code independence
or seed pairing, but still requires valid marginal Student intervals. Synthetic
coverage on the earlier grid does not establish physical coverage here.

## Result

- Maximum absolute difference in normalized density: **0.017975**.
- Maximum diagnostic difference half-width: **0.061293**.
- Maximum pointwise sum of absolute difference and half-width: **0.079145**.
- **356 / 560** nodes exceed the existing precision target **0.005**.

The earlier density tolerance is 0.01. Even the maximum mean difference exceeds
it in these laboratory coordinates. Increasing the number of seeds alone is
therefore not an evidence-based remedy. First inspect differences in mean front
location, time evolution and boundary setup without fitting away the mismatch.
These observations do not establish a solver defect: a fixed-coordinate mean
and a centered intrinsic shock shape are different estimands.

## Location and time diagnostic

The isotonic 50% crossing of the full mean profile is 0.130238 λ₁ for
MC3D and 0.204730 λ₁ for SPARTA, a separation of 0.074492 λ₁. Crossings
are descriptive: no translation is applied to the reported differences or bands.
The eight group-mean crossings range from 0.09722 to 0.18003 for MC3D and
0.13684 to 0.28349 for SPARTA.

Between the first and second 32-window halves, the ensemble crossing moves
from 0.121069 to 0.139093 for MC3D and from 0.184394 to 0.225237 for SPARTA.
Maximum absolute half-profile changes are 0.005298 and 0.010342 respectively.
These are observed changes, not significance tests or proof of systematic drift.
The report retains all eight paired group half-profile changes; windows are
not independent replicates. A stationary noise realization can also show a
half difference. No stationarity pass/fail is inferred from these numbers.

The mean front separation makes a location contribution plausible, but does
not quantify its causal share or establish equal intrinsic shapes. Before a
new campaign, audit identical physical initial/boundary states and the time
normalization in both solvers. Any shape-only translation diagnostic must stay
separate from the laboratory-frame acceptance result.

## Historical input equivalence

A separate read-only audit checked all 24 training run pairs, verifying 129
historical files against recorded hashes before and after reading. It compared
MC3D normalized inputs converted to SI and its effective run manifest against
SPARTA's saved reference and actual input deck. Matching quantities were:

- Reservoir density, temperature and all velocity components; initial left/right
  states and discontinuity position.
- Domain dimensions, cell counts and model-particle weight.
- Periodic transverse boundaries and declared two-sided reservoir injection.
- All 64 physical sampling-window endpoints, warm-up and production step counts.

Sampling spans 0.0019540413973 to 0.0058621241919 seconds in both codes.
No mismatch in these declared inputs explains the observed mean-front offset.
This is not proof of equivalent boundary algorithms or realized initial samples.
MC3D adaptive integration and SPARTA fixed-step averaging also have different
quadrature; equal endpoints alone do not establish equal temporal sampling.
Collision kernels and experimental conditions are outside this input audit.

Local reproduction (from primary checkout, output must not exist):
`python3 output/audit_historical_setup.py`. The one-off script and report remain
local, like the earlier storage audit; the report includes all inspected hashes.
Report: `output/historical-setup-audit-v1.json`, SHA-256
`5a73621562b6255fa507a79ca501dca5336735a469699cc0984023052153a073`.
Next inspect boundary sampling and initialization in the pinned implementations,
and quantify temporal-discretization sensitivity before choosing production.

## Remaining acceptance requirements

The full scope remains PHYSICS_PLAN.md stages 9–11. The synthetic band is only
one component. Stage 9 still needs stationarity, grid/time/population/domain
convergence, inter-code agreement and profile/thickness comparison with the
experimental primary source. Gradient thickness is not qualified by a density
band. The isothermal wedge and dynamic wall campaigns remain subsequent stages;
this diagnostic does not replace them or change their tolerances.

Before launching a fresh campaign, resolve the laboratory-frame mismatch and
freeze the intended estimand, physical setup, independent seeds, sampling
schedule, precision and stopping rules. Do not infer campaign size by scaling
this historical variance alone. Preserve all prior incomplete verdicts.

## Reproduction

From this worktree with evidence in the primary checkout:

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 tools/audit_fixed_profile_precision.py \
  --root ../.. --output ../../output/fixed-profile-physical-precision-new.json
python3 -m unittest discover -s tests/python
```

The output must not exist. Local evidence:
`output/fixed-profile-physical-precision-v2.json`. The report records manifest,
bank and audit source hashes; inputs are checked again after analysis. It retains
all group means, means, widths and differences for inspection. No physical solver
run or model change is part of this audit.

## Pinned boundary implementation findings

Source inspection used MC3D revision
`398078186333600d9a704ce0c4da447468eac8a6` recorded by the historical run
manifest, and the local SPARTA source checkout at
`95b9abaa8bd548991cc3c3f1c58b34722f7ade74`.

MC3D `HyperFreeBoundary::GenerateInflow` calls `AccumulateInflow` with minimum
batch one. The fractional expected count is carried between steps. SPARTA
`FixEmitFace::perform_task`, with the deck's default count mode, adds a uniform
random number and truncates each step's expected count. For a prescribed sequence
of expected counts, MC3D's cumulative rounding discrepancy stays below one
particle per face cell; SPARTA's increments have rounding variance f(1-f),
where f is the fractional expected count. Neither statement establishes equal
finite-particle fluctuations, even with equal mean prescribed reservoir flux.
Adaptive steps and the evolving system prevent using this observation alone as
a quantitative explanation of front variance or mean bias.

Both implementations sample an entry time uniformly within the step and advance
only the remaining fraction. MC3D's hyperfree velocity sampler uses an unbounded
flux-weighted normal distribution plus Gaussian tangential velocities. The
pinned SPARTA rejection sampler restricts the normalized normal thermal deviate
beta to [-3,3]. Its tangential draw is not subject to this cutoff. For inward
drift s in units sqrt(2 kT/m), the ideal incoming density is proportional to
(beta+s) exp(-beta²), beta >= -s. Integrating this law gives excluded normal
flux probability 0.0000221085 upstream and 0.000266525 downstream for these
reservoirs. These are tail probabilities, not observed density or shock errors;
SPARTA still samples its prescribed number of incoming particles from the
truncated velocity law. Energy and momentum effects require separate evaluation.

No solver change is justified solely by this inspection. Next use an isolated
reservoir test to measure count, momentum and energy injection against analytic
flux moments under the actual boundary implementations, then assess time-step
sensitivity. Do not replace MC3D's rounding or alter SPARTA merely to make their
fluctuations agree. The initial-step MC3D sampler uses independent Gaussian
velocities without batch temperature normalization; SPARTA initialization still
needs a corresponding distribution and population audit.

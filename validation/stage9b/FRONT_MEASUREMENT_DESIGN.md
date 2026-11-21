# Independent front measurement: storage audit and decision

This follows the [shared-noise control](FRONT_SHARED_NOISE_CONTROL.md).
**Existing stored outputs do not establish an independent front measurement.**
Keep the front-dependent plateau correction disabled. This decision does not
complete physical acceptance and does not authorize a new solver campaign.

## What is actually stored

Read-only audit of all 24 training runs for each solver:

| Output | Observed content | Independence limitation |
| --- | --- | --- |
| MC3D: 64 cell-statistics windows per run | Integrated particle count and velocity moments, cell volume, derived density/temperature/stress/heat flux | No per-particle contributions, cohort labels or independent measurement replicates |
| MC3D: 60 timed `.dat` snapshots per run | `x;y;z;N;ro;T;vx;vy;vz;E`, one row per cell | These are cell snapshots, not particle lists; sampling times differ from integrated windows |
| MC3D conservation log | Whole-domain stocks and boundary transfers | Same evolving particle system; not an independent local front sensor |
| SPARTA production profiles | `dump ... grid`, window-averaged density and velocity/thermal moments | Grid aggregates of the same evolving population, not independent particle measurements |

The audit verified 202 input/output hashes before and after inspection. It used
the hashed output inventories, sampled MC3D window 0/63 and lexicographically
first/last timed snapshot headers in every run, and inspected each SPARTA dump
declaration and first production profile. It does not claim to parse every
historical snapshot. No held-out arrays or physical profile fits were used.

The MC3D source agrees with the sampled format: `src/simulation_runner.cpp`
writes cell snapshots and `src/statistical_diagnostics.cpp` writes integrated
moments. `src/particle.h` stores position and velocity without a stable particle
identity or measurement-cohort field in this revision. Source inspection is
supporting evidence; historical conclusions rely on the hashed outputs.

Temperature or velocity from the same cells may be useful alternative
observables, but their shared particles and moment calculations preclude an
independence claim. Disjoint spatial nodes, alternating windows and independent
random weights applied to the same aggregate are likewise not independent
realizations of the physical density noise.

## Two different uncertainty targets

1. **Measurement randomization conditional on a fixed particle state.** Two
   independent sampling/weight streams can have independent randomization
   errors conditional on that state. Both still share the finite-particle
   realization, its fluctuations and dynamics. This only tests an additional
   measurement layer; it does not reproduce the independent physical-noise
   control in the synthetic experiment.
2. **Finite-particle/dynamical uncertainty.** Independently seeded trajectories
   can estimate the distribution of the complete estimator, including front
   motion. They do not share exactly the same latent front history. They support
   ensemble uncertainty evaluation, not a paired same-front subtraction of noise.

A single random partition into two cohorts also has complementary sampling
constraints. Collision interactions can create dependence between tagged
cohorts. Neither tagging nor raw particle dumps alone establishes the required
independence. A physical independence claim needs its own validated model.

## Chosen next step

Do not add particle dumps or cohort instrumentation just to assert independence.
First build a small estimator-qualification experiment that retains the coupled
front/plateau noise. Evaluate the **complete uncorrected estimator**, including
its front extraction, mean fitting and run offsets, against known truth in
synthetic independent runs. Keep comparisons paired within each synthetic trial
and train any fitted quantities on separate groups.

Freeze the estimand before implementation: the laboratory-frame ensemble mean
profile and a translation-invariant shape of a latent shock are different
objects when the front moves. For the next baseline use the **fixed-coordinate
ensemble mean density**, with no run/window alignment and no front-dependent
correction. Define its known synthetic truth by averaging over the prescribed
shift distribution, not by using the unshifted profile. This baseline can
quantify coupled-noise bias and uncertainty without pretending that front and
plateau errors are independent. It does not replace eventual shock-shape or
experimental validation.

Freeze the first experiment's shift/noise laws, number of independent training
and evaluation groups, estimator, uncertainty procedure and reported diagnostics
before running it. Include common noise of both signs and persistent run-level
offsets. Report Monte Carlo precision for bias and coverage; do not tune the
procedure on evaluation outcomes or promote an inconclusive screen to acceptance.
No new physical campaign is needed for that step.

## If additional solver measurements become necessary

A separate instrumentation proposal must specify: stable particle identities
and explicit inlet/deletion rules; assignment/weight RNG streams independent of
collision and injection RNG; exact sampling times and weights; cohort-resolved
zeroth moments with population/volume normalization; and manifest provenance
for the measurement design. Prefer aggregate cohort moments to large particle
trajectories when they suffice for the stated estimand.

Required checks include unchanged solver trajectories when measurement is
turned on, reconstruction of unsplit totals, controlled-state tests of
conditional randomization errors, and explicit tests for dependence remaining
through shared states/collisions. Passing these checks still does not by itself
qualify independent finite-particle noise. Do not launch a costly campaign
until its measurement target and independence assumptions are reviewable.

## Local evidence

Audit script: `output/audit_front_measurement_storage.py` (one-off read-only
inspection, not part of the production toolchain). Run from the primary checkout:

```sh
python3 output/audit_front_measurement_storage.py
```

The script refuses to overwrite its report. Report:
`output/front-measurement-storage-v1.json`, SHA-256
`8eb91abd1b839c035a19abeaa7e5a8680cff5ec87adb8e7a3634bbf3071934fe`.
The report includes the audit script hash, all inspected input hashes, sampled
headers and per-run inventories. No solver source or numerical thresholds change
in this documentation-only step; existing tests are not rerun for prose edits.

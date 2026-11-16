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
`output/fixed-profile-physical-precision-v1.json`. The report records manifest,
bank and audit source hashes; inputs are checked again after analysis. It retains
all group means, means, widths and differences for inspection. No physical solver
run or model change is part of this audit.

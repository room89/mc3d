# Independent C/E precision pilot

Freeze TIMESTEP_PILOT_PROTOCOL.json before execution. The 12 new calculations
use the historical engineering-shock physical setup and pinned binaries, with
three new seeds and the existing C/E time-step limits. This is a precision and
cost pilot, not an experimental-condition campaign or acceptance sample.

Historical C/E comparisons have only three seeds, and the larger historical C
sample differs materially from that subset. Both solvers pass isolated reservoir
moment checks at the tested precision. New C/E data will measure reproducibility
and between-run variability before selecting a larger independent comparison.
No parameter tuning or profile alignment is permitted. Same seed numbers across
codes/cases do not justify pairing their random paths.

The 0.01 density tolerance and 0.005 precision target are retained. Three seeds
cannot qualify the physical simultaneous band or resolve all convergence gates.
Report run-level means and half changes, not windows as independent replicates.
Any sample-size projection must state its distributional assumptions and sensitivity;
never declare sufficient precision solely by scaling historical variance.

Two variant runners execute concurrently. Historical sequential cost for all
12 cases was about 52 minutes, excluding preflight and analysis; concurrency may
reduce elapsed time but contention can offset the gain. Actual timing is recorded.
No automatic solver restart or expansion beyond these 12 runs is authorized by
this protocol. The user's standing authorization permits starting this prepared
pilot without another confirmation.

Evidence directory: `output/shock-timestep-pilot-v2`, with a copied protocol,
frozen tools, binary hashes, logs and independent C/E subdirectories. Existing
`shock_validation.py` records execution separately from physical gates. Its legacy
analysis is retained for traceability and cannot promote this pilot to acceptance.

Reproduction uses the copied runtime tools, once for each variant:

```sh
python3 runtime-tools/shock_validation.py --binary PATH_TO_PINNED_MC3D \
  --sparta-bin PATH_TO_PINNED_SPARTA --output-dir C --variants C \
  --seeds 210001 210002 210003 --load-mode full --no-time-budget
```

Replace C with E for the second runner. Paths are supplied from the campaign
manifest; the binary hashes must match the frozen protocol. The next stage
requires stationarity, all numerical refinements, shock thickness and primary
experimental comparison; none is waived by completing this pilot.

The first attempt in `output/shock-timestep-pilot-v1` terminated with return
code 3 in both runners because the sandbox denied `ps` for RSS monitoring.
After inspecting the terminal exit records, a new directory was launched with
process-monitoring access. The original evidence is retained; protocol, seeds,
binaries and physical parameters are unchanged. This is an explicit infrastructure
retry, not an automatic solver restart or selection based on physical results.

## Analysis after both runners terminate

```sh
OPENBLAS_NUM_THREADS=1 python3 tools/analyze_timestep_pilot.py \
  ../../output/shock-timestep-pilot-v2 \
  --output ../../output/shock-timestep-pilot-analysis-v1.json
```

The analyzer rejects an active runner, an execution error, incomplete cases,
changed binary identities or changed saved SPARTA/MC3D outputs. It verifies
runtime source hashes and snapshots campaign file hashes before reading physical
profiles, then checks them again after analysis. Output must be outside the
campaign and must not already exist. It retains each run's mean and half change,
between-run SD and the four fixed-coordinate comparisons. It does not produce
an acceptance verdict or infer a confidence band from three seeds.

## Completed pilot: results

All 12 planned calculations completed successfully. Both variant runners exit
with code 2 because physical acceptance remains incomplete, not because a solver
failed. The independent analyzer completed and verified 2,187 file hashes before
and after reading. No runs or windows were excluded.

| Comparison | Maximum absolute normalized-density difference |
| --- | ---: |
| MC3D C − SPARTA C | 0.0331965 |
| MC3D E − SPARTA E | 0.0285119 |
| MC3D C − MC3D E | 0.0084791 |
| SPARTA C − SPARTA E | 0.0061330 |

| Case | Maximum between-run density SD | Range of per-run mean-front positions, λ₁ |
| --- | ---: | ---: |
| MC3D C | 0.0379359 | 0.020961–0.309820 |
| SPARTA C | 0.0053066 | 0.263956–0.291803 |
| MC3D E | 0.0124070 | 0.121194–0.211112 |
| SPARTA E | 0.0050926 | 0.289242–0.306958 |

The inter-code mean discrepancies exceed 0.01 in both cases. This does not
establish a statistically significant solver bias from three seeds; it does
show that the earlier three-seed apparent agreement is not reproducible at the
required precision. Do not select favorable seeds or align away this difference.
Between-run SD is not a simultaneous uncertainty bound. Half-profile changes
remain appreciable and are retained per run in the report.

Next diagnose the relation between initial conserved quantities, front movement
and run offsets, and determine whether the laboratory-frame transient target
can support the intended physical comparison. A larger production sample or a
change of estimand needs its own justification; reducing the step alone did not
remove the observed mismatch. Thickness and experimental acceptance remain open.

C's runner elapsed 23.44 wall minutes; E's elapsed 84.28 wall minutes. Recorded
run/preflight durations were 1,352.41 and 2,284.67 seconds respectively. These
clocks disagree materially for E; do not treat wall time as CPU cost or infer
a cause such as machine sleep without further evidence.

Local report: `output/shock-timestep-pilot-analysis-v1.json`, SHA-256
`71f1182b45968c0839e8955230937cbdf1816623603b6a10d81e0c2c81dc2371`. Large per-run arrays and raw outputs remain local.

A second complete extraction from the original campaign using the worktree
analyzer matched the full report JSON exactly (`output/shock-timestep-pilot-analysis-replay-v1.json`).

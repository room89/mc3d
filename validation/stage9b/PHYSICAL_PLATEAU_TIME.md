# Physical time and downstream group diagnostic

This read-only follow-up to [the mean diagnostic](PHYSICAL_PLATEAU_MEAN.md)
checks sampling schedules and decomposes downstream averages by original run
group. No new model is fitted and physical acceptance remains incomplete.

## Evidence and method

All 72 runs have 64 windows spanning 0.0019540413973 to 0.0058621241919 seconds.
MC3D uses configured integration bounds converted with its SI time scale;
these bounds agree with the completed output manifest. SPARTA uses verified
averaging cycles (258 steps per window), timestep and dump-step headers.
SPARTA averages discrete step-end samples, so matching cycle bounds does not
prove identical continuous-time sampling kernels between solvers.

The tool verifies 2,497 historical files against hashes linked to the bank
manifest, including all SPARTA dump files, and checks them again after reading.
Bank and training-model provenance are retained. MC3D profile values are taken
from the previously verified banks; this is not a fresh extraction of every
MC3D cell file. No historical evidence is rewritten.

Window midpoints are affine in window number (normalized maximum discrepancy
4.45e-16). Replacing window index with SI time in the existing linear term
therefore cannot improve that fit. Eight consecutive blocks of eight windows
and original three-run groups are reported without recentering individual runs.
Residuals subtract the fixed full-training surface from the preceding step;
training residuals here are in-sample. Removing one group from the *summary*
measures its influence; it does not refit the surface.

## Descriptive result

Values below are last-half minus first-half normalized density after subtracting
the training surface. The ranges exclude each group in turn from the average.

| Solver / split | Half-change | Positive groups | Leave-one-group-out average range |
| --- | ---: | ---: | ---: |
| MC3D training | +0.000717 | 7 / 8 | +0.000470 to +0.000909 |
| MC3D held-out | +0.000796 | 3 / 4 | +0.000079 to +0.001226 |
| SPARTA training | -0.000634 | 3 / 8 | -0.000867 to -0.000363 |
| SPARTA held-out | -0.000009 | 2 / 4 | -0.000356 to +0.000342 |

MC3D's positive aggregate half-change is not solely caused by one group, but
held-out group g11 has substantial influence. Training block residuals fall
to -0.001536 in block four and rise to +0.001296 in block six, then fall again;
the linear trend does not describe this evolution. SPARTA's held-out temporal
change is small in aggregate but changes sign across groups. The earlier
SPARTA mean-transfer offset is not removed by this time-coordinate audit.

These are descriptive patterns in previously inspected historical engineering
data. They neither establish a deterministic transient nor distinguish it from
correlated fluctuations. They do not establish a solver defect, stationary
noise, confidence coverage or experimental acceptance.

## Next decision

Before selecting a richer mean model, compare training-only temporal structure
with front motion and conserved/flux quantities to identify a plausible common
transient. Freeze any candidate and its evaluation criteria before collecting
fresh runs. The existing held-out bank is now exploratory evidence, not a new
independent qualification set. Do not tune polynomial degree or discard early
windows to make its residuals pass a stationarity gate.

## Reproduction

Run from this worktree with the primary evidence checkout supplied explicitly:

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 tools/diagnose_plateau_time.py \
  --root ../.. --output ../../output/physical-plateau-time-new.json
python3 -m unittest discover -s tests/python
```

The output must not exist. Local report `output/physical-plateau-time-v2.json`
SHA-256: `496428b0af57bcad2599ca6d37dcf97ac6b3be2bb4a220d828118a7361b13678`.
Complete report recomputation matched exactly. Tests cover physical-time labels,
retained offsets, group influence, nonfinite data and invalid window schedules.

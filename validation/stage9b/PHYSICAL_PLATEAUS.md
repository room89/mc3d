# Historical physical plateau diagnostic

Physical acceptance is incomplete. This read-only diagnostic tests whether the
exact constant-plateau assumption of the synthetic covariance prototype can be
transferred to existing MC3D/SPARTA profiles. It does not estimate a spatial
covariance or certify a physical noise model.

## Evidence and fixed protocol

Use the previously verified fixed-grid fluctuation banks. Each solver supplies
36 historical runs of 64 consecutive windows: training groups g00–g07 (24 runs)
and held-out groups g08–g11 (12 runs). Their seed identities are disjoint within
each solver. Restore observations as `fluctuations + mean_profile`; using the
centered arrays alone would hide differences between the split means. These
historical held-out groups have already been inspected in earlier research;
they are not fresh independent validation. The campaign also predates the
matched reference-domain production matrix.

For each solver, locate the 50% crossing of the training ensemble mean. Freeze
physical-grid regions at [-12,-8] and [8,12] collision mean free paths relative
to that center. Require at least one mean free path of clearance from bank edges;
if unavailable, fail rather than move the region. Apply the identical masks to
held-out runs. No per-window alignment, cyclic time wrapping, clipping or
per-run mean subtraction is performed.

Inspect each window's 10–90% transition for overlap with the regions. Compare
first/last 32-window means and first/second spatial-half means, using group means
as replicates for descriptive summaries. Normal intervals use only eight
training/four held-out groups and are not multiplicity-adjusted; they cannot
certify stationarity or formal equivalence.

For covariance transfer, use the regional spatial-average time series and one
training scalar mean. Report its second moment and adjacent-window product
about that fixed mean, separately by time half and split. These quantities retain
mean mismatch, drift and persistent run effects. They are NOT isolated noise
variance or the separable temporal factor expected by the profile fitter.

## Observed limitations

All **4608 windows** have resolved 10–90% transitions and **zero region overlap**.
The transition envelopes stay approximately within [-2.7, 3.1] collision mean
free paths, well inside the selected regions. Geometric clearance alone does
not establish constant plateau means.

| Solver / region | Held-out minus training mean | Held-out/training full second moment |
| --- | ---: | ---: |
| MC3D upstream | -0.00000086 | 1.061 |
| MC3D downstream | -0.00004744 | 1.125 |
| SPARTA upstream | -0.00001195 | 1.007 |
| SPARTA downstream | -0.00066435 | 1.060 |

Upstream spatial-half changes are about +0.00093 to +0.00101 normalized density
in both solvers and both splits, with descriptive intervals excluding zero.
This indicates residual spatial structure rather than an exactly constant mean.
MC3D downstream training first-to-last mean change is +0.00081183; its descriptive
interval is [0.00021390, 0.00140977]. SPARTA upstream changes are -0.00008262
(training) and -0.00013410 (held-out), with descriptive intervals excluding zero.
SPARTA downstream mean transfer has interval [-0.00113315, -0.00019556]. Other
intervals overlap zero; that is not proof of stability.

These small normalized-density effects are **not new physical tolerances or a
verdict against either solver**. They prevent treating these samples as the
exact known-mean plateaus assumed in the synthetic prototype. The second-moment
ratios also include uncertainty from the fitted training mean and cannot alone
prove stationary noise transfer. No candidate region was moved or discarded
based on held-out results.

Next research requirement: separate spatial mean structure and temporal drift
from fluctuations using training runs, then evaluate fixed predictions on
separate runs with uncertainty of the learned mean included. Selecting plateau
edges until held-out diagnostics pass would invalidate that comparison. Existing
historical diagnostics can guide design, but eventual qualification needs fresh
independent evidence and a validated spatial noise model.

## Reproduce

```sh
OPENBLAS_NUM_THREADS=1 python tools/diagnose_physical_plateaus.py \
  --manifest output/shock-physical-fluctuation-banks/manifest.json \
  --output output/physical-plateau-diagnosis.json
```

The tool verifies bank hashes, dimensions and identities through the existing
loader, and rechecks bank hashes after analysis. The report records source hashes,
all split identities, group statistics, physical masks and limitations. It does
not re-extract raw solver trajectories or launch simulations. Detailed output
remains local; original banks are never modified.

# Fixed-grid simultaneous band screen

The [frozen protocol](SIMULTANEOUS_BAND_PROTOCOL.json) extends the
[fixed-profile screen](FIXED_PROFILE_COVERAGE.md) with a conservative band.
It was written before execution. No physical data or solver settings change.

## Method

Keep the same arithmetic group-mean estimator and illustrative noise laws.
At each of 113 nodes use a two-sided Student interval with per-node error
0.05/113 and seven degrees of freedom. The Bonferroni union bound does not
require spatial independence, but it does require valid marginal intervals.
Here the Student approximation for nonlinear group means is checked only
under the prescribed synthetic laws, not assumed validated for physical data.

The critical value is 6.206423865825787. The t7 tail is calculated by the
trigonometric integral after substituting t=sqrt(7)*tan(theta); tests compare
it with independent numerical quadrature and the known pointwise critical
value. Half-widths are 2.62470 times the previous pointwise half-widths.
This band targets the 113 fixed nodes, not every point of a continuous curve.

Use 2,048 independent trials with a new seed. Cases share random draws within
a trial. Whole-profile coverage is one Bernoulli result per independent trial.
The predeclared exploratory criterion is Wilson lower bound >=0.93 and
coverage MCSE <=0.005 in every case. The 0.93 threshold is a screening tolerance,
not a claim that 93% satisfies nominal 95% qualification. No criteria are tuned
following evaluation. No training sample is needed for this fixed estimator.

## Results

| Noise case | Whole-grid coverage | Wilson 95% interval | MCSE | Mean half-width |
| --- | ---: | ---: | ---: | ---: |
| Independent nodes | 99.27% | 98.80–99.56% | 0.188 pp | 0.005439 |
| Same-sign common | 99.32% | 98.86–99.59% | 0.182 pp | 0.005552 |
| Opposite-sign common | 99.02% | 98.50–99.37% | 0.217 pp | 0.005527 |

All frozen screen criteria are met. Coverage is conservative in this model
family, at a substantial width cost. Widths are normalized density units.
Passing this screen does not qualify the noise law for physical outputs,
intrinsic shock shape, gradient thickness or experimental agreement. No
post-result adjustment to the critical value or width was made.

## Execution and reproduction

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 tools/screen_simultaneous_band.py \
  --workers 12 --output output/simultaneous-band-new.json
python3 -m unittest discover -s tests/python
```

The output must not exist. Per-trial streams derive from the frozen seed and
trial index, so scheduling does not change numerical results. The initial run
used eight worker processes and the independent replay used all 12 logical
CPUs; BLAS threads remain one per worker to avoid oversubscription.

Local report: `output/simultaneous-band-v1.json`, SHA-256
`643621699c3fdf7f4c859e9835f796c32b5530e8a414e1ca5b009766a9abd533`.
The replay is `output/simultaneous-band-replay-v1.json`; its full JSON matches
the initial report exactly. All 223 Python tests
passed. Tests cover the Student tail, family-error allocation, Wilson limits
and trial stream reproducibility. Generated trials remain local.

## Next decision before physical computation

Retain this band as a conservative fixed-grid reference; do not spend compute
optimizing its width on the evaluation trials. Specify the next physical
campaign's estimand, independent groups/seeds, sampling schedule and required
precision in advance. If the target is intrinsic shock shape rather than the
laboratory-frame mean, the current band is insufficient and a separate
estimator qualification is required. Historical engineering outputs do not
become a fresh qualification sample because this synthetic screen passed.

No multi-hour physical run was launched for this synthetic check.

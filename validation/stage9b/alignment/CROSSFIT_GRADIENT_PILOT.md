# Independent gradient-selection pilot

Conclusion: **not qualified; do not integrate into physical acceptance**. None of the 12 cases passes both fixed pilot gates. The solver and existing acceptance reports are unchanged.

Each case uses 32 paired repetitions (random seeds 1200000–1200031), three simulation seeds and 64 windows per seed. Noise is 3% of the density jump. Shapes are logistic, asymmetric Gompertz and a broad flat peak. Temporal AR coefficients are 0 and 0.8; spatial correlation scales are 0 and 0.2. The flat-peak CDF is 1/2 + 5/8 (z - z^5/5), z clipped to [-1,1], z=x/2.5; its exact gradient thickness is 4.

For each held-out seed, the other two seeds construct the alignment template and select its derivative maximum. The held-out mean profile is differentiated there without another maximization. The three slopes are averaged before inversion. Alignment still depends on each measured window: this is independent peak selection, not complete independence of all preprocessing.

Gates fixed before this pilot: absolute mean relative bias plus two Monte Carlo standard errors at most 1.5%; standard deviation no more than 1.25 times the template-aligned maximum estimator. These are pilot screening gates, not calibrated physical confidence intervals. No execution errors occurred.

All bias columns below are percentages. Oracle means exact front alignment; known peak additionally uses the analytic maximum location. SD ratio compares independent selection with the template-aligned maximum on the same repetitions.

| Shape | AR | Spatial | Independent bias | Template max bias | Oracle max bias | Oracle known peak bias | Independent bound | SD ratio |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| logistic | 0.0 | 0.0 | +0.959 | -0.836 | -0.809 | +0.879 | 1.539 | 1.339 |
| logistic | 0.0 | 0.2 | +1.328 | -0.094 | +0.229 | +1.590 | 2.221 | 1.495 |
| logistic | 0.8 | 0.0 | +0.835 | -0.773 | -0.782 | +0.638 | 1.488 | 1.389 |
| logistic | 0.8 | 0.2 | +1.432 | -0.069 | +0.187 | +1.489 | 2.225 | 1.116 |
| gompertz | 0.0 | 0.0 | +1.135 | +0.027 | +0.021 | +0.685 | 1.569 | 1.176 |
| gompertz | 0.0 | 0.2 | +2.178 | +0.602 | +0.859 | +1.160 | 2.842 | 1.249 |
| gompertz | 0.8 | 0.0 | +1.347 | -0.024 | -0.003 | +0.536 | 1.800 | 1.252 |
| gompertz | 0.8 | 0.2 | +2.435 | +0.546 | +0.765 | +1.099 | 3.095 | 1.298 |
| flat_peak | 0.0 | 0.0 | +0.933 | -2.653 | -2.722 | +0.813 | 1.579 | 1.920 |
| flat_peak | 0.0 | 0.2 | +1.488 | -3.101 | -3.007 | +1.523 | 2.386 | 1.659 |
| flat_peak | 0.8 | 0.0 | +0.830 | -2.757 | -2.700 | +0.573 | 1.495 | 1.727 |
| flat_peak | 0.8 | 0.2 | +1.471 | -3.106 | -2.973 | +1.422 | 2.251 | 1.111 |

The candidate trades negative maximum-selection bias for positive bias and larger variability. Positive bias also occurs in the exact-front/known-peak control, so template peak error alone cannot explain the result. Inversion of a noisy slope and finite-grid differentiation are possible contributors; this pilot does not separate them. Results must not be compared numerically with the earlier qualification as if paired: its random seeds and repetition counts differ.

Next experiment: separately measure mean slope bias and variance before inversion, first with analytic noiseless profiles, then known front/known peak, then estimated peak. Increase repetitions with a fresh fixed seed range to distinguish Monte Carlo fluctuation from bias. Only after this decomposition should a bias correction or alternative derivative fit be chosen and checked on held-out synthetic cases. Width 10–90 failures remain a separate issue. Any later adoption also requires rebuilding templates and selecting peaks inside resampling and checking interval coverage.

Reproduce from the worktree with a NumPy environment and a new output path:

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python tools/experiment_crossfit_gradient.py \
  --repeats 32 --output output/crossfit-gradient-pilot.json
```

Aggregate results and the full local artifact SHA-256 are in [crossfit-gradient-pilot.json](crossfit-gradient-pilot.json). Full repetitions and fold diagnostics remain in `output/crossfit-gradient-pilot.json` in the primary checkout. No new solver runs were needed.

Follow-up: the [128-repeat decomposition](GRADIENT_BIAS_DECOMPOSITION.md)
does not reproduce the uniformly positive bias of this smaller pilot. Treat
the values above as results of this batch, not calibrated correction constants.

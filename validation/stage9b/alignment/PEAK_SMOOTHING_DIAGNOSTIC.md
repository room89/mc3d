# Peak smoothing and training alignment diagnostics

The controlled Gompertz experiment identifies noisy crossing50 alignment of training windows as a source of peak-location and thickness bias. Smoothing reduces noisy peak-selection error; removing it is not supported. This is diagnostic evidence, not qualification of a replacement estimator or physical acceptance.

## Fixed-front smoothing experiment

Protocol fixed before execution: 12 cases, 256 paired repetitions each, fresh random seeds 3600000–3600255, three seeds × 64 windows, 3% density noise. Four training-profile variants share the same exact front alignment and the same held-out raw derivative (local cubic, full width 1): raw, Gaussian smoothing (sigma 0.25, truncated at 3 sigma), monotone approximation, and smoothing then monotone approximation. Each peak is selected only on the other two seeds. Slopes are averaged before inversion. Eight noiseless controls per shape are also included. All 3072 noisy repetitions completed without errors.

Thickness biases and paired changes below are percentages / percentage points. ± denotes one Monte Carlo standard error of the paired change. Scenarios share random numbers and must not be counted as independent confirmations.

| Shape | AR | Spatial | Raw bias | Smoothed bias | Known-peak bias | Smoothed + monotone minus raw |
|---|---:|---:|---:|---:|---:|---:|
| logistic | 0.0 | 0.0 | +0.988 | +0.326 | +0.081 | -0.662 ± 0.128 |
| logistic | 0.0 | 0.2 | +1.141 | +0.523 | +0.147 | -0.618 ± 0.100 |
| logistic | 0.8 | 0.0 | +1.393 | +0.428 | +0.075 | -0.966 ± 0.113 |
| logistic | 0.8 | 0.2 | +1.678 | +0.865 | +0.206 | -0.813 ± 0.103 |
| gompertz | 0.0 | 0.0 | +0.727 | +0.236 | +0.149 | -0.491 ± 0.071 |
| gompertz | 0.0 | 0.2 | +0.754 | +0.326 | +0.186 | -0.427 ± 0.045 |
| gompertz | 0.8 | 0.0 | +0.871 | +0.227 | +0.141 | -0.644 ± 0.068 |
| gompertz | 0.8 | 0.2 | +0.922 | +0.394 | +0.221 | -0.528 ± 0.055 |
| flat_peak | 0.0 | 0.0 | +0.073 | -0.094 | +0.016 | -0.166 ± 0.133 |
| flat_peak | 0.0 | 0.2 | +0.181 | -0.056 | +0.082 | -0.237 ± 0.116 |
| flat_peak | 0.8 | 0.0 | +0.526 | +0.233 | +0.011 | -0.293 ± 0.139 |
| flat_peak | 0.8 | 0.2 | +0.776 | +0.328 | +0.141 | -0.448 ± 0.118 |

Monotone approximation did not change the selected peaks or measured slopes in this batch, with or without smoothing. This is an observed result for these profiles, not a general guarantee. Smoothing reduces the positive bias by 0.427–0.966 percentage points for logistic/Gompertz. Broad-peak effects are smaller and noisier.

## Training alignment experiment

After the smoothing result, a focused follow-up uses the same 256 paired repetitions in each of the four Gompertz scenarios (1024 repetitions, no errors). Held-out measurements remain aligned at the analytic front. All branches use the same smoothing and monotone operation. The branches differ only in training alignment and final centering:

1. Exact training fronts, no recentering.
2. Exact training fronts, recenter the smoothed template at its 50% crossing.
3. Training windows aligned at their measured 50% crossings, followed by the identical smoothing and recentering, using the existing template builder.

The exact-training branch matches the earlier fixed-front experiment on each repetition. This is a paired controlled comparison, not a comparison between unrelated batches.

| AR | Spatial | Exact + recentered bias | Crossing + recentered bias | Effect of crossing alignment ± SE | Peak RMSE, exact → crossing |
|---:|---:|---:|---:|---:|---:|
| 0.0 | 0.0 | +0.237 | +0.522 | +0.285 ± 0.041 | 0.0412 → 0.0962 |
| 0.0 | 0.2 | +0.320 | +1.197 | +0.878 ± 0.063 | 0.0574 → 0.1561 |
| 0.8 | 0.0 | +0.225 | +0.539 | +0.314 ± 0.040 | 0.0417 → 0.0937 |
| 0.8 | 0.2 | +0.396 | +1.271 | +0.875 ± 0.068 | 0.0598 → 0.1543 |

Crossing-based training alignment adds 0.285–0.878 percentage points to thickness bias. Mean peak error becomes +0.084–0.150 lambda, versus +0.021–0.027 lambda with exact training alignment and the same recentering. Recentring alone changes thickness by −0.007 to +0.002 percentage points in the noisy cases. In the noiseless control, replacing exact training fronts with crossing estimates changes thickness by only −0.000276 percentage points. These contrasts support noisy initial training alignment as a material cause for this asymmetric shape. They do not establish that it is the only error in the full pipeline, which additionally estimates held-out translations.

## Next correction to test

Keep smoothing. Refine training-window translations against a template using only the two training seeds, then rebuild the template; never involve the held-out seed. Fix the number of refinement passes before testing and retain translation-only fits, coverage checks and explicit failures. The initial crossing can remain a starting guess, but should no longer be the final training alignment.

First test the revised template on the same controlled exact-held-out-front comparison to check whether peak error decreases. Then test the complete alignment/measurement pipeline on a fresh validation seed range, all three shapes, and existing bias/variability gates. Rebuilding templates during resampling, uncertainty coverage and width 10–90 qualification remain separate requirements. No iterative correction has been implemented or qualified in this change.

## Reproduction

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python tools/experiment_peak_smoothing.py \
  --repeats 256 --output output/peak-smoothing-experiment.json
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python tools/experiment_template_training.py \
  --repeats 256 --output output/template-training-experiment.json
```

Use new output paths. Aggregate JSON files beside this report preserve the protocol, noiseless results, method summaries, paired contrasts and SHA-256 hashes of full local artifacts. Full artifacts under `output/` retain per-repetition slopes and fold peak errors. No solver simulations or changes to physical-acceptance reports were needed.

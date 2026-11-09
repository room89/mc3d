# Refined training-template estimator

One training-refinement pass plus independent peak selection with a 1.5-lambda derivative fit passes the fresh 3%-noise metric screen (12/12 cases, 256 repetitions per case). This qualifies only the tested metric-bias/variability screen, not confidence-interval coverage or physical acceptance. The additional 1.5%-noise screen also passes 12/12 cases (128 fresh repetitions each, seed start 9600000). Full-procedure coverage and physical precision remain unqualified; see the completed reanalysis below.

## Fixed procedure

For each held-out seed, construct the original smoothed monotone template from the other two seeds. Fit training-window translations against that template once; reconstruct the average directly from original measured samples and rebuild/recenter the smoothed template. Align the held-out windows against this refined template. Choose the derivative peak from the independent template and measure the held-out mean derivative there. Average three slopes before inversion. Width 10–90 is measured on the mean aligned profile. Only translations are fitted.

Exactly one refinement pass is supported; no convergence tuning or shape/width scaling is applied. The candidate remains opt-in (`refine_passes=1`, or reanalysis `--alignment template_refined`). The original template mode and crossing50 defaults retain their definitions.

## Selection and fresh validation

1. One pass, derivative width 1: 64 repetitions, seeds 4800000–4800063. Gradient screen passes 10/12 cases; logistic variability/precision prevents qualification.
2. Width increased to 1.5 before a fresh 64-repeat screen, seeds 6000000–6000063. Gradient passes 12/12; an additional width check reaches 0.526% in one case, above the existing 0.5% limit. This batch is not a joint success.
3. Parameters frozen; 256 fresh repetitions, seeds 7200000–7200255. Joint width/gradient screen passes 12/12. All cases have zero errors.

Each case has three seeds × 64 windows. Noise is 3% of the density jump; temporal AR coefficients are 0/0.8 and spatial scales 0/0.2. Shapes are logistic, Gompertz and a flat peak. Scenarios share random numbers; do not pool them as independent confirmations. Gates are unchanged: absolute mean bias plus two Monte Carlo SE at most 0.5% for width, 1.5% for gradient; gradient SD at most 1.25 times the original independent-peak estimator (unrefined template, width 1) on the same repetitions. These Monte Carlo SE are not physical uncertainty intervals.

| Shape | AR | Spatial | Width bound (%) | Gradient bound (%) | Gradient SD ratio |
|---|---:|---:|---:|---:|---:|
| logistic | 0.0 | 0.0 | 0.180 | 0.309 | 0.574 |
| logistic | 0.0 | 0.2 | 0.277 | 0.388 | 0.847 |
| logistic | 0.8 | 0.0 | 0.116 | 0.402 | 0.616 |
| logistic | 0.8 | 0.2 | 0.291 | 0.683 | 0.941 |
| gompertz | 0.0 | 0.0 | 0.110 | 0.377 | 0.560 |
| gompertz | 0.0 | 0.2 | 0.206 | 0.195 | 0.748 |
| gompertz | 0.8 | 0.0 | 0.096 | 0.420 | 0.591 |
| gompertz | 0.8 | 0.2 | 0.208 | 0.350 | 0.772 |
| flat_peak | 0.0 | 0.0 | 0.170 | 0.204 | 0.703 |
| flat_peak | 0.0 | 0.2 | 0.263 | 0.494 | 0.990 |
| flat_peak | 0.8 | 0.0 | 0.158 | 0.225 | 0.711 |
| flat_peak | 0.8 | 0.2 | 0.265 | 0.482 | 1.013 |

## Full-procedure uncertainty

`shock_refined_estimator.py` resamples circular blocks of 8/16 raw windows independently within each original seed, then rebuilds all templates, training shifts, held-out shifts and peaks. Original seed identities are retained to avoid a resampled copy of a held-out seed entering its training fold. Between-seed variation is retained through a conservative seed-t envelope. Its cross-fold dependence and finite-sample coverage remain under calibration. Processing sensitivity compares derivative widths 1 and 2 against the fixed 1.5 estimate.

Initial diagnostic C result: width 4.652915 with statistical radius 0.049844; gradient thickness 4.096279 with statistical radius 0.063290 and processing sensitivity 0.053865. Density simultaneous radius is 0.011542. These are unqualified uncertainty estimates; they do not pass the physical precision requirement. Full-matrix reanalysis and coverage calibration are separate jobs.

The implementation caches only immutable source coordinates/density and each window’s initial 50% crossing. Resampled templates and fitted shifts are never cached. Batched translation fitting matches the scalar implementation in regression checks, and cached processing reproduces the physical C point metrics and density exactly.

## Reproduction

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python tools/qualify_refined_template.py \
  --fit-width 1.5 --scope metrics --seed-start 7200000 --repeats 256 --jobs 3 \
  --output output/refined-template-validation.json
python tools/reanalyze_shock_mean.py output/shock-validation-restart/series \
  --alignment template_refined --bootstrap-repeats 128 --jobs 3 \
  --output output/refined-full-reanalysis.json
python tools/calibrate_refined_intervals.py --replicates 32 --bootstrap-repeats 128 \
  --jobs 3 --output output/refined-interval-coverage-pilot.json
```

Reanalysis saves a `.partial.json` checkpoint after each complete variant; calibration saves per-repetition `.parts` checkpoints. Add `--resume` only for an unfinished output. Input and algorithm hashes prevent resuming a different procedure. Existing final outputs are never overwritten.

Aggregate screen JSON files beside this report contain full-artifact hashes. Raw repetitions remain in the permanent project `output/` directory. The original physical report is unchanged.

## Completed physical reanalysis and coverage pilot

All ten saved variants (MC3D A–G and SPARTA C/D/E), including first/second halves,
have been reprocessed. All nine convergence/cross-code comparisons fail the
unchanged precision requirements. All ten half comparisons remain unconfirmed.
For C versus SPARTA C:

| Metric | Difference | Combined uncertainty including processing | Precision limit |
|---|---:|---:|---:|
| Width 10–90 | +0.009745 | 0.122078 | 0.069648 |
| Gradient thickness | −0.002737 | 0.274108 | 0.061485 |
| Density maximum difference | 0.00351 | 0.023434 | 0.005000 |

The gradient means differ by about 0.07%, but this does not establish agreement
at the required confidence and precision. The separate plateau diagnostics find
acceptable mean plateaus and boundary distances in all variants. Most plateau
half comparisons still lack precision. These diagnostics do not silently replace
the original every-window stationarity rule.

The full-procedure coverage pilot on correlated Gompertz has 29/32 width hits,
32/32 gradient hits and 32/32 simultaneous density hits. Width's 95% Wilson
interval is approximately 76–97%; this small pilot cannot certify 95% coverage.
Block-only metric intervals give 26/32 and 30/32 hits. Dropping the between-seed
envelope to obtain smaller physical intervals is therefore not justified.

A separate exploratory independent-group method uses 12 independent groups of
three seeds. It combines group profiles and slopes and resamples whole groups;
no group participates in another group's templates. A 32-repeat correlated-noise
pilot gives width 30/32 and gradient 31/32 for both logistic and Gompertz; density
is 32/32 and 31/32 respectively. This does not yet qualify group inference or
justify a large physical campaign. The group t multiplier uses df=11 (2.201),
from the [NIST Student-t table](https://www.itl.nist.gov/div898/handbook/eda/section3/eda3672.htm).

The experimental Alsmeyer source remains unavailable for verified digitization.
The user confirmed that they do not have an accessible copy or dataset.
The overall physical acceptance status is still `incomplete`.

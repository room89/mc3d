# Cross-fitted template alignment

Implementation uses independent seeds to construct a monotone smoothed template,
then fits only translation with Huber loss over its 10–90% transition.
The default remains crossing50; enable the candidate with `--alignment template`.
No solver changes or new physical runs are involved.

The noisy alignment qualification finished: 11 of 16 cases pass; the candidate
is not qualified. The separate noiseless/edge-case checks and final plots remain
outstanding. A follow-up independent gradient-selection pilot passes 0 of 12
cases; see [the pilot report](CROSSFIT_GRADIENT_PILOT.md). Neither candidate is
ready for physical acceptance.

Commands (NumPy environment):

```sh
python tools/qualify_shock_alignment.py --output output/alignment-qualification.json
python tools/reanalyze_shock_mean.py output/shock-validation-restart/series \
  --alignment template --output output/template-profile-reanalysis.json
```

The reanalysis keeps the original validation report unchanged and remains
experimental. Its intervals do not fully include template-estimation uncertainty.
Input hashes, per-seed templates and per-window shifts are retained in the new
report. Local artifacts use permanent project output directories.

The follow-up [slope/reciprocal decomposition](GRADIENT_BIAS_DECOMPOSITION.md)
completed 128 fresh repetitions per case. Reciprocal nonlinearity contributes
only 0.022–0.109 percentage points for the template candidate; no bias correction
or physical-acceptance change was applied.

The [controlled smoothing/training study](PEAK_SMOOTHING_DIAGNOSTIC.md)
finds that noisy crossing50 alignment of training windows adds 0.285–0.878
percentage points of thickness bias on Gompertz. Smoothing helps; the next
candidate should refine training translations without using the held-out seed.

The [one-pass refined estimator](REFINED_TEMPLATE_ESTIMATOR.md) passes the
fresh 256-repeat, 3%-noise metric screen in all 12 cases. This is distinct from
the unrefined candidate above. Coverage and physical precision are still open.

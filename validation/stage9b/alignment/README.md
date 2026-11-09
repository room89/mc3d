# Cross-fitted template alignment

Implementation uses independent seeds to construct a monotone smoothed template,
then fits only translation with Huber loss over its 10–90% transition.
The default remains crossing50; enable the candidate with `--alignment template`.
No solver changes or new physical runs are involved.

Synthetic qualification is running. Do not regard this candidate as qualified
until all fixed cases finish and the separate noiseless/edge-case checks are
reviewed. The complete planned qualification and final plots are not yet delivered.

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

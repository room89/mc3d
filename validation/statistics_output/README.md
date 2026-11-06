# Statistics output performance

The stage-9b run sampled before this change spent approximately 97% of main-thread
samples in `StatisticalDiagnostics::Write`. Legacy zero snapshot interval triggered
serialization of all 64 windows on every step, including unchanged windows.

This change caches the elapsed exposure and completion state of each written
window. It retains per-step accumulation, window-boundary flushes, interrupted-run
flushes and final output. Dynamic surface temperatures are still refreshed at
output time. Generated normal-shock setups use ten snapshot intervals per setup
duration; the stage-9b extension of the duration inherits that explicit interval.
Legacy configs with zero snapshot interval keep their snapshot behavior.

Reproduce the benchmark with Release binaries:

```sh
python3 tools/benchmark_statistics_output.py \
  --baseline /path/to/old/MC3dSolver \
  --candidate build/MC3dSolver \
  --output-dir /tmp/statistics-output-benchmark
```

The benchmark uses a 400-cell shock, seed 123, 64 windows, normalized end time
0.02 and three repeats per variant. Resolution diagnostics are disabled to isolate
statistics/snapshot cost. It compares every final statistics table, conservation
file and final snapshot byte for byte. It is a short performance/regression check,
not physical acceptance or a full-series runtime prediction.

One worker is used for this exact comparison: two preliminary baseline-only runs
with four workers had identical initial snapshots but different final snapshots
at seed 123. Both manifests reported full mode and four workers throughout.
The cause of that existing parallel variability has not been diagnosed here.

## Interrupted series audit

The original `/tmp/shock-validation-series` is retained without edits. Its report
records successful execution (not physical acceptance) for A-123, A-456, A-789 and
B-123. B-456's manifest still says running, but no solver or series runner remains.
Its log ends during step 9299, around normalized time 0.8518, with no error or
shutdown record. No MC3dSolver crash report was found in the user's DiagnosticReports.
The available evidence does not establish the termination cause.

The completed-run ledger accounts for 16144.92 seconds; the incomplete run's last
log records another 1414.77 seconds. These are a lower bound of 17559.69 seconds
(4.88 hours), not a completed-run budget reconciliation. The original series occupies
13 GiB; the filesystem had about 20 GiB available during inspection. This is not
proof that disk space caused termination.

A new binary must use a fresh series directory. The existing runner deliberately
rejects resume with a changed binary hash. Preserve original results and account
for the incomplete run before allocating the remaining 24-hour total budget.

## Results

Release build, local Apple Silicon machine; medians of three serial repeats:

| Variant | Seconds | Speedup |
| --- | ---: | ---: |
| Original binary, every-step snapshots | 37.384 | 1.00× |
| Cached statistics, every-step snapshots | 1.617 | 23.11× |
| Cached statistics, snapshot interval 0.002 | 1.275 | 29.33× |

All 131 compared files matched byte for byte across all nine runs. Exact input,
binary hashes and individual timings are in `benchmark.json`. Validation also
passed all 180 CTests and 31 Python unit tests. The integration regression covers
unchanged completed-window files, overlapping windows, dynamic surface output,
dense/sparse cadence equivalence, interruption and zero-duration output.

All output hashes listed in the four completed manifests were rechecked against
the original binary hash and matched; see `original_series_hash_check.json`.

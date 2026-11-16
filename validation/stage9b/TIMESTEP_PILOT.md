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

Evidence directory: `output/shock-timestep-pilot-v1`, with a copied protocol,
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

# Population sensitivity control

The completed domain diagnostic did not establish stationary fixed-coordinate
profiles. Test particle-number sensitivity independently using variant E,
longitudinal domain multiplier 2, warmup 20 original-domain transit times,
64 equal sampling windows over [20,28], and population multiplier 4.
This changes 160 to 640 model particles per upstream cell and divides model
particle weight by four. Grid, reservoirs, physical domain, gas model and
physical time schedule remain unchanged. Use the same pinned MC3D and SPARTA
binaries as DOMAIN_DIAGNOSTIC.md, seeds 210001, 210002 and 210003, with a new
output directory. Same seed labels do not imply paired random trajectories.

Compare with the completed multiplier-1 doubled-domain campaign. Retain all
runs and windows, with no automatic retry or favorable-window selection.
Report original-coordinate mean density differences, run-level front and
width descriptions, half-run changes and full-domain density inventories.
Do not fit alignment as an acceptance correction. A fourfold particle count
suggests a factor-two noise reduction only under an independent-particle
approximation; it is not a promised variance reduction or a pass criterion.
Three seeds provide exploratory evidence, not calibrated uncertainty.

The runner option is `--population-multiplier 4`. Its default is 1, and changing
it on resume is rejected before existing evidence is rewritten. This control
does not itself close numerical convergence, experimental comparison or
stages 9–11.

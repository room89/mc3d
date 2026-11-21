# PR 59 archive and split

Original snapshot: `2c58c95701e45c7e3b196e7e86c441fc8496b594`. Base: `b5725aeec7ecf23e0a764c2133b9b3381aed51ff`.
The original branch is retained; running computations continue from its unchanged
checkout. This split does not change numerical implementations or acceptance gates.

[PR59_FILE_DISPOSITION.tsv](PR59_FILE_DISPOSITION.tsv) accounts for all 208 changed
files. `retain` preserves original bytes; `rewrite` replaces cumulative prose;
`archive` excludes historical output or superseded diagnostics from the new PRs.
Stage numbers identify the seven sequential review branches. Tests for archived
standalone diagnostics are archived with their implementations.

Recreate the local full tracked-source archive from the pinned commit:

```sh
git archive --format=tar.gz --prefix=mc3d-pr59-snapshot/ \
  --output=pr59-original-snapshot.tar.gz 2c58c95701e45c7e3b196e7e86c441fc8496b594
shasum -a 256 pr59-original-snapshot.tar.gz
```

Local archive: 2503292 bytes; SHA-256
`3be7fbaa45837cd09690672e2efb101a4ddbf3801b938411abe72cf96d79b6f4`.
All 208 changed files were read back from this archive and checked against the
individual SHA-256 values in the inventory. Archive compression may differ across
Git/zlib versions; individual file hashes identify the source independently.
Raw simulation output is not tracked and is not included in this source archive.

Primary digitizations and reference inputs are retained in stage 4. Shared
`experiment_crossfit_gradient.py` is retained because calibrations import it.
The last two stages isolate global-search and likelihood research; neither is a
qualified production estimator. Detailed archived evidence is available from the
[pinned source tree](https://github.com/room89/mc3d/tree/2c58c95701e45c7e3b196e7e86c441fc8496b594).

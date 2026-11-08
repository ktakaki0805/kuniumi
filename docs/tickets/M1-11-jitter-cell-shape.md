# M1-11 Default jitter and the spread of cell areas

- Status: draft
- Milestone: M1
- Depends on: M1-06
- Crates: kuniumi-sphere (only if the default changes)
- Design: 球面グリッドとデータ構造
- Decisions: D-14

## Goal

The owner decides whether the default jitter of the grid's points stays at
0.35 or is lowered, knowing what it does to the spread of cell areas and to
the share of hexagonal cells; the code and the documents agree with the
decision.

## Scope

In:
- With the owner, decide `points::DEFAULT_JITTER` (now 0.35, chosen in M1-04
  for the closest pair and to break up the spiral's rows).
- If it changes: the new value with its reasons in the documentation of
  `DEFAULT_JITTER` (closest pair measured again with
  `examples/point_spacing.rs`), and the fingerprints that depend on it
  recorded again as a world-breaking change (`points`, `delaunay` and `grid`
  tests), listed in `CHANGELOG.md`.
- If it stays: the area ratio of M1-06 (about 3.4 for 650,000 cells) recorded
  as accepted.

Out:
- Other point sets (icosahedral, random), and any change to the algorithms.

## Implementation notes

Found in M1-06. Its acceptance criterion recorded the ratio of the largest to
the smallest cell area and expected it below 3. With the default jitter the
ratio, an extreme-value statistic, grows with the number of cells (seed 1;
values in units of the mean area, measured with a throwaway program):

| Cells | Jitter 0.25 | Jitter 0.30 | Jitter 0.35 |
| --- | --- | --- | --- |
| 1,000 | 1.93 | 2.22 | 2.56 |
| 10,000 | 2.14 | 2.52 | 2.98 |
| 100,000 | 2.20 | 2.64 | 3.25 |
| 650,000 | 2.34 | 2.82 | 3.43 (seeds 2, 3: 3.39, 3.45) |

At 650,000 cells and jitter 0.35 the smallest cell is 0.48 of the mean area
and the largest 1.63; the 1st and 99th percentiles are 0.67 and 1.35, the
0.1th and 99.9th 0.59 and 1.45. The plain spiral (jitter 0) gives 0.94 to
1.05.

The jitter also changes the topology unevenly with latitude: the plain
spiral is a nearly hexagonal lattice whose shape varies with latitude, and
the same jitter breaks it more where the lattice is less regular. Share of
cells with six neighbours by 10° band from 90° N to 90° S (650,000 cells,
seed 1):

- jitter 0: 0.91 0.95 0.99 0.98 1.00 0.97 1.00 1.00 1.00 1.00 1.00 1.00 0.97
  1.00 0.98 0.99 0.95 0.91
- jitter 0.25: 0.68 0.60 0.68 0.70 0.70 0.43 0.60 0.83 0.92 0.92 0.84 0.60
  0.43 0.70 0.70 0.69 0.59 0.69
- jitter 0.35: 0.54 0.50 0.55 0.56 0.55 0.42 0.50 0.62 0.68 0.68 0.62 0.50
  0.43 0.55 0.56 0.55 0.51 0.56

The ring-size map of M1-06 (`out/M1-06/ring_size_equirect_650000.png`)
shows these bands; the area map shows none. Whether the bands matter for
drainage networks (the reason for the jitter) is part of the decision.

## Acceptance criteria

- [ ] The owner's decision and its reasons are recorded in this ticket's
      Result and in the documentation of `DEFAULT_JITTER`.
- [ ] If the default changes: the area ratio and the six-neighbour shares by
      latitude band for 650,000 cells are measured again and recorded, every
      changed fingerprint is listed under *World-breaking* in `CHANGELOG.md`.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-sphere --release grid -- --nocapture
cargo run --release -p kuniumi-sphere --example point_spacing -- 1000
```

## Result

_(filled in when the ticket is implemented)_

# M2-13 Age steps around very slow ridge cells

- Status: draft
- Milestone: M2
- Depends on: M2-11
- Crates: kuniumi-sim (stages/tectonics)
- Design: 地殻から標高へ / 海洋底：年代から水深; プレートテクトニクス / 境界ごとの処理
- Decisions: none

## Goal

Decide whether the remaining steps of ocean-floor age next to very slow
ridge cells are smoothed, so that the sea floor of M2-05 has no step of
kilometres beside a ridge.

## Scope

In:
- The finding of the review of M2-11: the time-nearest growth of M2-11
  removes most age seams (pairs of neighbouring oceanic cells of one plate
  more than 20 Myr apart: 82 at 100,000 cells, seed 1, against 1,955 with
  M2-04), but a ridge cell that spreads very slowly is still the first to
  reach its own neighbours when no faster front is near, as in clusters
  of slow cells near a relative pole. At 100,000 cells, seed 1, ridge cell
  39194 (0.096 cm/yr, 13° N 80° W) has a neighbour of ranked age 105.9 Myr,
  a step of about 3.2 km in M2-05's sea floor; 20 of the 82 seams touch a
  ridge cell, and the largest step is 106 Myr. At 650,000 cells, seeds
  1–10: no seam in six worlds; 3–110 seams in the others (seed 5: 110, 15
  next to a ridge cell), with largest steps of 25–75 Myr.
- Options to decide with the owner: a floor on the half spreading rate of
  a ridge cell used for the times (for example the plate's mean rate, or a
  fixed minimum); rates smoothed along each ridge; or keep the steps
  until M6.

Out:
- The time-evolving tectonics of M6.

## Implementation notes

To be refined: the chosen option and its effect on the seams, the mean
age and the land fraction of M2-11 (which stays within 20–40 %).

## Acceptance criteria

- [ ] Seams and the largest step next to a ridge cell reported at
      100,000 cells (seed 1) and 650,000 cells (seeds 1–10), within the
      bound the owner sets.
- [ ] Maps of `crust_age` of the affected seeds viewed and described.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo run --release -p kuniumi-cli -- generate --seed 5 --out out/M2-13
```

## Result

_(filled in when the ticket is implemented)_

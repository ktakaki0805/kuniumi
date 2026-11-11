# M2-12 Ocean floor of plates without a ridge

- Status: draft
- Milestone: M2
- Depends on: M2-11
- Crates: kuniumi-sim (stages/tectonics)
- Design: 地殻から標高へ / 海洋底：年代から水深; プレートテクトニクス / プレート速度と分裂
- Decisions: none

## Goal

Decide what age the ocean floor of a plate without a divergent boundary
gets, so that no large uniform patch of the oldest floor remains in the
one-shot planet.

## Scope

In:
- The finding of M2-11: the ages come from fronts that grow from the
  divergent cells of each plate, ranked onto Earth's area–age
  distribution; a plate whose boundaries are all convergent or transform
  has no front, and its ocean floor is 180 Myr old everywhere, a flat
  abyss about 6 km deep with a step at the plate boundary. Measured with
  Earth's parameters (12 plates): at 650,000 cells, seeds 1–10, such plates
  hold 10.3–14.7 % of the ocean floor in four worlds (seeds 3, 6, 9, 10),
  whose mean ocean-floor age becomes 72–78 Myr instead of 60 Myr; at
  10,000 cells up to 20.5 % (seed 4).
- Options to decide with the owner: rank the floor of such plates too, by
  the time from another source (for example their transform boundaries,
  or the boundaries of any type at the plate's mean rate); give them a
  ridge (a divergent boundary where the relative motion is slowest);
  or keep them at 180 Myr.

Out:
- The time-evolving tectonics of M6, in which every ocean floor has a
  ridge of origin.

## Implementation notes

To be refined: the chosen option, how it enters `ridge_times` or
`earth_area_age`, and its effect on the land fraction of M2-11 (which
stays within 20–40 %).

## Acceptance criteria

- [ ] With Earth's parameters, seeds 1–10 at 650,000 cells: the share of
      the ocean floor at exactly 180 Myr and the mean age reported, as the
      owner decides.
- [ ] Maps of `crust_age` of the affected seeds viewed and described.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo run --release -p kuniumi-cli -- generate --seed 3 --out out/M2-12
```

## Result

_(filled in when the ticket is implemented)_

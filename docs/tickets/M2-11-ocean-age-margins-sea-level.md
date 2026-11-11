# M2-11 Ocean-floor age, passive margins and sea level

- Status: draft
- Milestone: M2
- Depends on: M2-04
- Crates: kuniumi-sim (stages/tectonics)
- Design: 地殻から標高へ; プレートテクトニクス / プレート速度と分裂; 検証方法
- Decisions: none

## Goal

Decide, before M2-05, how the one-shot planet reaches Earth's land fraction
of about 29 % from 40 % continental crust: the ocean floor of M2-04 is much
older than Earth's, and no continental margin is submerged.

## Scope

In:
- The finding of M2-04 (its Implementation notes, "Prototype
  measurements"): with the age of M2-04 (distance to the nearest ridge of
  the plate ÷ its half spreading rate, capped at 180 Myr) the area-weighted
  mean age of the ocean floor is 79–144 Myr (Earth: about 60 Myr) and
  23–68 % of it sits at the cap; with M2-05's formulas (Parsons & Sclater
  depth, Airy isostasy, 1.335e18 m³ of water) the sea level falls 0.7–1.5 km
  below the datum and the land fraction is 40.0 % for every seed (seeds 1–5
  at 100,000 cells, 1–10 at 10,000).
- The age is discontinuous inside plates: cells take their age from the
  nearest ridge cell, whose rate may differ greatly from its neighbours'
  (a ridge's rate falls to zero at its relative pole). Measured in the
  review of M2-04 at 100,000 cells, seed 1: 1,955 pairs of neighbouring
  oceanic cells of the same plate differ by more than 20 Myr, seams
  thousands of kilometres long reach 127–159 Myr, and next to a ridge cell
  spreading at 0.10 cm/yr the age jumps from 0 to 180 Myr in one step, a
  step of about 3.7 km in M2-05's sea floor. Taking for each cell the ridge
  it is reached from first by growth at each ridge cell's own half rate
  (the "time-nearest" variant of M2-04's prototype) would largely remove these
  seams; smoothing the rates along the ridges is another option.
- Even Earth's triangular area–age distribution would give a mean depth of
  4.9 km against the 4.4 km that the water fills on 60 % of the surface, so
  the design's thinned continental margins (地殻から標高へ / 大陸: 25–30 km,
  submerged as continental shelves) are needed as well.
- Options to decide with the owner: thin passive margins (continental cells
  near oceanic crust of the same plate) and continental rifts; faster
  ridges, for example plate speeds that grow with the length of their
  subduction boundaries (design, プレート速度と分裂), which changes M2-02
  and is world-breaking; another age model; or a combination.

Out:
- The time-evolving tectonics of M6, which produces margins and ages by
  itself.

## Implementation notes

To be refined: the chosen option, its parameters, and how the land
fraction and the mean ocean-floor age are measured (with M2-05's formulas,
or after M2-05 if the owner reorders the tickets).

## Acceptance criteria

- [ ] With Earth's parameters, seeds 1–10: land fraction within the range
      that M2-05's criterion asks (20–40 %), measured as the owner decides.
- [ ] Maps of the changed fields viewed and described.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo run --release -p kuniumi-cli -- generate --seed 1 --cells 100000 --out out/M2-11
```

## Result

_(filled in when the ticket is implemented)_

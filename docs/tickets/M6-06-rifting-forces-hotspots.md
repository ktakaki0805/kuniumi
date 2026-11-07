# M6-06 Rifting, plate speeds and hotspots

- Status: draft
- Milestone: M6
- Depends on: M6-04, M6-05
- Crates: kuniumi-sim (stages/tectonics)
- Design: プレートテクトニクス / プレート速度と分裂
- Decisions: none

## Goal

Plates evolve on their own: speeds follow slab pull, large old continental
plates split (preferring old sutures, giving Wilson cycles), and fixed
hotspots leave age-ordered island chains.

## Scope

In:
- Plate speed proportional to the length of its subducting boundary, reduced
  for plates with much continental crust; capped to 1–10 cm/yr.
- Rifting: probability grows with plate area and time since the last split; the
  rift line follows old sutures where they exist, otherwise a seeded great
  circle segment; the two halves get diverging rotations; rifting thins the
  crust (25–30 km margins become shelves).
- Hotspots: a few points fixed in the mantle frame (stream
  `"tectonics/hotspots"`), uplifting and thickening the crust passing over them.

## Implementation notes

To be refined: the speed law's constants, the rift probability model, hotspot
count and strength.

## Acceptance criteria

- [ ] Plate speeds stay within 1–10 cm/yr over a full run.
- [ ] A long run splits the supercontinent and opens an ocean with passive
      margins; at least one rift reuses an old suture in a seeded test.
- [ ] A hotspot under a moving oceanic plate produces an age-ordered chain.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-sim rifting
```

## Result

_(filled in when the ticket is implemented)_

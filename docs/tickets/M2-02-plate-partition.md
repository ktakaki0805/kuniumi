# M2-02 Plate partition and Euler poles

- Status: ready
- Milestone: M2
- Depends on: M2-01
- Crates: kuniumi-sim (stages/tectonics), kuniumi-world
- Design: プレートテクトニクス / 表現; プレートテクトニクス / プレート速度と分裂
- Decisions: none

## Goal

The surface is divided into irregular, contiguous plates, each rotating about
its own Euler pole at Earth-like speeds. This is the one-shot tectonics that
M6 later replaces with time evolution.

## Scope

In:
- `Settings::plate_count` (default 12; valid 2–64).
- Partition by seeded multi-source growth: choose seed cells, then grow plates
  from a priority queue whose priorities carry random noise per cell (stream
  `"tectonics/partition"`), which gives irregular boundaries. The heap is keyed
  by `(priority via total_cmp, cell index)`, so ties are deterministic.
- Euler pole: a random unit vector; angular speed chosen so that the plate's
  maximum surface speed lies between 1 and 10 cm/yr (3.2 × 10⁻¹⁰ to
  3.2 × 10⁻⁹ m/s).
- Outputs: `plate_id` (class), `plate_velocity_east` and
  `plate_velocity_north` (m/s), from `ω × r` scaled by the radius.

Out:
- Continents and boundaries (M2-03).

## Acceptance criteria

- [ ] Every cell belongs to exactly one plate and every plate is connected
      (graph search test).
- [ ] Plate sizes are irregular: the ratio of the largest to the smallest plate
      area is between 3 and 50 for seeds 1–20 (distribution reported).
- [ ] All surface speeds are within the range above, for planets of Earth's
      radius and of half Earth's radius.
- [ ] Golden fingerprint of `plate_id` for 10,000 cells, seed 1.
- [ ] `out/M2-02/plate_id.png` viewed and described.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo run -p kuniumi-cli -- generate --seed 1 --cells 100000 --out out/M2-02
```

## Result

_(filled in when the ticket is implemented)_

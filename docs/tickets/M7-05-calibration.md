# M7-05 Calibration against the validation metrics

- Status: draft
- Milestone: M7
- Depends on: M7-04
- Crates: kuniumi-sim, kuniumi-cli, `tools/` (Python, optional)
- Design: 検証方法; 侵食と地形発達 / 調整の目安
- Decisions: D-11, D-22

## Goal

Tune the default parameters so that, over many seeds, Earth's values fall inside
the distribution of every validation metric.

## Scope

In:
- Use `kuniumi validate` (M5-06) on 20–50 seeds with Earth-like parameters.
- Targets from the design: land fraction about 29%, land mean about 0.8 km and
  ocean mean about 3.7 km, highest peak 8–9 km, main orogens 2–4 km on average,
  deepest point about 11 km, Hack exponent about 0.6, dry land about 30%
  concentrated at 15–35° and on cold-current west coasts, clear rain shadows.
- Parameter changes recorded with before and after metrics; defaults updated
  (world-breaking change recorded).
- `docs/roadmap.md`: M7 marked as met.

## Implementation notes

To be refined: which parameters to tune and in what order; whether a Python
sweep script under `tools/` helps (D-11).

From M3-05 (the owner kept the default erosion settings and asked that
the issues be recorded here): with the defaults (1,000 steps of 10⁴ yr,
`K` = 7.5 × 10⁻¹⁰, `m` = 0.45, `K_d` = 5 × 10⁻³ m²/yr) the highest peaks
of the Earth worlds of seeds 1–5 at 650,000 cells reach 14,955, 7,949,
12,521, 5,185 and 8,606 m (design: 8–9 km); the steady-state slope
`U / (K Qᵐ)` steepens as the cells shrink, so the peaks depend on the
resolution (5,960 m for seed 1 at 100,000 cells). The interiors are planed
to the sea level (M7-03). The default `K_d` is sub-grid at 28 km between
cell centres (M3-04: an effective coarse-grid value of about 5,300 m²/yr
at 100,000 cells visibly widens the valleys); the owner chose no separate
calibration ticket, so it belongs here.

## Acceptance criteria

- [ ] Every metric's distribution contains Earth's value (table in the Result).
- [ ] The terrain-only test of the design (Earth's major orogens' uplift with
      erosion only) compared with real topography and described.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo run --release -p kuniumi-cli -- validate --seeds 20 --out out/M7-05
```

## Result

_(filled in when the ticket is implemented)_

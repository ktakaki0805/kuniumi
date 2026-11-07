# M4-01 Insolation by latitude and month

- Status: draft
- Milestone: M4
- Depends on: M2-01
- Crates: kuniumi-sim (stages/climate), kuniumi-world
- Design: 気候モデル / 構成要素
- Decisions: D-12

## Goal

Daily-mean top-of-atmosphere insolation (W/m²) for every cell and each of the
12 months, from stellar flux, obliquity and the orbit, for any planet in the v1
range.

## Scope

In:
- Solar declination and planet–star distance through the orbit; add
  `perihelion_longitude_rad` to `PlanetParams` with an Earth default.
- Daily-mean insolation
  `Q = (S/π)(ā/r)² (h₀ sinφ sinδ + cosφ cosδ sin h₀)` with
  `cos h₀ = −tanφ tanδ`, clamped for polar day and night.
- Monthly means from several sub-steps per month.

Out:
- Tidally locked insolation (v2; that regime stays rejected).

## Implementation notes

To be refined: month boundaries in orbital longitude versus time (Kepler's
equation), number of sub-steps; reference Berger (1978).

## Acceptance criteria

- [ ] Earth (S = 1361 W/m², obliquity 23.44°): annual-mean insolation ≈
      415.5 W/m² at the equator and ≈ 172.3 W/m² at the poles (within 0.5%;
      values from an independent numerical integration); the global annual mean
      equals S/4 × (1 − e²)^−½ within 1e-6 relative.
- [ ] Polar night gives exactly 0; with obliquity 60° or 90° the poles receive
      more annual-mean insolation than the equator.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-sim insolation
```

## Result

_(filled in when the ticket is implemented)_

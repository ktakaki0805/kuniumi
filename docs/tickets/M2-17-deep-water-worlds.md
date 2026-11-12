# M2-17 Sea level of deep water worlds

- Status: draft
- Milestone: M2
- Depends on: M2-08
- Crates: kuniumi-world or kuniumi-sim
- Design: 地殻から標高へ / 海面：海水量から決める
- Decisions: D-12

## Goal

No valid planet gets a sea level far from physical because its ocean is
deep compared with its radius.

## Scope

In:
- One of two ways, chosen with the owner when this ticket is refined:
  - a check across parameters in `PlanetParams::validate`: the ocean's
    mean depth `V / (4πR²)`, or its volume against the planet's
    `(4/3)πR³`, within a bound (as `ParamError::Inconsistent`, M2-08); or
  - a sea level from a spherical shell for the water above the highest
    land in `stages::crust::sea_level_m`: depth
    `(R_top³ + 3(V − C)/(4π))^(1/3) − R_top` above the brim, instead of
    the thin layer `(V − C) / (Σ a_i R²)` (world-breaking only for water
    worlds whose layer is thick enough to change the bits).

Out:
- Ocean volumes and radii themselves (bounded by M2-08).

## Implementation notes

Found in the review of M2-08 (R1-2): the bounds allow 1e20 m³ of water
(`MAX_OCEAN_VOLUME_M3`) on a 1,500 km planet (`MIN_RADIUS_M`), seven times
the planet's volume of 1.41e19 m³. The Crust stage's water-world formula
then gives a sea level of +3,533,314 m (CLI, 20,000 cells, seed 1), while a
spherical shell holding the same volume is 1,509 km deep (Python,
`(R³ + 3V/4π)^(1/3) − R`). At 13,000 km the two agree within 0.4 %
(47,087 m against 46,918 m). Later stages (climate) read the sea level.

## Acceptance criteria

- [ ] To be made measurable when the ticket is refined.
- [ ] `mise run ci` passes.

## Verification

```sh
mise exec -- cargo test -p kuniumi-world
mise exec -- cargo test -p kuniumi-sim crust
```

## Result

_(filled in when the ticket is implemented)_

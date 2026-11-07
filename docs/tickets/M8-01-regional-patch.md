# M8-01 Regional patch extraction and boundary conditions

- Status: draft
- Milestone: M8
- Depends on: M7-05
- Crates: kuniumi-sim (new `regional` module), kuniumi-world
- Design: 高解像度化と出力 / 3つの方法
- Decisions: none

## Goal

Cut a region of about 1,000 km square out of a finished world onto a fine local
grid, with everything the regional simulation needs from the global result.

## Scope

In:
- A local grid (regular in a Lambert azimuthal equal-area projection centred on
  the region, 100 m–1 km spacing) and its mapping to and from the sphere.
- Interpolated inputs: uplift rate, erodibility, precipitation and temperature
  (monthly), initial elevation.
- Boundary conditions from the global world: river inflow points with their
  discharge, outlet (coast) heights fixed, sea level.
- A `RegionSpec` (centre, size, spacing) in `Settings`.

## Implementation notes

To be refined: interpolation (barycentric on the Delaunay triangles), how many
global cells a region spans, memory limits for 1,000 × 1,000 km at 100 m.

## Acceptance criteria

- [ ] Projection round trips within 1e-9 m.
- [ ] Interpolated fields match the global values at the global cell centres.
- [ ] Every global river entering the region appears as an inflow with matching
      discharge.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-sim regional
```

## Result

_(filled in when the ticket is implemented)_

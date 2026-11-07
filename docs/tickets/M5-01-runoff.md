# M5-01 Evapotranspiration and runoff

- Status: draft
- Milestone: M5
- Depends on: M4-07, M3-05
- Crates: kuniumi-sim (stages/hydrology)
- Design: 水文（川・湖）
- Decisions: none

## Goal

Annual and monthly runoff (precipitation minus actual evapotranspiration) for
every land cell, and river discharge accumulated from it on the final terrain.

## Scope

In:
- Potential evapotranspiration from temperature (method chosen and documented
  here, e.g. Thornthwaite or Hargreaves-type with insolation).
- Annual actual evapotranspiration from the Budyko curve (aridity index
  PET/P), with monthly runoff distributed in proportion to the monthly water
  surplus (documented scheme).
- Discharge by accumulating runoff × area over the drainage graph of M3-01 on
  the final terrain (m³/s).
- Outputs: `runoff_annual`, `runoff_m01` … `runoff_m12` (m/s),
  `evapotranspiration_annual` (m/s), `discharge_annual` (m³/s).

Out:
- Lakes and closed basins (M5-02).

## Implementation notes

To be refined: the PET method, the Budyko parameter, units of the monthly
scheme. Sanity check from the design: 1 mm/yr of runoff over 1 km² is
1,000 m³/yr, so discharge (m³/s) ≈ runoff (mm/yr) × area (km²) / 31,600.

## Acceptance criteria

- [ ] Runoff is never negative annually and never exceeds precipitation.
- [ ] Water balance per cell: P − ET − runoff = 0 annually (1e-12 relative).
- [ ] An Amazon-like test basin (6 × 10⁶ km², 1,000 mm/yr runoff) gives about
      2 × 10⁵ m³/s at its outlet.
- [ ] Runoff maps viewed and described.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-sim hydrology
```

## Result

_(filled in when the ticket is implemented)_

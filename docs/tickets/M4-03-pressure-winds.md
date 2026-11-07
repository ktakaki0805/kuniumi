# M4-03 Pressure belts and surface winds

- Status: draft
- Milestone: M4
- Depends on: M4-02
- Crates: kuniumi-sim (stages/climate)
- Design: 気候モデル / 構成要素
- Decisions: D-12

## Goal

Monthly sea-level pressure and surface wind vectors, with circulation cells
whose number and width depend on the rotation rate, plus thermal lows and highs
over land.

## Scope

In:
- Zonal pressure pattern: equatorial low following the thermal equator,
  subtropical highs, subpolar lows, polar highs; the number and width of the
  cells scale with the rotation rate (scaling and source documented).
- Thermal anomalies from the land–sea temperature contrast (monsoon lows and
  highs).
- Surface wind from geostrophic balance, turned 20–30° towards low pressure by
  friction; near the equator the wind follows the pressure gradient directly.
- Outputs: monthly pressure (Pa) and wind east and north (m/s).

## Implementation notes

To be refined: functional form of the zonal profile, magnitudes, smoothing,
handling of the Coriolis parameter near the equator.

## Acceptance criteria

- [ ] Earth: easterly trade winds within about ±30°, westerlies at 35–60° in
      both hemispheres; the ITCZ shifts with the seasons.
- [ ] A ten times longer day widens the Hadley cells.
- [ ] Wind maps viewed and described.
- [ ] `mise run ci` passes.

## Verification

Wind maps for January and July rendered to `out/M4-03/` and viewed.

## Result

_(filled in when the ticket is implemented)_

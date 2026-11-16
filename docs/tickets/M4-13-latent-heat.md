# M4-13 Latent heat: evaporation paid by the energy balance

- Status: draft
- Milestone: M4
- Depends on: M4-06
- Crates: kuniumi-sim (stages/climate)
- Design: 気候モデル / 構成要素; 気候モデル / 計算の流れ
- Decisions: none

## Goal

Evaporation takes its latent heat from the surface's energy balance, and
condensation gives it back to the atmosphere, so that the hydrological cycle
of a hot planet is bounded by the energy it receives.

## Scope

In:
- The moisture model of M4-06 evaporates by the bulk formula and recycles a
  fixed share of the land's precipitation without any energy cost: on
  Earth's world its latent heat flux is about 83 W/m² (2.88 mm/day; Earth's
  is about 80–88 W/m²), but at 2,500 W/m² of sunlight (M4-06,
  `extreme_planets_stay_finite`) the planet rains 112 mm/day, about
  3,200 W/m² of latent heat against a few hundred absorbed; the hottest
  corner of the v1 range (2,500 W/m², e = 0.99, obliquity 90°, 8 h,
  0.1 bar; monthly means up to 12,029 K) 930 mm/day; and the 8 h, 0.1 bar
  planet 9 mm/day.
- The latent heat `L E` as a sink of the surface in the energy balance of
  M4-02 and M4-05, and `L P` released where it rains (or carried with the
  vapour), with the energy and water budgets both conserved; or, if the
  coupling is too large a change, an energy limit on evaporation (for
  example the Priestley–Taylor or equilibrium evaporation of the absorbed
  energy).
- The effect on Earth's calibration (`A`, `D` of M4-05, `U_g` of M4-06): the
  energy balance's diffusion `D` already carries the moist heat transport.

Out:
- Soil water and evapotranspiration of land (M5-01).

## Implementation notes

To be refined, and decided with the owner (it couples the climate's parts
and changes their calibration): which form of the coupling, its iteration
with the loop of M4-07, the recalibration, the extreme planets.

## Acceptance criteria

- [ ] Energy and water conserved per sub-step.
- [ ] The latent heat flux of the 2,500 W/m² planet and of the hottest
      corner below their absorbed sunlight;
      Earth's criteria of M4-02, M4-05 and M4-06 still met.
- [ ] `mise run ci` passes.

## Verification

Precipitation and temperature maps of Earth and of the hot planet rendered to
`out/M4-13/`.

## Result

_(filled in when the ticket is implemented)_

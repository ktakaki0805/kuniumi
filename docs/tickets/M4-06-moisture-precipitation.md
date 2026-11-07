# M4-06 Moisture, precipitation and rain shadows

- Status: draft
- Milestone: M4
- Depends on: M4-05
- Crates: kuniumi-sim (stages/climate)
- Design: 気候モデル / 構成要素
- Decisions: none

## Goal

Monthly precipitation (m/s) from moisture evaporated over the oceans, carried by
the winds and rained out by saturation, convergence and orographic lift, with
rain shadows downwind of mountains.

## Scope

In:
- Column water vapour advected by the surface wind (semi-Lagrangian or upwind,
  deterministic).
- Evaporation from the ocean (bulk formula with sea-surface temperature and
  wind speed) and recycling over land (a fraction of precipitation
  re-evaporates).
- Precipitation where vapour exceeds a fraction of saturation
  (Clausius–Clapeyron with `det::exp`), plus convergence and orographic terms
  (`∝ max(0, u·∇h)`), removing the vapour that rains out.
- Outputs: monthly precipitation and vapour.

## Implementation notes

To be refined: constants, sub-stepping, the orographic coefficient. The linear
model of Smith & Barstad (2004) is reserved for the regional stage (M8-02).

## Acceptance criteria

- [ ] Water is conserved: evaporation − precipitation − storage change sums to
      about zero globally per year (documented tolerance).
- [ ] Earth-like: wet equator, dry subtropics, wet mid-latitude west coasts and
      clear rain shadows behind mountain ranges.
- [ ] Precipitation maps viewed and described.
- [ ] `mise run ci` passes.

## Verification

Precipitation maps for January and July rendered to `out/M4-06/` and viewed.

## Result

_(filled in when the ticket is implemented)_

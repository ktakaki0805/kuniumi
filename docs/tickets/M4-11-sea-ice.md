# M4-11 Sea ice in the energy balance

- Status: draft
- Milestone: M4
- Depends on: M4-05
- Crates: kuniumi-sim (stages/climate)
- Design: 気候モデル / 構成要素
- Decisions: none

## Goal

A frozen ocean insulates its mixed layer, so polar oceans get the large
seasonal range of the real Arctic instead of that of open water.

## Scope

In:
- In M4-02 sea ice only raises the albedo: the cell keeps the heat capacity
  of a 75 m mixed layer, so a frozen ocean swings by a few kelvin a year
  (the Arctic swings by about 30 K).
- A smaller effective heat capacity under sea ice (North, Mengel & Short
  1983 used about 9 times that of land, as Williams & Kasting 1997 report),
  or a thermodynamic sea-ice layer, with energy conserved through freezing
  and melting.

Out:
- Sea-ice motion.

## Implementation notes

To be refined: the formulation (switching heat capacity versus an ice
thickness with latent heat), its source, the coupling with the sea-surface
temperature of M4-05.

## Acceptance criteria

- [ ] Earth: the seasonal range of the ice-covered polar ocean within a
      documented range of observations.
- [ ] Energy (including latent heat) conserved per step.
- [ ] `mise run ci` passes.

## Verification

Seasonal-range maps rendered to `out/M4-11/`.

## Result

_(filled in when the ticket is implemented)_

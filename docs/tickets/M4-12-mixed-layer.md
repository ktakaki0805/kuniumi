# M4-12 Sea-surface fronts: a mixed layer under the atmosphere

- Status: draft
- Milestone: M4
- Depends on: M4-05
- Crates: kuniumi-sim (stages/climate)
- Design: 気候モデル / 構成要素
- Decisions: none

## Goal

The sea surface keeps the sharp patterns of its currents and upwelling: cold
coastal bands along the eastern shores of the oceans and warm boundary
currents, with subtropical east–west contrasts of several kelvin as on Earth.

## Scope

In:
- M4-05 adds the heat of the currents and of coastal upwelling to the
  single-layer energy balance of M4-02, whose diffusion (the atmosphere's
  eddies) acts on the sea surface at the grid's scale: the anomalies spread
  over the diffusion length (about 3,000 km), so the subtropical west coasts
  of the Earth-parameter world are only 0.8–1.3 K colder than the east coasts
  (observed 4–8 K), and the cooling of a coastal cell shows as broad cool
  water (M4-05, Limitations).
- A mixed-layer temperature of its own under the atmosphere's, coupled by an
  air–sea heat flux, advected by the currents and cooled by upwelling, with
  the atmosphere diffusing only its own temperature; energy conserved.
- Possibly with it: the Ekman overturning's heat transport and equatorial
  upwelling (the cold tongue), which the barotropic currents of M4-04 lack.

Out:
- Sea ice (M4-11).

## Implementation notes

To be refined, and decided with the owner (it changes the structure of the
energy balance of M4-02 and its calibration): the air–sea coupling and its
source, the heat capacities of the two layers, the recalibration of `A` and
`D`, the effect on the seasonal cycle (M4-02's criterion 2), the cost.

## Acceptance criteria

- [ ] Earth: the subtropical west coasts colder than the east coasts by a
      documented margin closer to the observed 4–8 K; the criteria of M4-02
      and M4-05 still met.
- [ ] Energy conserved per sub-step.
- [ ] `mise run ci` passes.

## Verification

Sea-surface temperature maps rendered to `out/M4-12/` and compared with those
of `out/M4-05/`.

## Result

_(filled in when the ticket is implemented)_

# M4-05 Sea-surface temperature transport and upwelling

- Status: draft
- Milestone: M4
- Depends on: M4-04
- Crates: kuniumi-sim (stages/climate)
- Design: 気候モデル / 構成要素
- Decisions: none

## Goal

Ocean currents move heat: warm western boundary currents warm high latitudes,
and equatorward eastern boundary currents with coastal upwelling cool the
subtropical west coasts of continents.

## Scope

In:
- Advection of sea-surface temperature by the currents (upwind, deterministic),
  coupled to the energy balance of M4-02 as a source term.
- Coastal upwelling where alongshore wind drives offshore Ekman transport,
  cooling the sea surface.
- Output: monthly sea-surface temperature (K).

## Implementation notes

To be refined: coupling to M4-02 (operator splitting), upwelling strength.

## Acceptance criteria

- [ ] Earth-like: west coasts at 15–30° latitude have colder sea surfaces than
      east coasts at the same latitude; at high latitudes the eastern sides of
      oceans are warmer than the western sides.
- [ ] The advection step conserves heat in a closed basin.
- [ ] `mise run ci` passes.

## Verification

Sea-surface temperature maps rendered to `out/M4-05/` and viewed.

## Result

_(filled in when the ticket is implemented)_

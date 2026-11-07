# M8-03 Projections and game-facing exports

- Status: draft
- Milestone: M8
- Depends on: M8-02
- Crates: kuniumi-io, kuniumi-cli
- Design: 高解像度化と出力 / ゲームへの出力
- Decisions: D-07

## Goal

Exports a game can load directly: world and regional maps in suitable
projections, the layers the design lists, and optionally a hex board.

## Scope

In (to be confirmed with the owner before refining; the game engine and map use
are an open question in `docs/decisions.md`):
- Projections: equirectangular and an equal-area world projection; Lambert
  azimuthal equal-area for regions; optionally a Goldberg polyhedron hex board.
- Layers: elevation (16-bit or more), land/sea, rivers (polylines with
  discharge), lakes (fresh/salt), monthly temperature and precipitation,
  Köppen, biome, NPP, plates and rock types.
- Optional resource hints from plate history (arc metals, craton iron and gold,
  sedimentary-basin oil and coal) as a layer.
- `docs/roadmap.md`: M8 marked as met (with M8-02).

## Implementation notes

Ask the owner which engine and which map uses to target before refining this
ticket; the answer may also add engine bindings (D-07) as new tickets.

## Acceptance criteria

- [ ] To be written when the ticket is refined.
- [ ] `mise run ci` passes.

## Verification

To be written when the ticket is refined.

## Result

_(filled in when the ticket is implemented)_

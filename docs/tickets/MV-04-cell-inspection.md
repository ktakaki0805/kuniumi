# MV-04 Cell inspection in the viewer

- Status: draft
- Milestone: MV
- Depends on: MV-02
- Crates: kuniumi-viewer
- Design: 実装ロードマップ / 技術スタック
- Decisions: D-10, D-19

## Goal

Hovering over or clicking the globe of `kuniumi-viewer` shows the cell under
the cursor: its id, latitude and longitude, and the values of all fields,
in the window itself.

## Scope

In (split from MV-02 when it was refined):
- Picking: the ray through the cursor, intersected with the drawn surface
  (relief included), and the cell found with `SphereGrid::nearest` through a
  `Locator`.
- An on-screen panel with the cell id, latitude and longitude (degrees at
  the boundary), every field's value with its unit (class fields as their
  id), and the current layer with its legend (colour map and range).
- The same text printed to the terminal on click.

Out:
- Browser build (MV-03); editing worlds.

## Implementation notes

To be refined: how text is drawn — winit and wgpu draw none. Options: an
8×8 bitmap font ported from a public-domain source (recorded in
`THIRD_PARTY.md`) drawn as textured quads, no new crate; or a text or GUI
crate (`glyphon`, `egui`), which needs the owner's approval under D-10 and
may tie the `wgpu` version. Whether the picking ray accounts for the relief
exactly (iterate: intersect the sphere of the cell's drawn radius, find the
cell again) or by rendering cell ids into an offscreen buffer.

## Acceptance criteria

- [ ] Picking finds the cell under the cursor for any view and relief
      (tests: rays through known cells at several exaggerations).
- [ ] The panel shows the cell id, latitude, longitude and every field's
      value with its unit; a screenshot with the panel is saved and
      described.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo run --release -p kuniumi-viewer -- out/MV-02/650k
```

## Result

_(filled in when the ticket is implemented)_

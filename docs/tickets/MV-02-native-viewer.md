# MV-02 Native globe viewer

- Status: draft
- Milestone: MV
- Depends on: MV-01, M2-06
- Crates: new `kuniumi-viewer`
- Design: 実装ロードマップ / 技術スタック
- Decisions: D-10, D-19

## Goal

`kuniumi-viewer <world dir>` opens a window with a rotatable 3D globe showing
any field of the world, with relief, sea level and per-cell inspection.

## Scope

In:
- New crate `kuniumi-viewer` using `wgpu` and `winit`; versions listed in the
  D-10 table; the whole dependency tree passes cargo-deny (extend the licence
  list only with a reason).
- Mesh: each Voronoi cell as a triangle fan from its centre to its Voronoi
  vertices, coloured per cell with the colour maps of `kuniumi-io`.
- Orbit camera (drag to rotate, wheel to zoom), layer selection (keyboard is
  enough), adjustable relief exaggeration, sea-level surface.
- Hover or click shows the cell id, latitude and longitude, and the values of
  all fields (picking through `SphereGrid::nearest`).
- `f32` for GPU data with a crate-level `#![allow(clippy::disallowed_types)]`,
  permitted in this crate only; nothing computed here feeds back into
  generation.

Out:
- Browser build (MV-03); editing worlds; animation of time steps.

## Implementation notes

To be refined: `wgpu`/`winit` versions; how the viewer fits `mise run ci` —
its build time, and the minimum Rust version of `wgpu`, which may be newer than
the MSRV of D-24 (exclude the viewer from the `msrv` task or raise the MSRV;
ask the owner); shader structure; handling about 4 million triangles for
650,000 cells.

## Acceptance criteria

- [ ] A 650,000-cell world renders at interactive frame rates on the owner's
      machine (measured and recorded).
- [ ] Every field of the world can be displayed; class fields use the
      categorical palette.
- [ ] Screenshots of three layers saved and described.
- [ ] `cargo tree -p kuniumi-sim` shows no `wgpu` or `winit`.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo run --release -p kuniumi-viewer -- out/M2-06
```

## Result

_(filled in when the ticket is implemented)_

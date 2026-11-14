# MV-02 Native globe viewer

- Status: done
- Milestone: MV
- Depends on: MV-01, M2-06
- Crates: new `kuniumi-viewer`
- Design: 実装ロードマップ / 技術スタック
- Decisions: D-10, D-19, D-24

## Goal

`kuniumi-viewer <world dir>` opens a window with a rotatable 3D globe showing
any field of the world, with relief and the sea level. (Per-cell inspection
moved to MV-04 when this ticket was refined.)

## Scope

In:
- New crate `kuniumi-viewer` (a binary) using `wgpu` and `winit`, already in
  the D-10 table; the whole dependency tree passes cargo-deny with the current
  licence list.
- Mesh: each Voronoi cell as a triangle fan from its centre to its Voronoi
  vertices, coloured per cell with the colour maps of `kuniumi-io` (the same
  colours as the world directory's PNG maps).
- Orbit camera (drag to rotate, wheel to zoom), layer selection with the
  keyboard, adjustable relief exaggeration, a sea-level surface.
- Headless screenshots to PNG and a frame-time benchmark, for the evidence.
- `f32` for GPU data with a crate-level `#![allow(clippy::disallowed_types)]`,
  permitted in this crate only; nothing computed here feeds back into
  generation.

Out:
- Hover or click inspection of a cell (MV-04); browser build (MV-03); editing
  worlds; animation of time steps.

## Implementation notes

### Size and split (approved by the owner)

The whole draft scope came to an estimated 1,900 changed lines: about 1,100
for the globe (mesh, camera, layers, relief, sea, GPU code, screenshots,
benchmark, tests) and 500 for inspection, most of it an on-screen text
overlay (winit and wgpu draw no text; MV-03's browser build has no terminal
either), plus documents. The split keeps this ticket to the globe and moves
inspection to the new draft MV-04, on which MV-03 then depends as well.

### Dependencies and the toolchain

- `wgpu = "30.0.2"` (default features) and `winit = "0.30.13"` with
  `default-features = false, features = ["rwh_06", "x11", "wayland",
  "wayland-dlopen", "wayland-csd-adwaita-notitle"]`. winit's default
  `wayland-csd-adwaita` draws the title text of client-side decorations with
  `ab_glyph`, whose `ttf-parser` cargo-deny reports as unmaintained
  (RUSTSEC advisory); the `notitle` variant keeps the decorations without
  the font stack. Measured in a scratch crate: 250 packages, cargo-deny
  passes with the current `deny.toml` (advisories, bans, licences and
  sources ok; duplicate-version warnings only), and a cold build takes 15 s
  in debug, 20 s in release and 11 s for clippy on the owner's machine.
  Neither crate is on the generation path, so neither is pinned to an exact
  version (D-10 rule 2).
- `wgpu` 30 needs Rust 1.87 (`wgpu`, `wgpu-core`, `wgpu-hal`, `naga` all
  declare 1.87); the workspace MSRV is 1.85 (D-24). Chosen by the owner
  when the plan was approved: the viewer declares its own
  `rust-version = "1.87"`, the `msrv` task checks `--workspace --exclude kuniumi-viewer` (checked in
  a scratch copy of the workspace: passes with Rust 1.85.0 and the lock
  file that includes the viewer), and D-24 records that the viewer's
  minimum version follows `wgpu`. The library crates and the CLI keep 1.85.
- No other crate: arguments are parsed by hand (a handful of flags), blocking
  on wgpu's futures uses `std::task::Waker::noop` (stable since 1.85), PNG
  files go through `kuniumi_io::png`, vertex data is converted to bytes with
  `to_ne_bytes`.

### Crate layout (`crates/kuniumi-viewer/src/`)

- `main.rs`: `#![allow(clippy::disallowed_types)]` with the reason; parses
  the arguments and runs one of three modes: interactive window, headless
  `--screenshot`, `--benchmark` (in the window or offscreen).
- `args.rs`: `args::parse(impl IntoIterator<Item = String>) ->
  Result<Parsed, String>` (`Parsed::Run(Args)` or `Parsed::Help`) for
  `kuniumi-viewer <world dir> [--layer <name>] [--relief <factor>]
  [--no-sea] [--lat-deg <φ>] [--lon-deg <λ>] [--distance <radii>]
  [--size <w>x<h>] [--screenshot <png> | --benchmark <frames>
  [--offscreen]]`; degrees are converted with `kuniumi_math::units` at this
  boundary.
- `mesh.rs`: `GlobeMesh::new(grid: &SphereGrid, height_above_sea_m:
  &[f64]) -> GlobeMesh` with `vertices: Vec<[f32; 4]>` (unit direction and
  height above the sea level in metres), `indices: Vec<u32>` and
  `cell_vertex_offsets: Vec<u32>`; `GlobeMesh::vertex_colours(&self,
  cell_colours: &[Rgb8]) -> Vec<[u8; 4]>`. Cell `v` with ring `n_0 …
  n_{k−1}` gets vertex `c_v` (its centre, height `h_v`) and one vertex per
  Delaunay triangle `t_m = (v, n_m, n_{m+1})`: its circumcentre
  `w/|w|`, `w = (b − a) × (c − a)` (the Voronoi vertex, as in
  `kuniumi_sphere::voronoi`), at the height `(h_a + h_b + h_c)/3`, both
  computed in `f64` on the triangle rotated so that its smallest index comes
  first, so the three cells that share a Voronoi vertex get the same bits
  and the surface has no cracks; triangles `(c_v, t_m, t_{m+1})`,
  counter-clockwise from outside. For `n` cells: `7n − 12` vertices and
  `6n − 12` triangles (650,000 cells: 4.55 million vertices, 3.9 million
  triangles, 73 MB of positions, 18 MB of colours, 47 MB of indices; each
  buffer below wgpu's default `max_buffer_size` of 256 MiB, and the
  adapter's limit is checked with a clear error). The heights are those of
  `surface_height` minus the sea level, or zero for a world without it.
- `layers.rs`: one layer per field in the world's order; the colours of a
  layer are `world_dir::map_colours(field, world, metrics)` with
  `ElevationMetrics::of_world`, exactly the colours of `<name>.png`; class
  fields therefore use the categorical palette. `Layers::next`/`previous`
  cycle with wrap-around; `Legend` (colour map and range) is printed to the
  terminal on every change.
- `camera.rs`: `OrbitCamera { lat_rad, lon_rad, distance }` (distance from
  the planet's centre in planet radii, the planet drawn with radius 1),
  looking at the centre with +Z (north) up; `view_proj(aspect) -> [[f32;
  4]; 4]` (right-handed look-at and a perspective of 45° vertical field of
  view with depth in [0, 1], computed in `f64` with `kuniumi_math::det` and
  converted at the end); `drag(dx_px, dy_px, height_px)` rotates by an angle
  proportional to the height above the surface, latitude clamped to ±89°;
  `zoom(steps, min_distance)` scales the height above the surface by
  `1.1^−steps`, clamped between just above the highest drawn relief and 10
  radii.
- `gpu.rs` and `shader.wgsl`: one vertex buffer of positions and heights,
  one of colours (`Unorm8x4`, rewritten on a layer change), one index buffer,
  one uniform buffer (view-projection, relief scale `e/R` with the
  planet radius `R` in metres, light direction, sea colour). The vertex
  shader places a vertex at radius `1 + e·h/R`; the sea pass draws the same
  mesh with `h = 0`, translucent blue, depth test without depth writes, so
  every terrain point below the sea level lies radially below the sea
  surface (the terrain triangle's point is `S·q` with `S ≤ 1` and `q` on
  the sea triangle); at `e = 0` the sea pass is skipped. Fragment shading:
  the colour decoded from sRGB, Lambert from a light above-left of the
  camera plus ambient, from the face normal `cross(dpdx(p), dpdy(p))`;
  depth `Depth32Float`, 4× MSAA, an sRGB target. Screenshots render into an
  offscreen texture of the same pipeline and are read back row by row.
- Keys: Left/Right previous/next layer, Up/Down relief steps of 0, 1, 2, 5,
  10, 20, 50, 100 (default 20), S sea on/off, R reset view, Esc or Q quit.
  The window title names the world, the layer and its unit, the relief and
  the sea; the terminal prints the legend of each layer shown; redraws
  happen on input (`ControlFlow::Wait`), every frame only in the benchmark.
- Benchmark: `--benchmark N` spins the globe by 0.5° of longitude per frame
  for `N` frames and prints the mean, median, 95th percentile and maximum
  time (`std::time::Instant`, allowed here: the viewer is not generation)
  between frames presented in the window (`Fifo`), or with `--offscreen`
  the time to draw each frame without a window until the GPU has finished
  it; the load time (read, mesh and layers) is printed at start.

### Tests (no GPU needed, so `mise run ci` runs them anywhere)

- `mesh`: vertex and index counts `7n − 12` and `3(6n − 12)` on the
  octahedron and Fibonacci grids of 1,000 and 20,000 cells; every Voronoi
  vertex equidistant from its three centres (in `f64` before rounding,
  within 1e-12) and its three copies bit-identical; every triangle
  counter-clockwise from outside; the solid angles of all triangles (from
  the `f32` vertices, in `f64`) sum to 4π within 1e-9, so the fans cover
  the sphere once; heights: centres equal the cell's, Voronoi vertices the
  mean of three (by hand on the octahedron); `vertex_colours` repeats each
  cell's colour over its vertex range.
- `layers`: on a generated 1,000-cell world (stages grid, tectonics, crust),
  every field is a layer, each layer's colours equal
  `map_colours(...).colours(values)` cell by cell, class fields use
  `Colormap::Categorical` and `surface_height` the elevation map at the sea
  level; cycling visits every layer and wraps.
- `camera`: the view-projection maps the centre of the planet to the centre
  of the screen and the point under the camera in front of it; north is up
  at latitude 0; a drag to the right turns the globe to the right; the
  latitude and zoom clamps hold; values checked by hand.
- `args`: every flag, defaults, and the error messages of bad input.
- `shader`: `shader.wgsl` parses and validates with `wgpu::naga` (front end
  and validator), so a broken shader fails `cargo test` without a GPU.

## Acceptance criteria

- [x] A 650,000-cell world (`kuniumi generate --seed 1 --cells 650000 --out
      out/MV-02/650k`) renders at interactive frame rates on the owner's
      machine: `kuniumi-viewer out/MV-02/650k --benchmark 600` (relief and
      sea on) reports a mean frame time of at most 33 ms with vsync; the
      uncapped rate (`--no-vsync`), the load time and the machine are
      recorded. (Approved wording. The uncapped rate is measured with
      `--offscreen` instead of `--no-vsync`; see Deviations in the Result.)
- [x] Every field of the world can be displayed and cycled with the
      keyboard, in the colours of its PNG map; class fields use the
      categorical palette (tests `layers`).
- [x] The mesh covers the sphere once without cracks (tests `mesh`); the
      camera tests pass.
- [x] Screenshots of three layers (`surface_height` with relief and sea,
      `plate_id`, `crust_age`) saved under `out/MV-02/` with `--screenshot`
      and described.
- [x] `cargo tree -p kuniumi-sim`, `-p kuniumi-io` and `-p kuniumi-cli` show
      no `wgpu` or `winit`.
- [x] `mise run ci` passes.

## Verification

```sh
cargo run --release -p kuniumi-cli -- generate --seed 1 --cells 650000 --out out/MV-02/650k
cargo run --release -p kuniumi-viewer -- out/MV-02/650k
cargo run --release -p kuniumi-viewer -- out/MV-02/650k --benchmark 600
cargo run --release -p kuniumi-viewer -- out/MV-02/650k --benchmark 600 --offscreen
cargo run --release -p kuniumi-viewer -- out/MV-02/650k --layer plate_id --screenshot out/MV-02/plate_id.png
cargo tree -p kuniumi-sim | grep -E "wgpu|winit"   # nothing
```

## Result

What was done:
- New crate `kuniumi-viewer` (binary `kuniumi-viewer`, `rust-version =
  "1.87"`): `main.rs` (the crate-level `f32` allowance with its reason,
  `Scene`, the three modes), `args.rs`, `mesh.rs`, `layers.rs`,
  `camera.rs`, `view.rs`, `gpu.rs` (`Renderer`, `Offscreen`,
  `request_device`, `block_on`), `shader.wgsl`, `app.rs` (the winit event
  loop, `press` for the keys, `FrameStats`). Dependencies `wgpu` 30.0.2 and
  `winit` 0.30.13 (without the Wayland title font); `Cargo.lock` gains 232
  packages and changes none that was there.
- `mise.toml`: the `msrv` task excludes `kuniumi-viewer`. `docs/decisions.md`:
  D-06 (the crate moves into the table), D-24 (the viewer's exception).
  `docs/design-ja.md` (プロジェクト方針, 言語・構成: `viewer`), `README.md`
  (layout and usage), `CLAUDE.md` (workspace block), `CHANGELOG.md` (Added,
  Changed). Split: MV-04 (draft) and MV-03's dependency, roadmap row.

Deviations from the plan:
- `--no-vsync` became `--offscreen`. The window offers `Mailbox` and `Fifo`
  on this machine, but on Wayland winit paces redraws by the compositor's
  frame callbacks, so a window never exceeds the display's 60 Hz: with
  `AutoNoVsync` the benchmark still measured 16.3–16.7 ms. The drawing
  time without that pacing is measured offscreen instead: each frame drawn
  into a texture and timed until the GPU has finished it.
- The window title does not show the last frame time (the benchmark prints
  frame times); on GNOME/Wayland the title is not drawn in the window
  anyway (`wayland-csd-adwaita-notitle`), only in the overview.
- Exact float comparisons in the tests use `to_bits` (clippy's
  `float_cmp`), as the other crates' tests do.
- Added beyond the plan: `press` (the key handling as a function of the
  scene and the view) with a test, so the controls are checked without a
  window; the mesh tolerance for 4π is 1e-9 rather than 1e-5. After review
  round 1: `--relief` lies in [0, 100] (the largest key step; at ×1000
  Earth's trenches are drawn through the centre, inside out, seen in a
  trial screenshot), and the camera's largest distance is 10 radii of the
  highest drawn point instead of 10 planet radii, so no relief empties the
  distance range (R1-1).

Evidence:

1. Frame rate (this machine: NVIDIA GeForce RTX 3060 Ti, Vulkan, Ubuntu
   GNOME on Wayland, 60 Hz display; release build; world
   `out/MV-02/650k`, `kuniumi generate --seed 1 --cells 650000`, 650,000
   cells, 3,899,988 triangles, relief ×20, sea on):
   - `--benchmark 600` in the window, 1280×800: mean 16.74 and 16.77 ms
     (59.7 and 59.6 fps), median 16.66 ms, p95 16.92 and 16.87 ms, max
     33.40 and 33.36 ms (one missed refresh); at 1853×1013 (the largest
     window the screen allows): mean 16.76 ms, p95 16.87 ms. Interactive,
     at the display's rate.
   - `--benchmark 600 --offscreen`, 1280×800: mean 0.99–1.00 ms (about
     1,000 fps), median 0.98–1.00 ms, p95 1.02–1.07 ms, max 2.14–2.19 ms
     (three runs); 2560×1600: mean 1.12 ms; relief 0 and no sea (one pass):
     0.57 ms.
   - Load: `world_dir::read` 0.21–0.22 s, mesh and layers 0.12–0.13 s
     (three runs); a whole `--screenshot` run, including GPU setup, upload
     and read-back, 1.5–1.6 s and 609 MB peak resident memory.
   - After the round-1 fixes: window mean 16.75 ms, offscreen mean
     1.00 ms, read 0.21 s, mesh and layers 0.12 s.
2. Layers (`layers` tests, 5): on a generated 1,000-cell world (seed 3,
   stages grid, tectonics, crust, 13 fields) every field is a layer;
   for each, `Layers::colours` equals `map_colours(field, world,
   ElevationMetrics::of_world(world)).colours(values, range)` cell by cell
   and the layer's bounds equal the resolved range; class fields
   (`neighbor_count`, `plate_id`, `crust_type`, …) use
   `Colormap::Categorical`; `surface_height` uses the elevation map at the
   world's sea level from 11 km below to 9 km above it (the axis of
   `hypsometry.png`); other quantities viridis. Cycling visits every
   layer and wraps; an unknown `--layer` lists the fields; a world without
   `surface_height` starts on its first field, one without fields is grey;
   `metrics` follows `world_dir::write`'s checks (NaN height, too many
   bins). Keys (`app::tests::keys`): Right, Left, Up (20 → 50), Down
   (→ 20 → 10), S twice, R after a drag and zoom restores the camera,
   other keys do nothing, Q and Esc quit, and a higher relief pushes the
   camera out. In a screenshot of `plate_id` without the sea
   (`out/MV-02/650k_plate_id_nosea.png`), a throwaway Python check (PNG
   decoded with the standard library, sRGB to linear) finds 99.06 % of the
   360,244 globe pixels to be a shade in [0.35, 1] (the shader's ambient
   to full light) of one ColorBrewer Paired colour within two 8-bit levels
   (the round-2 reviewer's own decoder, within one level: 98.87 %);
   the rest are multisampled edges between two plates or against the
   background. Interactive use (drag, wheel, keys in the window) could not
   be driven from this session (no input tool); the window was opened and
   presented 600 frames in each window benchmark.
3. Mesh and camera: `mesh` tests (6): `7n − 12` vertices and `3(6n − 12)`
   indices on the octahedron and Fibonacci grids of 1,000 and 20,000
   cells; every Voronoi vertex unit within 1e-15 and equidistant from its
   three centres within 1e-12, the same `f64` bits from each of its three
   cells and the same `f32` bits in the mesh (2n − 4 triangles, three
   copies each); every triangle counter-clockwise from outside; the solid
   angles of all triangles sum to 4π within 1e-9; the octahedron by hand
   (face centres `(±1, ±1, ±1)/√3` at the mean height); non-finite heights
   drawn at 0; `vertex_colours` per cell range. A mutation that computes
   each Voronoi vertex in its own cell's rotation fails
   `voronoi_vertices_are_equidistant_and_shared_bit_for_bit` (it passed
   before the test checked the `f64` bits of all three rotations with
   heights that are not exact in binary). `camera` tests (4): values by
   hand (centre at the screen centre, the point below at depth 2/3, east
   at x = (1 + √2)/3, north up), drag direction and speed, latitude clamp
   at ±89°, zoom by 1.1 per step and its clamps, also for a relief 20
   radii high; `tests::scene_heights_lie_above_the_sea_level` (main.rs,
   review R1-2): every centre vertex holds `(surface_height − sea level)`
   as `f32` bit for bit, the top is the largest of them, the layer
   choice, a radius that is not finite and positive refused (R2-2), and a
   world without `surface_height` is flat; `view` tests (2): the
   relief steps, `e/R`, the top radius, no sea pass at `e = 0` or without
   heights; `gpu` tests (3): `shader.wgsl` parses and validates with
   `wgpu::naga`, its uniform block is the 112 bytes written, the byte
   layout (matrix columns, light, sea colour and opacity, relief scale),
   sRGB decoding; `args` (3) and `FrameStats` (1). Test-first:
   against stubs (commit "viewer crate with the tests … against stubs"),
   17 of the 24 tests then present failed (the argument parser, the shader
   check, the uniform layout and the frame statistics were implemented).
4. Screenshots (`--screenshot`, 1280×800, camera at latitude 20°,
   longitude 0°, 3 radii from the centre unless stated), viewed:
   - `out/MV-02/650k_surface_height.png` (relief ×20, sea on): the globe
     on the dark grey background, lit from the upper left; three continents
     in green (west, east, south) and a fourth at the right limb, a
     yellow-green collision belt across the lower left joining
     the western continent to the southern one, another on the northern
     coast of the eastern continent; light-blue shelves around every coast,
     deep blue ocean under the translucent sea, a darker band at the limb.
   - `out/MV-02/650k_surface_height_close.png` (relief ×50, latitude 5°,
     longitude −25°, 1.35 radii): the collision plateau in yellow-green
     with shaded facets, a steep dark wall down into the trench beside it,
     the shelf in light blue at the coast of the green lowland: the relief
     reads as 3D. `out/MV-02/650k_relief100_close.png`, the same view at
     ×100 (the largest relief): taller walls, the same shapes, nothing
     inverted.
   - `out/MV-02/650k_plate_id.png` (sea on): the plates in the categorical
     palette (orange, brown, blue, green, light green, light blue, yellow
     at the limb); where the sea covers a plate its colour is tinted blue
     (the orange plate's ocean part shows as mauve), the continents keep
     their colours and their raised edges catch the light;
     `out/MV-02/650k_plate_id_nosea.png`: the same plates untinted, the
     orange plate covering the western continent and the ocean north of
     it.
   - `out/MV-02/650k_crust_age.png` (sea on): continents yellow (the
     largest age, 180 Myr), the ocean floor in viridis from dark purple
     (young, south-east) to teal and green (old, between the continents),
     the margins tan where the shallow sea covers yellow crust.
   - The three layer screenshots, rendered again after the round-1 fixes:
     byte-identical.
   - `out/MV-02/20k_surface.png` (20,000 cells): the same style with the
     hexagonal cells visible.
5. `cargo tree -p <crate> -e all` shows no line with `wgpu` or `winit` for
   `kuniumi-sim`, `kuniumi-io`, `kuniumi-cli`, `kuniumi-math`,
   `kuniumi-sphere`, `kuniumi-world` and `kuniumi-config`.
6. `mise run ci` passes (log `out/MV-02/ci.log`): fmt, clippy `-D
   warnings`, doc, 1,133 tests passed in debug and release together (10
   ignored: the existing timing tests), MSRV 1.85.0 without the viewer,
   cargo-deny (advisories, bans, licences, sources ok; duplicate-version
   warnings only, as allowed by `deny.toml`).

Review:
- Round 1 (APPROVE): the reviewer re-ran `mise run ci`, regenerated the
  650,000-cell world, the benchmarks (window 16.70 ms, offscreen 0.99 ms)
  and the screenshots (byte-identical, plus an odd size and full-screen
  renders without background pixels, so no cracks), the stub commit (17 of
  24 failing) and mutations. R1-1 minor, `--relief 20000` panicked in
  `f64::clamp` (top above 10 radii) (fixed: `--relief` in [0, 100], the
  distance range scales with the top, test). R1-2 minor, no test checked
  that the mesh heights are `surface_height − sea level` (fixed: a test of
  `Scene::new`). R1-3 nit, the sea colour in the uniform block was
  untested (fixed). R1-4 minor, the first criterion had been reworded after
  approval (fixed: the approved wording restored with a note pointing to
  the deviation; to be relayed to the owner when the ticket lands). R1-5 nit, the viewer's
  `rust-version = "1.87"` is not checked by any task (D-24 now says so).
- Round 2 (APPROVE): the reviewer re-ran `mise run ci`, built the viewer
  from scratch in a copy, re-ran the benchmarks (window 16.72 ms,
  offscreen 1.00 ms), re-rendered all six 650,000-cell screenshots
  (byte-identical) and its own palette check (98.87 % within one 8-bit
  level), and confirmed R1-1 to R1-5 resolved with mutations. R2-1 minor,
  no test accepted `--relief 100`, and `MAX_RELIEF` repeated the last
  relief step by hand (fixed: derived from `RELIEF_STEPS`, test; the
  half-open range now fails `every_flag`). R2-2 nit, a hand-edited
  `world.json` with a negative radius made the new distance clamp panic
  (fixed: `Scene::new` refuses a radius that is not finite and positive,
  test; the edited 20,000-cell world now exits with that message). R2-3
  nit, the Result said R1-4 was relayed before it was (fixed: "to be
  relayed when the ticket lands").
- Round 3 (APPROVE): the reviewer re-ran `mise run ci`, regenerated the
  650,000-cell world (all 30 files identical), the benchmarks (window
  16.66 ms, offscreen 0.99–1.02 ms), all seven screenshots
  (byte-identical), its palette check (98.88 %), renders at the pole and
  across the antimeridian (no cracks), and checked that the viewer adds no
  feature to the crates on the generation path (libm, robust and the I/O
  crates keep the same features). R2-1 to R2-3 resolved. R3-1 minor: the
  ×100 cap holds for Earth but not for the smallest valid planets: at
  1,500 km with Earth's ocean volume the sea floor lies 50.7 km deep and
  the drawn radius falls below 0 above ×29.6, so ×50 and ×100 turn the
  globe inside out (no crash; Down recovers); the `MAX_RELIEF` comment also
  calls Earth's trenches 7 km deep (about 11 km). Recorded in the draft
  follow-up MV-05 rather than fixed after the last round, so that no
  unreviewed code lands.

Follow-ups: MV-04 (draft, cell inspection, split from this ticket); MV-05
(draft, the relief exaggeration bounded by each world's deepest point,
review R3-1). The roadmap's MV criterion is not met yet: inspection
(MV-04) and the browser (MV-03) remain.

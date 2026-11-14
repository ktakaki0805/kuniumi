# MV-04 Cell picking in the viewer

- Status: done
- Milestone: MV
- Depends on: MV-02
- Crates: kuniumi-viewer
- Design: 実装ロードマップ / 技術スタック
- Decisions: D-10, D-19

## Goal

Hovering over the globe of `kuniumi-viewer` finds the cell under the cursor
exactly as drawn (relief included) and highlights it; clicking prints the
cell's id, latitude, longitude, the value of every field with its unit and
the legend of the layer shown to the terminal. (The on-screen panel with the
same text moves to MV-06; see "Size and split".)

## Scope

In:
- Picking: the ray through the cursor, intersected with the drawn triangles
  of the globe (relief included, front faces only, as the GPU draws them),
  the cells to test found by a walk along the ray that starts at
  `SphereGrid::nearest`.
- The cell under the cursor highlighted on the globe while hovering (its
  vertex colours inverted), updated when the cursor moves, the view changes
  or the layer changes; no highlight when the cursor is off the globe or
  outside the window.
- A click (press and release of the left button with at most 3 px of
  movement) prints the cell's text to the terminal: id, latitude and
  longitude in degrees (at the boundary), every field's value with its unit
  (class fields as their id), the layer shown with its legend.
- `--cursor <x>,<y>` with `--screenshot`: the same pick at that pixel, the
  cell highlighted in the image, a crosshair drawn around the pixel, and the
  text printed, so the evidence does not need a mouse.

Out:
- The on-screen panel and text drawing (MV-06, proposed); browser build
  (MV-03); editing worlds; the relief bound for small worlds (MV-05).

## Implementation notes

### Size and split (proposed to the owner)

The whole draft scope comes to an estimated 1,300 changed lines of code and
tests (picking 450, text drawing and the panel 600, the rest 250) plus the
documents, well above one commit. The split keeps picking, the highlight and
the terminal text here (about 750 lines) and moves the on-screen panel to
the new draft MV-06 (a ported public-domain 8×8 bitmap font drawn as
textured quads, no new crate; about 600 lines), on which MV-03 then depends
as well (the browser has no terminal).

### Picking (`pick.rs`)

The drawn surface. A vertex `(x, y, z, h)` of `GlobeMesh` is drawn at
`q = (x, y, z)·(1 + σ·h)` with `σ` the frame's `relief_scale_per_m` (the
`f32` handed to the GPU); the picker computes the same in `f64`. The fans of
all cells tile the unit sphere once (MV-02's mesh tests: shared vertices
bit-identical, every triangle counter-clockwise, solid angles summing to
4π), and the radial projection of a straight triangle is the spherical
triangle of its vertex directions, so while every drawn radius `1 + σ·h`
is positive the surface is a radial graph: every direction from the centre
meets it once.

The ray. `cursor_ndc(x_px, y_px, width_px, height_px) = [2x/W − 1,
1 − 2y/H]` (pixel coordinates from the top-left corner as winit reports
them; the centre of pixel `(i, j)` is `(i + ½, j + ½)`).
`OrbitCamera::ray(aspect, ndc) -> (Vec3, Vec3)`: origin the eye, direction
`normalise(f + s·ndc_x·a/k + u·ndc_y/k)` with the camera's forward, right
and up vectors `f, s, u`, `a = W/H` and `k = 1/tan(FOV_Y/2)`, the inverse of
`view_proj`'s perspective.

The search. `Picker::new(mesh: &GlobeMesh) -> Picker` stores per cell `v`
the cosine of `ρ_v + 10⁻⁵` rad, where `ρ_v` is the largest angle between
the cell's drawn centre direction `ĉ_v` and its drawn Voronoi vertex
directions (the fan lies in this spherical cap; the margin covers the `f32`
rounding of the vertices, about 10⁻⁷ rad, and the rounding of the arc's
plane), the smallest `cos ρ_v`, and the lowest and highest vertex heights
`h_min ≤ 0 ≤ h_max`. `Picker::pick(&mut self, grid, mesh, σ, origin,
dir) -> Option<Hit>` (`Hit { cell, distance }`, planet radii):
1. Shell radii `r_top = 1 + σ·h_max` and `r_bot = (1 + σ·h_min)·min_v
   cos ρ_v`, both widened by 10⁻⁶ relative: no point of a flat triangle
   lies above its highest vertex, but it sags below its vertices, and a
   point `x` of cell `v`'s fan satisfies `|x| ≥ x·ĉ_v ≥ minᵢ |qᵢ| cos ρ_v`
   over its vertices `qᵢ` (review R1-1). Entry `t0` into the sphere
   `r_top` (none: no hit); `t1` the entry into the sphere `r_bot`, or the
   exit from `r_top` if the ray misses it. The first hit lies in
   `[t0, t1]`: at `t1` the ray is on or under every drawn point.
2. The directions of the ray over `[t0, t1]` sweep the great-circle arc
   from `a = x(t0)/|x(t0)|` to `b = x(t1)/|x(t1)|`, shorter than π (the
   segment stays at radius `≥ r_bot > 0`).
3. Candidates: starting at `grid.nearest(a, hint)` (`hint` the last start;
   the walk from a cell nearby takes a step or two), a walk over neighbours
   (depth-first, with a stack; the order does not change the result)
   through the cells whose cap comes within `ρ_v + 10⁻⁵` of the arc (an
   endpoint within the cap, or the point of the great circle nearest `ĉ_v`
   inside the arc with `√(1 − (ĉ_v·n)²) ≥ cos(ρ_v + 10⁻⁵)`, `n` the unit
   normal of the plane of the ray and the centre; for an arc shorter than
   10⁻⁷ rad, the endpoint test alone). The cells whose drawn fans the arc crosses form a
   connected chain, each inside the set, so the walk reaches every cell
   that can be hit.
4. Each candidate's triangles are tested with the Möller–Trumbore algorithm
   (T. Möller and B. Trumbore 1997, "Fast, minimum storage ray-triangle
   intersection", *J. Graphics Tools* 2(1), 21–28), front faces only
   (`det > 0` for triangles counter-clockwise from outside, as the GPU
   culls back faces), barycentric bounds inclusive, `t > 0`. The hit is the
   smallest `t`, ties to the lower cell index.
5. If `r_bot < 0.01` (a small world at a high exaggeration, turned inside
   out until MV-05), every triangle is tested instead: slow (all 3.9
   million triangles of a 650,000-cell world) but still what is drawn.

The cell is the one whose triangle is hit, not the nearest centre, so the
pick agrees with the drawn fans up to the `f32` rounding of the GPU. Under
the translucent sea the pick is the terrain seen through it. The `Locator`
named in the draft is not used: it borrows the grid, which the window's
state owns, and the start cell of the previous pick is a better start for a
moving cursor (the first pick of a session walks about `√n` steps, some
20 µs for 650,000 cells).

Expected cost: a ray looking down tests about ten cells; a grazing ray at
the limb sweeps an arc of about 37° at ×20 on Earth (some 150 cells long,
about 500 candidates, 3,000 triangles) and about 80° at ×100. Measured in
release on the 650,000-cell world and recorded.

### The cell's text (`inspect.rs`)

`cell_text(world: &World, layers: &Layers, cell: u32) -> Vec<String>`:

```
cell 123456 of 650000, lat -12.3456 deg, lon 98.7654 deg
  cell_area             1.93302e-5 sr
  neighbor_count        6
  surface_height        -4321.57 m
  ...
layer 13/13, surface_height (m): elevation, blue below the sea level at ...
```

Latitude and longitude of the cell's centre from `coords::lat_lon`,
divided by `RADIANS_PER_DEGREE`, 4 decimals. Field names padded to the
longest; quantities with 6 significant digits, trailing zeros trimmed,
scientific notation below 10⁻³ and from 10⁶ (`format_value`); class fields
as their integer id; a value that is not finite as Rust prints it (`NaN`,
`inf`). The last line is the layer's index and `Layer::legend()`. The same
lines are printed on a click and, in MV-06, drawn in the panel.

### The highlight and the window

- `Layers::cell_colour(world, cell) -> Rgb8` (the layer's colour of one
  cell) and `highlight(colour) = 255 − colour` per channel.
  `Renderer::set_cell_colour(vertices, colour)` rewrites one cell's range
  of the colour buffer (`GlobeMesh::cell_vertices`; `write_buffer` at
  `4·first_vertex`, aligned).
  The window restores the previous cell's colours and writes the new ones
  when the picked cell changes, and writes the highlight again after a
  layer change.
- `app.rs`: the cursor position is kept; a pick runs after `CursorMoved`,
  a drag, a zoom, a relief step, a reset, a resize and a layer change, and a
  redraw is requested when the picked cell changes; `CursorLeft` clears it.
  `is_click(press, release)`: at most 3 px apart. The benchmark ignores the
  cursor.
- `args.rs`: `--cursor <x>,<y>` (pixel indices inside `--size`), only with
  `--screenshot`. `main.rs`: the pick at the pixel centre, the highlight
  written through `Renderer::set_cell_colour` before drawing (the write the
  window uses), the text printed with the time of the pick, and after
  read-back a crosshair (four white arms from 3 to 12 px from the pixel,
  outlined in black from 2 to 13 px, the pixel and its eight neighbours left
  as drawn) painted into the pixels by `crosshair(pixels, width, height,
  x, y)`.
- `README.md` and `USAGE`: hover, click and `--cursor`.

### Tests (no GPU needed)

- `camera`: for several cameras and aspects, points `origin + t·dir`
  (`t` = 0.5, 2, 5) of the ray through an NDC point project back to it with
  `view_proj` within 10⁻¹², `dir` is a unit vector and the screen centre's
  ray is the forward direction; `cursor_ndc` by hand (corners, centre).
- `pick`, octahedron by hand: camera on +X at 3 radii, heights with
  `σ·h₀ = 0.25`: the ray to the centre hits cell 0 at distance
  `3 − 1.25 = 1.75`; rays pointing away or passing beside the globe miss.
- `pick`, rays straight down: on a generated 1,000-cell world (seed 3,
  stages grid, tectonics, crust) at exaggerations 0, 1, 20 and 100, for
  every cell and the centroid direction `ĝ` of each of its triangles, the
  ray from `3ĝ` towards the centre hits that cell at direction `ĝ` (within
  10⁻¹², on the near side): a radial ray meets the radial graph once, so
  the expected cell needs no algorithm; also through `OrbitCamera` at the
  latitude and longitude of `ĝ` and the screen centre.
- `pick` against brute force: the same world and a variant with 20 km
  spikes on every seventh cell (strong occlusion), at the four
  exaggerations, cameras at several latitudes, longitudes and distances
  (down to just above the relief, so rays graze), 21 × 21 cursor positions
  each: the pick equals the smallest front-facing hit over every triangle
  of the mesh, computed in the test by an independent formula (the ray's
  intersection with the triangle's plane and a point-in-triangle test by
  signed volumes): the same distance within 10⁻¹² and the same cell, or,
  on a shared edge, a cell whose triangle the reference also hits at that
  distance; misses agree.
- `pick`, inside out (`r_bot < 0.01`, synthetic heights): agrees with
  brute force.
- Added in review: limb rays inside the sag of the flat triangles (×0 and
  ×1), oblique rays on the flat octahedron, rays grazing the corners of
  cells, back faces in front of the surface, the tie rule by hand.
- `inspect`: `format_value` by hand (0, −4321.5734, 1.9330183543·10⁻⁵,
  6.371·10⁶, 123456, 999999.5 → `1e6`, 0.001, NaN); `cell_text` of the
  1,000-cell world: the header with latitude and longitude computed in the
  test with std's `to_degrees`, one line per field in the world's order
  with its unit, class ids as integers, the legend line.
- `layers`: `cell_colour` equals `colours()[cell]`; `highlight` by hand.
- `app`: `is_click` at 0, 3 and 3.1 px.
- `args`: `--cursor` parsed, and refused without `--screenshot`, outside
  the image or malformed.
- `main`: `crosshair` by hand on a small buffer (arms, outline, untouched
  centre, clipped at the edges).

## Acceptance criteria

- [x] Picking finds the cell under the cursor for any view and relief:
      the `pick` tests above pass (radial rays find their known cells at
      exaggerations 0, 1, 20 and 100; every ray of the camera sweeps
      agrees with brute force over all triangles, also with occluding
      spikes and inside out).
- [x] A pick on the 650,000-cell world takes at most 1 ms in release for a
      view straight down and for a grazing view at the limb at ×20 and
      ×100 (`--cursor` prints the time; numbers recorded).
- [x] The clicked cell's text lists the cell id, latitude, longitude, every
      field's value with its unit (class fields as ids) and the layer's
      legend (`inspect` tests); screenshots with `--cursor` of a
      20,000-cell world and of the 650,000-cell world (close view, relief
      ×20 and ×100, one near the limb) show the highlighted cell under the
      crosshair, saved under `out/MV-04/` and described with the printed
      text.
- [x] `mise run ci` passes.

## Verification

```sh
cargo run --release -p kuniumi-viewer -- out/MV-02/650k
cargo run --release -p kuniumi-viewer -- out/MV-02/20k --distance 1.3 --cursor 640,400 --screenshot out/MV-04/20k_centre.png
cargo run --release -p kuniumi-viewer -- out/MV-02/650k --relief 100 --cursor 640,80 --screenshot out/MV-04/650k_limb_x100.png
```

## Result

What was done:
- `pick.rs` (new): `Picker` (per cell the cosine of its cap's radius, the
  extreme vertex heights, the hint and the visit marks of the walk),
  `Picker::pick` and `Hit { cell, distance }`, as planned: the shell
  `[t0, t1]`, the arc of the ray's directions, the walk over neighbours from
  `SphereGrid::nearest` through the caps that reach the arc, Möller–Trumbore
  on front faces, every triangle below `MIN_BOTTOM_RADIUS`. The module
  documentation gives the argument.
- `camera.rs`: `OrbitCamera::ray` and `cursor_ndc`. `inspect.rs` (new):
  `format_value`, `cell_text`. `layers.rs`: `Layers::cell_colour`,
  `highlight`. `mesh.rs`: `GlobeMesh::cell_vertices`. `gpu.rs`:
  `Renderer::set_cell_colour`. `app.rs`: the cursor kept, a pick after
  every cursor move, drag, zoom, key and resize, the highlight moved when
  the cell changes and painted again after a layer change, `CursorLeft`,
  `is_click`, the text printed on a click. `args.rs`: `--cursor`.
  `main.rs`: `Scene::picker`, `Scene::pick`, the screenshot with
  `--cursor` (pick at the pixel centre, its time printed, highlight,
  `crosshair` after read-back). `README.md`, `CHANGELOG.md`; the split:
  MV-06 (draft), MV-03 depends on it, the roadmap row.

Deviations from the plan:
- The cap margin is 10⁻⁵ rad instead of 10⁻⁶, and the walk is replaced by
  the test of every triangle below a drawn radius of 0.01 rather than 0:
  working out the precision of the arc's plane (`Arc`'s documentation),
  the plane computed from `o × d` tilts by up to `2ε|o|/|o × d|`, which
  stays inside 10⁻⁵ for rays that reach a radius of at least 0.01; arcs
  shorter than 10⁻⁷ rad are tested by their endpoints. 10⁻⁵ rad is 64 m on
  Earth against cells of 30 km, so the walk costs the same.
- The crosshair's white arms run from 3 to 12 px with the black outline
  from 2 to 13 px, so the pixel and its eight neighbours stay as drawn (the
  plan said "2 to 12 px").
- The lower bound of the drawn surface is `(1 + σ·h_min)·min_v cos ρ_v`,
  not `1 + σ·h_min` as planned: flat triangles sag below their vertices
  (by about 7·10⁻⁴ radii on a 1,000-cell grid, 1.2·10⁻⁶ at 650,000 cells),
  so with the planned bound a ray grazing the limb could be cut off before
  its hit and return no cell (review R1-1; the notes above now give the
  bound and its proof).
- `Hit` has no `point` (nothing uses it; the distance gives it); the walk
  is depth-first with a stack (the plan said breadth-first; the order does
  not change the result); the GPU write is
  `Renderer::set_cell_colour(vertices, colour)` (planned:
  `set_cell_colours(first_vertex, colours)`), and the screenshot writes the
  highlight through it rather than into the colours before the upload
  (review R1-2), so it runs the window's path.
- The camera sweeps use 21 × 21 screen points instead of 41 × 41, to keep
  the brute-force tests at about a second in debug; the corner-grazing
  and limb tests target the rays the sweeps rarely produce.
- Where the ray meets only the translucent sea surface, at the rim of the
  limb above a sunken ocean floor, there is no cell: the pick is the
  terrain, as planned, and the rim has none under it.
- Added: `GlobeMesh::cell_vertices`, `Offscreen::renderer`; the tests
  `back_faces_in_front_are_not_picked` and
  `rays_grazing_the_corners_of_cells_agree_with_brute_force`, written after
  mutations of the picker survived the planned tests (below), and in
  review `limb_rays_in_the_sag_of_flat_triangles_agree_with_brute_force`,
  `oblique_rays_on_the_flat_octahedron_agree_with_brute_force` and
  `ties_go_to_the_lower_cell`.
- Size: about 830 lines of code and documentation and 750 of tests in the
  viewer, against the estimate of 750 in all; the excess is the documented
  argument of the search and the brute-force tests.

Evidence:

1. Picking (`pick` tests, 10, plus `camera` and `main` tests):
   - `octahedron_by_hand`: cell 0 at distance 1.75 (centre drawn at 1.25,
     σ = 2⁻¹²), cell 5 at 2 from below, misses away from and past the
     globe, three rays against brute force.
   - `back_faces_in_front_are_not_picked`: an octahedron cell drawn through
     the centre to x ≈ −1.4, its faces turned inwards: seen from
     (−3, 0, 0), cell 1 at distance 2, as the GPU culls them.
   - `radial_rays_find_their_cell`: a generated 1,000-cell world (seed 3,
     stages grid, tectonics, crust) at ×0, ×1, ×20 and ×100: for every cell
     and every one of its 5,988 triangles in all, the ray from 3ĝ towards
     the centre (ĝ the direction of the triangle's centroid) picks that
     cell, on the near side in the direction ĝ (within 10⁻¹²); every tenth
     cell also through `OrbitCamera` and the screen centre.
   - `camera_rays_agree_with_brute_force` and
     `occluding_spikes_agree_with_brute_force` (the same world, and with
     20 km more on every seventh cell, spikes of 0.31 radii at ×100): 4
     exaggerations × 4 cameras (latitude 20°, 3 radii; −60°, 1.5 radii; 5°,
     0.02 above the top of the relief; 89°, 2 radii) × 21 × 21 screen points
     = 7,056 rays each, every pick equal to the nearest front-facing hit
     over all triangles computed by the test's own formulas (plane
     distance, signed volumes of the edges): the same distance within
     10⁻¹² relative on a triangle of the picked cell (an edge counts within
     10⁻⁹ of the edges' scale), misses only where no triangle is hit;
     between a third and all of the rays hit.
   - `limb_rays_in_the_sag_of_flat_triangles_agree_with_brute_force`
     (review R1-1): the 1,000-cell world at ×0 and ×1, rays from 3 radii
     passing 1 to 0.997 radii from the centre (16 steps), from 40 origins
     in 8 azimuths each, 5,120 rays per exaggeration, inside the sag of the
     flat triangles at the limb. `oblique_rays_on_the_flat_octahedron_…`:
     the octahedron at σ = 0, rays passing 0 to 1 radius from the centre,
     4,400 rays. Both failed before the fix (the reviewer found 28 of 4,000
     limb rays and 192 of 960 octahedron rays without a pick) and pass
     after it. `ties_go_to_the_lower_cell` (review R1-4): `closer` by hand.
   - `rays_grazing_the_corners_of_cells_agree_with_brute_force`: a
     500-cell Fibonacci grid with heights from −5.0 to +6.7 km at the four
     exaggerations: for every triangle, rays at 60° from the vertical
     towards a point 10⁻⁴ (barycentric) from its outer corner, moving
     across the corner in both senses (2,988 triangles, 5,976 rays per
     exaggeration).
   - `inside_out_relief_tests_every_triangle`: every fifth cell 100 km deep
     at ×100 (drawn radius −0.57): 242 rays against brute force.
   - `camera::rays_project_back_to_their_screen_point` (3 cameras, 3
     aspects, 4 screen points, points at 0.5, 2 and 5 along the ray project
     back within 10⁻¹²), `camera::cursor_positions_by_hand`,
     `tests::the_scene_picks_what_the_view_shows` (`Scene::pick` with the
     σ of the view: the screen centre above a triangle of cell 42 picks it;
     the corner of the view picks nothing).
   - Mutations of `pick.rs`, each run against the whole viewer suite in
     release: without the arc's interior test, with the arc reduced to its
     start, with only the first triangle hit of a cell, with each cap from
     one vertex only, without the every-triangle fallback, with caps 10⁻²
     rad too small: caught by the planned tests. Caps 3·10⁻³ and 10⁻³ rad
     too small and a two-sided triangle test survived them; the two added
     tests catch them (corner-grazing rays and back faces in front, the
     latter also by `inside_out_relief_tests_every_triangle`). Using the
     exit of the top shell instead of the entry into the bottom one only
     lengthens the arc (more cells tested, the same picks) and is not
     caught, as expected.
   - Test-first: commit "tests against stubs": 17 tests failed (the 15 new
     ones and 2 existing ones that build a `Scene`, whose picker was a
     stub). In review, the two sag tests were run before the fix and
     failed.
2. Time of a pick on the 650,000-cell world (`out/MV-02/650k`, release, 13th
   Gen Intel Core i7-13700KF; the first pick of a process, from hint 0,
   printed by `--cursor`): straight down (`--cursor 640,400`, ×20)
   0.066 ms; at the top limb at ×20 (`--cursor 640,y` for y = 65 to 110,
   grazing over the north pole) 0.050 to 0.104 ms; at ×100 (y = 75 to 110,
   over the sunken polar sea onto the raised far side) 0.120 to 0.252 ms;
   a close view at ×100 0.091 ms. Misses over the limb take longer (the
   arc spans the whole chord through the shell): at ×100, y = 30 to 72,
   0.079 to 0.243 ms here and 0.446 ms in the reviewer's run; at ×20,
   y = 58 to 64, 0.072 to 0.098 ms. All below 1 ms; unchanged after the
   round-1 fixes (0.061 ms down, 0.109 and 0.205 ms at the limb). Building the picker adds
   about 20 ms to the load (mesh, layers and picker 0.14 to 0.15 s; mesh
   and layers 0.12 to 0.13 s in MV-02).
3. The cell's text (`inspect` tests: `values_by_hand`, 19 values by hand;
   `the_text_of_a_cell`, cells 0, 1, 517 and 999 of the 1,000-cell world
   against latitude and longitude from std's `to_degrees`, every field in
   order with its unit, class ids as integers, the legend line, and a
   world without fields). `layers::cell_colours_and_the_highlight`,
   `app::clicks_and_drags`, `args` (`--cursor` accepted; refused without
   `--screenshot`, with `--benchmark`, outside the image, malformed),
   `tests::crosshair_by_hand`. Screenshots (`--screenshot --cursor`,
   1280×800, viewed; text excerpts as printed):
   - `out/MV-04/20k_centre.png` (20,000 cells, 1.3 radii, ×20): coast in
     the upper right, shelf and deep sea to the left; the cell under the
     crosshair (640, 400) is a purple hexagon (the inverse of the
     light-blue shelf, tinted by the sea), the crosshair inside it near its
     left edge. Printed: cell 6532, lat 19.8269, lon 0.5382 (under the
     camera at lat 20, lon 0), surface_height −1689.49 m (917 m under the
     sea level of −772 m: shelf), 13 field lines, the legend.
   - `out/MV-04/20k_plates_x100.png` (`plate_id`, 1.25 radii, ×100,
     cursor 500, 250): the light-orange plate (seen through the sea as
     mauve) beside a blue one; the picked hexagon dark blue (the inverse of
     plate 6's light orange) with the crosshair inside; printed plate_id 6,
     neighbor_count 6, surface_height −3809.78 m.
   - `out/MV-04/20k_wall_x100.png` and its 4× crop `…_crop.png` (×100, no
     sea, the foot of a shelf wall, cursor 470, 545): the picked cell
     orange (inverse of deep blue), partly hidden behind the wall's facets
     in front of it; the crosshair's open centre shows orange, so the
     drawn cell at that pixel is the picked one; printed surface_height
     −5696.17 m (the deep floor, not the shelf in front).
   - `out/MV-04/650k_close_x100.png` (650,000 cells, 1.1 radii, lat 5, lon
     −25, ×100, cursor 700, 300): a flat green continent beside the shelf;
     the picked cell a purple pentagon with the crosshair inside; printed
     neighbor_count 5, surface_height 300 m.
   - `out/MV-04/650k_down.png` (default view, cursor 640, 400): cell
     214129, lat 19.9256, lon 0.0805, picked 2.0055 radii from the camera
     (surface_height −2449.47 m, 1,751 m under the sea level, drawn at
     ×20 at radius 1 − 0.0055).
   - `out/MV-04/650k_limb_x20.png` (cursor 640, 66) and its 10× crop: the
     crosshair at the top limb, on the inner edge of the dark rim of the
     sea, where the ocean floor seen through the sea begins; cell
     898, lat 85.7150, lon 1.1028 (over the pole), −5998.45 m. At y = 62 to
     64 the pixel shows only the rim: "no cell".
   - `out/MV-04/650k_limb_x100.png` (cursor 640, 80) and its 10× crop: at
     ×100 the ocean floor sinks, the rim widens, and beyond the north pole
     a raised continent shows above the horizon; the crosshair sits on its
     silhouette: cell 27975, lat 65.9841, lon 179.7340 (the far side),
     crust_type 1, surface_height 629.765 m, 3.232 radii from the camera
     (more than the 3 to the centre).
   - Agreement of the pick with the GPU's image (throwaway script, PNGs
     decoded with Python's standard library): for random pixels, a
     render with `--cursor` against one without: every picked cell's
     highlight changed the pixel, and every pixel without a cell was the
     background (48, 48, 48): 20,000 cells, 1.6 radii, ×100, no sea: 30
     hits of 30; 650,000 cells, ×20: 11 hits, 19 background; 650,000
     cells, 2 radii, ×100: 31 hits, 9 background. After the round-1 fixes
     (the highlight now written by `Renderer::set_cell_colour`): 20,000
     cells at ×0, no sea: 25 of 25; 650,000 cells, 2 radii, ×100: 20 hits,
     5 background. All seven screenshots above rendered again: byte for
     byte identical.
   - The window: `--benchmark 120` opens it and presents at 60.1 fps (mean
     16.65 ms) with the new code; hover and click could not be driven from
     this session (no input tool).
4. `mise run ci` passes (log `out/MV-04/ci.log`, after the round-1
   fixes): fmt, clippy `-D warnings`, doc, 1,169 tests passed in debug and
   release together (10 ignored: the existing timing tests; the viewer's
   44 tests take 2.7 s in debug and 0.6 s in release), MSRV 1.85.0 without
   the viewer, cargo-deny (advisories, bans, licences, sources ok;
   duplicate-version warnings only).

Review:
- Round 1 (CHANGES_REQUESTED): the reviewer re-ran `mise run ci`, the
  timings and screenshots, and probed the picker with scratch tests
  against the brute-force reference. R1-1 blocker: the planned lower
  bound `1 + σ·h_min` of the drawn surface ignores the sag of flat
  triangles below their vertices, so rays grazing the limb at ×0 (or in a
  world without heights) could return no cell for a visible face (28 of
  4,000 limb rays on the 1,000-cell world, 192 of 960 oblique rays on the
  octahedron) (fixed: the bound `(1 + σ·h_min)·min_v cos ρ_v` with its
  proof, two tests that failed before the fix). R1-2 minor: the hover
  highlight's GPU write was never run (fixed: the `--cursor` screenshot
  writes the highlight through it; images byte-identical). R1-3 minor:
  the ticket's notes disagreed with the code in several names and
  numbers, and the verification command was not a limb view (fixed: notes
  updated, deviations listed, the planned 10⁻¹² tolerance restored). R1-4
  nit: the tie rule untested (fixed: a test by hand). R1-5 nit: the
  changelog's pick time omitted misses, up to 0.446 ms (fixed: "under
  0.5 ms").
- Round 2 (APPROVE): the reviewer re-ran `mise run ci`, re-rendered three
  screenshots (byte-identical), confirmed R1-1 to R1-5 resolved with
  mutations (the old bound and a bound made too high fail the new tests;
  the reversed tie rule fails), probed 4,000 limb rays per case on the
  1,000-cell world at all eight relief steps, flat grids from 4 to 20,000
  cells and the octahedron, about 2 million oblique rays and a deep-pit
  world just above the fallback radius, all without a difference from an
  exhaustive search, and measured picks of 0.051 to 0.265 ms on the
  650,000-cell world. R2-1 nit: the Result counted 11 `pick` tests
  instead of 10 (fixed). R2-2 nit: a test comment gave the octahedron's
  sag as 1/√3 (the picker's bound); the drawn fans sag to about 0.81
  radii (fixed).

Follow-ups: MV-06 (draft, the on-screen panel, split from this ticket when
it was refined; MV-03 now depends on it). The roadmap's MV criterion is not
met yet: the panel (MV-06), the browser (MV-03) and the relief bound
(MV-05) remain.

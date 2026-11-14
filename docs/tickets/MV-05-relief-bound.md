# MV-05 Relief exaggeration bounded by the deepest point

- Status: done
- Milestone: MV
- Depends on: MV-02
- Crates: kuniumi-viewer
- Design: 実装ロードマップ / 技術スタック
- Decisions: D-12, D-19

## Goal

No relief exaggeration the viewer offers turns a valid world inside out:
the drawn radius `1 + e·h/R` of the deepest point stays clearly above 0 on
every planet within the bounds of D-12.

## Scope

In (review R3-1 of MV-02):
- The exaggeration of `kuniumi-viewer` is capped by its world, not only by
  the fixed ×100 of `RELIEF_STEPS`: the deepest height above the sea level
  is stored beside the highest (`view::Relief`), and the Up key, `--relief`
  and the start view keep `1 + e·h_min/R` at or above a floor of 0.5 planet
  radii (chosen when the ticket was refined: a cap, not a floor in the
  vertex shader; see below).
- The doc comment of `MAX_RELIEF` (`args.rs`) gives Earth's trenches as
  7 km deep; Earth's deepest point is about 11 km below the sea level
  (7.2 km is the lowest cell of MV-02's generated world).

Out:
- Everything else of the viewer (MV-03, MV-07); a bound for the
  highest point (a high exaggeration only moves the camera out).

## Implementation notes

Found by the third review of MV-02: a world of the smallest radius
(1,500 km, Earth's density, Earth's ocean volume, seed 1, 100,000 cells)
has an ocean 46 km deep on average and its lowest cell 50.7 km below the
sea level, so the drawn radius falls below 0 for `e > 29.6`; at ×100 the
inverted sea floor fills the view around the sea sphere, at ×20 it shrinks
to a ball of 0.3–0.5 radii under it. The radius of the deepest point at
×100 measured on 20,000-cell worlds: 0.45 at 3,000 km, 0.78 at 4,500 km,
0.89 at Earth's size, 0.97 at 13,000 km. No crash; Down recovers the view.
The near and far planes come from the highest point only, so they need the
deepest one too if the floor is in the shader.

Reproduced when the ticket was refined: `out/MV-05/small.toml`
(`radius_m = 1500000.0`, `surface_gravity_m_s2 = 2.3089`, i.e.
9.80665 m/s² × 1,500/6,371, Earth's bulk density; everything else Earth's),
`kuniumi generate --seed 1 --cells 100000 --planet out/MV-05/small.toml`:
sea level 42,681.6 m above the reference, every cell under the sea, from
37.1 to 50.7 km deep (`metrics.txt`: `lowest_m = -50681.6`). The viewer of
`main` at `--relief 100` draws the inverted sea floor around the sea sphere
over the whole view (`out/MV-05/before_small_x100.png`). Earth (seed 1,
100,000 cells): lowest cell 7,283 m below the sea level.

### Size

About 250 changed lines in `kuniumi-viewer` (half of them tests) plus the
documents; one part, no split.

### The choice: a cap per world, not a floor in the shader

- A floor in the vertex shader, radius `max(1 + e·h/R, r_min)`, keeps ×100
  on every world but flattens every point deeper than `(1 − r_min)·R/e` onto
  one sphere: on the 1,500 km world at ×100 that is everything deeper than
  7.5 km, i.e. its whole sea floor (37–51 km deep), drawn as a smooth ball
  of radius 0.5. The picker (MV-04), which finds what is drawn, would need
  the same floor, and the near and far planes the deepest point.
- A cap per world keeps the exaggeration uniform, so the globe is the world
  scaled, and nothing downstream changes: the shader, the picker (its
  neighbour walk then always applies, since its lower bound `r_bot` is at
  least 0.5 · min cos ρ, far above `MIN_BOTTOM_RADIUS` = 0.01), the near
  and far planes (every drawn vertex lies between radius 0.5 and `top`)
  and the sea pass.

### The cap

`view::MIN_DRAWN_RADIUS = 0.5` planet radii: no vertex is drawn closer to
the centre than half the radius. With the depth `d = −h_min > 0` (m) of the
deepest vertex below the sea level and the planet radius `R` (m), the
largest exaggeration of a world is

    e_max = min(MAX_RELIEF, (1 − MIN_DRAWN_RADIUS)·R/d·(1 − δ)),  δ = 10⁻⁶,

and `MAX_RELIEF` (100) for a world without heights or without a point below
the sea level. By hand: the 1,500 km world, `d` = 50,681.6 m: 0.5 ×
1,500,000/50,681.6 = 14.80 (instead of 100); Earth's generated worlds,
`d` ≈ 7.3 km: 437, so 100; Earth's Challenger Deep, `d` = 10,935 m: 291, so
100: Earth keeps ×100. A 3,000 km world with the 16.5 km deep floor of the
MV-02 review: 90.9. The cap holds for any depth, even one beyond the radius
(then `e_max` < 1).

The margin `δ`: the GPU computes the radius in `f32` from
`σ = fl32(e/R)` and `h = fl32(h_v)`; the roundings of `σ`, of `h`, of the
product and of the sum (at most 2⁻²⁴ relative each, or one rounding with a
fused multiply-add) move the radius near 0.5 by less than 1.2·10⁻⁷, and
`δ` raises the exact radius at `e_max` to 0.5 + 5·10⁻⁷, so the drawn radius
is at least 0.5.

### API (`crates/kuniumi-viewer/src/`)

- `mesh.rs`: `GlobeMesh::min_height_m: f64`, the lowest height above the
  sea level of any vertex (m), at most 0, computed as `max_height_m` is
  (the Voronoi vertices' means of three lie within their cells' heights).
- `view.rs`:
  - `pub const MAX_RELIEF: f64 = RELIEF_STEPS[RELIEF_STEPS.len() − 1]`,
    moved here from `args.rs`, its comment corrected: at ×1000 Earth's
    deepest point, about 11 km below the sea level on a radius of
    6,371 km, would be drawn beyond the centre (1 − 1000 × 10.9/6,371 =
    −0.72), inside out.
  - `pub const MIN_DRAWN_RADIUS: f64 = 0.5`; the margin `δ` as a constant.
  - `Relief::min_height_m: f64` (at most 0), from the mesh in `Scene::new`.
  - `Relief::max_relief(&self) -> f64`: `e_max` above.
  - `relief_up(relief: f64, max: f64) -> f64`: the smallest step above
    `relief`, lowered to `max`; `relief` itself if no step lies above it or
    `relief ≥ max`. `relief_down` is unchanged (from 14.80 down to 10).
  - `relief_text(relief: f64) -> String`: the exaggeration rounded to two
    decimals for the title and the terminal (`20`, `14.8`).
- `main.rs`: `start_view(args: &Args, relief: &Relief) -> View`: the camera
  of the arguments, the relief `min(args.relief, relief.max_relief())`, the
  distance clamped as now. `run` prints a line when it lowered the requested
  relief (`relief x100 lowered to x14.8, the most for this world: its
  deepest point, 50.7 km below the sea level, stays at half the radius`).
- `app.rs`: `press` passes `scene.relief.max_relief()` to `relief_up`; the
  terminal line and the window title use `relief_text`; at the cap the
  terminal line says so.
- `args.rs`: uses `view::MAX_RELIEF`; the usage of `--relief` adds that
  each world may lower it.
- `pick.rs`: step 5 of the module documentation no longer names "a small
  world drawn inside out until MV-05"; the search of every triangle stays
  for meshes whose cells reach a hemisphere (coarse grids).
- `README.md`: one clause on the cap.

### Tests (written first)

- `view::tests`: `max_relief` by hand: `R` = 1,500 km with `h_min` =
  −50 km gives 15·(1 − δ); Earth's radius with `h_min` = −10,935 m gives
  100; `h_min` = 0 or no heights give 100; `relief_up` with a cap: (10,
  14.8) → 14.8, (14.8, 14.8) → 14.8, (5, 14.8) → 10, (0, 0.5) → 0.5, and
  with the cap at 100 as before; `relief_text`.
- `mesh::tests`: `min_height_m` on the octahedron by hand, and no vertex
  lower.
- `app::tests`, one test per world: the 1,500 km world (seed 1, 100,000
  cells, `g` = 2.3089 m/s²; about 3 s to generate in debug) and Earth
  (seed 1, 100,000 cells). For every start relief `--relief` gives (0, 0.5,
  every step, 14, 14.79, 15, 37, 99.5, the default and 100), parsed by
  `args::parse` and passed through `start_view`, and every relief that up to
  ten presses of Up and of Down reach from it through `press`: the drawn
  radius of every vertex, in `f32` as the shader computes it (`1 + σ·h`
  with separate roundings, and with `mul_add`) and in `f64` from the `f32`
  values as the picker does, is at least 0.5. The small world's largest
  reached relief is its `max_relief`, below 100 and within 10⁻⁴ of
  0.5·R/d with `d` its lowest `surface_height` below the sea level computed
  in the test (14.80); Earth's is 100.
- The existing `keys` and `frames` tests follow the new signatures.

## Acceptance criteria

- [x] `Relief::max_relief` matches the hand values: a 1,500 km planet whose
      deepest point lies 50 km below the sea level gets 15·(1 − 10⁻⁶);
      Earth with the Challenger Deep (10,935 m) gets 100; a world without
      heights or without a point below the sea level gets 100 (test).
- [x] On the 1,500 km world above and on Earth (seed 1, 100,000 cells each),
      the drawn radius of every vertex stays at or above 0.5 for every
      exaggeration the keys and `--relief` reach, computed in `f32` as the
      shader does and in `f64` as the picker does (test); the small world's
      largest reachable exaggeration is its cap, 14.80 (within 10⁻⁴ of the
      hand value), Earth's is 100.
- [x] Screenshots of the small world at `--relief 100`, with and without the
      sea, saved under `out/MV-05/` and described, beside the inverted view
      of `main`'s viewer (`out/MV-05/before_small_x100.png`); Earth at
      `--relief 100` renders byte-identical to `main`'s viewer.
- [x] The doc comment of `MAX_RELIEF` gives Earth's deepest point as about
      11 km below the sea level.
- [x] `mise run ci` passes.

## Verification

```sh
printf '[planet]\nradius_m = 1500000.0\nsurface_gravity_m_s2 = 2.3089\n' > out/MV-05/small.toml
cargo run --release -p kuniumi-cli -- generate --seed 1 --cells 100000 --planet out/MV-05/small.toml --out out/MV-05/small
cargo run --release -p kuniumi-viewer -- out/MV-05/small --relief 100 --screenshot out/MV-05/small.png
cargo run --release -p kuniumi-viewer -- out/MV-05/small --relief 100 --no-sea --screenshot out/MV-05/small_nosea.png
cargo run --release -p kuniumi-cli -- generate --seed 1 --cells 100000 --out out/MV-05/earth
cargo run --release -p kuniumi-viewer -- out/MV-05/earth --relief 100 --screenshot out/MV-05/earth_x100.png
```

## Result

The plan was approved by the owner as refined (a cap per world, floor 0.5).

What was done:
- `view.rs`: `MAX_RELIEF` (moved from `args.rs`, its comment corrected to
  Earth's deepest point, about 11 km), `MIN_DRAWN_RADIUS = 0.5`,
  `RELIEF_MARGIN = 10⁻⁶`, `Relief::min_height_m`, `Relief::max_relief`
  (the cap, its `f32` argument in the doc comment), `relief_up(relief,
  max)`, `relief_text`.
- `mesh.rs`: `GlobeMesh::min_height_m`. `main.rs`: `start_view` (the start
  relief lowered to the cap), `relief_note` (the line printed when it is
  lowered), `Scene::new` fills `min_height_m`. `app.rs`: Up stops at the
  cap and says so on the terminal (`view::relief_line`, after review
  R1-1); the title and the terminal show the relief with `relief_text`.
  `args.rs`: the usage of `--relief` and of the keys. `pick.rs`: step 5 of
  the module documentation. `README.md`, `CHANGELOG.md` (Fixed).

Deviations from the plan:
- `relief_note` is a function of its own (tested) rather than a `println!`
  inside `run`, and it stays silent when the lowering is below the two
  decimals the title shows (`--relief 14.799` on the small world, whose cap
  is 14.798).

Evidence:

1. `Relief::max_relief` by hand (`view::tests::relief_limits_by_hand`):
   1,500 km and 50 km deep: 15·(1 − 10⁻⁶) within 10⁻¹², drawn radius
   0.5 + 5·10⁻⁷; Earth with the Challenger Deep (10,935 m): 100 exactly;
   32 km deep on Earth's radius: 99.546875·(1 − 10⁻⁶); 3,000 km deep on
   1,500 km: 0.25·(1 − 10⁻⁶); no point below the sea level, or no heights:
   100. `relief_steps`: `relief_up` with a cap (10 → 14.8, 14.8 → 14.8,
   5 → 10, 0 → 0.5 under a cap of 0.5, Up never lowers 20 under 14.8),
   the old steps under 100 unchanged; `relief_texts`;
   `mesh::tests` (`min_height_m` on the octahedron by hand, −600 m; −1 m
   for a world under the sea, 0 for land; −5,000.1 m on the Fibonacci grids
   and no vertex lower); `tests::scene_heights_lie_above_the_sea_level`
   (`Relief::min_height_m` is the lowest `surface_height − sea level`, 0
   without heights); `tests::relief_note_names_the_deepest_point`;
   `view::tests::relief_lines` (review R1-1: the line of Up and Down names
   the cap at the small world's cap, not below it nor at ×100).
2. The two worlds (`app::tests`), each generated in the test (seed 1,
   100,000 cells; all 63 tests of the viewer take 3.5 s in debug):
   - `relief_keeps_the_deepest_point_of_a_small_world_above_the_floor`: the
     1,500 km planet (`g` = 2.3089 m/s², Earth's ocean volume), deepest
     cell 50,681.6 m below the sea level, as in review R3-1 of MV-02. Its
     cap is 14.798266, within 10⁻⁴ of 0.5 × 1,500,000/50,681.6 computed
     in the test from the world's fields, and below 100. The reliefs
     reached from the default and from `--relief` 0, 0.5, 14, 14.79, 15,
     37, 99.5 and every step, through `start_view`, then ten presses of Up,
     ten of Down and ten of Up through `press`, are 0, 0.5, 1, 2, 5, 10,
     14, 14.79 and the cap, the largest. At each, every one of the 699,988
     vertices is drawn at a radius of at least 0.5: in `f32` with separate
     roundings, in `f32` with `mul_add` and in `f64` from the `f32` values
     with the `σ` of the frame. At the cap the lowest lies below 0.50001
     (0.5 + 5·10⁻⁷ by hand); at ×15 the deepest cell would lie below 0.5.
   - `relief_on_earth_reaches_the_largest_step_above_the_floor`: cap 100;
     the reliefs reached are 0, 0.5, 1, 2, 5, 10, 14, 14.79, 15, 20, 37,
     50, 99.5 and 100, the largest; every vertex at least 0.5 at each; at
     ×100 the lowest drawn radius lies in [0.85, 0.92] (1 − 100 ×
     7,283.2/6,371,000 = 0.886 by hand).
   - Test-first: against stubs (commit "tests of the relief cap against
     stubs"), 7 of the 61 tests then present failed; the Earth test passed,
     as it should, since Earth needs no cap. Mutations: `start_view`
     without the cap, or Up with the cap at 100, fail the small world's
     test; `RELIEF_MARGIN = 0` fails it too, with a fused `f32` radius of
     0.49999997 at vertex 117,032 (50,681.563 m below the sea level), and
     fails `relief_limits_by_hand`, so the margin is needed.
3. Screenshots (release build, 1280×800, camera at latitude 20°, 3 radii;
   `out/MV-05/small.toml` and the worlds as in Verification), viewed:
   - `out/MV-05/before_small_x100.png` (`main`'s viewer, `--relief 100`):
     the inverted sea floor fills the whole view with dark blue terrain
     around a grey ring and the translucent sea sphere, with splinters of
     the floor crossing the sphere at the lower right;
     `before_small_x100_nosea.png`: the same without the sea, a grey disc
     cut out of the inverted floor.
   - `out/MV-05/small.png` (this branch, `--relief 100`): the terminal
     says `relief x100 lowered to x14.8, the most for this world: its
     deepest point, 50.7 km below the sea level on a radius of 1500 km,
     stays at 0.5 of the radius`. The translucent sea sphere on the grey
     background, and inside it a darker ball about 0.54 of its diameter
     (370 against 680 px; the floor of 37–51 km lies at radii 0.50–0.63 at
     ×14.8, most of it near 0.54), with raised plateaus and scarps showing
     through. Nothing is inverted and nothing reaches outside the sea.
   - `out/MV-05/small_nosea.png`: the same ball without the sea, in dark
     blue lit from the upper left: two lines of steep scarps across the
     upper half, two raised rounded plateaus at the lower left and right,
     a wall of terraces along the right limb; a closed, outward-facing
     surface.
   - `out/MV-05/small_default.png` (no `--relief`, the default ×20 lowered
     to ×14.8): byte-identical to `small.png`.
   - `out/MV-05/earth_x100.png` (Earth, `--relief 100`): byte-identical to
     `out/MV-05/before_earth_x100.png` from `main`'s viewer; continents in
     green, shelves in light blue with tall walls, the floor under the
     translucent sea, as in MV-02.
   - In passing: `--cursor 640,400` on the small world at `--relief 100`
     picks cell 33,215 at 2.4918 radii from the camera (the floor at
     radius 0.508) in 0.032 ms with the neighbour walk; `main`'s viewer
     picked the inverted surface behind the centre, 5.23 radii away, by
     testing every triangle in 7.9 ms. A window benchmark of 60 frames
     at `--relief 50` printed the lowering note and ran.
4. The doc comment of `MAX_RELIEF` (now in `view.rs`) gives Earth's
   deepest point as about 11 km below the sea level (the Challenger Deep,
   10,935 m), drawn at 1 − 1000 × 10.9/6,371 = −0.72 at ×1000.
5. `mise run ci` passes (log `out/MV-05/ci.log`, exit 0): fmt, clippy
   `-D warnings`, doc, 1,205 tests passed in debug and release together
   (10 ignored: the existing timing tests), MSRV 1.85.0 without the viewer,
   cargo-deny (advisories, bans, licences, sources ok). After the round-1
   fixes (`out/MV-05/ci2.log`, exit 0): 1,207 passed, 10 ignored; the
   screenshots `small.png` and `earth_x100.png` rendered again are
   byte-identical.

Review:
- Round 1 (APPROVE): the reviewer re-ran `mise run ci` (1,205 passed, 10
  ignored), rebuilt the stub commit (7 of 61 failing), regenerated both
  worlds (same metrics), built `main`'s viewer itself (its Earth ×100 is
  md5-identical to the branch's, and its render of the small world at
  `--relief 14.798266` is byte-identical to the branch's at `--relief
  100`), looked at the screenshots, re-checked the f32 bound (1.19·10⁻⁷ <
  5·10⁻⁷) and ran mutations (the margin at 0, 10⁻⁷ and negative, the
  `present` check, `relief_up` ignoring `max`, `min_height_m` scaled or
  zeroed: all caught). R1-1 minor, no test covered the terminal line at
  the cap (fixed: `view::relief_line`, test `relief_lines`; the mutation
  that never names the cap now fails it). R1-2 nit, `(1 − MIN_DRAWN_RADIUS)`
  and `MIN_DRAWN_RADIUS` are both 0.5, so a test cannot tell them apart
  (not changed: the reviewer suggests nothing while the floor is a
  constant). R1-3 nit, two edited lines of `README.md` and `CHANGELOG.md`
  were not rewrapped (fixed).
- Round 2 (APPROVE): the reviewer re-ran `mise run ci` (1,207 passed, 10
  ignored), regenerated both worlds (every file byte-identical to
  `out/MV-05/`), rebuilt `main`'s viewer and the branch's, re-rendered and
  looked at the screenshots (the same md5 sums), reproduced the picks, and
  tried two extreme worlds: 1,500 km with 10¹¹ km³ of water (deepest point
  1,516 km down) is lowered to ×0.49 and drawn as a closed ball at about
  0.5, and a world without ocean keeps ×100 without a note. R1-1 and R1-3
  resolved, R1-2 accepted. R2-1 nit, no test checks that `press` hands the
  world's own cap to `relief_line` (not changed: one line, the function
  itself is tested, and a change now would land unreviewed code). R2-2
  nit, a Result line edited for R1-1 was not rewrapped (fixed when the
  Result was finalised).

Follow-ups: none.

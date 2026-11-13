# M2-17 Sea level of deep water worlds

- Status: done
- Milestone: M2
- Depends on: M2-08
- Crates: kuniumi-sim (stages::crust); a doc comment of kuniumi-world
- Design: 地殻から標高へ / 海面：海水量から決める
- Decisions: D-12

## Goal

No valid planet gets a sea level far from physical because its ocean is
deep compared with its radius.

## Scope

In:
- Of the two ways the draft offered, the owner chose B at the refinement
  ("B: 球殻として計算する (Recommended)", approving the plan below):
  - ~~A. a check across parameters in `PlanetParams::validate`: the
    ocean's mean depth `V / (4πR²)`, or its volume against the planet's
    `(4/3)πR³`, within a bound~~ (not chosen; outline below);
  - **B.** a sea level from a spherical shell for the water above the
    highest land in `stages::crust::sea_level_m`: depth
    `(R_top³ + 3(V − C)/A)^(1/3) − R_top` above the brim (`A = Σ a_i`,
    about 4π), instead of the thin layer `(V − C) / (A R²)`.
- `stages::crust::water_volume_m3` uses the same shell above the brim, so
  that it stays the inverse of `sea_level_m` for every volume.
- The design's 海面 section states the shell above the brim; the doc
  comment of `PlanetParams::MAX_OCEAN_VOLUME_M3` gives the new sea level.

Out:
- Ocean volumes and radii themselves (bounded by M2-08).
- The thin-layer bathtub below the brim, `Σ max(0, z − h_i) a_i R²`
  (design, 海面): its error against spherical shells is about
  `2|h|/R` of the volume (at most about 1 % for heights within ±8 km at
  the smallest radius), and replacing it would change every world.
- Later stages' use of the sea level and of the radius (a sea surface far
  above the datum puts the climate on a larger sphere than `radius_m`):
  not in M2.

## Implementation notes

Found in the review of M2-08 (R1-2): the bounds allow 1e20 m³ of water
(`MAX_OCEAN_VOLUME_M3`) on a 1,500 km planet (`MIN_RADIUS_M`), seven times
the planet's volume of 1.41e19 m³. The Crust stage's water-world formula
then gives a sea level of +3,533,312 m (CLI, 20,000 cells, seed 1, gravity
2.5 m/s²), while a spherical shell holding the same volume is about
1,509 km deep. Later stages (climate) read the sea level.

### Numbers measured at the refinement

`kuniumi generate --seed 1`, 20,000 cells (`out/M2-17/*.toml`), the heights
and cell areas read back from the `.kfield` files and the levels computed
in Python with `math.fsum` (thin: the stage's formula; shell: option B;
spherical: every cell's column as a spherical shell,
`Σ a_i ((R + z)³ − (R + h_i)³)/3 = V`, by bisection):

| Planet | `h_max` | Brim `C` | Thin | Shell (B) | Spherical |
| --- | --- | --- | --- | --- | --- |
| 1,500 km, 1e20 m³ | 4,221.8 m | 2.173e17 m³ | +3,533,312 m | +1,508,308 m | +1,508,308 m (0.22 m lower) |
| 1,500 km, Earth's 1.335e18 m³ | 4,221.8 m | 2.173e17 m³ | +43,751 m | +42,545 m | +42,544 m |
| 13,000 km, 1e20 m³ | 5,603.0 m | 1.824e19 m³ | +44,103 m | +43,957 m | +43,955 m |

The shell above the brim agrees with the fully spherical model within
1.5e-7 relative at the smallest radius; the thin layer is 134 % too high
there with the most water, and 2.8 % (1.2 km) too high with Earth's water.

### Option B (chosen)

With `A = Σ a_i` (sr, in cell order), `R = radius_m`, `h_max` the highest
height, `R_top = R + h_max` and `C = W(h_max)` the water the basins hold
when filled to the brim (the thin bathtub, unchanged):

- `water_volume_m3(z)`:
  - `z ≤ h_max`: `Σ_{h_i < z} (z − h_i) a_i R²`, unchanged (bit for bit);
  - `z > h_max`: `C + (A/3)((R_top + d)³ − R_top³)` with `d = z − h_max`,
    computed without cancellation as `C + A d (R_top² + R_top d + d²/3)`.
- `sea_level_m(V)`:
  - `V ≤ 0` and `0 < V < C`: unchanged (bit for bit);
  - `V ≥ C`: `z = h_max + d` with `q = (V − C) / (A/3)` (m³/sr; dividing
    by `A/3` rather than multiplying by 3 keeps `V = f64::MAX` finite),
    `b = (R_top³ + q)^(1/3)` and `d = q / (R_top² + R_top b + b²)`, the
    stable form of `b − R_top` (from `b³ − R_top³ = (b − R_top)(b² +
    R_top b + R_top²)`): `d ≥ 0`, exactly 0 at `V = C` (so the brim stays
    `h_max` bit for bit), and accurate to a few ulps even when `q` is tiny
    against `R_top³`.
- The cube root is `det::pow(x, 1.0 / 3.0)` (pure-Rust `libm`, the same
  bits on every platform; it only enters the denominator, so its error of
  about 1e-15 relative reaches `d` scaled down). No new `det` function, so
  `kuniumi-math` and its golden machinery are untouched.
- `A` is the grid's own area sum (not `4π`), as in the thin layer now, so
  that for a thin shell `d → (V − C)/(A R_top²)`.
- Both functions require `R_top > 0` (the highest height lies above the
  planet's centre) when they use the shell, and assert it (`# Panics`): a
  NaN must never reach an output (conventions, rule 9). With the bounds of
  `PlanetParams::validate` (`R ≥ 1,500 km`) and heights within ±16 km it
  always holds. The result is finite for every finite volume.
- Doc comments of both functions and of `MAX_OCEAN_VOLUME_M3` (which now
  names this ticket as open) are updated; the design's 海面 section gains a
  paragraph: the water above the highest land is a spherical shell on the
  sphere of radius `R + h_max`, with the formula and the 1,500 km example.

Tests (`stages::crust::tests`), written first:

- `sea_level_on_the_octahedron` (by hand): `R` = 1,000 m, heights −4,000 …
  1,000 m, `R_top` = 2,000 m. A water world with `z` = 1,100 m holds
  `C + (A/3)(2,100³ − 2,000³)`, which is 15,022 units of `a R²` (12,500 at
  the brim plus `2 × 1.261e9 / 1e6` = 2,522); the thin layer would give
  1,420.3 m. Both the sea level of 15,022 units and the volume at 1,100 m
  are checked by hand.
- `no_water_and_water_worlds`: the water worlds (2C, 10C, 1e25 m³ and
  `f64::MAX`) against the test's own formula, `cbrt(R_top³ + 3(V − C)/A)
  − R_top` with std `cbrt` (the naive form, independent of the stage's
  stable one), within 1e-9 relative; every cell under water; the reference
  volume at the level, extended in the test by the shell above the brim,
  holds `V` within 1e-9 relative. Filled to the brim stays `h_max` bit for
  bit; a flat planet's layer is a shell too.
- `volume_and_sea_level_are_monotonic`: the sweep of levels continues
  through the brim to 100 km above it, the volumes up to 200 C; strictly
  growing.
- A new `deep_water_worlds_against_spherical_shells`: the smallest planet
  with the most water (1,500 km, 1e20 m³) and the largest (13,000 km,
  1e20 m³) through the Crust stage (1,000 cells, seed 3): the stage's sea
  level lies within 1e-5 relative of the fully spherical model written in
  the test (`Σ a_i ((R + z)³ − (R + h_i)³)/3 = V`, bisection with std
  functions), and the old thin layer is more than 100 % too high at
  1,500 km. Also Earth's water on the smallest planet. (Implemented with
  a derived bound instead of 1e-5; see Result.)
- Everything else unchanged: `sea_level_holds_the_volume` (only volumes
  below the brim), `land_fraction_of_earth_like_worlds`, the stage test,
  the golden fingerprints and `GOLDEN_SEA_LEVEL_M` of the CLI (Earth with
  Earth's water, not a water world).

World-breaking: yes, but only for water worlds (`V ≥ C`, every cell under
water): their sea level is lower (0.33 % at 13,000 km with the most
water, a factor 2.3 at 1,500 km). Worlds with land keep every bit; no
golden value changes.

Size: about 250 changed lines (code 60, tests 120, design, changelog,
ticket).

### Option A (outline, not chosen)

`PlanetParams::MAX_OCEAN_MEAN_DEPTH_M` (a bound for the owner to choose;
for example 50 km, which keeps M2-08's largest case, a layer of 47.1 km on
13,000 km) and `ocean_mean_depth_m() = V / (4πR²)`; `validate` checks it
after the bulk density, as `ParamError::Inconsistent { name:
"ocean_mean_depth_m", value, fields: ["radius_m", "ocean_volume_m3"] }`,
and `kuniumi-config` names both keys as written. With 50 km a 1,500 km
planet holds at most 1.41e18 m³ (1.06 Earth oceans), and the thin layer
is still about 3 % (1.4 km) too high there; parameter sets that passed
before are rejected (*Changed*, not world-breaking); D-12 and the design's
プロジェクト方針 list the third check.

## Acceptance criteria

- [x] The water above the highest land is a spherical shell: the
      octahedron by hand (15,022 units at `R` = 1 km give 1,100 m, and
      back) and the 10,000-cell worlds' water worlds against the test's
      std-`cbrt` formula within 1e-9 relative, with `water_volume_m3` the
      inverse of `sea_level_m` within 1e-9 relative
      (`stages::crust::tests`).
- [x] Deep water worlds are near physical: 1,500 km and 13,000 km with
      1e20 m³, and 1,500 km with Earth's water, through the Crust stage
      within the error of the bathtub's prisms below the brim of the fully
      spherical model, `|z − z_s| ≤ |C_s − C| / (A R_top²)` with
      `|C_s − C| ≤ (2m + m²) C`, `m = max |h_i| / R`, and within 1e-3
      relative (`deep_water_worlds_against_spherical_shells`); the CLI's
      1,500 km world of the table gives about +1,508,308 m.
- [x] Volume and sea level stay strictly monotonic through the brim, and
      the sea level is finite for every finite volume (`f64::MAX`).
- [x] Worlds with land are unchanged bit for bit: no golden value or
      fingerprint changes, and the tests of volumes below the brim pass
      unchanged.
- [x] The design's 海面 section and the doc comments state the shell;
      `CHANGELOG.md` lists the change under *World-breaking*.
- [x] `mise run ci` passes.

## Verification

```sh
mise exec -- cargo test -p kuniumi-sim crust
mise exec -- cargo run --release -p kuniumi-cli -- generate --seed 1 \
  --planet out/M2-17/small_max.toml --out out/M2-17/small_max
mise run ci
```

Images: `out/M2-17/small_max/hypsometry.png` and `surface_height.png`
(every bar far below the sea-level line, now at about +1,508 km).

## Result

What was done:

- `stages::crust::sea_level_m`: for a water world (`V > 0`, `V ≥ C`) the
  water above the highest height is a spherical shell on the sphere of
  radius `R_top = R + h_max`: `z = h_max + q / (R_top² + R_top b + b²)`
  with `q = (V − C)/(A/3)` and `b = det::pow(R_top³ + q, 1/3)`. Dry worlds
  and worlds with land (bisection) are unchanged; the bisection and the
  brim now call the private `bathtub_m3`, the former loop of
  `water_volume_m3`, bit for bit.
- `stages::crust::water_volume_m3`: unchanged up to the highest height;
  above it `C + A d (R_top² + R_top d + d²/3)`, the same shell, so it is
  the inverse of `sea_level_m` for every volume.
- Both panic, through the private `shell_base_m`, for a water world whose
  highest height lies at or below the planet's centre (`R + h_max ≤ 0`),
  which the bounds of `PlanetParams::validate` exclude, so that no NaN
  can reach an output.
- Docs: the module doc and both functions' doc comments (the formula, the
  stable form, the old thin layer and the 1,500 km example), the doc
  comment of `PlanetParams::MAX_OCEAN_VOLUME_M3` (the new sea level), a
  paragraph and the formula in the design's 海面 section, `CHANGELOG.md`
  (*Changed* and *World-breaking*).

Deviations from the plan, and why:

- Criterion 2: the plan's tolerance against the fully spherical model,
  1e-5 relative, was too tight. The stage keeps the design's prisms below
  the brim, whose error reaches the sea level: at 1,000 cells, seed 3,
  1,500 km with Earth's water, the stage is 8.96 m (2.1e-4) above the
  spherical model (0.83 m at 20,000 cells, seed 1, the refinement's
  table). The test now asserts the bound that this error gives,
  `|z − z_s| ≤ |C_s − C| / (A R_top²)` (both models hold the water above
  the brim in the same shell, whose slope is at least `A R_top²`), with
  `|C_s − C| ≤ (2m + m²) C` checked too, plus 1e-3 relative.
- The octahedron got a second water world by hand (50,500 units at
  2,000 m) besides the planned 15,022 units at 1,100 m.

Evidence:

- Criterion 1, the shell: `sea_level_on_the_octahedron` (15,022 units →
  1,100 m and 50,500 units → 2,000 m within 1e-9 m, and the volumes at
  1,000, 1,100 and 2,000 m by hand within 1e-12 relative; the thin layer
  would give 1,420.3 m); `no_water_and_water_worlds` (2C, 10C and 1e25 m³
  on the 10,000-cell world of seed 1 against
  `reference_water_world_level_m`, std `cbrt` in the direct form, within
  1e-9 relative; `water_volume_m3` and `reference_shell_volume_m3` at the
  level hold `V` within 1e-9 relative; each level is below the old thin
  layer; 1e25 m³ gives a sea surface within 1e-3 of the radius of a ball of
  that volume, 1.337e8 m; a flat planet against the closed form within
  1e-6 m and against `spherical_volume_m3` within 1e-9 relative).
- Criterion 2, deep water worlds:
  `deep_water_worlds_against_spherical_shells` (1,000 cells, seed 3):

  | Planet | Stage | Spherical | Off | Bound | Old thin layer |
  | --- | --- | --- | --- | --- | --- |
  | 1,500 km, 1e20 m³ | 1,508,353.707 m | 1,508,351.350 m | 2.357 m | 9.408 m | 3,533,487.004 m |
  | 1,500 km, 1.335e18 m³ | 42,719.241 m | 42,710.280 m | 8.961 m | 9.408 m | 43,926.458 m |
  | 13,000 km, 1e20 m³ | 44,137.729 m | 44,136.557 m | 1.172 m | 1.179 m | 44,285.341 m |

  CLI, 20,000 cells, seed 1 (`out/M2-17/small_max.toml`, 1,500 km,
  2.5 m/s², 1e11 km³): `sea_level_m` = 1,508,307.904 m, equal to the
  refinement's Python value from the `.kfield` files to the millimetre
  (1,508,307.904 m; the fully spherical model 1,508,307.687 m); before,
  3,533,311.836 m. `out/M2-17/large_max.toml` (13,000 km, 20 m/s²):
  43,956.660 m (Python 43,956.660 m), before 44,103.042 m.
- Criterion 3, monotonic and finite: `volume_and_sea_level_are_monotonic`
  (the volume every 100 m up to 100 km above the brim, and the sea level
  for volumes up to 200 C, strictly growing, on two worlds);
  `no_water_and_water_worlds` (`V = f64::MAX` gives a finite 3.50e102 m,
  against the reference within 1e-9 relative); the M2-08 corner test
  (`the_bounds_of_the_planet_parameters_generate_finite_worlds`, all 128
  corners) passes unchanged.
- Criterion 4, worlds with land unchanged: `git diff main` touches no
  golden value or fingerprint; `generate_golden_fingerprints` (CLI, the
  13 fields and `GOLDEN_SEA_LEVEL_M`) passes; `sea_level_holds_the_volume`,
  `land_fraction_of_earth_like_worlds` and the stage test pass unchanged.
- Criterion 5: design 海面 section, doc comments, `CHANGELOG.md`.
- Criterion 6: `mise run ci` passes (989 tests over debug and release,
  clippy, rustdoc, MSRV 1.85.0, cargo-deny ok).
- Mutation: with the shell based at the datum (`R`) instead of
  `R + h_max`, the octahedron, the water-world and the deep-water tests
  fail; before the implementation, the thin layer failed them too.

Images looked at: `out/M2-17/small_max/hypsometry.png` (every bar at the
far left, about 1,504 to 1,516 km below the black sea-level line at the
far right, with Earth's red ocean-mode line next to it: a water world
whose highest point lies 1,504,086 m deep), `out/M2-17/small_max/surface_height.png`
(uniform deep blue, every cell far below the sea level, the plate pattern
faintly visible), `out/M2-17/large_max/hypsometry.png` (the bars from
−52 km to −38 km below the sea level: the ocean floor's spread, the
continents' spike at −43.7 km; `highest_m` = −38,353.6 m, the shell's depth
above the brim of the Python check).

Review: 1 round.

- Round 1 (APPROVE): the reviewer re-ran `mise run ci` (989 tests), redid
  the octahedron by hand, checked the `(2m + m²) C` bound algebraically and
  the ratio of offset to bound, `(R_top / (R + z))²`, evaluated the shell
  formula with 60-digit decimals on its own CLI worlds (2e-15 relative),
  regenerated both CLI worlds (the same sea levels; `surface_height.kfield`
  byte-identical to the old code's, `world.json` differing only in
  `sea_level_m`), and ran three mutations (the shell based at `R`, `d²/2`
  for `d²/3`, the naive `b − R_top`), each caught. R1-1 nit: "finite for
  every finite volume" fails for `R_top` below about 2e-162 m (0/0 at
  `V = C`) and is wrong above 5.6e102 m, both excluded by `validate`
  (fixed: the doc states the planet-scale range, 1e-100 m to 1e100 m,
  for which it holds). R1-2 nit: "the water above the highest land" for
  worlds without land (fixed: "the highest height" in the module doc and
  in `MAX_OCEAN_VOLUME_M3`). `mise run ci` passes after the fixes
  (989 tests).

Follow-ups: none.

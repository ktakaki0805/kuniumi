# M2-05 Elevation from crust and sea level

- Status: done
- Milestone: M2
- Depends on: M2-04, M2-11
- Crates: kuniumi-sim (stages/crust, pipeline), kuniumi-cli (tests), tools
- Design: 地殻から標高へ
- Decisions: none

## Goal

The Crust stage turns crust into relief and solves the sea level from the ocean
volume, giving the first real elevation map.

## Scope

In:
- Ocean depth from age after Parsons & Sclater (1977), with age converted to
  Myr only inside the formula: `d = 2500 + 350·√t` for young and
  `d = 6400 − 3200·exp(−t/62.8)` for old ocean floor (m). The two published
  curves overlap and switching at 70 Myr would leave a 78 m step, so switch
  where they intersect (t ≈ 26.3 Myr, d ≈ 4294 m), computed once by bisection.
- Continental elevation by Airy isostasy:
  `h = (T − T_ref)·(ρ_m − ρ_c)/ρ_m + h_ref` with ρ_c = 2800 kg/m³,
  ρ_m = 3300 kg/m³, T_ref = 35 km, h_ref = 300 m. Document that gravity
  cancels in Airy isostasy.
- Trenches deepened to about 7–10 km.
- Sea level `z` by bisection so that `Σ max(0, z − h_i) · a_i · R² =
  ocean_volume_m3`; ocean volume 0 (no ocean) and more water than any basin
  holds (water world) behave as documented.
- Outputs: `surface_height` (m) and `World::sea_level_m`.

Out:
- Metrics and histograms (M2-06).

## Implementation notes

Gaps found while refining (the ticket was `ready`), settled inside the
ticket's scope: the trench model, what "as documented" means for no water
and for a water world, how the stage reads its inputs, the reference values
of the formula tests, and the crates touched (the CLI tests change because
`generate` now runs the Crust stage and stops at Climate).

### Heights (`stages::crust`)

All heights are relative to the datum of the reference radius, on which
35 km continental crust stands at +300 m; the sea level is solved afterwards
on the same datum.

- Oceanic cells (`crust_type` 0): `h = −d(t)` with Parsons & Sclater (1977,
  "An analysis of the variation of ocean floor bathymetry and heat flow
  with age", JGR 82(5), 803–827): `d = 2500 + 350 √t` (m) for `t < t*`,
  `d = 6400 − 3200 e^(−t/62.8)` (m) for `t ≥ t*`, `t = age_s / SECONDS_PER_MYR`
  (Myr only inside the formula). `t*` by bisection on `[20, 70]` Myr (the
  first curve is below the second at 20 Myr and above it at 70 Myr; their
  difference grows strictly in between, so the root is unique) until no
  double lies between the ends; `t*` is the upper end. Computed once per
  run by `OceanDepth::new()`. Reference (`tools/crust_reference.py`,
  mpmath): `t* = 26.268197500086…` Myr, `d* = 4293.8378393156…` m.
- Trench cells (`trench` 1, always oceanic): `h = −max(d(t), 8 km)`,
  `TRENCH_DEPTH_M` = 8 km, the middle of the ticket's 7–10 km and the value
  of M2-11's prototype; deeper than the ocean floor of any age (the old
  curve never reaches 6.4 km), so every trench lies at −8 km, 7.0–7.3 km
  below the sea level of Earth's water. Most of Earth's major trenches are 7–10 km deep
  (the deepest, the Mariana trench, about 10.9 km); an age-dependent trench
  is left to M6-03.
- Continental cells (`crust_type` 1): Airy isostasy,
  `h = (T − T_ref)(ρ_m − ρ_c)/ρ_m + h_ref` with `ρ_c` = 2,800 kg/m³,
  `ρ_m` = 3,300 kg/m³, `T_ref` = 35 km (`CONTINENTAL_THICKNESS_M`),
  `h_ref` = 300 m. From the balance of the columns at the depth of the
  deepest root, `ρ_c g (h + T_ref + b) = ρ_c g T_ref + ρ_m g b` with root
  `b = T − T_ref − h`: `g` cancels, so `surface_gravity_m_s2` is not an
  input (the stub's doc said otherwise; corrected). 15 km → −2,730.30 m,
  35 km → +300 m, 70 km → +5,603.03 m.

### Sea level

`V(z) = Σ max(0, z − h_i) a_i R²` (m³), summed in cell order; a bathtub
model: every cell below `z` holds water, connected to the open ocean or not
(enclosed basins and the pockets of M2-10; inland water is M5-02).

- `V = 0` (no ocean): `z` = the lowest height, where the first water would
  collect; no cell is under water (depth `z − h_i` ≤ 0 everywhere).
- `0 < V < V(h_max)`: bisection on `[h_min, h_max]` keeping
  `V(lo) < V ≤ V(hi)`, halving until no double lies between the ends (at
  most 100 halvings); `z` = the upper end, the lowest level found that holds
  at least `V`.
- `V ≥ V(h_max)` (water world): every cell under water, `z` in closed form,
  `h_max + (V − V(h_max)) / (Σ a_i R²)`: the rest of the water as a layer
  over the whole planet. Finite for every valid volume with Earth's radius
  (bounds on the radius: M2-08).

### API

```rust
pub const SURFACE_HEIGHT: &str = "surface_height";       // m
pub const CRUST_DENSITY_KG_M3: f64 = 2800.0;
pub const MANTLE_DENSITY_KG_M3: f64 = 3300.0;
pub const REFERENCE_HEIGHT_M: f64 = 300.0;
pub const TRENCH_DEPTH_M: f64 = 8000.0;
pub struct OceanDepth { /* switch age */ }
impl OceanDepth {
    pub fn new() -> Self;                  // bisection, also Default
    pub fn switch_age_s(&self) -> f64;
    pub fn depth_m(&self, age_s: f64) -> f64;
}
pub fn continental_height_m(thickness_m: f64) -> f64;
pub fn surface_heights_m(continental: &[bool], thickness_m: &[f64],
                         age_s: &[f64], trench: &[bool]) -> Vec<f64>;
pub fn water_volume_m3(grid: &SphereGrid, heights_m: &[f64],
                       radius_m: f64, level_m: f64) -> f64;
pub fn sea_level_m(grid: &SphereGrid, heights_m: &[f64],
                   radius_m: f64, volume_m3: f64) -> f64;
pub fn run(world: &mut World, progress: &mut dyn FnMut(f64)) -> Result<(), SimError>;
```

`run` reads `crust_type`, `crust_thickness`, `crust_age` and `trench` from
the world (it panics, documented, if the Tectonics stage has not added
them: the runner always runs it first), adds `surface_height` (m) and sets
`World::sea_level_m`; no progress of its own. `pipeline::implementation`
returns it for `StageId::Crust`, so `generate` stops at Climate. Larger and
smaller values by explicit comparisons (determinism rule 2).

### Tests

In `stages/crust/tests.rs` (std `sqrt`/`exp` allowed as references):
formula values by hand and from `tools/crust_reference.py` (depth(0) =
2,500 m exactly, 1, 4, 16, 25 Myr exactly, 30, 60, 100, 180 Myr within
1e-9 m; `t*` within 1e-9 Myr; the two curves within 1e-6 m at `t*` and the
depth at the doubles next to it within 1e-6 m; strictly increasing on
0–180 Myr in steps of 0.01 Myr and at `t* ± 1e-6` Myr; Airy at 15, 35,
70 km); heights by cell kind; the stage's heights on real worlds against a
reference written with std functions (1e-6 m); sea level on the octahedron
by hand; `|V(z) − V| / V < 1e-9` against a volume summed in the test, on
real worlds for volumes from 1e12 m³ to just below the capacity;
monotonicity of `V(z)` and of `z(V)`; `V = 0`, `V` = the capacity, ten
times it and `f64::MAX`; the stage's fields, sea level and independence of
gravity; land fraction of Earth's parameters for seeds 1–10 at 10,000 cells
(20–40 %, the numbers reported, also at 100,000 and 650,000 cells by a
scratch tool). The CLI tests: the new field and stage, the error naming
`climate`, the range of `surface_height` and the sea level on the
1,000-cell world, and the golden values (`surface_height` and the sea level
in `world.json`, 10,000 cells, seed 1) recorded after the tests above pass
and the maps are viewed.

`surface_height.png` is drawn with viridis, as every quantity of a world
directory (M2-01, Out). For the visual check a scratch tool also draws the
heights with `Colormap::Elevation` at the solved sea level; drawing
`surface_height` that way in world directories is the follow-up draft
M2-14.

## Acceptance criteria

- [x] Formula tests: depth(0) = 2500 m; depth is continuous across the switch
      (difference below 1e-6 m) and increasing in t; thickness 35 km → 300 m,
      70 km → 5603 m.
- [x] Bisection: `|volume(z) − V| / V < 1e-9`; monotonicity test; V = 0 and a
      huge V behave as documented.
- [x] Earth defaults, seeds 1–10 at 10,000 cells: land fraction between
      20% and 40% (also reported at 100,000 and 650,000 cells).
- [x] `surface_height.png` (650,000 cells and the 10,000-cell golden case,
      seed 1) viewed and described (continents, ridges and trenches
      visible).
- [x] Golden fingerprints for 10,000 cells, seed 1: `surface_height` and
      the sea level.
- [x] `mise run ci` passes.

## Verification

```sh
cargo run --release -p kuniumi-cli -- generate --seed 1 --out out/M2-05/650k
cargo run --release -p kuniumi-cli -- generate --seed 1 --cells 10000 --out out/M2-05/10k
mise exec -- cargo test -p kuniumi-sim --release stages::crust -- --nocapture
```

## Result

What was done:
- `kuniumi-sim`, `stages::crust` (was a doc stub): `OceanDepth` (Parsons &
  Sclater depth, the switch `t*` by bisection once per run),
  `continental_height_m` (Airy, gravity cancels, documented with the
  balance of the columns), `surface_heights_m`, `water_volume_m3`,
  `sea_level_m` (no water, bisection, water world in closed form), `run`,
  and the constants `SURFACE_HEIGHT`, `CRUST_DENSITY_KG_M3`,
  `MANTLE_DENSITY_KG_M3`, `REFERENCE_HEIGHT_M` and `TRENCH_DEPTH_M`;
  tests in `stages/crust/tests.rs` (15). `pipeline::implementation` runs
  the stage for `StageId::Crust`; the documentation of `generate`,
  `implementation` and `StageId::Crust` and the tests of `lib.rs` and
  `pipeline` follow (the run stops at Climate).
- `kuniumi-cli`: tests only (the new field and stage, the error naming
  `climate`, the heights of the 1,000-cell world by kind of cell and the
  sea level from `world.json`, the golden values).
- `tools/crust_reference.py` (mpmath) and its line in `tools/README.md`;
  `CHANGELOG.md` (Added, Changed); the follow-up draft M2-14 in the README.

Deviations from the plan: none in the method. The constants of the
formula tests are the nearest doubles of the mpmath values (clippy rejects
literals with more digits).

Evidence:

1. Formula tests: `ocean_depth_by_hand` (depth(0) = 2,500 m bit for bit;
   1, 4, 16 and 25 Myr exactly 2,850, 3,200, 3,900 and 4,250 m; 30, 60,
   100 and 180 Myr within 1e-9 m of `tools/crust_reference.py`);
   `the_depth_curves_switch_where_they_meet` (`t*` within 1e-9 Myr of
   26.268197500086…; the two curves, with std functions, within 1e-6 m at
   `t*` and at the depth 4,293.84 m; one second either side of the switch
   the depth comes from the young and from the old curve and differs by
   less than 1e-6 m); `depth_increases_with_age` (strictly, every 0.01 Myr
   from 0 to 180 Myr, and at `t* − 1e-6`, `t*`, `t* + 1e-6` Myr; within
   [2,500, 6,400) m); `airy_heights_by_hand` (35 km → 300 m bit for bit,
   70 km → 184,900/33 = 5,603.03 m, 15 km → −2,730.30 m, 36 km →
   +151.5 m, within 1e-9 m); `heights_by_kind_of_cell` (continental cells
   ignore the trench flag, oceanic cells their thickness, trenches −8 km
   at 0, 60 and 180 Myr); `heights_match_a_reference_on_real_worlds`
   (seeds 1–10 at 10,000 cells: every cell within 1e-6 m of the formulas
   written with std `sqrt`/`exp` and their own bisection; continental
   cells, trenches and floor on both curves all occur).
2. Bisection: `sea_level_on_the_octahedron` (by hand: `z` = −2,000, 250,
   1,000, 1,100 and −3,999 m within 1e-9 m for `Σ (z − h_i)` = 3,000,
   9,000, 12,500, 13,100 and 1 m; `V = 0` gives −4,000 m bit for bit; the
   volume by hand at five levels within 1e-12 relative);
   `sea_level_holds_the_volume` (seeds 1–3 at 10,000 cells, `V` from 1e12
   m³ to `(1 − 1e-9)` of the capacity and Earth's 1.335e18 m³:
   `|V(z) − V| / V < 1e-9` against `reference_volume_m3`, a sorted sum with
   Neumaier's compensation; the stage's own sea level too);
   `volume_and_sea_level_are_monotonic` (seeds 1–2: `V(z)` exactly 0 up to
   the lowest height and strictly growing in 1,000 steps to 100 m above
   the highest; `z(V)` strictly growing for 200 volumes up to twice the
   capacity); `no_water_and_water_worlds` (`V = 0`: `z` the lowest height
   bit for bit, volume exactly 0, no cell below it; `V` = the capacity: `z`
   the highest height bit for bit; 2 and 10 times the capacity and 1e25 m³:
   every cell under water, `z = h_max + (V − C)/(A R²)` within 1e-12
   relative, volume within 1e-9; `f64::MAX`: finite, `V/(A R²)` within
   1e-9; a flat planet: a layer);
   `the_stage_adds_surface_height_and_sets_the_sea_level` also runs the
   stage with no water (the lowest height), 1e21 m³ (a water world), Mars's
   radius with Earth's water (a water world, 9 km deep on average) and with
   the water scaled to Mars's area (partly flooded, within 1e-9), each bit
   for bit equal to `sea_level_m` with the world's own radius and volume
   (it fails when `run` hard-codes Earth's radius or volume); the two panics
   (`heights_need_one_value_of_each_per_cell`,
   `the_sea_level_needs_a_height_per_cell`).
3. Land fraction (`land_fraction_of_earth_like_worlds`, Earth's
   parameters, seeds 1–10 at 10,000 cells): 28.22–30.11 % (29.61, 28.75,
   29.41, 30.11, 28.98, 28.96, 29.85, 30.06, 28.38, 28.22), mean 29.23 %;
   sea level −716 to −982 m. A scratch tool on the real pipeline:
   100,000 cells, seeds 1–5, 28.23–29.99 % (mean 29.46 %); 650,000 cells,
   seeds 1–10, 28.74–31.42 % (mean 29.98 %), sea level −688 to −891 m. At
   650,000 cells: mean land elevation above the sea 0.97–1.20 km (Earth
   0.84 km), mode +0.95 to +1.15 km (Earth about +0.1 km; M2-11, finding
   5), mean ocean depth 3.67–3.82 km (Earth 3.7 km), submerged continental
   crust 8.6–11.3 % of the surface, trenches 7.1–7.3 km below the sea,
   highest point 3.3–6.5 km above it.
4. Maps (`kuniumi generate --seed 1 --out out/M2-05/650k` and `--cells
   10000 --out out/M2-05/10k`; elevation maps at the solved sea level by
   the scratch tool in `out/M2-05/elevation/`), viewed:
   - `650k/surface_height.png` (viridis, −8 km to +5.6 km): continents
     uniform green (+300 m), fading to blue over the few hundred
     kilometres of their thinned margins; arcs as lighter green bands
     along the active coasts (the south-west coast of the eastern
     continent, about 105° E to 155° E, 0° to 30° S; the north coast of
     the north-western continent); a yellow collision plateau (+5.6 km)
     near 35° S 130° W; trenches as dark purple lines one cell wide along
     the convergent boundaries (from about 137° W 56° N to 115° W 22° N
     and around the central basin near 50° W–0°, 0°–25° S); ridges as
     pale lines in the blue ocean floor (the meridional ridge from about
     55° W 20° S to 30° W 75° S, one near 110° E–135° E, 60° S), the floor
     darkening (deeper) away from them, darkest in the north.
   - `elevation/elevation_650000_1.png` (land and sea): continents above
     the sea inside rings of light-blue shelves on every passive coast,
     none along the arcs (tan) with a dark trench just offshore; the
     collision plateau white and brown; pale ridge lines in the deep blue;
     an age step at a plate boundary in the northern ocean (floors of
     different plates); 30.15 % land.
   - `10k/surface_height.png` and `elevation/elevation_10000_1.png` (the
     golden case, 2° cells): the same features: trench lines one or two
     cells wide, arcs on the coasts behind them, the collision plateau near
     12° S 5° E, broad pale ridge bands (south of the central continent
     at 45°–65° S, and at 40°–70° W, 45°–70° N), shelves round the passive
     coasts and two submerged rifts inside the central continents (near
     12° N 0° and 40° S 80° E); 29.61 % land.
5. Golden values (`tests/cli.rs::generate_golden_fingerprints`, 10,000
   cells, seed 1), recorded after the tests above passed and the maps were
   viewed: `surface_height` `0xE587_6747_8054_1DE7` and the sea level
   −781.7637678179066 m (`GOLDEN_SEA_LEVEL_M`, compared bit for bit with
   `world.json`); the twelve earlier values unchanged.
6. `mise run ci` passes: fmt, clippy `-D warnings`, the tests in debug and
   release, MSRV 1.85.0, cargo-deny. `cargo doc -p kuniumi-sim` gives no
   warnings. Test-first: against stubs (commit "tests of the Crust stage
   against stubs") 14 tests failed: the 12 of `stages::crust` that check
   values, `generate_returns_the_partial_world_when_allowed` and the CLI
   test `generate_writes_the_world_so_far_then_exits_with_code_2`.

Timing (release, this machine): the Crust stage takes 0.03 s at 650,000
cells; `kuniumi generate --seed 1` (650,000 cells) 5.1 s and 204 MB (5.0 s
and 203 MB after M2-11); at 100,000 cells 1.1 s and 50 MB.

Review:
- Round 1 (APPROVE): the reviewer re-ran `mise run ci`, the reference
  script (and an independent mpmath root), the land fractions at 10,000,
  100,000 and 650,000 cells, the volume of the solved sea level from the
  written `.kfield` files (1.5e-15 and 7.5e-15 relative), the stub commit,
  the maps and the timing, and mutated the densities, the trench, the
  switch and `R²` (each caught). R1-1 minor, no test checked that `run`
  passes the world's radius and volume (fixed: four more worlds in
  `the_stage_adds_surface_height_and_sets_the_sea_level`; hard-coding
  either value now fails). R1-2 minor, the Crust stage narrows the radii
  that give finite worlds (an infinite sea level below about 1e-145 m with
  Earth's water) and M2-08 did not say so (fixed: M2-08's notes, and its
  criterion now runs through the Crust stage). R1-3 minor, M2-11's
  changelog entry gave other land fractions for "the formulas of M2-05"
  (fixed: "without its deepened trenches"). R1-4 nit, 14 tests failed
  against the stubs, not 13 (fixed). R1-5 nit, the bracket width of
  `MAX_HALVINGS` for ±16 km is 2.5e-26 m (fixed). R1-6 nit, the ticket
  depends on M2-11 too (fixed, here and in the README).

Follow-ups: M2-14 (draft, elevation maps of `surface_height` in world
directories).

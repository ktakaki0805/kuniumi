# Changelog

All notable changes to kuniumi are recorded here, following
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Nothing is released
before 1.0, so everything sits under *Unreleased*.

Besides Added, Changed, Fixed and Removed, each version has a **World-breaking**
section: every change after which an existing seed produces a different world
(`docs/conventions.md`, "World-breaking changes").

## [Unreleased]

### Added

- Cargo workspace: `kuniumi-math`, `kuniumi-sphere`, `kuniumi-world`,
  `kuniumi-sim`, `kuniumi-io`, `kuniumi-cli`.
- Deterministic random numbers (xoshiro256\*\* with named streams), `Vec3`, SI
  unit helpers.
- `SphereGrid` data structure with full validation and an octahedron test grid.
- `PlanetParams` with Earth defaults; tidally locked and stagnant-lid regimes
  reserved for v2.
- PNG writer and the `.kfield` field format with fingerprints.
- CLI commands `version`, `params`, `generate` (stub) and `test-image`.
- Determinism lints (`clippy.toml`) and cargo-deny configuration.
- Development environment in `mise.toml`: pinned Rust toolchains (development
  and MSRV 1.85.0), cargo-deny, and the tasks `format`, `check` and `ci`.
- CI workflow running the same mise tasks, with a weekly cross-platform
  determinism check.
- Design document, conventions, decisions, roadmap and tickets.
- Claude Code skills `implement-ticket` and `implement-all-tickets`, and the
  subagents `ticket-reviewer` and `ticket-implementer`.
- `kuniumi_math::det`: deterministic `sin`, `cos`, `sin_cos`, `tan`, `atan`,
  `atan2`, `asin`, `acos`, `exp`, `ln` and `pow` over the pure-Rust `libm`
  0.2.16 (pinned exactly), with accuracy, reference golden, special-value and
  output-fingerprint tests (M1-01).
- `tools/det_reference.py` and `tools/det_error.py`, the first development
  scripts (mpmath reference values and error measurement of `det`); Python
  3.14.10 pinned in `mise.toml`, packages in `tools/requirements.txt` (M1-01).
- `kuniumi_sphere::coords`: `from_lat_lon`, `lat_lon` (longitude in (−π, π],
  0 at the poles), the local frame `east` and `north` (defined at the poles as
  the limit along longitude 0), and `great_circle_rad`; and
  `Vec3::rotate_about` (Rodrigues' formula, right-hand rule, length kept
  within 1e-15 relative), all built on `det` and checked by output
  fingerprints; the example `rotate_length` measures the rotation's length
  change (M1-02).
- `kuniumi_sphere::predicates::orient3d` and `Orientation`: the exact sign of
  `det[b − a; c − a; d − a]` (`Positive` for a right-handed tetrahedron, the
  opposite of Shewchuk's sign), over Shewchuk's adaptive predicates from the
  `robust` crate 1.2.0 (pinned exactly); exact for coordinates that are zero or
  of magnitude in [2⁻²⁰⁰, 2²⁰⁰], checked against exact integer evaluations
  (M1-03).
- `kuniumi_sphere::points::fibonacci`: `n` jittered Fibonacci points on the
  unit sphere, the cell centres of the global grid: the golden-angle spiral
  `z = 1 − (2i + 1)/n`, `λ = i · π(3 − √5)`, each point moved in its tangent
  plane by an offset uniform over a disc of radius `jitter × √(4π/n)` and
  normalised; default jitter 0.35, clamped to [0, 0.45]; random numbers from
  the stream `grid/points` (two per point); `jitter = 0` gives the plain
  spiral bit for bit. Also `mean_spacing_rad`, and the example
  `point_spacing`, which measures the closest pair over many seeds (M1-04).
- `kuniumi_sphere::delaunay::triangulate`: the spherical Delaunay
  triangulation of unit vectors as their 3D convex hull, built in-house by
  randomised incremental construction with conflict lists (insertion order from
  the stream `grid/insertion`, Fisher–Yates). Returns `Triangulation`:
  `2n − 4` triangles counter-clockwise from outside and the triangle across
  each edge, in a canonical order (smallest index first, sorted), so the result
  depends only on the points, not on the seed. Every decision uses the exact
  `orient3d`; exact ties (co-circular points) are broken by the point indices
  with Simulation of Simplicity (Edelsbrunner and Mücke 1990), so degenerate
  inputs such as the cube give a valid triangulation. Rejects too few or too
  many points, non-unit points, duplicates, points that are not hull vertices,
  points closer than 1e-12 and points within a closed hemisphere
  (`DelaunayError`). 650,000 points take about 2.5 s in release (M1-05).
- `SphereGrid::fibonacci(cells, seed)`, the global grid: the jittered
  Fibonacci points (stream `grid/points`) and their Voronoi cells; and
  `SphereGrid::from_points(points, seed)` for any point set: the dual of the
  Delaunay triangulation, with circumcentres `normalize((b − a) × (c − a))` as
  Voronoi vertices, neighbour rings from walking the triangle fan around each
  centre (counter-clockwise, no angles), cell areas as sums of signed solid
  angles (Van Oosterom and Strackee 1983) of the fan triangles split at the
  midpoints of the Voronoi edges, and per neighbour pair the Voronoi edge
  length and the centre distance in radians, bitwise symmetric
  (`SphereGrid::edge_lengths_rad`, `SphereGrid::centre_distances_rad`).
  `BuildError` (including `Degenerate` for a cell without positive area and
  `Invalid` should a built grid fail `validate()`, which `from_points` always
  runs) and `MAX_CELLS` (CSR offsets in `u32`). 650,000 cells take about
  2.8 s in release; their areas sum to 4π within 1.5e-14 relative (M1-06).
- `SphereGrid::nearest(p, hint)`: the cell nearest to a point, the one with
  the largest computed `p · c`, ties to the lower index (the brute-force
  answer, whatever the hint), by a steepest-ascent walk on the neighbour
  graph from the hint and a search through near ties (within 1e-12
  relative) around its end, which also settles ties between co-circular
  cells that are not neighbours. `SphereGrid::locator()` builds a `Locator`
  of about `n/4` latitude/longitude buckets that supplies start cells for
  queries without a hint (`Locator::start`, `Locator::nearest`). A 2048 ×
  1024 map queried row by row on 650,000 cells takes about 0.09 s in
  release (M1-07).
- `kuniumi_io::colormap`: `Colormap::{Viridis, Diverging { centre },
  Elevation { sea_level }, Categorical}` and `Range::{Auto, Fixed}` (auto:
  the finite minimum and maximum). Viridis from 52 control points of the CC0
  table (within 0.62 of an 8-bit unit of the full table); diverging
  ColorBrewer RdBu (11 classes, blue low, white at the centre, the same scale
  on both sides); elevation with a hard break at sea level (blues below,
  greens and browns above, each side scaled separately); categorical
  ColorBrewer Paired (12 classes). Non-finite values are magenta
  (`NO_DATA`). Linear interpolation in `f64`, rounded half away from zero
  (M1-08).
- `kuniumi_io::render`: `equirect` and `orthographic` (globe) maps of a
  per-cell field, and `equirect_cells` and `orthographic_cells`, the cell
  under each pixel; one `SphereGrid::nearest` query per pixel from the
  previous pixel's cell. A 2048 × 1024 map of 650,000 cells renders in about
  0.1 s in release. The example `render_maps` writes test maps of every
  colour map (M1-08).
- `tools/colormap_tables.py` (the colour tables from their sources and test
  reference values) and `THIRD_PARTY.md` (M1-08).
- `kuniumi_world::Settings { cells }`, the generation settings, separate from
  the planet parameters: 650,000 cells by default, 1,000 to 715,827,884
  (`MAX_CELLS` of the grid); `Settings::validate` names an out-of-range
  setting with `ParamError::OutOfRange` (M1-09).
- The Grid stage, `kuniumi_sim::stages::grid::run`: `SphereGrid::fibonacci`
  with the settings' number of cells, random streams `grid/points` and
  `grid/insertion`; a golden test checks that 10,000 cells, seed 1 give
  exactly the fingerprints recorded in M1-06 (M1-09).
- `kuniumi grid --seed N --out DIR [--cells N] [--width PX]`: builds the
  grid and writes `cell_area.png` (equirectangular, area relative to the
  mean, diverging map), `cell_area_globe.png` (the same on a globe centred
  on 45° N, 180° E), `neighbors.png` (ring size, categorical) and
  `stats.txt` (cell count, area sum and its error against 4π, smallest,
  mean and largest area and their ratio, mean spacing, ring-size histogram,
  neighbour entries against `6n − 12`, build time) (M1-09).
- Conformance tests of the floating-point operations that rule 2 of the
  determinism conventions allows directly (`kuniumi-math/src/exact_tests.rs`:
  rounding, `%`, sign operations, comparisons, `min`/`max`/`clamp`, `sqrt`,
  `mul_add`, integer conversions, decimal formatting, at their hard cases
  with hand values), and a check that the reduction `(i · g) % 2π` in
  `points` is exact for every point of a 650,000-point spiral (M1-10).
- Example `jitter_survey` in `kuniumi-sphere`: for a grid size, a seed and
  a list of jitters, the spread of cell areas, the closest centres, the
  explicit diffusion limit, the ring sizes and hexagon shares by latitude,
  and the path bias of steepest-descent drainage on a uniform slope, with
  tests of the path bias against hand values (M1-11).
- A test of the example `jitter_survey`: its steepest-descent receiver takes
  the longer of two steps when that step is the steeper, which a receiver
  chosen by the drop over the squared step length fails; the example's
  documentation lists all six of its tests (M1-12).
- The stage runner `kuniumi_sim::pipeline::run`: validates the inputs, runs
  the Grid stage, then every later stage of `StageId::PIPELINE` in order from
  a table of implementations (`pipeline::implementation`, `StageFn`), with
  progress 0, the stage's own reports and 1 per stage, a debug check that
  every output is finite and every class field holds class ids, and the
  stage's name appended to `World::stages`;
  it stops at the first missing stage with `SimError::NotImplemented`, or
  with `Settings::allow_partial` returns the world built so far.
  `kuniumi_sim::next_stage` names the stage a partial world stopped at;
  `SimError::Field` reports a stage output the world rejects (M2-01).
- `World::add_field`, `World::add_class_field` and `World::field_mut`;
  `FieldKind::{Quantity, Class}` and `Field::kind`; `FieldError` (invalid
  name, duplicate, wrong length, a class value that is not an integer in
  [0, 2⁵³)); field names are `[a-z][a-z0-9_]*` and not a Windows device
  name (`con`, `nul`, `com1`, …). `World` records its
  `settings` and the `stages` that have run (M2-01).
- The Grid stage adds the fields `cell_area` (sr) and `neighbor_count`
  (class) (`stages::grid::add_fields`) (M2-01).
- `kuniumi_io::world_dir::write(dir, &world)`: a `.kfield` file and a
  2048 × 1024 equirectangular PNG per field (categorical palette for class
  fields, viridis otherwise), then `world.json` (format version 1: generator
  version, seed, settings, planet parameters, sea level, stages, fields with
  units and kinds), through serde mirror structs (`WorldJson`) with `serde`
  1.0 and `serde_json` 1.0 (feature `float_roundtrip`, so every `f64` parses
  back bit for bit); the same world gives the same bytes; an existing
  `world.json` is removed before the fields are written, so a failed rewrite
  leaves none behind (M2-01).
- The Tectonics stage (one-shot plates), `kuniumi_sim::stages::tectonics`:
  `partition` divides the cells into `Settings::plate_count` connected
  plates by growth from seed cells at least `0.4 √(4π/K)` rad apart, with a
  growth rate per plate uniform in [1, 3) and a slowness per cell
  `exp(3 (u − ½))` (multi-source Dijkstra, heap ordered by priority, cell
  and plate; stream `tectonics/partition`); `euler_poles` gives each plate a
  uniformly random axis and an angular speed at which its fastest cell
  moves at a log-uniform 1 to 10 cm/yr (stream `tectonics/poles`);
  `EulerPole::velocity_m_s` is `ω × r`. The stage adds the fields
  `plate_id` (class), `plate_velocity_east` and `plate_velocity_north`
  (m/s) (M2-02).
- `Settings::plate_count` (default 12, valid 2 to 64, rejected by name
  otherwise), recorded in `world.json` as `settings.plate_count` (M2-02).
- Continents and plate boundaries in the Tectonics stage
  (`kuniumi_sim::stages::tectonics`): `continents` grows one continental
  blob per plate (random seed cell, growth rate uniform in [1, 3), slowness
  per cell `exp(3 (u − ½))`, confined to its plate; stream
  `tectonics/continents`) on plates taken in random order until they cover
  twice the continental area, and stops when the area reaches
  `PlanetParams::continental_crust_fraction` (within one cell);
  `boundaries` classifies every cell next to another plate as convergent,
  divergent or transform (`BoundaryType`) from the relative velocity of the
  two plates across a boundary normal summed along the boundary within
  0.08 rad (`NORMAL_WINDOW_RAD`), with a 45° threshold, and gives its
  convergence rate (m/s, positive when the plates approach) and the plate
  across. The stage adds the fields `crust_type` (class: 0 oceanic, 1
  continental), `boundary_type` (class: 0 none, 1 convergent, 2 divergent,
  3 transform) and `convergence_rate` (m/s) (M2-03).
- The one-shot crust in the Tectonics stage
  (`kuniumi_sim::stages::tectonics::crust`): `nearest_sources` (multi-source
  Dijkstra inside each plate, edges of centre distance × radius, heap
  ordered by distance, cell and source); ocean-floor age = distance to the
  nearest divergent cell of the plate ÷ its half spreading rate, capped at
  180 Myr (continental crust 180 Myr); oceanic crust 7 km; continental crust
  35 km, thickened at collisions by up to 35 km (`min(1, c / 5 cm/yr)`,
  full to 200 km from the front, cosine taper to 0 at 400 km, both plates)
  and into arcs of up to 50 km peaking 200 km inland of subduction fronts;
  rock uplift `0.1 c` (collisions) and `0.03 c` (arcs) times the same
  profiles, at most 10 mm/yr, on continental crust only; trenches on the
  oceanic side of subduction boundaries (the older floor where ocean meets
  ocean). The stage adds the fields `crust_thickness` (m), `crust_age` (s),
  `uplift_rate` (m/s) and `trench` (class) (M2-04).
- The Crust stage (`kuniumi_sim::stages::crust`), the first relief: ocean
  floor at the depth of its age after Parsons & Sclater (1977),
  `2500 + 350 √t` and `6400 − 3200 e^(−t/62.8)` m (`t` in Myr), switched
  where the two curves meet (26.27 Myr, 4,294 m; `OceanDepth`, found by
  bisection); trench cells at 8 km below the datum (`TRENCH_DEPTH_M`);
  continental crust by Airy isostasy, `(T − 35 km)(ρ_m − ρ_c)/ρ_m + 300 m`
  with 2,800 and 3,300 kg/m³ (`continental_height_m`; gravity cancels, so
  the surface gravity is not used); the sea level that holds
  `ocean_volume_m3` (`sea_level_m`, bisection on the volume
  `water_volume_m3`; with no water the lowest height, with more water
  than the basins hold a layer over the whole planet). The stage adds the
  field `surface_height` (m) and sets `World::sea_level_m`. With Earth's
  parameters the land fraction is 28.2–30.1 % at 10,000 cells (seeds
  1–10) and 28.7–31.4 % at 650,000 cells; it takes 0.03 s at 650,000
  cells (M2-05).
- `tools/crust_reference.py`: mpmath reference values of the depth curve
  and Airy isostasy for the tests of the Crust stage (M2-05).
- `kuniumi_sim::metrics`, the first validation metrics of the design
  (検証方法), area-weighted over the cells with heights relative to the
  sea level (a cell below it is ocean, every other cell land): the land
  fraction, the mean land elevation and ocean depth, the lowest and
  highest elevation, percentiles (`percentiles_m`, the lower quantile, at
  `PERCENTILES`), the elevation histogram in 100 m bins from the sea level
  (`Hypsometry`, `HYPSOMETRY_BIN_M`) with its land and ocean modes, the
  valley between them and whether the two peaks are separate (the valley
  holds at most `VALLEY_RATIO` = ½ of the smaller mode), all in
  `ElevationMetrics`, and Earth's approximate values (`metrics::earth`).
  With Earth's parameters at 650,000 cells (seeds 1–10) every histogram
  has two separate peaks: land modes +0.95 to +1.15 km (Earth about
  +0.1 km) and ocean modes −5.15 to −5.35 km (Earth about −4.5 km) (M2-06).
- `kuniumi_io::chart::histogram`: bar charts drawn into RGB pixel buffers,
  without text, with grid lines, axes and vertical markers (M2-06).
- `tools/hypsometry_reference.py`: the metrics of a world directory with
  exact rational arithmetic, for the tests of `metrics.txt` (M2-06).
- `kuniumi-config`: settings files in TOML (D-15), through `toml` 1.2.1 and
  `serde` (not pinned). `load` reads `[planet]` (`PlanetParams`) and
  `[generation]` (`Settings`): each quantity under its SI key or one
  human-unit key (`radius_km`, `rotation_period_h`, `obliquity_deg`,
  `orbital_period_days`, `ocean_volume_km3`, `surface_pressure_bar`), not
  both; missing keys take Earth's values and the default settings; values
  are checked by `validate`. `ConfigError` names the key and the number as
  written, or the line and column for malformed TOML, unknown tables and
  keys and wrong types. `write` writes every field with its SI key, and
  `load` reads it back bit for bit (M2-07).
- `kuniumi_math::units`: `hours_to_s`, `days_to_s`, `km_to_m`, `km3_to_m3`,
  `bar_to_pa` and their constants (M2-07).
- `kuniumi generate --planet FILE` reads a settings file; `--cells` replaces
  the file's `cells` and no longer has a default of its own (M2-07).
- The mise task `doc`: `cargo doc --workspace --no-deps --locked` with
  warnings as errors (`RUSTDOCFLAGS="-D warnings"`). `mise run check`, and
  so `mise run ci`, runs it, and so does the `lint` job of the CI workflow
  (M2-09).

### Changed

- `SphereGrid::from_parts` takes the edge lengths and centre distances as
  two more arrays; `validate()` checks them (lengths, ranges, bitwise
  symmetry, centre distances against the centres; `GridError::EdgeLength`,
  `GridError::CentreDistance`), rejects non-finite areas and positions and
  coordinates outside `orient3d`'s exact range, decides the counter-clockwise
  order of the rings exactly with `orient3d`, and checks the area sum
  against 4π within 1e-12 relative (was 1e-9), summed with Neumaier's
  compensation. `SphereGrid::octahedron()` gives each cell the double nearest
  to 2π/3 (one unit in the last place above `4.0 * PI / 6.0`) and the hand
  values of the new arrays (M1-06).
- `kuniumi_io::png` compresses through the `png` crate 0.18 (default
  deflate level, adaptive filters) instead of the in-house stored-deflate
  encoder: a 2048 × 1024 map takes 0.4–1.0 MB instead of 6.3 MB. The pixel
  type `png::Rgb8` is now `kuniumi_io::Rgb8` (M1-08).
- `kuniumi_sim::generate(params, &settings, seed, progress)` takes the
  settings, validates them before any stage, runs the Grid stage (progress
  0 and 1), and returns `SimError::NotImplemented { stage: Tectonics }`
  instead of stopping before the grid; new `SimError::Grid` for a grid that
  cannot be built, and `SimError` reports its source (M1-09).
- The CLI parses its arguments with `clap` 4.6 (derive, subcommands) instead
  of a hand-written parser; `version`, `params`, `generate`, `test-image`,
  `help`, `--help` and `--version` keep working, every subcommand has its own
  `--help`, and every error, usage errors included, exits with code 2.
  `generate` takes `--cells` (default 650,000) and now prints the Grid
  stage's progress before the error that the tectonics stage is not
  implemented yet (M1-09).
- Determinism rule 2 (`docs/conventions.md`) names every floating-point
  operation library code may use directly: those with a uniquely defined
  result (correctly rounded arithmetic, `sqrt`, `mul_add`, integer to
  `f64`; exact sign operations, comparisons, classifications, bit
  reinterpretation, rounding, `fract`, `%`, `f64` to integer; `min`, `max`
  and `clamp` only where the sign of a zero result cannot reach the output;
  decimal formatting for printed output), and that a NaN, whose sign and
  payload are not defined, must not reach an operation that exposes its
  bits or an output. `clippy.toml` also bans
  `f64::to_degrees`, `f64::to_radians` and `f64::abs_sub`, whose precision
  std leaves unspecified. No library code changes (M1-10).
- The documentation of `points::DEFAULT_JITTER` records the owner's
  decision to keep 0.35, with its measured effect on drainage, cell areas
  (largest to smallest area about 3.4, accepted), ring sizes and explicit
  time steps. No values change (M1-11).
- `kuniumi generate` writes the world built so far into `--out` (one
  `wrote` line per file) before it exits with code 2 naming the first stage
  that is not implemented yet; for now the grid's two fields and
  `world.json` (M2-01).
- `kuniumi_sim::generate` runs through the stage runner; `Settings` has the
  new flag `allow_partial` (default `false`, which keeps the old behaviour),
  `Field` the new member `kind`, and `World` the members `settings` and
  `stages`; `.kfield` files are unchanged and are read back as quantities
  (M2-01).
- `kuniumi_sim::generate` and `kuniumi generate` run the Tectonics stage
  after the Grid stage and now stop at the Crust stage
  (``stage `crust` is not implemented yet``); the world they write has the
  three plate fields too (M2-02).
- The world written by `kuniumi generate` has the three fields of the
  continents and boundaries after the plate fields; the plate fields and
  their golden values are unchanged (M2-03).
- The one-shot crust reaches Earth's land fraction (M2-11). Ocean-floor
  age: the first arrival of fronts that grow from every ridge cell at its
  own half spreading rate (`ridge_times`, which removes most of the seams
  where neighbouring ridge cells spread at very different rates; steps
  remain next to very slow ridge cells, M2-13), ranked onto Earth's
  triangular area–age distribution (`earth_area_age`, Parsons 1982: at
  most 180 Myr, mean 60 Myr; plates without a ridge stay at 180 Myr); the
  older side of an ocean–ocean subduction is decided by these times.
  Continental margins: crust thins linearly from 15 km at the coastline to
  35 km at 600 km inland (`margin_distances`, `MARGIN_COAST_THICKNESS_M`,
  `MARGIN_WIDTH_M`), from every coast on an ocean of at least 0.5 % of the
  surface (`MARGIN_MIN_OCEAN_FRACTION`) except at convergent fronts, and
  from continental rifts; continental crust is now 15–70 km thick. With
  the formulas of M2-05, without its deepened trenches, the land fraction
  for Earth's parameters is 27.9–29.7 % at 10,000 cells (seeds 1–10) and
  28.7–31.4 % at 650,000 cells, instead of 40.0 % (with the trenches:
  M2-05).
- `kuniumi_sim::generate` and `kuniumi generate` run the Crust stage after
  the Tectonics stage and now stop at the Climate stage
  (``stage `climate` is not implemented yet``); the world they write has
  the field `surface_height`, and `world.json` the solved sea level instead
  of 0 m. The earlier fields and their golden values are unchanged (M2-05).
- World directories (`world_dir::write`, `kuniumi generate`) of a world
  with `surface_height` hold `metrics.txt` (the metrics of
  `kuniumi_sim::metrics` with Earth's values, the axes of the chart and the
  histogram) and `hypsometry.png` (the histogram as a 1040 × 540 bar
  chart), written before `world.json`; existing ones are removed before
  the first field is written. A non-finite `surface_height` is rejected
  before anything is written. `kuniumi-io` now depends on `kuniumi-sim`
  (M2-06).
- `kuniumi params` prints the settings file of Earth with the default
  settings (TOML, SI keys) in place of its aligned listing (M2-07).
- `PlanetParams::validate` bounds the planet to the Earth-like range of
  D-12 instead of accepting any positive finite value (any non-negative
  ocean volume): radius 1,500 to 13,000 km, surface gravity 1 to 40 m/s²,
  rotation period 8 to 240 h, orbital period 90 to 3,650 days, stellar flux
  400 to 2,500 W/m², ocean volume 0 to 10¹¹ km³, surface pressure 0.1 to
  10 bar (the `MIN_` and `MAX_` constants of `PlanetParams`), and the bulk
  density of the radius and gravity (`PlanetParams::bulk_density_kg_m3`,
  `3g / (4πGR)`) 3,000 to 16,000 kg/m³, reported as the new
  `ParamError::Inconsistent` and, from settings files, the new
  `ConfigError::Inconsistent` naming both keys as written. Parameters
  outside these bounds that were accepted before, such as Mars's surface
  pressure (0.006 bar) or a 13,000 km planet with Earth's gravity, are now
  rejected; worlds within them are unchanged. A compile-time assertion
  keeps every planet within the bounds slower than break-up. New
  `kuniumi_math::constants::GRAVITATIONAL_CONSTANT_M3_KG_S2` (CODATA 2018)
  (M2-08).
- `continents` fills the oceanic pockets that its growth leaves inside the
  continents: every enclosed oceanic region smaller than 0.1 % of the
  surface (`MIN_ENCLOSED_OCEAN_FRACTION`), except the largest ocean, becomes
  continental, and as much area is trimmed from the coast that grew last
  (cells on a coast whose continental neighbours on their plate form one
  run), so the continental fraction stays within one cell above the
  target, each plate's continent stays connected and no two kept seas are
  joined (a trimmed cell's oceanic neighbours lie in one region). Mostly
  holes of one to a few cells, at ocean depth in M2-05: 12–42 cells at
  10,000 cells and 43–205 at 100,000 for Earth's 40 %. Larger enclosed seas
  are kept. The ocean regions of `margin_distances` come from the same
  helper, unchanged (M2-10).
- The ocean floor of a plate without a divergent cell is no longer
  180 Myr old everywhere: it is ranked after all floor that a ridge
  reaches, by its distance from its plate's boundary (the new
  `ridgeless_distances`; `earth_area_age` takes these distances as a fourth
  argument), so it takes the oldest part of Earth's area–age distribution,
  older away from the boundary; the trench rule is unchanged. No ocean
  floor is at 180 Myr any more and the mean age is 60 Myr in every world:
  at 650,000 cells (seeds 1–10) such plates hold 10.3–14.7 % of the ocean
  floor in seeds 3, 6, 9 and 10, now 111–122 Myr old at their boundary and
  close to 180 Myr inside; in these four worlds the mean ocean-floor age
  falls from 72–78 Myr to 60.0 Myr, the sea level rises from between −819
  and −892 m to between −689 and −726 m, and the land fraction falls from
  29.4–31.4 % to 28.5–30.8 % (M2-12).
- The fronts of `ridge_times` grow at least at the new
  `MIN_SPREADING_RATE_M_S`, 1 cm/yr of full spreading (0.5 cm/yr half
  rate), about the rate of Earth's ultraslow ridges: a ridge cell that
  spreads more slowly, as near the relative pole of its plates, no longer
  leaves a step of up to about 100 Myr in age (kilometres in the sea floor)
  next to it. At 100,000 cells, seed 1, the pairs of neighbouring oceanic
  cells of one plate more than 20 Myr apart fall from 81 to 0 and the
  largest step next to a ridge cell from 105.9 to 11.5 Myr; at 650,000
  cells (seeds 1–10) from 0–110 pairs to none, and from 1.7–57.1 Myr to
  1.6–4.9 Myr. The convergence rates and boundary types are unchanged; the
  land fraction changes by less than 0.01 points and the sea level by less
  than 1 m (M2-13).
- `surface_height.png` of a world directory shows land and sea: drawn with
  the elevation map at the world's sea level (blue below it, green to
  brown above it, a hard break at the coast) instead of viridis, over the
  x axis of `hypsometry.png` shifted by the sea level (−11 km to +9 km
  above it for Earth-like worlds; the field's own range only if a bound
  overflows), so the chart is the map's legend and an elevation has the
  same colour in every such world. `world_dir::colormap(field)` becomes
  `world_dir::map_colours(field, world, metrics)`, which returns the range
  too; the maps of the other fields, the fields themselves, `world.json`
  and `metrics.txt` are unchanged (M2-14).

### Removed

- The hand-written argument parser and usage text of the CLI (M1-09).

### Fixed

- The two rustdoc warnings of the workspace: the public documentation of
  `kuniumi_io::colormap` no longer links to the private `VIRIDIS`, and that
  of `kuniumi_sphere::points` no longer spells out the target of its `TAU`
  link, which resolves to `core::f64::consts::TAU` as before (M2-09).

### World-breaking

- M2-13: `crust_age` changes in every world where a ridge cell spreading
  more slowly than 1 cm/yr is the first to reach one of its neighbours
  (seeds 1–8 at 650,000 cells; not seed 9, which has no such ridge cell,
  nor seed 10, where faster fronts reach every neighbour of its slow ridge
  cells first; the golden world, 9 of whose 413 ridge cells are below the
  floor), and with it `surface_height`, the sea level and, in some worlds,
  `trench` (the side that sinks where ocean meets ocean follows the times:
  for example 112, 72 and 66 cells in seeds 3, 5 and 6 at 650,000 cells,
  12, 15 and 2 cells in seeds 1, 3 and 5 at 100,000 cells; none in the
  golden world); the plates, continents, boundaries, thickness and uplift
  are unchanged.
- M2-12: `crust_age` changes in every world with a plate without a
  divergent cell (4 of seeds 1–10 at 650,000 cells, 4 at 10,000 cells, the
  golden world among them), and with it `surface_height` and the sea
  level; the ages of the other worlds, the plates, continents, boundaries,
  thickness, uplift and `trench` are unchanged.
- M2-11: `crust_thickness`, `crust_age` and `trench` change for every seed
  (thinned margins, the time-nearest ridge and Earth's area–age
  distribution); the plates, continents, boundaries and `uplift_rate` are
  unchanged.
- M2-10: `crust_type` changes for nearly every seed (the oceanic pockets
  filled, the coast that grew last trimmed), and with it `crust_thickness`,
  `crust_age`, `uplift_rate`, `surface_height` and the sea level, and
  `trench` for some seeds (12 of seeds 1–20 at 10,000 cells, up to 13
  cells); the plates, poles and boundaries are unchanged. For Earth's
  water the land fraction falls by up to 0.2 percentage points and the sea
  level moves by a few metres.

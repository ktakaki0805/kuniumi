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

### Removed

- The hand-written argument parser and usage text of the CLI (M1-09).

### World-breaking

- None yet.

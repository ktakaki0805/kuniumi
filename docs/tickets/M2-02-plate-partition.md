# M2-02 Plate partition and Euler poles

- Status: done
- Milestone: M2
- Depends on: M2-01
- Crates: kuniumi-sim (stages/tectonics, pipeline), kuniumi-world, kuniumi-io (world.json settings), kuniumi-cli (tests; one line of main.rs)
- Design: プレートテクトニクス / 表現; プレートテクトニクス / プレート速度と分裂
- Decisions: none

## Goal

The surface is divided into irregular, contiguous plates, each rotating about
its own Euler pole at Earth-like speeds. This is the one-shot tectonics that
M6 later replaces with time evolution.

## Scope

In:
- `Settings::plate_count` (default 12; valid 2–64).
- Partition by seeded multi-source growth: choose seed cells, then grow plates
  from a priority queue whose priorities carry random noise per cell (stream
  `"tectonics/partition"`), which gives irregular boundaries. The heap is keyed
  by `(priority via total_cmp, cell index)`, so ties are deterministic.
- Euler pole: a random unit vector; angular speed chosen so that the plate's
  maximum surface speed lies between 1 and 10 cm/yr (3.2 × 10⁻¹⁰ to
  3.2 × 10⁻⁹ m/s).
- Outputs: `plate_id` (class), `plate_velocity_east` and
  `plate_velocity_north` (m/s), from `ω × r` scaled by the radius.

Out:
- Continents and boundaries (M2-03).

## Implementation notes

Gaps found while refining (the ticket was `ready`). The partition method was
put to the owner (below); the rest are implementation choices inside the
ticket's crates.

### Question for the owner: the partition method

The method as written (seed cells drawn uniformly, then growth with noise per
cell only) cannot meet the size criterion (largest/smallest plate area in
[3, 50] for seeds 1–20). Measured with a throwaway prototype on the real grid,
12 plates, `SphereGrid::fibonacci(cells, seed)` for seeds 1–20:

| Method | 10,000 cells | 100,000 cells | Boundaries |
| --- | --- | --- | --- |
| Literal, moderate noise (slowness `1 + u`) | 2.2–32.4 | 2.5–19.7 | Voronoi polygons, straight edges, cell-scale roughness |
| Literal, strong noise (`exp(5(u − ½))`), seed spacing 0.3 | 3.0–43.8 | 2.1–9.8 | fractal, coastline-like; depends strongly on resolution |
| A: per-plate growth rates in [1, 3), seed spacing 0.4, noise `exp(3(u − ½))` | 5.3–41.4 (median 10) | 4.0–30.1 (median 9) | curved, irregular; a few large and some small plates |

(`u` uniform in [0, 1) per cell; seed spacing as a fraction of
`√(4π / plate_count)` rad.) Earth's 12 largest plates span a ratio of about
20 (Pacific about 2.6 sr, the twelfth, Amur, about 0.13 sr; Bird 2003, "An
updated digital model of plate boundaries", G³ 4(3), 1027).

Recommended: A, the plan below. Alternatives: keep the literal method and
relax the lower bound of the criterion to 2; or the literal method with strong
noise, accepting fractal boundaries and a borderline criterion.

### Decision (owner, 2026-11-09)

The owner chose A: growth rates per plate uniform in [1, 3), seed spacing at
least `0.4 √(4π/K)` rad, slowness per cell `exp(3(u − ½))`, the fastest
speed of each plate log-uniform in 1–10 cm/yr, and the speed criterion read
as "the fastest cell of every plate moves at 1–10 cm/yr and no cell
faster". This approves the plan below.

### `kuniumi-world`

```rust
pub struct Settings { pub cells: u32, pub plate_count: u32, pub allow_partial: bool }
impl Settings {
    pub const DEFAULT_PLATE_COUNT: u32 = 12;
    pub const MIN_PLATE_COUNT: u32 = 2;
    pub const MAX_PLATE_COUNT: u32 = 64;
}
```

`validate` checks `cells`, then `plate_count`
(`ParamError::OutOfRange { name: "plate_count", .. }`).

### `kuniumi-io`

`SettingsJson` gains `plate_count` (between `cells` and `allow_partial`, the
order of `Settings`), so `world.json` records it. `FORMAT_VERSION` stays 1:
nothing reads `world.json` yet (MV-01 is the first reader), and the settings
grow with the stages.

### `kuniumi-sim`: `stages::tectonics`

The Tectonics stage becomes a `StageFn` in `pipeline::implementation`; the
pipeline now stops at `crust`.

```rust
pub const PARTITION_STREAM: &str = "tectonics/partition"; // seeds, rates, noise
pub const POLES_STREAM: &str = "tectonics/poles";         // axes, speeds
pub const PLATE_ID: &str = "plate_id";                    // class
pub const PLATE_VELOCITY_EAST: &str = "plate_velocity_east";   // m/s
pub const PLATE_VELOCITY_NORTH: &str = "plate_velocity_north"; // m/s
pub const MIN_PLATE_SPEED_M_S: f64 = 0.01 / SECONDS_PER_YEAR;  // 1 cm/yr, 3.169e-10
pub const MAX_PLATE_SPEED_M_S: f64 = 0.10 / SECONDS_PER_YEAR;  // 10 cm/yr, 3.169e-9

pub fn partition(grid: &SphereGrid, plate_count: u32, seed: u64) -> Vec<u32>;
pub struct EulerPole { pub axis: Vec3, pub angular_speed_rad_s: f64 }
impl EulerPole { pub fn velocity_m_s(&self, position: Vec3, radius_m: f64) -> Vec3; }
pub fn euler_poles(grid: &SphereGrid, plates: &[u32], plate_count: u32,
                   radius_m: f64, seed: u64) -> Vec<EulerPole>;
pub fn run(world: &mut World, progress: &mut dyn FnMut(f64)) -> Result<(), SimError>;
```

Partition (`K = plate_count`, `n` cells; all draws from `PARTITION_STREAM`
in this order):

1. Seed cells: draw `rng.below(n)`; accept a cell that is not a seed yet and
   lies at least `d_min = 0.4 √(4π/K)` rad (great circle) from every seed
   accepted so far; plate `k` is the `k`-th accepted seed. After 10,000
   rejections in a row the spacing is dropped (a fallback for tiny grids, so
   the loop always ends).
2. Growth rates `r_k = 1 + 2 u_k` in plate order (in [1, 3)).
3. Slowness per cell `w_i = exp(3 (u_i − ½))` in cell order (in
   [e^−1.5, e^1.5), a factor of 20; `det::exp`).
4. Multi-source Dijkstra: each seed enters the heap with priority 0. Pop the
   smallest entry; skip it if its cell has a plate; otherwise give the cell
   that plate and push every neighbour `j` without a plate with priority
   `p + d_ij (w_i + w_j) / 2 / r_k` (`d_ij` the centre distance). The heap is
   a `BinaryHeap` of `Reverse((priority via total_cmp, cell, plate))`: with the
   plate as a last key no two different entries compare equal, so the order of
   pops is the sorted order whatever the heap's internals.

Every cell is reached (the grid is connected) and is given the plate of the
neighbour that pushed it, so every plate is connected.

Euler poles (from `POLES_STREAM`, plate by plate, in this order): the
target maximum speed, log-uniform, `v_k = v_min exp(u ln(v_max / v_min))`,
clamped to `[v_min, v_max]` (Earth's plate speeds spread over a decade with
more slow plates than fast ones); then the axis, uniform on the sphere
(`z = 2u − 1`, `λ = 2πv`, `det::sin_cos`); the lever
`s_k = max_i |axis × p_i|` over the plate's cells (the sine of the largest
angle from the axis); while `s_k < 1e-6` (only possible for a plate of one
cell next to its axis) the axis is drawn again; `ω_k = v_k / (R s_k)`. So
the fastest cell of each plate moves at `v_k` and no cell faster.

Velocity of cell `i` on plate `k`: `v = (ω_k axis_k) × (R p_i)` (m/s),
`east = v · coords::east(p_i)`, `north = v · coords::north(p_i)`. The stage
adds `plate_id`, `plate_velocity_east`, `plate_velocity_north`, in that order,
and reports no progress of its own.

### Tests

- `stages::tectonics::tests`: the streams are under `"tectonics/"`;
  hand-computed velocities (axis +Z on the equator and at 60° N, axis +X at
  the north pole); for 20 seeds at 10,000 cells and for 1,000 cells with 2,
  12 and 64 plates: every cell has one plate id in `0..K`, every plate is
  non-empty and connected (breadth-first search inside the plate from one of
  its cells reaches all of them); seeds keep their spacing; the area ratio
  of seeds 1–20 (printed); speeds for `R` = 6,371 km and 3,185.5 km; the
  velocity of a cell equals the rate at which `Vec3::rotate_about` about its
  pole moves it (centred difference, right-hand rule); same seed gives the
  same bits, another seed differs.
- `Settings`: `plate_count` bounds (2 and 64 valid, 0, 1, 65 rejected by name).
- `pipeline`, `lib.rs`: the real table now has the Tectonics stage; `generate`
  stops at `crust`.
- `world_dir`: `plate_count` in `world.json`, the layouts updated.
- CLI: `generate` prints `[tectonics] 0%` / `100%`, writes the three new
  fields and stops at `crust`; golden fingerprints of `plate_id` and both
  velocity fields for 10,000 cells, seed 1, recorded after the invariant
  tests and the visual check.

## Acceptance criteria

- [x] Every cell belongs to exactly one plate and every plate is connected
      (graph search test), for seeds 1–20 at 10,000 cells and for 1,000 cells
      with 2, 12 and 64 plates.
- [x] Plate sizes are irregular: the ratio of the largest to the smallest plate
      area is between 3 and 50 for seeds 1–20 at 10,000 cells (distribution
      reported, and also measured at 100,000 cells).
- [x] For planets of Earth's radius and of half Earth's radius, the fastest
      cell of every plate moves at 1 to 10 cm/yr and no cell faster than
      10 cm/yr (within 1e-12 relative).
- [x] Golden fingerprint of `plate_id` for 10,000 cells, seed 1.
- [x] `out/M2-02/plate_id.png` viewed and described.
- [x] `mise run ci` passes.

## Verification

```sh
cargo run -p kuniumi-cli -- generate --seed 1 --cells 100000 --out out/M2-02
```

## Result

What was done:
- `kuniumi-world`: `Settings::plate_count` with `DEFAULT_PLATE_COUNT` (12),
  `MIN_PLATE_COUNT` (2) and `MAX_PLATE_COUNT` (64); `validate` checks it
  after `cells` (`OutOfRange { name: "plate_count" }`).
- `kuniumi-sim`: `stages::tectonics` (the stage `run`, `partition`,
  `EulerPole`, `euler_poles`, the field names, the streams
  `tectonics/partition` and `tectonics/poles`, the speed bounds and the
  partition constants `SEED_SPACING`, `MAX_GROWTH_RATE`, `SLOWNESS_SPREAD`),
  tests in `stages/tectonics/tests.rs`; `pipeline::implementation` returns
  it for `StageId::Tectonics`, so `generate` now stops at `crust`.
- `kuniumi-io`: `SettingsJson::plate_count`; the module example and the
  layout tests of `world.json`.
- `kuniumi-cli`: tests only (stderr, the files written, `world.json`, golden
  values); the binary itself is unchanged apart from `..Settings::default()`.
- `CHANGELOG.md` (Added, Changed).

Deviations from the plan: none in behaviour. The fallback of the seed
choice is documented on `partition` itself (MAX_REJECTIONS = 10,000 in a
row, as planned); the per-plate draws of the poles are the speed first, then
the axis (the notes above were corrected to this order).

Evidence:

1. Partition and connectivity:
   `every_cell_has_one_plate_and_plates_are_connected` checks, for seeds
   1–20 at 10,000 cells (12 plates), for 1,000 cells with 2, 12 and 64
   plates (seeds 1–3) and for the octahedron with 1 to 6 plates, that every
   cell has one id in `0..K`, every plate has a cell, and a breadth-first
   search through neighbours of the same plate from the plate's first cell
   reaches all its cells. `seed_cells_keep_their_spacing` (10,000 and 1,000
   cells, 2 to 64 plates, five seeds each: every pair of seeds at least
   `0.4 √(4π/K)` rad apart); `seed_cells_end_on_a_full_grid` (as many
   plates as cells on the octahedron and on 1,000 cells, whose closest
   neighbours are nearer than the spacing, so only the fallback can finish:
   every cell becomes one seed and one plate);
   `the_same_seed_gives_the_same_plates` (bits of plates and poles equal for
   the same seed, different for another).
2. Plate sizes (`plate_sizes_are_irregular`, area-weighted, 12 plates),
   largest / smallest area for seeds 1–20 at 10,000 cells, sorted:
   7.7, 8.8, 9.3, 9.9, 10.0, 10.1, 10.1, 10.9, 11.0, 11.8, 12.7, 13.2, 14.3,
   15.7, 16.0, 17.3, 17.7, 19.3, 19.4, 23.1 (all in [3, 50]; median 11.8).
   Measured at 100,000 cells with the same functions (throwaway program in
   the scratchpad): 3.8, 4.9, 4.9, 5.1, 6.9, 7.1, 7.6, 8.1, 9.3, 10.3, 12.1,
   13.1, 13.7, 14.5, 15.4, 16.5, 16.5, 17.8, 25.2, 34.7 (median 11). Seed 1
   at 100,000 cells: areas 2.31, 1.62, 1.57, 1.27, 1.16, 1.15, 0.94, 0.92,
   0.73, 0.44, 0.29, 0.16 sr (Earth: Pacific about 2.6 sr, the twelfth
   plate about 0.13 sr, Bird 2003).
3. Speeds (`plate_speeds_stay_within_1_to_10_cm_per_year`): for seeds 1–20
   at 10,000 cells and radii 6,371 km and 3,185.5 km, the fastest cell of
   each of the 12 plates, from the east and north components, lies in
   [1, 10] cm/yr within 1e-12 relative, so no cell is faster than 10 cm/yr.
   Over the 240 plates (Earth's radius) the fastest speeds run from 1.013 to
   9.862 cm/yr, median 3.29 cm/yr; 114 of 240 lie below √10 cm/yr, the
   median of the log-uniform distribution (`plate_speeds_are_log_uniform`
   accepts 35 % to 65 %). Seed 1 at 100,000 cells: angular speeds 0.10 to
   0.77 °/Myr (Earth's large plates: about 0.1 to 1 °/Myr).
   `the_speed_bounds_are_1_and_10_cm_per_year` (3.1688e-10 and
   3.1688e-9 m/s for a Julian year); `velocities_of_hand_computed_rotations`
   (axis +Z at 0° N 0° E: east `ωR`; at 60° N: east `ωR/2`; axis +X at the
   north pole: east `−ωR`; at 0° N 90° E: north `ωR`);
   `velocities_follow_the_right_hand_rule` (every seventh cell: the velocity
   equals the centred difference of `Vec3::rotate_about` by ±1e-6 rad about
   its pole, within 1e-8 of `ωR`); `a_plate_of_one_cell_gets_a_finite_pole`
   (six one-cell plates on the octahedron, 50 seeds: unit axes, lever at
   least 1e-6, finite `ω`, speed in range);
   `an_axis_through_a_plate_of_one_cell_is_drawn_again` (seeds 1–3: a grid
   of 100 Fibonacci points whose cell 0 is moved exactly onto the first
   axis drawn for plate 0, cell 0 alone as plate 0: the pole returned has
   another axis, a lever of at least 1e-6, a finite `ω` and a speed in
   range; it fails when the redraw is removed);
   `euler_poles_rejects_a_partition_of_another_grid` (a plates array of the
   wrong length panics); `the_stage_adds_its_fields`
   (names, units, kinds and order of the three fields, values equal to
   `partition`, `euler_poles` and the velocities, all finite, no progress
   of its own).
4. Golden fingerprints (`tests/cli.rs::generate_golden_fingerprints`, 10,000
   cells, seed 1, through `kuniumi generate` and
   `kuniumi_io::field_file::fingerprint`): `plate_id`
   `0x5E9F_DB0A_B6B7_D521`, `plate_velocity_east` `0x8FFF_DF16_93AF_9B25`,
   `plate_velocity_north` `0x1A26_5DF6_79D2_2BB9`, recorded after the tests
   above passed and the maps below were looked at; `cell_area` and
   `neighbor_count` unchanged. `generate_writes_the_world_so_far_then_exits_with_code_2`
   (1,000 cells: stdout and stderr byte for byte, now with
   `[tectonics] 0%`, `[tectonics] 100%` and ``stage `crust` is not
   implemented yet``; 2048 × 1024 maps of the five fields; `plate_count` 12
   and the stages `grid`, `tectonics` in `world.json`; 12 non-empty plates;
   no speed above 10 cm/yr); `generate_writes_the_same_bytes_twice` (the 11
   files byte-identical).
5. Maps (`kuniumi generate --seed 1 --cells 100000 --out out/M2-02`,
   release, exit code 2 naming `crust`), viewed:
   - `plate_id.png`: 12 plates in the 12 colours of the palette; a north
     polar cap (orange) reaching down to about 40° to 60° N, a large plate in
     the south centre (green) and a band across the south (purple) are the
     largest; a dark blue oval near 0° N 20° W between red and green, a
     brown plate in the north-east and a thin south polar cap (lavender) are
     the smallest. Boundaries are curved arcs, rough at the scale of a few
     cells, without the straight edges of a Voronoi partition; slower plates
     sit in bays of faster ones. Pink and purple continue across the
     antimeridian; no plate has a detached piece.
   - `plate_velocity_east.png` and `plate_velocity_north.png` (viridis over
     each field's range): smooth gradients inside every plate, as a rigid
     rotation gives, and jumps exactly along the outlines of `plate_id.png`;
     over the south polar cap both components swing with longitude, as the
     east/north frame turns about the pole.
   - `out/M2-02/10k/plate_id.png` (the golden case, 10,000 cells, seed 1):
     12 contiguous plates, polar caps light green and lavender, small plates
     light blue, red and orange; steps of single cells (about 2°) along the
     boundaries; yellow and dark blue continue across the antimeridian.
   - `out/M2-02/650k/plate_id.png` (650,000 cells, seed 1): 12 plates with
     curved boundaries; relative to the plates' size they are smoother than
     at 100,000 cells, since the slowness varies per cell.
6. `mise run ci` passes: fmt, clippy `-D warnings`, the tests in debug and
   release (336 in debug, 327 in release, which skips the debug-only
   tests), MSRV 1.85.0, cargo-deny (advisories, bans, licences, sources
   ok; warnings only for allowed licences not in use). `cargo doc` gives no
   warnings for `kuniumi-sim`.

Timing (release, this machine): `kuniumi generate --seed 1 --cells 650000`
takes 3.9 s and 175 MB (the same as after M2-01: the grid dominates); the
partition takes 0.33 s and the poles 0.004 s. At 100,000 cells `generate`
takes 0.76 s and 43 MB.

Review:
- Round 1 (APPROVE): the reviewer re-ran `mise run ci` and the maps
  (byte-identical), wrote its own partition and poles from this ticket's
  notes (bit-identical for 10,000 cells × seeds 1–20, 1,000 cells with 2
  and 64 plates, 30,000 cells with 33 plates, radii from 100 km to
  6,371 km), checked connectivity with union-find up to 650,000 cells with
  64 plates, and recomputed both size distributions. R1-1 minor, no test
  reached the redraw of an axis through a one-cell plate (fixed:
  `an_axis_through_a_plate_of_one_cell_is_drawn_again`, built on the
  reviewer's construction; it fails when the redraw is removed). R1-2 nit,
  the log-uniform speeds had no source (fixed: stated as the owner's
  modelling choice over the design's range). R1-3 nit, a test missing from
  the module's list (fixed). R1-4 nit, `euler_poles` did not check the
  length of `plates` (fixed: assertion, `# Panics`, test). R1-5 nit, two
  numbers of this Result (fixed: 0.94 sr; tests counted per run).
- Round 2 (APPROVE): R1-1 to R1-5 confirmed resolved (the redraw mutation
  fails only the new test); the reviewer also mutated the growth (dropping
  the skip of assigned cells), the rates and spacing, the lever (minimum
  instead of maximum), the sign of `ω × r` and the speed distribution, each
  caught by a test, and recomputed the golden fingerprints from its own run.
  R2-1 minor, radii at the extremes of `PlanetParams::validate` (5e-324 m,
  `f64::MAX`) give non-finite or inaccurate velocities; not reachable from
  the CLI and caught by the runner's debug assertion (recorded as the draft
  M2-08, as the reviewer suggested, since the bounds of the planet
  parameters are the owner's decision). R2-2 nit, the `Crates:` line and the
  median speed (fixed: one line of `main.rs` named; median 3.29 cm/yr, the
  mean of the 120th and 121st values). R2-3 nit, a private intra-doc link in
  `kuniumi-io/src/colormap.rs` that predates this ticket (recorded as the
  draft M2-09).

Follow-ups: M2-08 (draft, physical bounds for the planet parameters), M2-09
(draft, rustdoc warnings and whether `mise run ci` checks them).

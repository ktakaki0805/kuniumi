# M3-03 Implicit stream-power erosion with uplift

- Status: done
- Milestone: M3
- Depends on: M3-02
- Crates: kuniumi-world (`ErosionSettings`), kuniumi-sim (stages/erosion)
- Design: 侵食と地形発達 / 基本方程式; 侵食と地形発達 / 1ステップの処理; 侵食と地形発達 / 調整の目安
- Decisions: none

## Goal

River incision balanced against tectonic uplift, solved implicitly so large time
steps stay stable: the core of realistic relief.

## Scope

In:
- `∂h/∂t = U − K Q^m S` (n = 1), with Q the discharge (m³/s) from the drainage
  graph, solved per step in stack order (Braun & Willett 2013):
  `h_i ← (h_i + Δt·U_i + F·h_r) / (1 + F)`, `F = Δt·K·Q_iᵐ / L_ir`.
- Units: with E in m/s and Q in m³/s, K is in m^(1−3m)·s^(m−1). Document the
  default K, m = 0.45, Δt and how they were chosen.
- Ocean cells are a fixed base level (not eroded here).
- `ErosionSettings` (K, m, Δt in s, step count) in `kuniumi-world`, with its
  defaults and `validate`, taken by the solver and the driver as an argument;
  and a uniform precipitation placeholder of 1 m/yr (in m/s) until the
  climate stage exists.
- A step driver: rebuild drainage, route depressions, solve, repeat; progress
  reported per step.

Out:
- Diffusion (M3-04); deposition and isostatic rebound (M7).
- `Settings::erosion` and its records in `world.json` and settings files
  (M3-07, split off this ticket by the owner).
- Wiring the Erosion stage into the pipeline, its world fields and its golden
  fingerprint in the CLI tests (M3-05).

## Implementation notes

Gaps found while refining (the ticket was `ready`): where `ErosionSettings`
lives; the API; what the base level of a river is where it reaches the sea
or a lake; what happens in the closed depressions that M3-02 routes across;
how the rounding of a nearly converged cell is kept from making a flat; the
defaults; the validity of the settings; how the criteria are measured; and
how the images are made.

### Split (owner, 2026-11-15)

`ErosionSettings` as a field of `Settings` is not reached by `kuniumi-sim`
alone: `world.json` records every field of `Settings` (`kuniumi-io`, whose
`SettingsJson` denies unknown fields, so a new key means a new format
version), and settings files write and read every field (`kuniumi-config`).
With the solver that would have been about 1,800 changed lines in five
crates, the settings plumbing independent of the solver. The owner chose to
split: this ticket defines `ErosionSettings` in `kuniumi-world` with its
defaults and `validate`, and the solver and the driver take it as an
argument; M3-07 (draft, after M3-03, before M3-05) adds `Settings::erosion`
and records it in `world.json` and in settings files.

### API

`kuniumi-world` (`settings.rs`, re-exported from the crate root):

```rust
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ErosionSettings {
    pub erodibility: f64,        // K, m^(1−3m) s^(m−1)
    pub discharge_exponent: f64, // m
    pub time_step_s: f64,        // Δt (s)
    pub steps: u32,
}
impl ErosionSettings {
    pub const DEFAULT_ERODIBILITY: f64 = 7.5e-10;
    pub const DEFAULT_DISCHARGE_EXPONENT: f64 = 0.45;
    pub const DEFAULT_TIME_STEP_S: f64 = 1.0e4 * SECONDS_PER_YEAR; // 3.15576e11 s
    pub const DEFAULT_STEPS: u32 = 1_000;
    pub fn validate(&self) -> Result<(), ParamError>;
}
impl Default for ErosionSettings { /* the constants */ }
```

The field `erodibility` has no unit suffix because its unit,
m^(1−3m) s^(m−1), depends on `discharge_exponent`; the doc comment states
it (docs/conventions.md, "Units").

`kuniumi-sim`, module `stages::erosion` (tests in `stages/erosion/tests.rs`,
so `cargo test -p kuniumi-sim erosion` selects them):

```rust
/// 1 m/yr in m/s (1 / SECONDS_PER_YEAR), until the Climate stage exists.
pub const PLACEHOLDER_PRECIPITATION_M_S: f64;

/// The network one incision step follows.
pub struct Network<'a> {
    pub receivers: &'a [u32],            // Drainage::receivers
    pub stack: &'a [u32],                // Drainage::stack
    pub receiver_distance_m: &'a [f64],  // L_i (m), > 0 except at base levels
    pub discharge_m3_s: &'a [f64],       // Q_i (m³/s), ≥ 0
}
/// One implicit step on `height_m`, in place.
pub fn incise(network: &Network<'_>, uplift_rate_m_s: &[f64], sea_level_m: f64,
              settings: &ErosionSettings, height_m: &mut [f64]);

/// L_i = d_{i r(i)} R (m), from SphereGrid::centre_distances_rad; 0 at base levels.
pub fn receiver_distances_m(grid: &SphereGrid, drainage: &Drainage, radius_m: f64) -> Vec<f64>;

/// What erosion depends on besides the heights.
pub struct ErosionInputs<'a> {
    pub grid: &'a SphereGrid,
    pub radius_m: f64,
    pub sea_level_m: f64,
    pub uplift_rate_m_s: &'a [f64],
    pub precipitation_m_s: &'a [f64],
}
/// `settings.steps` steps; progress (k + 1) / steps after step k.
pub fn erode(inputs: &ErosionInputs<'_>, settings: &ErosionSettings,
             height_m: &mut [f64], progress: &mut dyn FnMut(f64));
```

`incise` takes slices rather than a `Drainage` so that the synthetic 1D
channel of the first criterion needs no grid. Panics, documented: lengths
that differ from the number of cells, a stack that is not every cell once,
a non-finite height, uplift or sea level, a discharge that is negative or
not finite, a receiver distance that is not finite and positive at a cell
that is not a base level, a precipitation that is negative or not finite,
invalid settings (`ErosionSettings::validate`). The private `next_up` of
`drainage::depressions` becomes `pub(crate)` for the solver.

### Algorithm

One step of `erode` (design, 侵食と地形発達 / 1ステップの処理, steps 1–4):

1. `h_fill = fill_depressions(grid, h, sea_level)` (M3-02).
2. `drainage = Drainage::build(grid, &h_fill, sea_level, R)` (M3-01): routing
   on the filled surface, erosion on the original one.
3. `Q = drainage.accumulate(P_i a_i R²)` (m³/s; the per-cell input is
   computed once per call of `erode`).
4. `L = receiver_distances_m(grid, &drainage, R)`.
5. `incise` on `h`.
6. `progress((k + 1) / steps)`.

`incise` (Braun & Willett 2013, "A very efficient O(n), implicit and
parallel method to solve the stream power equation governing fluvial
incision and landscape evolution", *Geomorphology* 180–181, 170–179, the
implicit update for n = 1), along the stack, each cell after its receiver,
with `W` the water surface after the step (m):

- An ocean cell (`h_i ≤ sea_level` at the start of the step) is a fixed base
  level: unchanged, uplift included; `W_i = sea_level`.
- `a = h_i + Δt U_i`.
- A land base level (a pit, only on land without ocean after the fill):
  `h_i = a`, `W_i = a`.
- Otherwise, with `r` the receiver and `b = W_r`, the water surface the
  river reaches: if `a ≤ b` the cell lies under the water of its receiver (a
  lake of a closed depression) and is not eroded: `h_i = a`, `W_i = b`. Else
  `F = Δt K Q_iᵐ / L_i` (`Q_iᵐ` by `det::pow`; `F = 0` if `Δt K` or `Q_iᵐ`
  is 0, so that no `0 · ∞` occurs) and
  `h_i = b + (a − b) / (1 + F)`, the ticket's update written relative to the
  base, which equals `(a + F b) / (1 + F)` in exact arithmetic and keeps the
  cell above `b` in floating point; if rounding gives `h_i ≤ b`, then
  `h_i = next_up(b)`; if it gives `h_i > a`, then `h_i = a`. `W_i = h_i`.

Choices in this, each documented in the code:

- Rivers grade to the water surface: a river reaching the sea erodes towards
  the sea level, not towards the sea floor of its ocean receiver (which lies
  kilometres deeper: grading to it would turn the coast into deep sea in one
  step), and a river entering a lake towards the lake's surface. This is the
  "fixed base level" of the ticket.
- Closed depressions (M3-02 routes the water across them on the filled
  surface): a cell at or below the water surface of its receiver is under a
  lake and is not eroded, only uplifted; the B&W update would otherwise
  raise it (deposition, which is M7) or cut the lake floor. A lake drains as
  its outlet is cut down.
- With these, every land cell that is eroded ends strictly above the water
  surface of its receiver and not above `a`: erosion never takes a cell
  below its receiver and never raises one. With `U = 0`, heights never rise
  and land stays land, so no land cell loses its last strictly lower
  neighbour (no new pit) and no closed depression appears; with uniform `U`
  no closed depression appears, while the rounding of `h + Δt U` could, in
  principle, merge two heights one ULP apart into a flat where the sum
  crosses a power of two (stated in the docs; the tests measure it). Uplift
  that differs between neighbours can dam a valley, which is tectonics,
  not the scheme.
- Uplift on ocean cells is ignored (they are fixed); the sea level is not
  changed, so the ocean below it keeps its volume.

The arithmetic is `+ − × ÷`, comparisons, `det::pow` and `next_up` on bit
patterns (docs/conventions.md, "Determinism", rule 2); the order is the
stack's. Time O(n) per step besides the fill (O(n log n)).

### Defaults and how they were chosen

- `m = 0.45`, the ticket's value, in the design's 0.4–0.5 (基本方程式).
- `K = 7.5e-10 m^(1−3m) s^(m−1)`: the design tunes K by results (調整の目安:
  major orogens 2–4 km high on average, the highest peaks 8–9 km). At
  steady state `S = U / (K Qᵐ)`; for a river crossing a belt on the default
  grid (650,000 cells, 28.0 km spacing `√(4π/n) R`, 7.85 × 10⁸ m² per cell),
  whose drainage area grows as Hack's law, `A_k = a k^(1/0.6)` at the k-th
  cell from the divide, with 1 m/yr of runoff, the divide stands at
  `Σ_k U L / (K (P A_k)ᵐ)`: 3.0–3.8 km for 5–8 cells (140–220 km) at the
  4 mm/yr of a collision front of M2-04 (0.1 × 4 cm/yr), and 7.5–9.4 km at
  the 10 mm/yr maximum uplift of M2-04. The same K written for drainage area
  instead of discharge (`E = K_A Aᵐ S`, `K = K_A / Pᵐ`) is
  `K_A = 1.0e-5 m^0.1 yr⁻¹`. Calibration against the metrics is M3-05 and M7.
- `Δt = 10⁴ yr` (3.15576 × 10¹¹ s), the lower end of the design's 1–10万年
  (基本方程式): on the default grid `F` runs from 0.036 for one cell's runoff
  to 1.5 for 10⁵ m³/s, so knickpoints move without much smearing; larger
  steps stay stable (the second criterion).
- `steps = 1,000`, the lower end of the design's 1,000–3,000 (侵食と地形発達),
  10 Myr in all, the few Myr per round of the climate–erosion alternation
  (気候との交互反復) over two or three rounds.
- `PLACEHOLDER_PRECIPITATION_M_S = 1 / SECONDS_PER_YEAR` (1 m/yr), the
  ticket's value.

`validate`: `erodibility` finite and positive, `discharge_exponent` finite
and non-negative, `time_step_s` finite and positive (`ParamError::OutOfRange`
naming the field, in field order); every `steps` is valid (0 leaves the
surface unchanged). These are the bounds within which the law is defined
and the solver stays finite; physical bounds for user input, if any, come
with the settings files (M3-07).

### Tests (`stages/erosion/tests.rs`; `ErosionSettings` in `kuniumi-world`)

- `one_step_by_hand`: a small network given as slices, `m = 0.5` and
  `Q ∈ {4, 16, 64} m³/s` (`Q^0.5` exact), `Δt K / L` chosen so that `F` is 1,
  4 or 8 and every division exact: an ocean cell unchanged (uplift
  ignored); a coastal cell grading to the sea level, not to the −4,000 m
  floor of its ocean receiver; a cell draining to a land cell; a lake cell
  below its receiver's water surface, only uplifted; a river entering the
  lake grading to its surface; a cell with `F` so large that the update
  rounds to its base (`next_up(b)`); a cell with `Q = 0` (not eroded); a
  pit (uplift only) and a cell draining into it (graded to the pit); a cell
  exactly at the sea level with uplift (ocean, fixed) and a cell draining
  into it (graded to the sea level). Values by hand, with the ticket's
  `(a + F b)/(1 + F)`.
- `extreme_settings_stay_finite`: `Δt K` underflowing to 0 while `Qᵐ`
  overflows, and the reverse: no `0 · ∞`, the cell not eroded; `F`
  overflowing to ∞, with `a − b` finite and with `a − b` overflowing too
  (no `∞ / ∞`): the cells end one ULP above their bases.
- `rounding_never_raises_an_eroded_cell`: `a = 1 + 2⁻⁵²`, `b = −2⁻⁵³`,
  `F = 1e-30`, where `b + (a − b)/(1 + F)` rounds to `1 + 2⁻⁵¹ > a`: the
  cell stays at `a`.
- `steady_state_of_a_1d_channel` (criterion 1): 50 cells, cell 0 an ocean
  outlet 100 m deep (sea level 0), `r(i) = i − 1`, `L = 1 km`,
  `Q_i = (50 − i)` m³/s (the cells upstream of `i`, 1 m³/s each),
  `U = 1 mm/yr`, default K and m,
  `Δt` 10⁴, 10⁵ and 10⁶ yr, from land at 1 m: steps until no height changes
  (at most a set cap); then every `(h_i − W_r) / L` within 1e-6 relative of
  `U / (K Q_iᵐ)`, the reference computed in the test with std `powf`.
- `large_time_steps_are_stable` (criterion 2): the Earth worlds of seeds
  1–3 at 10,000 cells, the dry Earth world of seed 1 (no water: the sea
  level is the lowest height, shared by the trench cells, so cells lie
  exactly at it), a smooth random surface (plane waves, sea level −150 m)
  and heights uniform in [−1,000, 1,000) m (sea level −900 m, full of
  closed depressions); `U = 0`, a uniform 1 mm/yr and, on the Earth
  worlds, their own `uplift_rate`; 10 steps each of `Δt` 10⁴, 10⁵ and
  10⁶ yr, the drainage of each step rebuilt in the test as `erode` builds
  it. After every step: all
  heights finite; ocean cells unchanged bit for bit; every land cell with
  `a > W_r` strictly above `W_r` and not above `a`, every other land cell
  exactly `a`; with `U = 0` or uniform: no land cell that had a strictly
  lower neighbour has none (no new pit; asserted), and no new closed
  depression (cells whose exact spill height, from Priority-Flood without
  ε written independently in the test, exceeds their height); land stays
  land; with the Earth uplift the pits it dams are counted and reported.
- `land_volume_decreases_without_uplift` (criterion 3): same surfaces,
  `U = 0`, 20 steps of the default `Δt` and of 10⁶ yr: the land volume
  `Σ_land (h − sea level) a R²` (index order) falls strictly at every step,
  no height rises, and the land cells stay the same.
- `steady_state_of_an_island`: a cap of land 0.6 rad in radius at 10 m on
  the 10,000-cell grid, ocean 1 km deep around it, uniform `U = 1 mm/yr`,
  steps of 10⁶ yr until no height changes by more than 1e-12 m, on planets
  of radius 1,500, 6,371 and 13,000 km: every land cell's slope to the water
  surface of its receiver within 1e-6 relative of `U / (K Qᵐ)`, with the
  distance from the great-circle distance of the centres (std `acos`) and
  `Q` from a walk down the receivers of every cell, both computed in the
  test; `receiver_distances_m` within 2e-12 rad × R of that distance at
  every cell (the driver end to end).
- `erode_reports_progress_per_step`: fractions `k/steps` for k = 1..steps;
  `steps = 0` reports nothing and leaves the heights unchanged bit for bit.
  `erode_steps_compose`: five steps at once equal five single steps, bit
  for bit.
- `receiver_distances_by_hand`: on the octahedron, `π/2 · R` for every
  non-base cell, 0 at the base levels.
- Panics: each condition listed under API (`the_unbroken_chain_runs`
  checks that the chain they break runs); those of `erode` on heights,
  uplift and the sea level with 0 steps
  (`erode_without_steps_accepts_the_unchanged_world` checks the unbroken
  case).
- `kuniumi-world`: the defaults, `validate` (each field rejected when NaN,
  infinite or out of its domain, in field order; 0 steps accepted).
- Golden fingerprint: FNV-1a 64 over the heights' bits after 10 steps of
  the default settings on the Earth world of seed 1 at 10,000 cells with
  its `uplift_rate` and the placeholder precipitation, recorded only after
  the tests above pass and the images are viewed.

### Images and timing

A scratch tool outside the repository (a Cargo package in the scratchpad
with path dependencies on the workspace crates, as in M3-01 and M3-02)
generates the Earth world of seed 1 at 100,000 cells through the Crust
stage, runs `erode` with the default settings, its `uplift_rate` and the
placeholder precipitation, and draws into `out/M3-03/`: `surface_height`
before and after (equirectangular, the elevation colours of the world
maps, M2-14, with hill shading), close-ups of a mountain belt and of a
continental interior (the centre of an orthographic view), and the
logarithm of the discharge after erosion; after 100 steps (1 Myr) as well
as after the default 1,000. It times
one step in release at 100,000 and 650,000 cells (Earth world of seed 1,
best of 5), for M3-05's budget.

## Acceptance criteria

- [x] Steady state on a synthetic 1D channel (50 cells, one ocean outlet,
      `L = 1 km`, `Q_i` the number of cells upstream × 1 m³/s, uniform
      `U = 1 mm/yr`, default K and m), for `Δt` of 10⁴, 10⁵ and 10⁶ yr: the
      slope of every land cell to its base converges to `S = U / (K Qᵐ)`
      within 1e-6 relative.
- [x] Unconditional stability: time steps of 10⁴, 10⁵ and 10⁶ years (in s),
      on the surfaces of the Tests notes, never create NaN or infinity; no
      cell drops below its receiver (every eroded cell ends strictly above
      the water surface of its receiver, the others are only uplifted); with
      no uplift or uniform uplift, no new pits and no new closed
      depressions.
- [x] Without uplift, total land volume decreases strictly at every step
      (20 steps of 10⁴ and of 10⁶ yr on the same surfaces), no height rises
      and no land cell becomes ocean.
- [x] Before and after images of `surface_height` on a 100,000-cell planet
      (Earth, seed 1, default settings) viewed and described (dendritic
      valleys appear).
- [x] `mise run ci` passes.

## Verification

```sh
mise exec -- cargo test -p kuniumi-sim erosion
mise exec -- cargo test -p kuniumi-sim --release erosion
mise exec -- cargo test -p kuniumi-world
```

Images in `out/M3-03/` from the scratch tool (see Images and timing).

## Result

What was done:
- `kuniumi-world`: `ErosionSettings` (`erodibility`, `discharge_exponent`,
  `time_step_s`, `steps`), its documented defaults, `Default` and
  `validate`, re-exported from the crate root; 3 tests.
- `kuniumi-sim`: in `stages::erosion`, `PLACEHOLDER_PRECIPITATION_M_S`,
  `Network`, `incise`, `receiver_distances_m`, `ErosionInputs` and `erode`,
  as planned; 30 tests in `stages/erosion/tests.rs`. `drainage::next_up` is
  now `pub(crate)` for the solver; the documentation of `drainage` points to
  `erode`.
- The split agreed with the owner: M3-07 (draft, after M3-03, before M3-05)
  adds `Settings::erosion` and its records; M3-05 now depends on it.
- `CHANGELOG.md` (Added); a note from this ticket in M3-05 (item 5 below).

Deviations from the plan: none in the code. In the tests: `one_step_by_hand`
uses `F` = 1, 4 and 8 (not 1, 3 and 7) so that every division is exact;
the reference spill heights of the stability test come from Priority-Flood
without ε written in the test rather than from the Planchon–Darboux
iteration (exact as well, and much faster); the stability and volume tests
run on two more surfaces, uniform random heights full of closed
depressions (sea level −900 m, review round 1, R1-4), so that the lake
rule is exercised, and the dry Earth world, whose cells at the sea level
are fixed (R1-3); the island test runs at three radii with distances and
discharges computed in the test (R1-1, R1-2); three tests were added
(`extreme_settings_stay_finite`, `rounding_never_raises_an_eroded_cell`,
the second because the clamp to `a` survived every other test as a
mutant, and `erode_steps_compose`). The images show the surface after
100 steps (1 Myr) as well as after the default 1,000 steps (10 Myr).
After review round 1 `erode` checks the heights, the uplift and the sea
level before the first step (R1-5), and `incise` takes the limit `b⁺`
where `F` overflows instead of computing `∞ / ∞` (R1-6); the golden value
is unchanged.

Evidence:

1. Criterion 1 (`steady_state_of_a_1d_channel`): from land at 1 m, the
   chain is steady (no height changes by more than 1e-12 m) after 219, 40
   and 16 steps of 10⁴, 10⁵ and 10⁶ yr; the divide stands at 602.728 m in
   each case and the worst relative error of the 49 slopes against
   `U / (K Qᵐ)` (std `powf`) is 8.0e-14, 1.0e-14 and 2.2e-15 (bound 1e-6).
   The fixed point of the update is the steady state exactly:
   `(h_i − b)(1 + F) = h_i − b + Δt U` gives `(h_i − b)/L = U/(K Qᵐ)`.
2. Criterion 2 (`large_time_steps_are_stable`): on the Earth worlds of seeds
   1–3, the dry Earth world of seed 1, the smooth random surface and the
   uniform random heights, without uplift, with 1 mm/yr everywhere and
   (Earth) with their own `uplift_rate`, 10 steps each of 10⁴, 10⁵ and
   10⁶ yr (480 steps): every height finite; every ocean cell unchanged bit
   for bit (on the dry world, the cells exactly at the sea level); every
   land cell with `a > W_r` strictly above the water surface of its
   receiver and not above `a`, every other one exactly `a` (cells under
   water in 0–16,526 cell-steps per run, most on the uniform random
   heights at 10⁴ yr, 1,001 on the dry world, 412 on the smooth surface);
   land stays land; with no or uniform uplift no land cell lost its last
   strictly lower neighbour and no cell came to lie in a closed
   depression. The Earth uplift dammed one pit (seed 3, 10⁶ yr).
   `extreme_settings_stay_finite`: `Δt K = 0` with `Qᵐ = ∞`, `Δt K = ∞`
   with `Qᵐ = 0`, `F = ∞`, and `F = ∞` with `a − b = ∞` stay finite.
3. Criterion 3 (`land_volume_decreases_without_uplift`): over 20 steps the
   land volume falls strictly at every step, no height rises and the land
   cells stay the same, on all six surfaces; Earth seed 1: 1.5333e17 m³ to
   1.2076e17 m³ with steps of 10⁴ yr and to 3.1e6 m³ with steps of 10⁶ yr
   (the other surfaces alike: the dry world 2.5177e18 m³ to 2.2055e18 m³
   and 1.1e8 m³, the uniform random heights 4.6040e17 m³ to 2.9645e17 m³
   and 6.7e6 m³).
4. By hand and end to end: `one_step_by_hand` (thirteen cells: the ocean
   fixed, a coast graded to the sea level rather than to its −4,000 m
   receiver, 51; a land receiver, 57; a lake floor only uplifted, 41; a
   river entering the lake graded to its surface, 58; a rounding to the
   base, one ULP above 57; no discharge, 301; a pit, 503, and its donor,
   551.5; a cell exactly at the sea level with uplift, fixed at 0, and its
   donor, 10; a pit sinking below the sea level, −5, and its donor graded
   to the sea level, 10), `rounding_never_raises_an_eroded_cell`,
   `receiver_distances_by_hand`, and `steady_state_of_an_island` (867 land
   cells, rain varying with z, `P (1 + z/2)`, steady after 34, 36 and 39
   steps of 10⁶ yr at radii of 1,500, 6,371 and 13,000 km, highest 1,418,
   1,638 and 1,759 m, every slope within 1.2e-12 relative of `U / (K Qᵐ)`
   with the distance and the discharge computed in the test,
   `receiver_distances_m` within 1.6e-14 rad × R of the great-circle
   distance).
   `erode_reports_progress_per_step`, `erode_steps_compose`; 17 panic
   tests; `kuniumi-world`: defaults, the whole domain accepted, every bad
   value rejected by name in field order.
5. Criterion 4, images (`out/M3-03/`, a scratch tool outside the
   repository, `<scratchpad>/m3_03_maps`, a Cargo package with path
   dependencies on the workspace crates; Earth world of seed 1 at 100,000
   cells, default settings with its `uplift_rate` and 1 m/yr of rain;
   elevation colours of the world maps, −11 to +9 km about the sea level,
   with hill shading on land), viewed:
   - `before_100000.png`: flat continents (+1,017 m above the sea level,
     the uniform interiors of M2-05) with smooth plateaus along the
     collision belts and arcs (up to 6,320 m).
   - `after_100_100000.png` (1 Myr): dendritic valley networks cover every
     continent, cut back from the coasts into the plateaus, branching
     inland; the belts are dissected into ridges and valleys; the highest
     point has risen to 8,564 m where uplift outpaces erosion on the
     steepening belts. Close-ups: `interior_after_100_100000.png` (3.0° N,
     89.3° E, 33 cells from the sea), branching dark valleys between
     lighter interfluves; `belt_after_100_100000.png` (13.4° S, 175.6° W),
     the plateau cut into a ridge with spurs.
   - `after_1000_100000.png` (10 Myr, the default run): the interiors,
     which have no uplift, are cut down to the sea level (median land
     elevation 1,017 m → 0 m, mean 1,096 m → 79 m, 99th percentile
     4,811 m → 2,274 m); only the uplifted belts keep their relief
     (highest 5,960 m), dissected (`belt_after_1000_100000.png`).
     `change_100000.png`: about −1,000 m over every interior, up to −2 km
     along the belts.
   - `discharge_after_{100,1000}_100000.png` and their close-ups: tree-like
     networks on every continent, tributaries joining into trunk rivers
     that reach the coast, the same at 1 and 10 Myr (the network forms in
     the first knickpoint sweep and stays). On some continents the trunk
     rivers run nearly parallel to each other, cut into the formerly flat
     interiors along the ε-routing of M3-02 (M3-06).
   The planing of the interiors follows from the physics so far: with no
   uplift the base level is the sea, and a headwater cell relaxes towards
   its receiver with the time `L / (K Qᵐ)` ≈ 0.3 Myr at 100,000 cells;
   isostatic rebound (M7-03), deposition (M7-02) and harder old crust
   (M7-04), which slow it, come later. Noted in M3-05, which runs the
   default settings.
6. Golden value (`golden_fingerprint_of_an_eroded_earth_world`, the Earth
   world of seed 1 at 10,000 cells after 10 default steps):
   `0x73A6_193E_F4E0_5D72`, recorded after 1–5; the same in debug and
   release. The golden values of M3-01 and M3-02 are unchanged.
7. Timing (release, this machine, one step of `erode` on the Earth world
   of seed 1, best of 5): 7.8 ms at 100,000 cells (8.6 ms after 100
   steps), 80.9 ms at 650,000 cells (93.4 ms after 100 steps): about
   90 s for the default 1,000 steps, the design's 0.1 s per step.
8. Test-first: against the stubs (commit "tests of the incision and
   ErosionSettings against stubs") 23 of the 25 tests of `erosion` failed
   (`the_unbroken_chain_runs` and `erode_steps_compose` pass with a no-op)
   and `erosion_settings_reject_values_out_of_their_domain_by_name`.
   Mutations, each caught by the final tests (release; 23 mutants, six
   added in review round 1 and two in round 2): the ocean's floor as the base (8
   tests), no lake rule (2), the lake's floor as its surface (2), no
   `next_up` guard (4), no clamp to `a` (1), no guard against `0 · ∞` (1),
   the ticket's `(a + F b)/(1 + F)` form (1, the golden value: the forms
   agree in exact arithmetic), a pit not uplifted (1), ocean uplifted (5),
   edge lengths for distances (3), the area for the discharge (3), routing
   on the unfilled surface (2), progress `k/steps` (1), `m = 0.5` (3), no
   stack-order check (1); no branch for `F = ∞` (1), the runoff scaled by
   `R · 6.371e6` instead of `R²` (1), the distance of the first neighbour
   for every receiver (2), ocean only strictly below the sea level (2), a
   pit's water surface at the sea level (1), `erode` not checking its
   heights (1); the rain of cell 0 for every cell (1), a sunken pit's
   water surface below the sea level (1).
9. `mise run ci` passes: fmt, clippy `-D warnings`, rustdoc, the browser
   checks, the tests in debug and release, MSRV 1.85.0, cargo-deny.

Review:
- Round 1 (APPROVE): five minor findings and two nits, all fixed. R1-1,
  nothing tested how the results scale with the radius (a runoff of
  `P a R · 6.371e6` passed every test): the island test now runs at 1,500,
  6,371 and 13,000 km with the distance and the discharge computed in the
  test. R1-2, `receiver_distances_by_hand` could not tell neighbours apart
  on the octahedron: the island test checks `receiver_distances_m`
  against the great-circle distance on the 10,000-cell grid. R1-3, no
  cell lay exactly at the sea level: `one_step_by_hand` has one, with
  uplift, and the dry Earth world joins the stability and volume tests.
  R1-4, the uniform random surface had 12 depression cells, not "full of
  closed depressions": its sea level is now −900 m (16,526 under-water
  cell-steps at 10⁴ yr). R1-5, `erode` did not check heights, uplift and
  sea level with 0 steps: it does now, before the first step, with three
  panic tests. R1-6 (nit), `∞ / ∞` where both `a − b` and `F` overflow:
  `incise` takes `b⁺` where `F` is infinite, tested. R1-7 (nit), an
  unwrapped doc line in `drainage.rs`, rewrapped. The reviewer re-ran CI,
  the tests, the timing and the images (byte-identical) and probed the
  island steady state at four radii and a dry world.
- Round 2 (APPROVE): R1-1 to R1-7 resolved; the reviewer re-ran CI, the
  tests, the timing, the images (byte-identical), the 21 mutants and 13 of
  its own. R2-1 minor, no test used rain that varies between cells (the
  rain of cell 0 everywhere passed): fixed, the island test's rain varies
  with z. R2-2 nit, a pit sinking below the sea level (negative uplift)
  was untested and its water surface undocumented: fixed, two cells in
  `one_step_by_hand` and a sentence in the documentation of `incise`.
  R2-3 nit, an unwrapped line in `docs/tickets/README.md`: rewrapped.

Follow-ups: M3-07 (split off, draft).

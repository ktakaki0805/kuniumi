# M4-07 Climate stage and monthly outputs

- Status: done
- Milestone: M4
- Depends on: M4-06, M4-09, M4-14
- Crates: kuniumi-sim; kuniumi-io, kuniumi-cli, kuniumi-viewer (tests)
- Design: 気候モデル / 計算の流れ; 実装ロードマップ / 計算コストの目安
- Decisions: D-08, D-22

## Goal

The Climate stage runs in the pipeline before Erosion and writes the monthly
temperature and precipitation, the winds and the currents of every cell.
Erosion keeps its uniform placeholder precipitation for now (decided with
the owner during the implementation, below; ticket M4-19).

## Scope

In:
- `stages::climate::run`, a `StageFn` in the pipeline's table in place of
  `Implementation::Skipped`: the models of M4-02 to M4-06 composed once, in
  the order of the design and as they were calibrated (below, "The chain"),
  each run to its own periodic state: insolation → energy balance →
  pressure and winds → currents (annual) → sea-surface temperature with
  currents and upwelling → moisture and precipitation. A model that does not
  reach its periodic state keeps the months of its last year (as the models
  already do), documented with the stage.
- Fields (below, "Fields"): `temperature_m01` … `temperature_m12` (K),
  `precipitation_m01` … `precipitation_m12` (m/s), their annual means
  `temperature_annual` and `precipitation_annual`, the annual-mean surface
  wind `wind_east`, `wind_north` (m/s) and the currents `current_east`,
  `current_north` (m/s): 30 fields.
- The Erosion stage keeps the uniform placeholder of 1 m/yr: the
  climate's precipitation is exactly 0 on 13 to 35 % of the land, whose
  uplift would never erode (below, "Erosion"; chosen with the owner). The
  values of every existing field stay bit for bit.
- Maps: the world directory's map of every field (as for every field
  today); the annual means and two contrasting months (January and July) of
  the temperature and the precipitation, the winds, the currents and the
  rivers viewed at 650,000 cells.
- Golden fingerprints for 10,000 cells, seed 1: the stage's 30 fields, in
  the stage's test and in the CLI's golden world (whose 16 other fields
  keep their values).
- Timing of `kuniumi generate` at 650,000 cells, per stage.
- The tests that generate worlds kept affordable in debug builds (below,
  "Test time").

Out:
- Warnings for a climate that has not converged (split into ticket M4-16
  with the owner's approval). Decided with the owner in M4-02 for this ticket: a planet
  whose energy balance has not converged within its 100 years
  (`energy_balance::Temperature::converged`) keeps the last year's
  temperatures and is reported as a warning, with the last change, saying
  that the climate may be in transition (near an ice–albedo bifurcation the
  last year can lie mid-way into a snowball: Earth at 1,296 W/m² changed by
  11.5 K in its 100th year); the moisture's `converged`, `unsettled_steps`
  and `failed_solves` and the currents' solves likewise. The temperatures
  have no fixed bound at the corners of the v1 range (monthly means from
  about 95 K to about 11,600 K, M4-02).
- The coupled periodic state of temperature, winds and currents (follow-up
  M4-17; measured below, "The coupled loop").
- Colour scales shared by the months, a square-root scale for the
  precipitation (follow-up M4-18).
- Erosion from the climate's precipitation (follow-up M4-19, after the dry
  land has rain, M4-08 or a background precipitation, and relief has a
  limit).
- Climate–erosion alternation (M7-01), latent heat (M4-13), sea ice
  (M4-11), the solver of the periodic state (M4-15), Earth's real
  topography (M4-08).

## Implementation notes

The ticket was `draft`. It left open the periodicity criterion and the
maximum number of model years, which fields are monthly and which annual
only, and whether wind and currents are two scalar fields each. A scratch
tool outside the repository (`<scratchpad>/m4_07_tool`, a Cargo package with
path dependencies on `main`, the public API of the climate modules) measured
the chain and the coupled loop (release unless said otherwise, this
machine, the Earth-parameter world of seed 1 through the Crust stage).

### Chosen with the owner

Approved by the owner as recommended when the ticket was refined:

1. **The chain**: the models composed once, as M4-05 and
   M4-06 calibrated and tested them (winds and currents of the energy
   balance without currents), each to its own periodic state; not the
   coupled loop, which does not settle and leaves the calibration
   (measurements below). The coupled periodic state, with its
   recalibration, becomes follow-up M4-17.
2. **Split**: the warnings of an unconverged climate need a
   place in the world that survives the world directory (`World::warnings`,
   `world.json` format version 5, read back) and the CLI's report: an
   independent change of `kuniumi-world`, `kuniumi-io` and `kuniumi-cli`
   (about 250 lines with tests): ticket M4-16 (draft, after M4-07), so
   that M4-07 stays at about 950 lines (about 1,200 with it).
3. **Fields**: monthly temperature and precipitation (the inputs of Köppen,
   biomes and runoff, M5-01, M5-04, M5-05); wind and currents annual only,
   two scalar fields each (as `plate_velocity_east`/`_north`); no separate
   sea-surface temperature, which is the temperature fields at the ocean
   cells, bit for bit (below).
4. **Follow-ups**: M4-16 (warnings), M4-17 (coupled periodic state),
   M4-18 (climate maps on shared scales), all `draft`, in
   `docs/tickets/README.md`.

### Measurements of the refinement

The chain below, Earth, seed 1:

| Cells | Build | Energy balance | Circulation, currents | Sea surface | Moisture | Total |
| --- | --- | --- | --- | --- | --- | --- |
| 1,000 | release | 0.23 s (12 years) | 0.00 s | 0.18 s (9) | 0.04 s (3) | 0.45 s |
| 1,000 | debug | 3.2 s | 0.03 s | 2.5 s | 0.5 s | 6.3 s |
| 10,000 | release | 2.4 s (10) | 0.07 s | 2.7 s (9) | 0.5 s (3) | 5.7 s |
| 10,000 | debug | 22.5 s | 0.6 s | 25.7 s | 5.2 s | 54 s |
| 650,000 | release | about 4.3 min (M4-09) | | about 5 min (M4-09) | 81 s (M4-14) | about 11 min |

The 128 corners of the v1 box of `the_bounds_of_the_planet_parameters_
generate_finite_worlds` (1,000 cells, seed 3, release): 29 s for the chain
of all of them; every value finite. Not periodic: the energy balance of 4
(corners 56–59: 8 h, 0.1 bar, the longest orbit, 2,500 W/m² and the most
water, at either radius and gravity; last change 4.7 to 14.6 K after 100
years, the thin, fast planets of M4-02), the sea surface of 2 of them (57,
59: 7.4 and 6.6 K), the moisture of 11 after its 10 years (all without
water, whose only ocean is the cells at the sea level, at the largest
radius: the cold ones by 2 to 4 × 10⁻⁶ kg/m² in a column of under
10⁻³ kg/m², the hot ones by up to 55 kg/m² as their atmosphere fills). No
sub-step left unsettled, no failed solve, every basin's currents
converged.

`mise run test-debug` on `main` takes 2 min 19 s (the `kuniumi-sim` tests
103 s).

### The coupled loop (measured, not chosen)

The design's loop updates the winds from each month's temperature and the
currents once a year from the annual winds, "年ごとの差が十分小さくなるまで".
M4-05 and M4-06 left it to this ticket; they calibrated `A`, `D` of the sea
surface and the gust speed with the winds and currents of the energy
balance without currents. The scratch tool iterated the chain as a fixed
point: winds and currents from the sea surface's periodic temperatures,
the sea surface again with them, and so on:

Sea surface 1 is the calibrated chain (winds and currents of the energy
balance without currents); sea surface `k + 1` has the winds and currents
of sea surface `k`:

| | 10,000 cells, sea surface 1, 2, 3, 4 | 40,000 cells, 1, 2, 3, 4 |
| --- | --- | --- |
| largest change of a monthly sea-level temperature from the one before (the energy balance before 1) | 10.3, 3.2, 2.7, 2.5 K | 8.5, 5.4, 7.5, 9.9 K |
| RMS change | 2.0, 0.75, 0.22, 0.13 K | 1.5, 1.5, 0.50, 0.47 K |
| largest change of a monthly wind, of its winds from those before | 9.8, 7.7, 7.2, 7.1 m/s | 15.5, 26.0, 31.2, 57.3 m/s |
| mean surface temperature | 288.98, 289.31, 289.28, 289.25 K | 287.99, 288.85, 289.00, 289.07 K |
| `P₂` coefficient (target −28.0 K on 40,000 cells) | −24.6, −23.6, −23.8, −23.9 K | −27.97, −25.5, −25.0, −24.6 K |
| July thermal equator of its winds | 6.4, 7.3, 7.1, 6.8° | 3.3, 5.2, 6.2, 7.2° |
| global precipitation (with the energy balance's winds: 2.54 and 2.88) | 2.67, 2.71, 2.70, 2.69 mm/day | 3.17, 3.43, 3.53, 3.63 mm/day |

The loop does not settle in the largest change (local oscillations of
several kelvin), and on the 40,000-cell world of the calibrations it moves
away from them: 1.1 K warmer, a flatter `P₂`, 26 % more precipitation, a
thermal equator that swings further each iteration and winds that change
more each iteration. Making it converge
(damping, a smoother circulation) and recalibrating `A`, `D` and `U_g` with
it is its own work: follow-up M4-17. This ticket composes the models once,
as calibrated.

### The chain (`stages::climate::climate`)

For the planet `params`, the grid, the Crust stage's `surface_height` `h`
(m) and the sea level `z` (m), with the insolation of
`insolation::monthly_insolation`:

1. `T₀ = energy_balance::monthly_temperature(params, grid, h, z)`: the
   energy balance without currents (M4-02) to its periodic state
   (`PERIODIC_TOLERANCE_K` = 10⁻³ K, at most `MAX_YEARS` = 100 years).
2. `W = circulation::monthly_circulation(params, grid, T₀)`: monthly
   pressure and surface winds (M4-03).
3. `C = ocean::ocean_currents(params, grid, h, z, W, T₀)`: the currents of
   the annual-mean wind stress (M4-04), solved once.
4. `S = sea_surface::monthly_sea_surface_temperature(params, grid, h, z, W,
   T₀, C)`: the energy balance with the heat of the currents and of the
   monthly coastal upwelling (M4-05) to its periodic state (the same
   tolerance and years).
5. `M = moisture::monthly_moisture(params, grid, h, z, W, S.temperature)`:
   vapour, evaporation and precipitation (M4-06) to their periodic state
   (`PERIODIC_TOLERANCE_KG_M2` = 10⁻⁶ kg/m², at most 10 years).

Exactly the composition of the tests of M4-05 and M4-06
(`test_worlds::sea_model`), so their calibrations and golden values hold
for the stage. The periodicity criterion of each model is its own; the
stage adds none and no outer iteration. The climate is that of the relief
of the Crust stage, before erosion (the climate–erosion alternation is
M7-01).

```rust
/// The climate of a planet: the models composed once (module documentation).
pub struct Climate {
    /// The energy balance without currents (M4-02): the winds' and currents' temperatures.
    pub energy_balance: Temperature,
    /// The monthly pressure and winds of `energy_balance` (M4-03).
    pub circulation: Circulation,
    /// The currents of their annual-mean stress (M4-04).
    pub currents: OceanCurrents,
    /// The sea surface with currents and upwelling (M4-05): the climate's temperatures.
    pub sea_surface: SeaSurface,
    /// The moisture of `circulation` over `sea_surface` (M4-06).
    pub moisture: Moisture,
}

/// The climate of the surface `heights_m` (m, on the datum of `sea_level_m`)
/// of the planet `params` (valid) on `grid`; `progress` receives the
/// fraction done after each model.
pub fn climate(params: &PlanetParams, grid: &SphereGrid, heights_m: &[f64],
               sea_level_m: f64, progress: &mut dyn FnMut(f64)) -> Climate;

/// Runs the Climate stage on `world` (a `StageFn`): `climate` of its
/// `surface_height`, then adds the fields below in their order.
pub fn run(world: &mut World, progress: &mut dyn FnMut(f64)) -> Result<(), SimError>;

pub const TEMPERATURE_ANNUAL: &str = "temperature_annual";
pub const TEMPERATURE_MONTHS: [&str; MONTHS] = ["temperature_m01", …, "temperature_m12"];
pub const PRECIPITATION_ANNUAL: &str = "precipitation_annual";
pub const PRECIPITATION_MONTHS: [&str; MONTHS] = ["precipitation_m01", …];
pub const WIND_EAST: &str = "wind_east";
pub const WIND_NORTH: &str = "wind_north";
pub const CURRENT_EAST: &str = "current_east";
pub const CURRENT_NORTH: &str = "current_north";
```

In `stages/climate.rs` (the module root, which only declares the models
now), tests in `stages/climate/stage_tests.rs`. `climate` is public so that
the Earth calibration (M4-08) can run it on Earth's relief outside the
pipeline.

Progress: `climate` reports 0.40 after the energy balance, 0.42 after the
circulation and the currents, 0.87 after the sea surface (their shares of
the time at 650,000 cells); the runner reports 0 and 1. Finer reports
(per model year) would need a callback in the periodic-state drivers: not
in this ticket.

### Fields

Added by `run` in this order; months in calendar order, `k` = 1 … 12:

| Field | Unit | Values |
| --- | --- | --- |
| `temperature_annual` | K | `(1/12) Σ_k T_k`, summed in month order |
| `temperature_mNN` | K | `S.temperature.surface_k[k]`: the monthly-mean surface temperature at the cell's height (lapse rate 6.5 K/km above the sea level); over the ocean the sea-surface temperature |
| `precipitation_annual` | m/s | `(1/12) Σ_k P_k` |
| `precipitation_mNN` | m/s | `M.precipitation_m_s[k]`, as liquid water |
| `wind_east`, `wind_north` | m/s | `(1/12) Σ_k` of the monthly surface wind of `W` |
| `current_east`, `current_north` | m/s | `C.current_east_m_s`, `C.current_north_m_s`, 0 on land |

Annual means first so that the viewer's first climate layers are the
overview. The months are equal (twelfths of the orbital period,
M4-01), so the mean of the monthly means is the annual mean.

- Sea-surface temperature: no field of its own. An ocean cell (height at
  most the sea level) has no height above the sea, so its surface
  temperature is its sea-level temperature, the mixed layer's, bit for bit
  (M4-05); the temperature fields at the ocean cells are the sea-surface
  temperature, as the stage documents.
- Wind and currents: two scalar fields each, east and north, as the
  plates' velocity. The currents are annual by construction (the design:
  "海流は年平均の風で年1回"); the winds are written as their annual mean,
  so the monsoons show in the monthly precipitation, not in the winds.
  Monthly winds would add 24 fields (125 MB at 650,000 cells).
- Not written: the energy balance without currents (an intermediate), the
  pressure, the evaporation, the vapour, the upwelling (diagnostics of the
  models; later tickets can add what they need).
- Size at 650,000 cells: 30 fields of 5.2 MB, 156 MB more `.kfield` files,
  and 30 maps; the world has 16 fields today.

### Erosion (changed with the owner during the implementation)

The plan had `erosion::run` read `precipitation_annual` instead of
`PLACEHOLDER_PRECIPITATION_M_S`; it was implemented and measured, then
withdrawn: the Erosion stage keeps the placeholder, its code and tests
unchanged but for documentation, until ticket M4-19 (draft) makes it
follow the climate once the dry land has rain (M4-08's calibration or a
background precipitation) and relief has a limit. Chosen by the owner
among: the placeholder for now (chosen); the climate's rain with the
runaway peaks below as a known limitation; a floor under the rain of
erosion; a background precipitation in the moisture model first.

Found when the plan's erosion ran with the default 1,000 erosion steps (10 Myr), the Earth-parameter worlds
of seeds 1–3 eroded under the climate's `precipitation_annual` reach peaks
of 44 to 101 km on 10,000 cells and 39 to 106 km on 40,000 (`main`, under
1 m/yr everywhere: 2.3 to 5.6 and 4.0 to 9.8 km); 1,000 cells, seed 3:
38.9 km against 3.4. The moisture model of M4-06 gives exactly 0 on 26 to
35 % of the land on 10,000 cells and 13 to 33 % on 40,000 (below
50 mm/yr: 39 to 54 % and 27 to 47 %; Earth's hyper-arid land is about
7.5 %), uplifting cells among them (13 to 45 and 84 to 120), and a cell
without rain has no discharge, so its uplift (up to 10 mm/yr) is never
eroded: 10 Myr × 10 mm/yr = 100 km. With a floor under the rain of erosion
(scratch build, not chosen): 0.1 m/yr gives 6.6 to 27 km, 0.3 m/yr 3.8 to
17 km.

### Pipeline

- `pipeline::implementation(StageId::Climate)` is
  `Implementation::Stage(stages::climate::run)`. No stage of the real
  table is skipped any more; `Implementation::Skipped` stays for tables of
  tests (its own tests use it). Docs of `pipeline`, `lib.rs`, `generate`,
  `next_stage`, `stages/climate.rs`, `erosion.rs` and `Settings` that say
  that Climate is skipped until M4-07 are updated.
- A world generated now has the stages `grid`, `tectonics`, `crust`,
  `climate`, `erosion` and stops at Hydrology.

### Test time

Every `kuniumi_sim::generate` now runs the climate: 6 s at 1,000 cells and
54 s at 10,000 in debug builds. So:

- `the_bounds_of_the_planet_parameters_generate_finite_worlds` (Earth and
  128 corners at 1,000 cells, now with the climate): split into two tests
  (the corners at the lower and at the upper radius, so they run in
  parallel); release builds run every corner, debug builds 16 (those
  whose index halved is a multiple of 8, 8 at each radius) and Earth, as
  debug builds already run a tenth of the erosion steps.
- Tests that need eroded worlds but not the climate (the drainage metrics'
  `eroded_worlds`, 4 worlds of 10,000 cells) use a crate-internal helper,
  `generate_without_climate` (`#[cfg(test)]`): the pipeline's table with
  the Climate stage skipped. Erosion does not use the climate, so their
  worlds, and their asserted values, stay bit for bit those of today.
- `kuniumi-viewer`'s 100,000-cell scene test (`app.rs`, `scene_of`, relief
  only, erosion steps 0) runs the pipeline through the Crust stage
  (`pipeline::run` with the stages after Crust missing); its 1,000-cell
  tests keep `generate`.
- `kuniumi-io`'s `reads_back_a_generated_world_bit_exactly` generates
  2,000 cells instead of 20,000 (about 2 minutes of climate in debug; the
  round trip does not depend on the size); its field counts (16) become 46.
- The CLI's golden worlds (10,000 cells, four runs in parallel) and its
  1,000–1,500-cell runs keep their sizes.

### Tests

New, `stages/climate/stage_tests.rs` (the stage run once per process on the
10,000-cell Earth world of seed 1 through the Crust stage, shared):

- `the_stage_composes_the_models`: every field bit for bit against the
  models composed by hand from `test_worlds` (`earth_10k`,
  `earth_sea_10k`, and `monthly_moisture` of their winds and sea surface),
  the annual means recomputed from the months; the fields' names, units,
  kinds and order; `world.stages` gains `climate`, the sea level and the
  other fields unchanged.
- `the_climate_is_earth_like`: every value finite, every temperature
  positive, every precipitation non-negative, the currents 0 on land; the
  area-weighted global annual means: surface temperature within 1.5 K of
  the 288 K of the calibration on 40,000 cells (the scratch tool: 288.98 K
  on this world) and precipitation 0.8 to 1.2 m/yr (0.93 m/yr);
  January warmer than July south of 30° S and colder north of 30° N.
- `climate_reports_its_progress`: 0.40, 0.42, 0.87 in order from `climate`;
  the stage's reports within [0, 1] and increasing.
- `existing_outputs_leave_the_world_unchanged` (one of the 30 names present
  → `SimError::Field`, nothing added, no progress) and
  `the_stage_needs_the_heights` (panic).
- `golden_fingerprints_of_the_climate_stage`: FNV-1a of the 30 fields,
  recorded after the tests above pass and the maps are viewed.

Changed:

- `lib.rs` and `pipeline/tests.rs`: the stages that run (Climate now), the
  field list of a generated world, `next_stage` without the skipped
  Climate, `the_real_table`.
- `kuniumi-cli/tests/cli.rs`: `GOLDEN_FIELDS` with the 30 climate fields
  after the Crust stage's, the 16 others unchanged; the 0-step world keeps
  the climate's fingerprints too (the climate is that of the Crust stage's
  relief); `GOLDEN_CRUST_SURFACE_HEIGHT`, the sea level and the metrics
  unchanged; the progress lines and files of the climate.
- `kuniumi-io` read tests and `kuniumi-viewer` tests: field counts and
  stage lists (the panel's test of 60 fields: 46 generated and 14 copies).

### Images and timing

`cargo run --release -p kuniumi-cli -- generate --seed 1 --out out/M4-07`
(650,000 cells): the maps of `temperature_annual`, `temperature_m01`,
`temperature_m07`, `precipitation_annual`, `precipitation_m01`,
`precipitation_m07`, `wind_east`, `wind_north`, `current_east`,
`current_north`, and `rivers.png`, viewed and described
(each map over its own range, which the description gives); the time of
each stage from the progress lines' time stamps, and the total.

### World-breaking

Every world gains the stage `climate` and its 30 fields (and a world
directory their 60 files); every existing field keeps its values bit for
bit.

### Size

About 850 changed lines without the tickets: the stage and its
documentation 280, its tests 340, the pipeline, the Erosion stage's
documentation and the other crates' tests 200, documents (`CHANGELOG.md`,
`design-ja.md` 計算コストの目安 with the measured time) 50.

### Design document

`実装ロードマップ / 計算コストの目安`: the measured time of the climate stage
and of `kuniumi generate` at 650,000 cells. `気候モデル / 計算の流れ` keeps the
loop as the design; a sentence says that the stage composes the models
once (each to its own periodic state) until M4-17 couples them.

## Acceptance criteria

- [x] `the_stage_composes_the_models`: the 30 fields equal the models
      composed by hand bit for bit, with the annual means of the months,
      in the order, units and kinds of "Fields".
- [x] `the_climate_is_earth_like`: on the 10,000-cell Earth world of seed 1,
      every value finite, temperatures positive, precipitation
      non-negative, currents 0 on land; global annual-mean surface
      temperature within 1.5 K of 288 K and precipitation 0.8–1.2 m/yr;
      January warmer than July south of 30° S and colder north of 30° N.
- [x] The values of every field that existed before keep their golden
      fingerprints (the CLI's 16, the Erosion stage's own), the Erosion
      stage keeping the uniform 1 m/yr (owner's choice).
- [x] The 129 planets of the bounds test generate finite worlds through
      the Climate and Erosion stages (all in release builds, Earth and
      16 corners in debug builds).
- [x] `kuniumi generate --seed 1` (650,000 cells) writes the 30 climate
      fields with their maps; the maps listed in "Images and timing" viewed
      and described.
- [x] Golden fingerprints for 10,000 cells, seed 1: the climate stage's 30
      fields, in the stage's test and in the CLI's golden world (new fields:
      world-breaking, in `CHANGELOG.md`).
- [x] Timing of `kuniumi generate` at 650,000 cells recorded per stage
      (release).
- [x] `mise run ci` passes.

## Verification

```sh
cargo test --release -p kuniumi-sim climate
cargo test -p kuniumi-cli
cargo run --release -p kuniumi-cli -- generate --seed 1 --out out/M4-07
```

## Result

What was done:
- `kuniumi-sim::stages::climate` (`climate.rs`): the Climate stage, `run`
  (a `StageFn`), and `climate` with its result `Climate` (the energy
  balance without currents, its circulation, the currents, the sea surface
  with currents and upwelling, the moisture), the models of M4-02 to M4-06
  composed once as calibrated; the field names `TEMPERATURE_ANNUAL`,
  `TEMPERATURE_MONTHS`, `PRECIPITATION_ANNUAL`, `PRECIPITATION_MONTHS`,
  `WIND_EAST`, `WIND_NORTH`, `CURRENT_EAST`, `CURRENT_NORTH`, and
  `PROGRESS`. The module documentation gained "The Climate stage" (the
  chain, the coupled loop not taken, the unconverged states, erosion),
  "Fields" and "Cost". 7 tests in `climate/stage_tests.rs`.
- `pipeline::implementation(Climate)` is the stage; the documentation of
  `pipeline`, `lib.rs` (`generate`, `next_stage`) and `erosion.rs` (which
  keeps `PLACEHOLDER_PRECIPITATION_M_S`, until M4-19). `lib.rs`: the
  test helper `generate_without_climate`, the test
  `the_climate_changes_no_other_field`, the bounds test split in two
  (`…_at_the_lower_radius`, `…_at_the_upper_radius`) with 16 corners in
  debug builds; the pipeline's and the drainage metrics' tests follow.
- Other crates' tests: `kuniumi-cli` (`WORLD_FIELDS`, `GOLDEN_FIELDS` with
  the 30 climate fields, the climate's progress lines and `world.json`
  lines, the 0-step world's climate), `kuniumi-io` (the round trip on
  2,000 cells instead of 20,000, 46 fields), `kuniumi-viewer` (the
  100,000-cell scene through the Crust stage, the panel's 60 fields with
  14 copies).
- `docs/design-ja.md` (気候モデル / 計算の流れ: the stage's composition,
  no SST field, annual winds, erosion's rain; 実装ロードマップ /
  計算コストの目安: the measured cost), `README.md` (what `generate`
  writes), `CHANGELOG.md` (Added, Changed, World-breaking); follow-ups
  M4-16, M4-17, M4-18, M4-19 (draft) in `docs/tickets/README.md`.

Deviations from the plan:
- Erosion (changed with the owner during the implementation, above,
  "Erosion"): the Erosion stage keeps the uniform 1 m/yr; the plan's
  `erosion::run` from `precipitation_annual` was implemented, measured
  (runaway peaks of 39 to 106 km) and withdrawn. So the criterion on the
  erosion's discharge became "every existing field keeps its values", and
  the world-breaking change is only the new stage and fields.
- `PROGRESS` is the measured shares at 650,000 cells, 0.35, 0.41 and 0.88,
  not the planned 0.40, 0.42 and 0.87 (from the estimates).
- Tests beyond the plan: `the_field_names_are_the_ticket_s`,
  `the_climate_changes_no_other_field` (for the owner's choice on erosion)
  and `the_orbits_at_their_bounds_generate_finite_worlds` (review round 1); `climate_reports_its_progress` checks the stage's reports
  equal `PROGRESS`, not only that they increase.
- Images: besides the world directory's maps (each over its own range),
  the scratch tool `<scratchpad>/m4_07_maps` (a Cargo package with path
  dependencies, the real reader and colour maps) drew the temperatures on
  one range for all months, the precipitation on a square-root scale and
  the winds and currents on a diverging scale (`out/M4-07/scaled/`), to
  describe them; the world directory's precipitation maps are nearly
  uniformly dark (M4-18).

Evidence:

1. `the_stage_composes_the_models` (debug and release): the 30 fields
   equal bit for bit the models composed by hand (`test_worlds::earth_10k`,
   `earth_sea_10k`, their circulation and currents, `monthly_moisture`),
   the annual means the sums of the months divided by 12; the names,
   units (K, m/s) and kinds in the order of "Fields"; the earlier fields
   and the sea level unchanged; over the ocean `surface_k` is
   `sea_level_k` bit for bit.
2. `the_climate_is_earth_like` (10,000 cells, seed 1): every value finite,
   temperatures positive, precipitation non-negative, currents +0 on land;
   global annual-mean surface temperature 288.981 K, precipitation
   0.9263 m/yr; south of 30° S January 284.24 K against July 281.33 K,
   north of 30° N January 276.44 K against July 280.98 K.
3. Existing fields unchanged: the CLI's `GOLDEN_FIELDS` of the 16 fields,
   `GOLDEN_CRUST_SURFACE_HEIGHT`, the sea level and the reference metrics
   pass unchanged; the Erosion stage's own golden values
   (`stage_tests.rs`, `diffusion_tests.rs`, `tests.rs`) unchanged;
   `the_climate_changes_no_other_field` (1,000 cells, seed 2, 5 steps):
   with and without the Climate stage every other field and the sea level
   bit for bit.
4. The bounds tests: release, Earth and the 128 corners through the
   Climate and Erosion stages, all finite; debug, Earth and 16 corners.
   `the_orbits_at_their_bounds_generate_finite_worlds` (added in review
   round 1; debug and release): obliquity 90° with e = 0.99 at 2,500 W/m²,
   8 h, 0.1 bar and at 400 W/m², 240 h, 10 bar; obliquity 180°; obliquity
   0 with e = 0.99 and the perihelion at 270°; continental fractions 0
   and 1: all finite.
   (The corners' convergence, measured in the refinement: 4 energy
   balances and 2 sea surfaces not periodic after 100 years, 11
   moistures after 10; reported by M4-16.)
5. `kuniumi generate --seed 1 --out out/M4-07` (650,000 cells) writes the
   46 fields (97 files, 352 MB) and exits with code 2 at Hydrology.
   Viewed:
   - `temperature_annual.png` (242.7–306.3 K): a warm equatorial band,
     cold poles; continents cooler than the ocean at their latitude, the
     ranges along their convergent margins as cold lines (the lapse rate),
     the large southern continent (40–70° S) the coldest land. Zonal means
     from 260.6 K (80–90° S) to 300.4 K (0–10° S) to 261.7 K (80–90° N);
     global 288.23 K, land 284.88 K, ocean 289.68 K.
   - `temperature_m01.png` (244.9–306.2 K) and `temperature_m07.png`
     (234.7–306.3 K), and on one range (234.0–306.5 K,
     `scaled/temperature_m01_shared.png`, `_m07_shared.png`): the two
     northern continents (20–55° N) colder than the ocean in January and
     about as warm as it in July; the large southern continent dark in
     July (its winter) and lighter in January; the oceans nearly the same
     in both. Zonal means 40–50° N 279.6 K in January against 284.2 K in
     July, 40–50° S 286.0 K against 280.5 K.
   - `precipitation_annual.png`, `precipitation_m01.png`,
     `precipitation_m07.png` (0 to 165.5, 124.9 and 136.9 mm/day, linear):
     nearly uniformly dark, a faint equatorial band and bright windward
     coasts only. On √(mm/day) over 0–16 mm/day
     (`scaled/precipitation_*_sqrt.png`), as M4-06's maps: the intense
     equatorial band over the oceans and across the equatorial coasts; dry
     subtropical belts at 15–35° over the oceans and across the bodies of
     the continents there; storm tracks at 40–60°, wettest on the west
     coasts and windward ranges; the band north of the equator in July
     (5–15° N, crossing the western continent's north) and south of it over
     the continents in January. Global 0.944 m/yr (land 0.829, ocean
     0.993); zonal means 11.1 mm/day at 0–10° S and 7.9 at 0–10° N,
     0.1–0.5 mm/day at 20–40° in both hemispheres; 33.7 % of the land
     with no rain in the year.
   - `wind_east.png` (−6.24 to 6.28 m/s) and on ±8 m/s
     (`scaled/wind_east_div.png`): westerlies at 35–65° in both
     hemispheres, easterly trades within 30°, polar easterlies, weak
     westerly patches near the equator; `wind_north`
     (`scaled/wind_north_div.png`): equatorward trades (southward at
     10–30° N, northward at 10–30° S), poleward westerlies.
   - `current_east.png`, `current_north.png` (−0.17 to 0.18 m/s) and on
     ±0.1 m/s: eastward drift under the westerlies, westward under the
     trades, the strongest currents along the western boundaries of the
     basins: poleward at 15–35° on the east coasts of the continents
     (subtropical gyres), equatorward at 45–55° (subpolar gyres), 0 on
     land.
   - `rivers.png`: dendritic networks over every continent, as before the
     climate (the rivers follow the uniform 1 m/yr; the fields are those of
     `main` bit for bit, criterion 3).
6. Golden fingerprints (10,000 cells, seed 1): 
   `golden_fingerprints_of_the_climate_stage` (30 fields, FNV-1a of the
   values, debug and release alike) and the CLI's `GOLDEN_FIELDS` (46
   `.kfield` fingerprints; the 30 new, the 16 others unchanged), recorded
   after the tests above passed and the maps were viewed.
7. Timing (release, this machine, `kuniumi generate --seed 1`, 650,000
   cells, from the time stamps of the progress lines): grid 2.9 s,
   tectonics 1.1 s, crust 0.1 s, climate 686.1 s (the energy balance
   240.8 s, the circulation and the currents 40.0 s, the sea surface
   322.9 s, the moisture 82.4 s), erosion 108.9 s, writing the world
   6.3 s; in all 805.5 s, 13.4 minutes. The Climate stage on 10,000 cells
   5.7 s (debug 54 s), on 1,000 cells 0.45 s (debug 6.3 s).
8. `mise run ci` passes after the fixes of review round 1: fmt, clippy
   `-D warnings`, rustdoc, the browser build's lints, 2,005 tests summed
   over debug and release, MSRV 1.85.0, cargo-deny (advisories, bans,
   licences and sources ok; duplicate-version warnings only); 5 min 14 s.
   The debug tests of `kuniumi-sim` take 125 s (103 s on `main`), the
   CLI's 92 s (17 s), `kuniumi-io`'s 11 s (8 s), `kuniumi-viewer`'s 15 s.

Review:
- Round 1 (APPROVE; the reviewer re-ran `mise run ci`, regenerated the
  650,000-cell world, whose 46 `.kfield` files and grid were byte-identical
  to these, timed it within 5 %, reproduced the statistics of the maps,
  checked the stage's golden values against the CLI's world, and ran 27
  orbital-extreme planets and the 128 corners with probes). R1-1 minor: the
  module's "Cost" gave the energy balance's 10 years at 650,000 cells, it
  takes 11: corrected, with the sea surface's 12 and the moisture's 3.
  R1-2 minor: no committed test ran the chain at the bounds of the orbit
  (obliquity, eccentricity, perihelion) or of the continental fraction:
  `the_orbits_at_their_bounds_generate_finite_worlds` added (6 planets).
  R1-3 nit: the Scope's rivers map was neither listed nor described: added
  to "Images and timing" and the Result. R1-4 nit: the field documentation
  now says that the heights and the land are those of the Crust stage's
  relief. `mise run ci` passes again; the loop ended on the approval.

Follow-ups: M4-16 (warnings of an unconverged climate, split with the
owner's approval), M4-17 (the coupled periodic state of temperature, winds
and currents), M4-18 (climate maps on shared scales), M4-19 (erosion from
the climate's precipitation), all draft.

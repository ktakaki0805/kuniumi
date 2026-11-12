# M2-08 Physical bounds for the planet parameters

- Status: done
- Milestone: M2
- Depends on: none
- Crates: kuniumi-world, kuniumi-math (the gravitational constant),
  kuniumi-config (the error of the new check); tests and one doc comment
  of kuniumi-sim; a test of kuniumi-cli
- Design: プロジェクト方針
- Decisions: D-12

## Goal

`PlanetParams::validate` rejects values outside the range of Earth-like
planets that v1 generates, so that no stage meets a radius, gravity or period
for which its arithmetic overflows or underflows.

## Scope

In:
- Lower and upper bounds for `radius_m`, `surface_gravity_m_s2`,
  `rotation_period_s`, `orbital_period_s`, `stellar_flux_w_m2`,
  `ocean_volume_m3` and `surface_pressure_pa` (now any positive finite
  value, or any non-negative one for the ocean volume), chosen with the owner
  for the Earth-like range of D-12 (the table below); `ParamError::OutOfRange`
  names the parameter as now.
- The bounds as public associated constants of `PlanetParams`, like
  `Settings::MIN_CELLS`, each with a doc comment giving its value in SI and
  in the settings file's human unit, and its reason.
- Two checks across parameters, chosen by the owner when the bounds were
  approved: the bulk density that the gravity and the radius imply lies
  in a range, and the rotation is slower than break-up (below).
- The doc comment of each bounded field names its bounds; D-12 and the
  design's 対象惑星 row record the bounds and the two checks.
- Tests that use values outside the new bounds move inside them:
  `kuniumi-config`'s awkward and random planets, and the water world of
  `stages::crust`'s tests (1e21 m³). `stages::crust::sea_level_m`'s note on
  the radius points to the bounds.

Out:
- Settings files (M2-07), which check the same bounds through `validate`.
- Other checks across parameters, such as the ocean's depth relative to
  the radius (1e20 m³ on a 1,500 km planet is several times the planet's
  volume, but finite): draft M2-17.
- The bounds of `obliquity_rad`, `orbital_eccentricity` and
  `continental_crust_fraction`, unchanged.
- The allowed range in the error message (`ParamError::OutOfRange` keeps its
  fields).

## Implementation notes

Found by the review of M2-02 (R2-1): with `radius_m = 5e-324` the Tectonics
stage's `radius_m * lever` underflows to 0, so the angular speed and every
velocity are infinite (caught only by the runner's debug assertion); with
`radius_m = f64::MAX` the angular speed is subnormal and the fastest speed is
off by about 1e-7 relative. Radii from 1e-310 m to 1e300 m were fine
through the Tectonics stage. The Crust stage (M2-05) narrows this: its sea
level divides the ocean volume by the planet's area `4πR²`, so with Earth's
water a radius below about 1e-145 m gives an infinite sea level (and `R²`
underflows to 0 below about 1.5e-154 m), caught only by the runner's debug
assertion; measured in the review of M2-05 (R1-2), 1e-100 m still gave a
finite 1.06e217 m. `stages::crust::sea_level_m` documents a radius of 1 m
or more.

### Bounds

Proposed in the refinement and approved by the owner, with the checks
across parameters below.

Each bound is a round number in the settings file's human unit
(`radius_km`, `rotation_period_h`, `orbital_period_days`,
`ocean_volume_km3`, `surface_pressure_bar`) whose conversion to SI is exact,
so that writing the bound in either key gives the same value.

| Field | Earth | Lower | Upper |
| --- | --- | --- | --- |
| `radius_m` | 6,371 km | 1.5e6 m (1,500 km, 0.24 R⊕) | 1.3e7 m (13,000 km, 2.04 R⊕) |
| `surface_gravity_m_s2` | 9.81 | 1.0 m/s² (0.10 g⊕) | 40.0 m/s² (4.1 g⊕) |
| `rotation_period_s` | 23.93 h | 28,800 s (8 h) | 864,000 s (240 h, 10 days) |
| `orbital_period_s` | 365.26 days | 7,776,000 s (90 days) | 315,360,000 s (3,650 days) |
| `stellar_flux_w_m2` | 1,361 | 400 W/m² (0.29 S⊕) | 2,500 W/m² (1.84 S⊕) |
| `ocean_volume_m3` | 1.335e9 km³ | 0 (a dry planet; −0 too) | 1e20 m³ (1e11 km³, 75 Earth oceans) |
| `surface_pressure_pa` | 1.013 bar | 1e4 Pa (0.1 bar) | 1e6 Pa (10 bar) |

Reasons:

- Radius: the Moon (1,737 km), Mercury (2,440 km) and Mars (3,390 km) lie
  inside. Above 1.5 to 2 Earth radii, planets mostly have thick volatile
  envelopes rather than a rocky surface: the radius gap of Fulton et al.
  (2017, AJ 154, 109) and Rogers (2015, ApJ 801, 41).
- Gravity: the Moon (1.62 m/s²) and Mars (3.72 m/s²) lie inside; an
  Earth-composition planet of 13,000 km has about 33 m/s² (mass–radius
  relation of Zeng et al. 2016, ApJ 819, 127, `R = 1.00 M^(1/3.7)` for
  Earth's core mass fraction).
- Rotation, lower: three times Earth's rate. At Earth's density the
  equatorial centrifugal acceleration is then 3.1 % of gravity (Earth:
  0.35 %), so the spherical grid stays a fair approximation; break-up is at
  1.4 h.
- Rotation, upper: a tenth of Earth's rate. Slower rotators approach the
  day–night circulation of tidally locked planets (v2), and with the
  shortest year (90 days) a year still holds at least 8 solar days, so the
  daily-mean insolation of the climate model (design, 気候モデル) applies.
- Orbit, lower: at Earth's flux a 90-day year belongs to a star of about
  0.57 solar masses (`L ∝ M⁴`, so `P = M^2.5` years); planets in the
  habitable zones of smaller stars are expected to be tidally locked
  (Kasting et al. 1993, Icarus 101, 108), the v2 regime.
- Orbit, upper: at the lowest flux the habitable zone of an F0 star
  (6 L☉, 1.6 M☉) lies at 4.5 AU, a year of about 8 years.
- Flux: the optimistic habitable zone of a Sun-like star spans 0.32 S⊕
  (early Mars) to 1.78 S⊕ (recent Venus) (Kopparapu et al. 2013, ApJ 765,
  131), widened to round numbers. Mars receives 586 W/m²; Venus, 2,601 W/m²,
  lies outside.
- Ocean: 0 is a dry planet; 1e20 m³ drowns every planet of the box (a layer
  47 km deep on a 13,000 km planet, whose relief from the deepest trench to
  the highest mountain is about 20 km).
- Pressure: below about 0.1 bar water boils below 46 °C, so an ocean cannot
  be temperate (the triple point is 612 Pa); above 10 bar the atmosphere is
  far from the Earth calibration of the climate model (M4-08).

### Checks across parameters

The owner's answer to the bounds above: "提案の範囲＋組み合わせの検査", that
is, the bounds plus a check of the bulk density from gravity and radius and
of the rotation against break-up.

- **Bulk density.** `ρ = 3 g / (4π G R)` (kg/m³), from
  `g = G M / R²` and `M = (4/3) π R³ ρ`, with the gravitational constant
  `G = 6.674 30e-11 m³ kg⁻¹ s⁻²` (CODATA 2018), must lie in
  `[MIN_BULK_DENSITY_KG_M3, MAX_BULK_DENSITY_KG_M3]`. Earth: 5,506 kg/m³
  (5,513 from its mass: standard gravity, 9.80665 m/s², lies below
  `GM/R²` = 9.82 m/s² because it includes the centrifugal acceleration).
  Without it the box holds 275 kg/m³ (13,000 km, 1 m/s²) and 95,000 kg/m³
  (1,500 km, 40 m/s²).

  Bounds, 3,000 to 16,000 kg/m³ (the owner's answer: "3,000〜16,000
  kg/m³ (Recommended)", with break-up checked at compile time):
  - Lower: rocky bodies of little self-compression are about as dense as
    silicate rock: the Moon 3,344 kg/m³, Io 3,528, Mars 3,934. With it the
    smallest planet (1,500 km) needs at least 1.26 m/s².
  - Upper: the densest rocky planet of the box, pure iron at 13,000 km,
    about 15,900 kg/m³ (Zeng et al. 2016, `R = 0.86 M^(1/3.7)` for a core
    mass fraction of 1; an Earth-like composition gives 9,000 kg/m³ there,
    5,500 at Earth's radius). A constant bound rejects no possible rocky
    planet of the box, but admits small planets denser than iron could be
    at their size (uncompressed iron is 7,870 kg/m³); at 13,000 km the
    gravity bound already caps the density at 11,000 kg/m³.
- **Break-up.** The rotation must be slower than break-up: the centrifugal
  acceleration at the equator below gravity, `ω² R < g` with
  `ω = 2π / P_rot`, that is `P_rot > 2π √(R / g)`. The bounds already
  ensure it: the slowest break-up of the box (13,000 km, 1 m/s²) is 6.29 h,
  shorter than the 8 h lower bound, and with the density check break-up is
  at most 1.9 h (`P = √(3π / (G ρ))` at 3,000 kg/m³). A run-time check
  could never fail, so it is a compile-time assertion over the constants,
  `4π² · MAX_RADIUS_M < MIN_SURFACE_GRAVITY_M_S2 · MIN_ROTATION_PERIOD_S²`
  (const float arithmetic, stable since Rust 1.82, below the MSRV), plus a
  test; changing a bound so that break-up becomes possible stops the
  build.

### API

```rust
impl PlanetParams {
    pub const MIN_RADIUS_M: f64 = 1.5e6;
    pub const MAX_RADIUS_M: f64 = 1.3e7;
    pub const MIN_SURFACE_GRAVITY_M_S2: f64 = 1.0;
    pub const MAX_SURFACE_GRAVITY_M_S2: f64 = 40.0;
    pub const MIN_ROTATION_PERIOD_S: f64 = 28_800.0;
    pub const MAX_ROTATION_PERIOD_S: f64 = 864_000.0;
    pub const MIN_ORBITAL_PERIOD_S: f64 = 7_776_000.0;
    pub const MAX_ORBITAL_PERIOD_S: f64 = 315_360_000.0;
    pub const MIN_STELLAR_FLUX_W_M2: f64 = 400.0;
    pub const MAX_STELLAR_FLUX_W_M2: f64 = 2_500.0;
    pub const MAX_OCEAN_VOLUME_M3: f64 = 1e20; // the lower bound is 0
    pub const MIN_SURFACE_PRESSURE_PA: f64 = 1e4;
    pub const MAX_SURFACE_PRESSURE_PA: f64 = 1e6;
    pub const MIN_BULK_DENSITY_KG_M3: f64 = 3_000.0;
    pub const MAX_BULK_DENSITY_KG_M3: f64 = 16_000.0;

    /// 3 g / (4π G R) (kg/m³).
    pub fn bulk_density_kg_m3(&self) -> f64;
}

// kuniumi-math
pub mod constants {
    /// Newtonian constant of gravitation (m³ kg⁻¹ s⁻²), CODATA 2018.
    pub const GRAVITATIONAL_CONSTANT_M3_KG_S2: f64 = 6.674_30e-11;
}

pub enum ParamError {
    OutOfRange { .. },   // unchanged
    /// Parameters, each within its bounds, that together describe an
    /// impossible planet.
    Inconsistent {
        /// The derived quantity, with its unit: "bulk_density_kg_m3".
        name: &'static str,
        /// Its value.
        value: f64,
        /// The fields it comes from, in field order:
        /// ["radius_m", "surface_gravity_m_s2"].
        fields: [&'static str; 2],
    },
    Unsupported { .. },  // unchanged
}

// kuniumi-config
pub enum ConfigError {
    // …
    /// Keys whose values, each in range, do not fit together.
    Inconsistent {
        /// The keys as written, with their table, in field order
        /// (["planet.radius_km", "planet.surface_gravity_m_s2"]).
        keys: [&'static str; 2],
        /// What PlanetParams::validate reported.
        error: ParamError,
    },
}
```

`ParamError::Inconsistent` displays as "bulk_density_kg_m3 from radius_m
and surface_gravity_m_s2 is out of range: 275", in the form of
`OutOfRange` ("… is out of range: v"); `ConfigError` names the keys as
written instead of the fields.

`validate` checks every field in field order with the existing `within`
helper (`v.is_finite() && lo <= v && v <= hi`, inclusive; NaN and ±∞ fail),
and returns the first failure as `OutOfRange { name, value }`, as now;
then the bulk density (`Inconsistent`); then the regimes, as now. The ocean
volume's lower bound stays 0, so −0 passes as before. The `positive` helper
goes.

### Tests

- `kuniumi-world` (`params::tests`): for each of the seven fields, the
  lower and upper bound pass; the next double below the lower bound and
  above the upper bound (from the bits, since `f64::next_up` needs Rust
  1.86, above the MSRV), NaN, +∞ and −∞ fail as `OutOfRange` naming that
  field, with every other field Earth's. The bounds hold Earth's value, and
  the human-unit bounds of the table convert to the constants exactly
  (`kuniumi_math::units`). The bulk density: Earth's is 5,506 kg/m³
  (hand calculation above, to 1e-3 relative); gravities giving exactly the
  density bounds pass at the smallest, Earth's and the largest radius, the
  next double beyond fails as `Inconsistent` naming both fields; a
  per-field failure is reported before the density. Break-up: the
  compile-time assertion, and a test that the slowest break-up period of
  the box, `2π √(MAX_RADIUS_M / MIN_SURFACE_GRAVITY_M_S2)`, is below
  `MIN_ROTATION_PERIOD_S`.
- `kuniumi-sim`: the Earth parameters and the extremes of the valid
  region generate through the Crust stage (`allow_partial`, 1,000 cells,
  one fixed seed): each of the 2⁷ = 128 corners of the box (every bounded
  field at its lower or upper bound) with the gravity moved, where the
  corner's density is out of range, to the nearest gravity whose density
  is in range. The sea level and every field value are finite, checked
  explicitly so that release builds check it too. If the 128 runs take
  more than a few seconds in debug, the corners of the two fields the
  stages read (radius, ocean volume) plus the all-lower and all-upper
  corners are enough, and the Result says so.
- `kuniumi-config`: the awkward planets keep subnormals, −0, 17-digit values
  and exponents where the bounds allow them (ocean volume, obliquity,
  eccentricity, continental fraction) and put the bounded fields at their
  bounds and at awkward values inside; the random planets draw inside the
  bounds and the density range. A file with each bound in its human-unit
  key loads to the constant bit for bit; one just outside fails naming its
  key; a file whose radius and gravity give a density out of range fails
  with `ConfigError::Inconsistent` naming both keys as written
  (`radius_km` when the file used it).
- The existing tests (`rejects_negative_radius`,
  `generate_validates_params_first`, the NaN gravity of the pipeline tests)
  keep passing unchanged.

No golden value changes: Earth and the other test planets lie inside the
bounds, and no stage changes. Parameter sets that passed before and lie
outside the bounds are now rejected; the changelog lists this under
*Changed*. `World-breaking: no`.

## Acceptance criteria

- [x] Each bounded parameter: its bounds accepted, values just outside and
      non-finite values rejected by name
      (`params::tests` in `kuniumi-world`).
- [x] A radius and gravity whose bulk density `3g / (4πGR)` is outside
      the density bounds are rejected as `Inconsistent` naming both
      fields; the density bounds themselves pass (`params::tests`).
- [x] Break-up cannot happen inside the bounds: a compile-time assertion
      over the constants and a test.
- [x] The human-unit bounds of the table load from a settings file to the
      SI constants bit for bit, a value just outside fails naming its key,
      and an inconsistent radius and gravity fail naming both keys
      (`kuniumi-config` tests).
- [x] The Earth parameters and the extremes of the valid region (the 128
      corners of the box, gravity moved into the density range) generate
      finite worlds through the Crust stage (`kuniumi-sim` test, finite
      sea level and field values).
- [x] No golden value changes; D-12 and the design's 対象惑星 row list the
      bounds and the two checks.
- [x] `mise run ci` passes.

## Verification

```sh
mise exec -- cargo test -p kuniumi-world
mise exec -- cargo test -p kuniumi-config
mise exec -- cargo test -p kuniumi-sim bounds
mise run ci
```

## Result

What was done:

- `PlanetParams` has the bounds as `MIN_`/`MAX_` constants (each doc
  comment gives the value in SI and in the human unit, and its reason),
  `MIN_BULK_DENSITY_KG_M3` = 3,000, `MAX_BULK_DENSITY_KG_M3` = 16,000 and
  `bulk_density_kg_m3()` = `3g / (4πGR)`. `validate` checks every float
  field within its bounds in field order (`OutOfRange`), then the density
  (the new `ParamError::Inconsistent { name: "bulk_density_kg_m3", value,
  fields: ["radius_m", "surface_gravity_m_s2"] }`), then the regimes. A
  `const` assertion, `4π² · MAX_RADIUS_M < MIN_SURFACE_GRAVITY_M_S2 ·
  MIN_ROTATION_PERIOD_S²` (5.1e8 < 8.3e8), keeps every planet within the
  bounds slower than break-up.
- `kuniumi_math::constants::GRAVITATIONAL_CONSTANT_M3_KG_S2` =
  6.674 30e-11 (CODATA 2018).
- `kuniumi-config`: `ConfigError::Inconsistent { keys, error }` names both
  keys as written (the SI key for a value left at Earth's), displayed as
  "`planet.radius_km` and `planet.surface_gravity_m_s2` do not fit
  together: bulk_density_kg_m3 = 2698.254998805064 is out of range".
- D-12, the design's プロジェクト方針 (a paragraph under the table, pointed
  to from the 対象惑星 row) and the changelog (*Changed*) list the bounds.

Deviations from the plan, and why:

- The gravity's lower bound, 1 m/s², never decides alone: with the density
  bound the lowest gravity of a valid planet is 1.258 m/s² (1,500 km at
  3,000 kg/m³). The bound stays as approved (it names the gravity for
  values far below); the test checks that 1 m/s² passes its own check and
  fails as `Inconsistent` at the smallest radius, where the density is
  highest (2,385 kg/m³).
- Tests that used values outside the bounds: `kuniumi-config`'s Mars-like
  file had Mars's 0.006 bar (now 0.7 bar, with Mars's gravity, 3.721 m/s²,
  added), the largest TOML integer moved from `radius_m` to
  `ocean_volume_m3`; the Crust stage's gravity test used Mars's gravity at
  Earth's radius (2,089 kg/m³), now twice Earth's gravity (11,012 kg/m³);
  its water world 1e21 m³ is now `MAX_OCEAN_VOLUME_M3`. Each of those
  planets is now asserted valid. One CLI case was added for an
  inconsistent file.
- Size: about 1,400 changed lines, most of them tests and this ticket;
  the checks across parameters (the owner's choice in the refinement) add
  the new error variants and their tests.

Evidence:

- Criterion 1, the field bounds: `params::tests::each_bounded_parameter_accepts_its_bounds_and_rejects_the_values_beyond`
  (both bounds pass, −0 for the ocean; the next double beyond each bound,
  NaN, +∞ and −∞ fail as `OutOfRange` naming the field, with the value's
  bits); `the_bounds_hold_earth_and_are_round_in_human_units` (Earth inside
  every bound and the density range; the nine bounds with a human unit
  equal `units::km_to_m(1_500.0)` and so on bit for bit).
- Criterion 2, the density: `the_bulk_density_comes_from_gravity_and_radius`
  (Earth 5,505.78 kg/m³ against Python within 1e-12, and within 0.2 % of
  5,513 kg/m³ from Earth's mass); `a_density_out_of_range_is_inconsistent`
  (the extreme gravities at each density bound at 1,500 km, 6,371 km and,
  for 3,000 kg/m³, 13,000 km pass; the next double beyond fails as
  `Inconsistent` naming both fields; at 13,000 km the gravity bound comes
  first, 40 m/s² being 11,006 kg/m³);
  `checks_run_fields_first_then_the_density_then_the_regimes`;
  `inconsistent_parameters_are_displayed_with_their_fields`.
- Criterion 3, break-up: the `const` assertion in `params.rs`, and
  `every_planet_within_the_bounds_rotates_slower_than_break_up` (slowest
  break-up 22,654.3 s = 6.29 h < 28,800 s, against Python; at 3,000 kg/m³
  break-up is 6,861 s).
- Criterion 4, settings files:
  `the_bounds_in_human_units_load_to_the_si_bounds` (12 files with a bound
  under its human-unit key, or the SI key for gravity and flux, load to
  the constant bit for bit; 13 files just beyond fail as `Invalid` naming
  the key as written, with the SI value);
  `a_radius_and_gravity_that_do_not_fit_name_both_keys` (three files, the
  keys as written, the density against Python within 1e-12, the full
  message; a field out of its bounds is reported first);
  `generate_reports_a_bad_planet_file_and_writes_nothing` in the CLI tests.
  The awkward and 10,000 random planets of
  `random_and_awkward_planets_round_trip_bit_exactly` lie inside the
  bounds and the density range.
- Criterion 5, finite worlds: `kuniumi-sim`
  `tests::the_bounds_of_the_planet_parameters_generate_finite_worlds`:
  Earth and all 128 corners of the box (gravity moved just inside the
  density range where needed: 1 → 10.9 m/s² at 13,000 km, 1 → 1.26 and
  40 → 6.71 m/s² at 1,500 km), 1,000 cells, seed 3, through the Crust
  stage; finite sea level and field values checked explicitly. 2.45 s in
  debug, so all 128 corners are kept. Mutation: with `MIN_RADIUS_M` =
  1e-150 m and the density bounds opened, the release run fails at the
  infinite sea level.
- Criterion 6: no golden value or fingerprint changed (`git diff main`
  touches none); D-12 and the design list the bounds and both checks.
- Criterion 7: `mise run ci` passes (921 tests over debug and release,
  clippy, MSRV 1.85.0, cargo-deny ok), before and after the review fixes.

Extreme planets through the CLI (20,000 cells, seed 1,
`out/M2-08/*.toml`):

| Planet | Sea level | Land |
| --- | --- | --- |
| 1,500 km, 2.5 m/s², no water | −8,000 m (the lowest height) | 100 % |
| 1,500 km, 2.5 m/s², 1.6e7 km³ | −4,110 m | 53.3 % |
| 13,000 km, 20 m/s², Earth's water | −3,979 m | 51.3 % |
| 13,000 km, 40 m/s², 1e11 km³, every other field at a bound | +44,105 m | 0 % (mean depth 47,087 m) |

Images looked at: `out/M2-08/large_max/surface_height.png` (plates,
continents and trench lines, no no-data colour; flat continental
interiors, since the 600 km margins are narrow on a 13,000 km planet),
`out/M2-08/large_max/hypsometry.png` (every bar far left of the sea-level
line at +44 km: a water world, as the 47 km layer of the hand estimate
says), `out/M2-08/small_dry/surface_height.png` (the same plate layout,
the grid being on the unit sphere, with wide margin gradients and
mountain belts on the small planet), and
`out/M2-08/large_earth_water/hypsometry.png` (two modes with the sea level
between them).

Review: 1 round.

- Round 1 (APPROVE): the reviewer re-ran `mise run ci` (and clippy and the
  MSRV build in clean target directories), recomputed every reference
  value and the numbers of the doc comments in Python, mutated
  `MIN_ROTATION_PERIOD_S` to 20,000 s (the build fails with the
  assertion's message on 1.97.0 and 1.85.0) and the radius bound (the
  corner test fails in release), and regenerated the extreme planets.
  R1-1 minor, the ticket's API text gave another form of the
  `Inconsistent` message (fixed: the text now gives the implemented form,
  that of `OutOfRange`). R1-2 minor, 1e20 m³ on a 1,500 km planet gives a
  thin-layer sea level of +3,533 km where a spherical shell would be
  1,509 km deep, and no ticket recorded it (fixed: draft M2-17, and the
  `MAX_OCEAN_VOLUME_M3` doc names the smallest planet). R1-3 nit, Earth's
  density from its mass is 5,513 kg/m³, not 5,514 (fixed). R1-4 nit, an
  unclear sentence on the gravity's lower bound (fixed). R1-5 nit, the
  `Crates:` line missed a doc comment of `kuniumi-sim` (fixed).

Follow-ups: M2-17 (draft, sea level of deep water worlds).

//! Tests of the Erosion stage in the pipeline (ticket M3-05):
//!
//! - `the_stage_erodes_and_routes_the_final_surface`: [`run`] equals
//!   [`erode`] of the Crust stage's surface roughened by [`roughen`] with
//!   the stream `erosion/roughness` of the world's seed and an amplitude of
//!   1 m (ticket M3-06), followed by the routing of the final surface
//!   composed by hand from the drainage primitives;
//! - `the_rivers_conserve_area_and_water`: drainage area and discharge of
//!   the base levels against the planet's area and rain, discharge against
//!   the rate times the area, depression depths;
//! - `the_stage_uses_every_erosion_setting` (review round 1, R1-1);
//! - `zero_steps_keep_the_surface_of_the_crust_stage` (without steps, no
//!   roughness either);
//! - `existing_outputs_leave_the_world_unchanged` and the panics, of the
//!   stage and of [`Rivers::new`];
//! - `golden_fingerprints_of_the_erosion_stage`.
#![allow(clippy::disallowed_methods)]

use core::f64::consts::PI;

use kuniumi_math::Rng;
use kuniumi_math::rng::fnv1a64;
use kuniumi_math::units::SECONDS_PER_YEAR;
use kuniumi_world::{ErosionSettings, FieldError, FieldKind, World};

use super::tests::{dry_world, worlds};
use super::*;
use crate::drainage::depression_depth_m;
use crate::stages::crust::SURFACE_HEIGHT;
use crate::stages::tectonics::UPLIFT_RATE;

/// `world`, a world of the Crust stage, with `steps` steps of the default
/// erosion settings.
fn with_steps(world: &World, steps: u32) -> World {
    let mut world = world.clone();
    world.settings.erosion = ErosionSettings {
        steps,
        ..ErosionSettings::default()
    };
    world
}

/// The values of the field `name` of `world`.
fn values<'a>(world: &'a World, name: &str) -> &'a [f64] {
    &world
        .field(name)
        .unwrap_or_else(|| panic!("no field `{name}`"))
        .values
}

/// The bits of `values`, to compare them exactly.
fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|v| v.to_bits()).collect()
}

/// Runs the stage on `world`, returning its progress reports.
fn run_stage(world: &mut World) -> Vec<f64> {
    let mut reports = Vec::new();
    run(world, &mut |f| reports.push(f)).expect("the stage runs");
    reports
}

/// What the stage should produce from the world `before` of the Crust
/// stage, composed by hand from the primitives: the eroded heights (of the
/// surface roughened by 1 m from the stream `erosion/roughness` of the
/// world's seed, when there are steps), the graph of the final surface
/// routed across its depressions, the filled surface, and the three
/// outputs.
struct Expected {
    height_m: Vec<f64>,
    drainage: Drainage,
    filled_m: Vec<f64>,
    drainage_area_m2: Vec<f64>,
    discharge_m3_s: Vec<f64>,
    depression_depth_m: Vec<f64>,
}

fn expected(before: &World, settings: &ErosionSettings) -> Expected {
    let grid = &before.grid;
    let (r, z) = (before.params.radius_m, before.sea_level_m);
    let rain = vec![PLACEHOLDER_PRECIPITATION_M_S; grid.len()];
    let inputs = ErosionInputs {
        grid,
        radius_m: r,
        sea_level_m: z,
        uplift_rate_m_s: values(before, UPLIFT_RATE),
        precipitation_m_s: &rain,
    };
    let mut height_m = values(before, SURFACE_HEIGHT).to_vec();
    if settings.steps > 0 {
        let mut rng = Rng::for_stream(before.seed, "erosion/roughness");
        roughen(&mut height_m, z, 1.0, &mut rng);
    }
    erode(&inputs, settings, &mut height_m, &mut |_| {});
    let filled_m = fill_depressions(grid, &height_m, z);
    let drainage = Drainage::build(grid, &filled_m, z, r);
    let area: Vec<f64> = (0..grid.len() as u32)
        .map(|i| grid.area_sr(i) * r * r)
        .collect();
    let water: Vec<f64> = (0..grid.len() as u32)
        .map(|i| PLACEHOLDER_PRECIPITATION_M_S * grid.area_sr(i) * r * r)
        .collect();
    Expected {
        drainage_area_m2: drainage.accumulate(&area),
        discharge_m3_s: drainage.accumulate(&water),
        depression_depth_m: depression_depth_m(&height_m, &filled_m),
        height_m,
        drainage,
        filled_m,
    }
}

/// `Σ values` with Neumaier's compensated summation.
fn compensated_sum(values: impl IntoIterator<Item = f64>) -> f64 {
    let (mut sum, mut c) = (0.0_f64, 0.0_f64);
    for v in values {
        let t = sum + v;
        c += if sum.abs() >= v.abs() {
            (sum - t) + v
        } else {
            (v - t) + sum
        };
        sum = t;
    }
    sum + c
}

/// The Earth worlds of seeds 1–3 and the dry Earth world of seed 1, at
/// 10,000 cells, through the Crust stage.
fn crust_worlds() -> Vec<&'static World> {
    worlds().iter().chain([dry_world()]).collect()
}

#[test]
fn the_stage_erodes_and_routes_the_final_surface() {
    for before in crust_worlds() {
        let mut world = with_steps(before, 3);
        let reports = run_stage(&mut world);
        assert_eq!(reports, [1.0 / 3.0, 2.0 / 3.0, 1.0], "seed {}", world.seed);
        let want = expected(before, &world.settings.erosion);
        assert!(
            bits(values(&world, SURFACE_HEIGHT)) == bits(&want.height_m),
            "seed {}",
            world.seed
        );
        for (name, expected) in [
            (DRAINAGE_AREA, &want.drainage_area_m2),
            (DISCHARGE, &want.discharge_m3_s),
            (DEPRESSION_DEPTH, &want.depression_depth_m),
        ] {
            assert!(
                bits(values(&world, name)) == bits(expected),
                "{name}, seed {}",
                world.seed
            );
        }
        // The heights were eroded: the stage did something.
        assert!(
            bits(values(&world, SURFACE_HEIGHT)) != bits(values(before, SURFACE_HEIGHT)),
            "seed {}",
            world.seed
        );
        // And roughened first: not the erosion of the bare surface.
        let mut bare = values(before, SURFACE_HEIGHT).to_vec();
        let rain = vec![PLACEHOLDER_PRECIPITATION_M_S; bare.len()];
        let inputs = ErosionInputs {
            grid: &before.grid,
            radius_m: before.params.radius_m,
            sea_level_m: before.sea_level_m,
            uplift_rate_m_s: values(before, UPLIFT_RATE),
            precipitation_m_s: &rain,
        };
        erode(&inputs, &world.settings.erosion, &mut bare, &mut |_| {});
        assert!(
            bits(values(&world, SURFACE_HEIGHT)) != bits(&bare),
            "seed {}",
            world.seed
        );
        // The three outputs come after the fields of the earlier stages,
        // which keep their order and, but for the heights, their values;
        // the sea level is unchanged.
        let fields: Vec<_> = world
            .fields
            .iter()
            .map(|f| (f.name.as_str(), f.unit.as_str(), f.kind))
            .collect();
        let mut expected_fields: Vec<_> = before
            .fields
            .iter()
            .map(|f| (f.name.as_str(), f.unit.as_str(), f.kind))
            .collect();
        expected_fields.extend([
            (DRAINAGE_AREA, "m^2", FieldKind::Quantity),
            (DISCHARGE, "m^3/s", FieldKind::Quantity),
            (DEPRESSION_DEPTH, "m", FieldKind::Quantity),
        ]);
        assert_eq!(fields, expected_fields, "seed {}", world.seed);
        for (a, b) in before.fields.iter().zip(&world.fields) {
            if a.name != SURFACE_HEIGHT {
                assert!(bits(&a.values) == bits(&b.values), "{}", a.name);
            }
        }
        assert_eq!(world.sea_level_m.to_bits(), before.sea_level_m.to_bits());
        assert_eq!(world.stages, before.stages, "the runner records the stage");
    }
}

#[test]
fn the_rivers_conserve_area_and_water() {
    for before in crust_worlds() {
        let mut world = with_steps(before, 10);
        run_stage(&mut world);
        let want = expected(before, &world.settings.erosion);
        let (grid, r, z) = (&world.grid, world.params.radius_m, world.sea_level_m);
        let p = PLACEHOLDER_PRECIPITATION_M_S;
        let height = values(&world, SURFACE_HEIGHT);
        let area = values(&world, DRAINAGE_AREA);
        let discharge = values(&world, DISCHARGE);
        let depth = values(&world, DEPRESSION_DEPTH);
        let base = want.drainage.base_levels();
        let planet_m2 = 4.0 * PI * r * r;
        let total_area = compensated_sum(base.iter().map(|&i| area[i as usize]));
        let total_water = compensated_sum(base.iter().map(|&i| discharge[i as usize]));
        let relative = |x: f64, y: f64| (x - y).abs() / y;
        assert!(
            relative(total_area, planet_m2) <= 1e-12,
            "seed {}: {total_area:e} against {planet_m2:e}",
            world.seed
        );
        assert!(
            relative(total_water, p * planet_m2) <= 1e-12,
            "seed {}: {total_water:e} against {:e}",
            world.seed,
            p * planet_m2
        );
        let mut deepest = 0.0_f64;
        for i in 0..grid.len() {
            let own = grid.area_sr(i as u32) * r * r;
            assert!(area[i] >= own, "seed {}, cell {i}", world.seed);
            assert!(
                relative(discharge[i], p * area[i]) <= 1e-12,
                "seed {}, cell {i}: {} against {}",
                world.seed,
                discharge[i],
                p * area[i]
            );
            assert!(depth[i] >= 0.0, "seed {}, cell {i}", world.seed);
            if height[i] <= z {
                assert_eq!(depth[i].to_bits(), 0, "seed {}, cell {i}", world.seed);
            }
            deepest = deepest.max(depth[i]);
        }
        // Routed across the depressions: no pit. (The dry world has cells
        // exactly at its sea level, the trenches, which are outlets.)
        assert!(want.drainage.pits().is_empty(), "seed {}", world.seed);
        assert!(height.iter().any(|&h| h <= z), "seed {}", world.seed);
        let ocean = before.params.ocean_volume_m3 > 0.0;
        // The filled surface lies nowhere below the surface.
        assert!(
            want.filled_m.iter().zip(height).all(|(f, h)| f >= h),
            "seed {}",
            world.seed
        );
        println!(
            "seed {} (ocean: {ocean}): {} base levels, area {:.15e} of {:.15e} m², deepest \
             depression {deepest:.3} m",
            world.seed,
            base.len(),
            total_area,
            planet_m2
        );
    }
}

/// A change of one erosion setting.
type Change = fn(&mut ErosionSettings);

#[test]
fn the_stage_uses_every_erosion_setting() {
    // Each of the four rates of `world.settings.erosion` changed alone from
    // its default (K 2e-9 instead of 7.5e-10, m 0.5 instead of 0.45, K_d
    // 0.5 m²/yr instead of 5e-3, Δt 2 × 10⁴ yr instead of 10⁴): the stage
    // erodes as `erode` does with those settings, bit for bit, and not as
    // with the defaults.
    let before = &worlds()[0];
    let defaults = ErosionSettings {
        steps: 3,
        ..ErosionSettings::default()
    };
    let with_defaults = expected(before, &defaults);
    let changes: [(&str, Change); 4] = [
        ("erodibility", |s| s.erodibility = 2.0e-9),
        ("discharge_exponent", |s| s.discharge_exponent = 0.5),
        ("diffusivity_m2_s", |s| {
            s.diffusivity_m2_s = 0.5 / SECONDS_PER_YEAR;
        }),
        ("time_step_s", |s| s.time_step_s = 2.0e4 * SECONDS_PER_YEAR),
    ];
    for (name, change) in changes {
        let mut settings = defaults;
        change(&mut settings);
        assert_eq!(settings.validate(), Ok(()), "{name}");
        let mut world = before.clone();
        world.settings.erosion = settings;
        run_stage(&mut world);
        let want = expected(before, &settings);
        let height = values(&world, SURFACE_HEIGHT);
        assert!(bits(height) == bits(&want.height_m), "{name}");
        assert!(bits(height) != bits(&with_defaults.height_m), "{name}");
        assert!(
            bits(values(&world, DRAINAGE_AREA)) == bits(&want.drainage_area_m2),
            "{name}"
        );
        assert!(
            bits(values(&world, DISCHARGE)) == bits(&want.discharge_m3_s),
            "{name}"
        );
    }
}

#[test]
fn zero_steps_keep_the_surface_of_the_crust_stage() {
    let before = &worlds()[0];
    let mut world = with_steps(before, 0);
    assert!(run_stage(&mut world).is_empty(), "no step, no report");
    assert!(bits(values(&world, SURFACE_HEIGHT)) == bits(values(before, SURFACE_HEIGHT)));
    let want = expected(before, &world.settings.erosion);
    assert!(bits(values(&world, DRAINAGE_AREA)) == bits(&want.drainage_area_m2));
    assert!(bits(values(&world, DISCHARGE)) == bits(&want.discharge_m3_s));
    assert!(bits(values(&world, DEPRESSION_DEPTH)) == bits(&want.depression_depth_m));
}

#[test]
fn existing_outputs_leave_the_world_unchanged() {
    for name in [DRAINAGE_AREA, DISCHARGE, DEPRESSION_DEPTH] {
        let mut world = with_steps(&worlds()[0], 2);
        let n = world.grid.len();
        world
            .add_field(name, "m", vec![0.0; n])
            .expect("a new field");
        let before = world.clone();
        let result = run(&mut world, &mut |f| panic!("progress {f} before the check"));
        assert_eq!(
            result,
            Err(SimError::Field(FieldError::Duplicate { name: name.into() }))
        );
        assert_eq!(world, before, "{name}");
    }
}

#[test]
#[should_panic(expected = "the Erosion stage needs the field `uplift_rate` of the Tectonics stage")]
fn the_stage_needs_the_uplift() {
    let mut world = with_steps(&worlds()[0], 1);
    world.fields.retain(|f| f.name != UPLIFT_RATE);
    let _ = run(&mut world, &mut |_| {});
}

#[test]
#[should_panic(expected = "the Erosion stage needs the field `surface_height` of the Crust stage")]
fn the_stage_needs_the_heights() {
    let mut world = with_steps(&worlds()[0], 1);
    world.fields.retain(|f| f.name != SURFACE_HEIGHT);
    let _ = run(&mut world, &mut |_| {});
}

/// [`Rivers::new`] of the world of seed 1 with the precipitation `rain`.
fn rivers_with_rain(rain: &[f64]) -> Rivers {
    let world = &worlds()[0];
    let inputs = ErosionInputs {
        grid: &world.grid,
        radius_m: world.params.radius_m,
        sea_level_m: world.sea_level_m,
        uplift_rate_m_s: values(world, UPLIFT_RATE),
        precipitation_m_s: rain,
    };
    Rivers::new(&inputs, values(world, SURFACE_HEIGHT))
}

#[test]
#[should_panic(
    expected = "the rivers need one precipitation rate per cell: 3 values for 10000 cells"
)]
fn rivers_need_one_rate_per_cell() {
    let _ = rivers_with_rain(&[0.0; 3]);
}

#[test]
fn rivers_need_finite_non_negative_rates() {
    for bad in [-1.0, -5e-324, f64::NAN, f64::INFINITY] {
        let mut rain = vec![PLACEHOLDER_PRECIPITATION_M_S; worlds()[0].grid.len()];
        rain[17] = bad;
        let message = std::panic::catch_unwind(|| rivers_with_rain(&rain))
            .expect_err("a bad rate panics")
            .downcast::<String>()
            .expect("a formatted message");
        assert_eq!(
            *message,
            format!("the rivers need finite, non-negative precipitation: cell 17 has {bad}")
        );
    }
    // Zero and -0 are fine: no rain.
    let mut rain = vec![0.0; worlds()[0].grid.len()];
    rain[3] = -0.0;
    let rivers = rivers_with_rain(&rain);
    assert!(rivers.discharge_m3_s.iter().all(|&q| q == 0.0));
}

// ---------------------------------------------------------------------------
// Golden values.

#[test]
fn golden_fingerprints_of_the_erosion_stage() {
    // The Earth world of seed 1 at 10,000 cells after the Erosion stage with
    // 10 steps of the default settings (with diffusion): FNV-1a over the
    // bits of each field's values (f64, little-endian). A change means that
    // the same world erodes or drains differently. Recorded after the other
    // tests passed and the maps of out/M3-05 were viewed (ticket M3-05,
    // Result), and again by ticket M3-06 (the land roughened by up to 1 m
    // before the first step, world-breaking) after the tests of
    // `roughness_tests.rs` and of this file passed and the maps of
    // out/M3-06 were viewed; `kuniumi generate` writes the same values (the
    // CLI's golden test, through `.kfield` fingerprints).
    let mut world = with_steps(&worlds()[0], 10);
    run_stage(&mut world);
    let fingerprints: Vec<(&str, u64)> = GOLDEN_FINGERPRINTS
        .iter()
        .map(|&(name, _)| {
            let bytes: Vec<u8> = values(&world, name)
                .iter()
                .flat_map(|v| v.to_bits().to_le_bytes())
                .collect();
            (name, fnv1a64(&bytes))
        })
        .collect();
    for (name, fingerprint) in &fingerprints {
        println!("{name}: {fingerprint:#018X}");
    }
    assert_eq!(fingerprints, GOLDEN_FINGERPRINTS);
}

/// The fingerprints of `golden_fingerprints_of_the_erosion_stage`. Until
/// ticket M3-06 the heights' equalled that of
/// `diffusion_tests::golden_fingerprint_of_an_eroded_and_diffused_earth_world`
/// (ticket M3-04), which runs `erode` on the same world with the same
/// settings outside the pipeline; the stage now erodes the roughened
/// surface instead.
const GOLDEN_FINGERPRINTS: [(&str, u64); 4] = [
    (SURFACE_HEIGHT, 0xFDE1_B5E3_689C_E2AE),
    (DRAINAGE_AREA, 0x381B_B3AE_6288_1477),
    (DISCHARGE, 0x423E_19A2_C079_1331),
    (DEPRESSION_DEPTH, 0xBF73_131B_22D5_2440),
];

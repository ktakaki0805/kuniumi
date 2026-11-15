//! Tests of the hillslope diffusion of the Erosion stage (ticket M3-04):
//!
//! - `diffuse_by_hand`: one land cell on the octahedron relaxes towards
//!   the sea level, not towards the ocean floor;
//! - `land_stays_land_and_the_ocean_stays`: [`diffuse`] on Earth worlds for
//!   steps up to 10⁸ times the explicit limit (where the solver's tolerance
//!   would put some coastal cells at the sea level without the guard of
//!   `diffuse`);
//! - `zero_diffusivity_changes_nothing`;
//! - `erode_diffuses_after_incising`, `erode_steps_compose_with_diffusion`:
//!   a step of [`erode`] is the incision followed by the diffusion;
//! - the panics of [`diffuse`] and the diffusion check of [`erode`];
//! - `golden_fingerprint_of_an_eroded_and_diffused_earth_world`.
#![allow(clippy::disallowed_methods)]

use core::f64::consts::PI;

use kuniumi_math::rng::fnv1a64;
use kuniumi_sphere::SphereGrid;

use super::tests::{dry_world, rain, worlds};
use super::*;
use crate::stages::crust::SURFACE_HEIGHT;
use crate::stages::tectonics::UPLIFT_RATE;

/// The largest `Σ_j (ℓ_ij / d_ij) / a_i` of `grid`, from the grid: the
/// diffusion number of a step (its size in units of the explicit limit)
/// is `N = Δt K_d / R²` times it, and forward
/// Euler is stable for `N ≤ 1`.
fn largest_rate(grid: &SphereGrid) -> f64 {
    (0..grid.len() as u32)
        .map(|i| {
            (0..grid.neighbors(i).len())
                .map(|k| grid.edge_lengths_rad(i)[k] / grid.centre_distances_rad(i)[k])
                .sum::<f64>()
                / grid.area_sr(i)
        })
        .fold(0.0, f64::max)
}

/// The default settings with the diffusivity that gives diffusion number
/// `number` on `grid` at `radius_m` with the default time step.
fn settings_for(grid: &SphereGrid, radius_m: f64, number: f64) -> ErosionSettings {
    let defaults = ErosionSettings::default();
    ErosionSettings {
        diffusivity_m2_s: number / largest_rate(grid) * radius_m * radius_m / defaults.time_step_s,
        ..defaults
    }
}

#[test]
fn diffuse_by_hand() {
    // The octahedron on the unit sphere, sea level 0, K_d = 1 m²/s and
    // Δt = c s: the land cell +Z at 100 m among ocean cells 4 km deep (and
    // −Z exactly at the sea level, so ocean). Its four neighbours are
    // ocean, seen at the sea level: a h' + 4 c w (h' − 0) = a h, so
    // h' = a h / (a + 4 c w), w = acos(1/3) / (π/2), a = 2π/3. With the
    // ocean floor as the boundary it would be (a h − 16,000 c w)/(a + 4cw),
    // below the sea for c ≥ 0.03.
    let grid = SphereGrid::octahedron();
    let laplacian = Laplacian::new(&grid);
    let w = 1.230_959_417_340_774_7 / (PI / 2.0);
    let a = 2.0 * PI / 3.0;
    for c in [1e-2, 1.0, 1e2, 1e6] {
        let settings = ErosionSettings {
            diffusivity_m2_s: 1.0,
            time_step_s: c,
            ..ErosionSettings::default()
        };
        let start = [-4000.0, -4000.0, -4000.0, -4000.0, 100.0, 0.0];
        let mut h = start;
        let solve = diffuse(&laplacian, 1.0, 0.0, &settings, &mut h);
        assert!(solve.converged(), "c = {c}: {solve:?}");
        // Within 1e-12 of the starting height: the step adds the solved
        // change to it, so where the cell falls nearly to the sea level
        // the result keeps the absolute, not the relative, precision.
        let expected = a * 100.0 / (a + 4.0 * c * w);
        assert!(
            (h[4] - expected).abs() <= 1e-12 * 100.0,
            "c = {c}: {} for {expected}",
            h[4]
        );
        assert!(h[4] > 0.0);
        for i in [0, 1, 2, 3, 5] {
            assert_eq!(h[i].to_bits(), start[i].to_bits(), "c = {c}: ocean {i}");
        }
    }
}

#[test]
fn land_stays_land_and_the_ocean_stays() {
    // The Earth worlds of seeds 1–3 and the dry Earth world of seed 1
    // (whose lowest cells lie exactly at the sea level), three steps of
    // diffusion numbers 10⁻² to 10⁸: ocean cells unchanged bit for bit,
    // land cells above the sea level and not above the highest land cell
    // of the step before (the maximum principle with the sea level as the
    // boundary value).
    for world in worlds().iter().chain([dry_world()]) {
        let grid = &world.grid;
        let laplacian = Laplacian::new(grid);
        let radius_m = world.params.radius_m;
        let sea = world.sea_level_m;
        let start = &world.field(SURFACE_HEIGHT).expect("heights").values;
        for number in [1e-2, 1.0, 1e2, 1e4, 1e6, 1e8] {
            let settings = settings_for(grid, radius_m, number);
            let label = format!("seed {}, N = {number:e}", world.seed);
            let mut h = start.clone();
            for step in 0..3 {
                let before = h.clone();
                let top = before.iter().fold(sea, |m, &v| m.max(v));
                let solve = diffuse(&laplacian, radius_m, sea, &settings, &mut h);
                assert!(solve.converged(), "{label}: {solve:?}");
                for i in 0..h.len() {
                    if before[i] <= sea {
                        assert_eq!(h[i].to_bits(), before[i].to_bits(), "{label}: ocean {i}");
                    } else {
                        assert!(
                            sea < h[i] && h[i] <= top + 1e-6,
                            "{label}, step {step}: land cell {i} at {} (sea {sea}, top {top})",
                            h[i]
                        );
                    }
                }
            }
            let lowered = h
                .iter()
                .zip(start)
                .map(|(a, b)| b - a)
                .fold(0.0_f64, f64::max);
            println!("{label}: largest lowering in 3 steps {lowered:.3e} m");
        }
    }
}

#[test]
fn zero_diffusivity_changes_nothing() {
    let world = &worlds()[0];
    let laplacian = Laplacian::new(&world.grid);
    let start = &world.field(SURFACE_HEIGHT).expect("heights").values;
    let settings = ErosionSettings {
        diffusivity_m2_s: 0.0,
        ..ErosionSettings::default()
    };
    let mut h = start.clone();
    let solve = diffuse(
        &laplacian,
        world.params.radius_m,
        world.sea_level_m,
        &settings,
        &mut h,
    );
    assert_eq!(solve.iterations, 0);
    assert!(h.iter().zip(start).all(|(a, b)| a.to_bits() == b.to_bits()));
}

/// One step of [`erode`] written out: the routing, discharge and incision
/// of M3-03, then [`diffuse`].
fn step_by_hand(inputs: &ErosionInputs<'_>, settings: &ErosionSettings, height_m: &mut [f64]) {
    let grid = inputs.grid;
    let radius_m = inputs.radius_m;
    let inflow: Vec<f64> = (0..grid.len() as u32)
        .map(|i| inputs.precipitation_m_s[i as usize] * grid.area_sr(i) * radius_m * radius_m)
        .collect();
    let filled = fill_depressions(grid, height_m, inputs.sea_level_m);
    let drainage = Drainage::build(grid, &filled, inputs.sea_level_m, radius_m);
    let discharge = drainage.accumulate(&inflow);
    let distance = receiver_distances_m(grid, &drainage, radius_m);
    let network = Network {
        receivers: drainage.receivers(),
        stack: drainage.stack(),
        receiver_distance_m: &distance,
        discharge_m3_s: &discharge,
    };
    incise(
        &network,
        inputs.uplift_rate_m_s,
        inputs.sea_level_m,
        settings,
        height_m,
    );
    diffuse(
        &Laplacian::new(grid),
        radius_m,
        inputs.sea_level_m,
        settings,
        height_m,
    );
}

#[test]
fn erode_diffuses_after_incising() {
    // The default settings, and a diffusivity 10⁶ times larger (so that
    // the diffusion is felt), and none: one step of `erode` equals the
    // step written out, bit for bit.
    let world = &worlds()[0];
    let uplift = &world.field(UPLIFT_RATE).expect("uplift").values;
    let precipitation = rain(&world.grid);
    let inputs = ErosionInputs {
        grid: &world.grid,
        radius_m: world.params.radius_m,
        sea_level_m: world.sea_level_m,
        uplift_rate_m_s: uplift,
        precipitation_m_s: &precipitation,
    };
    let start = &world.field(SURFACE_HEIGHT).expect("heights").values;
    let defaults = ErosionSettings {
        steps: 1,
        ..ErosionSettings::default()
    };
    let mut changes = Vec::new();
    for diffusivity_m2_s in [
        defaults.diffusivity_m2_s,
        1e6 * defaults.diffusivity_m2_s,
        0.0,
    ] {
        let settings = ErosionSettings {
            diffusivity_m2_s,
            ..defaults
        };
        let mut by_erode = start.clone();
        erode(&inputs, &settings, &mut by_erode, &mut |_| {});
        let mut by_hand = start.clone();
        step_by_hand(&inputs, &settings, &mut by_hand);
        assert!(
            by_erode
                .iter()
                .zip(&by_hand)
                .all(|(a, b)| a.to_bits() == b.to_bits()),
            "K_d = {diffusivity_m2_s}"
        );
        changes.push(by_erode);
    }
    // What the diffusion changed in this step, against the incision alone.
    for (k, name) in [(0, "default"), (1, "10⁶ × default")] {
        let largest = changes[k]
            .iter()
            .zip(&changes[2])
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max);
        println!("diffusivity {name}: largest change by the diffusion in one step {largest:.3e} m");
        assert!(largest > 0.0, "{name}");
    }
}

#[test]
fn erode_steps_compose_with_diffusion() {
    // Five steps at once equal five single steps, bit for bit, with a
    // diffusivity large enough to be felt.
    let world = &worlds()[1];
    let uplift = &world.field(UPLIFT_RATE).expect("uplift").values;
    let precipitation = rain(&world.grid);
    let inputs = ErosionInputs {
        grid: &world.grid,
        radius_m: world.params.radius_m,
        sea_level_m: world.sea_level_m,
        uplift_rate_m_s: uplift,
        precipitation_m_s: &precipitation,
    };
    let settings = settings_for(&world.grid, world.params.radius_m, 1.0);
    let start = &world.field(SURFACE_HEIGHT).expect("heights").values;
    let mut at_once = start.clone();
    erode(
        &inputs,
        &ErosionSettings {
            steps: 5,
            ..settings
        },
        &mut at_once,
        &mut |_| {},
    );
    let mut one_by_one = start.clone();
    for _ in 0..5 {
        erode(
            &inputs,
            &ErosionSettings {
                steps: 1,
                ..settings
            },
            &mut one_by_one,
            &mut |_| {},
        );
    }
    assert!(
        at_once
            .iter()
            .zip(&one_by_one)
            .all(|(a, b)| a.to_bits() == b.to_bits())
    );
}

// ---------------------------------------------------------------------------
// Panics.

/// Runs `diffuse` on the octahedron after `change` has altered its valid
/// inputs (radius, sea level, settings, heights).
fn diffuse_octahedron(
    change: impl FnOnce(&mut f64, &mut f64, &mut ErosionSettings, &mut Vec<f64>),
) {
    let laplacian = Laplacian::new(&SphereGrid::octahedron());
    let (mut radius_m, mut sea_level_m) = (1.0, 0.0);
    let mut settings = ErosionSettings {
        diffusivity_m2_s: 1.0,
        time_step_s: 1.0,
        ..ErosionSettings::default()
    };
    let mut h = vec![-10.0, 10.0, 20.0, 30.0, 40.0, 50.0];
    change(&mut radius_m, &mut sea_level_m, &mut settings, &mut h);
    diffuse(&laplacian, radius_m, sea_level_m, &settings, &mut h);
}

#[test]
fn the_unbroken_diffusion_runs() {
    diffuse_octahedron(|_, _, _, _| {});
}

#[test]
#[should_panic(expected = "diffuse needs one height per cell")]
fn diffuse_needs_one_height_per_cell() {
    diffuse_octahedron(|_, _, _, h| h.push(0.0));
}

#[test]
#[should_panic(expected = "diffuse needs a finite positive radius")]
fn diffuse_needs_a_finite_positive_radius() {
    diffuse_octahedron(|r, _, _, _| *r = f64::INFINITY);
}

#[test]
#[should_panic(expected = "diffuse needs a finite sea level")]
fn diffuse_needs_a_finite_sea_level() {
    diffuse_octahedron(|_, s, _, _| *s = f64::NAN);
}

#[test]
#[should_panic(expected = "diffuse needs finite heights")]
fn diffuse_needs_finite_heights() {
    diffuse_octahedron(|_, _, _, h| h[0] = f64::NEG_INFINITY);
}

#[test]
#[should_panic(expected = "diffuse needs valid settings")]
fn diffuse_needs_valid_settings() {
    diffuse_octahedron(|_, _, s, _| s.diffusivity_m2_s = -1.0);
}

#[test]
#[should_panic(expected = "diffuse needs a finite step Δt K_d / R²")]
fn diffuse_needs_a_finite_step() {
    // Valid settings (`Δt K_d` = 1 m²) on a radius whose square underflows
    // to 0. A product `Δt K_d` that overflows is not valid (ticket M3-07),
    // so it is rejected as such.
    diffuse_octahedron(|r, _, _, _| *r = 1e-200);
}

#[test]
#[should_panic(
    expected = "diffuse needs valid settings: diffusion_step_m2 from diffusivity_m2_s \
                           and time_step_s is out of range: inf"
)]
fn diffuse_needs_a_diffusion_step_that_does_not_overflow() {
    diffuse_octahedron(|_, _, s, _| {
        s.diffusivity_m2_s = 1e200;
        s.time_step_s = 1e200;
    });
}

#[test]
#[should_panic(expected = "erode needs a finite step Δt K_d / R²")]
fn erode_needs_a_finite_step() {
    // Checked before the first step, so with 0 steps as well.
    let world = &worlds()[0];
    let n = world.grid.len();
    let (no_uplift, precipitation) = (vec![0.0; n], rain(&world.grid));
    // The default settings, valid (`Δt K_d` = 50 m²), on a radius whose
    // square underflows to 0; a product `Δt K_d` that overflows is not
    // valid (ticket M3-07).
    let inputs = ErosionInputs {
        grid: &world.grid,
        radius_m: 1e-200,
        sea_level_m: world.sea_level_m,
        uplift_rate_m_s: &no_uplift,
        precipitation_m_s: &precipitation,
    };
    let settings = ErosionSettings {
        steps: 0,
        ..ErosionSettings::default()
    };
    let mut h = world.field(SURFACE_HEIGHT).expect("heights").values.clone();
    erode(&inputs, &settings, &mut h, &mut |_| {});
}

// ---------------------------------------------------------------------------
// Golden value.

#[test]
fn golden_fingerprint_of_an_eroded_and_diffused_earth_world() {
    // The Earth world of seed 1 at 10,000 cells after 10 steps of the
    // default settings, hillslope diffusion included, with its own uplift
    // and the placeholder precipitation: FNV-1a over the bits of the
    // heights (f64, little-endian). A change means that the same world
    // erodes differently. Recorded after the other tests passed and the
    // images of out/M3-04 were viewed (ticket M3-04, Result).
    let world = &worlds()[0];
    let uplift = &world.field(UPLIFT_RATE).expect("uplift").values;
    let precipitation = rain(&world.grid);
    let inputs = ErosionInputs {
        grid: &world.grid,
        radius_m: world.params.radius_m,
        sea_level_m: world.sea_level_m,
        uplift_rate_m_s: uplift,
        precipitation_m_s: &precipitation,
    };
    let settings = ErosionSettings {
        steps: 10,
        ..ErosionSettings::default()
    };
    let mut h = world.field(SURFACE_HEIGHT).expect("heights").values.clone();
    erode(&inputs, &settings, &mut h, &mut |_| {});
    let bytes: Vec<u8> = h.iter().flat_map(|v| v.to_bits().to_le_bytes()).collect();
    let fingerprint = fnv1a64(&bytes);
    println!("fingerprint {fingerprint:#018X}");
    assert_eq!(fingerprint, GOLDEN_FINGERPRINT);
}

/// The fingerprint of `golden_fingerprint_of_an_eroded_and_diffused_earth_world`.
const GOLDEN_FINGERPRINT: u64 = 0xB715_B888_F60A_A864;

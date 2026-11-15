//! Tests of the river incision and its driver (ticket M3-03):
//!
//! - `one_step_by_hand` and `extreme_settings_stay_finite`: [`incise`] on
//!   small networks given as slices, values by hand;
//! - `steady_state_of_a_1d_channel`: the slope of a chain of cells converges
//!   to `U / (K Qᵐ)` for time steps of 10⁴, 10⁵ and 10⁶ years;
//! - `large_time_steps_are_stable`: [`erode`] on the Earth worlds of seeds
//!   1–3 and on random surfaces, without uplift, with uniform uplift and
//!   with the Earth worlds' own: no NaN, no cell below its receiver, no new
//!   pit and no new closed depression;
//! - `land_volume_decreases_without_uplift`;
//! - `steady_state_of_an_island`: the driver end to end on the sphere;
//! - `erode_reports_progress_per_step`, `erode_steps_compose`,
//!   `receiver_distances_by_hand`;
//! - the panics of [`incise`] and [`erode`];
//! - `golden_fingerprint_of_an_eroded_earth_world`.
//!
//! Since ticket M3-04 a step of [`erode`] also diffuses the land; these
//! tests run it without diffusion (`settings`), and `diffusion_tests.rs`
//! tests the whole step.
#![allow(clippy::disallowed_methods)]

use core::cmp::Ordering;
use core::f64::consts::PI;
use std::collections::BinaryHeap;
use std::sync::OnceLock;

use kuniumi_math::rng::fnv1a64;
use kuniumi_math::units::SECONDS_PER_YEAR;
use kuniumi_math::{Rng, Vec3};
use kuniumi_sphere::SphereGrid;
use kuniumi_world::{PlanetParams, Settings, World};

use super::*;
use crate::stages::crust::SURFACE_HEIGHT;
use crate::stages::tectonics::UPLIFT_RATE;

/// Earth's radius (m), `PlanetParams::earth().radius_m`.
const EARTH_RADIUS_M: f64 = 6_371_000.0;

/// The cells of the grids of the tests on the sphere.
const CELLS: u32 = 10_000;

/// An uplift rate of 1 mm/yr, in m/s.
const ONE_MM_PER_YEAR_M_S: f64 = 1.0e-3 / SECONDS_PER_YEAR;

/// The time steps of the stability criterion, in years.
const STEP_YEARS: [f64; 3] = [1.0e4, 1.0e5, 1.0e6];

/// The 10,000-cell grid of seed 1, built once.
pub(super) fn grid() -> &'static SphereGrid {
    static GRID: OnceLock<SphereGrid> = OnceLock::new();
    GRID.get_or_init(|| SphereGrid::fibonacci(CELLS, 1).expect("grid"))
}

/// Earth's worlds at 10,000 cells for seeds 1–3, through the Crust stage,
/// generated once.
pub(super) fn worlds() -> &'static [World] {
    static WORLDS: OnceLock<Vec<World>> = OnceLock::new();
    WORLDS.get_or_init(|| {
        let settings = Settings {
            cells: CELLS,
            allow_partial: true,
            ..Settings::default()
        };
        (1..=3)
            .map(|seed| {
                let world = crate::generate_through_crust(&PlanetParams::earth(), &settings, seed)
                    .expect("a partial world");
                assert_eq!(world.stages, ["grid", "tectonics", "crust"], "seed {seed}");
                world
            })
            .collect()
    })
}

/// Earth's world of seed 1 at 10,000 cells without water
/// (`ocean_volume_m3 = 0`), through the Crust stage, generated once: the
/// sea level is the lowest height, which the trench cells share, so
/// several cells lie exactly at the sea level.
pub(super) fn dry_world() -> &'static World {
    static WORLD: OnceLock<World> = OnceLock::new();
    WORLD.get_or_init(|| {
        let settings = Settings {
            cells: CELLS,
            allow_partial: true,
            ..Settings::default()
        };
        let params = PlanetParams {
            ocean_volume_m3: 0.0,
            ..PlanetParams::earth()
        };
        let world = crate::generate_through_crust(&params, &settings, 1).expect("a partial world");
        assert_eq!(world.stages, ["grid", "tectonics", "crust"]);
        world
    })
}

/// A smooth random surface: a sum of 24 plane waves of random direction,
/// wavenumber 2–12 per radian and amplitude 1,000 m over the wavenumber
/// (as in the tests of ticket M3-02).
fn waves_m(grid: &SphereGrid) -> Vec<f64> {
    let mut rng = Rng::for_stream(5, "test/erosion/waves");
    let waves: Vec<(Vec3, f64, f64, f64)> = (0..24)
        .map(|_| {
            let z = rng.range_f64(-1.0, 1.0);
            let lon = rng.range_f64(-PI, PI);
            let r = (1.0 - z * z).sqrt();
            let dir = Vec3::new(r * lon.cos(), r * lon.sin(), z);
            let k = rng.range_f64(2.0, 12.0);
            (dir, k, rng.range_f64(0.0, 2.0 * PI), 1000.0 / k)
        })
        .collect();
    grid.positions()
        .iter()
        .map(|&p| {
            waves
                .iter()
                .map(|&(w, k, phase, a)| a * (k * p.dot(w) + phase).sin())
                .sum()
        })
        .collect()
}

/// A surface to erode.
struct Surface<'a> {
    name: String,
    grid: &'a SphereGrid,
    radius_m: f64,
    sea_level_m: f64,
    heights_m: Vec<f64>,
    /// The `uplift_rate` of an Earth world.
    own_uplift_m_s: Option<Vec<f64>>,
}

/// The surfaces of the stability and volume tests: the Earth worlds of
/// seeds 1–3 and the dry Earth world of seed 1 (with their uplift); a
/// smooth random surface with sea level −150 m; heights uniform in
/// [−1,000, 1,000) m with sea level −900 m, full of closed depressions.
fn surfaces() -> Vec<Surface<'static>> {
    let mut out: Vec<Surface<'static>> = worlds()
        .iter()
        .chain([dry_world()])
        .map(|w| Surface {
            name: if w.params.ocean_volume_m3 > 0.0 {
                format!("Earth, seed {}", w.seed)
            } else {
                format!("dry Earth, seed {}", w.seed)
            },
            grid: &w.grid,
            radius_m: w.params.radius_m,
            sea_level_m: w.sea_level_m,
            heights_m: w.field(SURFACE_HEIGHT).expect("heights").values.clone(),
            own_uplift_m_s: Some(w.field(UPLIFT_RATE).expect("uplift").values.clone()),
        })
        .collect();
    let grid = grid();
    let mut rng = Rng::for_stream(1, "test/erosion/uniform");
    let uniform = (0..grid.len())
        .map(|_| rng.range_f64(-1000.0, 1000.0))
        .collect();
    out.push(Surface {
        name: "smooth random".into(),
        grid,
        radius_m: EARTH_RADIUS_M,
        sea_level_m: -150.0,
        heights_m: waves_m(grid),
        own_uplift_m_s: None,
    });
    out.push(Surface {
        name: "uniform random".into(),
        grid,
        radius_m: EARTH_RADIUS_M,
        sea_level_m: -900.0,
        heights_m: uniform,
        own_uplift_m_s: None,
    });
    out
}

/// The default settings with a time step of `years` and `steps` steps,
/// without hillslope diffusion (`diffusivity_m2_s` 0): these tests check
/// the incision's own guarantees, such as "no height rises", which the
/// diffusion of each step (ticket M3-04) does not keep; the tests of the
/// whole step are in `diffusion_tests.rs`.
fn settings(years: f64, steps: u32) -> ErosionSettings {
    ErosionSettings {
        diffusivity_m2_s: 0.0,
        time_step_s: years * SECONDS_PER_YEAR,
        steps,
        ..ErosionSettings::default()
    }
}

/// The placeholder precipitation on every cell of `grid`.
pub(super) fn rain(grid: &SphereGrid) -> Vec<f64> {
    vec![PLACEHOLDER_PRECIPITATION_M_S; grid.len()]
}

/// The drainage graph that a step of [`erode`] follows from the surface
/// `height_m`: built on the filled surface.
fn routing(grid: &SphereGrid, height_m: &[f64], sea_level_m: f64, radius_m: f64) -> Drainage {
    let filled = fill_depressions(grid, height_m, sea_level_m);
    Drainage::build(grid, &filled, sea_level_m, radius_m)
}

/// The water surface after a step along the graph `d` the step followed:
/// the sea level on the cells that were ocean, the height itself on a
/// base level, `max(h_i, W_r)` elsewhere, in stack order.
fn water_surface_m(d: &Drainage, height_m: &[f64], sea_level_m: f64, ocean: &[bool]) -> Vec<f64> {
    let mut w = vec![0.0; height_m.len()];
    for &i in d.stack() {
        let i = i as usize;
        let r = d.receiver(i as u32) as usize;
        w[i] = if ocean[i] {
            sea_level_m
        } else if r == i || height_m[i] > w[r] {
            height_m[i]
        } else {
            w[r]
        };
    }
    w
}

/// The spill height of every cell (m): the lowest, over the paths to the
/// ocean, of the highest height on the path. Priority-Flood without ε
/// (Barnes, Lehman & Mulla 2014, Algorithm 1), written here independently
/// of the code under test; every surface of these tests has ocean.
fn spill_heights_m(grid: &SphereGrid, height_m: &[f64], sea_level_m: f64) -> Vec<f64> {
    /// A min-heap entry: the spill height and the cell.
    struct Open(f64, u32);
    impl PartialEq for Open {
        fn eq(&self, other: &Self) -> bool {
            self.cmp(other) == Ordering::Equal
        }
    }
    impl Eq for Open {}
    impl PartialOrd for Open {
        fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
            Some(self.cmp(other))
        }
    }
    impl Ord for Open {
        fn cmp(&self, other: &Self) -> Ordering {
            other.0.total_cmp(&self.0).then(other.1.cmp(&self.1))
        }
    }
    let mut spill = height_m.to_vec();
    let mut closed: Vec<bool> = height_m.iter().map(|&h| h <= sea_level_m).collect();
    assert!(closed.contains(&true), "the surface has ocean");
    let mut open: BinaryHeap<Open> = (0..grid.len() as u32)
        .filter(|&i| closed[i as usize])
        .map(|i| Open(height_m[i as usize], i))
        .collect();
    while let Some(Open(z, c)) = open.pop() {
        for &j in grid.neighbors(c) {
            let ju = j as usize;
            if !closed[ju] {
                closed[ju] = true;
                spill[ju] = if height_m[ju] < z { z } else { height_m[ju] };
                open.push(Open(spill[ju], j));
            }
        }
    }
    spill
}

/// Whether land cell `i` has a neighbour strictly lower than itself.
fn has_lower_neighbour(grid: &SphereGrid, height_m: &[f64], i: u32) -> bool {
    grid.neighbors(i)
        .iter()
        .any(|&j| height_m[j as usize] < height_m[i as usize])
}

/// The volume of land above the sea level (m³), `Σ (h − sea level) a R²`
/// over the land cells in index order.
fn land_volume_m3(grid: &SphereGrid, height_m: &[f64], sea_level_m: f64, radius_m: f64) -> f64 {
    (0..grid.len() as u32)
        .filter(|&i| height_m[i as usize] > sea_level_m)
        .map(|i| (height_m[i as usize] - sea_level_m) * grid.area_sr(i) * radius_m * radius_m)
        .sum()
}

/// `|a − b| / |b|`.
fn relative(a: f64, b: f64) -> f64 {
    (a - b).abs() / b.abs()
}

/// The least double above `x > 0`.
fn up(x: f64) -> f64 {
    assert!(x > 0.0);
    f64::from_bits(x.to_bits() + 1)
}

// ---------------------------------------------------------------------------
// By hand.

#[test]
fn one_step_by_hand() {
    // K = 0.25, m = 0.5 and Δt = 1 s, so that F = 0.25 · √Q / L, with
    // √Q exact for Q = 4, 16, 64 and every division below exact. Sea
    // level 0. The update is the ticket's (a + F b) / (1 + F), with
    // a = h + Δt U and b the water surface of the receiver.
    let settings = ErosionSettings {
        erodibility: 0.25,
        discharge_exponent: 0.5,
        diffusivity_m2_s: 0.0,
        time_step_s: 1.0,
        steps: 1,
    };
    // 0: ocean, 4 km deep, uplift ignored.
    // 1: coast, h 100, U 2, Q 4, L 0.5: F = 1, b = sea level 0 (not the
    //    sea floor): (102 + 0) / 2 = 51.
    // 2: h 81, Q 16, L 0.25 into 1: F = 4, b = 51: (81 + 204) / 5 = 57.
    // 3: a lake floor at 40 into 2, U 1: a = 41 ≤ 57, only uplifted.
    // 4: a river entering the lake, h 66, Q 64, L 0.25 into 3: F = 8,
    //    b = the lake's surface 57, not its floor 41: (66 + 456) / 9 = 58.
    // 5: h 200 into 2 with L = 1e-300: F = 5e299, the update rounds to
    //    the base, 57, and the cell is kept one ULP above it.
    // 6: h 300 into 2 with Q = 0, U 1: no erosion, 301.
    // 7: a pit (land base level), h 500, U 3: 503.
    // 8: a cell exactly at the sea level, U 5: ocean, so fixed at 0.
    // 9: h 600 into the pit 7, Q 4, L 0.5: F = 1, b = 503:
    //    (600 + 503) / 2 = 551.5.
    // 10: h 20 into 8, Q 4, L 0.5: F = 1, b = the sea level 0: 10.
    // 11: a pit at 2 m sinking at 7 m/s: a = −5, below the sea level, so
    //     it becomes ocean and its water surface is the sea level 0.
    // 12: h 20 into 11, Q 4, L 0.5: F = 1, b = 0 (not −5): 10.
    let receivers = [0, 0, 1, 2, 3, 2, 2, 7, 8, 7, 8, 11, 11];
    let stack = [0, 7, 8, 11, 1, 2, 3, 5, 6, 9, 10, 12, 4];
    let distance_m = [
        0.0, 0.5, 0.25, 1.0, 0.25, 1e-300, 1.0, 0.0, 0.0, 0.5, 0.5, 0.0, 0.5,
    ];
    let discharge_m3_s = [
        100.0, 4.0, 16.0, 64.0, 64.0, 4.0, 0.0, 4.0, 8.0, 4.0, 4.0, 8.0, 4.0,
    ];
    let uplift_m_s = [
        1e-3, 2.0, 0.0, 1.0, 0.0, 0.0, 1.0, 3.0, 5.0, 0.0, 0.0, -7.0, 0.0,
    ];
    let mut h = [
        -4000.0, 100.0, 81.0, 40.0, 66.0, 200.0, 300.0, 500.0, 0.0, 600.0, 20.0, 2.0, 20.0,
    ];
    let network = Network {
        receivers: &receivers,
        stack: &stack,
        receiver_distance_m: &distance_m,
        discharge_m3_s: &discharge_m3_s,
    };
    incise(&network, &uplift_m_s, 0.0, &settings, &mut h);
    let expected = [
        -4000.0,
        51.0,
        57.0,
        41.0,
        58.0,
        up(57.0),
        301.0,
        503.0,
        0.0,
        551.5,
        10.0,
        -5.0,
        10.0,
    ];
    for (i, (&got, &want)) in h.iter().zip(&expected).enumerate() {
        assert_eq!(got.to_bits(), want.to_bits(), "cell {i}: {got} for {want}");
    }
}

#[test]
fn extreme_settings_stay_finite() {
    // Δt K underflows to 0 while Qᵐ overflows to ∞, and Δt K overflows to
    // ∞ while Qᵐ is 0: neither may give 0 · ∞ = NaN. In both cases the
    // cell is not eroded.
    for (erodibility, time_step_s, discharge_exponent, q) in
        [(1e-200, 1e-200, 2.0, 1e200), (1e200, 1e200, 1.0, 0.0)]
    {
        let settings = ErosionSettings {
            erodibility,
            discharge_exponent,
            diffusivity_m2_s: 0.0,
            time_step_s,
            steps: 1,
        };
        let mut h = [-10.0, 10.0];
        let network = Network {
            receivers: &[0, 0],
            stack: &[0, 1],
            receiver_distance_m: &[0.0, 1.0],
            discharge_m3_s: &[q, q],
        };
        incise(&network, &[0.0, 0.0], 0.0, &settings, &mut h);
        assert_eq!(
            h.map(f64::to_bits),
            [-10.0, 10.0].map(f64::to_bits),
            "{settings:?}"
        );
    }
    // Δt K Qᵐ / L overflows to ∞: the cell ends one ULP above its base.
    let settings = ErosionSettings {
        erodibility: 1e200,
        discharge_exponent: 1.0,
        diffusivity_m2_s: 0.0,
        time_step_s: 1e200,
        steps: 1,
    };
    let mut h = [-10.0, 10.0, 20.0];
    let network = Network {
        receivers: &[0, 0, 1],
        stack: &[0, 1, 2],
        receiver_distance_m: &[0.0, 1.0, 1.0],
        discharge_m3_s: &[1.0, 1.0, 1.0],
    };
    incise(&network, &[0.0; 3], 0.0, &settings, &mut h);
    assert_eq!(h[0].to_bits(), (-10.0_f64).to_bits());
    assert_eq!(
        h[1].to_bits(),
        5e-324_f64.to_bits(),
        "one ULP above the sea level 0"
    );
    assert_eq!(h[2].to_bits(), 1e-323_f64.to_bits(), "one ULP above cell 1");
    // a − b overflows to ∞ as well as F: still no ∞ / ∞; the cell ends one
    // ULP above its base, the sea level −1.7e308 m.
    let sea_level_m = -1.7e308;
    let mut h = [sea_level_m, 1.7e308];
    let network = Network {
        receivers: &[0, 0],
        stack: &[0, 1],
        receiver_distance_m: &[0.0, 1.0],
        discharge_m3_s: &[1.0, 1.0],
    };
    incise(&network, &[0.0; 2], sea_level_m, &settings, &mut h);
    assert_eq!(h[0].to_bits(), sea_level_m.to_bits());
    assert_eq!(
        h[1].to_bits(),
        sea_level_m.to_bits() - 1,
        "the next double above"
    );
}

#[test]
fn rounding_never_raises_an_eroded_cell() {
    // a = 1 + 2⁻⁵², b = −2⁻⁵³ (the sea level of the ocean receiver) and
    // F = 1e-30, so that 1 + F = 1: a − b = 1 + 1.5 ULP rounds (ties to
    // even) up to 1 + 2⁻⁵¹, and b plus that, 1 + 1.5 ULP again, rounds up
    // to 1 + 2⁻⁵¹ > a. Erosion must not raise the cell: it stays at a.
    let settings = ErosionSettings {
        erodibility: 1e-30,
        discharge_exponent: 1.0,
        diffusivity_m2_s: 0.0,
        time_step_s: 1.0,
        steps: 1,
    };
    let a = 1.0 + f64::EPSILON;
    let b = -f64::EPSILON / 2.0;
    let mut h = [-1.0, a];
    let network = Network {
        receivers: &[0, 0],
        stack: &[0, 1],
        receiver_distance_m: &[0.0, 1.0],
        discharge_m3_s: &[1.0, 1.0],
    };
    // The rounding the clamp guards against, step by step.
    let above = 1.0 + 2.0 * f64::EPSILON;
    assert_eq!((a - b).to_bits(), above.to_bits());
    assert_eq!((b + (a - b) / (1.0 + 1e-30)).to_bits(), above.to_bits());
    incise(&network, &[0.0, 0.0], b, &settings, &mut h);
    assert_eq!(h[1].to_bits(), a.to_bits());
}

// ---------------------------------------------------------------------------
// Criterion 1: steady state of a 1D channel.

#[test]
fn steady_state_of_a_1d_channel() {
    // 50 cells: cell 0 an ocean outlet 100 m deep (sea level 0), cell i
    // draining to i − 1 over 1 km, with the discharge of the cells
    // upstream of it, 1 m³/s each (Q_i = 50 − i), uniform uplift 1 mm/yr
    // and the default K and m. At steady state the slope of every land
    // cell to its base (the sea level for cell 1) is U / (K Qᵐ), computed
    // here with std `powf`.
    const N: usize = 50;
    let receivers: Vec<u32> = (0..N as u32).map(|i| i.saturating_sub(1)).collect();
    let stack: Vec<u32> = (0..N as u32).collect();
    let mut distance_m = vec![1000.0; N];
    distance_m[0] = 0.0;
    let discharge_m3_s: Vec<f64> = (0..N).map(|i| (N - i) as f64).collect();
    let uplift_m_s = vec![ONE_MM_PER_YEAR_M_S; N];
    let network = Network {
        receivers: &receivers,
        stack: &stack,
        receiver_distance_m: &distance_m,
        discharge_m3_s: &discharge_m3_s,
    };
    for years in STEP_YEARS {
        let settings = settings(years, 1);
        let mut h = vec![1.0; N];
        h[0] = -100.0;
        let mut steps = 0;
        let mut change_m = f64::INFINITY;
        while change_m > 1e-12 && steps < 1_000_000 {
            let before = h.clone();
            incise(&network, &uplift_m_s, 0.0, &settings, &mut h);
            change_m = h
                .iter()
                .zip(&before)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0, f64::max);
            steps += 1;
        }
        assert!(
            change_m <= 1e-12,
            "{years} yr: not converged after {steps} steps"
        );
        assert_eq!(h[0].to_bits(), (-100.0_f64).to_bits(), "the ocean is fixed");
        let mut worst = 0.0_f64;
        for i in 1..N {
            let base_m = if i == 1 { 0.0 } else { h[i - 1] };
            let slope = (h[i] - base_m) / 1000.0;
            let expected = ONE_MM_PER_YEAR_M_S
                / (ErosionSettings::DEFAULT_ERODIBILITY
                    * discharge_m3_s[i].powf(ErosionSettings::DEFAULT_DISCHARGE_EXPONENT));
            let error = relative(slope, expected);
            assert!(
                error < 1e-6,
                "{years} yr, cell {i}: slope {slope} for {expected}"
            );
            worst = worst.max(error);
        }
        println!(
            "Δt = {years} yr: steady after {steps} steps, divide at {:.3} m, worst slope error \
             {worst:.2e}",
            h[N - 1]
        );
    }
}

// ---------------------------------------------------------------------------
// Criterion 2: stability for large time steps.

#[test]
fn large_time_steps_are_stable() {
    for surface in surfaces() {
        let grid = surface.grid;
        let n = grid.len();
        let sea = surface.sea_level_m;
        let mut uplifts: Vec<(&str, Vec<f64>)> = vec![
            ("no uplift", vec![0.0; n]),
            ("uniform uplift", vec![ONE_MM_PER_YEAR_M_S; n]),
        ];
        if let Some(own) = &surface.own_uplift_m_s {
            uplifts.push(("own uplift", own.clone()));
        }
        let precipitation = rain(grid);
        for (uplift_name, uplift) in &uplifts {
            let uniform = *uplift_name != "own uplift";
            let inputs = ErosionInputs {
                grid,
                radius_m: surface.radius_m,
                sea_level_m: sea,
                uplift_rate_m_s: uplift,
                precipitation_m_s: &precipitation,
            };
            for years in STEP_YEARS {
                let settings = settings(years, 1);
                let label = format!("{}, {uplift_name}, {years} yr", surface.name);
                let mut h = surface.heights_m.clone();
                let mut spill_before = spill_heights_m(grid, &h, sea);
                let (mut eroded, mut submerged, mut dammed_pits) = (0, 0, 0);
                for step in 0..10 {
                    let before = h.clone();
                    let d = routing(grid, &before, sea, surface.radius_m);
                    erode(&inputs, &settings, &mut h, &mut |_| {});
                    let ocean: Vec<bool> = before.iter().map(|&v| v <= sea).collect();
                    assert!(h.iter().all(|v| v.is_finite()), "{label}, step {step}");
                    let w = water_surface_m(&d, &h, sea, &ocean);
                    for i in 0..n {
                        if ocean[i] {
                            assert_eq!(h[i].to_bits(), before[i].to_bits(), "{label}: ocean {i}");
                            continue;
                        }
                        // Uplift is never negative here, so land stays land.
                        assert!(h[i] > sea, "{label}, step {step}: cell {i} drowned");
                        let a = before[i] + settings.time_step_s * uplift[i];
                        let r = d.receiver(i as u32) as usize;
                        if r != i && a > w[r] {
                            eroded += 1;
                            assert!(
                                w[r] < h[i] && h[i] <= a,
                                "{label}, step {step}: cell {i} at {} for base {} and {a}",
                                h[i],
                                w[r]
                            );
                        } else {
                            submerged += usize::from(r != i);
                            assert_eq!(h[i].to_bits(), a.to_bits(), "{label}, step {step}: {i}");
                        }
                    }
                    let spill_after = spill_heights_m(grid, &h, sea);
                    for i in 0..n as u32 {
                        let iu = i as usize;
                        if ocean[iu] {
                            continue;
                        }
                        let new_pit = has_lower_neighbour(grid, &before, i)
                            && !has_lower_neighbour(grid, &h, i);
                        let new_depression =
                            spill_before[iu] <= before[iu] && spill_after[iu] > h[iu];
                        if uniform {
                            assert!(!new_pit, "{label}, step {step}: new pit {i}");
                            assert!(!new_depression, "{label}, step {step}: new depression {i}");
                        } else {
                            dammed_pits += usize::from(new_pit);
                        }
                    }
                    spill_before = spill_after;
                }
                println!(
                    "{label}: {eroded} cells eroded, {submerged} under water, {dammed_pits} \
                     pits dammed by uplift"
                );
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Criterion 3: the land volume without uplift.

#[test]
fn land_volume_decreases_without_uplift() {
    for surface in surfaces() {
        let grid = surface.grid;
        let sea = surface.sea_level_m;
        let no_uplift = vec![0.0; grid.len()];
        let precipitation = rain(grid);
        let inputs = ErosionInputs {
            grid,
            radius_m: surface.radius_m,
            sea_level_m: sea,
            uplift_rate_m_s: &no_uplift,
            precipitation_m_s: &precipitation,
        };
        for years in [1.0e4, 1.0e6] {
            let settings = settings(years, 1);
            let mut h = surface.heights_m.clone();
            let start_m3 = land_volume_m3(grid, &h, sea, surface.radius_m);
            let mut volume_m3 = start_m3;
            for step in 0..20 {
                let before = h.clone();
                erode(&inputs, &settings, &mut h, &mut |_| {});
                for (i, (&now, &was)) in h.iter().zip(&before).enumerate() {
                    assert!(
                        now <= was,
                        "{}, {years} yr, step {step}: {i} rose",
                        surface.name
                    );
                    assert_eq!(now > sea, was > sea, "{}: cell {i}", surface.name);
                }
                let next_m3 = land_volume_m3(grid, &h, sea, surface.radius_m);
                assert!(
                    next_m3 < volume_m3,
                    "{}, {years} yr, step {step}: {next_m3} m³ after {volume_m3} m³",
                    surface.name
                );
                volume_m3 = next_m3;
            }
            println!(
                "{}, {years} yr: land volume {:.4e} m³ → {:.4e} m³ in 20 steps",
                surface.name, start_m3, volume_m3
            );
        }
    }
}

// ---------------------------------------------------------------------------
// The driver end to end.

#[test]
fn steady_state_of_an_island() {
    // A cap of land 0.6 rad in radius at 10 m, ocean 1 km deep around it,
    // uniform uplift of 1 mm/yr and steps of 10⁶ yr, on planets of radius
    // 1,500, 6,371 and 13,000 km (the bounds of `PlanetParams` and Earth):
    // once no height changes any more, every land cell's slope to the
    // water surface of its receiver is U / (K Qᵐ). The test takes the
    // distance to the receiver from the great-circle distance of the two
    // centres (std `acos`) and the discharge from a walk down the
    // receivers of every cell, not from the code under test; it also
    // checks `receiver_distances_m` against that distance on this grid,
    // where the neighbours of a cell are at different distances. The rain
    // varies with z, P_i = P (1 + z_i / 2), so that each cell's own rate
    // must reach the discharge.
    let grid = grid();
    let n = grid.len();
    let centre = Vec3::new(1.0, 2.0, 3.0).normalized().expect("a direction");
    let uplift = vec![ONE_MM_PER_YEAR_M_S; n];
    let precipitation: Vec<f64> = grid
        .positions()
        .iter()
        .map(|p| PLACEHOLDER_PRECIPITATION_M_S * (1.0 + 0.5 * p.z))
        .collect();
    for radius_m in [1.5e6, EARTH_RADIUS_M, 1.3e7] {
        let mut h: Vec<f64> = grid
            .positions()
            .iter()
            .map(|p| {
                if p.dot(centre) > 0.6_f64.cos() {
                    10.0
                } else {
                    -1000.0
                }
            })
            .collect();
        let inputs = ErosionInputs {
            grid,
            radius_m,
            sea_level_m: 0.0,
            uplift_rate_m_s: &uplift,
            precipitation_m_s: &precipitation,
        };
        let settings = settings(1.0e6, 1);
        let mut steps = 0;
        let mut change_m = f64::INFINITY;
        while change_m > 1e-12 && steps < 5_000 {
            let before = h.clone();
            erode(&inputs, &settings, &mut h, &mut |_| {});
            change_m = h
                .iter()
                .zip(&before)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0, f64::max);
            steps += 1;
        }
        assert!(
            change_m <= 1e-12,
            "R = {radius_m} m: not steady after {steps} steps: {change_m} m"
        );
        let d = routing(grid, &h, 0.0, radius_m);
        // Q: each cell's runoff added to every cell on its way down.
        let mut q = vec![0.0; n];
        for j in 0..n as u32 {
            let runoff = precipitation[j as usize] * grid.area_sr(j) * radius_m * radius_m;
            let mut c = j;
            loop {
                q[c as usize] += runoff;
                let r = d.receiver(c);
                if r == c {
                    break;
                }
                c = r;
            }
        }
        let l = receiver_distances_m(grid, &d, radius_m);
        let ocean: Vec<bool> = h.iter().map(|&v| v <= 0.0).collect();
        let w = water_surface_m(&d, &h, 0.0, &ocean);
        let (mut land, mut worst, mut worst_distance) = (0, 0.0_f64, 0.0_f64);
        for i in 0..n {
            let r = d.receiver(i as u32) as usize;
            if r == i {
                assert_eq!(l[i].to_bits(), 0.0_f64.to_bits(), "base level {i}");
                continue;
            }
            let great_circle_m =
                grid.position(i as u32).dot(grid.position(r as u32)).acos() * radius_m;
            let distance_error = (l[i] - great_circle_m).abs() / radius_m;
            assert!(
                distance_error < 2e-12,
                "cell {i}: {} m for {great_circle_m} m",
                l[i]
            );
            worst_distance = worst_distance.max(distance_error);
            if ocean[i] {
                continue;
            }
            land += 1;
            let slope = (h[i] - w[r]) / great_circle_m;
            let expected = ONE_MM_PER_YEAR_M_S
                / (ErosionSettings::DEFAULT_ERODIBILITY
                    * q[i].powf(ErosionSettings::DEFAULT_DISCHARGE_EXPONENT));
            let error = relative(slope, expected);
            assert!(
                error < 1e-6,
                "R = {radius_m} m, cell {i}: slope {slope} for {expected}"
            );
            worst = worst.max(error);
        }
        assert!(land > 0 && d.pits().is_empty());
        let top = h.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
        println!(
            "island, R = {radius_m} m: steady after {steps} steps, {land} land cells, highest \
             {top:.3} m, worst slope error {worst:.2e}, worst distance error {worst_distance:.2e} rad"
        );
    }
}

#[test]
fn erode_reports_progress_per_step() {
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
    let start = world.field(SURFACE_HEIGHT).expect("heights").values.clone();
    let mut reports = Vec::new();
    let mut h = start.clone();
    erode(&inputs, &settings(1.0e4, 5), &mut h, &mut |f| {
        reports.push(f)
    });
    let expected: Vec<f64> = (1..=5).map(|k| f64::from(k) / 5.0).collect();
    assert_eq!(reports, expected);
    assert_ne!(h, start, "five steps change the surface");

    reports.clear();
    let mut h = start.clone();
    erode(&inputs, &settings(1.0e4, 0), &mut h, &mut |f| {
        reports.push(f)
    });
    assert!(reports.is_empty());
    assert!(
        h.iter()
            .zip(&start)
            .all(|(a, b)| a.to_bits() == b.to_bits())
    );
}

#[test]
fn erode_steps_compose() {
    // Five steps at once equal five single steps, bit for bit.
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
    let start = &world.field(SURFACE_HEIGHT).expect("heights").values;
    let mut at_once = start.clone();
    erode(&inputs, &settings(1.0e5, 5), &mut at_once, &mut |_| {});
    let mut one_by_one = start.clone();
    for _ in 0..5 {
        erode(&inputs, &settings(1.0e5, 1), &mut one_by_one, &mut |_| {});
    }
    assert!(
        at_once
            .iter()
            .zip(&one_by_one)
            .all(|(a, b)| a.to_bits() == b.to_bits())
    );
}

#[test]
fn receiver_distances_by_hand() {
    // The octahedron: every centre distance is π/2. Cell 5 is the ocean,
    // cell 4 a pit enclosed by cells 0–3 (the example of the drainage
    // module, without the fill).
    let grid = SphereGrid::octahedron();
    let height_m = [300.0, 200.0, 400.0, 250.0, 10.0, -100.0];
    let d = Drainage::build(&grid, &height_m, 0.0, EARTH_RADIUS_M);
    let l = receiver_distances_m(&grid, &d, EARTH_RADIUS_M);
    assert_eq!(l.len(), 6);
    for i in 0..6_u32 {
        let r = d.receiver(i);
        if r == i {
            assert_eq!(l[i as usize].to_bits(), 0.0_f64.to_bits(), "base level {i}");
        } else {
            let k = grid
                .neighbors(i)
                .iter()
                .position(|&j| j == r)
                .expect("a neighbour");
            let stored = grid.centre_distances_rad(i)[k] * EARTH_RADIUS_M;
            assert_eq!(l[i as usize].to_bits(), stored.to_bits(), "cell {i}");
            assert!(
                relative(l[i as usize], PI / 2.0 * EARTH_RADIUS_M) < 1e-12,
                "cell {i}"
            );
        }
    }
    assert_eq!(d.base_levels(), [4, 5]);
}

// ---------------------------------------------------------------------------
// Panics.

/// Runs [`incise`] on a three-cell chain (ocean, then two land cells)
/// after `change` has modified it.
fn incise_chain(change: impl FnOnce(&mut ChainParts)) {
    let mut parts = ChainParts {
        receivers: vec![0, 0, 1],
        stack: vec![0, 1, 2],
        distance_m: vec![0.0, 1000.0, 1000.0],
        discharge_m3_s: vec![3.0, 2.0, 1.0],
        uplift_m_s: vec![0.0; 3],
        height_m: vec![-10.0, 10.0, 20.0],
        sea_level_m: 0.0,
        settings: ErosionSettings::default(),
    };
    change(&mut parts);
    let network = Network {
        receivers: &parts.receivers,
        stack: &parts.stack,
        receiver_distance_m: &parts.distance_m,
        discharge_m3_s: &parts.discharge_m3_s,
    };
    incise(
        &network,
        &parts.uplift_m_s,
        parts.sea_level_m,
        &parts.settings,
        &mut parts.height_m,
    );
}

/// The inputs of [`incise`], owned.
struct ChainParts {
    receivers: Vec<u32>,
    stack: Vec<u32>,
    distance_m: Vec<f64>,
    discharge_m3_s: Vec<f64>,
    uplift_m_s: Vec<f64>,
    height_m: Vec<f64>,
    sea_level_m: f64,
    settings: ErosionSettings,
}

#[test]
fn the_unbroken_chain_runs() {
    incise_chain(|_| {});
}

#[test]
#[should_panic(expected = "incise needs one uplift rate per cell")]
fn incise_needs_one_uplift_rate_per_cell() {
    incise_chain(|p| {
        p.uplift_m_s.pop();
    });
}

#[test]
#[should_panic(expected = "incise needs one receiver per cell")]
fn incise_needs_one_receiver_per_cell() {
    incise_chain(|p| p.receivers.push(0));
}

#[test]
#[should_panic(expected = "incise needs a stack of every cell once")]
fn incise_needs_every_cell_once_in_the_stack() {
    incise_chain(|p| p.stack = vec![0, 1, 1]);
}

#[test]
#[should_panic(expected = "incise needs each cell after its receiver in the stack")]
fn incise_needs_each_cell_after_its_receiver() {
    incise_chain(|p| p.stack = vec![0, 2, 1]);
}

#[test]
#[should_panic(expected = "incise needs finite heights")]
fn incise_needs_finite_heights() {
    incise_chain(|p| p.height_m[2] = f64::NAN);
}

#[test]
#[should_panic(expected = "incise needs finite uplift rates")]
fn incise_needs_finite_uplift_rates() {
    incise_chain(|p| p.uplift_m_s[1] = f64::INFINITY);
}

#[test]
#[should_panic(expected = "incise needs a finite sea level")]
fn incise_needs_a_finite_sea_level() {
    incise_chain(|p| p.sea_level_m = f64::NAN);
}

#[test]
#[should_panic(expected = "incise needs finite, non-negative discharges")]
fn incise_needs_non_negative_discharges() {
    incise_chain(|p| p.discharge_m3_s[2] = -1.0);
}

#[test]
#[should_panic(expected = "incise needs a finite, positive distance to the receiver")]
fn incise_needs_positive_receiver_distances() {
    incise_chain(|p| p.distance_m[1] = 0.0);
}

#[test]
#[should_panic(expected = "incise needs valid settings")]
fn incise_needs_valid_settings() {
    incise_chain(|p| p.settings.erodibility = 0.0);
}

/// Runs one step of [`erode`] on the Earth world of seed 1 after `change`
/// has modified its inputs (radius, uplift, precipitation).
fn erode_world(change: impl FnOnce(&mut f64, &mut Vec<f64>, &mut Vec<f64>)) {
    let world = &worlds()[0];
    let mut radius_m = world.params.radius_m;
    let mut uplift = world.field(UPLIFT_RATE).expect("uplift").values.clone();
    let mut precipitation = rain(&world.grid);
    change(&mut radius_m, &mut uplift, &mut precipitation);
    let inputs = ErosionInputs {
        grid: &world.grid,
        radius_m,
        sea_level_m: world.sea_level_m,
        uplift_rate_m_s: &uplift,
        precipitation_m_s: &precipitation,
    };
    let mut h = world.field(SURFACE_HEIGHT).expect("heights").values.clone();
    erode(&inputs, &settings(1.0e4, 1), &mut h, &mut |_| {});
}

/// Runs [`erode`] with 0 steps on the Earth world of seed 1 after `change`
/// has modified its heights, uplift and sea level: the inputs are checked
/// before the first step.
fn erode_world_without_steps(change: impl FnOnce(&mut Vec<f64>, &mut Vec<f64>, &mut f64)) {
    let world = &worlds()[0];
    let mut h = world.field(SURFACE_HEIGHT).expect("heights").values.clone();
    let mut uplift = world.field(UPLIFT_RATE).expect("uplift").values.clone();
    let mut sea_level_m = world.sea_level_m;
    change(&mut h, &mut uplift, &mut sea_level_m);
    let precipitation = rain(&world.grid);
    let inputs = ErosionInputs {
        grid: &world.grid,
        radius_m: world.params.radius_m,
        sea_level_m,
        uplift_rate_m_s: &uplift,
        precipitation_m_s: &precipitation,
    };
    erode(&inputs, &settings(1.0e4, 0), &mut h, &mut |_| {});
}

#[test]
fn erode_without_steps_accepts_the_unchanged_world() {
    erode_world_without_steps(|_, _, _| {});
}

#[test]
#[should_panic(expected = "erode needs finite heights")]
fn erode_needs_finite_heights() {
    erode_world_without_steps(|h, _, _| h[5] = f64::NAN);
}

#[test]
#[should_panic(expected = "erode needs finite uplift rates")]
fn erode_needs_finite_uplift_rates() {
    erode_world_without_steps(|_, u, _| u[5] = f64::NEG_INFINITY);
}

#[test]
#[should_panic(expected = "erode needs a finite sea level")]
fn erode_needs_a_finite_sea_level() {
    erode_world_without_steps(|_, _, sea| *sea = f64::INFINITY);
}

#[test]
#[should_panic(expected = "erode needs one precipitation rate per cell")]
fn erode_needs_one_precipitation_rate_per_cell() {
    erode_world(|_, _, p| p.truncate(10));
}

#[test]
#[should_panic(expected = "erode needs finite, non-negative precipitation")]
fn erode_needs_non_negative_precipitation() {
    erode_world(|_, _, p| p[17] = -1e-9);
}

#[test]
#[should_panic(expected = "erode needs a finite positive radius")]
fn erode_needs_a_finite_positive_radius() {
    erode_world(|r, _, _| *r = 0.0);
}

#[test]
#[should_panic(expected = "erode needs one height per cell")]
fn erode_needs_one_height_per_cell() {
    let world = &worlds()[0];
    let uplift = vec![0.0; world.grid.len()];
    let precipitation = rain(&world.grid);
    let inputs = ErosionInputs {
        grid: &world.grid,
        radius_m: world.params.radius_m,
        sea_level_m: world.sea_level_m,
        uplift_rate_m_s: &uplift,
        precipitation_m_s: &precipitation,
    };
    erode(&inputs, &settings(1.0e4, 1), &mut [0.0; 3], &mut |_| {});
}

// ---------------------------------------------------------------------------
// Golden value.

#[test]
fn golden_fingerprint_of_an_eroded_earth_world() {
    // The Earth world of seed 1 at 10,000 cells after 10 steps of the
    // default settings without diffusion (`settings`; the incision of
    // ticket M3-03 alone), with its own uplift and the placeholder
    // precipitation: FNV-1a over the bits of the heights (f64,
    // little-endian). A change means that the same world erodes
    // differently. Recorded after the other tests passed and the images of
    // out/M3-03 were viewed (ticket M3-03, Result); unchanged by ticket
    // M3-04, whose default step diffuses as well
    // (`diffusion_tests::golden_fingerprint_of_an_eroded_and_diffused_earth_world`).
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
    let mut h = world.field(SURFACE_HEIGHT).expect("heights").values.clone();
    erode(&inputs, &settings(1.0e4, 10), &mut h, &mut |_| {});
    let bytes: Vec<u8> = h.iter().flat_map(|v| v.to_bits().to_le_bytes()).collect();
    let fingerprint = fnv1a64(&bytes);
    let top = h.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
    println!("fingerprint {fingerprint:#018X}, highest {top:.3} m");
    assert_eq!(fingerprint, GOLDEN_FINGERPRINT);
}

/// The fingerprint of `golden_fingerprint_of_an_eroded_earth_world`.
const GOLDEN_FINGERPRINT: u64 = 0x73A6_193E_F4E0_5D72;

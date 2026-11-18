//! Tests of the sea-surface temperature (ticket M4-05).
//!
//! Reference values are closed forms worked out in the tests (the Ekman
//! transport, the upwelling of a straight coast, the conservation of heat),
//! the energy balance of ticket M4-02 and North & Coakley's observed
//! sea-level temperature (NCC81 Table 1). The criteria of the ticket:
//!
//! - criterion 1: `earth_west_coasts_are_colder_and_eastern_high_latitudes_warmer`;
//! - criterion 2: `the_advection_conserves_heat_in_a_closed_basin`,
//!   `the_budget_closes_in_every_substep`,
//!   `upwelling_cools_its_cell_and_returns_the_heat_by_area`;
//! - criterion 3: `earth_stays_calibrated`,
//!   `earth_converges_to_a_periodic_state`,
//!   `without_currents_it_is_the_energy_balance`,
//!   `extreme_planets_stay_finite`;
//! - criterion 4: `golden_fingerprint_of_earth_sea_surface`.
//!
//! Ticket M4-09 (the multigrid): `a_substep_matches_the_jacobi_solver`.
//!
//! Ticket M4-15 (the balance step): criterion 1,
//! `slow_sea_surfaces_reach_a_stable_periodic_state`; criterion 3,
//! `earth_stays_within_the_tolerance_of_its_periodic_state`.
//!
//! The Earth tests share one run of the Earth-parameter world of seed 1
//! with its M4-02 temperatures, M4-03 winds and M4-04 currents: 40,000
//! cells in release builds, the world of the calibration (about half a
//! minute with the energy balance's own run), and 10,000 cells in debug
//! builds, where the world differs and criterion 1 and the calibration
//! have the thresholds stated in the tests.
#![allow(clippy::disallowed_methods)]

use core::f64::consts::PI;
use std::sync::OnceLock;

use kuniumi_math::Rng;
use kuniumi_math::rng::fnv1a64;
use kuniumi_world::{PlanetParams, Settings};

use super::*;
use crate::stages::climate::circulation::monthly_circulation;
use crate::stages::climate::energy_balance::monthly_temperature;
use crate::stages::climate::ocean::{OceanModel, ocean_currents};
use crate::stages::climate::test_worlds::{
    earth_10k, earth_reference, earth_sea_10k, earth_sea_reference, sea_model,
};
use crate::stages::crust::SURFACE_HEIGHT;

/// One degree in radians.
const DEG: f64 = PI / 180.0;

/// The 10,000-cell grid of seed 1, built once.
fn grid_10k() -> &'static SphereGrid {
    static GRID: OnceLock<SphereGrid> = OnceLock::new();
    GRID.get_or_init(|| SphereGrid::fibonacci(10_000, 1).expect("grid"))
}

/// Latitude and longitude of cell `i` in degrees.
fn lat_lon_deg(grid: &SphereGrid, i: usize) -> (f64, f64) {
    let (lat, lon) = coords::lat_lon(grid.position(i as u32));
    (lat / DEG, lon / DEG)
}

/// Whether ocean cell `i` has a land neighbour within 60° of due west
/// (`west`) and within 60° of due east (`east`), in its local frame (as in
/// the tests of ticket M4-04).
fn coast_sides(grid: &SphereGrid, ocean: &[bool], i: usize) -> (bool, bool) {
    let p = grid.position(i as u32);
    let e = coords::east(p);
    let (mut west, mut east) = (false, false);
    for &j in grid.neighbors(i as u32) {
        if ocean[j as usize] {
            continue;
        }
        let d = grid.position(j) - p;
        let c = d.dot(e) / d.length();
        west |= c < -0.5;
        east |= c > 0.5;
    }
    (west, east)
}

/// The annual mean of twelve monthly fields.
fn annual(months: &[Vec<f64>; MONTHS]) -> Vec<f64> {
    let n = months[0].len();
    (0..n)
        .map(|i| months.iter().map(|m| m[i]).sum::<f64>() / MONTHS as f64)
        .collect()
}

/// Asserts that `actual` and `expected` are the same double, bit for bit.
fn assert_same(actual: f64, expected: f64) {
    assert_eq!(
        actual.to_bits(),
        expected.to_bits(),
        "{actual} against {expected}"
    );
}

/// The CSR entry of cell `j` that points back to `i`.
fn reverse_entry(grid: &SphereGrid, offsets: &[usize], i: usize, j: usize) -> usize {
    let s = grid
        .neighbors(j as u32)
        .iter()
        .position(|&b| b as usize == i)
        .expect("symmetric neighbours");
    offsets[j] + s
}

/// The CSR offsets of a grid.
fn offsets(grid: &SphereGrid) -> Vec<usize> {
    let mut out = vec![0];
    for i in 0..grid.len() {
        out.push(out[i] + grid.neighbors(i as u32).len());
    }
    out
}

/// Heights of a mask: −1,000 m in the ocean, +1,000 m on land (sea level 0).
fn heights_of(ocean: &[bool]) -> Vec<f64> {
    ocean
        .iter()
        .map(|&o| if o { -1_000.0 } else { 1_000.0 })
        .collect()
}

/// A zonal stress `τ_x(φ)` (Pa).
fn zonal_stress(grid: &SphereGrid, tau: impl Fn(f64) -> f64) -> WindStress {
    WindStress {
        east_pa: grid.positions().iter().map(|p| tau(p.z.asin())).collect(),
        north_pa: vec![0.0; grid.len()],
    }
}

/// A uniform stress (Pa, east and north).
fn uniform_stress(grid: &SphereGrid, east: f64, north: f64) -> WindStress {
    WindStress {
        east_pa: vec![east; grid.len()],
        north_pa: vec![north; grid.len()],
    }
}

/// No upwelling in any month.
fn no_upwelling(n: usize) -> [Vec<f64>; MONTHS] {
    core::array::from_fn(|_| vec![0.0; n])
}

// ---------------------------------------------------------------------------
// Parameters, the Ekman transport and the upwelling.

#[test]
fn the_parameters() {
    // κ = 0.3, r_E = 1/(2 days), T_fz = −1.8 °C, ρ_w c_w = 1025 × 3990; the
    // energy balance with currents is M4-02's but for A and D⊕, and D scales
    // with the rotation period as there (× 4 at twice Earth's).
    let p = SeaSurfaceParams::earth();
    assert_same(p.upwelling_contrast, 0.3);
    assert_same(p.ekman_friction_rate_per_s, 1.0 / (2.0 * 86_400.0));
    assert!((p.freezing_temperature_k - 271.35).abs() < 1e-12);
    assert_same(p.volumetric_heat_capacity_j_m3_k, 4_089_750.0);
    let earth = PlanetParams::earth();
    let alone = EnergyBalanceParams::earth();
    let with = SeaSurfaceParams::energy_balance(&earth);
    assert_eq!(
        with,
        EnergyBalanceParams {
            olr_at_freezing_w_m2: SeaSurfaceParams::EARTH_OLR_AT_FREEZING_W_M2,
            diffusion_w_m2_k: SeaSurfaceParams::EARTH_DIFFUSION_W_M2_K,
            ..alone
        }
    );
    let slow = PlanetParams {
        rotation_period_s: 2.0 * earth.rotation_period_s,
        ..earth.clone()
    };
    assert_same(
        SeaSurfaceParams::energy_balance(&slow).diffusion_w_m2_k,
        4.0 * SeaSurfaceParams::EARTH_DIFFUSION_W_M2_K,
    );
    // The energy balance alone keeps its own values.
    assert_same(
        EnergyBalanceParams::for_planet(&slow).diffusion_w_m2_k,
        4.0 * 0.583,
    );
}

#[test]
fn the_ekman_transport_is_right_of_the_stress() {
    // M = (r_E T − f k × T) / (r_E² + f²), T = τ/ρ₀, checked against the
    // same balance written with 3-D vectors (k the radial unit vector, the
    // stress a tangent vector) at every cell, and its limits: away from the
    // equator it is T/|f| at 90° to the right of the stress in the north
    // and to the left in the south, turned back towards the stress by
    // atan(r_E/|f|) (at 45°, f = 1.0313e-4 s⁻¹: 3.21°, and
    // |M| = 0.1/1025 / √(r_E² + f²) = 0.94455 m²/s for 0.1 Pa); at the
    // equator downwind, T/r_E = 0.1/1025 × 172,800 = 16.859 m²/s.
    let g = grid_10k();
    let params = SeaSurfaceParams::earth();
    let ocean = OceanParams::earth();
    let stress = WindStress {
        east_pa: g.positions().iter().map(|p| 0.1 * (1.0 + p.x)).collect(),
        north_pa: g.positions().iter().map(|p| 0.05 * p.y).collect(),
    };
    let (east, north) = ekman_transport_m2_s(&params, &ocean, g, &stress);
    let r = params.ekman_friction_rate_per_s;
    for i in 0..g.len() {
        let p = g.position(i as u32);
        let t = (coords::east(p) * stress.east_pa[i] + coords::north(p) * stress.north_pa[i])
            / ocean.seawater_density_kg_m3;
        let f = 2.0 * ocean.rotation_rate_rad_s * p.z;
        let m = (t * r - p.cross(t) * f) / (r * r + f * f);
        let scale = m.length().max(1e-6);
        assert!((east[i] - m.dot(coords::east(p))).abs() <= 1e-12 * scale);
        assert!((north[i] - m.dot(coords::north(p))).abs() <= 1e-12 * scale);
    }
    // The limits under a uniform eastward stress of 0.1 Pa.
    let stress = uniform_stress(g, 0.1, 0.0);
    let (east, north) = ekman_transport_m2_s(&params, &ocean, g, &stress);
    let t = 0.1 / 1_025.0;
    let near = |lat: f64| {
        (0..g.len())
            .min_by(|&a, &b| {
                (lat_lon_deg(g, a).0 - lat)
                    .abs()
                    .total_cmp(&(lat_lon_deg(g, b).0 - lat).abs())
            })
            .expect("a cell")
    };
    let north45 = near(45.0);
    let f45 = 2.0 * ocean.rotation_rate_rad_s * (45.0 * DEG).sin();
    assert!((f45 - 1.0313e-4).abs() < 1e-8);
    assert!((t / (r * r + f45 * f45).sqrt() - 0.94455).abs() < 1e-4);
    // Southward (to the right of an eastward stress) in the north, with a
    // little downwind; northward in the south.
    let speed = east[north45].hypot(north[north45]);
    assert!((speed - 0.94455).abs() < 0.03, "{speed}");
    let turn = (east[north45] / -north[north45]).atan() / DEG;
    assert!(north[north45] < 0.0 && (turn - 3.21).abs() < 0.2, "{turn}");
    let south45 = near(-45.0);
    assert!(north[south45] > 0.0 && east[south45] > 0.0);
    assert!((east[south45].hypot(north[south45]) - 0.94455).abs() < 0.03);
    // At the equator, the four equatorial cells of the octahedron.
    assert!((t / r - 16.859).abs() < 1e-3);
    let octahedron = SphereGrid::octahedron();
    let (east, north) = ekman_transport_m2_s(
        &params,
        &ocean,
        &octahedron,
        &uniform_stress(&octahedron, 0.1, 0.0),
    );
    let mut equatorial = 0;
    for (i, p) in octahedron.positions().iter().enumerate() {
        if p.z == 0.0 {
            equatorial += 1;
            assert!((east[i] - t / r).abs() <= 1e-12 * (t / r));
            assert_same(north[i].abs(), 0.0);
        }
    }
    assert_eq!(equatorial, 4);
}

/// The ocean west of a continent between 0° and 60° E from 60° S to 60° N.
fn continent(grid: &SphereGrid) -> Vec<bool> {
    (0..grid.len())
        .map(|i| {
            let (lat, lon) = lat_lon_deg(grid, i);
            !(lon > 0.0 && lon < 60.0 && lat.abs() < 60.0)
        })
        .collect()
}

#[test]
fn a_straight_coast_upwells_its_offshore_transport() {
    // A southward stress of 0.1 Pa along the meridional west coast of a
    // continent (0° E): in the north its Ekman transport, f T/(r_E² + f²)
    // across the coast (westward, offshore), must come from below, so the
    // upwelling of the coastal cells from 15° to 45° N adds up to
    // R ∫ f T/(r_E² + f²) dφ over the latitudes of their stretch of coast
    // (Simpson's rule here), from half a cell's spacing below the first to
    // half a spacing above the last (each cell takes the coast around it;
    // about 4.4 Sv on 10,000 cells); in the south the transport is onshore
    // there and offshore at the east coast (60° E). A northward stress
    // reverses both. No cell without a land neighbour upwells, nor any land
    // cell.
    let g = grid_10k();
    let params = SeaSurfaceParams::earth();
    let ocean_params = OceanParams::earth();
    let ocean = continent(g);
    let r = params.ekman_friction_rate_per_s;
    let t = 0.1 / 1_025.0;
    let crosswind = |phi: f64| {
        let f = 2.0 * ocean_params.rotation_rate_rad_s * phi.sin();
        f.abs() * t / (r * r + f * f)
    };
    let band = |w: &[f64], lat0: f64, lat1: f64, lon0: f64, lon1: f64| -> (f64, f64, f64) {
        // The upwelling of the band's coastal cells, and the latitudes of
        // its stretch of coast: the first and last cell's, widened by half
        // their mean spacing.
        let (mut sum, mut lo, mut hi, mut count) = (0.0, 90.0_f64, -90.0_f64, 0);
        for i in 0..g.len() {
            let (lat, lon) = lat_lon_deg(g, i);
            if !ocean[i] || lat < lat0 || lat > lat1 || lon < lon0 || lon > lon1 {
                continue;
            }
            if g.neighbors(i as u32).iter().any(|&j| !ocean[j as usize]) {
                sum += w[i];
                lo = lo.min(lat);
                hi = hi.max(lat);
                count += 1;
            }
        }
        let half = 0.5 * (hi - lo) / (count - 1) as f64;
        (sum, lo - half, hi + half)
    };
    let integral = |lat0: f64, lat1: f64| {
        let steps = 1_000;
        let h = (lat1 - lat0) * DEG / steps as f64;
        let mut sum = crosswind(lat0 * DEG) + crosswind(lat1 * DEG);
        for k in 1..steps {
            let weight = if k % 2 == 1 { 4.0 } else { 2.0 };
            sum += weight * crosswind(lat0 * DEG + k as f64 * h);
        }
        ocean_params.radius_m * sum * h / 3.0
    };
    let southward = coastal_upwelling_m3_s(
        &params,
        &ocean_params,
        g,
        &ocean,
        &uniform_stress(g, 0.0, -0.1),
    );
    let northward = coastal_upwelling_m3_s(
        &params,
        &ocean_params,
        g,
        &ocean,
        &uniform_stress(g, 0.0, 0.1),
    );
    for w in [&southward, &northward] {
        for i in 0..g.len() {
            let coastal = ocean[i] && g.neighbors(i as u32).iter().any(|&j| !ocean[j as usize]);
            assert!(w[i] >= 0.0 && w[i].is_finite());
            if !coastal {
                assert_eq!(w[i].to_bits(), 0.0_f64.to_bits(), "cell {i}");
            }
        }
    }
    let assert_coast = |w: &[f64], lat0: f64, lat1: f64, lon0: f64, lon1: f64| {
        let (sum, lo, hi) = band(w, lat0, lat1, lon0, lon1);
        let exact = integral(lo, hi);
        println!(
            "{lat0}..{lat1}, {lon0}..{lon1}: {:.4} Sv against {:.4} Sv ({lo:.2}..{hi:.2})",
            sum / 1e6,
            exact / 1e6
        );
        assert!((sum - exact).abs() <= 0.05 * exact);
        exact
    };
    let west = assert_coast(&southward, 15.0, 45.0, -10.0, 1.0);
    assert!((west / 1e6 - 4.4).abs() < 0.3, "{}", west / 1e6);
    assert_coast(&southward, -45.0, -15.0, 59.0, 70.0);
    assert_coast(&northward, -45.0, -15.0, -10.0, 1.0);
    assert_coast(&northward, 15.0, 45.0, 59.0, 70.0);
    // The onshore sides: at most 2 % of the offshore transport (notches of
    // the Voronoi coast).
    for (w, lat0, lat1, lon0, lon1) in [
        (&southward, -45.0, -15.0, -10.0, 1.0),
        (&southward, 15.0, 45.0, 59.0, 70.0),
        (&northward, 15.0, 45.0, -10.0, 1.0),
        (&northward, -45.0, -15.0, 59.0, 70.0),
    ] {
        let (sum, _, _) = band(w, lat0, lat1, lon0, lon1);
        println!(
            "onshore {lat0}..{lat1}, {lon0}..{lon1}: {:.4} Sv",
            sum / 1e6
        );
        assert!(sum <= 0.02 * west);
    }
}

// ---------------------------------------------------------------------------
// Conservation.

/// The closed basin of ticket M4-04 (10–50° N, 40° W to 40° E), on a planet
/// without radiation or diffusion (`A = B = D = 0`), its currents under the
/// zonal stress `−0.1 cos(π (φ − 10°)/40°)` Pa, and no upwelling: the
/// advection alone.
fn closed_basin() -> (Vec<bool>, SeaSurfaceModel, OceanCurrents) {
    let g = grid_10k();
    let ocean: Vec<bool> = (0..g.len())
        .map(|i| {
            let (lat, lon) = lat_lon_deg(g, i);
            lat > 10.0 && lat < 50.0 && lon > -40.0 && lon < 40.0
        })
        .collect();
    let stress = zonal_stress(g, |phi| {
        -0.1 * (PI * (phi - 10.0 * DEG) / (40.0 * DEG)).cos()
    });
    let currents = OceanModel::new(&OceanParams::earth(), g, &ocean).currents(&stress);
    let earth = PlanetParams::earth();
    let silent = EnergyBalanceParams {
        olr_at_freezing_w_m2: 0.0,
        olr_slope_w_m2_k: 0.0,
        diffusion_w_m2_k: 0.0,
        ..EnergyBalanceParams::earth()
    };
    let balance = EnergyBalance::new(&silent, g, &heights_of(&ocean), 0.0, earth.orbital_period_s);
    let model = SeaSurfaceModel::new(
        balance,
        &SeaSurfaceParams::earth(),
        g,
        earth.radius_m,
        &currents.transport_m3_s,
        no_upwelling(g.len()),
    );
    (ocean, model, currents)
}

#[test]
fn the_advection_conserves_heat_in_a_closed_basin() {
    // Criterion 2. The advective flux of an edge is the same from both
    // sides with the opposite sign, bit for bit, so the fluxes of a random
    // temperature field sum to zero over the basin (to rounding, 10⁻¹² of
    // Σ |flux|); a uniform temperature gives each cell fluxes that sum to
    // zero within 10⁻¹² of its largest (its transports do, ticket M4-04);
    // and a year of implicit steps of the advection alone keeps the basin's
    // heat Σ a C T within the bound of the solver's tolerance and every
    // temperature within the start's range (the upwind scheme is monotone),
    // while moving the temperatures by kelvins.
    let g = grid_10k();
    let n = g.len();
    let (ocean, model, currents) = closed_basin();
    let offsets = offsets(g);
    let mut rng = Rng::for_stream(1, "test/sea_surface/basin");
    let random: Vec<f64> = (0..n)
        .map(|i| {
            if ocean[i] {
                rng.range_f64(270.0, 300.0)
            } else {
                280.0
            }
        })
        .collect();
    let flux = model.advective_heat_flux_w(&random);
    let (mut sum, mut magnitude, mut moving) = (0.0, 0.0, 0);
    for i in 0..n {
        for k in offsets[i]..offsets[i + 1] {
            let j = g.neighbors(i as u32)[k - offsets[i]] as usize;
            let back = reverse_entry(g, &offsets, i, j);
            assert_eq!(flux[back].abs().to_bits(), flux[k].abs().to_bits());
            assert!(flux[k] == 0.0 || flux[back].is_sign_negative() != flux[k].is_sign_negative());
            if !(ocean[i] && ocean[j]) {
                assert_same(flux[k].abs(), 0.0);
            }
            sum += flux[k];
            magnitude += flux[k].abs();
            moving += usize::from(flux[k] != 0.0);
        }
    }
    println!("basin: Σ flux {sum:e} of Σ |flux| {magnitude:e}, {moving} entries");
    assert!(moving > 1_000);
    assert!(sum.abs() <= 1e-12 * magnitude);
    let uniform = model.advective_heat_flux_w(&vec![288.0; n]);
    let mut worst: f64 = 0.0;
    for i in 0..n {
        let cell = &uniform[offsets[i]..offsets[i + 1]];
        let largest = cell.iter().fold(0.0_f64, |a, x| a.max(x.abs()));
        let total: f64 = cell.iter().sum();
        if largest > 0.0 {
            worst = worst.max(total.abs() / largest);
        }
    }
    println!("uniform: largest cell sum {worst:e} of its largest flux");
    assert!(worst <= 1e-12);
    // The advection alone, a year of sub-steps (no sunlight, no radiation).
    let areas: Vec<f64> = (0..n).map(|i| g.area_sr(i as u32)).collect();
    let capacity = model.balance().heat_capacity_j_m2_k().to_vec();
    let heat = |t: &[f64]| -> f64 {
        let mut h = 0.0;
        for i in 0..n {
            if ocean[i] {
                h += areas[i] * capacity[i] * t[i];
            }
        }
        h
    };
    let (lo, hi) = (0..n)
        .filter(|&i| ocean[i])
        .fold((f64::MAX, f64::MIN), |(a, b), i| {
            (a.min(random[i]), b.max(random[i]))
        });
    let start = heat(&random);
    let mut t = random.clone();
    let dark = vec![0.0; n];
    let mut drift: f64 = 0.0;
    for month in 0..MONTHS {
        for _ in 0..energy_balance::SUBSTEPS_PER_MONTH {
            let budget = model.substep(month, &mut t, &dark);
            assert!(budget.solve.converged());
            drift = drift.max((heat(&t) - start).abs() / start);
            for i in 0..n {
                if ocean[i] {
                    assert!(t[i] >= lo - 1e-9 && t[i] <= hi + 1e-9, "cell {i}: {}", t[i]);
                } else {
                    assert!((t[i] - 280.0).abs() <= 1e-9);
                }
            }
        }
    }
    let moved = (0..n).fold(0.0_f64, |a, i| a.max((t[i] - random[i]).abs()));
    println!("a year of advection: heat drift {drift:e}, largest change {moved:.2} K");
    assert!(moved > 1.0);
    assert!(drift <= 1e-10);
    // The currents are those of the basin.
    assert!(currents.transport_m3_s.iter().any(|&q| q > 1e6));
}

/// The largest `|storage − (absorbed − emitted)|` (W/m²) of a year of
/// sub-steps of `model` from its warm start, and the most heat upwelling
/// took in a sub-step; every solve converged.
fn year_budget(model: &SeaSurfaceModel, insolation: &[Vec<f64>; MONTHS]) -> (f64, f64) {
    let mut t = model.warm_start(insolation);
    let (mut worst, mut upwelled): (f64, f64) = (0.0, 0.0);
    for (month, ins) in insolation.iter().enumerate() {
        for _ in 0..energy_balance::SUBSTEPS_PER_MONTH {
            let b = model.substep(month, &mut t, ins);
            assert!(b.solve.converged());
            worst = worst.max((b.storage_w_m2 - (b.absorbed_w_m2 - b.emitted_w_m2)).abs());
            upwelled = upwelled.max(b.upwelling_w_m2);
        }
    }
    (worst, upwelled)
}

#[test]
fn the_budget_closes_in_every_substep() {
    // Criterion 2: Earth on 10,000 cells, a year from the warm start; the
    // advection and the upwelling with its return create no heat, so the
    // storage equals the sunlight absorbed less the radiation emitted in
    // every sub-step within 10⁻⁵ W/m², while upwelling takes (and the
    // return gives back) more than 1 W/m² over the planet: without the
    // return the budget would miss by that much.
    let run = earth_10k();
    let model = &earth_sea_10k().model;
    let insolation = monthly_insolation(&PlanetParams::earth(), &run.world.grid);
    let (worst, upwelled) = year_budget(model, &insolation);
    println!("largest imbalance {worst:e} W/m², upwelling up to {upwelled:.3} W/m²");
    assert!(worst <= 1e-5);
    assert!(upwelled > 1.0);
}

#[test]
fn upwelling_cools_its_cell_and_returns_the_heat_by_area() {
    // One sub-step on a planet without sunlight, radiation, diffusion or
    // currents (A = B = D = 0), 2,000 cells, land only at cell 0, a uniform
    // 290 K, and 1 Sv of upwelling in one neighbour `c` of cell 0, written
    // out here: u_c = κ ρ_w c_w W / R² = 0.3 × 1025 × 3990 × 10⁶ / R². With
    // m_i = a_i C_i / Δt and the share s_i = a_i / Σ_ocean a, the step is
    // m_c (T'_c − T₀) = −u_c (T'_c − T_fz) + s_c u_c (T'_c − T_fz), so
    // T'_c = T_fz + m_c (T₀ − T_fz) / (m_c + (1 − s_c) u_c), and every other
    // ocean cell warms by s_i u_c (T'_c − T_fz) / m_i; land keeps 290 K, and
    // the heat of the planet does not change.
    let grid = SphereGrid::fibonacci(2_000, 1).expect("grid");
    let n = grid.len();
    let ocean: Vec<bool> = (0..n).map(|i| i != 0).collect();
    let earth = PlanetParams::earth();
    let silent = EnergyBalanceParams {
        olr_at_freezing_w_m2: 0.0,
        olr_slope_w_m2_k: 0.0,
        diffusion_w_m2_k: 0.0,
        ..EnergyBalanceParams::earth()
    };
    let balance = EnergyBalance::new(
        &silent,
        &grid,
        &heights_of(&ocean),
        0.0,
        earth.orbital_period_s,
    );
    let c = grid.neighbors(0)[0] as usize;
    let mut upwelling = no_upwelling(n);
    upwelling[0][c] = 1e6;
    let entries: usize = (0..n).map(|i| grid.neighbors(i as u32).len()).sum();
    let model = SeaSurfaceModel::new(
        balance,
        &SeaSurfaceParams::earth(),
        &grid,
        earth.radius_m,
        &vec![0.0; entries],
        upwelling,
    );
    let r2 = earth.radius_m * earth.radius_m;
    let u_c = 0.3 * (1025.0 * 3990.0) * 1e6 / r2;
    let t_fz = 273.15 - 1.8;
    let dt = earth.orbital_period_s / 48.0;
    let ocean_capacity = EnergyBalanceParams::earth().ocean_heat_capacity_j_m2_k;
    let areas: Vec<f64> = (0..n).map(|i| grid.area_sr(i as u32)).collect();
    let ocean_area: f64 = areas[1..].iter().sum();
    let mass = |i: usize| areas[i] * ocean_capacity / dt;
    let share = |i: usize| areas[i] / ocean_area;
    let mut t = vec![290.0; n];
    let budget = model.substep(0, &mut t, &vec![0.0; n]);
    let expected_c = t_fz + mass(c) * (290.0 - t_fz) / (mass(c) + (1.0 - share(c)) * u_c);
    println!(
        "cell {c}: {:.9} K against {expected_c:.9} K; budget {budget:?}",
        t[c]
    );
    assert!((t[c] - expected_c).abs() <= 1e-9);
    assert!(t[c] < 290.0 - 0.1);
    let mut heat = 0.0;
    for i in 0..n {
        heat += areas[i] * model.balance().heat_capacity_j_m2_k()[i] * (t[i] - 290.0);
        if i == 0 {
            assert!((t[i] - 290.0).abs() <= 1e-9);
        } else if i != c {
            let expected = 290.0 + share(i) * u_c * (expected_c - t_fz) / mass(i);
            assert!(
                (t[i] - expected).abs() <= 1e-9,
                "cell {i}: {} against {expected}",
                t[i]
            );
        }
    }
    // The heat the cell lost is the heat the others gained: none is lost.
    let moved = (1.0 - share(c)) * u_c * (expected_c - t_fz) * dt;
    println!(
        "heat lost by the cell {:e} W; the planet's change {:e} of it",
        moved / dt * r2,
        heat / moved
    );
    assert!(heat.abs() <= 1e-9 * moved);
    let taken = u_c * (expected_c - t_fz) / areas.iter().sum::<f64>();
    assert!((budget.upwelling_w_m2 - taken).abs() <= 1e-9 * taken);
    assert!(budget.storage_w_m2.abs() <= 1e-9 * taken);
}

#[test]
fn without_currents_it_is_the_energy_balance() {
    // With no transport and no upwelling the sea surface is the energy
    // balance of ticket M4-02 solved by another method (BiCGSTAB for the
    // conjugate gradient): Earth's world on 2,000 cells, the same years and
    // the same periodic state within 10⁻⁴ K.
    let params = PlanetParams::earth();
    let settings = Settings {
        cells: 2_000,
        allow_partial: true,
        ..Settings::default()
    };
    let world = crate::generate_through_crust(&params, &settings, 1).expect("a partial world");
    let heights = &world.field(SURFACE_HEIGHT).expect("heights").values;
    let balance = EnergyBalance::new(
        &EnergyBalanceParams::earth(),
        &world.grid,
        heights,
        world.sea_level_m,
        params.orbital_period_s,
    );
    let insolation = monthly_insolation(&params, &world.grid);
    let alone = balance.periodic_state(&insolation);
    let entries: usize = (0..world.grid.len())
        .map(|i| world.grid.neighbors(i as u32).len())
        .sum();
    let model = SeaSurfaceModel::new(
        balance,
        &SeaSurfaceParams::earth(),
        &world.grid,
        params.radius_m,
        &vec![0.0; entries],
        no_upwelling(world.grid.len()),
    );
    let coupled = model.periodic_state(&insolation);
    let mut largest: f64 = 0.0;
    for (a, b) in alone.sea_level_k.iter().zip(&coupled.sea_level_k) {
        for (x, y) in a.iter().zip(b) {
            largest = largest.max((x - y).abs());
        }
    }
    println!(
        "{} and {} years, largest difference {largest:e} K",
        alone.years, coupled.years
    );
    assert_eq!(alone.years, coupled.years);
    assert!(largest <= 1e-4);
}

// ---------------------------------------------------------------------------
// Earth.

/// The area-weighted annual mean of `values` over the ocean cells of the
/// latitudes `[lo, hi)`: of the eastern-boundary cells (land within 60° of
/// due east, none within 60° of due west: the west coasts of continents)
/// and of the western-boundary cells (the reverse).
fn coast_means(grid: &SphereGrid, ocean: &[bool], values: &[f64], lo: f64, hi: f64) -> (f64, f64) {
    let mut sums = [[0.0; 2]; 2];
    for i in 0..grid.len() {
        let lat = lat_lon_deg(grid, i).0;
        if !ocean[i] || lat < lo || lat >= hi {
            continue;
        }
        let (west, east) = coast_sides(grid, ocean, i);
        let kind = if east && !west {
            0
        } else if west && !east {
            1
        } else {
            continue;
        };
        let a = grid.area_sr(i as u32);
        sums[kind][0] += a;
        sums[kind][1] += a * values[i];
    }
    (sums[0][1] / sums[0][0], sums[1][1] / sums[1][0])
}

#[test]
fn earth_west_coasts_are_colder_and_eastern_high_latitudes_warmer() {
    // Criterion 1 on the Earth-parameter world of seed 1: the annual-mean
    // sea-surface temperature of the eastern-boundary cells (the west
    // coasts of continents) against the western-boundary cells (the east
    // coasts). In 15–30° N and 15–30° S the west coasts are colder, by at
    // least 0.5 K on 40,000 cells and 0.2 K on 10,000 (prototype: 0.80 and
    // 1.33 K; 0.44 and 1.85 K); in 45–65° N and 45–65° S the eastern sides
    // of the oceans are warmer by at least 1 K (2.02 and 3.64 K; 3.92 and
    // 5.78 K). The energy balance alone, printed for comparison, has no
    // such contrast (prototype, 40,000 cells: +0.07, +0.17, −0.46, −0.37 K).
    let run = earth_reference();
    let sea = earth_sea_reference();
    let grid = &run.world.grid;
    let ocean = sea.model.balance().ocean();
    let coupled = annual(&sea.temperature.sea_level_k);
    let alone = annual(&run.temperature.sea_level_k);
    let subtropics = if cfg!(debug_assertions) { 0.2 } else { 0.5 };
    for (lo, hi, sign, least) in [
        (15.0, 30.0, -1.0, subtropics),
        (-30.0, -15.0, -1.0, subtropics),
        (45.0, 65.0, 1.0, 1.0),
        (-65.0, -45.0, 1.0, 1.0),
    ] {
        let (west_coasts, east_coasts) = coast_means(grid, ocean, &coupled, lo, hi);
        let (w0, e0) = coast_means(grid, ocean, &alone, lo, hi);
        println!(
            "{lo}..{hi}: west coasts {west_coasts:.2} K, east coasts {east_coasts:.2} K \
             (difference {:+.2} K; the energy balance alone {:+.2} K)",
            west_coasts - east_coasts,
            w0 - e0
        );
        assert!(sign * (west_coasts - east_coasts) >= least);
    }
}

/// The global mean of the annual-mean surface temperature, the `P₂`
/// coefficient `5 ⟨T P₂⟩` of the annual-mean sea-level temperature, and the
/// largest departure of a 10° band's zonal mean from North & Coakley's fit
/// `14.9 − 28.0 P₂(sin φ)` °C (NCC81 Table 1).
fn climate_summary(grid: &SphereGrid, t: &Temperature) -> (f64, f64, f64) {
    let sea_level = annual(&t.sea_level_k);
    let surface = annual(&t.surface_k);
    let p2 = |z: f64| 0.5 * (3.0 * z * z - 1.0);
    let (mut area, mut mean, mut coefficient) = (0.0, 0.0, 0.0);
    let mut bands = [[0.0; 3]; 18];
    for i in 0..grid.len() {
        let a = grid.area_sr(i as u32);
        let z = grid.position(i as u32).z;
        area += a;
        mean += a * surface[i];
        coefficient += a * sea_level[i] * p2(z);
        let b = (((lat_lon_deg(grid, i).0 + 90.0) / 10.0) as usize).min(17);
        bands[b][0] += a;
        bands[b][1] += a * sea_level[i];
        bands[b][2] += a * (ZERO_CELSIUS_K + 14.9 - 28.0 * p2(z));
    }
    let worst = bands
        .iter()
        .fold(0.0_f64, |w, b| w.max(((b[1] - b[2]) / b[0]).abs()));
    (mean / area, 5.0 * coefficient / area, worst)
}

#[test]
fn earth_stays_calibrated() {
    // Criterion 3 with the refitted A and D: on the 40,000-cell world of
    // the calibration (release) the P₂ coefficient −28.0 ± 0.5 K, the mean
    // surface temperature 288.0 ± 0.3 K and every 10° band within 5 K of
    // North & Coakley's fit, as the energy balance alone of ticket M4-02.
    // On the 10,000-cell world of debug builds (another world) within 2 K
    // and 0.5 K of the energy balance alone on the same world (prototype:
    // 1.3 K and 0.25 K).
    let run = earth_reference();
    let sea = earth_sea_reference();
    let grid = &run.world.grid;
    let (mean, p2, worst) = climate_summary(grid, &sea.temperature);
    let (mean0, p20, worst0) = climate_summary(grid, &run.temperature);
    println!(
        "with currents: mean surface {mean:.3} K, P₂ {p2:.3} K, worst band {worst:.2} K; \
         alone: {mean0:.3} K, {p20:.3} K, {worst0:.2} K"
    );
    if cfg!(debug_assertions) {
        assert!((p2 - p20).abs() <= 2.0);
        assert!((mean - mean0).abs() <= 0.5);
    } else {
        assert!((p2 + 28.0).abs() <= 0.5);
        assert!((mean - 288.0).abs() <= 0.3);
        assert!(worst <= 5.0);
    }
}

#[test]
fn earth_converges_to_a_periodic_state() {
    // Criterion 3: periodic within 40 years (prototype: 12 years on 40,000
    // cells, 9 on 10,000), every temperature finite and positive.
    let t = &earth_sea_reference().temperature;
    println!("{} years, last change {:e} K", t.years, t.change_k);
    assert!(t.converged);
    assert!(t.years <= 40);
    for month in t.sea_level_k.iter().chain(&t.surface_k) {
        assert!(month.iter().all(|v| v.is_finite() && *v > 0.0));
    }
}

#[test]
fn earth_stays_within_the_tolerance_of_its_periodic_state() {
    // Ticket M4-15, criterion 3 (as reworded with the owner): Earth's sea
    // surface (the shared run: 40,000 cells in release builds, 10,000 in
    // debug builds) within 3 × 10⁻⁴ K of the same model run to 10⁻⁷ K
    // (measured 2.11 × 10⁻⁴ and 9.97 × 10⁻⁵ K; the annual correction of
    // ticket M4-02 landed 8.9 × 10⁻⁴ and 3.0 × 10⁻⁴ K away).
    let sea = earth_sea_reference();
    let grid = &earth_reference().world.grid;
    let insolation = monthly_insolation(&PlanetParams::earth(), grid);
    let converged = sea.model.periodic_state_with(&insolation, 1e-7, 400);
    let mut difference = 0.0_f64;
    for (a, b) in sea
        .temperature
        .sea_level_k
        .iter()
        .zip(&converged.sea_level_k)
    {
        for (x, y) in a.iter().zip(b) {
            difference = difference.max((x - y).abs());
        }
    }
    println!(
        "{} cells: {} years; to 10⁻⁷ K in {} years, {difference:.2e} K away",
        grid.len(),
        sea.temperature.years,
        converged.years
    );
    assert!(converged.converged);
    assert!(difference <= 3e-4, "{difference:e}");
}

#[test]
fn slow_sea_surfaces_reach_a_stable_periodic_state() {
    // Ticket M4-15, criterion 1: the sea surface of Earth at 1,282 W/m²
    // (seed 1; 75 years with the annual correction of ticket M4-02) and of
    // 8 h at 0.1 bar with obliquity 60° (seed 1; 18 years), with the winds
    // and currents of their energy balances, on 10,000 cells in release
    // builds and 2,000 in debug builds; and in both builds the latter on
    // 2,000 cells, whose sea surface the balance step without its damping
    // swung between two years by 1.04 K at an upwelling coast in the ramp
    // of the albedo (the years without the step: periodic in 56 years).
    // Each periodic within 40 years, and 30 more years without the step
    // from the start of the last year move no monthly mean by more than
    // 10⁻² K and end with a change of at most the periodic tolerance.
    let cells = if cfg!(debug_assertions) {
        2_000
    } else {
        10_000
    };
    let earth = PlanetParams::earth();
    let thin = PlanetParams {
        rotation_period_s: 28_800.0,
        surface_pressure_pa: 1e4,
        obliquity_rad: 60.0 * DEG,
        ..earth.clone()
    };
    let mut planets = vec![
        (
            "Earth at 1,282 W/m²",
            PlanetParams {
                stellar_flux_w_m2: 1_282.0,
                ..earth.clone()
            },
            cells,
        ),
        ("8 h, 0.1 bar, 60°", thin.clone(), cells),
    ];
    if cells != 2_000 {
        planets.push(("8 h, 0.1 bar, 60°", thin, 2_000));
    }
    for (name, params, cells) in &planets {
        let settings = Settings {
            cells: *cells,
            allow_partial: true,
            ..Settings::default()
        };
        let world = crate::generate_through_crust(params, &settings, 1).expect("a partial world");
        let grid = &world.grid;
        let heights = &world.field(SURFACE_HEIGHT).expect("heights").values;
        let temperature = monthly_temperature(params, grid, heights, world.sea_level_m);
        assert!(temperature.converged, "{name}: the energy balance");
        let model = sea_model(params, grid, heights, world.sea_level_m, &temperature);
        let insolation = monthly_insolation(params, grid);
        let (t, start) = energy_balance::periodic_run(
            &model,
            &insolation,
            None,
            true,
            PERIODIC_TOLERANCE_K,
            MAX_YEARS,
        );
        let (bare, _) =
            energy_balance::periodic_run(&model, &insolation, Some(start), false, 0.0, 30);
        let mut moved = 0.0_f64;
        for (a, b) in t.sea_level_k.iter().zip(&bare.sea_level_k) {
            for (x, y) in a.iter().zip(b) {
                moved = moved.max((x - y).abs());
            }
        }
        println!(
            "{name} ({cells} cells): {} years, change {:.2e} K; 30 years without the step move \
             it by {moved:.2e} K, the last by {:.2e} K",
            t.years, t.change_k, bare.change_k
        );
        assert!(t.converged && t.years <= 40, "{name}: {} years", t.years);
        for month in t.sea_level_k.iter().chain(&t.surface_k) {
            assert!(month.iter().all(|v| v.is_finite() && *v > 0.0), "{name}");
        }
        assert!(moved <= 1e-2, "{name}: {moved:e}");
        assert!(
            bare.change_k <= PERIODIC_TOLERANCE_K,
            "{name}: {:e}",
            bare.change_k
        );
    }
}

#[test]
fn extreme_planets_stay_finite() {
    // Planets at the corners of the v1 range on 2,000 cells, each with the
    // temperatures of its energy balance, its winds and its currents: 8 h
    // at 0.1 bar (tropical winds of tens of m/s), 240 h at 10 bar, a water
    // world (no coast, so no upwelling), a world without ocean (nothing to
    // carry or upwell): every temperature finite and positive, the upwelling
    // finite and non-negative, and the run periodic or reported unconverged.
    let earth = PlanetParams::earth();
    let settings = Settings {
        cells: 2_000,
        allow_partial: true,
        ..Settings::default()
    };
    let planets = [
        (
            "8 h, 0.1 bar",
            PlanetParams {
                rotation_period_s: 28_800.0,
                surface_pressure_pa: 1e4,
                ..earth.clone()
            },
            None,
        ),
        (
            "240 h, 10 bar",
            PlanetParams {
                rotation_period_s: 864_000.0,
                surface_pressure_pa: 1e6,
                ..earth.clone()
            },
            None,
        ),
        ("water world", earth.clone(), Some(true)),
        ("no ocean", earth.clone(), Some(false)),
    ];
    for (name, params, all) in planets {
        let world = crate::generate_through_crust(&params, &settings, 1).expect("a partial world");
        let grid = &world.grid;
        let mut heights = world.field(SURFACE_HEIGHT).expect("heights").values.clone();
        if let Some(ocean) = all {
            heights = heights_of(&vec![ocean; grid.len()]);
        }
        let sea_level = if all.is_some() {
            0.0
        } else {
            world.sea_level_m
        };
        let temperature = monthly_temperature(&params, grid, &heights, sea_level);
        let circulation = monthly_circulation(&params, grid, &temperature);
        let currents = ocean_currents(
            &params,
            grid,
            &heights,
            sea_level,
            &circulation,
            &temperature,
        );
        let sea = monthly_sea_surface_temperature(
            &params,
            grid,
            &heights,
            sea_level,
            &circulation,
            &temperature,
            &currents,
        );
        let t = &sea.temperature;
        let (lo, hi) = t
            .surface_k
            .iter()
            .flatten()
            .fold((f64::MAX, f64::MIN), |(a, b), &x| (a.min(x), b.max(x)));
        let upwelled = sea
            .upwelling_m3_s
            .iter()
            .flatten()
            .fold(0.0_f64, |a, &w| a.max(w));
        println!(
            "{name}: {} years (converged {}), surface {lo:.1} to {hi:.1} K, \
             largest upwelling {:.3} Sv",
            t.years,
            t.converged,
            upwelled / 1e6
        );
        for month in t.sea_level_k.iter().chain(&t.surface_k) {
            assert!(month.iter().all(|v| v.is_finite() && *v > 0.0), "{name}");
        }
        for month in &sea.upwelling_m3_s {
            assert!(month.iter().all(|w| w.is_finite() && *w >= 0.0), "{name}");
        }
        assert_eq!(t.converged, t.change_k <= PERIODIC_TOLERANCE_K, "{name}");
        assert!(t.converged || t.years == MAX_YEARS, "{name}");
        if all.is_some() {
            assert_same(upwelled, 0.0);
        }
        // The budget of a year of sub-steps, with the upwelling and its
        // return in balance.
        let model = sea_model(&params, grid, &heights, sea_level, &temperature);
        let (worst, upwelling) = year_budget(&model, &monthly_insolation(&params, grid));
        println!("  largest imbalance {worst:e} W/m², upwelling up to {upwelling:.3} W/m²");
        assert!(worst <= 1e-5, "{name}");
    }
}

// ---------------------------------------------------------------------------
// Panics.

#[test]
#[should_panic(expected = "one transport per CSR entry")]
fn the_model_needs_one_transport_per_entry() {
    let (_, model, _) = closed_basin();
    let g = grid_10k();
    SeaSurfaceModel::new(
        model.balance().clone(),
        &SeaSurfaceParams::earth(),
        g,
        6_371_000.0,
        &[0.0; 3],
        no_upwelling(g.len()),
    );
}

#[test]
#[should_panic(expected = "antisymmetric transports")]
fn the_model_needs_antisymmetric_transports() {
    let (_, model, currents) = closed_basin();
    let g = grid_10k();
    let mut transport = currents.transport_m3_s.clone();
    let k = transport.iter().position(|&q| q > 0.0).expect("a flow");
    transport[k] *= 2.0;
    SeaSurfaceModel::new(
        model.balance().clone(),
        &SeaSurfaceParams::earth(),
        g,
        6_371_000.0,
        &transport,
        no_upwelling(g.len()),
    );
}

#[test]
#[should_panic(expected = "zero on land")]
fn the_model_needs_no_upwelling_on_land() {
    let (ocean, model, currents) = closed_basin();
    let g = grid_10k();
    let mut upwelling = no_upwelling(g.len());
    upwelling[3][ocean.iter().position(|&o| !o).expect("land")] = 1.0;
    SeaSurfaceModel::new(
        model.balance().clone(),
        &SeaSurfaceParams::earth(),
        g,
        6_371_000.0,
        &currents.transport_m3_s,
        upwelling,
    );
}

#[test]
#[should_panic(expected = "κ in [0, 1]")]
fn the_model_needs_a_contrast_between_zero_and_one() {
    let (_, model, currents) = closed_basin();
    let g = grid_10k();
    let params = SeaSurfaceParams {
        upwelling_contrast: 1.5,
        ..SeaSurfaceParams::earth()
    };
    SeaSurfaceModel::new(
        model.balance().clone(),
        &params,
        g,
        6_371_000.0,
        &currents.transport_m3_s,
        no_upwelling(g.len()),
    );
}

#[test]
#[should_panic(expected = "one ocean flag per cell")]
fn the_upwelling_needs_one_flag_per_cell() {
    let g = grid_10k();
    coastal_upwelling_m3_s(
        &SeaSurfaceParams::earth(),
        &OceanParams::earth(),
        g,
        &[true; 5],
        &uniform_stress(g, 0.1, 0.0),
    );
}

#[test]
#[should_panic(expected = "a finite stress")]
fn the_ekman_transport_needs_a_finite_stress() {
    let g = grid_10k();
    let mut stress = uniform_stress(g, 0.1, 0.0);
    stress.north_pa[7] = f64::NAN;
    ekman_transport_m2_s(
        &SeaSurfaceParams::earth(),
        &OceanParams::earth(),
        g,
        &stress,
    );
}

#[test]
#[should_panic(expected = "a month is 0 to 11")]
fn a_month_is_below_twelve() {
    let (_, model, _) = closed_basin();
    let n = grid_10k().len();
    model.substep(12, &mut vec![280.0; n], &vec![0.0; n]);
}

// ---------------------------------------------------------------------------
// Golden fingerprint.

#[test]
fn a_substep_matches_the_jacobi_solver() {
    // Ticket M4-09: a sub-step of July (10,000 cells, Earth, from the
    // periodic state) by BiCGSTAB with the multigrid of the annual-mean
    // matrix against the same step with the Jacobi preconditioner of
    // `sparse::bicgstab` that the model used before, its return vector
    // solved the same way: within 10⁻⁵ K (the bound of the energy balance's
    // `a_substep_matches_the_jacobi_solver`); and the multigrid needs fewer
    // iterations.
    let run = earth_10k();
    let model = &earth_sea_10k().model;
    let insolation = monthly_insolation(&PlanetParams::earth(), &run.world.grid);
    let start = earth_sea_10k().temperature.sea_level_k[6].clone();
    let mut multigrid = start.clone();
    let budget = model.substep(6, &mut multigrid, &insolation[6]);
    // The Jacobi step, as `step` and `solve_with_return` did it.
    let n = model.len();
    let p = model.balance().params();
    let (a, b) = (p.olr_at_freezing_w_m2, p.olr_slope_w_m2_k);
    let dt = model.balance().substep_s();
    let areas = model.balance().areas_sr();
    let capacities = model.balance().heat_capacity_j_m2_k();
    let u = &model.upwelling_coupling[6];
    let t_fz = model.params().freezing_temperature_k;
    let albedos = model.balance().albedos(&start);
    let matrix = model.matrix(&model.step_mass, u);
    let mut z = vec![0.0; n];
    crate::sparse::bicgstab(&matrix, &model.share, 1e-12, SOLVER_MAX_ITERATIONS, &mut z);
    let mut total = 0.0;
    for &x in u {
        total += x;
    }
    let mut rhs = Vec::with_capacity(n);
    for i in 0..n {
        let q = insolation[6][i] * (1.0 - albedos[i]);
        rhs.push(
            areas[i] * (capacities[i] * start[i] / dt + q - a + b * ZERO_CELSIUS_K) + u[i] * t_fz
                - model.share[i] * total * t_fz,
        );
    }
    let mut jacobi = start.clone();
    let solve = crate::sparse::bicgstab(
        &matrix,
        &rhs,
        SOLVER_TOLERANCE,
        SOLVER_MAX_ITERATIONS,
        &mut jacobi,
    );
    let (mut uy, mut denominator) = (0.0, 0.0);
    for i in 0..n {
        uy += u[i] * jacobi[i];
        denominator += model.step_mass[i] * z[i];
    }
    for i in 0..n {
        jacobi[i] += z[i] * (uy / denominator);
    }
    let worst = largest_difference_k(&multigrid, &jacobi);
    println!(
        "multigrid {} iterations, Jacobi {}; largest difference {worst:e} K",
        budget.solve.iterations, solve.iterations
    );
    assert!(budget.solve.converged() && solve.converged());
    assert!(worst <= 1e-5, "{worst:e}");
    assert!(budget.solve.iterations < solve.iterations);
}

#[test]
fn golden_fingerprint_of_earth_sea_surface() {
    // Earth's world of seed 1 on 10,000 cells through the Crust stage, the
    // periodic state of its energy balance (ticket M4-02), its winds and
    // currents, and the periodic state with currents and upwelling: FNV-1a
    // over the bits of the twelve monthly sea-level temperatures, the
    // twelve surface temperatures and the twelve upwellings (f64,
    // little-endian, cells in index order). `monthly_sea_surface_temperature`
    // gives the same bits. A change means that the same world has a
    // different sea surface. Recorded after the other tests passed and the
    // images of out/M4-05 were viewed (ticket M4-05, Result); recorded again
    // when the solver became the multigrid, after the other tests passed
    // and the images of out/M4-09 were viewed (ticket M4-09); and when the
    // balance step replaced the annual correction, after the other tests
    // passed and the images of out/M4-15 were viewed (ticket M4-15: 7 years
    // instead of 9, every monthly mean within 3.6 × 10⁻⁴ K of the old one).
    let run = earth_10k();
    let sea = earth_sea_10k();
    let t = &sea.temperature;
    let upwelling = sea.model.upwelling_m3_s();
    let bytes: Vec<u8> = t
        .sea_level_k
        .iter()
        .chain(&t.surface_k)
        .chain(upwelling)
        .flat_map(|m| m.iter().flat_map(|v| v.to_bits().to_le_bytes()))
        .collect();
    let fingerprint = fnv1a64(&bytes);
    println!("fingerprint {fingerprint:#018X}, {} years", t.years);
    let params = PlanetParams::earth();
    let grid = &run.world.grid;
    let heights = &run.world.field(SURFACE_HEIGHT).expect("heights").values;
    let circulation = monthly_circulation(&params, grid, &run.temperature);
    let currents = ocean_currents(
        &params,
        grid,
        heights,
        run.world.sea_level_m,
        &circulation,
        &run.temperature,
    );
    let again = monthly_sea_surface_temperature(
        &params,
        grid,
        heights,
        run.world.sea_level_m,
        &circulation,
        &run.temperature,
        &currents,
    );
    assert_eq!(&again.temperature, t);
    assert_eq!(&again.upwelling_m3_s, upwelling);
    assert_eq!(again.ocean, sea.model.balance().ocean());
    assert_eq!(fingerprint, GOLDEN_FINGERPRINT);
}

/// The fingerprint of `golden_fingerprint_of_earth_sea_surface`.
const GOLDEN_FINGERPRINT: u64 = 0xDB9F_73DB_B780_2E2E;

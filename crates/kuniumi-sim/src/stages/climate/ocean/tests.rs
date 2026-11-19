//! Tests of the ocean basins and the wind-driven gyres (ticket M4-04).
//!
//! Reference values are closed forms worked in the comments (the bulk
//! formula, the Bernoulli function, the curl of a solid-body stress, the
//! frictional zonal current `τ_x / (ρ₀ H r)` and its integral across a
//! channel, the Sverdrup transport) and the observed facts of the
//! criteria. The criteria of the ticket:
//!
//! - criterion 1: `earth_has_western_intensified_gyres`,
//!   `a_closed_basin_has_a_clockwise_western_intensified_gyre`;
//! - criterion 2: `no_flow_into_land`;
//! - criterion 3: `a_circumpolar_channel_carries_its_frictional_transport`,
//!   `the_curl_converges_on_strained_fields`,
//!   `open_bands_carry_the_frictional_zonal_current`,
//!   `zonal_wind_on_a_water_world_drives_frictional_flow`,
//!   `extreme_planets_stay_finite`;
//! - criterion 4: `golden_fingerprint_of_earth_currents`.
//!
//! Synthetic basins use the 10,000-cell grid of seed 1; the Earth tests use
//! the climate's shared Earth runs ([`crate::stages::climate::test_worlds`]):
//! 40,000 cells in release builds, 10,000 in debug builds, with the same
//! thresholds except in `open_bands_carry_the_frictional_zonal_current`
//! (3 % and 8 %: the exponential fitting slows zonal currents more on the
//! coarser grid).
#![allow(clippy::disallowed_methods)]

use std::sync::OnceLock;

use kuniumi_math::Rng;
use kuniumi_math::rng::fnv1a64;

use super::*;
use crate::stages::climate::circulation::{Cells, CirculationParams, monthly_circulation};
use crate::stages::climate::test_worlds::{Run, earth_10k, earth_reference};
use crate::stages::crust::SURFACE_HEIGHT;

/// One degree in radians.
const DEG: f64 = PI / 180.0;

/// The 10,000-cell grid of seed 1, built once.
fn grid_10k() -> &'static SphereGrid {
    static GRID: OnceLock<SphereGrid> = OnceLock::new();
    GRID.get_or_init(|| SphereGrid::fibonacci(10_000, 1).expect("grid"))
}

/// Latitude and longitude (degrees) of cell `i`.
fn lat_lon_deg(grid: &SphereGrid, i: usize) -> (f64, f64) {
    let p = grid.position(i as u32);
    (p.z.clamp(-1.0, 1.0).asin() / DEG, p.y.atan2(p.x) / DEG)
}

/// A zonal stress `τ_x(φ)` (Pa) with no northward part.
fn zonal_stress(grid: &SphereGrid, tau: impl Fn(f64) -> f64) -> WindStress {
    let n = grid.len();
    WindStress {
        east_pa: (0..n).map(|i| tau(lat_lon_deg(grid, i).0 * DEG)).collect(),
        north_pa: vec![0.0; n],
    }
}

/// Earth's parameters with the friction rate `r` (1/s).
fn with_friction(r: f64) -> OceanParams {
    OceanParams {
        friction_rate_per_s: r,
        ..OceanParams::earth()
    }
}

/// The area-weighted means of `values` (and of `other`) over the ocean
/// cells whose latitude lies in `[lo_deg, hi_deg)`, and the area of land
/// there over the area of the band.
fn band(
    grid: &SphereGrid,
    ocean: &[bool],
    values: &[f64],
    other: &[f64],
    lo_deg: f64,
    hi_deg: f64,
) -> (f64, f64, f64) {
    let (mut a_all, mut a_land, mut a_ocean, mut s1, mut s2) = (0.0, 0.0, 0.0, 0.0, 0.0);
    for i in 0..grid.len() {
        let lat = lat_lon_deg(grid, i).0;
        if lat < lo_deg || lat >= hi_deg {
            continue;
        }
        let a = grid.area_sr(i as u32);
        a_all += a;
        if !ocean[i] {
            a_land += a;
            continue;
        }
        a_ocean += a;
        s1 += a * values[i];
        s2 += a * other[i];
    }
    (s1 / a_ocean, s2 / a_ocean, a_land / a_all)
}

/// Whether ocean cell `i` has a land neighbour within 60° of due west
/// (`west`) and within 60° of due east (`east`), in its local frame.
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

/// Asserts that `x` and `y` have the same magnitude, bit for bit, and
/// opposite signs unless they are zero.
fn assert_antisymmetric(x: f64, y: f64) {
    assert_eq!(x.abs().to_bits(), y.abs().to_bits(), "{x} against {y}");
    assert!(
        x.abs().to_bits() == 0 || x.is_sign_negative() != y.is_sign_negative(),
        "{x} against {y}"
    );
}

/// The speed (m/s) of cell `i`.
fn speed(c: &OceanCurrents, i: usize) -> f64 {
    c.current_east_m_s[i].hypot(c.current_north_m_s[i])
}

// ---------------------------------------------------------------------------
// Parameters, labels, the Bernoulli function, the stress and its curl.

#[test]
fn the_parameters_scale_with_the_planet() {
    let earth = OceanParams::earth();
    assert_eq!(earth.friction_rate_per_s.to_bits(), 2.0e-6_f64.to_bits());
    assert_eq!(earth.layer_depth_m.to_bits(), 1_000.0_f64.to_bits());
    assert_eq!(
        earth.seawater_density_kg_m3.to_bits(),
        1_025.0_f64.to_bits()
    );
    assert_eq!(earth.drag_coefficient.to_bits(), 1.2e-3_f64.to_bits());
    assert_eq!(earth.gas_constant_j_kg_k.to_bits(), 287.05_f64.to_bits());
    assert_eq!(earth.radius_m.to_bits(), 6_371_000.0_f64.to_bits());
    // Ω = 2π / 86,164.0905 s = 7.292116e-5 rad/s.
    assert!((earth.rotation_rate_rad_s - 7.292_116e-5).abs() < 1e-11);
    let other = OceanParams::for_planet(&PlanetParams {
        rotation_period_s: 8.0 * 3_600.0,
        radius_m: 3_000_000.0,
        ..PlanetParams::earth()
    });
    assert!((other.rotation_rate_rad_s - 2.0 * PI / 28_800.0).abs() < 1e-18);
    assert_eq!(other.radius_m.to_bits(), 3_000_000.0_f64.to_bits());
    assert_eq!(
        other.friction_rate_per_s.to_bits(),
        earth.friction_rate_per_s.to_bits()
    );
}

#[test]
fn labels_are_connected_components_in_index_order() {
    // On the 10,000-cell grid: ocean in two polar caps (north of 60° N and
    // south of 60° S), a lake of the cells within 5° of (0°, 0°) and the
    // cells within 5° of (0°, 90°); land elsewhere, with an island of the
    // cells within 3° of (75° N, 0°) inside the northern ocean. Then the
    // basins are the four connected sets, labelled in the order of their
    // lowest cell index, and the land masses the main land and the island.
    let g = grid_10k();
    let n = g.len();
    let near = |i: usize, lat: f64, lon: f64, r: f64| {
        let q = coords::from_lat_lon(lat * DEG, lon * DEG);
        coords::great_circle_rad(g.position(i as u32), q) < r * DEG
    };
    let ocean: Vec<bool> = (0..n)
        .map(|i| {
            let lat = lat_lon_deg(g, i).0;
            ((lat.abs() > 60.0) && !near(i, 75.0, 0.0, 3.0))
                || near(i, 0.0, 0.0, 5.0)
                || near(i, 0.0, 90.0, 5.0)
        })
        .collect();
    let b = OceanBasins::new(g, &ocean);
    assert_eq!(b.basins(), 4);
    assert_eq!(b.land_masses(), 2);
    // Every ocean cell has a basin and no land mass, every land cell the
    // other way round; neighbours of one kind share their label.
    for i in 0..n {
        assert_eq!(b.basin(i as u32).is_some(), ocean[i]);
        assert_eq!(b.land_mass(i as u32).is_some(), !ocean[i]);
        for &j in g.neighbors(i as u32) {
            if ocean[i] == ocean[j as usize] {
                assert_eq!(b.basin(i as u32), b.basin(j));
                assert_eq!(b.land_mass(i as u32), b.land_mass(j));
            }
        }
    }
    // Labels in the order of the lowest index of each set.
    let mut first_basin = vec![usize::MAX; 4];
    let mut first_mass = vec![usize::MAX; 2];
    for i in 0..n {
        if let Some(k) = b.basin(i as u32) {
            first_basin[k as usize] = first_basin[k as usize].min(i);
        }
        if let Some(k) = b.land_mass(i as u32) {
            first_mass[k as usize] = first_mass[k as usize].min(i);
        }
    }
    assert!(
        first_basin.windows(2).all(|w| w[0] < w[1]),
        "{first_basin:?}"
    );
    assert!(first_mass.windows(2).all(|w| w[0] < w[1]), "{first_mass:?}");
    // The two caps, the two lakes, the main land and the island.
    let at = |lat: f64, lon: f64| {
        (0..n)
            .min_by(|&p, &q| {
                let x = coords::from_lat_lon(lat * DEG, lon * DEG);
                coords::great_circle_rad(g.position(p as u32), x)
                    .total_cmp(&coords::great_circle_rad(g.position(q as u32), x))
            })
            .expect("a cell") as u32
    };
    let caps = [b.basin(at(89.0, 0.0)), b.basin(at(-89.0, 0.0))];
    let lakes = [b.basin(at(0.0, 0.0)), b.basin(at(0.0, 90.0))];
    let mut all: Vec<u32> = caps
        .iter()
        .chain(&lakes)
        .map(|x| x.expect("ocean"))
        .collect();
    all.sort_unstable();
    assert_eq!(all, vec![0, 1, 2, 3]);
    assert_ne!(b.land_mass(at(75.0, 0.0)), b.land_mass(at(30.0, 0.0)));
    // A world without ocean, and one without land.
    let dry = OceanBasins::new(g, &vec![false; n]);
    assert_eq!((dry.basins(), dry.land_masses()), (0, 1));
    let wet = OceanBasins::new(g, &vec![true; n]);
    assert_eq!((wet.basins(), wet.land_masses()), (1, 0));
}

#[test]
fn the_bernoulli_function() {
    // B(x) = x / (eˣ − 1): B(0) = 1, B(1) = 1/(e − 1) = 0.581976706869326…,
    // B(−x) = B(x) + x, B(x) → x e^{−x} for large x, → −x for large −x,
    // and no overflow or NaN at ±800.
    assert_eq!(bernoulli(0.0).to_bits(), 1.0_f64.to_bits());
    assert!((bernoulli(1.0) - 0.581_976_706_869_326_4).abs() < 1e-15);
    assert!((bernoulli(-1.0) - 1.581_976_706_869_326_4).abs() < 1e-15);
    assert!((bernoulli(50.0) - 50.0 * (-50.0_f64).exp()).abs() < 1e-35);
    assert!((bernoulli(-50.0) - 50.0).abs() < 1e-13);
    assert_eq!(bernoulli(800.0).to_bits(), 0.0_f64.to_bits());
    assert_eq!(bernoulli(-800.0).to_bits(), 800.0_f64.to_bits());
    // Against the closed form with std's exp_m1, from 10⁻⁸ to 30, and the
    // identity, which the scheme's conservation rests on.
    let mut x: f64 = 1e-8;
    while x < 30.0 {
        for x in [x, -x] {
            let exact = x / x.exp_m1();
            let b = bernoulli(x);
            assert!(
                (b - exact).abs() <= 1e-12 * exact,
                "B({x}) = {b} against {exact}"
            );
            assert!((bernoulli(-x) - b - x).abs() <= 1e-14 * (1.0 + x.abs()));
        }
        x *= 1.37;
    }
}

#[test]
fn the_wind_stress_is_the_bulk_formula() {
    // Three cells; month k (1 to 12) blows (k, 0) m/s at cell 0, (3, 4)
    // m/s at cell 1 and nothing at cell 2, over 101,325 Pa at 288 K:
    // ρ_a = 101,325 / (287.05 × 288) = 1.225 650 kg/m³. Then
    // τ = ρ_a C_D |U| U: cell 1 every month 1.225 650 × 1.2e-3 × 5 × (3, 4)
    // = (0.022 062, 0.029 416) Pa; cell 0 the mean of ρ_a C_D k² over the
    // twelve months (Σ k² = 650), ρ_a C_D 650/12 = 0.079 667 Pa east.
    let params = OceanParams::earth();
    let n = 3;
    let circulation = Circulation {
        pressure_pa: core::array::from_fn(|_| vec![101_325.0; n]),
        wind_east_m_s: core::array::from_fn(|k| vec![(k + 1) as f64, 3.0, 0.0]),
        wind_north_m_s: core::array::from_fn(|_| vec![0.0, 4.0, 0.0]),
        thermal_equator_rad: [0.0; MONTHS],
        cells: Cells::new(&CirculationParams::earth()),
    };
    let temperature: [Vec<f64>; MONTHS] = core::array::from_fn(|_| vec![288.0; n]);
    let s = annual_wind_stress(&params, &circulation, &temperature);
    let rho: f64 = 101_325.0 / (287.05 * 288.0);
    assert!((rho - 1.225_650).abs() < 1e-6);
    assert!((s.east_pa[0] - rho * 1.2e-3 * 650.0 / 12.0).abs() < 1e-15);
    assert!((s.east_pa[0] - 0.079_667).abs() < 1e-6);
    assert_eq!(s.north_pa[0].to_bits(), 0.0_f64.to_bits());
    assert!((s.east_pa[1] - 0.022_062).abs() < 1e-6);
    assert!((s.north_pa[1] - 0.029_416).abs() < 1e-6);
    assert_eq!(s.east_pa[2].to_bits(), 0.0_f64.to_bits());
    assert_eq!(s.north_pa[2].to_bits(), 0.0_f64.to_bits());
}

#[test]
fn the_monthly_stress_averages_to_the_annual_bit_for_bit() {
    // Ticket M4-05: `monthly_wind_stress` gives each month's bulk formula,
    // and `annual_wind_stress`, now its mean, is bit for bit what ticket
    // M4-04 computed (each month's ρ_a C_D |U| U added into one sum in
    // month order, then divided by 12), written out here on the Earth run's
    // circulation and temperatures. Month 1 of cell 1 of the three-cell
    // case above: 1.225 650 × 1.2e-3 × 5 × (3, 4) Pa.
    let params = OceanParams::earth();
    let run = earth_10k();
    let grid = &run.world.grid;
    let circulation = monthly_circulation(&PlanetParams::earth(), grid, &run.temperature);
    let t = &run.temperature.sea_level_k;
    let n = grid.len();
    let months = monthly_wind_stress(&params, &circulation, t);
    let annual = annual_wind_stress(&params, &circulation, t);
    let (mut east, mut north) = (vec![0.0; n], vec![0.0; n]);
    for k in 0..MONTHS {
        for i in 0..n {
            let density = circulation.pressure_pa[k][i] / (params.gas_constant_j_kg_k * t[k][i]);
            let (u, v) = (
                circulation.wind_east_m_s[k][i],
                circulation.wind_north_m_s[k][i],
            );
            let factor = density * params.drag_coefficient * (u * u + v * v).sqrt();
            assert_eq!(months[k].east_pa[i].to_bits(), (factor * u).to_bits());
            assert_eq!(months[k].north_pa[i].to_bits(), (factor * v).to_bits());
            east[i] += factor * u;
            north[i] += factor * v;
        }
    }
    for i in 0..n {
        assert_eq!(annual.east_pa[i].to_bits(), (east[i] / 12.0).to_bits());
        assert_eq!(annual.north_pa[i].to_bits(), (north[i] / 12.0).to_bits());
    }
    let circulation = Circulation {
        pressure_pa: core::array::from_fn(|_| vec![101_325.0; 2]),
        wind_east_m_s: core::array::from_fn(|k| vec![(k + 1) as f64, 3.0]),
        wind_north_m_s: core::array::from_fn(|_| vec![0.0, 4.0]),
        thermal_equator_rad: [0.0; MONTHS],
        cells: Cells::new(&CirculationParams::earth()),
    };
    let temperature: [Vec<f64>; MONTHS] = core::array::from_fn(|_| vec![288.0; 2]);
    let months = monthly_wind_stress(&params, &circulation, &temperature);
    assert!((months[0].east_pa[1] - 0.022_062).abs() < 1e-6);
    assert!((months[0].north_pa[1] - 0.029_416).abs() < 1e-6);
    assert!((months[11].east_pa[0] - 1.225_650 * 1.2e-3 * 144.0).abs() < 1e-6);
}

#[test]
fn the_curl_satisfies_stokes() {
    // A solid-body stress τ = U cos φ towards the east (U = 0.1 Pa) has the
    // curl (1/(a cos φ)) · (−∂(U cos² φ)/∂φ) = 2U sin φ / a. The cells' curls
    // within 4 % of 2U/a (RMS within 1 %) on 10,000 cells; their sum over
    // the cap north of 30° N times the areas is the circulation around the
    // cap, eastward along 30° N: ∮ τ·dl = U cos 30° · 2π a cos 30°, within
    // 2 % (the cap's boundary is a staircase of cells); over the whole
    // sphere the curls cancel to rounding.
    let g = grid_10k();
    let n = g.len();
    let params = OceanParams::earth();
    let model = OceanModel::new(&params, g, &vec![true; n]);
    let u = 0.1;
    let stress = zonal_stress(g, |phi| u * phi.cos());
    let curl = model.stress_curl_pa_m(&stress);
    let a = params.radius_m;
    let scale = 2.0 * u / a;
    let (mut worst, mut sq) = (0.0_f64, 0.0);
    let (mut cap, mut total, mut magnitude) = (0.0, 0.0, 0.0);
    for (i, &curl) in curl.iter().enumerate() {
        let (lat, _) = lat_lon_deg(g, i);
        let exact = scale * (lat * DEG).sin();
        worst = worst.max((curl - exact).abs());
        sq += (curl - exact) * (curl - exact);
        let area_m2 = g.area_sr(i as u32) * a * a;
        if lat >= 30.0 {
            cap += curl * area_m2;
        }
        total += curl * area_m2;
        magnitude += (curl * area_m2).abs();
    }
    let rms = (sq / n as f64).sqrt();
    println!(
        "curl: largest error {:.4}, RMS {:.4} of 2U/a",
        worst / scale,
        rms / scale
    );
    assert!(worst <= 0.04 * scale, "largest error {}", worst / scale);
    assert!(rms <= 0.01 * scale, "RMS error {}", rms / scale);
    let circulation = u * (30.0 * DEG).cos() * 2.0 * PI * a * (30.0 * DEG).cos();
    println!("cap: {cap} against {circulation}");
    assert!((cap - circulation).abs() <= 0.02 * circulation);
    assert!(total.abs() <= 1e-14 * magnitude, "{total} of {magnitude}");
}

/// The largest and the RMS error of the cells' curls (Pa/m) of the stress
/// `field(p)` (Pa, east and north at the cell centre `p`) on `grid` against
/// `exact(p)`, relative to `peak`.
fn curl_errors(
    grid: &SphereGrid,
    field: &dyn Fn(Vec3) -> (f64, f64),
    exact: &dyn Fn(Vec3) -> f64,
    peak: f64,
) -> (f64, f64) {
    let model = OceanModel::new(&OceanParams::earth(), grid, &vec![true; grid.len()]);
    let (east_pa, north_pa) = grid.positions().iter().map(|&p| field(p)).unzip();
    let curl = model.stress_curl_pa_m(&WindStress { east_pa, north_pa });
    let (mut worst, mut sq) = (0.0_f64, 0.0);
    for (&p, &c) in grid.positions().iter().zip(&curl) {
        let e = c - exact(p);
        worst = worst.max(e.abs());
        sq += e * e;
    }
    (worst / peak, (sq / grid.len() as f64).sqrt() / peak)
}

#[test]
fn the_curl_converges_on_strained_fields() {
    // A solid-body stress has no strain, which hides an error of the edge
    // value (review R1-1 of ticket M4-04): with the stresses of the two
    // cells averaged, the error is of the order of the strain times the
    // offset of the edge's midpoint from the segment between the centres,
    // which does not shrink with the cells (5.3 % RMS and 25 % at most of
    // the peak at 10,000 cells for the first field below, about the same at
    // 40,000). Two strained fields, smooth at the poles too:
    // - the three-cell zonal stress `τ_x = −0.1 cos φ cos 6φ`, whose curl is
    //   `−(1/(a cos φ)) ∂(τ_x cos φ)/∂φ = −(0.1/a)(2 sin φ cos 6φ + 6 cos φ sin 6φ)`,
    //   below 6.1 × 0.1/a in magnitude;
    // - `τ = 0.1 r̂ × ∇(xy)` with the gradient on the unit sphere,
    //   `∇(xy) = (y, x, 0) − 2xy r̂`, whose curl is `0.1 ∇²(xy) / a =
    //   −0.6 xy / a` (a harmonic of degree 2), at most 0.3/a in magnitude.
    // With the stress interpolated to the Voronoi vertices: RMS within 3 %
    // and largest within 15 % of the peak at 10,000 cells, the RMS at
    // 40,000 cells at most 0.6 times that (first order).
    let a = 6_371_000.0;
    let three_cell = |p: Vec3| (three_cell_stress(p.z.clamp(-1.0, 1.0).asin()), 0.0);
    let three_cell_curl = |p: Vec3| {
        let lat = p.z.clamp(-1.0, 1.0).asin();
        -(0.1 / a) * (2.0 * lat.sin() * (6.0 * lat).cos() + 6.0 * lat.cos() * (6.0 * lat).sin())
    };
    let rotated = |p: Vec3| {
        let g = Vec3::new(p.y, p.x, 0.0) - p * (2.0 * p.x * p.y);
        let t = p.cross(g) * 0.1;
        (t.dot(coords::east(p)), t.dot(coords::north(p)))
    };
    let rotated_curl = |p: Vec3| -0.6 * p.x * p.y / a;
    let fine = SphereGrid::fibonacci(40_000, 1).expect("grid");
    type Case<'a> = (
        &'a str,
        &'a dyn Fn(Vec3) -> (f64, f64),
        &'a dyn Fn(Vec3) -> f64,
        f64,
    );
    let cases: [Case; 2] = [
        ("three-cell", &three_cell, &three_cell_curl, 6.1 * 0.1 / a),
        ("rotated gradient of xy", &rotated, &rotated_curl, 0.3 / a),
    ];
    for (name, field, exact, peak) in cases {
        let (worst, rms) = curl_errors(grid_10k(), field, exact, peak);
        let (worst_fine, rms_fine) = curl_errors(&fine, field, exact, peak);
        println!(
            "{name}: largest {worst:.4} and {worst_fine:.4}, RMS {rms:.4} and {rms_fine:.4} of the peak"
        );
        assert!(worst <= 0.15 && rms <= 0.03, "{name}: {worst}, {rms}");
        assert!(rms_fine <= 0.6 * rms, "{name}: {rms_fine} after {rms}");
    }
}

#[test]
fn the_beta_fluxes_are_exact_and_antisymmetric() {
    // G_ij = 2Ω (z_end − z_start) over the Voronoi vertices: the same bits
    // with the opposite sign from both sides of every edge, summing to zero
    // around every cell (up to rounding), and the positive ones of a cell
    // add up to 2Ω (z_max − z_min) over its vertices, which this test
    // recomputes from the cell centres (circumcentres of the fan), away
    // from the poles.
    let g = grid_10k();
    let n = g.len();
    let params = OceanParams::earth();
    let model = OceanModel::new(&params, g, &vec![true; n]);
    let two_omega = 2.0 * params.rotation_rate_rad_s;
    for i in 0..n {
        let ring = g.neighbors(i as u32);
        let d = ring.len();
        let range = model.range(i);
        // A cell around a pole has more than one extreme of z on its ring.
        let polar = g.position(i as u32).z.abs() > (85.0 * DEG).sin();
        let mut sum = 0.0;
        let mut positive = 0.0;
        for (m, k) in range.clone().enumerate() {
            let gij = model.beta_flux(k);
            let j = ring[m] as usize;
            let back = model
                .range(j)
                .find(|&q| model.neighbors[q] as usize == i)
                .expect("symmetric neighbours");
            assert_antisymmetric(gij, model.beta_flux(back));
            sum += gij;
            if gij > 0.0 {
                positive += gij;
            }
        }
        assert!(sum.abs() <= 1e-15 * two_omega, "cell {i}: {sum}");
        let p = g.position(i as u32);
        let zs: Vec<f64> = (0..d)
            .map(|m| {
                let (b, c) = (g.position(ring[m]), g.position(ring[(m + 1) % d]));
                let w = (b - p).cross(c - p);
                w.z / w.length()
            })
            .collect();
        let (lo, hi) = zs
            .iter()
            .fold((f64::MAX, f64::MIN), |(a, b), &z| (a.min(z), b.max(z)));
        if polar {
            continue;
        }
        assert!(
            (positive - two_omega * (hi - lo)).abs() <= 1e-12 * two_omega * (hi - lo),
            "cell {i}: {positive} against {}",
            two_omega * (hi - lo)
        );
    }
}

// ---------------------------------------------------------------------------
// Synthetic oceans.

/// `τ_x = −0.1 cos φ cos 6φ` Pa: easterlies at the equator, westerlies near
/// ±30°, easterlies near ±60°.
fn three_cell_stress(phi: f64) -> f64 {
    -0.1 * phi.cos() * (6.0 * phi).cos()
}

/// The 5° bands' mean `u` and `v` against `τ_x / (ρ₀ H r)` on a water
/// world: the largest difference of the means, the largest `|v|` mean and
/// the largest exact mean.
fn water_world_errors(params: &OceanParams) -> (f64, f64, f64) {
    let g = grid_10k();
    let n = g.len();
    let ocean = vec![true; n];
    let model = OceanModel::new(params, g, &ocean);
    let stress = zonal_stress(g, three_cell_stress);
    let c = model.currents(&stress);
    assert!(c.solves.iter().all(|s| s.converged()), "{:?}", c.solves);
    let exact: Vec<f64> = stress
        .east_pa
        .iter()
        .map(|t| {
            t / (params.seawater_density_kg_m3 * params.layer_depth_m * params.friction_rate_per_s)
        })
        .collect();
    let (mut diff, mut v_most, mut largest) = (0.0_f64, 0.0_f64, 0.0_f64);
    for b in 0..36 {
        let lo = -90.0 + 5.0 * b as f64;
        let (u, v, _) = band(
            g,
            &ocean,
            &c.current_east_m_s,
            &c.current_north_m_s,
            lo,
            lo + 5.0,
        );
        let (e, _, _) = band(g, &ocean, &exact, &exact, lo, lo + 5.0);
        diff = diff.max((u - e).abs());
        v_most = v_most.max(v.abs());
        largest = largest.max(e.abs());
    }
    (diff, v_most, largest)
}

#[test]
fn zonal_wind_on_a_water_world_drives_frictional_flow() {
    // On a planet without land the zonal mean of the momentum balance is
    // r u = τ_x / (ρ₀ H) (no zonal pressure gradient, no net meridional
    // flow across a latitude circle). With r = 10⁻⁵ s⁻¹ the Stommel layer
    // (r/β = 440 km at the equator) is resolved on 10,000 cells: every 5°
    // band's mean u within 3 % of the largest exact mean (prototype: 2.2 %),
    // its mean v below 1 % of it. With Earth's r the exponential fitting
    // slows the zonal flow by up to 17 % at this resolution (ticket M4-04):
    // within 20 %.
    let (diff, v, largest) = water_world_errors(&with_friction(1e-5));
    println!(
        "r = 1e-5: {:.4}, v {:.5} of {largest}",
        diff / largest,
        v / largest
    );
    assert!(diff <= 0.03 * largest);
    assert!(v <= 0.01 * largest);
    let (diff, v, largest) = water_world_errors(&OceanParams::earth());
    println!(
        "Earth's r: {:.4}, v {:.5} of {largest}",
        diff / largest,
        v / largest
    );
    assert!(diff <= 0.20 * largest);
    assert!(v <= 0.01 * largest);
}

/// The channel from 70° S to 30° S between two land caps, with
/// `τ_x = 0.1 |sin(π(φ + 30°)/40°)|` Pa, and the transport between the
/// caps (m³/s, eastward positive).
fn channel_transport(params: &OceanParams) -> f64 {
    let g = grid_10k();
    let n = g.len();
    let ocean: Vec<bool> = (0..n)
        .map(|i| {
            let lat = lat_lon_deg(g, i).0;
            lat > -70.0 && lat < -30.0
        })
        .collect();
    let model = OceanModel::new(params, g, &ocean);
    assert_eq!(model.basins().basins(), 1);
    assert_eq!(model.basins().land_masses(), 2);
    let stress = zonal_stress(g, |phi| {
        0.1 * (PI * (phi + 30.0 * DEG) / (40.0 * DEG)).sin().abs()
    });
    let c = model.currents(&stress);
    assert!(c.solves[0].converged());
    let south = (0..n)
        .min_by(|&p, &q| g.position(p as u32).z.total_cmp(&g.position(q as u32).z))
        .expect("a cell");
    let north = (0..n)
        .max_by(|&p, &q| g.position(p as u32).z.total_cmp(&g.position(q as u32).z))
        .expect("a cell");
    let value = |cell: usize| {
        let m = model.basins().land_mass(cell as u32).expect("land");
        c.coast_streamfunction_m2_s[0]
            .iter()
            .find(|(k, _)| *k == m)
            .expect("a coast of the basin")
            .1
    };
    // u = −∂ψ/∂y: an eastward current has ψ falling northwards.
    params.layer_depth_m * (value(south) - value(north))
}

#[test]
fn a_circumpolar_channel_carries_its_frictional_transport() {
    // The zonal mean of the momentum balance in a channel that is ocean all
    // around: r u = τ_x / (ρ₀ H), so the transport between the coasts is
    // ∫ τ_x / (ρ₀ r) dy = (0.1 a / (ρ₀ r)) · 40° · (2/π) (the mean of
    // |sin| over a half period): 27.62 Sv for r = 10⁻⁵ s⁻¹, 138.1 Sv for
    // Earth's r = 2 × 10⁻⁶ s⁻¹. It needs the island constants: with ψ = 0 on
    // both coasts there is none (the prototype's current then reversed next
    // to both coasts). Within 3 % for r = 10⁻⁵ (prototype: 0.1 %), within
    // 10 % with Earth's r (prototype: −7.3 %, the exponential fitting at
    // this resolution).
    let exact = |r: f64| 0.1 * 6_371_000.0 / (1_025.0 * r) * (40.0 * DEG) * (2.0 / PI);
    assert!((exact(1e-5) / 1e6 - 27.62).abs() < 0.01);
    assert!((exact(2e-6) / 1e6 - 138.1).abs() < 0.05);
    for (params, tolerance) in [(with_friction(1e-5), 0.03), (OceanParams::earth(), 0.10)] {
        let transport = channel_transport(&params);
        let e = exact(params.friction_rate_per_s);
        println!(
            "r = {}: {:.2} Sv against {:.2}",
            params.friction_rate_per_s,
            transport / 1e6,
            e / 1e6
        );
        assert!((transport - e).abs() <= tolerance * e);
    }
}

/// The closed basin from 10° to 50° (north, or south if `south`), 80° wide
/// (40° W to 40° E), with `τ_x = −0.1 cos(π(|φ| − 10°)/40°)` Pa (easterlies
/// at 10°, westerlies at 50°); the mask and the currents.
fn closed_basin(params: &OceanParams, south: bool) -> (Vec<bool>, OceanCurrents) {
    let g = grid_10k();
    let n = g.len();
    let ocean: Vec<bool> = (0..n)
        .map(|i| {
            let (lat, lon) = lat_lon_deg(g, i);
            let lat = if south { -lat } else { lat };
            lat > 10.0 && lat < 50.0 && lon > -40.0 && lon < 40.0
        })
        .collect();
    let model = OceanModel::new(params, g, &ocean);
    let stress = zonal_stress(g, |phi| {
        -0.1 * (PI * (phi.abs() - 10.0 * DEG) / (40.0 * DEG)).cos()
    });
    let c = model.currents(&stress);
    assert!(c.solves.iter().all(|s| s.converged()));
    (ocean, c)
}

/// The Sverdrup transport (m³/s) of the closed basin at latitude `phi`:
/// `−curl τ · L_x / (ρ₀ β)` with `L_x` = 80° · a cos φ, `β = 2Ω cos φ / a`
/// and, for `τ_x = −τ₀ cos(k(φ − φ₁))`,
/// `curl τ = −(1/(a cos φ)) ∂(τ_x cos φ)/∂φ = −(τ₀/a)(k sin(k(φ − φ₁)) + cos(k(φ − φ₁)) tan φ)`.
fn sverdrup_transport(phi: f64) -> f64 {
    let (a, omega, rho) = (6_371_000.0, 2.0 * PI / 86_164.090_5, 1_025.0);
    let (tau0, k, phi1) = (0.1, PI / (40.0 * DEG), 10.0 * DEG);
    let curl = -(tau0 / a) * (k * (k * (phi - phi1)).sin() + (k * (phi - phi1)).cos() * phi.tan());
    -curl * (80.0 * DEG * a * phi.cos()) / (rho * 2.0 * omega * phi.cos() / a)
}

#[test]
fn a_closed_basin_has_a_clockwise_western_intensified_gyre() {
    // Criterion 1 on an idealised subtropical basin: the curl of the
    // stress is negative, so ψ > 0 (clockwise) with its maximum next to the
    // western coast, a northward boundary current and a southward interior;
    // the southern mirror the other way round.
    let g = grid_10k();
    for south in [false, true] {
        let (ocean, c) = closed_basin(&OceanParams::earth(), south);
        let sign = if south { -1.0 } else { 1.0 };
        let psi: Vec<f64> = c.streamfunction_m2_s.iter().map(|p| sign * p).collect();
        let cells: Vec<usize> = (0..g.len()).filter(|&i| ocean[i]).collect();
        let top = cells
            .iter()
            .copied()
            .max_by(|&p, &q| psi[p].total_cmp(&psi[q]))
            .expect("ocean");
        let positive = cells.iter().filter(|&&i| psi[i] > 0.0).count();
        let lowest = cells.iter().fold(f64::MAX, |a, &i| a.min(psi[i]));
        println!(
            "south {south}: ψ max {:.2} Sv at {:?}, lowest {:.4} of it, positive {positive} of {}",
            psi[top] * 1_000.0 / 1e6,
            lat_lon_deg(g, top),
            lowest / psi[top],
            cells.len()
        );
        assert!(lowest >= -0.01 * psi[top]);
        assert!(positive as f64 >= 0.99 * cells.len() as f64);
        // The maximum in the western 15 % of the 80° width.
        assert!(lat_lon_deg(g, top).1 < -40.0 + 0.15 * 80.0);
        // The fastest current at a western-coast cell.
        let fastest = cells
            .iter()
            .copied()
            .max_by(|&p, &q| speed(&c, p).total_cmp(&speed(&c, q)))
            .expect("ocean");
        let (west, east) = coast_sides(g, &ocean, fastest);
        assert!(west && !east, "fastest at {:?}", lat_lon_deg(g, fastest));
        // The western-coast cells flow polewards... northward in the north,
        // southward in the south; the others the other way.
        let (mut vw, mut aw, mut vo, mut ao) = (0.0, 0.0, 0.0, 0.0);
        for &i in &cells {
            let a = g.area_sr(i as u32);
            let (west, east) = coast_sides(g, &ocean, i);
            if west && !east {
                vw += a * c.current_north_m_s[i];
                aw += a;
            } else {
                vo += a * c.current_north_m_s[i];
                ao += a;
            }
        }
        println!("  western coast v {:.4}, others {:.5}", vw / aw, vo / ao);
        assert!(sign * vw > 0.0 && sign * vo < 0.0);
    }
    // The Sverdrup interior: with little friction (r = 5 × 10⁻⁷ s⁻¹) the
    // largest transport lies between 85 % and 100 % of the largest
    // Sverdrup transport, 26.97 Sv at 28.6° N (by hand; the friction of
    // the meridional curvature of ψ takes the rest: prototype 24.0 Sv).
    let mut sverdrup: f64 = 0.0;
    let mut phi = 10.0 * DEG;
    while phi < 50.0 * DEG {
        sverdrup = sverdrup.max(sverdrup_transport(phi));
        phi += 0.01 * DEG;
    }
    assert!((sverdrup / 1e6 - 26.97).abs() < 0.01, "{}", sverdrup / 1e6);
    let params = with_friction(5e-7);
    let (_, c) = closed_basin(&params, false);
    let top = c.streamfunction_m2_s.iter().fold(0.0_f64, |a, &p| a.max(p)) * params.layer_depth_m;
    println!(
        "Sverdrup {:.2} Sv, model {:.2} Sv",
        sverdrup / 1e6,
        top / 1e6
    );
    assert!(top >= 0.85 * sverdrup && top <= sverdrup);
}

// ---------------------------------------------------------------------------
// Earth.

/// Earth's currents on the reference run (40,000 cells in release, 10,000
/// in debug) or on the 10,000-cell run, with the ocean mask and the stress.
struct EarthOcean {
    ocean: Vec<bool>,
    stress: WindStress,
    model: OceanModel,
    currents: OceanCurrents,
}

fn earth_ocean(run: &Run) -> EarthOcean {
    let params = PlanetParams::earth();
    let grid = &run.world.grid;
    let heights = &run.world.field(SURFACE_HEIGHT).expect("heights").values;
    let ocean: Vec<bool> = heights
        .iter()
        .map(|&h| h <= run.world.sea_level_m)
        .collect();
    let circulation = monthly_circulation(&params, grid, &run.temperature);
    let ocean_params = OceanParams::for_planet(&params);
    let stress = annual_wind_stress(&ocean_params, &circulation, &run.temperature.sea_level_k);
    let model = OceanModel::new(&ocean_params, grid, &ocean);
    let currents = model.currents(&stress);
    EarthOcean {
        ocean,
        stress,
        model,
        currents,
    }
}

fn earth_reference_ocean() -> &'static EarthOcean {
    static OCEAN: OnceLock<EarthOcean> = OnceLock::new();
    OCEAN.get_or_init(|| earth_ocean(earth_reference()))
}

#[test]
fn earth_has_western_intensified_gyres() {
    // Criterion 1 on the Earth-parameter world of seed 1: in 15–40° N the
    // western-boundary cells (land within 60° of due west, none within 60°
    // of due east) carry the northward return flow of a clockwise gyre and
    // the interior (no land among the cell's neighbours or theirs) its
    // southward Sverdrup flow, and the western cells are at least twice as
    // fast as the eastern-boundary ones and as the interior; in 15–40° S
    // the same with north and south swapped (prototype, 40,000 cells:
    // +0.0153 and −0.00033 m/s, speeds 0.0219, 0.0072, 0.0042 m/s; south
    // −0.0126 and +0.00028, 0.0200, 0.0058, 0.0034).
    let e = earth_reference_ocean();
    let g = &earth_reference().world.grid;
    let n = g.len();
    let coastal: Vec<bool> = (0..n)
        .map(|i| e.ocean[i] && g.neighbors(i as u32).iter().any(|&j| !e.ocean[j as usize]))
        .collect();
    let c = &e.currents;
    for (lo, hi, sign) in [(15.0, 40.0, 1.0), (-40.0, -15.0, -1.0)] {
        let mut sums = [[0.0; 3]; 3]; // west, east, interior: area, v, speed
        for i in 0..n {
            let lat = lat_lon_deg(g, i).0;
            if !e.ocean[i] || lat < lo || lat >= hi {
                continue;
            }
            let (west, east) = coast_sides(g, &e.ocean, i);
            let interior =
                !coastal[i] && g.neighbors(i as u32).iter().all(|&j| !coastal[j as usize]);
            let kind = if west && !east {
                0
            } else if east && !west {
                1
            } else if interior {
                2
            } else {
                continue;
            };
            let a = g.area_sr(i as u32);
            sums[kind][0] += a;
            sums[kind][1] += a * c.current_north_m_s[i];
            sums[kind][2] += a * speed(c, i);
        }
        let mean = |k: usize, q: usize| sums[k][q] / sums[k][0];
        println!(
            "{lo}..{hi}: west v {:.5} speed {:.5}; east v {:.5} speed {:.5}; interior v {:.5} speed {:.5}",
            mean(0, 1),
            mean(0, 2),
            mean(1, 1),
            mean(1, 2),
            mean(2, 1),
            mean(2, 2)
        );
        assert!(sign * mean(0, 1) > 0.0);
        assert!(sign * mean(2, 1) < 0.0);
        assert!(mean(0, 2) >= 2.0 * mean(1, 2));
        assert!(mean(0, 2) >= 2.0 * mean(2, 2));
    }
}

#[test]
fn open_bands_carry_the_frictional_zonal_current() {
    // In a 5° band that is ocean all around, r ū = τ̄_x / (ρ₀ H) in the
    // zonal mean: the difference at most 3 % of the largest such mean of
    // these bands on 40,000 cells (measured 1.2 %) and 8 % on 10,000
    // (measured 5.0 %): the exponential fitting slows zonal currents where
    // the boundary layer is not resolved, by up to 17 % on a 10,000-cell
    // water world (`zonal_wind_on_a_water_world_drives_frictional_flow`).
    let e = earth_reference_ocean();
    let g = &earth_reference().world.grid;
    let p = OceanParams::earth();
    let exact: Vec<f64> = e
        .stress
        .east_pa
        .iter()
        .map(|t| t / (p.seawater_density_kg_m3 * p.layer_depth_m * p.friction_rate_per_s))
        .collect();
    let mut open = Vec::new();
    for b in 0..36 {
        let lo = -90.0 + 5.0 * b as f64;
        let (u, _, land) = band(
            g,
            &e.ocean,
            &e.currents.current_east_m_s,
            &exact,
            lo,
            lo + 5.0,
        );
        let (x, _, _) = band(g, &e.ocean, &exact, &exact, lo, lo + 5.0);
        if land == 0.0 {
            open.push((lo, u, x));
        }
    }
    assert!(open.len() >= 2, "{open:?}");
    let largest = open.iter().fold(0.0_f64, |a, o| a.max(o.2.abs()));
    let tolerance = if g.len() >= 40_000 { 0.03 } else { 0.08 };
    for &(lo, u, x) in &open {
        println!("open band {lo}: u {u:.5} against {x:.5}");
        assert!(
            (u - x).abs() <= tolerance * largest,
            "band {lo}: {u} against {x}"
        );
    }
}

/// The checks of `no_flow_into_land` on one ocean: the transports through
/// coast edges and between land cells are exactly 0, those between ocean
/// cells antisymmetric bit for bit, every ocean cell's sum within 10⁻¹² of
/// the largest transport; every basin lists each land mass on its coast
/// once, in order, the reference at 0.
fn check_no_flow_into_land(grid: &SphereGrid, model: &OceanModel, c: &OceanCurrents) {
    let n = grid.len();
    let largest = c.transport_m3_s.iter().fold(0.0_f64, |a, q| a.max(q.abs()));
    assert!(largest > 0.0);
    let basins = model.basins();
    for i in 0..n {
        let range = model.range(i);
        let mut sum = 0.0;
        for k in range {
            let j = model.neighbors[k] as usize;
            let q = c.transport_m3_s[k];
            match (basins.basin(i as u32), basins.basin(j as u32)) {
                (Some(_), Some(_)) => {
                    let back = model
                        .range(j)
                        .find(|&b| model.neighbors[b] as usize == i)
                        .expect("symmetric");
                    assert_antisymmetric(q, c.transport_m3_s[back]);
                }
                _ => assert_eq!(q.to_bits(), 0.0_f64.to_bits(), "cell {i} to {j}"),
            }
            sum += q;
        }
        assert!(sum.abs() <= 1e-12 * largest, "cell {i}: {sum} of {largest}");
        if basins.basin(i as u32).is_none() {
            assert_eq!(c.current_east_m_s[i].to_bits(), 0.0_f64.to_bits());
            assert_eq!(c.current_north_m_s[i].to_bits(), 0.0_f64.to_bits());
            assert_eq!(c.streamfunction_m2_s[i].to_bits(), 0.0_f64.to_bits());
        }
    }
    for (b, coasts) in c.coast_streamfunction_m2_s.iter().enumerate() {
        let mut expected: Vec<u32> = (0..n)
            .filter(|&i| basins.basin(i as u32) == Some(b as u32))
            .flat_map(|i| {
                grid.neighbors(i as u32)
                    .iter()
                    .filter_map(|&j| basins.land_mass(j))
            })
            .collect();
        expected.sort_unstable();
        expected.dedup();
        assert_eq!(coasts.iter().map(|x| x.0).collect::<Vec<_>>(), expected);
        if !coasts.is_empty() {
            assert!(coasts.iter().any(|x| x.1.to_bits() == 0), "{coasts:?}");
        }
    }
}

#[test]
fn no_flow_into_land() {
    // Criterion 2: the transports of the Earth world and the closed basin
    // are exactly zero through every coast; along the straight western
    // coast of the closed basin (15–45°) the current lies within 20° (RMS)
    // of the meridian (prototype: 15.6° in the north, 16.5° in the south;
    // the Voronoi coast is a staircase, which the current follows on the
    // scale of a few cells).
    let e = earth_reference_ocean();
    check_no_flow_into_land(&earth_reference().world.grid, &e.model, &e.currents);
    let g = grid_10k();
    for south in [false, true] {
        let (ocean, c) = closed_basin(&OceanParams::earth(), south);
        let model = OceanModel::new(&OceanParams::earth(), g, &ocean);
        check_no_flow_into_land(g, &model, &c);
        let (mut sq, mut count) = (0.0, 0);
        for i in 0..g.len() {
            let lat = lat_lon_deg(g, i).0.abs();
            if !ocean[i] || !(15.0..=45.0).contains(&lat) {
                continue;
            }
            let (west, _) = coast_sides(g, &ocean, i);
            if west {
                let angle = c.current_east_m_s[i]
                    .abs()
                    .atan2(c.current_north_m_s[i].abs())
                    / DEG;
                sq += angle * angle;
                count += 1;
            }
        }
        let rms = (sq / count as f64).sqrt();
        println!("south {south}: {count} western-coast cells, RMS {rms:.1}°");
        assert!(count >= 10 && rms <= 20.0);
    }
}

#[test]
fn the_row_held_at_zero_is_satisfied() {
    // The reference land mass of the Earth world's largest basin, and the
    // cell held at zero on a water world, have no equation of their own;
    // the full system's columns and forcings sum to zero, so their row is
    // satisfied by the others' solution: its residual is within 10⁻⁶ of the
    // sum of the magnitudes of its terms (the solver stops at 10⁻¹⁰ of the
    // whole right-hand side).
    let e = earth_reference_ocean();
    let model = &e.model;
    let c = &e.currents;
    let stress = &e.stress;
    let circ = model.edge_circulations(stress);
    let basins = model.basins();
    let largest = (0..basins.basins())
        .max_by_key(|&b| {
            (0..model.len())
                .filter(|&i| basins.basin(i as u32) == Some(b))
                .count()
        })
        .expect("a basin");
    let system = &model.systems[largest as usize];
    let reference = system.reference.expect("a reference land mass");
    let (mut residual, mut magnitude) = (0.0, 0.0);
    for i in 0..model.len() {
        if basins.basin(i as u32) != Some(largest) {
            continue;
        }
        for k in model.range(i) {
            let j = model.neighbors[k] as usize;
            if basins.land_mass(j as u32) == Some(reference) {
                // The land side of the edge: c_ji (ψ_i − 0), and the
                // circulation of the edge seen from the land, −circ[k].
                let term = model.reverse_coupling(k) * c.streamfunction_m2_s[i];
                residual += term + circ[k];
                magnitude += term.abs() + circ[k].abs();
            }
        }
    }
    println!("reference row: residual {residual:e} of {magnitude:e}");
    assert!(residual.abs() <= 1e-6 * magnitude);

    let g = grid_10k();
    let water = OceanModel::new(&OceanParams::earth(), g, &vec![true; g.len()]);
    let stress = zonal_stress(g, three_cell_stress);
    let wc = water.currents(&stress);
    let circ = water.edge_circulations(&stress);
    let (mut residual, mut magnitude) = (0.0, 0.0);
    for k in water.range(0) {
        let j = water.neighbors[k] as usize;
        let term = water.coupling(k) * (wc.streamfunction_m2_s[j] - wc.streamfunction_m2_s[0]);
        residual += term - circ[k];
        magnitude += term.abs() + circ[k].abs();
    }
    assert_eq!(wc.streamfunction_m2_s[0].to_bits(), 0.0_f64.to_bits());
    println!("pinned cell: residual {residual:e} of {magnitude:e}");
    assert!(residual.abs() <= 1e-6 * magnitude);
}

#[test]
fn extreme_planets_stay_finite() {
    // Planets at the corners of the v1 range (8 h and 240 h, 1,500 and
    // 13,000 km) over a random land mask (each cell land with probability
    // 0.3), a water world, a world without ocean and one with a single
    // ocean cell, with random stresses of up to 1 Pa: every output finite,
    // a solve per basin.
    let g = SphereGrid::fibonacci(2_000, 2).expect("grid");
    let n = g.len();
    let mut rng = Rng::for_stream(5, "test/ocean/extreme");
    let random: Vec<bool> = (0..n).map(|_| rng.range_f64(0.0, 1.0) > 0.3).collect();
    let single: Vec<bool> = (0..n).map(|i| i == 17).collect();
    let stress = WindStress {
        east_pa: (0..n).map(|_| rng.range_f64(-1.0, 1.0)).collect(),
        north_pa: (0..n).map(|_| rng.range_f64(-1.0, 1.0)).collect(),
    };
    let earth = PlanetParams::earth();
    let planets = [
        PlanetParams {
            rotation_period_s: 8.0 * 3_600.0,
            ..earth.clone()
        },
        PlanetParams {
            rotation_period_s: 240.0 * 3_600.0,
            ..earth.clone()
        },
        PlanetParams {
            radius_m: 1_500_000.0,
            ..earth.clone()
        },
        PlanetParams {
            radius_m: 13_000_000.0,
            ..earth.clone()
        },
    ];
    for planet in &planets {
        let params = OceanParams::for_planet(planet);
        for ocean in [&random, &vec![true; n], &vec![false; n], &single] {
            let model = OceanModel::new(&params, &g, ocean);
            let c = model.currents(&stress);
            assert_eq!(c.solves.len(), model.basins().basins() as usize);
            assert!(c.solves.iter().all(|s| s.converged()), "{:?}", c.solves);
            assert!(
                c.streamfunction_m2_s
                    .iter()
                    .chain(&c.transport_m3_s)
                    .chain(&c.current_east_m_s)
                    .chain(&c.current_north_m_s)
                    .all(|x| x.is_finite())
            );
            assert!(
                c.coast_streamfunction_m2_s
                    .iter()
                    .flatten()
                    .all(|x| x.1.is_finite())
            );
        }
    }
    let dry = OceanModel::new(&OceanParams::earth(), &g, &vec![false; n]).currents(&stress);
    assert!(dry.solves.is_empty() && dry.transport_m3_s.iter().all(|q| q.to_bits() == 0));
}

#[test]
#[should_panic(expected = "one flag per cell")]
fn the_model_needs_one_flag_per_cell() {
    OceanModel::new(&OceanParams::earth(), &SphereGrid::octahedron(), &[true; 5]);
}

#[test]
#[should_panic(expected = "finite positive parameters")]
fn the_model_needs_positive_friction() {
    OceanModel::new(&with_friction(0.0), &SphereGrid::octahedron(), &[true; 6]);
}

#[test]
#[should_panic(expected = "one stress per cell")]
fn currents_need_one_stress_per_cell() {
    let model = OceanModel::new(&OceanParams::earth(), &SphereGrid::octahedron(), &[true; 6]);
    model.currents(&WindStress {
        east_pa: vec![0.0; 6],
        north_pa: vec![0.0; 5],
    });
}

#[test]
#[should_panic(expected = "finite stress")]
fn currents_need_a_finite_stress() {
    let model = OceanModel::new(&OceanParams::earth(), &SphereGrid::octahedron(), &[true; 6]);
    let mut east = vec![0.0; 6];
    east[3] = f64::NAN;
    model.currents(&WindStress {
        east_pa: east,
        north_pa: vec![0.0; 6],
    });
}

#[test]
#[should_panic(expected = "finite heights and a finite sea level")]
fn ocean_currents_need_a_finite_sea_level() {
    let run = earth_10k();
    let heights = &run.world.field(SURFACE_HEIGHT).expect("heights").values;
    let circulation =
        monthly_circulation(&PlanetParams::earth(), &run.world.grid, &run.temperature);
    ocean_currents(
        &PlanetParams::earth(),
        &run.world.grid,
        heights,
        f64::NAN,
        &circulation,
        &run.temperature,
    );
}

#[test]
#[should_panic(expected = "one height per cell")]
fn ocean_currents_need_one_height_per_cell() {
    let run = earth_10k();
    let circulation =
        monthly_circulation(&PlanetParams::earth(), &run.world.grid, &run.temperature);
    ocean_currents(
        &PlanetParams::earth(),
        &run.world.grid,
        &[0.0; 3],
        0.0,
        &circulation,
        &run.temperature,
    );
}

// ---------------------------------------------------------------------------
// Golden fingerprint.

#[test]
fn golden_fingerprint_of_earth_currents() {
    // Earth's world of seed 1 on 10,000 cells through the Crust stage, its
    // periodic temperatures (ticket M4-02), their circulation (ticket M4-03)
    // and the currents of its annual-mean stress: FNV-1a over the bits of
    // the streamfunction, the transports, the eastward and the northward
    // currents (f64, little-endian, cells and CSR entries in index order).
    // A change means that the same world has different currents. Recorded
    // after the other tests passed and the images of out/M4-04 were viewed
    // (ticket M4-04, Result); recorded again when the temperatures' solver
    // became the multigrid (ticket M4-09), when their periodic state came
    // from the balance step (ticket M4-15) and when the circulation's
    // thermal anomaly was smoothed over 1,000 km (ticket M4-17, after the
    // other tests passed and the images of out/M4-17 were viewed).
    let run = earth_10k();
    let heights = &run.world.field(SURFACE_HEIGHT).expect("heights").values;
    let circulation =
        monthly_circulation(&PlanetParams::earth(), &run.world.grid, &run.temperature);
    let c = ocean_currents(
        &PlanetParams::earth(),
        &run.world.grid,
        heights,
        run.world.sea_level_m,
        &circulation,
        &run.temperature,
    );
    let bytes: Vec<u8> = c
        .streamfunction_m2_s
        .iter()
        .chain(&c.transport_m3_s)
        .chain(&c.current_east_m_s)
        .chain(&c.current_north_m_s)
        .flat_map(|v| v.to_bits().to_le_bytes())
        .collect();
    let fingerprint = fnv1a64(&bytes);
    println!("fingerprint {fingerprint:#018X}");
    assert_eq!(fingerprint, GOLDEN_FINGERPRINT);
}

/// The fingerprint of `golden_fingerprint_of_earth_currents`.
const GOLDEN_FINGERPRINT: u64 = 0x8A1B_EC2A_584E_D38D;

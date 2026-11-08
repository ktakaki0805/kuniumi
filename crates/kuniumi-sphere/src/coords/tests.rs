//! Tests for `coords` (ticket M1-02).
//!
//! Reference values come from the definitions, not from the code under test:
//! the inputs of a round trip, the analytic east/north formulas evaluated with
//! std's `sin` and `cos`, and point pairs constructed at a known angle.
//!
//! - `round_trip_*`: `lat_lon(from_lat_lon(φ, λ))` returns `(φ, λ)` within
//!   1e-15 rad over 10⁵ random points, many of them next to a pole or to the
//!   antimeridian. Run with `--nocapture` to see the measured maxima.
//! - `lat_lon_*`: ranges, poles, the antimeridian and signed zeros.
//! - `frame_*`: `east` and `north` are unit, orthogonal to each other and to
//!   the position, right-handed, equal to the analytic formulas, point towards
//!   growing longitude and latitude, and are continuous across the
//!   antimeridian; their values at the poles.
//! - `great_circle_*`: symmetry, exact zero, right angles, antipodes, and
//!   pairs at 1e-9 rad from each other or from antipodal.
//! - `fingerprints_match_recorded`: recorded hashes of the output bits.

use core::f64::consts::{FRAC_PI_2, PI, TAU};

use kuniumi_math::rng::fnv1a64;
use kuniumi_math::{Rng, Vec3};

use super::{east, from_lat_lon, great_circle_rad, lat_lon, north};

/// The tolerance of the acceptance criteria, in radians or as a length.
const TOL: f64 = 1.0e-15;

/// Number of random round-trip points.
const ROUND_TRIP_POINTS: usize = 100_000;

/// Difference `a − b` of two longitudes, wrapped into (−π, π]: the meridians
/// −π and π are the same.
fn lon_diff(a: f64, b: f64) -> f64 {
    let d = a - b;
    if d > PI {
        d - TAU
    } else if d <= -PI {
        d + TAU
    } else {
        d
    }
}

/// A random sign applied to `x`.
fn random_sign(rng: &mut Rng, x: f64) -> f64 {
    if rng.next_u64() & 1 == 0 { x } else { -x }
}

/// A latitude within 1e-12 rad of a pole (either hemisphere).
fn near_pole_lat(rng: &mut Rng) -> f64 {
    let lat = FRAC_PI_2 - rng.range_f64(0.0, 1.0e-12);
    random_sign(rng, lat)
}

/// A longitude within 1e-12 rad of the antimeridian (either side).
fn near_antimeridian_lon(rng: &mut Rng) -> f64 {
    let lon = PI - rng.range_f64(0.0, 1.0e-12);
    random_sign(rng, lon)
}

/// The round-trip inputs `(φ, λ)`: 40 % uniform, 20 % next to a pole, 20 %
/// next to the antimeridian, 20 % next to both, then fixed edge cases.
fn round_trip_inputs() -> Vec<(f64, f64)> {
    let mut rng = Rng::from_seed(0x0102);
    let mut inputs = Vec::with_capacity(ROUND_TRIP_POINTS + 15);
    let fifth = ROUND_TRIP_POINTS / 5;
    for _ in 0..2 * fifth {
        let lat = rng.range_f64(-FRAC_PI_2, FRAC_PI_2);
        let lon = rng.range_f64(-PI, PI);
        inputs.push((lat, lon));
    }
    for _ in 0..fifth {
        let lat = near_pole_lat(&mut rng);
        let lon = rng.range_f64(-PI, PI);
        inputs.push((lat, lon));
    }
    for _ in 0..fifth {
        let lat = rng.range_f64(-FRAC_PI_2, FRAC_PI_2);
        let lon = near_antimeridian_lon(&mut rng);
        inputs.push((lat, lon));
    }
    for _ in 0..fifth {
        let lat = near_pole_lat(&mut rng);
        let lon = near_antimeridian_lon(&mut rng);
        inputs.push((lat, lon));
    }
    for lat in [0.0, FRAC_PI_2, -FRAC_PI_2] {
        for lon in [0.0, PI, -PI, FRAC_PI_2, -FRAC_PI_2] {
            inputs.push((lat, lon));
        }
    }
    inputs
}

/// A random unit vector, uniform on the sphere (rejection sampling).
fn random_unit(rng: &mut Rng) -> Vec3 {
    loop {
        let v = Vec3::new(
            rng.range_f64(-1.0, 1.0),
            rng.range_f64(-1.0, 1.0),
            rng.range_f64(-1.0, 1.0),
        );
        let len2 = v.length_squared();
        if len2 > 1.0e-6 && len2 <= 1.0 {
            return v.normalized().expect("non-zero");
        }
    }
}

/// A random unit vector perpendicular to the unit vector `a`: `a × u` for a
/// random `u` at least 30° away from `±a`, so that normalising does not
/// magnify the rounding errors of the cross product (`|a · t|` stays at a
/// few 1e-16).
fn random_tangent(rng: &mut Rng, a: Vec3) -> Vec3 {
    loop {
        let t = a.cross(random_unit(rng));
        if t.length_squared() > 0.25 {
            return t.normalized().expect("non-zero");
        }
    }
}

/// Unit vectors next to the poles whose `x² + y²` underflows or is tiny.
fn near_pole_vectors() -> Vec<Vec3> {
    vec![
        Vec3::new(1.0e-170, 0.0, 1.0),
        Vec3::new(0.0, -1.0e-170, 1.0),
        Vec3::new(-3.0e-200, 2.0e-200, -1.0),
        Vec3::new(1.0e-300, 1.0e-300, 1.0),
        Vec3::new(-5.0e-324, 0.0, -1.0),
        Vec3::new(1.0e-12, -1.0e-12, 1.0)
            .normalized()
            .expect("non-zero"),
        Vec3::new(1.0e-150, 1.0e-160, -1.0),
    ]
}

// ---------------------------------------------------------------------------
// Round trip

#[test]
fn round_trip_lat_lon() {
    let inputs = round_trip_inputs();
    assert_eq!(inputs.len(), ROUND_TRIP_POINTS + 15);
    let (mut max_lat, mut max_lon) = (0.0_f64, 0.0_f64);
    let (mut worst_lat, mut worst_lon) = ((0.0, 0.0), (0.0, 0.0));
    for &(lat, lon) in &inputs {
        let (lat2, lon2) = lat_lon(from_lat_lon(lat, lon));
        assert!(
            (-FRAC_PI_2..=FRAC_PI_2).contains(&lat2),
            "latitude {lat2} out of range for input ({lat}, {lon})"
        );
        assert!(
            lon2 > -PI && lon2 <= PI,
            "longitude {lon2} out of (−π, π] for input ({lat}, {lon})"
        );
        let e_lat = (lat2 - lat).abs();
        let e_lon = lon_diff(lon2, lon).abs();
        if e_lat > max_lat {
            max_lat = e_lat;
            worst_lat = (lat, lon);
        }
        if e_lon > max_lon {
            max_lon = e_lon;
            worst_lon = (lat, lon);
        }
    }
    println!(
        "round trip over {} points: max latitude error {max_lat:e} rad at {worst_lat:?}, \
         max longitude error {max_lon:e} rad at {worst_lon:?}",
        inputs.len()
    );
    assert!(
        max_lat <= TOL,
        "latitude error {max_lat:e} at {worst_lat:?}"
    );
    assert!(
        max_lon <= TOL,
        "longitude error {max_lon:e} at {worst_lon:?}"
    );
}

#[test]
fn round_trip_antimeridian_minus_pi_maps_to_pi() {
    // The double −π lies 1.2e-16 east of the real −π, on the antimeridian.
    for lat in [-1.0, 0.0, 0.5, FRAC_PI_2 - 1.0e-13] {
        let (_, lon) = lat_lon(from_lat_lon(lat, -PI));
        assert_eq!(lon.to_bits(), PI.to_bits(), "lat {lat}: {lon}");
        let (_, lon) = lat_lon(from_lat_lon(lat, PI));
        assert_eq!(lon.to_bits(), PI.to_bits(), "lat {lat}: {lon}");
    }
}

// ---------------------------------------------------------------------------
// lat_lon

#[test]
fn lat_lon_ranges_for_random_unit_vectors() {
    let mut rng = Rng::from_seed(0x0103);
    for _ in 0..100_000 {
        let p = random_unit(&mut rng);
        let (lat, lon) = lat_lon(p);
        assert!((-FRAC_PI_2..=FRAC_PI_2).contains(&lat), "{p:?}");
        assert!(lon > -PI && lon <= PI, "{p:?}");
        // Compare with the definition, evaluated with std's functions.
        let lat_ref = p.z.atan2((p.x * p.x + p.y * p.y).sqrt());
        let lon_ref = p.y.atan2(p.x);
        assert!((lat - lat_ref).abs() <= TOL, "{p:?}: {lat} vs {lat_ref}");
        assert!(
            lon_diff(lon, lon_ref).abs() <= TOL,
            "{p:?}: {lon} vs {lon_ref}"
        );
    }
}

#[test]
fn lat_lon_axes() {
    let cases = [
        (Vec3::X, 0.0, 0.0),
        (Vec3::Y, 0.0, FRAC_PI_2),
        (-Vec3::X, 0.0, PI),
        (-Vec3::Y, 0.0, -FRAC_PI_2),
        (Vec3::Z, FRAC_PI_2, 0.0),
        (-Vec3::Z, -FRAC_PI_2, 0.0),
    ];
    for (p, lat, lon) in cases {
        assert_eq!(lat_lon(p), (lat, lon), "{p:?}");
    }
}

#[test]
fn lat_lon_poles_have_longitude_positive_zero() {
    for p in [
        Vec3::new(0.0, 0.0, 1.0),
        Vec3::new(-0.0, 0.0, 1.0),
        Vec3::new(0.0, -0.0, 1.0),
        Vec3::new(-0.0, -0.0, 1.0),
        Vec3::new(0.0, 0.0, -1.0),
        Vec3::new(-0.0, 0.0, -1.0),
        Vec3::new(0.0, -0.0, -1.0),
        Vec3::new(-0.0, -0.0, -1.0),
    ] {
        let (lat, lon) = lat_lon(p);
        assert_eq!(
            lat.to_bits(),
            FRAC_PI_2.copysign(p.z).to_bits(),
            "{p:?}: {lat}"
        );
        assert_eq!(lon.to_bits(), 0.0_f64.to_bits(), "{p:?}: longitude {lon}");
    }
}

#[test]
fn lat_lon_antimeridian() {
    // atan2 returns −π for y = −0 or a tiny negative y with x < 0; the
    // antimeridian has the single longitude π.
    for p in [
        Vec3::new(-1.0, 0.0, 0.0),
        Vec3::new(-1.0, -0.0, 0.0),
        Vec3::new(-1.0, -1.0e-300, 0.0),
        Vec3::new(-1.0, -1.0e-17, 0.0),
        Vec3::new(-1.0e-12, -0.0, 1.0),
    ] {
        assert_eq!(lat_lon(p).1.to_bits(), PI.to_bits(), "{p:?}");
    }
    // Just west of the antimeridian the longitude stays negative.
    let (_, lon) = lat_lon(Vec3::new(-1.0, -1.0e-15, 0.0));
    assert!(lon > -PI && lon < -PI + 1.0e-14, "{lon}");
}

#[test]
fn lat_lon_near_poles() {
    for p in near_pole_vectors() {
        let (lat, lon) = lat_lon(p);
        // Where x² + y² underflows, std's reference gives ±π/2, which is also
        // the correctly rounded latitude there.
        let lat_ref = p.z.atan2((p.x * p.x + p.y * p.y).sqrt());
        assert!((lat - lat_ref).abs() <= TOL, "{p:?}: {lat} vs {lat_ref}");
        let lon_ref = p.y.atan2(p.x);
        if lon_ref > -PI {
            assert!((lon - lon_ref).abs() <= TOL, "{p:?}: {lon} vs {lon_ref}");
        }
    }
}

#[test]
fn lat_lon_and_east_depend_only_on_direction() {
    // Scales of 1e±200 make x² + y² underflow or overflow; the horizontal
    // length must be computed without squaring the raw components.
    let mut rng = Rng::from_seed(0x0104);
    for _ in 0..10_000 {
        let p = random_unit(&mut rng);
        let (lat, lon) = lat_lon(p);
        let e = east(p);
        for scale in [1.0e-200, 1.0e-3, 7.0, 1.0e200] {
            let (lat_s, lon_s) = lat_lon(p * scale);
            assert!((lat_s - lat).abs() <= TOL, "{p:?} × {scale}");
            assert!(lon_diff(lon_s, lon).abs() <= TOL, "{p:?} × {scale}");
            assert!((east(p * scale) - e).length() <= TOL, "{p:?} × {scale}");
        }
    }
}

// ---------------------------------------------------------------------------
// East and north

/// Largest deviation of the frame at `p` from the required properties:
/// unit `east` and `north`, orthogonal to each other and to `p`, and
/// `east × north == p`.
fn frame_deviation(p: Vec3) -> f64 {
    let e = east(p);
    let n = north(p);
    [
        (e.length() - 1.0).abs(),
        (n.length() - 1.0).abs(),
        e.dot(n).abs(),
        e.dot(p).abs(),
        n.dot(p).abs(),
        (e.cross(n) - p).length(),
    ]
    .into_iter()
    .fold(0.0, f64::max)
}

#[test]
fn frame_is_orthonormal_and_right_handed() {
    let mut points: Vec<Vec3> = round_trip_inputs()
        .into_iter()
        .map(|(lat, lon)| from_lat_lon(lat, lon))
        .collect();
    let mut rng = Rng::from_seed(0x0105);
    points.extend((0..100_000).map(|_| random_unit(&mut rng)));
    points.extend([Vec3::Z, -Vec3::Z, Vec3::X, -Vec3::X, Vec3::Y, -Vec3::Y]);
    points.extend(near_pole_vectors());
    let (mut max, mut worst) = (0.0_f64, Vec3::ZERO);
    for &p in &points {
        let d = frame_deviation(p);
        if d > max {
            max = d;
            worst = p;
        }
    }
    println!(
        "frame over {} points: max deviation {max:e} at {worst:?}",
        points.len()
    );
    assert!(max <= TOL, "frame deviation {max:e} at {worst:?}");
}

#[test]
fn frame_matches_analytic_formulas() {
    // east = (−sin λ, cos λ, 0), north = (−sin φ cos λ, −sin φ sin λ, cos φ).
    let (mut max_e, mut max_n) = (0.0_f64, 0.0_f64);
    for (lat, lon) in round_trip_inputs() {
        let p = from_lat_lon(lat, lon);
        let (sin_lat, cos_lat) = (lat.sin(), lat.cos());
        let (sin_lon, cos_lon) = (lon.sin(), lon.cos());
        let e_ref = Vec3::new(-sin_lon, cos_lon, 0.0);
        let n_ref = Vec3::new(-sin_lat * cos_lon, -sin_lat * sin_lon, cos_lat);
        max_e = max_e.max((east(p) - e_ref).length());
        max_n = max_n.max((north(p) - n_ref).length());
    }
    println!("analytic frame: max |Δeast| {max_e:e}, max |Δnorth| {max_n:e}");
    assert!(max_e <= TOL, "east deviates by {max_e:e}");
    assert!(max_n <= TOL, "north deviates by {max_n:e}");
}

#[test]
fn frame_points_towards_growing_longitude_and_latitude() {
    let mut rng = Rng::from_seed(0x0106);
    let step = 1.0e-6;
    for _ in 0..1_000 {
        let lat = rng.range_f64(-1.5, 1.5);
        let lon = rng.range_f64(-PI + 0.01, PI - 0.01);
        let p = from_lat_lon(lat, lon);
        // A step east: longitude grows by step / cos φ, latitude stays.
        let (lat_e, lon_e) = lat_lon((p + east(p) * step).normalized().expect("non-zero"));
        let d_lon = (lon_e - lon) * lat.cos() / step;
        assert!(
            (d_lon - 1.0).abs() < 1.0e-6,
            "east at ({lat}, {lon}): {d_lon}"
        );
        assert!((lat_e - lat).abs() < 1.0e-9, "east at ({lat}, {lon})");
        // A step north: latitude grows by step, longitude stays.
        let (lat_n, lon_n) = lat_lon((p + north(p) * step).normalized().expect("non-zero"));
        let d_lat = (lat_n - lat) / step;
        assert!(
            (d_lat - 1.0).abs() < 1.0e-6,
            "north at ({lat}, {lon}): {d_lat}"
        );
        assert!((lon_n - lon).abs() < 1.0e-9, "north at ({lat}, {lon})");
    }
}

#[test]
fn frame_is_continuous_across_the_antimeridian() {
    let mut rng = Rng::from_seed(0x0107);
    for _ in 0..10_000 {
        let lat = rng.range_f64(-1.5, 1.5);
        // Two points 1e-9 rad of longitude apart, on either side of ±π.
        let west = from_lat_lon(lat, PI - 0.5e-9);
        let east_side = from_lat_lon(lat, -PI + 0.5e-9);
        let de = (east(west) - east(east_side)).length();
        let dn = (north(west) - north(east_side)).length();
        assert!(de < 1.0e-8 && dn < 1.0e-8, "lat {lat}: {de:e}, {dn:e}");
    }
}

#[test]
fn frame_at_the_poles() {
    // The limit of the frame along the meridian of longitude 0.
    assert_eq!(east(Vec3::Z), Vec3::Y);
    assert_eq!(north(Vec3::Z), -Vec3::X);
    assert_eq!(east(-Vec3::Z), Vec3::Y);
    assert_eq!(north(-Vec3::Z), Vec3::X);
    for (lat, pole) in [
        (FRAC_PI_2 - 1.0e-9, Vec3::Z),
        (-FRAC_PI_2 + 1.0e-9, -Vec3::Z),
    ] {
        let near = from_lat_lon(lat, 0.0);
        assert!((east(near) - east(pole)).length() < 1.0e-8);
        assert!((north(near) - north(pole)).length() < 1.0e-8);
    }
}

#[test]
fn frame_on_the_axes() {
    assert_eq!(east(Vec3::X), Vec3::Y);
    assert_eq!(north(Vec3::X), Vec3::Z);
    assert_eq!(east(Vec3::Y), -Vec3::X);
    assert_eq!(north(Vec3::Y), Vec3::Z);
}

// ---------------------------------------------------------------------------
// Great-circle distance

#[test]
fn great_circle_is_symmetric_bitwise() {
    let mut rng = Rng::from_seed(0x0108);
    for _ in 0..100_000 {
        let a = random_unit(&mut rng);
        let b = random_unit(&mut rng);
        assert_eq!(
            great_circle_rad(a, b).to_bits(),
            great_circle_rad(b, a).to_bits(),
            "{a:?}, {b:?}"
        );
    }
}

#[test]
fn great_circle_is_zero_for_identical_inputs() {
    let mut rng = Rng::from_seed(0x0109);
    let mut points: Vec<Vec3> = (0..10_000).map(|_| random_unit(&mut rng)).collect();
    points.extend([Vec3::X, -Vec3::Y, Vec3::Z, Vec3::new(3.0, -4.0, 12.0)]);
    for p in points {
        assert_eq!(great_circle_rad(p, p).to_bits(), 0.0_f64.to_bits(), "{p:?}");
    }
}

#[test]
fn great_circle_right_angles_and_antipodes() {
    let axes = [Vec3::X, Vec3::Y, Vec3::Z, -Vec3::X, -Vec3::Y, -Vec3::Z];
    for &a in &axes {
        for &b in &axes {
            let expected = if a == b {
                0.0
            } else if a == -b {
                PI
            } else {
                FRAC_PI_2
            };
            let d = great_circle_rad(a, b);
            assert!((d - expected).abs() <= TOL, "{a:?}, {b:?}: {d}");
        }
    }
    let mut rng = Rng::from_seed(0x010A);
    let (mut max_right, mut max_anti) = (0.0_f64, 0.0_f64);
    for _ in 0..100_000 {
        let a = random_unit(&mut rng);
        let t = random_tangent(&mut rng, a);
        max_right = max_right.max((great_circle_rad(a, t) - FRAC_PI_2).abs());
        max_anti = max_anti.max((great_circle_rad(a, -a) - PI).abs());
    }
    println!("right angles: max error {max_right:e}; antipodes: max error {max_anti:e}");
    assert!(max_right <= TOL, "right angle error {max_right:e}");
    assert!(max_anti <= TOL, "antipode error {max_anti:e}");
}

#[test]
fn great_circle_tiny_separation() {
    // b = a cos θ + t sin θ lies at angle θ from a, up to the rounding of b's
    // components (about 1e-16 absolute, so 1e-7 relative to θ = 1e-9).
    let theta: f64 = 1.0e-9;
    let (sin_t, cos_t) = (theta.sin(), theta.cos());
    let mut rng = Rng::from_seed(0x010B);
    let mut max_rel = 0.0_f64;
    for _ in 0..100_000 {
        let a = random_unit(&mut rng);
        let t = random_tangent(&mut rng, a);
        let b = a * cos_t + t * sin_t;
        let rel = (great_circle_rad(a, b) - theta).abs() / theta;
        max_rel = max_rel.max(rel);
    }
    println!("separation 1e-9 rad: max relative error {max_rel:e}");
    assert!(max_rel < 1.0e-6, "relative error {max_rel:e}");
}

#[test]
fn great_circle_near_antipodal() {
    // b = −a cos θ + t sin θ lies at angle π − θ from a.
    let theta: f64 = 1.0e-9;
    let (sin_t, cos_t) = (theta.sin(), theta.cos());
    let mut rng = Rng::from_seed(0x010C);
    let mut max_abs = 0.0_f64;
    for _ in 0..100_000 {
        let a = random_unit(&mut rng);
        let t = random_tangent(&mut rng, a);
        let b = -a * cos_t + t * sin_t;
        max_abs = max_abs.max((great_circle_rad(a, b) - (PI - theta)).abs());
    }
    println!("separation π − 1e-9 rad: max absolute error {max_abs:e}");
    assert!(max_abs <= 2.0 * TOL, "absolute error {max_abs:e}");
}

#[test]
fn great_circle_matches_chord_formula() {
    // Independent reference: the chord c = |a − b| gives θ = 2 asin(c / 2),
    // accurate below π/2; above it, θ = π − 2 asin(|a + b| / 2).
    let mut rng = Rng::from_seed(0x010D);
    let mut max_abs = 0.0_f64;
    for _ in 0..100_000 {
        let a = random_unit(&mut rng);
        let b = random_unit(&mut rng);
        let reference = if a.dot(b) >= 0.0 {
            2.0 * ((a - b).length() / 2.0).asin()
        } else {
            PI - 2.0 * ((a + b).length() / 2.0).asin()
        };
        max_abs = max_abs.max((great_circle_rad(a, b) - reference).abs());
    }
    println!("random pairs against the chord formula: max error {max_abs:e}");
    assert!(max_abs <= 2.0 * TOL, "error {max_abs:e}");
}

#[test]
fn great_circle_depends_only_on_direction() {
    let mut rng = Rng::from_seed(0x010F);
    for _ in 0..10_000 {
        let a = random_unit(&mut rng);
        let b = random_unit(&mut rng);
        let d = great_circle_rad(a, b);
        for (sa, sb) in [(1.0e-3, 1.0), (7.0, 0.25), (1.0e3, 1.0e-3)] {
            let ds = great_circle_rad(a * sa, b * sb);
            assert!((ds - d).abs() <= 2.0 * TOL, "{a:?} × {sa}, {b:?} × {sb}");
        }
    }
}

#[test]
fn great_circle_matches_rotation_angle() {
    // A vector perpendicular to the axis turns by |θ| for θ in [−π, π].
    let mut rng = Rng::from_seed(0x010E);
    let mut max_abs = 0.0_f64;
    for _ in 0..100_000 {
        let k = random_unit(&mut rng);
        let v = random_tangent(&mut rng, k);
        let theta = rng.range_f64(-PI, PI);
        let d = great_circle_rad(v, v.rotate_about(k, theta));
        max_abs = max_abs.max((d - theta.abs()).abs());
    }
    println!("rotation angle against distance: max error {max_abs:e}");
    assert!(max_abs <= 2.0 * TOL, "error {max_abs:e}");
}

// ---------------------------------------------------------------------------
// Output fingerprints

/// FNV-1a over the little-endian bytes of the output bits.
fn fingerprint(outputs: impl Iterator<Item = f64>) -> u64 {
    let bytes: Vec<u8> = outputs.flat_map(|y| y.to_bits().to_le_bytes()).collect();
    fnv1a64(&bytes)
}

/// Components of a vector, for fingerprinting.
fn components(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}

fn fingerprints() -> Vec<(&'static str, u64)> {
    let mut angles = round_trip_inputs();
    angles.truncate(2_000);
    angles.extend(
        round_trip_inputs()
            .into_iter()
            .skip(ROUND_TRIP_POINTS - 2_000),
    );
    let mut rng = Rng::from_seed(0x01F0);
    let mut points: Vec<Vec3> = (0..4_000).map(|_| random_unit(&mut rng)).collect();
    points.extend([Vec3::Z, -Vec3::Z, -Vec3::X, Vec3::new(-1.0, -0.0, 0.0)]);
    points.extend(near_pole_vectors());
    let pairs: Vec<(Vec3, Vec3)> = points
        .iter()
        .zip(points.iter().rev())
        .map(|(&a, &b)| (a, b))
        .collect();
    vec![
        (
            "from_lat_lon",
            fingerprint(
                angles
                    .iter()
                    .flat_map(|&(lat, lon)| components(from_lat_lon(lat, lon))),
            ),
        ),
        (
            "lat_lon",
            fingerprint(points.iter().flat_map(|&p| {
                let (lat, lon) = lat_lon(p);
                [lat, lon]
            })),
        ),
        (
            "east",
            fingerprint(points.iter().flat_map(|&p| components(east(p)))),
        ),
        (
            "north",
            fingerprint(points.iter().flat_map(|&p| components(north(p)))),
        ),
        (
            "great_circle_rad",
            fingerprint(pairs.iter().map(|&(a, b)| great_circle_rad(a, b))),
        ),
    ]
}

/// Recorded on x86_64 Linux, in debug and in release, after every invariant
/// and reference test above had passed and the maps of the visual check
/// (ticket M1-02, Result) had been inspected. They are the determinism contract of
/// `coords`: a mismatch means some result changed in its last bits, which
/// changes every world that uses these functions (world-breaking).
const RECORDED_FINGERPRINTS: &[(&str, u64)] = &[
    ("from_lat_lon", 0x884E_F6EE_1EF7_55C3),
    ("lat_lon", 0x9F57_860A_2E20_7D8A),
    ("east", 0x8C4B_55B7_FC5C_E077),
    ("north", 0x508A_B3C9_B099_2304),
    ("great_circle_rad", 0xAFA3_B1F3_2AA3_CB79),
];

#[test]
fn fingerprints_match_recorded() {
    let actual = fingerprints();
    let names: Vec<&str> = actual.iter().map(|a| a.0).collect();
    let recorded_names: Vec<&str> = RECORDED_FINGERPRINTS.iter().map(|r| r.0).collect();
    assert_eq!(names, recorded_names);
    let mismatches: Vec<String> = actual
        .iter()
        .zip(RECORDED_FINGERPRINTS)
        .filter(|(a, r)| a.1 != r.1)
        .map(|(a, r)| format!("coords::{}: {:#018X}, recorded {:#018X}", a.0, a.1, r.1))
        .collect();
    assert!(
        mismatches.is_empty(),
        "output fingerprints changed (world-breaking):\n{}",
        mismatches.join("\n")
    );
}

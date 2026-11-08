//! A small `f64` 3-vector.
//!
//! Positions on the planet are unit vectors in the planet-centred frame:
//! right-handed, +Z along the rotation axis towards the north pole, +X through
//! longitude 0 on the equator (see docs/conventions.md). Latitude, longitude
//! and the local east/north frame are in `kuniumi_sphere::coords`; rotations
//! about an axis are [`Vec3::rotate_about`].

use core::ops::{Add, AddAssign, Div, Mul, Neg, Sub, SubAssign};

use crate::det;

/// Largest `| |axis|² − 1 |` that [`Vec3::rotate_about`] accepts as a unit
/// axis in debug builds (the unit tolerance of the cell grid).
const UNIT_AXIS_TOLERANCE: f64 = 1.0e-12;

/// A 3-vector of `f64`.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Vec3 {
    /// X component.
    pub x: f64,
    /// Y component.
    pub y: f64,
    /// Z component.
    pub z: f64,
}

impl Vec3 {
    /// The zero vector.
    pub const ZERO: Self = Self::new(0.0, 0.0, 0.0);
    /// Unit vector along +X (longitude 0 on the equator).
    pub const X: Self = Self::new(1.0, 0.0, 0.0);
    /// Unit vector along +Y (longitude +90° on the equator).
    pub const Y: Self = Self::new(0.0, 1.0, 0.0);
    /// Unit vector along +Z (the north pole).
    pub const Z: Self = Self::new(0.0, 0.0, 1.0);

    /// Creates a vector.
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// Dot product.
    pub fn dot(self, other: Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// Cross product.
    pub fn cross(self, other: Self) -> Self {
        Self::new(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )
    }

    /// Squared Euclidean length.
    pub fn length_squared(self) -> f64 {
        self.dot(self)
    }

    /// Euclidean length. `sqrt` is correctly rounded by IEEE 754, so this is
    /// deterministic.
    pub fn length(self) -> f64 {
        self.length_squared().sqrt()
    }

    /// Returns the vector scaled to unit length, or `None` for the zero vector.
    pub fn normalized(self) -> Option<Self> {
        let len = self.length();
        if len > 0.0 { Some(self / len) } else { None }
    }

    /// Rotates the vector about the unit vector `axis` by `angle_rad`.
    ///
    /// Right-hand rule: a positive angle turns the vector counter-clockwise
    /// as seen from the tip of `axis` looking towards the origin, so
    /// `Vec3::X.rotate_about(Vec3::Z, π/2)` is `Vec3::Y`.
    ///
    /// Rodrigues' rotation formula (O. Rodrigues, "Des lois géométriques qui
    /// régissent les déplacements d'un système solide dans l'espace",
    /// J. Math. Pures Appl. 5, 380–440, 1840), with `k` the axis and `θ` the
    /// angle, in the form that splits `v` into its parts along and across
    /// the axis:
    ///
    /// `R v = v∥ + (v⊥ cos θ + (k × v) sin θ / |k|)`, with
    /// `v∥ = k (k · v) / |k|²` and `v⊥ = v − v∥`.
    ///
    /// The divisions by `|k|²` and `|k|` make the result depend only on the
    /// direction of `k`. A normalised axis is off unit length by up to about
    /// 6e-16 (`|k|² − 1` within ±5.6 · 2⁻⁵³), and without them that error
    /// changes the length of the result: through the axial part by up to 6.4
    /// ULP of 1 for vectors near the axis at θ ≈ π, and through `k × v` by up
    /// to 1.4 ULP across the axis at θ ≈ ±π/2.
    ///
    /// The axial part must also reach the result unchanged by rounding, or
    /// `cos θ ≈ −1` flips the error left in `v⊥`. The projection is therefore
    /// applied twice: `r = v − k c₁` with `c₁ = (k · v) / |k|²`, then
    /// `c₂ = (k · r) / |k|²` and `v⊥ = r − k c₂`, and `k c₂` is added back
    /// with the axial part. With a single projection the residual of a few
    /// 1e-16 |v| changed the length by up to 4.7 ULP.
    ///
    /// `v = ±k` (and its multiples by powers of two) is returned exactly. Two
    /// other arrangements of the formula round worse:
    /// `v cos θ + (k × v) sin θ + k (k · v)(1 − cos θ)` cancels near θ = π,
    /// and `v + (k × v) sin θ − v⊥ (1 − cos θ)` loses `cos θ` when it is tiny
    /// (it turns `X` about `Z` by π/2 into `(1.1e-16, 1, 0)` instead of
    /// `(cos θ, 1, 0) = (6.1e-17, 1, 0)`).
    ///
    /// The sine and cosine come from one call to [`det::sin_cos`], and `√`
    /// is correctly rounded, so the result is the same on every platform.
    /// The largest relative change of length measured exactly by
    /// `examples/rotate_length.rs` is 2.6 ULP of 1 (5.7e-16), over 10⁷
    /// inputs in each of eight runs and 10⁶ in each of two more: random
    /// inputs as in the tests (three seeds); vectors within 1e-16 to 0.3 rad
    /// of the axis (near ±π and at any angle) or of the plane across it
    /// (near ±π/2 and ±π); vectors 45° from the axis near ±π; vectors across
    /// an axis with `|k|² − 1` of 4 to 12 · 2⁻⁵³ near ±π/2; lengths from
    /// 1e-100 to 1e100; and axes with zero components with angles up to
    /// 1e6. This is a sampled maximum, not a proven bound: a directed search
    /// in review (about 2 · 10¹⁰ evaluations) found 3.9 ULP (8.6e-16) for an
    /// axis within the tolerance below and 3.6 ULP (8.0e-16) for a normalised
    /// axis, which leaves a margin of about 14 % to 1e-15. An angle of ±0
    /// returns `v` unchanged.
    ///
    /// `axis` must be a unit vector: `| |axis|² − 1 |` at most 1e-12, checked
    /// by a debug assertion (a zero or arbitrary axis is a bug in the
    /// caller); within that tolerance its length does not matter.
    pub fn rotate_about(self, axis: Self, angle_rad: f64) -> Self {
        debug_assert!(
            (axis.length_squared() - 1.0).abs() <= UNIT_AXIS_TOLERANCE,
            "rotate_about requires a unit axis, got {axis:?}"
        );
        // Exactly the identity, so that a zero rotation never drifts.
        if angle_rad == 0.0 {
            return self;
        }
        let (sin, cos) = det::sin_cos(angle_rad);
        let k2 = axis.length_squared();
        let parallel = axis * (axis.dot(self) / k2);
        let rest = self - parallel;
        // The axial residual that the rounding of `parallel` leaves in `rest`.
        let correction = axis * (axis.dot(rest) / k2);
        let perpendicular = rest - correction;
        let across = axis.cross(self) * (sin / k2.sqrt());
        parallel + (correction + (perpendicular * cos + across))
    }
}

impl Add for Vec3 {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}

impl AddAssign for Vec3 {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl Sub for Vec3 {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

impl SubAssign for Vec3 {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl Mul<f64> for Vec3 {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        Self::new(self.x * rhs, self.y * rhs, self.z * rhs)
    }
}

impl Div<f64> for Vec3 {
    type Output = Self;
    fn div(self, rhs: f64) -> Self {
        Self::new(self.x / rhs, self.y / rhs, self.z / rhs)
    }
}

impl Neg for Vec3 {
    type Output = Self;
    fn neg(self) -> Self {
        Self::new(-self.x, -self.y, -self.z)
    }
}

// The rotation tests use std's `sin` and `cos` for an independent reference.
#[cfg(test)]
#[allow(clippy::disallowed_methods)]
mod tests {
    use core::f64::consts::{FRAC_PI_2, PI};

    use super::*;
    use crate::rng::{Rng, fnv1a64};

    /// The tolerance of the acceptance criteria.
    const TOL: f64 = 1.0e-15;

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

    /// Error bound for a single rotation of a unit vector, compared with an
    /// independent evaluation: `rotate_about` sums three terms of length up
    /// to 1, each after two to four roundings of at most 2⁻⁵³ |term| (dot
    /// product, projection, subtraction, products with sin and cos), so each
    /// component can be off by several 1e-16; the reference has errors of the
    /// same size. About 9 ULP of 1 leaves a margin over the measured maxima
    /// (printed by the tests, about 1.1e-15).
    const ROTATION_TOL: f64 = 2.0e-15;

    /// `|v|²` as an unevaluated sum `hi + lo`, exact to about 2⁻¹⁰⁰ relative:
    /// each square is split exactly with `mul_add` and the terms are added
    /// with Knuth's TwoSum.
    fn length_squared_dd(v: Vec3) -> (f64, f64) {
        let (mut hi, mut lo) = (0.0_f64, 0.0_f64);
        for c in [v.x, v.y, v.z] {
            let p = c * c;
            let p_err = c.mul_add(c, -p);
            let s = hi + p;
            let b = s - hi;
            let s_err = (hi - (s - b)) + (p - b);
            hi = s;
            lo += s_err + p_err;
        }
        (hi, lo)
    }

    /// `|r| / |v| − 1` to about 1e-30, so that the measurement does not add
    /// rounding errors of the size being measured. The high parts of `|r|²`
    /// and `|v|²` are within a factor 2 of each other, so their difference is
    /// exact (Sterbenz); `√(1 + x) − 1` is evaluated as `x / (√(1 + x) + 1)`.
    fn relative_length_change(r: Vec3, v: Vec3) -> f64 {
        let (r_hi, r_lo) = length_squared_dd(r);
        let (v_hi, v_lo) = length_squared_dd(v);
        let x = ((r_hi - v_hi) + (r_lo - v_lo)) / (v_hi + v_lo);
        x / ((1.0 + x).sqrt() + 1.0)
    }

    /// Rotation by the unit quaternion `(cos θ/2, k sin θ/2)`, evaluated with
    /// std's functions: `v + 2w (u × v) + 2 u × (u × v)` with `w = cos θ/2`,
    /// `u = k sin θ/2`. Independent of Rodrigues' formula in `rotate_about`.
    fn quaternion_rotate(v: Vec3, k: Vec3, angle: f64) -> Vec3 {
        let half = 0.5 * angle;
        let w = half.cos();
        let u = k * half.sin();
        let uv = u.cross(v);
        v + uv * (2.0 * w) + u.cross(uv) * 2.0
    }

    #[test]
    fn rotate_about_quarter_turn() {
        let r = Vec3::X.rotate_about(Vec3::Z, FRAC_PI_2);
        let err = (r - Vec3::Y).length();
        println!("X about Z by π/2: {r:?}, error {err:e}");
        assert!(err <= 1.0e-16, "{r:?}");
    }

    #[test]
    fn rotate_about_right_hand_rule() {
        // A positive angle turns counter-clockwise seen from the tip of the axis.
        let cases = [
            (Vec3::X, Vec3::Z, Vec3::Y),
            (Vec3::Y, Vec3::X, Vec3::Z),
            (Vec3::Z, Vec3::Y, Vec3::X),
            (Vec3::Y, Vec3::Z, -Vec3::X),
            (Vec3::Z, Vec3::X, -Vec3::Y),
            (Vec3::X, Vec3::Y, -Vec3::Z),
        ];
        for (v, axis, expected) in cases {
            let r = v.rotate_about(axis, FRAC_PI_2);
            assert!(
                (r - expected).length() <= 1.0e-16,
                "{v:?} about {axis:?}: {r:?}"
            );
            let back = v.rotate_about(-axis, FRAC_PI_2);
            assert!(
                (back + expected).length() <= 1.0e-16,
                "{v:?} about -{axis:?}"
            );
            let neg = v.rotate_about(axis, -FRAC_PI_2);
            assert!((neg + expected).length() <= 1.0e-16, "{v:?} by -π/2");
        }
        let half_turn = Vec3::X.rotate_about(Vec3::Z, PI);
        assert!((half_turn + Vec3::X).length() <= TOL, "{half_turn:?}");
    }

    /// Largest exact relative length change of `v.rotate_about(k, angle)`
    /// over `n` inputs drawn by `draw`, with the worst input.
    fn max_length_change(
        n: usize,
        rng: &mut Rng,
        mut draw: impl FnMut(&mut Rng, usize) -> (Vec3, Vec3, f64),
    ) -> (f64, (Vec3, Vec3, f64)) {
        let mut max = (0.0_f64, (Vec3::ZERO, Vec3::ZERO, 0.0));
        for i in 0..n {
            let (v, k, angle) = draw(rng, i);
            let change = relative_length_change(v.rotate_about(k, angle), v).abs();
            if change > max.0 {
                max = (change, (v, k, angle));
            }
        }
        max
    }

    /// A log-uniform number in `[lo, hi)`.
    fn log_uniform(rng: &mut Rng, lo: f64, hi: f64) -> f64 {
        rng.range_f64(lo.ln(), hi.ln()).exp()
    }

    /// The unit vector at angle `a` from the unit vector `k`: `k` turned
    /// towards a random tangent `t = k × u` (with `u` at least 30° from `±k`).
    fn near_axis(rng: &mut Rng, k: Vec3, a: f64) -> Vec3 {
        let t = loop {
            let t = k.cross(random_unit(rng));
            if t.length_squared() > 0.25 {
                break t.normalized().expect("non-zero");
            }
        };
        (k * a.cos() + t * a.sin()).normalized().expect("non-zero")
    }

    #[test]
    fn rotate_about_preserves_length() {
        // Criterion: the length changes by at most 1e-15 relative, measured
        // exactly. 10⁶ random unit axes (normalised, so off unit length by
        // up to about 1e-16), random vectors (half unit, half of length 1e-3
        // to 1e3) and angles in [−10π, 10π].
        let mut rng = Rng::from_seed(0x0201);
        let (max, worst) = max_length_change(1_000_000, &mut rng, |rng, i| {
            let k = random_unit(rng);
            let scale = if i % 2 == 0 {
                1.0
            } else {
                rng.range_f64(1.0e-3, 1.0e3)
            };
            (
                random_unit(rng) * scale,
                k,
                rng.range_f64(-10.0 * PI, 10.0 * PI),
            )
        });
        println!("rotation, 10⁶ random inputs: max relative length change {max:e} at {worst:?}");
        assert!(max <= TOL, "length changed by {max:e} at {worst:?}");
    }

    #[test]
    fn rotate_about_preserves_length_near_axis_and_half_turn() {
        // Where the axial component is most sensitive: vectors 1e-16 to 0.3
        // rad from the axis (either direction) and angles 1e-12 to 0.3 rad
        // from ±π, both log-uniform so that every scale is sampled; 10⁶
        // inputs, measured exactly.
        let mut rng = Rng::from_seed(0x0205);
        let (max, worst) = max_length_change(1_000_000, &mut rng, |rng, i| {
            let k = random_unit(rng);
            let distance = log_uniform(rng, 1.0e-16, 0.3);
            let v = near_axis(rng, k, distance);
            let v = if i % 4 < 2 { v } else { -v };
            let near_pi = PI - log_uniform(rng, 1.0e-12, 0.3);
            (v, k, if i % 2 == 0 { near_pi } else { -near_pi })
        });
        println!(
            "rotation near the axis, near ±π: max relative length change {max:e} at {worst:?}"
        );
        assert!(max <= TOL, "length changed by {max:e} at {worst:?}");
    }

    /// `|k|² − 1`, exactly (to about 1e-30).
    fn exact_length_squared_minus_one(k: Vec3) -> f64 {
        let (hi, lo) = length_squared_dd(k);
        (hi - 1.0) + lo
    }

    #[test]
    fn rotate_about_preserves_length_across_axis_with_skewed_axis() {
        // Where the length of the axis would matter most: vectors across the
        // axis with |v| just above 1 (the largest relative rounding) and
        // angles near ±π/2, where the `k × v` term carries the whole part
        // across the axis. Half the axes are normalised and skewed to
        // |k|² − 1 of 4 to 12 · 2⁻⁵³ (normalised vectors reach about ±5.6),
        // half are off unit length by up to 8e-13, near the tolerance of the
        // debug assertion; the result must not depend on it. 10⁶ inputs,
        // measured exactly. Without the division of `k × v` by |k| the first
        // half reaches 3.8 ULP (4.1 ULP over the 10⁷ inputs of
        // `examples/rotate_length.rs`) and the second half 4.0e-13.
        let mut rng = Rng::from_seed(0x0207);
        // A normalised axis, with its largest component stepped by up to 3
        // ULPs in the direction that makes |k|² − 1 largest (rejection
        // sampling among plain normalised axes would be too slow).
        let skewed_axis = |rng: &mut Rng| loop {
            let k = random_unit(rng);
            let mut best = (exact_length_squared_minus_one(k).abs(), k);
            let largest = if k.x.abs() >= k.y.abs().max(k.z.abs()) {
                0
            } else if k.y.abs() >= k.z.abs() {
                1
            } else {
                2
            };
            for step in [-3_i64, -2, -1, 1, 2, 3] {
                let mut c = [k.x, k.y, k.z];
                c[largest] = f64::from_bits(c[largest].to_bits().wrapping_add_signed(step));
                let candidate = Vec3::new(c[0], c[1], c[2]);
                let skew = exact_length_squared_minus_one(candidate).abs();
                if skew > best.0 {
                    best = (skew, candidate);
                }
            }
            if best.0 >= 4.0 * f64::EPSILON / 2.0 {
                return best.1;
            }
        };
        let (max, worst) = max_length_change(1_000_000, &mut rng, |rng, i| {
            let k = skewed_axis(rng);
            let v = near_axis(rng, k, FRAC_PI_2) * rng.range_f64(1.0, 1.05);
            let k = if i % 8 < 4 {
                k
            } else {
                let off = log_uniform(rng, 1.0e-15, 4.0e-13);
                k * (1.0 + if i % 2 == 0 { off } else { -off })
            };
            let offset = log_uniform(rng, 1.0e-6, 0.2);
            let angle = FRAC_PI_2 + if i % 2 == 0 { offset } else { -offset };
            (v, k, if i % 4 < 2 { angle } else { -angle })
        });
        println!(
            "rotation across a skewed axis, near ±π/2: max relative length change {max:e} \
             at {worst:?}"
        );
        assert!(max <= TOL, "length changed by {max:e} at {worst:?}");
    }

    #[test]
    fn rotate_about_review_cases() {
        // Inputs found in review (M1-02). Round 1: without the division by
        // |k|² the length changed by 5.2 ULP (and 6.4 ULP for k about k by π
        // below). Round 2: with a single projection, by 4.72, 4.72 and 4.51
        // ULP. Round 3: the largest changes found by search for the formula
        // without the division of k × v by |k|, 3.85 and 3.81 ULP. Round 4:
        // the largest found by directed search for the current formula, 3.88
        // ULP (axis off unit length by 5.1e-13) and 3.60 ULP (a normalised
        // axis). All confirmed with exact rational arithmetic.
        let v3 = |x: u64, y: u64, z: u64| {
            Vec3::new(f64::from_bits(x), f64::from_bits(y), f64::from_bits(z))
        };
        let cases = [
            (
                v3(
                    0x4050_950B_2A7B_9813,
                    0xC082_FF3E_B48E_55CC,
                    0x4070_34B7_C919_26BA,
                ),
                v3(
                    0xBFB0_EC02_4A04_A8FA,
                    0x3FEC_D09C_C690_A8F8,
                    0xBFDB_8313_7645_4E91,
                ),
                0xC02F_9FE4_067D_8ADB,
            ),
            (
                v3(
                    0x4050_628E_D086_5514,
                    0xC028_239B_98AE_B7B4,
                    0xC012_D5F6_FEBA_96AB,
                ),
                v3(
                    0x3FEF_5980_66D1_5197,
                    0xBFC1_D984_B7C3_66DB,
                    0xBFC2_74AB_D700_DB3C,
                ),
                0xC021_96E0_9405_D1EE,
            ),
            (
                v3(
                    0xBFE6_F782_0E14_811F,
                    0x3FE5_A8A9_AD8F_7643,
                    0x3FC4_F208_F8CE_7AAD,
                ),
                v3(
                    0xBFE6_F782_0DE4_DA64,
                    0x3FE5_A8A9_ADD9_0A3D,
                    0x3FC4_F208_F751_223F,
                ),
                0xC007_FC6C_3956_1264,
            ),
            (
                v3(
                    0xBFE7_4419_A504_39BD,
                    0x3FE5_265E_2365_D01B,
                    0xBFC7_CA75_E3D8_E5D3,
                ),
                v3(
                    0xBFE7_4419_A521_7705,
                    0x3FE5_265E_2305_52A9,
                    0xBFC7_CA75_E76B_D96F,
                ),
                0xC007_2E61_B3D5_FB7F,
            ),
            (
                v3(
                    0xBF2B_D2BE_FDAE_953E,
                    0xBFE7_06E6_36E0_A16F,
                    0x3FE7_0C56_5D0C_1A27,
                ),
                v3(
                    0xBFC6_14C3_5D81_118C,
                    0xBFE6_5263_0390_EA57,
                    0xBFE6_4126_E1FA_17CE,
                ),
                0x3FF9_21E1_C951_192F,
            ),
            (
                v3(
                    0xBFF0_9F2A_7C5D_1CD4,
                    0x3F12_C764_D1C0_7992,
                    0x3F4C_3BBA_18A2_3535,
                ),
                v3(
                    0xBFC4_468C_D455_B241,
                    0x3FEE_DA56_EC2B_C500,
                    0x3FCB_3FA7_2471_7B78,
                ),
                0x3FF9_4581_49A2_66B4,
            ),
            (
                v3(
                    0xBFE7_F697_FBB3_35E2,
                    0x3FE7_5943_14A2_4FC7,
                    0x3F97_6867_0C5E_3B3F,
                ),
                v3(
                    0xBFE5_FA23_A491_2373,
                    0xBFE6_B79F_34C8_EA3A,
                    0x3FC3_F910_1D23_13BB,
                ),
                0xBFF7_ED66_C3B9_420A,
            ),
            (
                v3(
                    0x3FC2_63DE_3331_89BC,
                    0x3FAA_2107_6B54_B695,
                    0xBFF0_59EE_8230_373D,
                ),
                v3(
                    0xBFAF_A023_A151_DABD,
                    0xBFEF_E254_70D5_5370,
                    0xBFAD_EBD1_67C5_1EDF,
                ),
                0x3FF9_23BA_E157_A194,
            ),
        ];
        for (v, k, angle) in cases {
            let angle = f64::from_bits(angle);
            let change = relative_length_change(v.rotate_about(k, angle), v).abs();
            println!("review case: relative length change {change:e}");
            assert!(
                change <= TOL,
                "{v:?} about {k:?}: length changed by {change:e}"
            );
        }
        let k = v3(
            0x3FE8_7636_CC78_1121,
            0x3FE4_2A90_73F5_1372,
            0x3FC1_68D9_7E53_F033,
        );
        assert_eq!(k.rotate_about(k, PI), k);
    }

    #[test]
    fn rotate_about_keeps_vectors_on_the_axis() {
        // v = ±k, and its multiples by powers of two, is exactly its own
        // projection (k · v / |k|² is then exactly ±1, ±4, ...), so it comes
        // back exactly whatever the angle; 10⁵ normalised axes.
        let mut rng = Rng::from_seed(0x0206);
        for i in 0..100_000 {
            let k = random_unit(&mut rng);
            let v = k * [1.0, -1.0, 4.0, -0.25][i % 4];
            for angle in [PI, -PI, FRAC_PI_2, rng.range_f64(-10.0 * PI, 10.0 * PI)] {
                assert_eq!(v.rotate_about(k, angle), v, "{v:?} about {k:?} by {angle}");
            }
        }
    }

    #[test]
    fn rotate_about_matches_quaternion_rotation() {
        let mut rng = Rng::from_seed(0x0202);
        let mut max_err = 0.0_f64;
        for _ in 0..100_000 {
            let k = random_unit(&mut rng);
            let v = random_unit(&mut rng);
            let angle = rng.range_f64(-10.0 * PI, 10.0 * PI);
            let err = (v.rotate_about(k, angle) - quaternion_rotate(v, k, angle)).length();
            max_err = max_err.max(err);
        }
        println!("rotation against quaternion rotation: max error {max_err:e}");
        assert!(max_err <= ROTATION_TOL, "error {max_err:e}");
    }

    #[test]
    fn rotate_about_fixes_axis_and_inverts() {
        let mut rng = Rng::from_seed(0x0203);
        let (mut max_axis, mut max_inverse) = (0.0_f64, 0.0_f64);
        for _ in 0..100_000 {
            let k = random_unit(&mut rng);
            let v = random_unit(&mut rng);
            let angle = rng.range_f64(-10.0 * PI, 10.0 * PI);
            max_axis = max_axis.max((k.rotate_about(k, angle) - k).length());
            let back = v.rotate_about(k, angle).rotate_about(k, -angle);
            max_inverse = max_inverse.max((back - v).length());
        }
        println!("rotation: axis moves by {max_axis:e}; θ then −θ errs by {max_inverse:e}");
        assert!(max_axis <= ROTATION_TOL, "axis moved by {max_axis:e}");
        // Two rotations, so twice the bound of one.
        assert!(
            max_inverse <= 2.0 * ROTATION_TOL,
            "inverse error {max_inverse:e}"
        );
    }

    #[test]
    fn rotate_about_zero_angle_is_identity() {
        let mut rng = Rng::from_seed(0x0204);
        for _ in 0..10_000 {
            let k = random_unit(&mut rng);
            let v = random_unit(&mut rng) * 3.0;
            for zero in [0.0, -0.0] {
                let r = v.rotate_about(k, zero);
                let bits = |a: Vec3| [a.x.to_bits(), a.y.to_bits(), a.z.to_bits()];
                assert_eq!(bits(r), bits(v), "{v:?} about {k:?}");
            }
        }
    }

    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "unit axis")]
    fn rotate_about_rejects_non_unit_axis() {
        let _ = Vec3::X.rotate_about(Vec3::new(0.0, 0.0, 2.0), 1.0);
    }

    /// Recorded on x86_64 Linux, in debug and in release, after the rotation
    /// tests above had passed and the rotated graticule of the visual check
    /// (ticket M1-02, Result) had been inspected: FNV-1a of the output bits of `rotate_about`
    /// over 10⁴ fixed inputs. A mismatch means the rotation changed in its
    /// last bits, which changes every world that uses it (world-breaking).
    const RECORDED_ROTATE_ABOUT_FINGERPRINT: u64 = 0xED69_192B_E30E_1591;

    #[test]
    fn rotate_about_fingerprint_matches_recorded() {
        let mut rng = Rng::from_seed(0x02F0);
        let mut bytes = Vec::new();
        for _ in 0..10_000 {
            let k = random_unit(&mut rng);
            let v = random_unit(&mut rng);
            let angle = rng.range_f64(-10.0 * PI, 10.0 * PI);
            let r = v.rotate_about(k, angle);
            for c in [r.x, r.y, r.z] {
                bytes.extend_from_slice(&c.to_bits().to_le_bytes());
            }
        }
        let actual = fnv1a64(&bytes);
        assert_eq!(
            actual, RECORDED_ROTATE_ABOUT_FINGERPRINT,
            "rotate_about fingerprint changed (world-breaking): {actual:#018X}"
        );
    }

    #[test]
    fn right_handed_axes() {
        assert_eq!(Vec3::X.cross(Vec3::Y), Vec3::Z);
        assert_eq!(Vec3::Y.cross(Vec3::Z), Vec3::X);
        assert_eq!(Vec3::Z.cross(Vec3::X), Vec3::Y);
    }

    #[test]
    fn normalize() {
        let v = Vec3::new(3.0, 0.0, 4.0).normalized().unwrap();
        assert_eq!(v, Vec3::new(0.6, 0.0, 0.8));
        assert!(Vec3::ZERO.normalized().is_none());
    }
}

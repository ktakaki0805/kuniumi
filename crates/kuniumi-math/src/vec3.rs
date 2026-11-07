//! A small `f64` 3-vector.
//!
//! Positions on the planet are unit vectors in the planet-centred frame:
//! right-handed, +Z along the rotation axis towards the north pole, +X through
//! longitude 0 on the equator (see docs/conventions.md).

use core::ops::{Add, AddAssign, Div, Mul, Neg, Sub, SubAssign};

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

#[cfg(test)]
mod tests {
    use super::*;

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

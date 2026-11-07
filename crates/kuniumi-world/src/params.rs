use core::fmt;

use kuniumi_math::units::deg_to_rad;

/// How the planet spins relative to its star.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RotationRegime {
    /// Rotates freely with a day/night cycle (Earth-like). Supported in v1.
    Free,
    /// Always shows the same face to its star. Planned for v2: the substellar
    /// point is fixed at longitude 0, latitude 0.
    TidallyLocked,
}

/// How the planet's lithosphere moves.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TectonicRegime {
    /// Mobile plates (Earth-like). Supported in v1.
    PlateTectonics,
    /// A single immobile lid shaped by volcanism and impacts. Planned for v2.
    StagnantLid,
}

/// Everything that describes a planet before generation. All values are SI.
#[derive(Clone, Debug, PartialEq)]
pub struct PlanetParams {
    /// Mean planetary radius (m).
    pub radius_m: f64,
    /// Surface gravitational acceleration (m/s²).
    pub surface_gravity_m_s2: f64,
    /// Sidereal rotation period (s). Ignored when tidally locked.
    pub rotation_period_s: f64,
    /// Axial tilt relative to the orbital plane (rad), in `[0, π]`.
    pub obliquity_rad: f64,
    /// Orbital period around the star (s).
    pub orbital_period_s: f64,
    /// Orbital eccentricity, in `[0, 1)`.
    pub orbital_eccentricity: f64,
    /// Mean stellar flux at the planet's orbit (W/m²).
    pub stellar_flux_w_m2: f64,
    /// Total volume of surface water (m³). Sea level is solved from this.
    pub ocean_volume_m3: f64,
    /// Mean surface air pressure (Pa).
    pub surface_pressure_pa: f64,
    /// Fraction of the surface underlain by continental crust, in `[0, 1]`.
    pub continental_crust_fraction: f64,
    /// Rotation regime.
    pub rotation: RotationRegime,
    /// Tectonic regime.
    pub tectonics: TectonicRegime,
}

impl PlanetParams {
    /// Approximate present-day Earth.
    pub fn earth() -> Self {
        Self {
            radius_m: 6_371_000.0,
            surface_gravity_m_s2: 9.806_65,
            rotation_period_s: 86_164.090_5,
            obliquity_rad: deg_to_rad(23.44),
            orbital_period_s: 31_558_149.8,
            orbital_eccentricity: 0.0167,
            stellar_flux_w_m2: 1361.0,
            ocean_volume_m3: 1.335e18,
            surface_pressure_pa: 101_325.0,
            continental_crust_fraction: 0.40,
            rotation: RotationRegime::Free,
            tectonics: TectonicRegime::PlateTectonics,
        }
    }

    /// Checks ranges, and rejects regimes that this version cannot generate.
    pub fn validate(&self) -> Result<(), ParamError> {
        fn positive(name: &'static str, v: f64) -> Result<(), ParamError> {
            if v.is_finite() && v > 0.0 {
                Ok(())
            } else {
                Err(ParamError::OutOfRange { name, value: v })
            }
        }
        fn within(name: &'static str, v: f64, lo: f64, hi: f64) -> Result<(), ParamError> {
            if v.is_finite() && v >= lo && v <= hi {
                Ok(())
            } else {
                Err(ParamError::OutOfRange { name, value: v })
            }
        }
        positive("radius_m", self.radius_m)?;
        positive("surface_gravity_m_s2", self.surface_gravity_m_s2)?;
        positive("rotation_period_s", self.rotation_period_s)?;
        within(
            "obliquity_rad",
            self.obliquity_rad,
            0.0,
            core::f64::consts::PI,
        )?;
        positive("orbital_period_s", self.orbital_period_s)?;
        within("orbital_eccentricity", self.orbital_eccentricity, 0.0, 0.99)?;
        positive("stellar_flux_w_m2", self.stellar_flux_w_m2)?;
        within("ocean_volume_m3", self.ocean_volume_m3, 0.0, f64::MAX)?;
        positive("surface_pressure_pa", self.surface_pressure_pa)?;
        within(
            "continental_crust_fraction",
            self.continental_crust_fraction,
            0.0,
            1.0,
        )?;
        if self.rotation == RotationRegime::TidallyLocked {
            return Err(ParamError::Unsupported {
                feature: "tidally locked rotation",
                planned: "v2",
            });
        }
        if self.tectonics == TectonicRegime::StagnantLid {
            return Err(ParamError::Unsupported {
                feature: "stagnant-lid tectonics",
                planned: "v2",
            });
        }
        Ok(())
    }
}

/// Why a set of [`PlanetParams`] cannot be generated.
#[derive(Clone, Debug, PartialEq)]
pub enum ParamError {
    /// A value is not finite or outside its valid range.
    OutOfRange {
        /// Field name.
        name: &'static str,
        /// The rejected value.
        value: f64,
    },
    /// The parameters ask for a feature this version does not implement yet.
    Unsupported {
        /// The missing feature.
        feature: &'static str,
        /// The release that is planned to add it.
        planned: &'static str,
    },
}

impl fmt::Display for ParamError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::OutOfRange { name, value } => write!(f, "{name} is out of range: {value}"),
            Self::Unsupported { feature, planned } => {
                write!(f, "{feature} is not supported yet (planned for {planned})")
            }
        }
    }
}

impl std::error::Error for ParamError {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn earth_is_valid() {
        assert_eq!(PlanetParams::earth().validate(), Ok(()));
    }

    #[test]
    fn rejects_negative_radius() {
        let p = PlanetParams {
            radius_m: -1.0,
            ..PlanetParams::earth()
        };
        assert!(matches!(
            p.validate(),
            Err(ParamError::OutOfRange {
                name: "radius_m",
                ..
            })
        ));
    }

    #[test]
    fn v2_regimes_are_reported_as_unsupported() {
        let p = PlanetParams {
            rotation: RotationRegime::TidallyLocked,
            ..PlanetParams::earth()
        };
        assert!(matches!(
            p.validate(),
            Err(ParamError::Unsupported { planned: "v2", .. })
        ));
        let p = PlanetParams {
            tectonics: TectonicRegime::StagnantLid,
            ..PlanetParams::earth()
        };
        assert!(matches!(
            p.validate(),
            Err(ParamError::Unsupported { planned: "v2", .. })
        ));
    }
}

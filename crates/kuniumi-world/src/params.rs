use core::f64::consts::PI;
use core::fmt;

use kuniumi_math::constants::GRAVITATIONAL_CONSTANT_M3_KG_S2;
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
///
/// v1 generates Earth-like planets (D-12): [`validate`](Self::validate)
/// bounds every quantity (the `MIN_` and `MAX_` constants, ticket M2-08)
/// and the bulk density that the radius and the gravity imply.
#[derive(Clone, Debug, PartialEq)]
pub struct PlanetParams {
    /// Mean planetary radius (m), in
    /// [`MIN_RADIUS_M`](Self::MIN_RADIUS_M)`..=`[`MAX_RADIUS_M`](Self::MAX_RADIUS_M)
    /// (1,500 to 13,000 km), with a
    /// [`bulk_density_kg_m3`](Self::bulk_density_kg_m3) in range.
    pub radius_m: f64,
    /// Surface gravitational acceleration (m/s²), in
    /// [`MIN_SURFACE_GRAVITY_M_S2`](Self::MIN_SURFACE_GRAVITY_M_S2)`..=`
    /// [`MAX_SURFACE_GRAVITY_M_S2`](Self::MAX_SURFACE_GRAVITY_M_S2) (1 to
    /// 40 m/s²), with a [`bulk_density_kg_m3`](Self::bulk_density_kg_m3) in
    /// range.
    pub surface_gravity_m_s2: f64,
    /// Sidereal rotation period (s), in
    /// [`MIN_ROTATION_PERIOD_S`](Self::MIN_ROTATION_PERIOD_S)`..=`
    /// [`MAX_ROTATION_PERIOD_S`](Self::MAX_ROTATION_PERIOD_S) (8 to 240 h).
    /// Ignored when tidally locked.
    pub rotation_period_s: f64,
    /// Axial tilt relative to the orbital plane (rad), in `[0, π]`.
    pub obliquity_rad: f64,
    /// Orbital period around the star (s), in
    /// [`MIN_ORBITAL_PERIOD_S`](Self::MIN_ORBITAL_PERIOD_S)`..=`
    /// [`MAX_ORBITAL_PERIOD_S`](Self::MAX_ORBITAL_PERIOD_S) (90 to 3,650
    /// days).
    pub orbital_period_s: f64,
    /// Orbital eccentricity, in `[0, 1)`.
    pub orbital_eccentricity: f64,
    /// Mean stellar flux at the planet's orbit (W/m²), in
    /// [`MIN_STELLAR_FLUX_W_M2`](Self::MIN_STELLAR_FLUX_W_M2)`..=`
    /// [`MAX_STELLAR_FLUX_W_M2`](Self::MAX_STELLAR_FLUX_W_M2) (400 to
    /// 2,500 W/m²).
    pub stellar_flux_w_m2: f64,
    /// Total volume of surface water (m³), in
    /// `0..=`[`MAX_OCEAN_VOLUME_M3`](Self::MAX_OCEAN_VOLUME_M3) (0 to
    /// 10¹¹ km³). Sea level is solved from this.
    pub ocean_volume_m3: f64,
    /// Mean surface air pressure (Pa), in
    /// [`MIN_SURFACE_PRESSURE_PA`](Self::MIN_SURFACE_PRESSURE_PA)`..=`
    /// [`MAX_SURFACE_PRESSURE_PA`](Self::MAX_SURFACE_PRESSURE_PA) (0.1 to
    /// 10 bar).
    pub surface_pressure_pa: f64,
    /// Fraction of the surface underlain by continental crust, in `[0, 1]`.
    pub continental_crust_fraction: f64,
    /// Rotation regime.
    pub rotation: RotationRegime,
    /// Tectonic regime.
    pub tectonics: TectonicRegime,
}

impl PlanetParams {
    /// The smallest radius: 1.5 × 10⁶ m (1,500 km, 0.24 Earth radii), a
    /// little below the Moon's 1,737 km, so that the Moon, Mercury and Mars
    /// lie inside (ticket M2-08).
    pub const MIN_RADIUS_M: f64 = 1.5e6;

    /// The largest radius: 1.3 × 10⁷ m (13,000 km, 2.04 Earth radii). Above
    /// 1.5 to 2 Earth radii planets mostly have thick volatile envelopes
    /// rather than a rocky surface (the radius gap of Fulton et al. 2017,
    /// AJ 154, 109; Rogers 2015, ApJ 801, 41).
    pub const MAX_RADIUS_M: f64 = 1.3e7;

    /// The lowest surface gravity: 1 m/s² (0.10 g), below the Moon's
    /// 1.62 m/s². With the density bounds no valid planet has less than
    /// 1.26 m/s² (1,500 km at 3,000 kg/m³), so only gravities far below
    /// are reported against this bound; the others fail the density
    /// check.
    pub const MIN_SURFACE_GRAVITY_M_S2: f64 = 1.0;

    /// The highest surface gravity: 40 m/s² (4.1 g). An Earth-composition
    /// planet of 13,000 km has about 33 m/s² (mass–radius relation of Zeng
    /// et al. 2016, ApJ 819, 127, `R = 1.00 M^(1/3.7)` in Earth units).
    pub const MAX_SURFACE_GRAVITY_M_S2: f64 = 40.0;

    /// The shortest rotation period: 28,800 s (8 h), three times Earth's
    /// rotation rate. At Earth's density the centrifugal acceleration at
    /// the equator is then 3.1 % of gravity (Earth: 0.35 %), so the
    /// spherical grid stays a fair approximation. Every planet within the
    /// bounds rotates slower than break-up (see [`validate`](Self::validate)).
    pub const MIN_ROTATION_PERIOD_S: f64 = 28_800.0;

    /// The longest rotation period: 864,000 s (240 h, 10 days), a tenth of
    /// Earth's rotation rate. Slower rotators approach the day–night
    /// circulation of tidally locked planets (v2); with the shortest year a
    /// year still holds at least 8 solar days.
    pub const MAX_ROTATION_PERIOD_S: f64 = 864_000.0;

    /// The shortest orbital period: 7,776,000 s (90 days). At Earth's flux
    /// a 90-day year belongs to a star of about 0.57 solar masses
    /// (`L ∝ M⁴`, so `P = M^2.5` years); planets in the habitable zones of
    /// smaller stars are expected to be tidally locked (Kasting et al.
    /// 1993, Icarus 101, 108), the v2 regime.
    pub const MIN_ORBITAL_PERIOD_S: f64 = 7_776_000.0;

    /// The longest orbital period: 315,360,000 s (3,650 days). At the
    /// lowest flux the habitable zone of an F0 star (6 solar luminosities,
    /// 1.6 solar masses) lies at 4.5 AU, a year of about 8 years.
    pub const MAX_ORBITAL_PERIOD_S: f64 = 315_360_000.0;

    /// The lowest stellar flux: 400 W/m² (0.29 of Earth's), outside the
    /// optimistic habitable zone of a Sun-like star, which starts at
    /// 0.32 (early Mars; Kopparapu et al. 2013, ApJ 765, 131). Mars receives
    /// 586 W/m².
    pub const MIN_STELLAR_FLUX_W_M2: f64 = 400.0;

    /// The highest stellar flux: 2,500 W/m² (1.84 of Earth's), outside the
    /// optimistic habitable zone, which ends at 1.78 (recent Venus;
    /// Kopparapu et al. 2013). Venus, at 2,601 W/m², lies outside.
    pub const MAX_STELLAR_FLUX_W_M2: f64 = 2_500.0;

    /// The largest ocean volume: 10²⁰ m³ (10¹¹ km³, 75 Earth oceans), which
    /// drowns every planet within the bounds: a layer 47 km deep on a
    /// 13,000 km planet, whose relief from the deepest trench to the
    /// highest mountain is about 20 km. On the smallest planet (1,500 km)
    /// it is seven times the planet's volume, and the thin-layer sea level
    /// of the Crust stage (3,533 km above the datum, where a spherical
    /// shell would be 1,509 km deep) is finite but far from physical
    /// (ticket M2-17). The smallest is 0, a dry planet.
    pub const MAX_OCEAN_VOLUME_M3: f64 = 1e20;

    /// The lowest surface pressure: 10⁴ Pa (0.1 bar). Below it water boils
    /// below 46 °C, so an ocean cannot be temperate (the triple point is
    /// 612 Pa).
    pub const MIN_SURFACE_PRESSURE_PA: f64 = 1e4;

    /// The highest surface pressure: 10⁶ Pa (10 bar), far from the Earth
    /// calibration of the climate model beyond it.
    pub const MAX_SURFACE_PRESSURE_PA: f64 = 1e6;

    /// The lowest bulk density ([`bulk_density_kg_m3`](Self::bulk_density_kg_m3)):
    /// 3,000 kg/m³. Rocky bodies of little self-compression are about as
    /// dense as silicate rock: the Moon 3,344 kg/m³, Io 3,528, Mars 3,934.
    pub const MIN_BULK_DENSITY_KG_M3: f64 = 3_000.0;

    /// The highest bulk density: 16,000 kg/m³, a little above the densest
    /// rocky planet within the bounds, pure iron at 13,000 km (about
    /// 15,900 kg/m³; Zeng et al. 2016, `R = 0.86 M^(1/3.7)` for a core
    /// mass fraction of 1). A constant bound rejects no possible rocky
    /// planet within the bounds, but admits small planets denser than iron
    /// could be at their size (uncompressed iron is 7,870 kg/m³).
    pub const MAX_BULK_DENSITY_KG_M3: f64 = 16_000.0;

    /// The mean density (kg/m³) of a planet with this radius and surface
    /// gravity: `ρ = 3 g / (4π G R)`, from `g = G M / R²` and
    /// `M = (4/3) π R³ ρ`, with `G` =
    /// [`GRAVITATIONAL_CONSTANT_M3_KG_S2`]. Earth's parameters give
    /// 5,506 kg/m³ (5,513 from Earth's mass: standard gravity lies below
    /// `G M / R²` because it includes the centrifugal acceleration).
    pub fn bulk_density_kg_m3(&self) -> f64 {
        3.0 * self.surface_gravity_m_s2
            / (4.0 * PI * GRAVITATIONAL_CONSTANT_M3_KG_S2 * self.radius_m)
    }

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

    /// Checks the parameters for a planet that this version can generate,
    /// and reports the first problem:
    ///
    /// 1. every float field within its bounds (the constants above, and
    ///    the ranges given with the fields), inclusive, in field order:
    ///    [`ParamError::OutOfRange`] with the field's name; NaN and ±∞ are
    ///    out of range;
    /// 2. the [`bulk_density_kg_m3`](Self::bulk_density_kg_m3) of the radius
    ///    and the gravity within
    ///    [`MIN_BULK_DENSITY_KG_M3`](Self::MIN_BULK_DENSITY_KG_M3)`..=`
    ///    [`MAX_BULK_DENSITY_KG_M3`](Self::MAX_BULK_DENSITY_KG_M3):
    ///    [`ParamError::Inconsistent`];
    /// 3. the regimes: [`ParamError::Unsupported`] for those of v2.
    ///
    /// The bounds also ensure that every planet rotates slower than
    /// break-up, where the centrifugal acceleration at the equator,
    /// `ω² R` with `ω = 2π / P`, would reach the surface gravity `g`: the
    /// longest break-up period, `2π √(R / g)` at the largest radius and the
    /// lowest gravity, is 6.29 h, shorter than the 8 h of
    /// [`MIN_ROTATION_PERIOD_S`](Self::MIN_ROTATION_PERIOD_S). A
    /// compile-time assertion over the constants keeps it so, so `validate`
    /// does not check it (ticket M2-08).
    pub fn validate(&self) -> Result<(), ParamError> {
        fn within(name: &'static str, v: f64, lo: f64, hi: f64) -> Result<(), ParamError> {
            if v.is_finite() && v >= lo && v <= hi {
                Ok(())
            } else {
                Err(ParamError::OutOfRange { name, value: v })
            }
        }
        within(
            "radius_m",
            self.radius_m,
            Self::MIN_RADIUS_M,
            Self::MAX_RADIUS_M,
        )?;
        within(
            "surface_gravity_m_s2",
            self.surface_gravity_m_s2,
            Self::MIN_SURFACE_GRAVITY_M_S2,
            Self::MAX_SURFACE_GRAVITY_M_S2,
        )?;
        within(
            "rotation_period_s",
            self.rotation_period_s,
            Self::MIN_ROTATION_PERIOD_S,
            Self::MAX_ROTATION_PERIOD_S,
        )?;
        within(
            "obliquity_rad",
            self.obliquity_rad,
            0.0,
            core::f64::consts::PI,
        )?;
        within(
            "orbital_period_s",
            self.orbital_period_s,
            Self::MIN_ORBITAL_PERIOD_S,
            Self::MAX_ORBITAL_PERIOD_S,
        )?;
        within("orbital_eccentricity", self.orbital_eccentricity, 0.0, 0.99)?;
        within(
            "stellar_flux_w_m2",
            self.stellar_flux_w_m2,
            Self::MIN_STELLAR_FLUX_W_M2,
            Self::MAX_STELLAR_FLUX_W_M2,
        )?;
        within(
            "ocean_volume_m3",
            self.ocean_volume_m3,
            0.0,
            Self::MAX_OCEAN_VOLUME_M3,
        )?;
        within(
            "surface_pressure_pa",
            self.surface_pressure_pa,
            Self::MIN_SURFACE_PRESSURE_PA,
            Self::MAX_SURFACE_PRESSURE_PA,
        )?;
        within(
            "continental_crust_fraction",
            self.continental_crust_fraction,
            0.0,
            1.0,
        )?;
        let density = self.bulk_density_kg_m3();
        if !(Self::MIN_BULK_DENSITY_KG_M3..=Self::MAX_BULK_DENSITY_KG_M3).contains(&density) {
            return Err(ParamError::Inconsistent {
                name: "bulk_density_kg_m3",
                value: density,
                fields: ["radius_m", "surface_gravity_m_s2"],
            });
        }
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

// Every planet within the bounds rotates slower than break-up (`ω² R < g`
// with `ω = 2π / P`, so `4π² R < g P²`), worst at the largest radius, the
// lowest gravity and the shortest period: 5.1e8 < 8.3e8. `validate` relies
// on it instead of checking it.
const _: () = assert!(
    4.0 * PI * PI * PlanetParams::MAX_RADIUS_M
        < PlanetParams::MIN_SURFACE_GRAVITY_M_S2
            * PlanetParams::MIN_ROTATION_PERIOD_S
            * PlanetParams::MIN_ROTATION_PERIOD_S,
    "the bounds of PlanetParams admit a rotation faster than break-up"
);

/// Why a set of [`PlanetParams`], or of [`Settings`](crate::Settings),
/// cannot be generated.
#[derive(Clone, Debug, PartialEq)]
pub enum ParamError {
    /// A value is not finite or outside its valid range.
    OutOfRange {
        /// Field name, of [`PlanetParams`] or of
        /// [`Settings`](crate::Settings).
        name: &'static str,
        /// The rejected value.
        value: f64,
    },
    /// Parameters, each within its bounds, that together describe a planet
    /// that cannot exist.
    Inconsistent {
        /// The derived quantity that is out of range, named with its unit
        /// (`"bulk_density_kg_m3"`).
        name: &'static str,
        /// Its value.
        value: f64,
        /// The fields of [`PlanetParams`] it is derived from, in field order
        /// (`["radius_m", "surface_gravity_m_s2"]`).
        fields: [&'static str; 2],
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
            Self::Inconsistent {
                name,
                value,
                fields: [first, second],
            } => write!(
                f,
                "{name} from {first} and {second} is out of range: {value}"
            ),
            Self::Unsupported { feature, planned } => {
                write!(f, "{feature} is not supported yet (planned for {planned})")
            }
        }
    }
}

impl std::error::Error for ParamError {}

#[cfg(test)]
mod tests {
    use kuniumi_math::units;

    use super::*;

    /// The next double above the finite `v` (`f64::next_up` needs Rust
    /// 1.86, above the MSRV).
    fn next_up(v: f64) -> f64 {
        if v == 0.0 {
            f64::from_bits(1)
        } else if v > 0.0 {
            f64::from_bits(v.to_bits() + 1)
        } else {
            f64::from_bits(v.to_bits() - 1)
        }
    }

    /// The next double below the finite `v`.
    fn next_down(v: f64) -> f64 {
        -next_up(-v)
    }

    fn earth() -> PlanetParams {
        PlanetParams::earth()
    }

    /// Earth with `radius_m`, and the gravity that keeps Earth's density
    /// (`g ∝ R` at a fixed density), so that only the radius is tested.
    fn with_radius(radius_m: f64) -> PlanetParams {
        let e = earth();
        PlanetParams {
            radius_m,
            surface_gravity_m_s2: e.surface_gravity_m_s2 * radius_m / e.radius_m,
            ..e
        }
    }

    /// Earth with `g` and the radius that keeps Earth's density, within the
    /// radius bounds.
    fn with_gravity(g: f64) -> PlanetParams {
        let e = earth();
        let radius_m = (e.radius_m * g / e.surface_gravity_m_s2)
            .clamp(PlanetParams::MIN_RADIUS_M, PlanetParams::MAX_RADIUS_M);
        PlanetParams {
            radius_m: if radius_m.is_finite() {
                radius_m
            } else {
                e.radius_m
            },
            surface_gravity_m_s2: g,
            ..e
        }
    }

    /// Each bounded field: its name, a planet with the field set to a value
    /// (and otherwise valid where it can be), and its bounds.
    #[allow(clippy::type_complexity)]
    fn bounded() -> [(&'static str, fn(f64) -> PlanetParams, f64, f64); 7] {
        [
            (
                "radius_m",
                with_radius,
                PlanetParams::MIN_RADIUS_M,
                PlanetParams::MAX_RADIUS_M,
            ),
            (
                "surface_gravity_m_s2",
                with_gravity,
                PlanetParams::MIN_SURFACE_GRAVITY_M_S2,
                PlanetParams::MAX_SURFACE_GRAVITY_M_S2,
            ),
            (
                "rotation_period_s",
                |v| PlanetParams {
                    rotation_period_s: v,
                    ..earth()
                },
                PlanetParams::MIN_ROTATION_PERIOD_S,
                PlanetParams::MAX_ROTATION_PERIOD_S,
            ),
            (
                "orbital_period_s",
                |v| PlanetParams {
                    orbital_period_s: v,
                    ..earth()
                },
                PlanetParams::MIN_ORBITAL_PERIOD_S,
                PlanetParams::MAX_ORBITAL_PERIOD_S,
            ),
            (
                "stellar_flux_w_m2",
                |v| PlanetParams {
                    stellar_flux_w_m2: v,
                    ..earth()
                },
                PlanetParams::MIN_STELLAR_FLUX_W_M2,
                PlanetParams::MAX_STELLAR_FLUX_W_M2,
            ),
            (
                "ocean_volume_m3",
                |v| PlanetParams {
                    ocean_volume_m3: v,
                    ..earth()
                },
                0.0,
                PlanetParams::MAX_OCEAN_VOLUME_M3,
            ),
            (
                "surface_pressure_pa",
                |v| PlanetParams {
                    surface_pressure_pa: v,
                    ..earth()
                },
                PlanetParams::MIN_SURFACE_PRESSURE_PA,
                PlanetParams::MAX_SURFACE_PRESSURE_PA,
            ),
        ]
    }

    #[test]
    fn earth_is_valid() {
        assert_eq!(PlanetParams::earth().validate(), Ok(()));
    }

    #[test]
    fn each_bounded_parameter_accepts_its_bounds_and_rejects_the_values_beyond() {
        for (name, planet, lo, hi) in bounded() {
            for v in [lo, hi] {
                let p = planet(v);
                if name == "surface_gravity_m_s2" && v.to_bits() == lo.to_bits() {
                    // 1 m/s² passes its own bounds, but at every radius
                    // gives a density below 3,000 kg/m³: at the smallest
                    // radius 1 / ((4/3) π G · 1.5e6) = 2,385 kg/m³ (hand
                    // calculation).
                    assert_eq!(p.radius_m.to_bits(), PlanetParams::MIN_RADIUS_M.to_bits());
                    assert!(
                        matches!(
                            p.validate(),
                            Err(ParamError::Inconsistent {
                                name: "bulk_density_kg_m3",
                                ..
                            })
                        ),
                        "{p:?}"
                    );
                } else {
                    assert_eq!(p.validate(), Ok(()), "{name} = {v:e}: {p:?}");
                }
            }
            if name == "ocean_volume_m3" {
                assert_eq!(planet(-0.0).validate(), Ok(()));
            }
            for v in [
                next_down(lo),
                next_up(hi),
                f64::NAN,
                f64::INFINITY,
                f64::NEG_INFINITY,
            ] {
                let p = planet(v);
                match p.validate() {
                    Err(ParamError::OutOfRange { name: n, value }) => {
                        assert_eq!(n, name, "{v:e}: {p:?}");
                        assert_eq!(value.to_bits(), v.to_bits(), "{name}");
                    }
                    other => panic!("{name} = {v:e}: {other:?}"),
                }
            }
        }
    }

    #[test]
    fn the_bounds_hold_earth_and_are_round_in_human_units() {
        let e = earth();
        let earth_values = [
            e.radius_m,
            e.surface_gravity_m_s2,
            e.rotation_period_s,
            e.orbital_period_s,
            e.stellar_flux_w_m2,
            e.ocean_volume_m3,
            e.surface_pressure_pa,
        ];
        for ((name, _, lo, hi), v) in bounded().into_iter().zip(earth_values) {
            assert!(lo < v && v < hi, "{name}: {v}");
        }
        let rho = e.bulk_density_kg_m3();
        assert!(
            PlanetParams::MIN_BULK_DENSITY_KG_M3 < rho
                && rho < PlanetParams::MAX_BULK_DENSITY_KG_M3
        );
        // Each bound is the exact conversion of a round number in the
        // settings file's human unit (docs/conventions.md, Units).
        for (si, human) in [
            (PlanetParams::MIN_RADIUS_M, units::km_to_m(1_500.0)),
            (PlanetParams::MAX_RADIUS_M, units::km_to_m(13_000.0)),
            (PlanetParams::MIN_ROTATION_PERIOD_S, units::hours_to_s(8.0)),
            (
                PlanetParams::MAX_ROTATION_PERIOD_S,
                units::hours_to_s(240.0),
            ),
            (PlanetParams::MIN_ORBITAL_PERIOD_S, units::days_to_s(90.0)),
            (
                PlanetParams::MAX_ORBITAL_PERIOD_S,
                units::days_to_s(3_650.0),
            ),
            (PlanetParams::MAX_OCEAN_VOLUME_M3, units::km3_to_m3(1e11)),
            (PlanetParams::MIN_SURFACE_PRESSURE_PA, units::bar_to_pa(0.1)),
            (
                PlanetParams::MAX_SURFACE_PRESSURE_PA,
                units::bar_to_pa(10.0),
            ),
        ] {
            assert_eq!(si.to_bits(), human.to_bits(), "{si:e}");
        }
    }

    #[test]
    fn the_bulk_density_comes_from_gravity_and_radius() {
        // Python: 3 * 9.80665 / (4 * math.pi * 6.6743e-11 * 6.371e6).
        let earth_density = 5_505.778_525_265_396;
        let rho = earth().bulk_density_kg_m3();
        assert!(
            ((rho - earth_density) / earth_density).abs() < 1e-12,
            "{rho}"
        );
        // Within 0.2 % of Earth's mean density from its mass, 5,513 kg/m³
        // (5.972e24 kg over a sphere of 6,371 km).
        assert!(((rho - 5_513.0) / 5_513.0).abs() < 2e-3, "{rho}");
        // Twice the gravity, twice the density; twice the radius, half.
        let p = PlanetParams {
            surface_gravity_m_s2: 2.0 * 9.806_65,
            ..earth()
        };
        assert_eq!(p.bulk_density_kg_m3().to_bits(), (2.0 * rho).to_bits());
        let p = PlanetParams {
            radius_m: 2.0 * 6_371_000.0,
            ..earth()
        };
        assert_eq!(p.bulk_density_kg_m3().to_bits(), (0.5 * rho).to_bits());
    }

    /// The gravity at `radius_m` on the density bound `bound`: the lowest
    /// gravity whose density is at least `bound` when `lower`, otherwise
    /// the highest whose density is at most `bound`.
    fn gravity_at_density(radius_m: f64, bound: f64, lower: bool) -> f64 {
        let density = |g: f64| {
            PlanetParams {
                radius_m,
                surface_gravity_m_s2: g,
                ..earth()
            }
            .bulk_density_kg_m3()
        };
        let mut g = bound * (4.0 / 3.0) * PI * GRAVITATIONAL_CONSTANT_M3_KG_S2 * radius_m;
        if lower {
            while density(g) < bound {
                g = next_up(g);
            }
            while density(next_down(g)) >= bound {
                g = next_down(g);
            }
        } else {
            while density(g) > bound {
                g = next_down(g);
            }
            while density(next_up(g)) <= bound {
                g = next_up(g);
            }
        }
        g
    }

    #[test]
    fn a_density_out_of_range_is_inconsistent() {
        let (lo, hi) = (
            PlanetParams::MIN_BULK_DENSITY_KG_M3,
            PlanetParams::MAX_BULK_DENSITY_KG_M3,
        );
        let earth_radius_m = earth().radius_m;
        for (radius_m, bound, lower) in [
            (PlanetParams::MIN_RADIUS_M, lo, true),
            (PlanetParams::MIN_RADIUS_M, hi, false),
            (earth_radius_m, lo, true),
            (earth_radius_m, hi, false),
            (PlanetParams::MAX_RADIUS_M, lo, true),
        ] {
            let g = gravity_at_density(radius_m, bound, lower);
            let at = PlanetParams {
                radius_m,
                surface_gravity_m_s2: g,
                ..earth()
            };
            assert_eq!(at.validate(), Ok(()), "{at:?}");
            let beyond = PlanetParams {
                surface_gravity_m_s2: if lower { next_down(g) } else { next_up(g) },
                ..at
            };
            assert_eq!(
                beyond.validate(),
                Err(ParamError::Inconsistent {
                    name: "bulk_density_kg_m3",
                    value: beyond.bulk_density_kg_m3(),
                    fields: ["radius_m", "surface_gravity_m_s2"],
                }),
                "{beyond:?}"
            );
        }
        // At the largest radius the gravity bound comes first: 40 m/s² is
        // 11,006 kg/m³ there (hand calculation).
        let p = PlanetParams {
            radius_m: PlanetParams::MAX_RADIUS_M,
            surface_gravity_m_s2: PlanetParams::MAX_SURFACE_GRAVITY_M_S2,
            ..earth()
        };
        assert_eq!(p.validate(), Ok(()));
        assert!((p.bulk_density_kg_m3() - 11_006.0).abs() < 1.0);
    }

    #[test]
    fn checks_run_fields_first_then_the_density_then_the_regimes() {
        // A gravity out of its bounds is named, though the density is out
        // of range too.
        let p = PlanetParams {
            radius_m: PlanetParams::MAX_RADIUS_M,
            surface_gravity_m_s2: 0.5,
            ..earth()
        };
        assert!(matches!(
            p.validate(),
            Err(ParamError::OutOfRange {
                name: "surface_gravity_m_s2",
                ..
            })
        ));
        // Earth's gravity on the largest planet: 2,698 kg/m³ (Python:
        // 3 * 9.80665 / (4 * math.pi * 6.6743e-11 * 1.3e7)), before the
        // unsupported regime.
        let p = PlanetParams {
            radius_m: PlanetParams::MAX_RADIUS_M,
            rotation: RotationRegime::TidallyLocked,
            ..earth()
        };
        match p.validate() {
            Err(ParamError::Inconsistent {
                name,
                value,
                fields,
            }) => {
                assert_eq!(name, "bulk_density_kg_m3");
                assert!((value - 2_698.254_998_805_064).abs() < 1e-9, "{value}");
                assert_eq!(fields, ["radius_m", "surface_gravity_m_s2"]);
            }
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn every_planet_within_the_bounds_rotates_slower_than_break_up() {
        // Break-up: `ω² R = g`, a period of 2π √(R / g), longest at the
        // largest radius and the lowest gravity. Python:
        // 2 * math.pi * math.sqrt(1.3e7 / 1.0) = 22,654.3 s (6.29 h).
        let slowest_break_up_s =
            2.0 * PI * (PlanetParams::MAX_RADIUS_M / PlanetParams::MIN_SURFACE_GRAVITY_M_S2).sqrt();
        assert!((slowest_break_up_s - 22_654.346_798_277_955).abs() < 1e-6);
        assert!(slowest_break_up_s < PlanetParams::MIN_ROTATION_PERIOD_S);
        // With the density bound, break-up is at most √(3π / (G ρ)):
        // 6,861 s (1.9 h) at 3,000 kg/m³ (hand calculation).
        let densest_break_up_s = (3.0 * PI
            / (GRAVITATIONAL_CONSTANT_M3_KG_S2 * PlanetParams::MIN_BULK_DENSITY_KG_M3))
            .sqrt();
        assert!(
            (densest_break_up_s - 6_861.0).abs() < 1.0,
            "{densest_break_up_s}"
        );
    }

    #[test]
    fn inconsistent_parameters_are_displayed_with_their_fields() {
        let e = ParamError::Inconsistent {
            name: "bulk_density_kg_m3",
            value: 275.0,
            fields: ["radius_m", "surface_gravity_m_s2"],
        };
        assert_eq!(
            e.to_string(),
            "bulk_density_kg_m3 from radius_m and surface_gravity_m_s2 is out of range: 275"
        );
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

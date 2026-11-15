use kuniumi_math::units::SECONDS_PER_YEAR;

use crate::ParamError;

/// How a planet is generated, as opposed to what it is ([`PlanetParams`]):
/// the resolution and the settings of each stage.
///
/// [`PlanetParams`]: crate::PlanetParams
#[derive(Clone, Debug, PartialEq)]
pub struct Settings {
    /// Number of cells of the global grid (design document,
    /// 球面グリッドとデータ構造; D-14), in
    /// [`MIN_CELLS`](Self::MIN_CELLS)`..=`[`MAX_CELLS`](Self::MAX_CELLS).
    /// The default, 650,000, puts the cell centres about 28 km apart on an
    /// Earth-sized planet.
    pub cells: u32,
    /// Number of tectonic plates the Tectonics stage divides the surface
    /// into, in [`MIN_PLATE_COUNT`](Self::MIN_PLATE_COUNT)`..=`
    /// [`MAX_PLATE_COUNT`](Self::MAX_PLATE_COUNT). The default, 12, is in
    /// the 10 to 15 plates the design starts from (design document,
    /// プレートテクトニクス / 初期条件の目安); Earth has seven large plates
    /// and many small ones.
    pub plate_count: u32,
    /// What `kuniumi_sim::generate` does at the first stage that is not
    /// implemented yet: with `false` (the default) it fails with
    /// `SimError::NotImplemented`; with `true` it returns the world built so
    /// far, whose [`World::stages`](crate::World::stages) end before that
    /// stage, so that the stages that exist can be looked at during
    /// development. No stage reads it, so it changes no generated value; it
    /// is recorded with the other settings ([`World::settings`]).
    ///
    /// [`World::settings`]: crate::World::settings
    pub allow_partial: bool,
}

impl Settings {
    /// The default number of cells: 650,000, the resolution of the design.
    pub const DEFAULT_CELLS: u32 = 650_000;

    /// The fewest cells a world may have: 1,000 (about 714 km apart on an
    /// Earth-sized planet), below which no stage is meaningful.
    pub const MIN_CELLS: u32 = 1_000;

    /// The most cells a world may have: [`kuniumi_sphere::MAX_CELLS`]
    /// (715,827,884), the largest grid whose neighbour lists fit `u32`
    /// offsets.
    pub const MAX_CELLS: u32 = {
        assert!(kuniumi_sphere::MAX_CELLS <= u32::MAX as usize);
        kuniumi_sphere::MAX_CELLS as u32
    };

    /// The default number of plates: 12.
    pub const DEFAULT_PLATE_COUNT: u32 = 12;

    /// The fewest plates: 2, the fewest that have a boundary.
    pub const MIN_PLATE_COUNT: u32 = 2;

    /// The most plates: 64, far below [`MIN_CELLS`](Self::MIN_CELLS), so
    /// that every plate can have a cell of its own.
    pub const MAX_PLATE_COUNT: u32 = 64;

    /// Checks every setting, in the order of the fields; the error names the
    /// first one out of range ([`ParamError::OutOfRange`] with the setting's
    /// name, such as `"cells"` or `"plate_count"`).
    pub fn validate(&self) -> Result<(), ParamError> {
        if !(Self::MIN_CELLS..=Self::MAX_CELLS).contains(&self.cells) {
            return Err(ParamError::OutOfRange {
                name: "cells",
                value: f64::from(self.cells),
            });
        }
        if !(Self::MIN_PLATE_COUNT..=Self::MAX_PLATE_COUNT).contains(&self.plate_count) {
            return Err(ParamError::OutOfRange {
                name: "plate_count",
                value: f64::from(self.plate_count),
            });
        }
        Ok(())
    }
}

impl Default for Settings {
    fn default() -> Self {
        Self {
            cells: Self::DEFAULT_CELLS,
            plate_count: Self::DEFAULT_PLATE_COUNT,
            allow_partial: false,
        }
    }
}

/// The settings of the Erosion stage (design document, 侵食と地形発達 /
/// 基本方程式): the stream-power law `E = K Qᵐ S` (n = 1) of the river
/// erosion (ticket M3-03), with the erosion rate `E` (m/s), the discharge
/// `Q` (m³/s) and the slope `S` to the receiver (m/m), solved implicitly
/// (Braun & Willett 2013), and the hillslope diffusion `K_d ∇²h` (ticket
/// M3-04), solved by backward Euler; [`steps`](Self::steps) steps of
/// [`time_step_s`](Self::time_step_s).
///
/// Not yet a field of [`Settings`]: ticket M3-07 adds it there and records
/// it in `world.json` and in settings files. Until then the solver of
/// `kuniumi_sim::stages::erosion` takes it as an argument.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ErosionSettings {
    /// The erodibility `K` of the stream-power law, in the SI unit
    /// m^(1−3m) s^(m−1) that makes `K Qᵐ S` an erosion rate in m/s for a
    /// discharge in m³/s; the unit depends on
    /// [`discharge_exponent`](Self::discharge_exponent), so the name has no
    /// suffix. Finite and positive. Default
    /// [`DEFAULT_ERODIBILITY`](Self::DEFAULT_ERODIBILITY).
    pub erodibility: f64,
    /// The exponent `m` of the discharge (dimensionless), finite and not
    /// negative. Default
    /// [`DEFAULT_DISCHARGE_EXPONENT`](Self::DEFAULT_DISCHARGE_EXPONENT).
    pub discharge_exponent: f64,
    /// The hillslope diffusivity `K_d` (m²/s) of the diffusion
    /// `∂h/∂t = K_d ∇²h`, finite and not negative; 0 turns the diffusion
    /// off. Default
    /// [`DEFAULT_DIFFUSIVITY_M2_S`](Self::DEFAULT_DIFFUSIVITY_M2_S).
    pub diffusivity_m2_s: f64,
    /// The time step `Δt` (s), finite and positive. Default
    /// [`DEFAULT_TIME_STEP_S`](Self::DEFAULT_TIME_STEP_S).
    pub time_step_s: f64,
    /// The number of time steps; 0 leaves the surface unchanged. Default
    /// [`DEFAULT_STEPS`](Self::DEFAULT_STEPS).
    pub steps: u32,
}

impl ErosionSettings {
    /// The default erodibility: 7.5 × 10⁻¹⁰ m^(1−3m) s^(m−1) (with
    /// `m` = 0.45: m^−0.35 s^−0.55).
    ///
    /// The design tunes `K` by its results, since its value depends on the
    /// units and the grid spacing: the major orogens 2–4 km high on
    /// average, the highest peaks 8–9 km (design document, 侵食と地形発達 /
    /// 調整の目安). At steady state the slope is `S = U / (K Qᵐ)`. Take a
    /// river crossing a belt that rises at the rate `U`, on the default
    /// grid (650,000 cells: `a` = 7.85 × 10⁸ m² per cell, `L` = 28.0 km
    /// between centres, `√(4π/n) R`), with 1 m/yr of runoff `P`, its
    /// drainage area growing as in Hack's law, `A_k = a k^(1/0.6)` at the
    /// k-th cell from the divide. Its divide stands at
    /// `Σ_k U L / (K (P A_k)ᵐ)` above the foreland: 3.0–3.8 km for 5–8
    /// cells (140–220 km) at 4 mm/yr, the uplift of a collision front
    /// converging at 4 cm/yr (`COLLISION_UPLIFT_RATIO` of the Tectonics
    /// stage, 0.1), and 7.5–9.4 km at 10 mm/yr, the fastest uplift of the
    /// Tectonics stage (ticket M3-03 has the numbers). The same law written
    /// for drainage area instead of discharge, `E = K_A Aᵐ S` with
    /// `K = K_A / Pᵐ`, has `K_A` = 1.0 × 10⁻⁵ m^0.1 yr⁻¹. Calibration
    /// against the metrics of the design is left to tickets M3-05 and M7.
    pub const DEFAULT_ERODIBILITY: f64 = 7.5e-10;

    /// The default discharge exponent: 0.45, in the 0.4–0.5 of the design
    /// document (侵食と地形発達 / 基本方程式).
    pub const DEFAULT_DISCHARGE_EXPONENT: f64 = 0.45;

    /// The default hillslope diffusivity: 5 × 10⁻³ m²/yr
    /// (1.584 × 10⁻¹⁰ m²/s), an Earth value for soil-mantled hillslopes:
    /// the mean of 49 ± 37 cm²/yr that McKean et al. (1993) computed for 34
    /// sites in California, Oregon and Washington from the rates of
    /// colluvium accumulation in hollows reported by Reneau (1988), as
    /// summarised by Dietrich et al. (2003, "Geomorphic transport laws for
    /// predicting landscape form and dynamics", AGU Geophysical Monograph
    /// 135, 103–132), who also cite 3.2 × 10⁻³ m²/yr for the Oregon Coast
    /// Range (Roering et al. 1999) as the coefficient of the nonlinear,
    /// critical-slope law, the linear diffusivity at low gradients; within
    /// the 4.4 × 10⁻⁴ to 3.6 × 10⁻² m²/yr of the field estimates
    /// summarised by Fernandes & Dietrich (1997), as quoted by Ruetenik et
    /// al. (2023, *Earth Surface Dynamics* 11, 865–880).
    ///
    /// On the default grid (cell centres 28 km apart) this acts below the
    /// grid scale: over the 10 Myr of the default run the diffusion length
    /// `√(K_d t)` is 220 m, so the diffusion changes the heights little
    /// (ticket M3-04). Calibration against the metrics of the design is
    /// left to tickets M3-05 and M7.
    pub const DEFAULT_DIFFUSIVITY_M2_S: f64 = 5.0e-3 / SECONDS_PER_YEAR;

    /// The default time step: 10,000 years (3.15576 × 10¹¹ s), the lower
    /// end of the 10,000 to 100,000 years that the design allows the
    /// implicit solve (侵食と地形発達 / 基本方程式). On the default grid, the
    /// factor `F = Δt K Qᵐ / L` of the solve then runs from 0.036 for the
    /// runoff of one cell to 1.5 for 10⁵ m³/s, so that knickpoints move
    /// up the rivers without much numerical smearing; longer steps stay
    /// stable.
    pub const DEFAULT_TIME_STEP_S: f64 = 1.0e4 * SECONDS_PER_YEAR;

    /// The default number of steps: 1,000, the lower end of the design's
    /// 1,000–3,000 (侵食と地形発達), so 10 Myr in all with the default time
    /// step: a few Myr for each of the two or three rounds of erosion
    /// between which the climate is computed again (侵食と地形発達 /
    /// 気候との交互反復).
    pub const DEFAULT_STEPS: u32 = 1_000;

    /// Checks every setting, in the order of the fields; the error names
    /// the first one that is not finite or out of its domain
    /// ([`ParamError::OutOfRange`] with the field's name, such as
    /// `"erodibility"`): the erodibility and the time step must be finite
    /// and positive, the discharge exponent and the diffusivity finite and
    /// not negative; every number of steps is valid. Within these the laws
    /// are defined and the solvers' results are finite, provided the step
    /// of the diffusion, `Δt K_d / R²` with the planet's radius `R`, does
    /// not overflow, which the solver checks (`kuniumi_sim`'s `erode`).
    pub fn validate(&self) -> Result<(), ParamError> {
        // Written so that NaN fails every check.
        let checks = [
            (
                "erodibility",
                self.erodibility,
                self.erodibility > 0.0 && self.erodibility.is_finite(),
            ),
            (
                "discharge_exponent",
                self.discharge_exponent,
                self.discharge_exponent >= 0.0 && self.discharge_exponent.is_finite(),
            ),
            (
                "diffusivity_m2_s",
                self.diffusivity_m2_s,
                self.diffusivity_m2_s >= 0.0 && self.diffusivity_m2_s.is_finite(),
            ),
            (
                "time_step_s",
                self.time_step_s,
                self.time_step_s > 0.0 && self.time_step_s.is_finite(),
            ),
        ];
        match checks.into_iter().find(|&(_, _, valid)| !valid) {
            Some((name, value, _)) => Err(ParamError::OutOfRange { name, value }),
            None => Ok(()),
        }
    }
}

impl Default for ErosionSettings {
    fn default() -> Self {
        Self {
            erodibility: Self::DEFAULT_ERODIBILITY,
            discharge_exponent: Self::DEFAULT_DISCHARGE_EXPONENT,
            diffusivity_m2_s: Self::DEFAULT_DIFFUSIVITY_M2_S,
            time_step_s: Self::DEFAULT_TIME_STEP_S,
            steps: Self::DEFAULT_STEPS,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn erosion_defaults_are_the_documented_values_and_valid() {
        let s = ErosionSettings::default();
        assert_eq!(s.erodibility.to_bits(), 7.5e-10_f64.to_bits());
        assert_eq!(s.discharge_exponent.to_bits(), 0.45_f64.to_bits());
        // 5 × 10⁻³ m²/yr over the Julian year, 1.584 × 10⁻¹⁰ m²/s.
        assert_eq!(
            s.diffusivity_m2_s.to_bits(),
            (5.0e-3 / 31_557_600.0_f64).to_bits()
        );
        assert!((s.diffusivity_m2_s - 1.584_404e-10).abs() < 1e-15);
        // 10,000 Julian years, an integer number of seconds.
        assert_eq!(s.time_step_s.to_bits(), 315_576_000_000.0_f64.to_bits());
        assert_eq!(s.steps, 1_000);
        assert_eq!(s.validate(), Ok(()));
    }

    #[test]
    fn erosion_settings_accept_their_whole_domain() {
        for (erodibility, discharge_exponent, diffusivity_m2_s, time_step_s, steps) in [
            (f64::MIN_POSITIVE, 0.0, 0.0, f64::MIN_POSITIVE, 0),
            (5e-324, 1.0, 5e-324, 5e-324, 1),
            (1.0, 0.5, 1.0, 1.0, u32::MAX),
            (f64::MAX, f64::MAX, f64::MAX, f64::MAX, 7),
        ] {
            let s = ErosionSettings {
                erodibility,
                discharge_exponent,
                diffusivity_m2_s,
                time_step_s,
                steps,
            };
            assert_eq!(s.validate(), Ok(()), "{s:?}");
        }
    }

    #[test]
    fn erosion_settings_reject_values_out_of_their_domain_by_name() {
        let bad_positive = [0.0, -0.0, -1.0, f64::INFINITY, f64::NEG_INFINITY, f64::NAN];
        let bad_non_negative = [-5e-324, -1.0, f64::INFINITY, f64::NEG_INFINITY, f64::NAN];
        let cases = bad_positive
            .iter()
            .map(|&v| {
                (
                    ErosionSettings {
                        erodibility: v,
                        ..ErosionSettings::default()
                    },
                    "erodibility",
                    v,
                )
            })
            .chain(bad_non_negative.iter().map(|&v| {
                (
                    ErosionSettings {
                        discharge_exponent: v,
                        ..ErosionSettings::default()
                    },
                    "discharge_exponent",
                    v,
                )
            }))
            .chain(bad_non_negative.iter().map(|&v| {
                (
                    ErosionSettings {
                        diffusivity_m2_s: v,
                        ..ErosionSettings::default()
                    },
                    "diffusivity_m2_s",
                    v,
                )
            }))
            .chain(bad_positive.iter().map(|&v| {
                (
                    ErosionSettings {
                        time_step_s: v,
                        ..ErosionSettings::default()
                    },
                    "time_step_s",
                    v,
                )
            }));
        for (s, name, value) in cases {
            match s.validate() {
                Err(ParamError::OutOfRange {
                    name: n,
                    value: got,
                }) => {
                    assert_eq!(n, name, "{s:?}");
                    assert_eq!(got.to_bits(), value.to_bits(), "{s:?}");
                }
                other => panic!("{s:?}: {other:?}"),
            }
        }
        // The fields are checked in order.
        let all_bad = ErosionSettings {
            erodibility: -1.0,
            discharge_exponent: -1.0,
            diffusivity_m2_s: -1.0,
            time_step_s: -1.0,
            steps: 0,
        };
        assert_eq!(
            all_bad.validate(),
            Err(ParamError::OutOfRange {
                name: "erodibility",
                value: -1.0,
            })
        );
        let two_bad = ErosionSettings {
            erodibility: 1.0,
            ..all_bad
        };
        assert_eq!(
            two_bad.validate(),
            Err(ParamError::OutOfRange {
                name: "discharge_exponent",
                value: -1.0,
            })
        );
        let three_bad = ErosionSettings {
            discharge_exponent: 1.0,
            ..two_bad
        };
        assert_eq!(
            three_bad.validate(),
            Err(ParamError::OutOfRange {
                name: "diffusivity_m2_s",
                value: -1.0,
            })
        );
    }

    #[test]
    fn default_has_650_000_cells_12_plates_and_is_valid() {
        let s = Settings::default();
        assert_eq!(s.cells, 650_000);
        assert_eq!(s.plate_count, 12);
        assert!(!s.allow_partial);
        assert_eq!(s.validate(), Ok(()));
    }

    #[test]
    fn plate_count_bounds_are_inclusive() {
        assert_eq!(
            (Settings::MIN_PLATE_COUNT, Settings::MAX_PLATE_COUNT),
            (2, 64)
        );
        for plate_count in [2, 3, 12, 63, 64] {
            let s = Settings {
                plate_count,
                ..Settings::default()
            };
            assert_eq!(s.validate(), Ok(()), "{plate_count}");
        }
    }

    #[test]
    fn rejects_plate_counts_out_of_range_by_name() {
        for plate_count in [0, 1, 65, 1_000, u32::MAX] {
            assert_eq!(
                Settings {
                    plate_count,
                    ..Settings::default()
                }
                .validate(),
                Err(ParamError::OutOfRange {
                    name: "plate_count",
                    value: f64::from(plate_count),
                }),
                "{plate_count}"
            );
        }
        // The cells are checked first.
        assert_eq!(
            Settings {
                cells: 999,
                plate_count: 0,
                allow_partial: false,
            }
            .validate(),
            Err(ParamError::OutOfRange {
                name: "cells",
                value: 999.0,
            })
        );
    }

    #[test]
    fn cell_bounds_are_inclusive() {
        // 715,827,884 is the largest n with 6n − 12 ≤ u32::MAX (M1-06).
        assert_eq!(Settings::MAX_CELLS, 715_827_884);
        for cells in [1_000, 1_001, 10_000, 715_827_883, 715_827_884] {
            for allow_partial in [false, true] {
                let s = Settings {
                    cells,
                    allow_partial,
                    ..Settings::default()
                };
                assert_eq!(s.validate(), Ok(()), "{cells}");
            }
        }
    }

    #[test]
    fn rejects_cells_out_of_range_by_name() {
        for cells in [0, 4, 999, 715_827_885, u32::MAX] {
            assert_eq!(
                Settings {
                    cells,
                    ..Settings::default()
                }
                .validate(),
                Err(ParamError::OutOfRange {
                    name: "cells",
                    value: f64::from(cells),
                }),
                "{cells}"
            );
        }
    }
}

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

#[cfg(test)]
mod tests {
    use super::*;

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

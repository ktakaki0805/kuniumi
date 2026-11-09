use crate::ParamError;

/// How a planet is generated, as opposed to what it is ([`PlanetParams`]):
/// the resolution and, in later milestones, the settings of each stage.
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

    /// Checks every setting; the error names the first one out of range
    /// ([`ParamError::OutOfRange`] with the setting's name, such as
    /// `"cells"`).
    pub fn validate(&self) -> Result<(), ParamError> {
        if !(Self::MIN_CELLS..=Self::MAX_CELLS).contains(&self.cells) {
            return Err(ParamError::OutOfRange {
                name: "cells",
                value: f64::from(self.cells),
            });
        }
        Ok(())
    }
}

impl Default for Settings {
    fn default() -> Self {
        Self {
            cells: Self::DEFAULT_CELLS,
            allow_partial: false,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_has_650_000_cells_and_is_valid() {
        let s = Settings::default();
        assert_eq!(s.cells, 650_000);
        assert!(!s.allow_partial);
        assert_eq!(s.validate(), Ok(()));
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

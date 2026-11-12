//! Crust to elevation (design: "地殻から標高へ"; ticket M2-05).
//!
//! The Crust stage turns the crust of the Tectonics stage into relief and
//! fills the basins with the planet's water:
//!
//! - ocean floor sinks as it ages, after Parsons & Sclater (1977)
//!   ([`OceanDepth`]); trenches lie deeper still, at [`TRENCH_DEPTH_M`];
//! - continental crust floats on the mantle by Airy isostasy, higher where
//!   it is thicker ([`continental_height_m`]);
//! - the sea level is the level at which the basins hold the planet's
//!   water, [`PlanetParams::ocean_volume_m3`] ([`sea_level_m`]). The land
//!   fraction is a result, never an input.
//!
//! Heights are relative to the datum of the reference radius, on which
//! continental crust of the reference thickness (35 km) stands at
//! [`REFERENCE_HEIGHT_M`] and the ocean floor at its depth below it; the sea
//! level is solved afterwards, on the same datum (docs/conventions.md,
//! "Elevation, sea level and time"). For Earth's water it lies 0.7–1.0 km
//! below the datum (ticket M2-05).
//!
//! Outputs ([`run`]): [`SURFACE_HEIGHT`] (m) and `World::sea_level_m` (m).
//! The stage draws no random numbers.
//!
//! [`PlanetParams::ocean_volume_m3`]: kuniumi_world::PlanetParams::ocean_volume_m3

use kuniumi_math::det;
use kuniumi_math::units::SECONDS_PER_MYR;
use kuniumi_sphere::SphereGrid;
use kuniumi_world::World;

use crate::SimError;
use crate::stages::tectonics::{
    CONTINENTAL_THICKNESS_M, CRUST_AGE, CRUST_THICKNESS, CRUST_TYPE, TRENCH,
};

/// The field of the height of the solid surface (m) above the datum of the
/// reference radius ([`surface_heights_m`]).
pub const SURFACE_HEIGHT: &str = "surface_height";

/// The density of continental crust (kg/m³): 2,800 kg/m³ (design,
/// 地殻から標高へ / 大陸).
pub const CRUST_DENSITY_KG_M3: f64 = 2_800.0;

/// The density of the mantle under the crust (kg/m³): 3,300 kg/m³ (design,
/// 地殻から標高へ / 大陸).
pub const MANTLE_DENSITY_KG_M3: f64 = 3_300.0;

/// The height (m) of continental crust of the reference thickness
/// [`CONTINENTAL_THICKNESS_M`] (35 km): +300 m, so that standard continental
/// crust stands about 300 m above the datum (design, 地殻から標高へ / 大陸).
pub const REFERENCE_HEIGHT_M: f64 = 300.0;

/// The depth (m) below the datum of the floor of a trench: 8 km. Trench
/// cells (the oceanic side of a subduction boundary, `trench` of the
/// Tectonics stage) lie at the deeper of this and the depth of their age,
/// and the ocean floor of any age is shallower than 6.4 km, so every trench
/// lies at −8 km, 7.0–7.3 km below the sea level of Earth's water (ticket
/// M2-05, seeds 1–10 at 10,000 and 650,000 cells). The middle of the
/// 7–10 km of the ticket, which most of Earth's major trenches reach (the
/// deepest, the Mariana trench, about 10.9 km), and the value of ticket
/// M2-11's prototype; a trench that depends on the age of the subducting
/// floor is left to the time-evolving tectonics (M6-03).
pub const TRENCH_DEPTH_M: f64 = 8_000.0;

/// The depth (m) of the ridge crest, `t = 0` in the young curve of Parsons
/// & Sclater (1977): 2,500 m.
const RIDGE_DEPTH_M: f64 = 2_500.0;

/// The growth of the depth with the square root of the age, in the young
/// curve: 350 m per √Myr.
const YOUNG_DEPTH_SLOPE_M: f64 = 350.0;

/// The depth (m) that old ocean floor approaches, in the old curve:
/// 6,400 m.
const OLD_DEPTH_LIMIT_M: f64 = 6_400.0;

/// The depth (m) by which old ocean floor is shallower than its limit at
/// age 0 of the old curve: 3,200 m.
const OLD_DEPTH_RANGE_M: f64 = 3_200.0;

/// The time constant of the old curve: 62.8 Myr.
const OLD_DEPTH_TIME_MYR: f64 = 62.8;

/// The ages (Myr) between which both curves of Parsons & Sclater (1977)
/// hold, 20 to 70 Myr: the young curve is below the old one at 20 Myr and
/// above it at 70 Myr.
const OVERLAP_MYR: (f64, f64) = (20.0, 70.0);

/// The most halvings of the bracket of the sea level ([`sea_level_m`]).
/// With heights within ±16 km the bracket is then below 2.6e-26 m wide, far
/// below any volume that matters; usually no double is left between its
/// ends long before.
const MAX_HALVINGS: u32 = 100;

/// The depth of the ocean floor (m) by its age, after Parsons & Sclater
/// (1977, "An analysis of the variation of ocean floor bathymetry and heat
/// flow with age", JGR 82(5), 803–827), with `t` the age in Myr:
///
/// - `d(t) = 2500 + 350 √t` (m) for `t < t*` (the young curve, which the
///   paper fits to ocean floor younger than 70 Myr);
/// - `d(t) = 6400 − 3200 e^(−t/62.8)` (m) for `t ≥ t*` (the old curve,
///   fitted to floor older than 20 Myr).
///
/// The two curves overlap from 20 to 70 Myr, and switching at 70 Myr would
/// leave a step of 78 m, so they are switched where they meet, at
/// `t* ≈ 26.27` Myr (depth about 4,294 m; design, 地殻から標高へ / 海洋底).
/// [`OceanDepth::new`] finds `t*` once, by bisection.
///
/// The age is converted to Myr only inside the formula; ages and depths are
/// in seconds and metres everywhere else.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OceanDepth {
    /// The age `t*` (Myr) at which the curves are switched.
    switch_myr: f64,
}

impl OceanDepth {
    /// The depth curve, with the age `t*` where the two curves meet found
    /// by bisection on [20, 70] Myr.
    ///
    /// The young curve minus the old one, `g(t) = 350 √t − 3900 +
    /// 3200 e^(−t/62.8)`, is negative at 20 Myr and positive at 70 Myr, and
    /// grows strictly in between (its slope is at least 0.3 m/Myr there), so
    /// it has one root. The bracket `[lo, hi]` with `g(lo) < 0 ≤ g(hi)` is
    /// halved until no double lies strictly between its ends (about 55
    /// halvings); `t*` is `hi`. At `t*` the two curves differ by about
    /// 1e-12 m.
    pub fn new() -> Self {
        let (mut lo, mut hi) = OVERLAP_MYR;
        loop {
            let mid = 0.5 * (lo + hi);
            if mid <= lo || mid >= hi {
                break;
            }
            if young_depth_m(mid) < old_depth_m(mid) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        Self { switch_myr: hi }
    }

    /// The age (s) at which the depth switches from the young curve to the
    /// old one: about 26.27 Myr (8.29 × 10¹⁴ s).
    pub fn switch_age_s(&self) -> f64 {
        self.switch_myr * SECONDS_PER_MYR
    }

    /// The depth (m, positive downwards) below the datum of ocean floor of
    /// age `age_s` (s, at least 0): the young curve below the switch age,
    /// the old one from it on. 2,500 m at the ridge, about 4,294 m at the
    /// switch, 6,218 m at 180 Myr; it grows strictly with the age and never
    /// reaches 6,400 m.
    pub fn depth_m(&self, age_s: f64) -> f64 {
        debug_assert!(age_s >= 0.0, "a negative age: {age_s} s");
        let t_myr = age_s / SECONDS_PER_MYR;
        if t_myr < self.switch_myr {
            young_depth_m(t_myr)
        } else {
            old_depth_m(t_myr)
        }
    }
}

impl Default for OceanDepth {
    /// [`OceanDepth::new`].
    fn default() -> Self {
        Self::new()
    }
}

/// The depth (m) of the young curve at the age `t_myr` (Myr).
fn young_depth_m(t_myr: f64) -> f64 {
    RIDGE_DEPTH_M + YOUNG_DEPTH_SLOPE_M * t_myr.sqrt()
}

/// The depth (m) of the old curve at the age `t_myr` (Myr).
fn old_depth_m(t_myr: f64) -> f64 {
    OLD_DEPTH_LIMIT_M - OLD_DEPTH_RANGE_M * det::exp(-t_myr / OLD_DEPTH_TIME_MYR)
}

/// The height (m) above the datum of continental crust `thickness_m` (m)
/// thick, by Airy isostasy:
///
/// `h = (T − T_ref)(ρ_m − ρ_c)/ρ_m + h_ref`,
///
/// with `T_ref` = [`CONTINENTAL_THICKNESS_M`] (35 km), `h_ref` =
/// [`REFERENCE_HEIGHT_M`] (300 m), `ρ_c` = [`CRUST_DENSITY_KG_M3`] and `ρ_m`
/// = [`MANTLE_DENSITY_KG_M3`]: each kilometre of crust raises the surface
/// by about 0.15 km. 15 km gives −2,730 m, 35 km +300 m, 70 km +5,603 m.
///
/// Crust floats on the mantle like ice on water (Airy 1855, "On the
/// computation of the effect of the attraction of mountain-masses, as
/// disturbing the apparent astronomical latitude of stations in geodetic
/// surveys", Phil. Trans. R. Soc. 145, 101–104): a column of crust
/// `T = T_ref + h' + b` thick, whose surface stands `h' = h − h_ref` above
/// that of the reference column, has a root that reaches `b` below the
/// reference column's base; at the depth of the root's base both columns
/// weigh the same,
/// `ρ_c g (h' + T_ref + b) = ρ_c g T_ref + ρ_m g b`, so `ρ_c h' = (ρ_m − ρ_c) b`
/// and `h' = (T − T_ref)(ρ_m − ρ_c)/ρ_m`. The surface gravity `g` cancels:
/// the height depends on the densities only, not on the planet's gravity
/// (design, 地殻から標高へ / 大陸).
pub fn continental_height_m(thickness_m: f64) -> f64 {
    (thickness_m - CONTINENTAL_THICKNESS_M) * (MANTLE_DENSITY_KG_M3 - CRUST_DENSITY_KG_M3)
        / MANTLE_DENSITY_KG_M3
        + REFERENCE_HEIGHT_M
}

/// The height (m) above the datum of the solid surface of every cell, from
/// its crust: `continental` (the crust type), `thickness_m` (m), `age_s`
/// (s) and `trench` of the Tectonics stage, one value per cell.
///
/// - Continental cells: [`continental_height_m`] of the thickness (their
///   age and trench flag are not used).
/// - Oceanic cells: minus the depth of their age ([`OceanDepth::depth_m`];
///   their thickness is not used: the depth curve is that of Earth's
///   oceanic crust, about 7 km thick).
/// - Oceanic trench cells: minus the larger of that depth and
///   [`TRENCH_DEPTH_M`], so −8 km.
///
/// # Panics
///
/// If the four slices are not of the same length.
pub fn surface_heights_m(
    continental: &[bool],
    thickness_m: &[f64],
    age_s: &[f64],
    trench: &[bool],
) -> Vec<f64> {
    let n = continental.len();
    assert!(
        thickness_m.len() == n && age_s.len() == n && trench.len() == n,
        "one crust type, thickness, age and trench flag per cell: {n}, {}, {}, {}",
        thickness_m.len(),
        age_s.len(),
        trench.len()
    );
    let depth = OceanDepth::new();
    (0..n)
        .map(|i| {
            if continental[i] {
                continental_height_m(thickness_m[i])
            } else {
                let d = depth.depth_m(age_s[i]);
                -if trench[i] && d < TRENCH_DEPTH_M {
                    TRENCH_DEPTH_M
                } else {
                    d
                }
            }
        })
        .collect()
}

/// The volume of water (m³) below the level `level_m` (m) over the heights
/// `heights_m` (m, one per cell of `grid`) on a planet of radius `radius_m`
/// (m): `Σ max(0, z − h_i) a_i R²`, with `a_i` the area of cell `i` (sr),
/// summed in cell order.
///
/// A bathtub model: every cell below the level holds water, whether it is
/// connected to the open ocean or not (enclosed basins, such as the seas
/// that continents enclose and continental rifts below the level, fill
/// too; lakes are milestone 5).
///
/// # Panics
///
/// If `heights_m` does not hold one value per cell of `grid`.
pub fn water_volume_m3(grid: &SphereGrid, heights_m: &[f64], radius_m: f64, level_m: f64) -> f64 {
    assert_eq!(
        heights_m.len(),
        grid.len(),
        "a height for each of the {} cells",
        grid.len()
    );
    let mut sum = 0.0;
    for (cell, &h) in (0..).zip(heights_m) {
        if h < level_m {
            sum += (level_m - h) * grid.area_sr(cell);
        }
    }
    sum * radius_m * radius_m
}

/// The sea level (m, on the datum of `heights_m`) at which the basins of
/// the surface `heights_m` (m, one per cell of `grid`) hold `volume_m3`
/// (m³, at least 0) of water on a planet of radius `radius_m` (m): the
/// level `z` with [`water_volume_m3`]`(z) = V` (design, 地殻から標高へ /
/// 海面).
///
/// With `h_min` and `h_max` the lowest and highest heights and
/// `C = water_volume_m3(h_max)` the water that the basins hold when filled
/// to the brim:
///
/// - `V = 0`, no ocean: `z = h_min`, the lowest height, where the first
///   water would collect. No cell lies below it, and the depth of water
///   `z − h_i` is 0 or less everywhere.
/// - `0 < V < C`: bisection on `[h_min, h_max]`, keeping
///   `water_volume_m3(lo) < V ≤ water_volume_m3(hi)`, halving until no
///   double lies strictly between the ends (at most 100 halvings); `z` is
///   `hi`, the lowest level found that holds at least `V`; it holds `V`
///   within 1e-9 relative for Earth-like volumes (ticket M2-05).
/// - `V ≥ C`, a water world: every cell lies under water, and the rest of
///   the water, `V − C`, covers the whole planet as a layer:
///   `z = h_max + (V − C) / (Σ a_i R²)`, in closed form. Finite for every
///   finite volume on a planet of radius 1 m or more, far below the
///   smallest radius that `PlanetParams::validate` accepts (1,500 km).
///
/// # Panics
///
/// If `heights_m` does not hold one value per cell of `grid`.
pub fn sea_level_m(grid: &SphereGrid, heights_m: &[f64], radius_m: f64, volume_m3: f64) -> f64 {
    assert_eq!(
        heights_m.len(),
        grid.len(),
        "a height for each of the {} cells",
        grid.len()
    );
    let (lowest, highest) = extremes(heights_m);
    if volume_m3 <= 0.0 {
        return lowest;
    }
    let capacity_m3 = water_volume_m3(grid, heights_m, radius_m, highest);
    if volume_m3 >= capacity_m3 {
        let area_sr: f64 = (0..grid.len() as u32).map(|cell| grid.area_sr(cell)).sum();
        return highest + (volume_m3 - capacity_m3) / (area_sr * radius_m * radius_m);
    }
    let (mut lo, mut hi) = (lowest, highest);
    for _ in 0..MAX_HALVINGS {
        let mid = 0.5 * (lo + hi);
        if mid <= lo || mid >= hi {
            break;
        }
        if water_volume_m3(grid, heights_m, radius_m, mid) < volume_m3 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    hi
}

/// The lowest and highest of `heights_m`, by comparisons in cell order (the
/// first of equal values is kept).
///
/// # Panics
///
/// If `heights_m` is empty.
fn extremes(heights_m: &[f64]) -> (f64, f64) {
    let (&first, rest) = heights_m.split_first().expect("at least one height");
    rest.iter().fold((first, first), |(lo, hi), &h| {
        (if h < lo { h } else { lo }, if h > hi { h } else { hi })
    })
}

/// Runs the Crust stage on `world` (a [`StageFn`](crate::pipeline::StageFn)):
/// reads the crust of the Tectonics stage ([`CRUST_TYPE`], [`CRUST_THICKNESS`],
/// [`CRUST_AGE`] and [`TRENCH`]), adds the field [`SURFACE_HEIGHT`] (m,
/// [`surface_heights_m`]) and sets `world.sea_level_m` to the level that
/// holds `world.params.ocean_volume_m3` on a planet of radius
/// `world.params.radius_m` ([`sea_level_m`]). It reports no progress of its
/// own. The surface gravity is not used (it cancels in Airy isostasy,
/// [`continental_height_m`]).
///
/// # Errors
///
/// [`SimError::Field`] if the world already has the field
/// [`SURFACE_HEIGHT`]; the world is then unchanged.
///
/// # Panics
///
/// If the world lacks one of the four fields of the Tectonics stage, which
/// the runner always runs first.
pub fn run(world: &mut World, _progress: &mut dyn FnMut(f64)) -> Result<(), SimError> {
    let heights = {
        let continental: Vec<bool> = input(world, CRUST_TYPE).iter().map(|&c| c != 0.0).collect();
        let trench: Vec<bool> = input(world, TRENCH).iter().map(|&t| t != 0.0).collect();
        surface_heights_m(
            &continental,
            input(world, CRUST_THICKNESS),
            input(world, CRUST_AGE),
            &trench,
        )
    };
    let sea_level = sea_level_m(
        &world.grid,
        &heights,
        world.params.radius_m,
        world.params.ocean_volume_m3,
    );
    world.add_field(SURFACE_HEIGHT, "m", heights)?;
    world.sea_level_m = sea_level;
    Ok(())
}

/// The values of the field `name` of `world`, an input of the stage.
///
/// # Panics
///
/// If the world has no such field.
fn input<'a>(world: &'a World, name: &str) -> &'a [f64] {
    match world.field(name) {
        Some(field) => &field.values,
        None => panic!("the Crust stage needs the field `{name}` of the Tectonics stage"),
    }
}

#[cfg(test)]
mod tests;

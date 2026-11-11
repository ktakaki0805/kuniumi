//! Plate tectonics (design: "プレートテクトニクス"; milestones 2 and 6).
//!
//! Milestone 2: one-shot plates, no time evolution. The surface is divided
//! into [`Settings::plate_count`] plates by competitive growth from seed
//! cells ([`partition`]), and each plate turns as a rigid cap about its own
//! Euler pole, at a speed chosen so that its fastest cell moves at 1 to
//! 10 cm/yr ([`euler_poles`]), the range of Earth's plates (design document,
//! プレートテクトニクス / プレート速度と分裂). Continental blobs grow inside
//! plates until they cover [`PlanetParams::continental_crust_fraction`] of
//! the surface ([`continents`]; about 40 % on Earth, design document,
//! プレートテクトニクス / 初期条件の目安), and every boundary cell is classified
//! as convergent, divergent or transform from the relative velocity of the
//! two plates ([`boundaries`]; the boundary types of 境界ごとの処理). Each
//! cell then gets its crust from its distance to the boundaries and coasts
//! of its own plate ([`crust`]): ocean floor ages away from the divergent
//! boundaries at their half spreading rates ([`ridge_times`]), with Earth's
//! distribution of area by age ([`earth_area_age`]); continental crust
//! thins towards the open ocean and continental rifts, where much of it
//! lies under the sea ([`margin_distances`]), thickens into plateaus at
//! collisions and into arcs inland of subduction zones, where it also
//! rises; the oceanic side of each subduction boundary is a trench.
//!
//! Milestone 6: kinematic simulation after Cortial et al. (2019). Rigid plates
//! of crust samples rotate about Euler poles in 1–2 Myr steps for 300–500 Myr,
//! starting from a supercontinent. Handles ocean–continent and ocean–ocean
//! subduction, continental collision and suturing, ridge spreading, rifting
//! along old sutures, and fixed hotspots. Samples are resampled onto the grid
//! about every 10 steps. Records the uplift rate of the final 20–50 Myr for
//! the erosion stage.
//!
//! Random streams: [`PARTITION_STREAM`] (seed cells, growth rates and the
//! slowness of each cell), [`POLES_STREAM`] (axes and speeds) and
//! [`CONTINENTS_STREAM`] (the slowness of each cell, the order of the
//! plates, the seed and rate of each blob), all under the stage's name
//! `"tectonics"` ([`StageId::name`](crate::StageId::name)).
//!
//! Outputs so far ([`run`]): [`PLATE_ID`] (class), [`PLATE_VELOCITY_EAST`]
//! and [`PLATE_VELOCITY_NORTH`] (m/s), [`CRUST_TYPE`] (class),
//! [`BOUNDARY_TYPE`] (class), [`CONVERGENCE_RATE`] (m/s),
//! [`CRUST_THICKNESS`] (m), [`CRUST_AGE`] (s), [`UPLIFT_RATE`] (m/s) and
//! [`TRENCH`] (class). Planned (M6): `orogeny_age` (s).
//!
//! [`Settings::plate_count`]: kuniumi_world::Settings::plate_count
//! [`PlanetParams::continental_crust_fraction`]: kuniumi_world::PlanetParams::continental_crust_fraction

use core::cmp::{Ordering, Reverse};
use std::collections::BinaryHeap;

use kuniumi_math::units::SECONDS_PER_YEAR;
use kuniumi_math::{Rng, Vec3, det};
use kuniumi_sphere::{SphereGrid, coords};
use kuniumi_world::World;

use crate::SimError;

mod boundaries;
mod continents;
mod crust;

pub use boundaries::{
    BOUNDARY_TYPE, Boundaries, BoundaryType, CONVERGENCE_RATE, NORMAL_WINDOW_RAD, boundaries,
};
pub use continents::{CONTINENT_PLATE_COVER, CONTINENTS_STREAM, CRUST_TYPE, continents};
pub use crust::{
    ARC_DISTANCE_M, ARC_HALF_WIDTH_M, ARC_THICKNESS_M, ARC_UPLIFT_RATIO, COLLISION_PLATEAU_M,
    COLLISION_REFERENCE_RATE_M_S, COLLISION_TAPER_M, COLLISION_THICKNESS_M, COLLISION_UPLIFT_RATIO,
    CONTINENTAL_THICKNESS_M, CRUST_AGE, CRUST_THICKNESS, Crust, MARGIN_COAST_THICKNESS_M,
    MARGIN_MIN_OCEAN_FRACTION, MARGIN_WIDTH_M, MAX_CRUST_AGE_S, MAX_UPLIFT_RATE_M_S, Nearest,
    OCEANIC_THICKNESS_M, RidgeTimes, TRENCH, UPLIFT_RATE, crust, earth_area_age, margin_distances,
    nearest_sources, ridge_times,
};

/// The random stream of the partition: the seed cells, then the growth rate
/// of each plate, then the slowness of each cell ([`partition`]).
pub const PARTITION_STREAM: &str = "tectonics/partition";

/// The random stream of the Euler poles: per plate, its fastest speed and
/// its axis ([`euler_poles`]).
pub const POLES_STREAM: &str = "tectonics/poles";

/// The class field of plate ids, `0..plate_count`.
pub const PLATE_ID: &str = "plate_id";

/// The field of the eastward component of each cell's plate velocity (m/s).
pub const PLATE_VELOCITY_EAST: &str = "plate_velocity_east";

/// The field of the northward component of each cell's plate velocity (m/s).
pub const PLATE_VELOCITY_NORTH: &str = "plate_velocity_north";

/// The slowest speed of the fastest cell of a plate: 1 cm/yr
/// (3.17 × 10⁻¹⁰ m/s), the low end of Earth's plate speeds (design
/// document, プレートテクトニクス / プレート速度と分裂).
pub const MIN_PLATE_SPEED_M_S: f64 = 0.01 / SECONDS_PER_YEAR;

/// The fastest speed of any cell: 10 cm/yr (3.17 × 10⁻⁹ m/s), the high end
/// of Earth's plate speeds.
pub const MAX_PLATE_SPEED_M_S: f64 = 0.10 / SECONDS_PER_YEAR;

/// The least distance between two seed cells, as a fraction of
/// `√(4π / plate_count)` rad, the side of a square of the mean plate area:
/// 0.4. Chosen with the owner in ticket M2-02, so that no plate starts
/// squeezed between others.
pub const SEED_SPACING: f64 = 0.4;

/// The fastest growth rate of a plate, relative to the slowest: rates are
/// uniform in `[1, 3)` (ticket M2-02). A plate's area grows roughly with the
/// square of its rate, which spreads the plate sizes.
pub const MAX_GROWTH_RATE: f64 = 3.0;

/// The spread of the slowness of the cells: `exp(3 (u − ½))` for `u`
/// uniform in `[0, 1)`, from `e^−1.5` to `e^1.5`, a factor of 20 (ticket
/// M2-02). It makes the plate boundaries irregular.
pub const SLOWNESS_SPREAD: f64 = 3.0;

/// Rejections of seed candidates in a row after which a seed no longer
/// needs [`SEED_SPACING`]; a fallback for grids too small for the spacing,
/// so that choosing seeds always ends.
const MAX_REJECTIONS: u32 = 10_000;

/// The least lever of a plate (the largest `|axis × p|` over its cells)
/// with which an axis is accepted: an axis that passes within about
/// 1e-6 rad of every cell of the plate (possible only for a plate of one
/// cell) is drawn again, so that the angular speed stays finite.
const MIN_LEVER: f64 = 1e-6;

/// The rotation of a rigid plate about an axis through the planet's centre.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EulerPole {
    /// The axis, a unit vector in the planet-centred frame. The plate turns
    /// counter-clockwise as seen from the tip of the axis looking towards
    /// the centre (right-hand rule, as [`Vec3::rotate_about`]).
    pub axis: Vec3,
    /// The angular speed (rad/s), positive.
    pub angular_speed_rad_s: f64,
}

impl EulerPole {
    /// The velocity (m/s) of the surface point at the unit vector
    /// `position` on a planet of radius `radius_m`: `ω × r` with
    /// `ω = angular_speed_rad_s · axis` and `r = radius_m · position`. It is
    /// tangent to the sphere, of length `ω R sin θ` for the angle `θ`
    /// between the axis and the point.
    pub fn velocity_m_s(&self, position: Vec3, radius_m: f64) -> Vec3 {
        (self.axis * self.angular_speed_rad_s).cross(position * radius_m)
    }
}

/// The heap entry of a growth: a cell that `source` (a plate in
/// [`partition`], a continental blob in [`continents`]) can reach at cost
/// `priority`. Ordered by the priority ([`f64::total_cmp`]), then the cell,
/// then the source, so that no two different entries compare equal and the
/// heap pops entries in one order whatever its internals.
#[derive(Clone, Copy, Debug)]
struct Entry {
    priority: f64,
    cell: u32,
    source: u32,
}

impl Ord for Entry {
    fn cmp(&self, other: &Self) -> Ordering {
        self.priority
            .total_cmp(&other.priority)
            .then(self.cell.cmp(&other.cell))
            .then(self.source.cmp(&other.source))
    }
}

impl PartialOrd for Entry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for Entry {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for Entry {}

/// Divides the cells of `grid` into `plate_count` plates; returns the plate
/// of each cell, in `0..plate_count`.
///
/// Every cell gets one plate and every plate is connected (each cell is
/// given the plate of the neighbour that reached it). Every plate has at
/// least its seed cell. The method (ticket M2-02, chosen with the owner),
/// with every random number from [`PARTITION_STREAM`] in this order:
///
/// 1. Seed cells: uniformly random cells ([`Rng::below`]), each at least
///    [`SEED_SPACING`]` · √(4π / K)` rad (great circle) from the seeds
///    before it; plate `k` grows from the `k`-th. After 10,000 rejected
///    candidates in a row the next free cell is taken whatever its
///    distance, so that the choice ends even on grids too small for the
///    spacing.
/// 2. A growth rate per plate, in plate order: `r_k = 1 + 2u`, uniform in
///    `[1, `[`MAX_GROWTH_RATE`]`)`.
/// 3. A slowness per cell, in cell order: `w_i = exp(3 (u − ½))`
///    ([`SLOWNESS_SPREAD`]).
/// 4. Growth from every seed at once, as Dijkstra's algorithm with several
///    sources: each seed enters a priority queue at priority 0; the entry of
///    least priority is taken; if its cell has no plate yet, the cell gets
///    the entry's plate `k` and every neighbour `j` without a plate enters
///    the queue at priority `p + d_ij (w_i + w_j) / 2 / r_k`, with `d_ij`
///    the distance between the centres (rad). Entries are ordered by
///    priority ([`f64::total_cmp`]), then cell, then plate, so no two
///    differ only in their order and the pops are the same whatever the
///    heap's internals.
///
/// The slowness of the cells roughens the boundaries; the rates bend them
/// into arcs around the slower plates and spread the plate sizes, a few
/// large and some small, as on Earth.
///
/// # Panics
///
/// If `plate_count` is 0 or larger than the number of cells.
pub fn partition(grid: &SphereGrid, plate_count: u32, seed: u64) -> Vec<u32> {
    let n = grid.len();
    assert!(
        plate_count >= 1 && plate_count as usize <= n,
        "cannot divide {n} cells into {plate_count} plates"
    );
    let mut rng = Rng::for_stream(seed, PARTITION_STREAM);
    let seeds = seed_cells(grid, plate_count, &mut rng);
    let rates: Vec<f64> = (0..plate_count)
        .map(|_| 1.0 + (MAX_GROWTH_RATE - 1.0) * rng.next_f64())
        .collect();
    let slowness: Vec<f64> = (0..n)
        .map(|_| det::exp(SLOWNESS_SPREAD * (rng.next_f64() - 0.5)))
        .collect();

    const NONE: u32 = u32::MAX;
    let mut plates = vec![NONE; n];
    let mut heap = BinaryHeap::with_capacity(2 * n);
    for (plate, &cell) in (0..).zip(&seeds) {
        heap.push(Reverse(Entry {
            priority: 0.0,
            cell,
            source: plate,
        }));
    }
    while let Some(Reverse(Entry {
        priority,
        cell,
        source: plate,
    })) = heap.pop()
    {
        if plates[cell as usize] != NONE {
            continue;
        }
        plates[cell as usize] = plate;
        let (w_i, rate) = (slowness[cell as usize], rates[plate as usize]);
        let distances_rad = grid.centre_distances_rad(cell);
        for (&j, &d_rad) in grid.neighbors(cell).iter().zip(distances_rad) {
            if plates[j as usize] == NONE {
                let step = d_rad * (w_i + slowness[j as usize]) * 0.5 / rate;
                heap.push(Reverse(Entry {
                    priority: priority + step,
                    cell: j,
                    source: plate,
                }));
            }
        }
    }
    debug_assert!(
        plates.iter().all(|&p| p < plate_count),
        "a cell was not reached; the grid is not connected"
    );
    plates
}

/// The seed cells of `plate_count` plates: cells drawn with
/// [`Rng::below`], each accepted if it is not a seed yet and lies at least
/// [`SEED_SPACING`]` · √(4π / plate_count)` rad (great circle) from every
/// seed accepted before it. After [`MAX_REJECTIONS`] rejections in a row the
/// next free cell is accepted whatever its distance, so the loop ends even
/// on grids too small for the spacing.
fn seed_cells(grid: &SphereGrid, plate_count: u32, rng: &mut Rng) -> Vec<u32> {
    let n = grid.len() as u64;
    let min_rad = SEED_SPACING * (4.0 * core::f64::consts::PI / f64::from(plate_count)).sqrt();
    let mut seeds: Vec<u32> = Vec::with_capacity(plate_count as usize);
    let mut rejections = 0;
    while seeds.len() < plate_count as usize {
        let cell = rng.below(n) as u32;
        let p = grid.position(cell);
        let crowded = rejections < MAX_REJECTIONS
            && seeds
                .iter()
                .any(|&s| coords::great_circle_rad(grid.position(s), p) < min_rad);
        if seeds.contains(&cell) || crowded {
            rejections += 1;
        } else {
            rejections = 0;
            seeds.push(cell);
        }
    }
    seeds
}

/// The Euler pole of each plate of the partition `plates` (the plate of
/// each cell of `grid`, in `0..plate_count`, every plate with a cell) on a
/// planet of radius `radius_m`.
///
/// For each plate in order, from [`POLES_STREAM`]:
///
/// 1. The speed of its fastest cell, log-uniform over
///    [[`MIN_PLATE_SPEED_M_S`], [`MAX_PLATE_SPEED_M_S`]]:
///    `v = v_min exp(u ln(v_max / v_min))`, clamped to the range against
///    rounding. The range is the design document's; the log-uniform
///    distribution is a modelling choice made with the owner (ticket M2-02),
///    which puts half the plates below √10 ≈ 3.2 cm/yr.
/// 2. The axis, uniform on the sphere (two numbers: `z = 2u − 1`,
///    `λ = 2πv`). Its lever `s = max |axis × p_i|` over the plate's cells
///    is the sine of the largest angle between the axis and a cell; while
///    `s` is below 1e-6 (an axis through a plate of one cell) the axis is
///    drawn again, so that `ω` stays finite.
/// 3. The angular speed `ω = v / (R s)`, so that the fastest cell moves at
///    `v` (up to rounding) and no cell faster.
///
/// # Panics
///
/// If `plates` does not hold one plate per cell of `grid`, a value of
/// `plates` is not below `plate_count`, or a plate has no cell.
pub fn euler_poles(
    grid: &SphereGrid,
    plates: &[u32],
    plate_count: u32,
    radius_m: f64,
    seed: u64,
) -> Vec<EulerPole> {
    assert_eq!(
        plates.len(),
        grid.len(),
        "a plate for each of the {} cells, got {}",
        grid.len(),
        plates.len()
    );
    let mut members = vec![Vec::new(); plate_count as usize];
    for (cell, &plate) in (0..).zip(plates) {
        members[plate as usize].push(grid.position(cell));
    }
    let mut rng = Rng::for_stream(seed, POLES_STREAM);
    let ln_ratio = det::ln(MAX_PLATE_SPEED_M_S / MIN_PLATE_SPEED_M_S);
    (0..)
        .zip(&members)
        .map(|(plate, positions): (u32, &Vec<Vec3>)| {
            assert!(!positions.is_empty(), "plate {plate} has no cell");
            let speed_m_s = (MIN_PLATE_SPEED_M_S * det::exp(ln_ratio * rng.next_f64()))
                .clamp(MIN_PLATE_SPEED_M_S, MAX_PLATE_SPEED_M_S);
            loop {
                let axis = random_axis(&mut rng);
                // The levers are non-negative, so the sign of a zero cannot
                // matter here; a zero lever is rejected.
                let lever = positions
                    .iter()
                    .map(|&p| axis.cross(p).length())
                    .fold(0.0, f64::max);
                if lever >= MIN_LEVER {
                    break EulerPole {
                        axis,
                        angular_speed_rad_s: speed_m_s / (radius_m * lever),
                    };
                }
            }
        })
        .collect()
}

/// A unit vector uniform on the sphere, from two random numbers: by
/// Archimedes' hat-box theorem `z = 2u − 1` uniform in `[−1, 1)` and the
/// longitude `λ = 2πv` give `(√(1 − z²) cos λ, √(1 − z²) sin λ, z)`.
/// `|z| ≤ 1`, so `1 − z²` is never negative; the length is 1 within a few
/// units of 1e-16.
fn random_axis(rng: &mut Rng) -> Vec3 {
    let z = 2.0 * rng.next_f64() - 1.0;
    let lon_rad = 2.0 * core::f64::consts::PI * rng.next_f64();
    let (sin, cos) = det::sin_cos(lon_rad);
    let horizontal = (1.0 - z * z).sqrt();
    Vec3::new(horizontal * cos, horizontal * sin, z)
}

/// The eastward and northward components (m/s) of the velocity of the
/// point at the unit vector `position` turning with `pole` on a planet of
/// radius `radius_m`: [`EulerPole::velocity_m_s`] projected on
/// [`coords::east`] and [`coords::north`]. The velocity is tangent, so the
/// two components hold all of it.
fn east_north_m_s(pole: &EulerPole, position: Vec3, radius_m: f64) -> (f64, f64) {
    let v = pole.velocity_m_s(position, radius_m);
    (
        v.dot(coords::east(position)),
        v.dot(coords::north(position)),
    )
}

/// The eastward and northward plate velocity (m/s) of every cell of
/// `grid`, given the plate of each cell and the poles of the plates.
fn velocities_m_s(
    grid: &SphereGrid,
    plates: &[u32],
    poles: &[EulerPole],
    radius_m: f64,
) -> (Vec<f64>, Vec<f64>) {
    (0..)
        .zip(plates)
        .map(|(cell, &plate)| east_north_m_s(&poles[plate as usize], grid.position(cell), radius_m))
        .unzip()
}

/// Runs the Tectonics stage on `world` (a [`StageFn`](crate::pipeline::StageFn)):
/// the plates of `world.settings.plate_count` ([`partition`]), their poles
/// for `world.params.radius_m` ([`euler_poles`]), the continents of
/// `world.params.continental_crust_fraction` ([`continents`]), the
/// boundaries ([`boundaries`]) and the crust ([`crust`]), all from
/// `world.seed`, and adds the fields [`PLATE_ID`], [`PLATE_VELOCITY_EAST`],
/// [`PLATE_VELOCITY_NORTH`], [`CRUST_TYPE`] (0 oceanic, 1 continental),
/// [`BOUNDARY_TYPE`] ([`BoundaryType::class`]), [`CONVERGENCE_RATE`],
/// [`CRUST_THICKNESS`], [`CRUST_AGE`], [`UPLIFT_RATE`] and [`TRENCH`] (0, or
/// 1 for a trench), in that order. It reports no progress of its own.
///
/// The world's settings and parameters must be valid (`Settings::validate`:
/// at least 1,000 cells and 2 to 64 plates; `PlanetParams::validate`: a
/// continental fraction in `[0, 1]`), as the runner ensures.
///
/// # Errors
///
/// [`SimError::Field`] if the world already has one of the fields.
pub fn run(world: &mut World, _progress: &mut dyn FnMut(f64)) -> Result<(), SimError> {
    let plate_count = world.settings.plate_count;
    let radius_m = world.params.radius_m;
    let plates = partition(&world.grid, plate_count, world.seed);
    let poles = euler_poles(&world.grid, &plates, plate_count, radius_m, world.seed);
    let (east, north) = velocities_m_s(&world.grid, &plates, &poles, radius_m);
    let continental = continents(
        &world.grid,
        &plates,
        plate_count,
        world.params.continental_crust_fraction,
        world.seed,
    );
    let bounds = boundaries(&world.grid, &plates, &poles, radius_m);
    let crust = crust(&world.grid, &plates, &continental, &bounds, radius_m);
    world.add_class_field(PLATE_ID, plates.iter().map(|&p| f64::from(p)).collect())?;
    world.add_field(PLATE_VELOCITY_EAST, "m/s", east)?;
    world.add_field(PLATE_VELOCITY_NORTH, "m/s", north)?;
    world.add_class_field(
        CRUST_TYPE,
        continental
            .iter()
            .map(|&c| f64::from(u8::from(c)))
            .collect(),
    )?;
    world.add_class_field(
        BOUNDARY_TYPE,
        bounds.kind.iter().map(|k| k.class()).collect(),
    )?;
    world.add_field(CONVERGENCE_RATE, "m/s", bounds.convergence_rate_m_s)?;
    world.add_field(CRUST_THICKNESS, "m", crust.thickness_m)?;
    world.add_field(CRUST_AGE, "s", crust.age_s)?;
    world.add_field(UPLIFT_RATE, "m/s", crust.uplift_rate_m_s)?;
    world.add_class_field(
        TRENCH,
        crust
            .trench
            .iter()
            .map(|&t| f64::from(u8::from(t)))
            .collect(),
    )?;
    Ok(())
}

#[cfg(test)]
mod tests;

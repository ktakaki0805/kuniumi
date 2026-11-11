//! Plate boundaries for the one-shot tectonics (ticket M2-03): the type of
//! each boundary cell (convergent, divergent or transform) and the rate at
//! which the plates approach there.

use kuniumi_math::Vec3;
use kuniumi_sphere::{SphereGrid, coords};

use super::EulerPole;

/// The class field of the boundary type of each cell
/// ([`BoundaryType::class`]): 0 none, 1 convergent, 2 divergent, 3
/// transform.
pub const BOUNDARY_TYPE: &str = "boundary_type";

/// The field of the convergence rate of each cell (m/s): the normal part of
/// the relative plate velocity across its boundary, positive when the
/// plates approach, negative when they separate, 0 away from boundaries
/// ([`Boundaries::convergence_rate_m_s`]).
pub const CONVERGENCE_RATE: &str = "convergence_rate";

/// The radius of the window over which the normal of a boundary is taken
/// (rad): 0.08, about 510 km on Earth and two to three cells at 10,000
/// cells (spacing 0.035 rad). Chosen with the owner in ticket M2-03: a
/// single Voronoi edge of a staircase boundary is often oblique to the
/// boundary as a whole, and the window averages the staircase out. In
/// radians, so that the types do not depend on the planet's radius.
///
/// Below about 5,000 cells (spacing above 0.05 rad) the window holds
/// little more than the cell itself, so the normal is the cell's own and
/// the types are noisier: on two-plate tests at 1,000 cells up to 2.5 % of
/// the boundary cells are misclassified, against none from 10,000 cells
/// up (measured in the review of ticket M2-03).
pub const NORMAL_WINDOW_RAD: f64 = 0.08;

/// The type of a cell's plate boundary.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u8)]
pub enum BoundaryType {
    /// No neighbour on another plate.
    None = 0,
    /// The plates approach: the normal part of their relative velocity is
    /// larger than the tangential part.
    Convergent = 1,
    /// The plates separate: the normal part is larger than the tangential
    /// part.
    Divergent = 2,
    /// The plates slide past each other: the tangential part is at least
    /// as large as the normal part.
    Transform = 3,
}

impl BoundaryType {
    /// The class id of the type in [`BOUNDARY_TYPE`].
    pub fn class(self) -> f64 {
        f64::from(self as u8)
    }
}

/// The plate boundaries of a partition ([`boundaries`]), one value per
/// cell.
#[derive(Clone, Debug, PartialEq)]
pub struct Boundaries {
    /// The type of the cell's boundary.
    pub kind: Vec<BoundaryType>,
    /// The convergence rate `c` (m/s) across the cell's boundary: positive
    /// when the plates approach, negative when they separate (the spreading
    /// rate is `−c`), and `+0` for [`BoundaryType::None`]. Transform cells
    /// keep their small normal part, `|c|` at most the tangential one.
    pub convergence_rate_m_s: Vec<f64>,
    /// The plate across the cell's boundary, or `None` for
    /// [`BoundaryType::None`].
    pub other_plate: Vec<Option<u32>>,
}

/// The raw normals of the boundary cells: for each cell, the other plates
/// among its neighbours in increasing order, each with
/// `g = Σ ℓ_ij u_ij` over its neighbours `j` on that plate.
struct RawNormals {
    /// CSR offsets into `entries`, one range per cell.
    offsets: Vec<usize>,
    /// `(other plate, g)`, sorted by plate within each cell.
    entries: Vec<(u32, Vec3)>,
}

impl RawNormals {
    fn new(grid: &SphereGrid, plates: &[u32]) -> Self {
        let mut offsets = Vec::with_capacity(grid.len() + 1);
        let mut entries: Vec<(u32, Vec3)> = Vec::new();
        offsets.push(0);
        for i in 0..grid.len() as u32 {
            let start = entries.len();
            let (a, p_i) = (plates[i as usize], grid.position(i));
            let lengths_rad = grid.edge_lengths_rad(i);
            for (&j, &length_rad) in grid.neighbors(i).iter().zip(lengths_rad) {
                let b = plates[j as usize];
                if b == a {
                    continue;
                }
                let u = unit_chord(p_i, grid.position(j));
                let cell_entries = &mut entries[start..];
                match cell_entries.iter_mut().find(|(plate, _)| *plate == b) {
                    Some((_, g)) => *g += u * length_rad,
                    None => entries.push((b, u * length_rad)),
                }
            }
            entries[start..].sort_unstable_by_key(|&(plate, _)| plate);
            offsets.push(entries.len());
        }
        Self { offsets, entries }
    }

    /// The other plates of cell `i` with their raw normals.
    fn of(&self, i: u32) -> &[(u32, Vec3)] {
        &self.entries[self.offsets[i as usize]..self.offsets[i as usize + 1]]
    }

    /// The raw normal of cell `k` for the pair `(a, b)`, oriented from `a`
    /// to `b`, or `None` if `k` is not in the band of the pair (not on `a`
    /// next to `b`, nor on `b` next to `a`).
    fn oriented(&self, plates: &[u32], k: u32, a: u32, b: u32) -> Option<Vec3> {
        let find = |other: u32| {
            self.of(k)
                .iter()
                .find(|&&(plate, _)| plate == other)
                .map(|&(_, g)| g)
        };
        match plates[k as usize] {
            plate if plate == a => find(b),
            plate if plate == b => find(a).map(|g| -g),
            _ => None,
        }
    }
}

/// The unit vector from `p` to `q` (distinct unit vectors): the direction
/// of the great circle from `p` to `q` at the midpoint of the arc, which is
/// the normal of the Voronoi edge between two neighbouring cells.
fn unit_chord(p: Vec3, q: Vec3) -> Vec3 {
    let chord = q - p;
    chord / chord.length()
}

/// The unit tangent at `p` along the projection of `v` on the tangent
/// plane, or `None` if that projection is shorter than 1e-9 of `|v|` (or
/// `v` is zero).
fn tangent_direction(p: Vec3, v: Vec3) -> Option<Vec3> {
    let tangent = v - p * v.dot(p);
    let length = tangent.length();
    (length > 1e-9 * v.length()).then(|| tangent / length)
}

/// The angular velocity vector `ω · axis` of a pole (rad/s).
fn angular_velocity(pole: &EulerPole) -> Vec3 {
    pole.axis * pole.angular_speed_rad_s
}

/// The boundary type, convergence rate and plate across of each cell of
/// `grid`, given the plate of each cell (`plates`, every plate with a pole
/// in `poles`) on a planet of radius `radius_m`.
///
/// A cell `i` of plate `a` is a boundary cell of the pair `(a, b)` if one
/// of its neighbours is on plate `b`; the band of `(a, b)` is the boundary
/// cells of that pair on both sides. The method (ticket M2-03, chosen with
/// the owner):
///
/// 1. Raw normal: `g_ib = Σ ℓ_ij u_ij` over the neighbours `j` on plate
///    `b`, `u_ij = (p_j − p_i)/|p_j − p_i|` the direction of the great
///    circle from `i` to `j` at the edge's midpoint (the normal of the
///    Voronoi edge) and `ℓ_ij` the edge length (rad). Up to the factor
///    `−2 a_i` it is the Green–Gauss gradient of plate `b`'s indicator on
///    the cell; it points from `a` to `b`.
/// 2. Window: `N = Σ g_k` over the cells `k` of the band reached from `i`
///    by a breadth-first search through neighbours in the band whose
///    centres lie within [`NORMAL_WINDOW_RAD`] of `p_i` (great circle),
///    each oriented from `a` to `b` (`g_kb` on `a`, `−g_ka` on `b`), summed
///    in the order of the search. The unit normal `n̂` is `N` projected on
///    the tangent plane at `p_i`, normalised. If that projection is shorter
///    than 1e-9 of `|N|` (a cell all but surrounded by `b`, whose normal is
///    undefined, such as a plate of one cell) the tangent direction from
///    `p_i` towards its first neighbour on `b` (in the order of
///    [`SphereGrid::neighbors`]) is used instead.
/// 3. The velocity of `b` relative to `a` at the cell,
///    `Δv = (ω_b − ω_a) × (R p_i)` (m/s, `ω = ω axis` of each pole); the
///    convergence rate `c = −Δv · n̂`, positive when the plates approach,
///    and the shear `s = Δv · (p_i × n̂)`.
/// 4. The type: convergent if `c > |s|`, divergent if `−c > |s|`, transform
///    otherwise. The threshold is 45° between the relative velocity and the
///    boundary normal, where the normal part stops dominating; no relative
///    motion at all counts as transform.
/// 5. A cell next to several other plates takes the pair with the largest
///    relative speed `|Δv|` ("the strongest of its edges"); between equal
///    speeds, the lower plate id.
///
/// Cells with no neighbour on another plate get [`BoundaryType::None`], a
/// rate of `+0` and no plate across.
///
/// Accuracy: on two plates split by great circles, the rate lies within
/// about 10 % of the relative speed `|Δv|` (mostly 3 to 6 %) of the exact
/// normal rate at 10,000 cells, and within 1 % at 100,000 cells; see
/// [`NORMAL_WINDOW_RAD`] for coarser grids. The poles scale with `1/R`
/// ([`euler_poles`](super::euler_poles)) and the window is in radians, so
/// the types and rates of a partition do not depend on the radius (up to
/// rounding).
///
/// # Panics
///
/// If `plates` does not hold one plate per cell of `grid`, or a plate has
/// no pole in `poles`.
pub fn boundaries(
    grid: &SphereGrid,
    plates: &[u32],
    poles: &[EulerPole],
    radius_m: f64,
) -> Boundaries {
    let n = grid.len();
    assert_eq!(
        plates.len(),
        n,
        "a plate for each of the {n} cells, got {}",
        plates.len()
    );
    let raw = RawNormals::new(grid, plates);
    let mut result = Boundaries {
        kind: vec![BoundaryType::None; n],
        convergence_rate_m_s: vec![0.0; n],
        other_plate: vec![None; n],
    };
    let mut visited = vec![false; n];
    let mut window: Vec<u32> = Vec::new();
    for i in 0..n as u32 {
        let (a, p) = (plates[i as usize], grid.position(i));
        let omega_a = angular_velocity(&poles[a as usize]);
        // (|Δv|, type, c, b) of the strongest pair so far.
        let mut strongest: Option<(f64, BoundaryType, f64, u32)> = None;
        for &(b, _) in raw.of(i) {
            // Breadth-first search through the band of (a, b) within the
            // window; `window` holds the cells reached, in search order.
            window.clear();
            window.push(i);
            visited[i as usize] = true;
            let mut sum = Vec3::ZERO;
            let mut next = 0;
            while let Some(&k) = window.get(next) {
                next += 1;
                sum += raw
                    .oriented(plates, k, a, b)
                    .expect("cells in the window are in the band");
                for &j in grid.neighbors(k) {
                    if !visited[j as usize]
                        && raw.oriented(plates, j, a, b).is_some()
                        && coords::great_circle_rad(p, grid.position(j)) <= NORMAL_WINDOW_RAD
                    {
                        visited[j as usize] = true;
                        window.push(j);
                    }
                }
            }
            for &k in &window {
                visited[k as usize] = false;
            }
            let normal = tangent_direction(p, sum).unwrap_or_else(|| {
                let first = grid
                    .neighbors(i)
                    .iter()
                    .find(|&&j| plates[j as usize] == b)
                    .expect("a boundary cell has a neighbour on the other plate");
                tangent_direction(p, unit_chord(p, grid.position(*first)))
                    .expect("neighbouring centres are neither equal nor opposite")
            });
            let dv = (angular_velocity(&poles[b as usize]) - omega_a).cross(p * radius_m);
            let c = -dv.dot(normal);
            let s = dv.dot(p.cross(normal));
            let kind = if c > s.abs() {
                BoundaryType::Convergent
            } else if -c > s.abs() {
                BoundaryType::Divergent
            } else {
                BoundaryType::Transform
            };
            let speed = dv.length();
            if strongest.is_none_or(|(best, ..)| speed > best) {
                strongest = Some((speed, kind, c, b));
            }
        }
        if let Some((_, kind, c, b)) = strongest {
            result.kind[i as usize] = kind;
            result.convergence_rate_m_s[i as usize] = c;
            result.other_plate[i as usize] = Some(b);
        }
    }
    result
}

#[cfg(test)]
mod tests;

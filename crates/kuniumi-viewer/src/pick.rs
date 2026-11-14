//! Picking: the cell under the cursor, found on the triangles of the globe
//! as they are drawn (ticket MV-04).
//!
//! # The drawn surface
//!
//! A vertex `(x, y, z, h)` of [`GlobeMesh`] is drawn at
//! `q = (x, y, z)·(1 + σ·h)`, with `σ` the frame's `relief_scale_per_m` (the
//! `f32` handed to the GPU, `shader.wgsl`); the picker computes the same in
//! `f64`. The fans of all cells tile the unit sphere once (the mesh's
//! shared vertices are bit-identical, its triangles counter-clockwise from
//! outside, their solid angles sum to 4π), and the radial projection of a
//! straight triangle onto the unit sphere is the spherical triangle of its
//! vertex directions. So while every drawn radius `1 + σ·h` is positive the
//! surface is a radial graph: every direction from the centre meets it once.
//!
//! # The search
//!
//! For a ray `x(t) = o + t·d` (`|d| = 1`) from outside the drawn surface:
//!
//! 1. The surface lies between the spheres `r_top = 1 + σ·h_max` and
//!    `r_bot = (1 + σ·h_min)·min_v cos ρ_v` (`h_min ≤ 0 ≤ h_max` the
//!    extreme vertex heights, `ρ_v` the radius of cell `v`'s cap, step 3;
//!    both widened by [`SHELL_MARGIN`] for the `f32` rounding of the vertex
//!    directions). No point of a triangle lies farther out than its
//!    farthest vertex, but a flat triangle sags below its vertices (review
//!    R1-1 of MV-04): a point `x` of cell `v`'s fan is a convex combination
//!    of drawn vertices `qᵢ` whose directions lie within `ρ_v` of `ĉ_v`, so
//!    `|x| ≥ x·ĉ_v ≥ minᵢ qᵢ·ĉ_v ≥ minᵢ |qᵢ| cos ρ_v` (for positive radii
//!    and `cos ρ_v > 0`). `t0` is the ray's entry into `r_top` (none: no
//!    hit), `t1` its entry into `r_bot`, or its exit from `r_top` if it
//!    misses `r_bot`. The first hit lies in `[t0, t1]`: at `t0` the ray is
//!    on or above the surface, at `t1` on or under it, so it has crossed
//!    the surface on the way.
//! 2. Over `[t0, t1]` the directions `x(t)/|x(t)|` sweep the great-circle
//!    arc from `a = x(t0)/|x(t0)|` to `b = x(t1)/|x(t1)|` in the plane of the
//!    ray and the centre, monotonically and through less than π (the segment
//!    stays at a radius of at least `r_bot > 0`).
//! 3. The fan of cell `v` lies in the spherical cap around its drawn centre
//!    direction `ĉ_v` of radius `ρ_v`, the largest angle to its drawn
//!    Voronoi vertex directions (a cap narrower than a hemisphere is convex).
//!    The candidates are found by a walk over neighbours, from the cell that
//!    [`SphereGrid::nearest`] gives for `a`, through the cells whose cap
//!    widened by [`CAP_MARGIN_RAD`] comes within reach of the arc (an
//!    endpoint inside it, or the point of the great circle nearest `ĉ_v`
//!    inside the arc and within the cap). The cells whose drawn fans the arc
//!    crosses form a connected chain (it passes from fan to fan across
//!    shared edges or vertices, between neighbours), and each of them passes
//!    the test, so the walk reaches every cell the ray can hit. The margin
//!    covers the difference between the exact Voronoi cell that `nearest`
//!    returns and the drawn fan (the `f32` rounding of the vertices, about
//!    10⁻⁷ rad) and the rounding of the arc's plane.
//! 4. Each candidate's triangles are tested with the Möller–Trumbore
//!    algorithm (T. Möller and B. Trumbore 1997, "Fast, minimum storage
//!    ray-triangle intersection", *Journal of Graphics Tools* 2(1), 21–28),
//!    front faces only (the GPU culls back faces; the first crossing of a
//!    closed surface from outside is a front face anyway), barycentric
//!    bounds inclusive, `t > 0`. The hit is the smallest `t`, ties to the
//!    lower cell index.
//! 5. Below `r_bot` = [`MIN_BOTTOM_RADIUS`], or where the bound of step 1
//!    does not hold (a vertex at a radius of 0 or less, a cap of a
//!    hemisphere or more), every triangle of the mesh is tested instead:
//!    slow, but still what is drawn. The viewer never draws a vertex below
//!    half the radius (`view::Relief::max_relief`, MV-05), so on its grids
//!    this applies only to an exaggeration that it does not offer.
//!
//! The cell is the one whose drawn triangle is hit, not the nearest centre,
//! so the pick agrees with the drawn fans up to the GPU's `f32` rounding.
//! Under the translucent sea the pick is the terrain seen through it.
//!
//! Nothing here reaches world data (D-19): the viewer only shows the cell.

use core::cmp::Ordering;

use kuniumi_math::{Vec3, det};
use kuniumi_sphere::SphereGrid;

use crate::mesh::GlobeMesh;

/// The margin by which each cell's cap is widened (rad): a hundred times the
/// `f32` rounding of the drawn vertex directions; 64 m on Earth, against a
/// cell spacing of some 30 km at 650,000 cells.
pub const CAP_MARGIN_RAD: f64 = 1.0e-5;

/// The relative margin by which the shell of the drawn surface is widened:
/// ten times the `f32` rounding of a vertex direction's length.
pub const SHELL_MARGIN: f64 = 1.0e-6;

/// The smallest lower bound `r_bot` of the drawn surface (planet radii,
/// module documentation, step 1) for which the walk is used; below it, every
/// triangle is tested (see [`Arc`] for the precision this keeps).
pub const MIN_BOTTOM_RADIUS: f64 = 0.01;

/// An arc whose endpoints are closer than this (the sine of their angle,
/// with a positive cosine) is tested by its endpoints alone: it lies within
/// 10⁻⁷ rad of them, well inside [`CAP_MARGIN_RAD`].
const SHORT_ARC_SIN: f64 = 1.0e-7;

/// Where a ray meets the drawn globe.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Hit {
    /// The cell whose drawn triangle the ray meets first.
    pub cell: u32,
    /// The distance from the ray's origin to the hit (planet radii).
    pub distance: f64,
}

/// What the picker derives from a mesh once, and the state of its walks.
#[derive(Clone, Debug)]
pub struct Picker {
    /// Per cell: the cosine of its cap's radius `ρ_v` widened by
    /// [`CAP_MARGIN_RAD`] (at most 0 for a cap of a hemisphere or more).
    cos_reach: Vec<f64>,
    /// The lowest vertex height (m), at most 0.
    min_height_m: f64,
    /// The highest vertex height (m), at least 0.
    max_height_m: f64,
    /// The smallest `cos ρ_v` over the cells (without the margin): the
    /// factor by which flat triangles may sag below their vertices.
    min_cos_rho: f64,
    /// The start cell of the last walk: the start of the next.
    hint: u32,
    /// Per cell: the number of the last walk that examined it.
    seen: Vec<u32>,
    /// The number of the current walk.
    walk: u32,
}

impl Picker {
    /// The picker of `mesh`: each cell's cap and the extreme vertex heights.
    pub fn new(mesh: &GlobeMesh) -> Self {
        let cells = mesh.cell_vertex_offsets.len() - 1;
        let mut cos_reach = Vec::with_capacity(cells);
        let mut min_cos_rho: f64 = 1.0;
        for range in mesh.cell_vertex_offsets.windows(2) {
            let (centre, end) = (range[0] as usize, range[1] as usize);
            let c = direction(mesh.vertices[centre]);
            let mut lowest = 1.0;
            for &v in &mesh.vertices[centre + 1..end] {
                let dot = c.dot(direction(v));
                if dot < lowest {
                    lowest = dot;
                }
            }
            if lowest < min_cos_rho {
                min_cos_rho = lowest;
            }
            // ρ = acos(lowest); rounding may put the dot product just
            // outside [−1, 1].
            let reach_rad = det::acos(lowest.clamp(-1.0, 1.0)) + CAP_MARGIN_RAD;
            cos_reach.push(det::cos(reach_rad));
        }
        let (mut min_height_m, mut max_height_m) = (0.0, 0.0);
        for v in &mesh.vertices {
            let h = f64::from(v[3]);
            if h < min_height_m {
                min_height_m = h;
            }
            if h > max_height_m {
                max_height_m = h;
            }
        }
        Self {
            cos_reach,
            min_height_m,
            max_height_m,
            min_cos_rho,
            hint: 0,
            seen: vec![0; cells],
            walk: 0,
        }
    }

    /// The cell whose drawn triangle the ray from `origin` (planet radii,
    /// outside the drawn surface) in the unit direction `dir` meets first,
    /// for the globe `mesh` of `grid` drawn with `relief_scale_per_m` (the
    /// frame's `σ`), or `None` if the ray misses the globe.
    ///
    /// # Panics
    ///
    /// If `mesh` is not the mesh of `grid` this picker was made for; in debug
    /// builds also if the ray is not finite.
    pub fn pick(
        &mut self,
        grid: &SphereGrid,
        mesh: &GlobeMesh,
        relief_scale_per_m: f32,
        origin: Vec3,
        dir: Vec3,
    ) -> Option<Hit> {
        let sigma = f64::from(relief_scale_per_m);
        let top = (1.0 + sigma * self.max_height_m) * (1.0 + SHELL_MARGIN);
        // The lowest vertex's radius times the largest sag (module
        // documentation, step 1); 0 where that bound does not hold.
        let lowest = 1.0 + sigma * self.min_height_m;
        let bottom = if lowest > 0.0 && self.min_cos_rho > 0.0 {
            lowest * self.min_cos_rho * (1.0 - SHELL_MARGIN)
        } else {
            0.0
        };
        let (entry, exit) = sphere(origin, dir, top)?;
        if exit <= 0.0 {
            return None;
        }
        if bottom < MIN_BOTTOM_RADIUS {
            let mut best = None;
            for cell in 0..self.cos_reach.len() as u32 {
                best = closer(best, cell, cell_hit(mesh, sigma, cell, origin, dir));
            }
            return best;
        }
        let t0 = entry.max(0.0);
        let t1 = match sphere(origin, dir, bottom) {
            Some((inner, _)) if inner >= t0 => inner,
            _ => exit,
        };
        let arc = Arc::new(origin, dir, t0, t1);
        let start = grid.nearest(arc.a, self.hint);
        self.hint = start;
        let walk = self.next_walk();
        self.seen[start as usize] = walk;
        let mut stack = vec![start];
        let mut best = None;
        while let Some(cell) = stack.pop() {
            best = closer(best, cell, cell_hit(mesh, sigma, cell, origin, dir));
            for &next in grid.neighbors(cell) {
                let k = next as usize;
                if self.seen[k] != walk {
                    self.seen[k] = walk;
                    let centre = direction(mesh.vertices[mesh.cell_vertex_offsets[k] as usize]);
                    if arc.reaches(centre, self.cos_reach[k]) {
                        stack.push(next);
                    }
                }
            }
        }
        best
    }

    /// The number of a new walk, every cell unseen by it.
    fn next_walk(&mut self) -> u32 {
        self.walk = self.walk.wrapping_add(1);
        if self.walk == 0 {
            self.seen.fill(0);
            self.walk = 1;
        }
        self.walk
    }
}

/// The great-circle arc swept by the directions of a ray segment.
///
/// Its plane's normal is `n = (o × d)/|o × d|`: `x(t) × x′(t) = o × d` for
/// `x(t) = o + t·d`, so the arc runs from `a` to `b` counter-clockwise
/// around `n`. The computed `o × d` is off by about `2ε|o|` (`ε = 2⁻⁵³`), so
/// the plane tilts by at most `2ε|o|/|o × d|`; an arc longer than
/// [`SHORT_ARC_SIN`] from a ray that reaches the radius `r ≥`
/// [`MIN_BOTTOM_RADIUS`] passes at `|o × d| ≥ 10⁻⁷·r` from the centre, so
/// for `|o| ≤ 20` the tilt stays below 5·10⁻⁶ rad, inside
/// [`CAP_MARGIN_RAD`].
struct Arc {
    /// The direction where the segment starts.
    a: Vec3,
    /// The direction where it ends.
    b: Vec3,
    /// The unit normal of its plane, `None` for an arc tested by its
    /// endpoints alone.
    normal: Option<Vec3>,
}

impl Arc {
    /// The arc of the ray `origin + t·dir` for `t` from `t0` to `t1`, which
    /// stays clear of the centre.
    fn new(origin: Vec3, dir: Vec3, t0: f64, t1: f64) -> Self {
        let a = origin + dir * t0;
        let b = origin + dir * t1;
        let (a, b) = (a / a.length(), b / b.length());
        let m = origin.cross(dir);
        let short = a.dot(b) > 0.0 && a.cross(b).length() < SHORT_ARC_SIN;
        let normal = (!short && m.length() > 0.0).then(|| m / m.length());
        Self { a, b, normal }
    }

    /// Whether the cap around the unit vector `c` with the radius whose
    /// cosine is `cos_reach` meets the arc: an endpoint in the cap, or the
    /// point of the great circle nearest `c` inside the arc (`c` in the
    /// wedge between `a` and `b`) and in the cap (`cos δ = √(1 − (c·n)²)`).
    fn reaches(&self, c: Vec3, cos_reach: f64) -> bool {
        if cos_reach <= 0.0 || c.dot(self.a) >= cos_reach || c.dot(self.b) >= cos_reach {
            return true;
        }
        let Some(n) = self.normal else {
            return false;
        };
        let s = c.dot(n);
        self.a.cross(c).dot(n) >= 0.0
            && c.cross(self.b).dot(n) >= 0.0
            && 1.0 - s * s >= cos_reach * cos_reach
    }
}

/// The distances along the ray `origin + t·dir` (`|dir| = 1`) at which it
/// enters and leaves the sphere of radius `radius` around the centre, if it
/// meets it: the roots of `t² + 2(o·d)t + |o|² − r² = 0`.
fn sphere(origin: Vec3, dir: Vec3, radius: f64) -> Option<(f64, f64)> {
    let b = origin.dot(dir);
    let disc = b * b - (origin.length_squared() - radius * radius);
    if disc.is_nan() || disc < 0.0 {
        return None;
    }
    let root = disc.sqrt();
    Some((-b - root, -b + root))
}

/// `hit` (the distance to a triangle of `cell`, if any) or `best`, whichever
/// is nearer; ties to the lower cell index.
fn closer(best: Option<Hit>, cell: u32, hit: Option<f64>) -> Option<Hit> {
    let Some(distance) = hit else {
        return best;
    };
    let new = Hit { cell, distance };
    match best {
        None => Some(new),
        Some(old) => match distance.total_cmp(&old.distance) {
            Ordering::Less => Some(new),
            Ordering::Equal if cell < old.cell => Some(new),
            _ => Some(old),
        },
    }
}

/// The unit direction of a mesh vertex.
fn direction(v: [f32; 4]) -> Vec3 {
    let p = Vec3::new(f64::from(v[0]), f64::from(v[1]), f64::from(v[2]));
    p / p.length()
}

/// A mesh vertex where the shader draws it: `(x, y, z)·(1 + σ·h)`.
fn drawn(v: [f32; 4], sigma: f64) -> Vec3 {
    Vec3::new(f64::from(v[0]), f64::from(v[1]), f64::from(v[2])) * (1.0 + sigma * f64::from(v[3]))
}

/// The distance along the ray to the nearest front face among the drawn
/// triangles of `cell`, if the ray meets one. Cell `v`'s triangles are
/// `indices[3(o_v − v)..3(o_{v+1} − v − 1)]` with `o` the vertex offsets: one
/// vertex per cell more than triangles.
fn cell_hit(mesh: &GlobeMesh, sigma: f64, cell: u32, origin: Vec3, dir: Vec3) -> Option<f64> {
    let v = cell as usize;
    let first = mesh.cell_vertex_offsets[v] as usize - v;
    let end = mesh.cell_vertex_offsets[v + 1] as usize - v - 1;
    let mut nearest: Option<f64> = None;
    for t in mesh.indices[3 * first..3 * end].chunks_exact(3) {
        let [p0, p1, p2] = [t[0], t[1], t[2]].map(|k| drawn(mesh.vertices[k as usize], sigma));
        if let Some(distance) = front_hit(origin, dir, p0, p1, p2) {
            if nearest.is_none_or(|n| distance < n) {
                nearest = Some(distance);
            }
        }
    }
    nearest
}

/// The distance `t > 0` at which the ray `origin + t·dir` meets the front of
/// the triangle `(p0, p1, p2)` (counter-clockwise seen from the front), by
/// the Möller–Trumbore algorithm (Möller and Trumbore 1997, *J. Graphics
/// Tools* 2(1), 21–28), with the barycentric bounds inclusive so that edges
/// and corners count; `None` for a back face, a ray in the triangle's plane
/// or a miss. `det = e₁·(d × e₂) = −d·(e₁ × e₂)` is positive for a front
/// face.
fn front_hit(origin: Vec3, dir: Vec3, p0: Vec3, p1: Vec3, p2: Vec3) -> Option<f64> {
    let e1 = p1 - p0;
    let e2 = p2 - p0;
    let p = dir.cross(e2);
    let det = e1.dot(p);
    if det.is_nan() || det <= 0.0 {
        return None;
    }
    let s = origin - p0;
    let u = s.dot(p);
    if u < 0.0 || u > det {
        return None;
    }
    let q = s.cross(e1);
    let v = dir.dot(q);
    if v < 0.0 || u + v > det {
        return None;
    }
    let t = e2.dot(q) / det;
    (t > 0.0).then_some(t)
}

#[cfg(test)]
mod tests {
    use kuniumi_sim::stages::crust::SURFACE_HEIGHT;
    use kuniumi_sphere::coords;
    use kuniumi_world::{PlanetParams, Settings, World};

    use super::*;
    use crate::camera::OrbitCamera;

    /// Earth's radius (m), as in `PlanetParams::earth`.
    const EARTH_RADIUS_M: f64 = 6_371_000.0;

    /// The exaggerations tested.
    const RELIEFS: [f64; 4] = [0.0, 1.0, 20.0, 100.0];

    fn world() -> World {
        let settings = Settings {
            cells: 1_000,
            allow_partial: true,
            ..Settings::default()
        };
        kuniumi_sim::generate(&PlanetParams::earth(), &settings, 3, &mut |_| {}).unwrap()
    }

    /// The heights above the sea level of `world`.
    fn above_sea(world: &World) -> Vec<f64> {
        let heights = &world.field(SURFACE_HEIGHT).unwrap().values;
        heights.iter().map(|h| h - world.sea_level_m).collect()
    }

    /// `σ = e/R` as the frame hands it to the GPU.
    fn sigma(relief: f64) -> f32 {
        (relief / EARTH_RADIUS_M) as f32
    }

    fn drawn(v: [f32; 4], sigma: f32) -> Vec3 {
        Vec3::new(f64::from(v[0]), f64::from(v[1]), f64::from(v[2]))
            * (1.0 + f64::from(sigma) * f64::from(v[3]))
    }

    /// The drawn triangles of `cell`.
    fn triangles(mesh: &GlobeMesh, sigma: f32, cell: u32) -> Vec<[Vec3; 3]> {
        let first = mesh.cell_vertex_offsets[cell as usize] as usize - cell as usize;
        let last = mesh.cell_vertex_offsets[cell as usize + 1] as usize - cell as usize - 1;
        mesh.indices[3 * first..3 * last]
            .chunks_exact(3)
            .map(|t| [t[0], t[1], t[2]].map(|k| drawn(mesh.vertices[k as usize], sigma)))
            .collect()
    }

    /// A front-facing hit of the reference: its distance, its cell, and
    /// whether it lies inside the triangle without the tolerance.
    #[derive(Clone, Copy, Debug)]
    struct Reference {
        distance: f64,
        cell: u32,
        strict: bool,
    }

    /// Every drawn triangle of `mesh` with its cell.
    fn all_triangles(mesh: &GlobeMesh, sigma: f32) -> Vec<([Vec3; 3], u32)> {
        (0..mesh.cell_vertex_offsets.len() as u32 - 1)
            .flat_map(|cell| {
                triangles(mesh, sigma, cell)
                    .into_iter()
                    .map(move |t| (t, cell))
            })
            .collect()
    }

    /// Every front-facing hit of the ray on the drawn `triangles`, by brute
    /// force and with formulas independent of the picker's: the line meets
    /// a triangle if the three signed volumes `d·((pᵢ − o) × (pⱼ − o))` of
    /// its edges agree in sign (within a tolerance relative to their scale,
    /// for hits on an edge), and the distance is that to the triangle's
    /// plane.
    fn reference(triangles: &[([Vec3; 3], u32)], o: Vec3, d: Vec3) -> Vec<Reference> {
        let mut hits = Vec::new();
        for &([p0, p1, p2], cell) in triangles {
            let normal = (p1 - p0).cross(p2 - p0);
            let facing = d.dot(normal);
            if facing >= 0.0 {
                continue;
            }
            let distance = normal.dot(p0 - o) / facing;
            if distance <= 0.0 {
                continue;
            }
            let (a, b, c) = (p0 - o, p1 - o, p2 - o);
            let volumes = [d.dot(a.cross(b)), d.dot(b.cross(c)), d.dot(c.cross(a))];
            let tolerance = 1e-9 * (a.length() * b.length() + b.length() * c.length());
            let inside = |sign: f64| volumes.iter().all(|&v| v * sign >= -tolerance);
            let strictly = |sign: f64| volumes.iter().all(|&v| v * sign >= 0.0);
            if inside(1.0) || inside(-1.0) {
                hits.push(Reference {
                    distance,
                    cell,
                    strict: strictly(1.0) || strictly(-1.0),
                });
            }
        }
        hits
    }

    /// Checks `hit` against the reference: a miss only where no triangle is
    /// hit strictly; otherwise the smallest distance (within 1e-12 relative)
    /// on a triangle of the picked cell.
    fn check(hit: Option<Hit>, triangles: &[([Vec3; 3], u32)], o: Vec3, d: Vec3, what: &str) {
        let all = reference(triangles, o, d);
        let nearest = |strict: bool| {
            all.iter()
                .filter(|r| r.strict || !strict)
                .map(|r| r.distance)
                .fold(f64::INFINITY, f64::min)
        };
        let (strict_min, loose_min) = (nearest(true), nearest(false));
        match hit {
            None => assert!(
                strict_min.is_infinite(),
                "{what}: missed a hit at {strict_min}"
            ),
            Some(hit) => {
                let tolerance = 1e-12 * hit.distance;
                assert!(
                    hit.distance <= strict_min + tolerance && hit.distance >= loose_min - tolerance,
                    "{what}: {hit:?}, reference from {loose_min} to {strict_min}"
                );
                assert!(
                    all.iter()
                        .any(|r| r.cell == hit.cell
                            && (r.distance - hit.distance).abs() <= tolerance),
                    "{what}: {hit:?} is not on a triangle of its cell"
                );
            }
        }
    }

    #[test]
    fn octahedron_by_hand() {
        // Cell 0 (+X) 1,024 m high, σ = 2⁻¹² per metre (exact in f32), so
        // its centre is drawn at radius 1.25; the others at 0.
        let grid = SphereGrid::octahedron();
        let mesh = GlobeMesh::new(&grid, &[1024.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        let sigma = 1.0 / 4096.0;
        let mut picker = Picker::new(&mesh);
        let origin = Vec3::new(3.0, 0.0, 0.0);
        let hit = picker.pick(&grid, &mesh, sigma, origin, -Vec3::X).unwrap();
        assert_eq!(hit.cell, 0);
        assert!((hit.distance - 1.75).abs() < 1e-12, "{hit:?}");
        // Straight down onto −Z from above it: radius 1, distance 2.
        let hit = picker
            .pick(&grid, &mesh, sigma, Vec3::new(0.0, 0.0, -3.0), Vec3::Z)
            .unwrap();
        assert_eq!(hit.cell, 5);
        assert!((hit.distance - 2.0).abs() < 1e-12, "{hit:?}");
        // Away from the globe, and past it: the line from (3, 0, 0) towards
        // (0, 1.5, 0) passes 4.5/√11.25 ≈ 1.34 from the centre, above the
        // highest point at 1.25.
        assert_eq!(picker.pick(&grid, &mesh, sigma, origin, Vec3::X), None);
        let past = Vec3::new(-3.0, 1.5, 0.0);
        assert_eq!(
            picker.pick(&grid, &mesh, sigma, origin, past / past.length()),
            None
        );
        let all = all_triangles(&mesh, sigma);
        for d in [
            -Vec3::X,
            Vec3::new(-1.0, 0.3, 0.1),
            Vec3::new(-1.0, -0.2, 0.4),
        ] {
            let d = d / d.length();
            let hit = picker.pick(&grid, &mesh, sigma, origin, d);
            check(hit, &all, origin, d, "octahedron");
        }
    }

    #[test]
    fn back_faces_in_front_are_not_picked() {
        // Cell 0 (+X) 9,830 m deep at σ = 2⁻¹² per metre: its centre drawn at
        // x = 1 − 2.4 ≈ −1.4, through the centre and beyond cell 1 (−X, at
        // x = −1); its corners at radius 1 − 0.8. Flipped through the plane of
        // its corners, its triangles face inwards: seen from (−3, 0, 0) they
        // are back faces, which the GPU culls, so the pick is cell 1, 2 away.
        let grid = SphereGrid::octahedron();
        let mesh = GlobeMesh::new(&grid, &[-9_830.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        let sigma = 1.0 / 4096.0;
        let mut picker = Picker::new(&mesh);
        let origin = Vec3::new(-3.0, 0.0, 0.0);
        let hit = picker.pick(&grid, &mesh, sigma, origin, Vec3::X).unwrap();
        assert_eq!(hit.cell, 1);
        assert!((hit.distance - 2.0).abs() < 1e-12, "{hit:?}");
    }

    #[test]
    fn rays_grazing_the_corners_of_cells_agree_with_brute_force() {
        // A ray that meets a cell only near one corner, moving across the
        // corner rather than into the cell, is found only if the cell's cap
        // reaches that corner: each ray comes in at 60° from the vertical
        // towards a point of the drawn triangle 10⁻⁴ (barycentric) from the
        // corner, moving perpendicular to the line from the cell's centre to
        // the corner, in both senses.
        let grid = SphereGrid::fibonacci(500, 11).unwrap();
        let heights: Vec<f64> = (0..grid.len())
            .map(|i| (i % 17) as f64 * 731.3 - 5_000.1)
            .collect();
        let mesh = GlobeMesh::new(&grid, &heights);
        let mut picker = Picker::new(&mesh);
        let (sin_in, cos_in) = (0.75_f64.sqrt(), 0.5);
        for relief in RELIEFS {
            let sigma = sigma(relief);
            let all = all_triangles(&mesh, sigma);
            for cell in 0..grid.len() as u32 {
                let c = drawn(
                    mesh.vertices[mesh.cell_vertices(cell).start as usize],
                    sigma,
                );
                for [p0, p1, p2] in triangles(&mesh, sigma, cell) {
                    let target = p0 * 1e-4 + p1 * (1.0 - 2e-4) + p2 * 1e-4;
                    let up = target / target.length();
                    let across = p1 / p1.length() - c / c.length();
                    let across = across - up * across.dot(up);
                    let tangent = up.cross(across);
                    let tangent = tangent / tangent.length();
                    for sense in [1.0, -1.0] {
                        let back = up * cos_in + tangent * (sense * sin_in);
                        let (o, d) = (target + back * 2.0, -back);
                        let hit = picker.pick(&grid, &mesh, sigma, o, d);
                        check(hit, &all, o, d, &format!("×{relief} corner of {cell}"));
                    }
                }
            }
        }
    }

    #[test]
    fn radial_rays_find_their_cell() {
        // A ray towards the centre meets the radial graph once, so the cell
        // of the direction it comes from is known without any search: the
        // centroid direction of each drawn triangle lies inside it.
        let world = world();
        let mesh = GlobeMesh::new(&world.grid, &above_sea(&world));
        let mut picker = Picker::new(&mesh);
        for relief in RELIEFS {
            let sigma = sigma(relief);
            for cell in 0..world.grid.len() as u32 {
                for (m, [p0, p1, p2]) in triangles(&mesh, sigma, cell).into_iter().enumerate() {
                    let g = p0 + p1 + p2;
                    let g = g / g.length();
                    let hit = picker.pick(&world.grid, &mesh, sigma, g * 3.0, -g);
                    assert_eq!(
                        hit.map(|h| h.cell),
                        Some(cell),
                        "×{relief}, cell {cell}, triangle {m}"
                    );
                    // On the near side, in the direction ĝ.
                    let point = g * (3.0 - hit.unwrap().distance);
                    assert!((point / point.length() - g).length() < 1e-12);
                    // Through the camera, for some of them.
                    if m == 0 && cell % 10 == 0 {
                        let (lat_rad, lon_rad) = coords::lat_lon(g);
                        let camera = OrbitCamera {
                            lat_rad,
                            lon_rad,
                            distance: 3.0,
                        };
                        let (o, d) = camera.ray(1.6, [0.0, 0.0]);
                        let hit = picker.pick(&world.grid, &mesh, sigma, o, d);
                        assert_eq!(hit.map(|h| h.cell), Some(cell), "×{relief}, camera {cell}");
                    }
                }
            }
        }
    }

    /// Sweeps 21 × 21 cursor positions of several cameras at every
    /// exaggeration and checks each pick against brute force.
    fn sweep(grid: &SphereGrid, heights: &[f64], what: &str) {
        let mesh = GlobeMesh::new(grid, heights);
        let mut picker = Picker::new(&mesh);
        let max_height_m = mesh
            .vertices
            .iter()
            .map(|v| f64::from(v[3]))
            .fold(0.0, f64::max);
        let mut hits = 0;
        let mut rays = 0;
        for relief in RELIEFS {
            let sigma = sigma(relief);
            let all = all_triangles(&mesh, sigma);
            let top = 1.0 + f64::from(sigma) * max_height_m;
            let cameras = [
                (20.0, 0.0, 3.0),
                (-60.0, 135.0, 1.5),
                (5.0, -100.0, top + 0.02),
                (89.0, 10.0, 2.0),
            ];
            for (lat_deg, lon_deg, distance) in cameras {
                let camera = OrbitCamera {
                    lat_rad: kuniumi_math::units::deg_to_rad(lat_deg),
                    lon_rad: kuniumi_math::units::deg_to_rad(lon_deg),
                    distance,
                };
                for i in 0..21 {
                    for j in 0..21 {
                        let ndc = [f64::from(i) / 10.0 - 1.0, f64::from(j) / 10.0 - 1.0];
                        let (o, d) = camera.ray(1.6, ndc);
                        let hit = picker.pick(grid, &mesh, sigma, o, d);
                        let label = format!("{what} ×{relief} {camera:?} {ndc:?}");
                        check(hit, &all, o, d, &label);
                        rays += 1;
                        hits += usize::from(hit.is_some());
                    }
                }
            }
        }
        // Both hits and misses were checked.
        assert!(
            hits > rays / 3 && hits < rays,
            "{what}: {hits} hits of {rays} rays"
        );
    }

    #[test]
    fn camera_rays_agree_with_brute_force() {
        let world = world();
        sweep(&world.grid, &above_sea(&world), "world");
    }

    #[test]
    fn occluding_spikes_agree_with_brute_force() {
        // Every seventh cell 20 km higher: at ×100 spikes of 0.3 radii that
        // hide the cells behind them from grazing rays.
        let world = world();
        let heights: Vec<f64> = above_sea(&world)
            .iter()
            .enumerate()
            .map(|(i, h)| if i % 7 == 0 { h + 20_000.0 } else { *h })
            .collect();
        sweep(&world.grid, &heights, "spikes");
    }

    /// Rays from `distance` radii that pass `closest` radii from the centre,
    /// from origins along `origins` directions, each in `azimuths`
    /// directions around the line to the centre.
    fn rays_past(distance: f64, closest: f64, origins: u32, azimuths: u32) -> Vec<(Vec3, Vec3)> {
        let points = SphereGrid::fibonacci(origins, 5).unwrap();
        let mut rays = Vec::new();
        for i in 0..origins {
            let o = points.position(i);
            let east = coords::east(o);
            let north = coords::north(o);
            let sin = closest / distance;
            let cos = (1.0 - sin * sin).sqrt();
            for k in 0..azimuths {
                let angle = core::f64::consts::TAU * f64::from(k) / f64::from(azimuths);
                let (sin_a, cos_a) = det::sin_cos(angle);
                let side = east * cos_a + north * sin_a;
                rays.push((o * distance, -o * cos + side * sin));
            }
        }
        rays
    }

    #[test]
    fn limb_rays_in_the_sag_of_flat_triangles_agree_with_brute_force() {
        // A flat triangle sags below its vertices, by about 7·10⁻⁴ radii on
        // a 1,000-cell grid (review R1-1): rays passing just inside the
        // vertices' sphere at the limb meet it beyond the vertices' lowest
        // radius. At ×0 every vertex lies at radius 1; at ×1 nearly so.
        let world = world();
        let mesh = GlobeMesh::new(&world.grid, &above_sea(&world));
        let mut picker = Picker::new(&mesh);
        for relief in [0.0, 1.0] {
            let sigma = sigma(relief);
            let all = all_triangles(&mesh, sigma);
            let mut hits = 0;
            for k in 0..16 {
                let closest = 1.0 - 2e-4 * f64::from(k);
                for (o, d) in rays_past(3.0, closest, 40, 8) {
                    let hit = picker.pick(&world.grid, &mesh, sigma, o, d);
                    check(hit, &all, o, d, &format!("×{relief} limb {closest}"));
                    hits += usize::from(hit.is_some());
                }
            }
            assert!(hits > 0, "×{relief}: no limb ray hit");
        }
    }

    #[test]
    fn oblique_rays_on_the_flat_octahedron_agree_with_brute_force() {
        // The octahedron's flat fan triangles sag to about 0.81 radii
        // between their vertices on the unit sphere (review R1-1); the
        // picker's bound for them is min cos ρ = 1/√3.
        let grid = SphereGrid::octahedron();
        let mesh = GlobeMesh::new(&grid, &[0.0; 6]);
        let mut picker = Picker::new(&mesh);
        let all = all_triangles(&mesh, 0.0);
        for k in 0..=10 {
            let closest = 0.1 * f64::from(k);
            for (o, d) in rays_past(3.0, closest, 50, 8) {
                let hit = picker.pick(&grid, &mesh, 0.0, o, d);
                check(
                    hit,
                    &all,
                    o,
                    d,
                    &format!("octahedron, {closest} from the centre"),
                );
            }
        }
    }

    #[test]
    fn ties_go_to_the_lower_cell() {
        let at = |cell, distance| Some(Hit { cell, distance });
        assert_eq!(closer(None, 7, None), None);
        assert_eq!(closer(None, 7, Some(2.0)), at(7, 2.0));
        assert_eq!(closer(at(5, 1.0), 3, Some(1.0)), at(3, 1.0));
        assert_eq!(closer(at(3, 1.0), 5, Some(1.0)), at(3, 1.0));
        assert_eq!(closer(at(3, 1.0), 9, Some(0.5)), at(9, 0.5));
        assert_eq!(closer(at(3, 1.0), 1, Some(1.5)), at(3, 1.0));
        assert_eq!(closer(at(3, 1.0), 1, None), at(3, 1.0));
    }

    #[test]
    fn inside_out_relief_tests_every_triangle() {
        // Every fifth cell 100 km deep: at ×100 drawn at radius
        // 1 − 100 · 10⁵/6.371·10⁶ < 0, through the centre.
        let world = world();
        let heights: Vec<f64> = (0..world.grid.len())
            .map(|i| if i % 5 == 0 { -100_000.0 } else { 1_000.0 })
            .collect();
        let mesh = GlobeMesh::new(&world.grid, &heights);
        let mut picker = Picker::new(&mesh);
        let sigma = sigma(100.0);
        assert!(1.0 + f64::from(sigma) * -100_000.0 < 0.0);
        let all = all_triangles(&mesh, sigma);
        for (lat_deg, lon_deg) in [(10.0, 20.0), (-45.0, -170.0)] {
            let camera = OrbitCamera {
                lat_rad: kuniumi_math::units::deg_to_rad(lat_deg),
                lon_rad: kuniumi_math::units::deg_to_rad(lon_deg),
                distance: 3.0,
            };
            for i in 0..11 {
                for j in 0..11 {
                    let ndc = [f64::from(i) / 5.0 - 1.0, f64::from(j) / 5.0 - 1.0];
                    let (o, d) = camera.ray(1.6, ndc);
                    let hit = picker.pick(&world.grid, &mesh, sigma, o, d);
                    check(hit, &all, o, d, &format!("inside out {ndc:?}"));
                }
            }
        }
    }
}

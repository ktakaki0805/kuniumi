//! The spherical Delaunay triangulation of the cell centres (design document,
//! 球面グリッドとデータ構造; D-14).
//!
//! [`triangulate`] returns the triangles of the 3D convex hull of points on
//! the unit sphere, with the neighbouring triangle across each edge. The
//! hull of points on a sphere is their spherical Delaunay triangulation: a
//! plane cuts the sphere in a circle, and the points of the sphere strictly
//! inside the circle (on the side away from the centre) lie strictly above
//! the plane. So a triangle of the points has an empty circumcircle exactly
//! when no point lies above its plane, that is when the plane supports the
//! hull (K. Q. Brown 1979, "Voronoi diagrams from convex hulls", *Inf.
//! Process. Lett.* 9(5), 223–228). The Voronoi cells, the grid's cells, are
//! the dual (ticket M1-06).
//!
//! # Algorithm
//!
//! Randomised incremental construction with conflict lists (M. de Berg,
//! O. Cheong, M. van Kreveld and M. Overmars 2008, *Computational Geometry:
//! Algorithms and Applications*, 3rd ed., Springer, chapter 11):
//!
//! 1. The indices are shuffled (Fisher–Yates, [`STREAM`]). The first four
//!    points form a tetrahedron; every other point is attached to one face
//!    of it that it sees (lies strictly above), in a linked list per face.
//! 2. The points are inserted in the shuffled order. The faces that point
//!    `p` sees form a connected region around its attached face; a
//!    depth-first search finds them. The edges between them and the faces
//!    `p` does not see form a simple cycle, the horizon; each horizon edge
//!    and `p` form a new face, and the visible faces are deleted.
//! 3. Each pending point of a deleted face is attached to the first new face
//!    it sees. If it sees none, it lies inside the new hull and is dropped:
//!    it sees the deleted face `f` and lies below every new face, so it lies
//!    in the cone from `p` over the horizon,
//!    `p + Σ aᵢ (hᵢ − p)` with `aᵢ ≥ 0`. Its height above `f`'s plane is
//!    `(1 − Σ aᵢ) h(p) + Σ aᵢ h(hᵢ)` with `h(p) > 0` and `h(hᵢ) ≤ 0`; it is
//!    positive only if `Σ aᵢ < 1`, which puts the point in the convex hull of
//!    `p` and the horizon.
//!
//! The book keeps every visible face of a pending point; here each point
//! keeps one and the others are found by the search. A point is reattached
//! only when its face is deleted, and then it sees a new face, so the number
//! of reattachments is at most the number of changes of the book's conflict
//! graph, expected O(n log n) over the insertion orders (the book's
//! analysis); each costs at most one test per new face, that is the degree
//! of the new vertex. With bounded degrees (at most 9 measured on the
//! jittered Fibonacci grid) the total is therefore O(n log n). Measured in
//! release (ticket M1-05): 4.5 to 4.8 visibility tests per `n log₂ n` from
//! 10⁴ to 2.6 × 10⁶ Fibonacci points, and about 2.5 s for 650,000 points.
//! Inputs with thousands of exactly co-circular points, whose triangulations
//! have vertices of very high degree, can be much slower: points on the
//! equator and on ±45° (k per circle, longitudes offset by half a step) and
//! the poles take 3.8 s for 48,002 points and 48 s for 192,002. Faces live
//! in flat arrays, freed slots are reused, and the conflict lists are linked
//! through per-point arrays: no allocation per face, no hashing.
//!
//! # Degenerate inputs
//!
//! Every visibility decision is the exact [`orient3d`], and an exact zero
//! (four points on one circle of the sphere, such as a face of a cube, or
//! four points on a great circle) is resolved by Simulation of Simplicity on
//! the point indices (`predicates::orient3d_sos`, H. Edelsbrunner and E. P.
//! Mücke 1990): every decision is that of one infinitesimally perturbed
//! configuration in general position. The result is the hull of the
//! perturbed points, unique whatever the insertion order: where several
//! triangulations of co-circular points are valid (each is Delaunay), the
//! indices decide which one is returned.
//!
//! # Output and determinism
//!
//! The triangles are counter-clockwise seen from outside, each rotated so
//! that its smallest index comes first, and sorted. The output is therefore
//! a function of the points alone: the seed decides only the insertion order
//! and so the running time. Everything else is integer bookkeeping and exact
//! predicates, so the result is the same on every platform.
//!
//! # Checks
//!
//! On the sphere every point is a vertex of the hull, and the hull contains
//! the centre when the points do not all lie in a closed hemisphere. The
//! input is checked before and the result after construction; see
//! [`DelaunayError`] for the order.

use core::cmp::Ordering;
use core::fmt;

use kuniumi_math::{Rng, Vec3};

use crate::grid::UNIT_TOLERANCE;
use crate::predicates::{Orientation, in_exact_range, orient3d, orient3d_sos};

/// The random stream of the insertion order: [`triangulate`] shuffles with
/// `Rng::for_stream(seed, STREAM)`. Stream names are part of the output
/// contract (docs/conventions.md), although the output of [`triangulate`]
/// does not depend on the order.
pub const STREAM: &str = "grid/insertion";

/// The smallest distance between two points that [`triangulate`] accepts:
/// 1e-12 as a chord on the unit sphere (the same in radians to within
/// 1e-36). Closer points are rejected with [`DelaunayError::TooClose`]:
/// their triangles would be too thin for circumcentres to be accurate.
pub const MIN_SEPARATION: f64 = 1.0e-12;

/// The largest number of points [`triangulate`] accepts: 2³⁰. Face slots
/// (at most `3n`) and visit marks (`2n + 1`) then fit in `u32`.
pub const MAX_POINTS: usize = 1 << 30;

/// Marks a free face slot, a missing link or a point without a face.
const NONE: u32 = u32::MAX;

/// A triangulation of the sphere: triangles of point indices and the
/// triangle across each edge.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Triangulation {
    triangles: Vec<[u32; 3]>,
    neighbors: Vec<[u32; 3]>,
}

impl Triangulation {
    /// The triangles: point indices, counter-clockwise seen from outside,
    /// the smallest index first; sorted, so the output is canonical.
    /// There are `2n − 4` of them for `n` points.
    pub fn triangles(&self) -> &[[u32; 3]] {
        &self.triangles
    }

    /// The neighbours: `neighbors()[t][k]` is the triangle across the edge
    /// from `triangles()[t][k]` to `triangles()[t][(k + 1) % 3]`, which runs
    /// the other way in that triangle.
    pub fn neighbors(&self) -> &[[u32; 3]] {
        &self.neighbors
    }

    /// The number of triangles.
    pub fn len(&self) -> usize {
        self.triangles.len()
    }

    /// Whether there are no triangles (never true for a result of
    /// [`triangulate`]).
    pub fn is_empty(&self) -> bool {
        self.triangles.is_empty()
    }

    /// The triangles and the neighbours, as described at
    /// [`triangles`](Self::triangles) and [`neighbors`](Self::neighbors).
    pub fn into_parts(self) -> (Vec<[u32; 3]>, Vec<[u32; 3]>) {
        (self.triangles, self.neighbors)
    }
}

/// Why [`triangulate`] cannot triangulate a point set.
///
/// The checks run in the order of the variants: the number of points, then
/// each point in index order, duplicates, and after construction the
/// vertices, the edge lengths and the triangles. Every reported index is
/// chosen by a rule that does not depend on the insertion order.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DelaunayError {
    /// Fewer than four points: no closed triangulation has fewer vertices.
    TooFew {
        /// The number of points.
        n: usize,
    },
    /// More than [`MAX_POINTS`] points.
    TooMany {
        /// The number of points.
        n: usize,
    },
    /// The first point that is not finite, has a non-zero coordinate of
    /// magnitude below 2⁻²⁰⁰ (outside the exact range of [`orient3d`];
    /// round such coordinates to zero), or is not a unit vector:
    /// `| |p|² − 1 | > 1e-12`, the grid's tolerance.
    BadPoint {
        /// Its index.
        point: u32,
    },
    /// Two points are equal (−0 and +0 are equal): `second` is the smallest
    /// index equal to an earlier point, `first` the smallest index equal to
    /// it.
    Duplicate {
        /// The earlier point.
        first: u32,
        /// The later point.
        second: u32,
    },
    /// A point is not a vertex of the hull: it lies inside or on the hull of
    /// the others. On the unit sphere this happens only when points nearly
    /// coincide: near a vertex the hull lies about `d s / 2` below the sphere
    /// at a distance `d` from it (`s` the edge length), so a point is lost
    /// only if that is below its error of length. For the default grid
    /// (`s` ≈ 0.0044) that means within about 1e-13 for unit vectors exact
    /// to rounding, and about 2e-10 for the full tolerance. The smallest such
    /// index.
    NotOnHull {
        /// Its index.
        point: u32,
    },
    /// Two points joined by an edge are closer than [`MIN_SEPARATION`]
    /// (`|a − b|² < 1e-24`). Each point's nearest neighbour is joined to it
    /// by a Delaunay edge, so if any two vertices are that close, at least
    /// one such pair is an edge and is reported. The smallest such edge
    /// `(a, b)`, `a < b`.
    TooClose {
        /// The smaller index.
        a: u32,
        /// The larger index.
        b: u32,
    },
    /// The centre of the sphere is not strictly below the plane of a
    /// triangle: the points lie in a closed hemisphere, so the hull does not
    /// enclose the centre and its triangles do not cover the sphere. The
    /// first such triangle of the canonical order.
    Hemisphere {
        /// The triangle, as it would appear in
        /// [`Triangulation::triangles`].
        triangle: [u32; 3],
    },
}

impl fmt::Display for DelaunayError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TooFew { n } => {
                write!(
                    f,
                    "a triangulation of the sphere needs at least 4 points, got {n}"
                )
            }
            Self::TooMany { n } => {
                write!(
                    f,
                    "at most {MAX_POINTS} points can be triangulated, got {n}"
                )
            }
            Self::BadPoint { point } => write!(
                f,
                "point {point} is not a unit vector with coordinates in the exact range"
            ),
            Self::Duplicate { first, second } => {
                write!(f, "point {second} is equal to point {first}")
            }
            Self::NotOnHull { point } => write!(
                f,
                "point {point} is not a vertex of the convex hull (nearly coincident points?)"
            ),
            Self::TooClose { a, b } => write!(
                f,
                "points {a} and {b} are closer than {MIN_SEPARATION:e} on the unit sphere"
            ),
            Self::Hemisphere { triangle } => write!(
                f,
                "the points lie in a closed hemisphere: the centre is not below triangle {triangle:?}"
            ),
        }
    }
}

impl std::error::Error for DelaunayError {}

/// The spherical Delaunay triangulation of `points`, unit vectors, as their
/// 3D convex hull; see the module documentation.
///
/// The insertion order is shuffled with `Rng::for_stream(seed, STREAM)`;
/// the result does not depend on it (nor, therefore, on `seed`), only the
/// running time does. Points on a common circle are triangulated by a rule
/// on their indices, so the result is always one valid triangulation, the
/// same on every run and platform.
///
/// # Errors
///
/// See [`DelaunayError`], whose variants are checked in order.
pub fn triangulate(points: &[Vec3], seed: u64) -> Result<Triangulation, DelaunayError> {
    let n = check_count(points.len())?;
    check_points(points)?;
    check_duplicates(points)?;
    let order = insertion_order(n, seed);
    let mut hull = Hull::new(points, &order);
    for &p in &order[4..] {
        hull.insert(p);
    }
    let triangulation = hull.into_triangulation();
    check_vertices(points.len(), &triangulation)?;
    check_separation(points, &triangulation)?;
    check_centre(points, &triangulation)?;
    Ok(triangulation)
}

/// The number of points as `u32`, if it is in [4, [`MAX_POINTS`]].
fn check_count(n: usize) -> Result<u32, DelaunayError> {
    if n < 4 {
        Err(DelaunayError::TooFew { n })
    } else if n > MAX_POINTS {
        Err(DelaunayError::TooMany { n })
    } else {
        Ok(n as u32)
    }
}

/// [`DelaunayError::BadPoint`] for the first point that is not a unit
/// vector with coordinates in `orient3d`'s exact range.
fn check_points(points: &[Vec3]) -> Result<(), DelaunayError> {
    for (i, p) in points.iter().enumerate() {
        // `in_exact_range` is false for NaN and infinities.
        let in_range = in_exact_range(p.x) && in_exact_range(p.y) && in_exact_range(p.z);
        if !in_range || (p.length_squared() - 1.0).abs() > UNIT_TOLERANCE {
            return Err(DelaunayError::BadPoint { point: i as u32 });
        }
    }
    Ok(())
}

/// The coordinates of `p` with −0 replaced by +0 (`−0 + 0 = +0`), so that
/// `total_cmp` orders equal points as equal.
fn duplicate_key(p: Vec3) -> [f64; 3] {
    [p.x + 0.0, p.y + 0.0, p.z + 0.0]
}

/// Lexicographic `total_cmp` of two keys.
fn compare_keys(a: &[f64; 3], b: &[f64; 3]) -> Ordering {
    a[0].total_cmp(&b[0])
        .then(a[1].total_cmp(&b[1]))
        .then(a[2].total_cmp(&b[2]))
}

/// [`DelaunayError::Duplicate`] if two points are equal. The indices are
/// sorted by coordinates, then by index, so equal points form runs in index
/// order; the pair reported is the first two indices of the run whose
/// second index is smallest.
fn check_duplicates(points: &[Vec3]) -> Result<(), DelaunayError> {
    let keys: Vec<[f64; 3]> = points.iter().map(|&p| duplicate_key(p)).collect();
    let mut order: Vec<u32> = (0..points.len() as u32).collect();
    order.sort_unstable_by(|&i, &j| {
        compare_keys(&keys[i as usize], &keys[j as usize]).then(i.cmp(&j))
    });
    let mut found: Option<(u32, u32)> = None;
    let mut start = 0;
    for k in 1..=order.len() {
        let run_ends = k == order.len()
            || compare_keys(&keys[order[k] as usize], &keys[order[start] as usize])
                != Ordering::Equal;
        if run_ends {
            if k - start >= 2 {
                let pair = (order[start], order[start + 1]);
                if found.is_none_or(|(_, second)| pair.1 < second) {
                    found = Some(pair);
                }
            }
            start = k;
        }
    }
    match found {
        Some((first, second)) => Err(DelaunayError::Duplicate { first, second }),
        None => Ok(()),
    }
}

/// The insertion order: `0..n` shuffled by Durstenfeld's Fisher–Yates with
/// `Rng::for_stream(seed, STREAM)`: for `i` from `n − 1` down to 1, swap
/// element `i` with element `below(i + 1)`.
fn insertion_order(n: u32, seed: u64) -> Vec<u32> {
    let mut rng = Rng::for_stream(seed, STREAM);
    let mut order: Vec<u32> = (0..n).collect();
    for i in (1..order.len()).rev() {
        let j = rng.below(i as u64 + 1) as usize;
        order.swap(i, j);
    }
    order
}

/// [`DelaunayError::NotOnHull`] for the smallest index that is not a vertex.
fn check_vertices(n: usize, t: &Triangulation) -> Result<(), DelaunayError> {
    let mut used = vec![false; n];
    for tri in &t.triangles {
        for &v in tri {
            used[v as usize] = true;
        }
    }
    match used.iter().position(|&u| !u) {
        Some(point) => Err(DelaunayError::NotOnHull {
            point: point as u32,
        }),
        None => Ok(()),
    }
}

/// [`DelaunayError::TooClose`] for the smallest edge shorter than
/// [`MIN_SEPARATION`].
fn check_separation(points: &[Vec3], t: &Triangulation) -> Result<(), DelaunayError> {
    let limit = MIN_SEPARATION * MIN_SEPARATION;
    let mut found: Option<(u32, u32)> = None;
    for tri in &t.triangles {
        for k in 0..3 {
            let (a, b) = (tri[k], tri[(k + 1) % 3]);
            // Each edge appears once in each direction; take it once.
            if a < b
                && (points[a as usize] - points[b as usize]).length_squared() < limit
                && found.is_none_or(|smallest| (a, b) < smallest)
            {
                found = Some((a, b));
            }
        }
    }
    match found {
        Some((a, b)) => Err(DelaunayError::TooClose { a, b }),
        None => Ok(()),
    }
}

/// [`DelaunayError::Hemisphere`] for the first triangle whose plane does not
/// have the centre strictly below it.
fn check_centre(points: &[Vec3], t: &Triangulation) -> Result<(), DelaunayError> {
    for tri in &t.triangles {
        let [a, b, c] = tri.map(|v| points[v as usize]);
        if orient3d(a, b, c, Vec3::ZERO) != Orientation::Negative {
            return Err(DelaunayError::Hemisphere { triangle: *tri });
        }
    }
    Ok(())
}

/// The convex hull under construction.
///
/// A face is a slot in the per-face arrays; its edge `k` runs from vertex
/// `k` to vertex `(k + 1) % 3`.
struct Hull<'a> {
    points: &'a [Vec3],
    /// Per face slot: the vertices, counter-clockwise seen from outside, or
    /// `[NONE; 3]` for a free slot.
    faces: Vec<[u32; 3]>,
    /// Per face slot: the face across each edge.
    adjacent: Vec<[u32; 3]>,
    /// Per face slot: the first pending point attached to it, or `NONE`.
    conflict_head: Vec<u32>,
    /// Per face slot: `2s` if insertion `s` found it visible, `2s + 1` if
    /// hidden (0 initially; insertions count from 1).
    stamp: Vec<u32>,
    /// Free face slots, reused last in, first out.
    free: Vec<u32>,
    /// Per point: the next point attached to the same face, or `NONE`.
    conflict_next: Vec<u32>,
    /// Per point: the face it is attached to; `NONE` once it is a vertex or
    /// has been dropped as lying inside the hull.
    conflict_face: Vec<u32>,
    /// Per point, meaningful for the horizon vertices of the current
    /// insertion only: the new face whose edge 0 starts there.
    new_face_at: Vec<u32>,
    /// The number of insertions so far.
    step: u32,
    /// Scratch: the depth-first search stack.
    stack: Vec<u32>,
    /// Scratch: the faces visible from the point being inserted.
    visible: Vec<u32>,
    /// Scratch: the faces created by the current insertion.
    new_faces: Vec<u32>,
}

impl<'a> Hull<'a> {
    /// The tetrahedron of `order[0..4]`, with every later point of `order`
    /// attached to the first of its faces that it sees.
    fn new(points: &'a [Vec3], order: &[u32]) -> Self {
        let n = points.len();
        let mut hull = Self {
            points,
            faces: Vec::with_capacity(2 * n),
            adjacent: Vec::with_capacity(2 * n),
            conflict_head: Vec::with_capacity(2 * n),
            stamp: Vec::with_capacity(2 * n),
            free: Vec::new(),
            conflict_next: vec![NONE; n],
            conflict_face: vec![NONE; n],
            new_face_at: vec![NONE; n],
            step: 0,
            stack: Vec::new(),
            visible: Vec::new(),
            new_faces: Vec::new(),
        };
        let (a, mut b, mut c, d) = (order[0], order[1], order[2], order[3]);
        // The perturbed orientation is never zero for distinct indices.
        if orient3d_sos(points, [a, b, c, d]) == Orientation::Positive {
            core::mem::swap(&mut b, &mut c);
        }
        // Now d lies below (a, b, c); the other faces are even permutations
        // of (a, b, c, d) with the left-out vertex last.
        for face in [[a, b, c], [a, c, d], [a, d, b], [b, d, c]] {
            hull.add_face(face);
        }
        for f in 0..4 {
            for k in 0..3 {
                let (u, v) = (hull.faces[f][k], hull.faces[f][(k + 1) % 3]);
                let g = (0..4)
                    .find(|&g| g != f && hull.edge_slot(g as u32, v, u).is_some())
                    .expect("every tetrahedron edge has a twin");
                hull.adjacent[f][k] = g as u32;
            }
        }
        for &q in &order[4..] {
            if let Some(f) = (0..4).find(|&f| hull.sees(f, q)) {
                hull.attach(q, f);
            }
        }
        hull
    }

    /// Whether `q` lies strictly above the plane of face `f` (perturbed).
    #[inline]
    fn sees(&self, f: u32, q: u32) -> bool {
        let [a, b, c] = self.faces[f as usize];
        orient3d_sos(self.points, [a, b, c, q]) == Orientation::Positive
    }

    /// The edge of face `f` running from `u` to `v`, if any.
    fn edge_slot(&self, f: u32, u: u32, v: u32) -> Option<usize> {
        let face = self.faces[f as usize];
        (0..3).find(|&k| face[k] == u && face[(k + 1) % 3] == v)
    }

    /// A new face with the given vertices and no links, in a free slot.
    fn add_face(&mut self, vertices: [u32; 3]) -> u32 {
        if let Some(f) = self.free.pop() {
            let i = f as usize;
            self.faces[i] = vertices;
            self.adjacent[i] = [NONE; 3];
            self.conflict_head[i] = NONE;
            f
        } else {
            let f = self.faces.len() as u32;
            self.faces.push(vertices);
            self.adjacent.push([NONE; 3]);
            self.conflict_head.push(NONE);
            self.stamp.push(0);
            f
        }
    }

    /// Attaches the pending point `q` to face `f`.
    fn attach(&mut self, q: u32, f: u32) {
        self.conflict_next[q as usize] = self.conflict_head[f as usize];
        self.conflict_head[f as usize] = q;
        self.conflict_face[q as usize] = f;
    }

    /// Adds point `p` to the hull, unless it has been dropped.
    fn insert(&mut self, p: u32) {
        let start = self.conflict_face[p as usize];
        if start == NONE {
            return;
        }
        self.step += 1;
        let seen = 2 * self.step;
        let hidden = seen + 1;
        let mut stack = core::mem::take(&mut self.stack);
        let mut visible = core::mem::take(&mut self.visible);
        let mut new_faces = core::mem::take(&mut self.new_faces);
        stack.clear();
        visible.clear();
        new_faces.clear();

        // The visible region, and a new face on each horizon edge.
        self.stamp[start as usize] = seen;
        stack.push(start);
        while let Some(f) = stack.pop() {
            visible.push(f);
            for k in 0..3 {
                let g = self.adjacent[f as usize][k];
                let mark = self.stamp[g as usize];
                if mark == seen {
                    continue;
                }
                if mark != hidden && self.sees(g, p) {
                    self.stamp[g as usize] = seen;
                    stack.push(g);
                    continue;
                }
                self.stamp[g as usize] = hidden;
                let face = self.faces[f as usize];
                let (u, v) = (face[k], face[(k + 1) % 3]);
                let nf = self.add_face([u, v, p]);
                self.adjacent[nf as usize][0] = g;
                let j = self
                    .edge_slot(g, v, u)
                    .expect("the hidden face shares the horizon edge");
                self.adjacent[g as usize][j] = nf;
                debug_assert_eq!(self.new_face_at[u as usize], NONE, "horizon is not simple");
                self.new_face_at[u as usize] = nf;
                new_faces.push(nf);
            }
        }

        // Link the new faces around p: the face (u, v, p) meets the face
        // starting at v across its edge 1 (v, p), which is that face's
        // edge 2 (p, v).
        for &nf in &new_faces {
            let v = self.faces[nf as usize][1];
            let next = self.new_face_at[v as usize];
            debug_assert_ne!(next, NONE, "horizon is not closed");
            self.adjacent[nf as usize][1] = next;
            self.adjacent[next as usize][2] = nf;
        }
        for &nf in &new_faces {
            self.new_face_at[self.faces[nf as usize][0] as usize] = NONE;
        }

        // Reattach the pending points of the visible faces and free them.
        for &f in &visible {
            let mut q = self.conflict_head[f as usize];
            while q != NONE {
                let next = self.conflict_next[q as usize];
                self.conflict_face[q as usize] = NONE;
                if q != p {
                    if let Some(&nf) = new_faces.iter().find(|&&nf| self.sees(nf, q)) {
                        self.attach(q, nf);
                    }
                }
                q = next;
            }
            let i = f as usize;
            self.conflict_head[i] = NONE;
            self.faces[i] = [NONE; 3];
            self.adjacent[i] = [NONE; 3];
            self.free.push(f);
        }
        self.conflict_face[p as usize] = NONE;

        self.stack = stack;
        self.visible = visible;
        self.new_faces = new_faces;
    }

    /// The live faces in canonical order: each rotated so that its smallest
    /// index comes first (a rotation keeps the orientation), sorted, with
    /// the links renumbered.
    fn into_triangulation(self) -> Triangulation {
        let mut live: Vec<([u32; 3], u32, usize)> = Vec::with_capacity(self.faces.len());
        for (f, face) in self.faces.iter().enumerate() {
            if face[0] != NONE {
                let r = (0..3).min_by_key(|&k| face[k]).unwrap_or(0);
                let rotated = [face[r], face[(r + 1) % 3], face[(r + 2) % 3]];
                live.push((rotated, f as u32, r));
            }
        }
        // The keys are distinct, so an unstable sort is deterministic.
        live.sort_unstable_by_key(|entry| entry.0);
        let mut renumber = vec![NONE; self.faces.len()];
        for (t, entry) in live.iter().enumerate() {
            renumber[entry.1 as usize] = t as u32;
        }
        let triangles = live.iter().map(|entry| entry.0).collect();
        let neighbors = live
            .iter()
            .map(|&(_, f, r)| {
                let links = self.adjacent[f as usize];
                [0, 1, 2].map(|k| renumber[links[(k + r) % 3] as usize])
            })
            .collect();
        Triangulation {
            triangles,
            neighbors,
        }
    }
}

// The tests use std's trigonometric functions as independent references.
#[cfg(test)]
#[allow(clippy::disallowed_methods)]
mod tests;

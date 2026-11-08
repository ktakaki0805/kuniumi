//! The Voronoi cells of the grid, as the dual of the spherical Delaunay
//! triangulation (design document, 球面グリッドとデータ構造; D-14).
//!
//! The Voronoi cell of a centre is the set of points of the sphere closer to
//! it than to any other centre. Its vertices are the circumcentres of the
//! Delaunay triangles around the centre, and its edges cross the Delaunay
//! edges: the cell of `v` has one edge per neighbour `n`, on the great
//! circle equidistant from `v` and `n`, between the circumcentres of the two
//! triangles that share the edge `(v, n)`.
//!
//! # Construction
//!
//! [`dual`] takes the triangles of [`Triangulation`] (counter-clockwise from
//! outside, smallest index first, canonical order) and computes:
//!
//! 1. per triangle `(a, b, c)`, the circumcentre on the sphere,
//!    `c_t = w / |w|` with `w = (b − a) × (c − a)`, the unit normal of the
//!    triangle's plane: every point of that plane is equidistant from `a`,
//!    `b` and `c`, and so is its intersection with the sphere. `w` points
//!    away from the centre of the sphere, which lies strictly below every
//!    triangle (checked exactly by [`triangulate`](crate::delaunay::triangulate)),
//!    so `c_t` is the circumcentre on the triangle's side, the Voronoi vertex.
//!    Three distinct points of a sphere are never collinear, so `w ≠ 0`.
//! 2. per cell `v`, the ring of neighbours, by walking the fan of triangles
//!    around `v` without any angle or trigonometry: start at the first
//!    triangle of the canonical order that contains `v`; in the current
//!    triangle, rotated to `(v, b, c)`, append `b` and step to the triangle
//!    across the edge `(c, v)`; after `deg(v)` steps (the number of
//!    triangles containing `v`) the walk is back at the start. This gives
//!    the neighbours `n_0 … n_{d−1}` and the triangles
//!    `t_m = (v, n_m, n_{m+1})`, both counter-clockwise seen from outside.
//! 3. per CSR entry `(v, n_m)`: the edge length, the great-circle distance
//!    between the circumcentres of `t_{m−1}` and `t_m` (the two triangles
//!    sharing the edge), and the centre distance, the great-circle distance
//!    between `v` and `n_m` ([`great_circle_rad`]). The arguments are passed
//!    in a fixed order (the smaller triangle or cell index first), so the
//!    entry of `(n_m, v)` is the same computation and the same bits.
//! 4. per cell, the area: the sum over its edges of the solid angles of the
//!    triangles from `p̂ = p/|p|` (the input is a unit vector only within
//!    1e-12, and the formula assumes unit vectors) to the edge, each split
//!    at the edge's midpoint `e = (c₁ + c₂)/|c₁ + c₂|`: for the edge to
//!    `n_m`, `Ω(p̂, c_{t_{m−1}}, e) + Ω(p̂, e, c_{t_m})`. Floating-point
//!    addition is commutative, so both cells of the edge get the same bits
//!    for `e`. The split keeps every triangle
//!    well conditioned: an edge can be nearly a semicircle (a cluster of
//!    centres and a centre nearly opposite it), and for two vertices
//!    `π − δ` apart both terms of the formula below are of order `δ` while
//!    their rounding errors are of order `ε`, so the angle would be off by
//!    about `2ε/δ` (review R4-1 of M1-06: 3e-12 to 7e-11 of 4π for
//!    `δ` = 1e-5 to 1e-7). The solid angle is by A. van Oosterom and
//!    J. Strackee (1983), "The solid angle of a plane triangle", *IEEE
//!    Trans. Biomed. Eng.* BME-30(2), 125–126, for unit vectors `a`, `b`,
//!    `c`:
//!
//!    `Ω = 2·atan2(a · (b × c), 1 + a·b + b·c + c·a)`,
//!
//!    the paper's formula with the signed triple product (the determinant
//!    `[a b c]`), so that `Ω` is negative for a clockwise triangle. A
//!    Voronoi cell is convex and contains its centre, so on the sphere every
//!    triangle of the fan is counter-clockwise (or, for co-circular centres
//!    with coincident circumcentres, empty). The input points, however, lie
//!    on the sphere only to within rounding (about 1e-16) or the unit
//!    tolerance (up to 5e-13 in length). For two centres `p`, `q` a
//!    distance `s` apart with `|q|² − |p|² = 2δ`, the plane equidistant from
//!    both, on which their Voronoi edge lies, passes `δ/s` from the centre of
//!    the sphere, which moves the edge by about `δ/s`; when that exceeds
//!    `s/2`, one of the two centres lies outside the fan of its Voronoi
//!    vertices, even with exactly computed circumcentres (their rounding
//!    adds to it). The signed solid angles still add up to the area of the
//!    polygon of the vertices, while their absolute values (the ticket's
//!    `|a · (b × c)|`) would count a part twice: with a centre 1e-7 rad from
//!    another and 4e-13 off unit length, the sum of the areas was 3.5e-8
//!    relative from 4π (review R1-1 of M1-06). A centre that is a hull
//!    vertex only barely can get a cell of zero or negative computed area;
//!    [`SphereGrid::from_points`](crate::SphereGrid::from_points) refuses
//!    such a set (review R3-1).
//!
//! # Why the areas sum to 4π
//!
//! Two neighbouring cells use the same two circumcentre values and the same
//! midpoint for their common edge, so the fans of all cells, counted with
//! the signs of their solid angles, cover the sphere once, with the rounded
//! circumcentres as well as with exact ones. The sum of the areas differs
//! from 4π only by the rounding of each solid angle and of the summation.
//! Measured for 1,000 to 650,000 cells (ticket M1-06, Result): within
//! 1.5e-16 relative with a compensated sum, within 1.5e-14 with a plain sum
//! in index order. A wrong ring or a circumcentre on the wrong side changes
//! it at once.
//!
//! # Determinism
//!
//! Everything is `+ − × ÷`, `sqrt` and [`det::atan2`], in a fixed order: the
//! triangles and the cells in index order, each ring from its fixed start.
//! The result is a function of the points and the triangulation, which
//! itself depends only on the points.

use kuniumi_math::{Vec3, det};

use crate::coords::great_circle_rad;
use crate::delaunay::Triangulation;

/// The arrays of a [`SphereGrid`](crate::SphereGrid) other than the
/// positions; see the documentation of its accessors.
pub(crate) struct Dual {
    /// CSR offsets, `n + 1` entries.
    pub(crate) neighbor_offsets: Vec<u32>,
    /// The rings, counter-clockwise from outside.
    pub(crate) neighbors: Vec<u32>,
    /// Per cell, the Voronoi area in steradians.
    pub(crate) areas_sr: Vec<f64>,
    /// Per CSR entry, the length of the Voronoi edge in radians.
    pub(crate) edge_lengths_rad: Vec<f64>,
    /// Per CSR entry, the distance between the two centres in radians.
    pub(crate) centre_distances_rad: Vec<f64>,
}

/// The Voronoi dual of `triangulation`, a valid triangulation of `points`
/// as returned by [`triangulate`](crate::delaunay::triangulate); see the
/// module documentation.
pub(crate) fn dual(points: &[Vec3], triangulation: &Triangulation) -> Dual {
    let n = points.len();
    let triangles = triangulation.triangles();
    let adjacent = triangulation.neighbors();
    let circumcentres: Vec<Vec3> = triangles
        .iter()
        .map(|&[a, b, c]| circumcentre(points[a as usize], points[b as usize], points[c as usize]))
        .collect();

    // The first triangle of each vertex and the number of triangles around it.
    let mut first = vec![u32::MAX; n];
    let mut degree = vec![0u32; n];
    for (t, triangle) in triangles.iter().enumerate() {
        for &v in triangle {
            let v = v as usize;
            degree[v] += 1;
            if first[v] == u32::MAX {
                first[v] = t as u32;
            }
        }
    }
    let mut neighbor_offsets = Vec::with_capacity(n + 1);
    let mut total = 0u32;
    neighbor_offsets.push(0);
    for &d in &degree {
        total += d;
        neighbor_offsets.push(total);
    }
    let entries = total as usize;

    let mut neighbors = Vec::with_capacity(entries);
    let mut edge_lengths_rad = Vec::with_capacity(entries);
    let mut centre_distances_rad = Vec::with_capacity(entries);
    let mut areas_sr = Vec::with_capacity(n);
    let mut fan: Vec<u32> = Vec::new();
    for v in 0..n as u32 {
        // Walk the fan: in t = (v, b, c), b is the next neighbour and the
        // next triangle lies across the edge (c, v).
        let start = neighbors.len();
        fan.clear();
        let mut t = first[v as usize];
        for _ in 0..degree[v as usize] {
            let triangle = triangles[t as usize];
            let k = slot(triangle, v);
            neighbors.push(triangle[(k + 1) % 3]);
            fan.push(t);
            t = adjacent[t as usize][(k + 2) % 3];
        }
        debug_assert_eq!(t, first[v as usize], "the fan of {v} does not close");

        let ring = &neighbors[start..];
        let d = fan.len();
        let p = points[v as usize];
        let p_unit = p / p.length();
        let mut area_sr = 0.0;
        for m in 0..d {
            // The edge to n_m lies between t_{m−1} and t_m.
            let (before, after) = (fan[(m + d - 1) % d], fan[m]);
            let (s, e) = (before.min(after), before.max(after));
            let (c_lo, c_hi) = (circumcentres[s as usize], circumcentres[e as usize]);
            edge_lengths_rad.push(great_circle_rad(c_lo, c_hi));
            let j = ring[m];
            let (lo, hi) = (v.min(j), v.max(j));
            centre_distances_rad.push(great_circle_rad(points[lo as usize], points[hi as usize]));
            // The fan triangle over this edge, split at the edge's midpoint
            // (the same bits for both cells of the edge: `+` commutes).
            let sum = c_lo + c_hi;
            let middle = sum / sum.length();
            area_sr += solid_angle_sr(p_unit, circumcentres[before as usize], middle);
            area_sr += solid_angle_sr(p_unit, middle, circumcentres[after as usize]);
        }
        areas_sr.push(area_sr);
    }
    Dual {
        neighbor_offsets,
        neighbors,
        areas_sr,
        edge_lengths_rad,
        centre_distances_rad,
    }
}

/// The position of `v` in `triangle`.
fn slot(triangle: [u32; 3], v: u32) -> usize {
    if triangle[0] == v {
        0
    } else if triangle[1] == v {
        1
    } else {
        debug_assert_eq!(triangle[2], v, "{v} is not a vertex of {triangle:?}");
        2
    }
}

/// The circumcentre on the unit sphere of the triangle `(a, b, c)`,
/// counter-clockwise seen from outside: `w / |w|` with
/// `w = (b − a) × (c − a)`.
fn circumcentre(a: Vec3, b: Vec3, c: Vec3) -> Vec3 {
    let w = (b - a).cross(c - a);
    w / w.length()
}

/// The signed solid angle of the triangle of unit vectors `(a, b, c)` seen
/// from the centre, in steradians, positive when the triangle is
/// counter-clockwise seen from outside (Van Oosterom and Strackee 1983):
/// `2·atan2(a · (b × c), 1 + a·b + b·c + c·a)`.
fn solid_angle_sr(a: Vec3, b: Vec3, c: Vec3) -> f64 {
    let triple = a.dot(b.cross(c));
    2.0 * det::atan2(triple, 1.0 + a.dot(b) + b.dot(c) + c.dot(a))
}

//! The globe as a triangle mesh: every Voronoi cell a fan of triangles from
//! its centre to its Voronoi vertices, so that each cell can take its own
//! colour.
//!
//! Cell `v` with the ring of neighbours `n_0 … n_{k−1}` (counter-clockwise
//! seen from outside) gets `1 + k` vertices: its centre, at the cell's
//! height, and one vertex per Delaunay triangle `t_m = (v, n_m, n_{m+1})`
//! around it, the triangle's circumcentre (the Voronoi vertex, as in
//! `kuniumi_sphere`'s `voronoi` module: `w/|w|` with `w = (b − a) × (c − a)`)
//! at the mean height of its three cells. Each Voronoi vertex belongs to
//! three cells; it is computed in `f64` on the triangle rotated so that its
//! smallest cell index comes first, so the three copies have the same bits
//! and the surface has no cracks. The fan's triangles `(v, t_m, t_{m+1})` are
//! counter-clockwise seen from outside.
//!
//! A grid of `n` cells has `6n − 12` neighbour entries (Euler's formula for a
//! triangulated sphere), so the mesh has `7n − 12` vertices and `6n − 12`
//! triangles: 4.55 million and 3.9 million for 650,000 cells.

use kuniumi_io::Rgb8;
use kuniumi_math::vec3::Vec3;
use kuniumi_sphere::SphereGrid;

/// The globe's vertices and triangles, ready for the GPU.
#[derive(Clone, Debug, PartialEq)]
pub struct GlobeMesh {
    /// Per vertex: the unit direction (x, y, z) and the height above the sea
    /// level (m). Cell `i` owns the vertices
    /// `cell_vertex_offsets[i]..cell_vertex_offsets[i + 1]`, its centre
    /// first.
    pub vertices: Vec<[f32; 4]>,
    /// Three vertex indices per triangle, counter-clockwise seen from
    /// outside.
    pub indices: Vec<u32>,
    /// `n + 1` offsets into `vertices`, one range per cell.
    pub cell_vertex_offsets: Vec<u32>,
    /// The largest height above the sea level of any vertex (m), at least 0.
    pub max_height_m: f64,
}

impl GlobeMesh {
    /// The mesh of `grid` with the heights above the sea level
    /// `height_above_sea_m` (m, one per cell; a non-finite height is drawn
    /// at 0).
    ///
    /// # Panics
    ///
    /// If `height_above_sea_m` does not hold one value per cell, or the mesh
    /// would need more than `u32::MAX` vertices.
    pub fn new(grid: &SphereGrid, height_above_sea_m: &[f64]) -> Self {
        let n = grid.len();
        assert_eq!(
            height_above_sea_m.len(),
            n,
            "the mesh needs one height per cell"
        );
        let heights: Vec<f64> = height_above_sea_m
            .iter()
            .map(|&h| if h.is_finite() { h } else { 0.0 })
            .collect();
        let entries: usize = (0..n as u32).map(|v| grid.neighbors(v).len()).sum();
        let vertex_count = n + entries;
        assert!(
            u32::try_from(vertex_count).is_ok(),
            "the mesh of {n} cells needs {vertex_count} vertices, more than u32 indices reach"
        );
        let mut vertices = Vec::with_capacity(vertex_count);
        let mut indices = Vec::with_capacity(3 * entries);
        let mut cell_vertex_offsets = Vec::with_capacity(n + 1);
        for v in 0..n as u32 {
            let ring = grid.neighbors(v);
            let k = ring.len() as u32;
            let centre = vertices.len() as u32;
            cell_vertex_offsets.push(centre);
            vertices.push(vertex(grid.position(v), heights[v as usize]));
            for (m, &b) in ring.iter().enumerate() {
                let c = ring[(m + 1) % ring.len()];
                let (p, h) = voronoi_vertex(grid, &heights, v, b, c);
                vertices.push(vertex(p, h));
            }
            for m in 0..k {
                indices.extend([centre, centre + 1 + m, centre + 1 + (m + 1) % k]);
            }
        }
        cell_vertex_offsets.push(vertices.len() as u32);
        let max_height_m = heights.iter().copied().fold(0.0, f64::max);
        Self {
            vertices,
            indices,
            cell_vertex_offsets,
            max_height_m,
        }
    }

    /// The vertices of `cell`: its centre first, then its Voronoi vertices.
    ///
    /// # Panics
    ///
    /// If `cell` is not a cell of the mesh.
    pub fn cell_vertices(&self, cell: u32) -> core::ops::Range<u32> {
        let v = cell as usize;
        self.cell_vertex_offsets[v]..self.cell_vertex_offsets[v + 1]
    }

    /// The colour of every vertex: each cell's colour from `cell_colours`
    /// over its vertex range, opaque, as `[r, g, b, 255]` (sRGB).
    ///
    /// # Panics
    ///
    /// If `cell_colours` does not hold one colour per cell.
    pub fn vertex_colours(&self, cell_colours: &[Rgb8]) -> Vec<[u8; 4]> {
        assert_eq!(
            cell_colours.len() + 1,
            self.cell_vertex_offsets.len(),
            "the mesh needs one colour per cell"
        );
        let mut colours = Vec::with_capacity(self.vertices.len());
        for (range, &[r, g, b]) in self.cell_vertex_offsets.windows(2).zip(cell_colours) {
            colours.resize(range[1] as usize, [r, g, b, 255]);
        }
        colours
    }
}

/// The Voronoi vertex of the Delaunay triangle `(a, b, c)` (cell indices,
/// counter-clockwise seen from outside) and its height, the mean of the
/// three cells' `heights`: both computed on the triangle rotated so that its
/// smallest index comes first, so every cell of the triangle gets the same
/// bits.
pub fn voronoi_vertex(grid: &SphereGrid, heights: &[f64], a: u32, b: u32, c: u32) -> (Vec3, f64) {
    // Rotating keeps the triangle counter-clockwise.
    let [a, b, c] = if a < b && a < c {
        [a, b, c]
    } else if b < c {
        [b, c, a]
    } else {
        [c, a, b]
    };
    let (pa, pb, pc) = (grid.position(a), grid.position(b), grid.position(c));
    let w = (pb - pa).cross(pc - pa);
    let height = (heights[a as usize] + heights[b as usize] + heights[c as usize]) / 3.0;
    (w / w.length(), height)
}

/// A vertex at the unit direction `p` and the height `h` (m), rounded to
/// `f32`.
fn vertex(p: Vec3, h: f64) -> [f32; 4] {
    [p.x as f32, p.y as f32, p.z as f32, h as f32]
}

// The tests use std's trigonometric functions as independent references.
#[cfg(test)]
#[allow(clippy::disallowed_methods)]
mod tests {
    use std::collections::BTreeMap;
    use std::f64::consts::PI;

    use super::*;

    fn grids() -> Vec<SphereGrid> {
        vec![
            SphereGrid::octahedron(),
            SphereGrid::fibonacci(1_000, 1).unwrap(),
            SphereGrid::fibonacci(20_000, 7).unwrap(),
        ]
    }

    /// Heights that differ between cells, some below the sea level, not
    /// exact in binary, so that the order of a sum shows in its bits.
    fn heights(grid: &SphereGrid) -> Vec<f64> {
        (0..grid.len())
            .map(|i| (i % 17) as f64 * 731.3 - 5_000.1)
            .collect()
    }

    fn direction(v: [f32; 4]) -> Vec3 {
        Vec3::new(f64::from(v[0]), f64::from(v[1]), f64::from(v[2]))
    }

    #[test]
    fn counts_follow_euler() {
        for grid in grids() {
            let n = grid.len();
            let mesh = GlobeMesh::new(&grid, &heights(&grid));
            assert_eq!(mesh.vertices.len(), 7 * n - 12);
            assert_eq!(mesh.indices.len(), 3 * (6 * n - 12));
            assert_eq!(mesh.cell_vertex_offsets.len(), n + 1);
            assert_eq!(mesh.cell_vertex_offsets[0], 0);
            assert_eq!(mesh.cell_vertex_offsets[n] as usize, mesh.vertices.len());
            for i in 0..n {
                let range = mesh.cell_vertex_offsets[i]..mesh.cell_vertex_offsets[i + 1];
                assert_eq!(range.len(), 1 + grid.neighbors(i as u32).len());
            }
            assert!(
                mesh.indices
                    .iter()
                    .all(|&k| (k as usize) < mesh.vertices.len())
            );
            let top = heights(&grid).into_iter().fold(0.0, f64::max);
            assert_eq!(mesh.max_height_m.to_bits(), top.to_bits());
        }
    }

    #[test]
    fn voronoi_vertices_are_equidistant_and_shared_bit_for_bit() {
        for grid in grids() {
            let h = heights(&grid);
            let mesh = GlobeMesh::new(&grid, &h);
            let mut copies: BTreeMap<[u32; 3], Vec<[u32; 4]>> = BTreeMap::new();
            for v in 0..grid.len() as u32 {
                let ring = grid.neighbors(v);
                let first = mesh.cell_vertex_offsets[v as usize] as usize;
                // The centre, at the cell's height.
                let p = grid.position(v);
                let centre = mesh.vertices[first];
                assert_eq!(
                    centre.map(f32::to_bits),
                    [p.x as f32, p.y as f32, p.z as f32, h[v as usize] as f32].map(f32::to_bits)
                );
                for (m, &b) in ring.iter().enumerate() {
                    let c = ring[(m + 1) % ring.len()];
                    let (q, height) = voronoi_vertex(&grid, &h, v, b, c);
                    // The same bits from each cell of the triangle, in f64.
                    for (x, y, z) in [(b, c, v), (c, v, b)] {
                        let (q2, height2) = voronoi_vertex(&grid, &h, x, y, z);
                        assert_eq!(
                            [q2.x, q2.y, q2.z, height2].map(f64::to_bits),
                            [q.x, q.y, q.z, height].map(f64::to_bits)
                        );
                    }
                    let (pa, pb, pc) = (grid.position(v), grid.position(b), grid.position(c));
                    assert!((q.length() - 1.0).abs() < 1e-15);
                    assert!((q.dot(pa) - q.dot(pb)).abs() < 1e-12);
                    assert!((q.dot(pa) - q.dot(pc)).abs() < 1e-12);
                    // On the triangle's side of the sphere.
                    assert!(q.dot(pa) > 0.0);
                    let mean = (h[v as usize] + h[b as usize] + h[c as usize]) / 3.0;
                    assert!((height - mean).abs() <= 1e-12 * mean.abs().max(1.0));
                    let vertex = mesh.vertices[first + 1 + m];
                    assert_eq!(
                        vertex.map(f32::to_bits),
                        [q.x as f32, q.y as f32, q.z as f32, height as f32].map(f32::to_bits)
                    );
                    let mut key = [v, b, c];
                    key.sort_unstable();
                    copies
                        .entry(key)
                        .or_default()
                        .push(vertex.map(f32::to_bits));
                }
            }
            // Every Delaunay triangle (2n − 4 of them) appears in its three
            // cells with the same bits.
            assert_eq!(copies.len(), 2 * grid.len() - 4);
            for (key, list) in &copies {
                assert_eq!(list.len(), 3, "{key:?}");
                assert!(
                    list.iter().all(|bits| *bits == list[0]),
                    "{key:?}: {list:?}"
                );
            }
        }
    }

    /// The signed solid angle of the triangle `(a, b, c)` seen from the
    /// centre (Van Oosterom and Strackee 1983), for any non-zero vectors.
    fn solid_angle_sr(a: Vec3, b: Vec3, c: Vec3) -> f64 {
        let (a, b, c) = (a / a.length(), b / b.length(), c / c.length());
        2.0 * a
            .dot(b.cross(c))
            .atan2(1.0 + a.dot(b) + b.dot(c) + c.dot(a))
    }

    #[test]
    fn triangles_are_counter_clockwise_and_cover_the_sphere_once() {
        for grid in grids() {
            let mesh = GlobeMesh::new(&grid, &heights(&grid));
            let mut total_sr = 0.0;
            for t in mesh.indices.chunks_exact(3) {
                let [a, b, c] = [t[0], t[1], t[2]].map(|k| direction(mesh.vertices[k as usize]));
                assert!(a.dot(b.cross(c)) > 0.0, "{t:?}");
                total_sr += solid_angle_sr(a, b, c);
            }
            // The copies of each vertex are equal, so the solid angles add up
            // to 4π up to the rounding of the sum alone.
            assert!(
                (total_sr - 4.0 * PI).abs() < 1e-9,
                "{}: {total_sr}",
                grid.len()
            );
        }
    }

    #[test]
    fn octahedron_by_hand() {
        // Cells ±X, ±Y, ±Z; every Voronoi vertex is a face centre
        // (±1, ±1, ±1)/√3 at the mean of its three cells' heights.
        let grid = SphereGrid::octahedron();
        let h = [100.0, -200.0, 300.0, -400.0, 500.0, -600.0];
        let mesh = GlobeMesh::new(&grid, &h);
        for v in 0..6u32 {
            let ring = grid.neighbors(v);
            assert_eq!(ring.len(), 4);
            let first = mesh.cell_vertex_offsets[v as usize] as usize;
            for (m, &b) in ring.iter().enumerate() {
                let c = ring[(m + 1) % 4];
                let face = grid.position(v) + grid.position(b) + grid.position(c);
                let expected = face / 3.0_f64.sqrt();
                let vertex = mesh.vertices[first + 1 + m];
                assert!((direction(vertex) - expected).length() < 1e-7, "{v} {m}");
                let mean = (h[v as usize] + h[b as usize] + h[c as usize]) / 3.0;
                assert_eq!(vertex[3].to_bits(), (mean as f32).to_bits());
            }
        }
        assert_eq!(mesh.max_height_m.to_bits(), 500.0_f64.to_bits());
        // A world below the sea level everywhere still has a top of 0.
        assert_eq!(
            (GlobeMesh::new(&grid, &[-1.0; 6]).max_height_m).to_bits(),
            0.0_f64.to_bits()
        );
        // Non-finite heights are drawn at 0.
        let odd = GlobeMesh::new(&grid, &[f64::NAN, 1.0, 1.0, 1.0, 1.0, f64::INFINITY]);
        assert_eq!(odd.vertices[0][3].to_bits(), 0.0_f32.to_bits());
        assert!(odd.vertices.iter().all(|v| v[3].is_finite()));
    }

    #[test]
    fn vertex_colours_repeat_each_cell_colour() {
        let grid = SphereGrid::fibonacci(1_000, 1).unwrap();
        let mesh = GlobeMesh::new(&grid, &vec![0.0; grid.len()]);
        let colours: Vec<Rgb8> = (0..grid.len())
            .map(|i| [i as u8, (i / 256) as u8, 7])
            .collect();
        let per_vertex = mesh.vertex_colours(&colours);
        assert_eq!(per_vertex.len(), mesh.vertices.len());
        for (i, colour) in colours.iter().enumerate() {
            let range = mesh.cell_vertices(i as u32);
            assert_eq!(range.len(), 1 + grid.neighbors(i as u32).len());
            for vertex_colour in &per_vertex[range.start as usize..range.end as usize] {
                assert_eq!(*vertex_colour, [colour[0], colour[1], colour[2], 255]);
            }
        }
    }

    #[test]
    #[should_panic(expected = "one height per cell")]
    fn rejects_heights_of_another_grid() {
        let _ = GlobeMesh::new(&SphereGrid::octahedron(), &[0.0; 5]);
    }
}

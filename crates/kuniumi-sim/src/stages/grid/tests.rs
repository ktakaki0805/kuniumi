//! Tests of the Grid stage:
//!
//! - `streams_are_under_the_stage_name`;
//! - `builds_the_fibonacci_grid_of_the_settings`;
//! - `golden_fingerprints`.

use kuniumi_math::rng::fnv1a64;
use kuniumi_sphere::{SphereGrid, delaunay, points};
use kuniumi_world::Settings;

use super::run;
use crate::StageId;

#[test]
fn streams_are_under_the_stage_name() {
    let prefix = format!("{}/", StageId::Grid.name());
    assert_eq!(prefix, "grid/");
    assert!(points::STREAM.starts_with(&prefix), "{}", points::STREAM);
    assert!(
        delaunay::STREAM.starts_with(&prefix),
        "{}",
        delaunay::STREAM
    );
}

#[test]
fn builds_the_fibonacci_grid_of_the_settings() {
    for (cells, seed) in [(1_000, 1), (1_000, 2), (2_500, 7)] {
        let grid = run(&Settings { cells }, seed).expect("valid settings");
        assert_eq!(grid.len(), cells as usize);
        assert_eq!(
            grid,
            SphereGrid::fibonacci(cells, seed).expect("valid grid"),
            "{cells} cells, seed {seed}"
        );
    }
    // The seed reaches the points.
    let one = run(&Settings { cells: 1_000 }, 1).expect("valid settings");
    let two = run(&Settings { cells: 1_000 }, 2).expect("valid settings");
    assert_ne!(one.positions(), two.positions());
}

// ---------------------------------------------------------------------------
// Golden fingerprints

/// FNV-1a of the positions (x, y, z) and then the areas, little-endian: the
/// bytes of `fingerprint_positions_areas` in the tests of
/// `kuniumi_sphere::SphereGrid` (ticket M1-06).
fn fingerprint_positions_areas(g: &SphereGrid) -> u64 {
    let mut bytes = Vec::new();
    for p in g.positions() {
        for x in [p.x, p.y, p.z] {
            bytes.extend_from_slice(&x.to_le_bytes());
        }
    }
    for i in 0..g.len() as u32 {
        bytes.extend_from_slice(&g.area_sr(i).to_le_bytes());
    }
    fnv1a64(&bytes)
}

/// FNV-1a of the CSR offsets, the neighbours, the edge lengths and the
/// centre distances, little-endian: the bytes of
/// `fingerprint_neighbours_edges` in the tests of `SphereGrid` (ticket
/// M1-06), rebuilt from the public accessors.
fn fingerprint_neighbours_edges(g: &SphereGrid) -> u64 {
    let cells = 0..g.len() as u32;
    let mut bytes = Vec::new();
    let mut offset = 0u32;
    bytes.extend_from_slice(&offset.to_le_bytes());
    for i in cells.clone() {
        offset += g.neighbors(i).len() as u32;
        bytes.extend_from_slice(&offset.to_le_bytes());
    }
    for i in cells.clone() {
        for &j in g.neighbors(i) {
            bytes.extend_from_slice(&j.to_le_bytes());
        }
    }
    for i in cells.clone() {
        for x in g.edge_lengths_rad(i) {
            bytes.extend_from_slice(&x.to_le_bytes());
        }
    }
    for i in cells {
        for x in g.centre_distances_rad(i) {
            bytes.extend_from_slice(&x.to_le_bytes());
        }
    }
    fnv1a64(&bytes)
}

/// The fingerprints recorded for `SphereGrid::fibonacci(10_000, 1)` by
/// ticket M1-06 (`RECORDED_POSITIONS_AREAS` and `RECORDED_NEIGHBOURS_EDGES`
/// in the tests of `SphereGrid`), after its invariant and reference tests
/// had passed and its images had been inspected. The Grid stage must build
/// exactly that grid for 10,000 cells and seed 1; a mismatch means every
/// world changed (world-breaking).
const RECORDED_POSITIONS_AREAS: u64 = 0xE4CB_7D8D_8425_FDEF;
/// As [`RECORDED_POSITIONS_AREAS`], for the neighbours and the edge data.
const RECORDED_NEIGHBOURS_EDGES: u64 = 0x01B8_F9CE_F488_A7C7;

#[test]
fn golden_fingerprints() {
    let grid = run(&Settings { cells: 10_000 }, 1).expect("valid settings");
    let positions_areas = fingerprint_positions_areas(&grid);
    let neighbours_edges = fingerprint_neighbours_edges(&grid);
    assert_eq!(
        (positions_areas, neighbours_edges),
        (RECORDED_POSITIONS_AREAS, RECORDED_NEIGHBOURS_EDGES),
        "Grid stage fingerprints changed (world-breaking): \
         {positions_areas:#018X}, {neighbours_edges:#018X}"
    );
}

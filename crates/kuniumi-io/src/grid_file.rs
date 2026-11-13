//! The `.kgrid` binary format for a [`SphereGrid`]: every array of the
//! grid, bit for bit, so that a world directory holds its own cells
//! (ticket MV-01; [`crate::world_dir`] writes it as `grid.kgrid`).
//!
//! Layout (all integers little-endian, as in [`crate::field_file`]), for a
//! grid of `n` cells with `m` neighbour entries (each neighbour pair is
//! listed twice, once in each cell's ring; `m = 6n − 12` for a grid that
//! tiles the sphere):
//!
//! | bytes      | content                                                   |
//! | ---------- | --------------------------------------------------------- |
//! | 8          | magic `b"KUNIGRD\0"`                                      |
//! | 4          | format version (`u32`, currently 1)                       |
//! | 8          | cell count `n` (`u64`)                                    |
//! | 8          | neighbour entry count `m` (`u64`)                         |
//! | 24 × n     | cell centres: `x`, `y`, `z` of each cell (`f64` bits)     |
//! | 4 × (n + 1)| CSR offsets (`u32`): cell `i`'s neighbours are entries `offsets[i]..offsets[i + 1]` |
//! | 4 × m      | neighbours (`u32`), counter-clockwise seen from outside   |
//! | 8 × n      | cell areas on the unit sphere (sr, `f64` bits)            |
//! | 8 × m      | Voronoi edge lengths per entry (rad, `f64` bits)          |
//! | 8 × m      | centre distances per entry (rad, `f64` bits)              |
//!
//! These are the arrays of [`SphereGrid::from_parts`], in its order, and
//! [`read_grid`] passes them to it, so a grid read back has passed
//! [`SphereGrid::validate`] and equals the one written bit for bit. A grid
//! of 650,000 cells takes 101,399,792 bytes, 156 bytes per cell (the
//! measurements of ticket MV-01 are in [`crate::world_dir`]).

use std::io::{self, BufWriter, Read, Write};

use kuniumi_math::Vec3;
use kuniumi_sphere::{MAX_CELLS, SphereGrid};

use crate::binary::Reader;

/// File magic.
pub const MAGIC: [u8; 8] = *b"KUNIGRD\0";

/// Current format version.
pub const FORMAT_VERSION: u32 = 1;

/// The name of the format in messages.
const FORMAT: &str = ".kgrid";

/// Writes `grid` in `.kgrid` format (buffered).
pub fn write_grid(w: impl Write, grid: &SphereGrid) -> io::Result<()> {
    let mut w = BufWriter::new(w);
    let n = grid.len() as u32;
    let cells = || 0..n;
    let entries: usize = cells().map(|i| grid.neighbors(i).len()).sum();
    w.write_all(&MAGIC)?;
    w.write_all(&FORMAT_VERSION.to_le_bytes())?;
    w.write_all(&(grid.len() as u64).to_le_bytes())?;
    w.write_all(&(entries as u64).to_le_bytes())?;
    for p in grid.positions() {
        for c in [p.x, p.y, p.z] {
            w.write_all(&c.to_bits().to_le_bytes())?;
        }
    }
    // The offsets are the running sums of the ring sizes, as stored: the
    // first is 0 and cell i's ring is offsets[i]..offsets[i + 1]
    // (`SphereGrid::validate`).
    let mut offset = 0u32;
    w.write_all(&offset.to_le_bytes())?;
    for i in cells() {
        offset += grid.neighbors(i).len() as u32;
        w.write_all(&offset.to_le_bytes())?;
    }
    for i in cells() {
        for j in grid.neighbors(i) {
            w.write_all(&j.to_le_bytes())?;
        }
    }
    for i in cells() {
        w.write_all(&grid.area_sr(i).to_bits().to_le_bytes())?;
    }
    for i in cells() {
        for v in grid.edge_lengths_rad(i) {
            w.write_all(&v.to_bits().to_le_bytes())?;
        }
    }
    for i in cells() {
        for v in grid.centre_distances_rad(i) {
            w.write_all(&v.to_bits().to_le_bytes())?;
        }
    }
    w.flush()
}

/// Reads a grid written by [`write_grid`].
///
/// # Errors
///
/// [`io::ErrorKind::InvalidData`] if the input does not start with
/// [`MAGIC`], if the header gives more than [`MAX_CELLS`] cells or more
/// than `u32::MAX` neighbour entries, or if the arrays do not form a valid
/// grid ([`SphereGrid::from_parts`], with the [`GridError`]'s message);
/// [`io::ErrorKind::Unsupported`] for a format version other than
/// [`FORMAT_VERSION`]; [`io::ErrorKind::UnexpectedEof`] if the input ends
/// early, naming the array. The arrays are read in chunks, so a header that
/// claims more than the input holds allocates no more than about twice the
/// input plus 64 KiB. Bytes after the last array are not read.
///
/// [`GridError`]: kuniumi_sphere::GridError
pub fn read_grid(r: impl Read) -> io::Result<SphereGrid> {
    let mut r = Reader::new(r, FORMAT);
    r.header(&MAGIC, FORMAT_VERSION)?;
    let n = r.u64("cell count")?;
    if n > MAX_CELLS as u64 {
        return Err(r.invalid(format!(
            "the {FORMAT} file gives {n} cells, more than a grid can have ({MAX_CELLS})"
        )));
    }
    let m = r.u64("neighbour entry count")?;
    if m > u64::from(u32::MAX) {
        return Err(r.invalid(format!(
            "the {FORMAT} file gives {m} neighbour entries, more than u32 offsets can address \
             ({})",
            u32::MAX
        )));
    }
    // Both fit `usize` on 32-bit targets too: 3 × MAX_CELLS < 2³².
    let (n, m) = (n as usize, m as usize);
    let positions = r
        .f64s(3 * n, "cell centres")?
        .chunks_exact(3)
        .map(|c| Vec3::new(c[0], c[1], c[2]))
        .collect();
    let offsets = r.u32s(n + 1, "offsets")?;
    let neighbors = r.u32s(m, "neighbours")?;
    let areas_sr = r.f64s(n, "areas")?;
    let edge_lengths_rad = r.f64s(m, "edge lengths")?;
    let centre_distances_rad = r.f64s(m, "centre distances")?;
    SphereGrid::from_parts(
        positions,
        offsets,
        neighbors,
        areas_sr,
        edge_lengths_rad,
        centre_distances_rad,
    )
    .map_err(|e| r.invalid(format!("the {FORMAT} file does not hold a valid grid: {e}")))
}

#[cfg(test)]
mod tests {
    use core::f64::consts::PI;

    use super::*;

    fn encode(grid: &SphereGrid) -> Vec<u8> {
        let mut bytes = Vec::new();
        write_grid(&mut bytes, grid).unwrap();
        bytes
    }

    /// Every array of `grid`, as bits, through the public accessors.
    fn arrays(grid: &SphereGrid) -> [Vec<u64>; 5] {
        let n = grid.len() as u32;
        let cells = || 0..n;
        [
            grid.positions()
                .iter()
                .flat_map(|p| [p.x, p.y, p.z])
                .map(f64::to_bits)
                .collect(),
            cells()
                .flat_map(|i| grid.neighbors(i).iter().map(|&j| u64::from(j)))
                .chain(cells().map(|i| grid.neighbors(i).len() as u64))
                .collect(),
            cells().map(|i| grid.area_sr(i).to_bits()).collect(),
            cells()
                .flat_map(|i| grid.edge_lengths_rad(i).iter().map(|v| v.to_bits()))
                .collect(),
            cells()
                .flat_map(|i| grid.centre_distances_rad(i).iter().map(|v| v.to_bits()))
                .collect(),
        ]
    }

    #[test]
    fn round_trips_bit_exactly() {
        for grid in [
            SphereGrid::octahedron(),
            SphereGrid::fibonacci(1_000, 1).unwrap(),
            SphereGrid::fibonacci(20_000, 7).unwrap(),
        ] {
            let bytes = encode(&grid);
            let n = grid.len();
            let m = 6 * n - 12;
            assert_eq!(
                bytes.len(),
                28 + 24 * n + 4 * (n + 1) + 4 * m + 8 * n + 16 * m
            );
            let back = read_grid(bytes.as_slice()).unwrap();
            assert_eq!(arrays(&back), arrays(&grid), "{n} cells");
            assert_eq!(back, grid);
            // Writing it again gives the same bytes.
            assert!(encode(&back) == bytes, "{n} cells");
        }
    }

    /// The bytes of the octahedron grid, written out by hand from the
    /// layout in the module documentation and the arrays documented at
    /// `SphereGrid::octahedron`.
    fn octahedron_bytes() -> Vec<u8> {
        let mut b = Vec::new();
        b.extend(b"KUNIGRD\0");
        b.extend(1u32.to_le_bytes());
        b.extend(6u64.to_le_bytes());
        b.extend(24u64.to_le_bytes());
        // +X, −X, +Y, −Y, +Z, −Z.
        let centres: [[f64; 3]; 6] = [
            [1.0, 0.0, 0.0],
            [-1.0, -0.0, -0.0],
            [0.0, 1.0, 0.0],
            [-0.0, -1.0, -0.0],
            [0.0, 0.0, 1.0],
            [-0.0, -0.0, -1.0],
        ];
        for c in centres.iter().flatten() {
            b.extend(c.to_bits().to_le_bytes());
        }
        for k in 0..=6u32 {
            b.extend((4 * k).to_le_bytes());
        }
        let rings: [[u32; 4]; 6] = [
            [2, 4, 3, 5],
            [2, 5, 3, 4],
            [4, 0, 5, 1],
            [4, 1, 5, 0],
            [0, 2, 1, 3],
            [0, 3, 1, 2],
        ];
        for j in rings.iter().flatten() {
            b.extend(j.to_le_bytes());
        }
        // 2π/3 sr, acos(1/3) rad and π/2 rad (mpmath, as in
        // `SphereGrid::octahedron`).
        for _ in 0..6 {
            b.extend(2.094_395_102_393_195_7f64.to_bits().to_le_bytes());
        }
        for _ in 0..24 {
            b.extend(1.230_959_417_340_774_7f64.to_bits().to_le_bytes());
        }
        for _ in 0..24 {
            b.extend((PI / 2.0).to_bits().to_le_bytes());
        }
        b
    }

    #[test]
    fn the_octahedron_has_the_documented_layout() {
        let bytes = octahedron_bytes();
        assert_eq!(bytes.len(), 728);
        assert!(encode(&SphereGrid::octahedron()) == bytes);
        assert_eq!(
            read_grid(bytes.as_slice()).unwrap(),
            SphereGrid::octahedron()
        );
    }

    #[test]
    fn every_proper_prefix_is_truncated() {
        let bytes = octahedron_bytes();
        let parts = [
            (0, "magic number"),
            (8, "format version"),
            (12, "cell count"),
            (20, "neighbour entry count"),
            (28, "cell centres"),
            (172, "offsets"),
            (200, "neighbours"),
            (296, "areas"),
            (344, "edge lengths"),
            (536, "centre distances"),
        ];
        for len in 0..bytes.len() {
            let err = read_grid(&bytes[..len]).unwrap_err();
            assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof, "{len}: {err}");
            let part = parts
                .iter()
                .rev()
                .find(|(start, _)| *start <= len)
                .unwrap()
                .1;
            let message = err.to_string();
            assert!(
                message.starts_with(&format!("truncated .kgrid file: it ends inside the {part}")),
                "{len}: {message}"
            );
        }
        // Inside an array, the message says how far it got.
        let err = read_grid(&bytes[..300]).unwrap_err();
        assert_eq!(
            err.to_string(),
            "truncated .kgrid file: it ends inside the areas (after 0 of 6)"
        );
        let err = read_grid(&bytes[..28 + 24 * 5 + 7]).unwrap_err();
        assert_eq!(
            err.to_string(),
            "truncated .kgrid file: it ends inside the cell centres (after 15 of 18)"
        );
    }

    #[test]
    fn rejects_other_files_and_versions() {
        let mut bytes = octahedron_bytes();
        bytes[..8].copy_from_slice(b"KUNIFLD\0");
        let err = read_grid(bytes.as_slice()).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        assert_eq!(
            err.to_string(),
            "not a .kgrid file: it does not start with the magic number KUNIGRD\\x00"
        );
        for version in [0u32, 2, u32::MAX] {
            let mut bytes = octahedron_bytes();
            bytes[8..12].copy_from_slice(&version.to_le_bytes());
            let err = read_grid(bytes.as_slice()).unwrap_err();
            assert_eq!(err.kind(), io::ErrorKind::Unsupported);
            assert_eq!(
                err.to_string(),
                format!(
                    ".kgrid format version {version} is not supported: this build reads only \
                     version 1"
                )
            );
        }
    }

    #[test]
    fn rejects_counts_beyond_the_limits_before_reading_the_arrays() {
        let header = |n: u64, m: u64| {
            let mut b = octahedron_bytes()[..12].to_vec();
            b.extend(n.to_le_bytes());
            b.extend(m.to_le_bytes());
            b
        };
        let max = MAX_CELLS as u64;
        let err = read_grid(header(max + 1, 0).as_slice()).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        assert_eq!(
            err.to_string(),
            "the .kgrid file gives 715827885 cells, more than a grid can have (715827884)"
        );
        let err = read_grid(header(6, 1 << 32).as_slice()).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        assert_eq!(
            err.to_string(),
            "the .kgrid file gives 4294967296 neighbour entries, more than u32 offsets can \
             address (4294967295)"
        );
        // The largest counts with nothing behind them are truncated,
        // without reserving gigabytes first.
        let err = read_grid(header(max, u64::from(u32::MAX)).as_slice()).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof);
        assert_eq!(
            err.to_string(),
            "truncated .kgrid file: it ends inside the cell centres (after 0 of 2147483652)"
        );
    }

    #[test]
    fn rejects_arrays_that_are_not_a_grid() {
        // The first centre, at 28, made 2 instead of 1; the first area, at
        // 296, made 1 sr instead of 2π/3; the last offset, at 196, made 23
        // instead of 24; a neighbour of cell 0, at 200, made 0.
        let cases: [(usize, Vec<u8>, &str); 4] = [
            (
                28,
                2.0f64.to_bits().to_le_bytes().to_vec(),
                "cell 0 centre is not a unit vector",
            ),
            (
                296,
                1.0f64.to_bits().to_le_bytes().to_vec(),
                "cell areas are not positive and finite or do not sum to 4π",
            ),
            (196, 23u32.to_le_bytes().to_vec(), "CSR offsets are invalid"),
            (
                200,
                0u32.to_le_bytes().to_vec(),
                "cell 0 has invalid neighbour 0",
            ),
        ];
        for (at, patch, message) in cases {
            let mut bytes = octahedron_bytes();
            bytes[at..at + patch.len()].copy_from_slice(&patch);
            let err = read_grid(bytes.as_slice()).unwrap_err();
            assert_eq!(err.kind(), io::ErrorKind::InvalidData, "{message}");
            assert_eq!(
                err.to_string(),
                format!("the .kgrid file does not hold a valid grid: {message}")
            );
        }
    }
}

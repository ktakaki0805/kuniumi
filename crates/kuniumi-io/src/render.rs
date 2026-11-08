//! Maps of per-cell fields: equirectangular and orthographic (globe) images.
//!
//! A map is a row-major buffer of pixels, top row first. Each pixel takes
//! the cell nearest to the point of the sphere at its centre
//! ([`SphereGrid::nearest`]), and the colour of that cell's value
//! ([`Colormap::colours`]): the cells are drawn as their Voronoi polygons,
//! with no interpolation between them.
//!
//! - Equirectangular ([`equirect`], [`equirect_cells`]): `W × H` pixels with
//!   `H = ⌊W/2⌋`; pixel `(x, y)` is centred on longitude
//!   `λ = −π + (x + ½)·2π/W` and latitude `φ = π/2 − (y + ½)·π/H`, so the
//!   left edge is the antimeridian, longitude 0 is in the middle and north
//!   is up.
//! - Orthographic ([`orthographic`], [`orthographic_cells`]): the globe seen
//!   from infinitely far above the view centre `c` (latitude `φ₀`, longitude
//!   `λ₀`), in an `S × S` square. Pixel `(x, y)` has the view coordinates
//!   `u = 2(x + ½)/S − 1` (east) and `v = 1 − 2(y + ½)/S` (north); where
//!   `u² + v² > 1` it shows [`BACKGROUND`], elsewhere the point
//!   `p = w·c + u·e + v·n`, `w = √(1 − u² − v²)`, with `e` and `n` the east
//!   and north directions at `c` ([`coords::east`], [`coords::north`]). For
//!   an odd `S` the centre pixel has `u = v = 0` exactly (`2(x + ½) = S`),
//!   so it shows the cell nearest `c`.
//! - At a pole, `λ₀` sets the orientation of the globe: as computed,
//!   `from_lat_lon(±π/2, λ₀)` lies about 6e-17 off the axis towards `λ₀`,
//!   so `e` and `n` are those of a view centred on the meridian `λ₀` just
//!   short of the pole. The meridian `λ₀` runs from the centre downwards in
//!   a view of the north pole and upwards in a view of the south pole.
//!
//! The search for each pixel's cell starts from the previous pixel's cell,
//! and at the start of a row from the first cell of the row above; the
//! result does not depend on where it starts (M1-07), so the maps are the
//! same on every platform: positions from [`coords`] (built on
//! `kuniumi_math::det`), `+ − × ÷` and `sqrt`.

use core::f64::consts::{FRAC_PI_2, PI, TAU};

use kuniumi_sphere::SphereGrid;
use kuniumi_sphere::coords;

use crate::colormap::{Colormap, Range, Rgb8};

/// The colour outside the disc of a globe: a neutral dark grey that no
/// colour map uses.
pub const BACKGROUND: Rgb8 = [48, 48, 48];

/// The cell under each pixel of an equirectangular map `width` pixels wide
/// and `width / 2` high (module documentation).
///
/// # Panics
///
/// If `width < 2`.
pub fn equirect_cells(grid: &SphereGrid, width: u32) -> Vec<u32> {
    assert!(
        width >= 2,
        "an equirectangular map needs a width of at least 2"
    );
    let height = width / 2;
    let (dlon_rad, dlat_rad) = (TAU / f64::from(width), PI / f64::from(height));
    let mut cells = Vec::with_capacity(width as usize * height as usize);
    let mut row_start = 0;
    for y in 0..height {
        let lat_rad = FRAC_PI_2 - (f64::from(y) + 0.5) * dlat_rad;
        let mut hint = row_start;
        for x in 0..width {
            let lon_rad = -PI + (f64::from(x) + 0.5) * dlon_rad;
            hint = grid.nearest(coords::from_lat_lon(lat_rad, lon_rad), hint);
            if x == 0 {
                row_start = hint;
            }
            cells.push(hint);
        }
    }
    cells
}

/// The cell under each pixel of a globe `size × size` pixels centred on
/// latitude `centre_lat_rad` and longitude `centre_lon_rad`, `None` outside
/// the disc (module documentation).
///
/// # Panics
///
/// If `size` is 0 or the view centre is not finite.
pub fn orthographic_cells(
    grid: &SphereGrid,
    size: u32,
    centre_lat_rad: f64,
    centre_lon_rad: f64,
) -> Vec<Option<u32>> {
    assert!(size >= 1, "a globe needs a size of at least 1");
    assert!(
        centre_lat_rad.is_finite() && centre_lon_rad.is_finite(),
        "a globe needs a finite view centre, got ({centre_lat_rad}, {centre_lon_rad})"
    );
    let centre = coords::from_lat_lon(centre_lat_rad, centre_lon_rad);
    let (east, north) = (coords::east(centre), coords::north(centre));
    let mut cells = Vec::with_capacity(size as usize * size as usize);
    let mut row_start = None;
    let mut hint = 0;
    for y in 0..size {
        if let Some(cell) = row_start {
            hint = cell;
        }
        row_start = None;
        for x in 0..size {
            let (u, v) = view_coordinates(x, y, size);
            let r2 = u * u + v * v;
            if r2 > 1.0 {
                cells.push(None);
                continue;
            }
            let w = (1.0 - r2).sqrt();
            let p = centre * w + east * u + north * v;
            hint = grid.nearest(p, hint);
            row_start.get_or_insert(hint);
            cells.push(Some(hint));
        }
    }
    cells
}

/// The view coordinates `(u, v)` of the centre of pixel `(x, y)` of a globe
/// `size` pixels square: `u = 2(x + ½)/S − 1` (east) and
/// `v = 1 − 2(y + ½)/S` (north). Each is rounded once after the exact
/// `2(x + ½)` and division, so for an odd `S` the centre pixel gets exactly
/// `(0, 0)`.
fn view_coordinates(x: u32, y: u32, size: u32) -> (f64, f64) {
    let s = f64::from(size);
    (
        2.0 * (f64::from(x) + 0.5) / s - 1.0,
        1.0 - 2.0 * (f64::from(y) + 0.5) / s,
    )
}

/// An equirectangular map of `values`, one per cell, `width` pixels wide and
/// `width / 2` high, coloured by `colormap` over `range`.
///
/// # Panics
///
/// If `values` does not hold one value per cell, if `width < 2`, or as
/// [`Range::resolve`].
pub fn equirect(
    grid: &SphereGrid,
    values: &[f64],
    colormap: Colormap,
    range: Range,
    width: u32,
) -> Vec<Rgb8> {
    let colours = cell_colours(grid, values, colormap, range);
    equirect_cells(grid, width)
        .into_iter()
        .map(|cell| colours[cell as usize])
        .collect()
}

/// A globe view of `values`, one per cell, `size × size` pixels centred on
/// latitude `centre_lat_rad` and longitude `centre_lon_rad`, coloured by
/// `colormap` over `range`, with [`BACKGROUND`] outside the disc.
///
/// # Panics
///
/// If `values` does not hold one value per cell, if `size` is 0, if the
/// view centre is not finite, or as [`Range::resolve`].
pub fn orthographic(
    grid: &SphereGrid,
    values: &[f64],
    colormap: Colormap,
    range: Range,
    size: u32,
    centre_lat_rad: f64,
    centre_lon_rad: f64,
) -> Vec<Rgb8> {
    let colours = cell_colours(grid, values, colormap, range);
    orthographic_cells(grid, size, centre_lat_rad, centre_lon_rad)
        .into_iter()
        .map(|cell| cell.map_or(BACKGROUND, |c| colours[c as usize]))
        .collect()
}

/// The colour of each cell, after checking that `values` has one value per
/// cell.
fn cell_colours(grid: &SphereGrid, values: &[f64], colormap: Colormap, range: Range) -> Vec<Rgb8> {
    assert_eq!(
        values.len(),
        grid.len(),
        "a map needs one value per cell: {} values for {} cells",
        values.len(),
        grid.len()
    );
    colormap.colours(values, range)
}

// The tests use std's trigonometric functions as independent references.
#[cfg(test)]
#[allow(clippy::disallowed_methods)]
mod tests;

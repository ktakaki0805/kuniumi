//! Drainage metrics (design: "検証方法"; ticket M3-08): Hack's law and the
//! exceedance distribution of the drainage areas, on the rivers of the
//! final surface. See the module documentation of [`crate::metrics`],
//! "Drainage", for the definitions and their sources.

use core::f64::consts::PI;

use kuniumi_math::det;
use kuniumi_sphere::SphereGrid;
use kuniumi_world::World;

use crate::stages::crust::SURFACE_HEIGHT;
use crate::stages::erosion::{DRAINAGE_AREA, receiver_distances_m, routed_drainage};

/// The smallest basin of the fit of Hack's law in [`DrainageMetrics::new`],
/// in mean cell areas `4πR²/n`: 100 (ticket M3-08), so that every basin is
/// at least about ten cell spacings across (78,000 km² at 650,000 cells on
/// Earth).
pub const HACK_MIN_CELLS: f64 = 100.0;

/// The largest surface `4πR²` (m²) of [`DrainageMetrics::radius_in_limits`]:
/// 10³⁰⁰ m², so that every drainage area, its logarithm and every value of
/// the exceedance ladder above it stay finite (Earth: 5.1 × 10¹⁴ m²).
const MAX_SURFACE_M2: f64 = 1e300;

/// The smallest cell area `a_i R²` (m²) of
/// [`DrainageMetrics::radius_in_limits`]: 10⁻³⁰⁰ m², so that every drainage
/// area and every value of the exceedance ladder below it are normal
/// doubles (Earth at 650,000 cells: about 7.8 × 10⁸ m²).
const MIN_CELL_AREA_M2: f64 = 1e-300;

/// The highest height (m) of [`DrainageMetrics::heights_in_limits`]:
/// 10³⁰⁰ m, so that filling the depressions, which raises a cell by at most
/// one ULP per cell above its spill height, never passes `f64::MAX`
/// (`drainage::fill_depressions`).
const MAX_HEIGHT_M: f64 = 1e300;

/// The drainage areas that [`exceedance`] takes (m²): from 10⁻³⁰⁰ to 10³⁰⁷,
/// so that its ladder, from below the smallest to above the largest, stays
/// within the normal doubles.
const EXCEEDANCE_AREA_RANGE_M2: (f64, f64) = (1e-300, 1e307);

/// The exponents of ten between which [`exceedance`] searches its first
/// value: `10^−301` lies below every area it takes, `10^308` above.
const LADDER_POWERS: (i32, i32) = (-301, 308);

/// The main-stream length `L_i` (m) of every cell of a drainage network:
/// the length of the longest flow path from the cell up to a source (a
/// cell without donors), the sum of the distances to the receivers along
/// it (module documentation of [`crate::metrics`], "Drainage").
///
/// The network is given as in `stages::erosion::Network`: `receivers`, the
/// receiver of each cell (the cell itself for a base level); `stack`,
/// every cell once, each after its receiver (`Drainage::stack`); and
/// `receiver_distance_m`, the distance (m) from each cell to its receiver
/// (`stages::erosion::receiver_distances_m`; not read at the base levels).
/// Every `L` starts at +0 and the stack is visited backwards, each cell
/// after all its donors: `L_r = max(L_r, L_i + ℓ_i)` for each cell `i` that
/// is not a base level, `r` its receiver. A source has `L = 0`, and so has
/// a base level without donors. Finite unless a sum of distances exceeds
/// `f64::MAX`.
///
/// # Panics
///
/// If the slices differ in length; if the stack does not hold every cell
/// once, or a cell comes before its receiver in it; if the distance of a
/// cell that is not a base level is not finite and positive.
pub fn main_stream_lengths_m(
    receivers: &[u32],
    stack: &[u32],
    receiver_distance_m: &[f64],
) -> Vec<f64> {
    let n = receivers.len();
    assert_eq!(
        stack.len(),
        n,
        "main_stream_lengths_m needs a stack of every cell once: {} entries for {n} cells",
        stack.len()
    );
    assert_eq!(
        receiver_distance_m.len(),
        n,
        "main_stream_lengths_m needs one distance per cell: {} values for {n} cells",
        receiver_distance_m.len()
    );
    let mut position = vec![usize::MAX; n];
    for (k, &i) in stack.iter().enumerate() {
        let i = i as usize;
        assert!(
            i < n && position[i] == usize::MAX,
            "main_stream_lengths_m needs a stack of every cell once: entry {k} is {i}"
        );
        position[i] = k;
    }
    for (i, &r) in receivers.iter().enumerate() {
        let r = r as usize;
        if r == i {
            continue;
        }
        assert!(
            r < n && position[r] < position[i],
            "main_stream_lengths_m needs each cell after its receiver in the stack: \
             cell {i} drains to {r}"
        );
        let l = receiver_distance_m[i];
        // Written so that NaN fails.
        assert!(
            l > 0.0 && l.is_finite(),
            "main_stream_lengths_m needs a finite, positive distance to the receiver: \
             cell {i} has {l}"
        );
    }
    let mut length_m = vec![0.0; n];
    for &i in stack.iter().rev() {
        let i = i as usize;
        let r = receivers[i] as usize;
        if r != i {
            let l = length_m[i] + receiver_distance_m[i];
            if l > length_m[r] {
                length_m[r] = l;
            }
        }
    }
    length_m
}

/// The least-squares fit of Hack's law `L = c A^h` (module documentation of
/// [`crate::metrics`], "Drainage").
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HackFit {
    /// Hack's exponent `h` (dimensionless); Earth about 0.6
    /// ([`super::earth::HACK_EXPONENT`]).
    pub exponent: f64,
    /// Hack's coefficient `c`, in m^(1−2h) (L in metres for A in square
    /// metres).
    pub coefficient: f64,
}

impl HackFit {
    /// The fit of `ln L` on `ln A` over the basins given, one pair per
    /// basin: their drainage areas `drainage_area_m2` (m²) and main-stream
    /// lengths `main_stream_length_m` (m).
    ///
    /// With `x = ln A` and `y = ln L` (`det::ln`): the means `x̄ = Σx / N`
    /// and `ȳ = Σy / N`, then `S_xx = Σ (x − x̄)²` and
    /// `S_xy = Σ (x − x̄)(y − ȳ)`, every sum in the order of the pairs from
    /// +0; `h = S_xy / S_xx` and `c = exp(ȳ − h x̄)` (`det::exp`). `None`
    /// with fewer than two pairs, when every `x` is the same (all areas
    /// equal: the rounded mean `x̄` need not equal them, so `S_xx` alone
    /// could be a rounding error rather than 0), when `S_xx = 0`, or when
    /// `h` or `c` is not finite: `c` overflows only for exponents far from
    /// any river network (below about −30 for areas near 10¹⁰ m²).
    ///
    /// # Panics
    ///
    /// If the slices differ in length, or an area or a length is not finite
    /// and positive.
    pub fn new(drainage_area_m2: &[f64], main_stream_length_m: &[f64]) -> Option<Self> {
        assert_eq!(
            main_stream_length_m.len(),
            drainage_area_m2.len(),
            "HackFit::new needs one length per area: {} lengths for {} areas",
            main_stream_length_m.len(),
            drainage_area_m2.len()
        );
        for (what, values) in [("area", drainage_area_m2), ("length", main_stream_length_m)] {
            // Written so that NaN fails.
            if let Some(k) = values.iter().position(|&v| !(v > 0.0 && v.is_finite())) {
                panic!(
                    "HackFit::new needs finite, positive values: the {what} of basin {k} is {}",
                    values[k]
                );
            }
        }
        if drainage_area_m2.len() < 2 {
            return None;
        }
        let x: Vec<f64> = drainage_area_m2.iter().map(|&a| det::ln(a)).collect();
        // The logarithms of positive finite areas are finite and never −0,
        // so equal bits are equal values.
        if x.iter().all(|xi| xi.to_bits() == x[0].to_bits()) {
            return None;
        }
        let y: Vec<f64> = main_stream_length_m.iter().map(|&l| det::ln(l)).collect();
        // Exact: fewer than 2⁵³ basins.
        let count = x.len() as f64;
        let x_mean = x.iter().fold(0.0, |sum, &v| sum + v) / count;
        let y_mean = y.iter().fold(0.0, |sum, &v| sum + v) / count;
        let (mut s_xx, mut s_xy) = (0.0, 0.0);
        for (&xi, &yi) in x.iter().zip(&y) {
            let dx = xi - x_mean;
            s_xx += dx * dx;
            s_xy += dx * (yi - y_mean);
        }
        if s_xx == 0.0 {
            return None;
        }
        let exponent = s_xy / s_xx;
        let coefficient = det::exp(y_mean - exponent * x_mean);
        (exponent.is_finite() && coefficient.is_finite()).then_some(Self {
            exponent,
            coefficient,
        })
    }
}

/// The share of the cells whose drainage area is at least a given area: one
/// value of the exceedance distribution ([`exceedance`]).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Exceedance {
    /// The mantissa `m` of the area, 1, 2 or 5.
    pub mantissa: u32,
    /// The power of ten `k` of the area.
    pub power_of_ten: i32,
    /// The area `a = m × 10^k` (m²), as [`exceedance`] computes it.
    pub drainage_area_m2: f64,
    /// The share of the weight of the cells whose drainage area is at
    /// least `a`, in [0, 1].
    pub share: f64,
}

/// The exceedance distribution of the drainage areas `drainage_area_m2`
/// (m²) of the cells given, weighted by `weight` (any unit, for example
/// their areas in sr): for the areas `a = m × 10^k`, `m` ∈ {1, 2, 5}, from
/// the largest not above the smallest drainage area to the smallest above
/// the largest, the share `Σ_{A_i ≥ a} w_i / Σ w_i`, both sums in cell order
/// from +0. The first share is 1 and the last 0; without cells the list is
/// empty.
///
/// `a` is `m · 10^k` for `k ≥ 0` and `m / 10^−k` for `k < 0`, with `10^j`
/// the product of `j` tens: exact for `0 ≤ j ≤ 22`, so for `|k| ≤ 22`
/// (areas from 10⁻²² to 5 × 10²² m², those of every planet) `a` is the
/// double nearest to `m × 10^k`; beyond, it is the result of the rounded
/// steps, the same on every platform.
///
/// # Panics
///
/// If the slices differ in length; if a weight is not finite and positive,
/// or the weights sum to infinity; if a drainage area is not finite and
/// between 10⁻³⁰⁰ and 10³⁰⁷ m².
pub fn exceedance(weight: &[f64], drainage_area_m2: &[f64]) -> Vec<Exceedance> {
    assert_eq!(
        drainage_area_m2.len(),
        weight.len(),
        "exceedance needs one drainage area per weight: {} areas for {} weights",
        drainage_area_m2.len(),
        weight.len()
    );
    // Written so that NaN fails.
    if let Some(i) = weight.iter().position(|&w| !(w > 0.0 && w.is_finite())) {
        panic!(
            "exceedance needs finite, positive weights: cell {i} has {}",
            weight[i]
        );
    }
    let total = weight.iter().fold(0.0, |sum, &w| sum + w);
    assert!(
        total.is_finite(),
        "exceedance needs weights whose sum is finite, got {total}"
    );
    let (lowest, highest) = EXCEEDANCE_AREA_RANGE_M2;
    if let Some(i) = drainage_area_m2
        .iter()
        .position(|&a| !(a >= lowest && a <= highest))
    {
        panic!(
            "exceedance needs drainage areas from {lowest:e} to {highest:e} m²: cell {i} has {:e}",
            drainage_area_m2[i]
        );
    }
    let Some((smallest, largest)) = drainage_area_m2.iter().fold(None, |range, &a| match range {
        None => Some((a, a)),
        Some((lo, hi)) => Some((if a < lo { a } else { lo }, if a > hi { a } else { hi })),
    }) else {
        return Vec::new();
    };
    let mut rung = Rung::at_or_below(smallest);
    let mut ladder = Vec::new();
    loop {
        let a = rung.value();
        let above =
            weight.iter().zip(drainage_area_m2).fold(
                0.0,
                |sum, (&w, &area)| if area >= a { sum + w } else { sum },
            );
        ladder.push(Exceedance {
            mantissa: rung.mantissa,
            power_of_ten: rung.power,
            drainage_area_m2: a,
            share: above / total,
        });
        if a > largest {
            return ladder;
        }
        rung = rung.next();
    }
}

/// A value `m × 10^k` of the ladder of [`exceedance`], `m` ∈ {1, 2, 5}.
#[derive(Clone, Copy, Debug)]
struct Rung {
    /// `m`.
    mantissa: u32,
    /// `k`.
    power: i32,
}

impl Rung {
    /// The largest value of the ladder not above `area_m2`, an area within
    /// [`EXCEEDANCE_AREA_RANGE_M2`]: the power `k` with
    /// `10^k ≤ area < 10^(k+1)` by bisection over [`LADDER_POWERS`], then
    /// the largest mantissa whose value is not above the area.
    fn at_or_below(area_m2: f64) -> Self {
        let (mut low, mut high) = LADDER_POWERS;
        // 10^low ≤ area < 10^high; the values grow strictly with k.
        while high - low > 1 {
            let mid = low + (high - low) / 2;
            if power_of_ten(mid) <= area_m2 {
                low = mid;
            } else {
                high = mid;
            }
        }
        let mut rung = Self {
            mantissa: 1,
            power: low,
        };
        for mantissa in [2, 5] {
            let next = Self { mantissa, ..rung };
            if next.value() <= area_m2 {
                rung = next;
            }
        }
        rung
    }

    /// The next value of the ladder: 1 → 2 → 5 → 10.
    fn next(self) -> Self {
        match self.mantissa {
            1 => Self {
                mantissa: 2,
                ..self
            },
            2 => Self {
                mantissa: 5,
                ..self
            },
            _ => Self {
                mantissa: 1,
                power: self.power + 1,
            },
        }
    }

    /// `m · 10^k` for `k ≥ 0`, `m / 10^−k` for `k < 0` ([`exceedance`]).
    fn value(self) -> f64 {
        let m = f64::from(self.mantissa);
        if self.power >= 0 {
            m * power_of_ten(self.power)
        } else {
            m / power_of_ten(-self.power)
        }
    }
}

/// `10^j` for `j ≥ 0` as the product of `j` tens, exact for `j ≤ 22`; for
/// `j < 0`, `1 / 10^−j`. Finite for `j ≤ 308`, normal for `j ≥ −307`.
fn power_of_ten(j: i32) -> f64 {
    let mut p = 1.0;
    for _ in 0..j.unsigned_abs() {
        p *= 10.0;
    }
    if j >= 0 { p } else { 1.0 / p }
}

/// The drainage metrics of a surface (module documentation of
/// [`crate::metrics`], "Drainage"): what the design's validation table asks
/// of the river networks.
#[derive(Clone, Debug, PartialEq)]
pub struct DrainageMetrics {
    /// The smallest basin of the fit, in mean cell areas `4πR²/n`
    /// ([`HACK_MIN_CELLS`] for [`new`](Self::new)).
    pub hack_min_cells: f64,
    /// The smallest basin of the fit (m²): `hack_min_cells · 4πR²/n`.
    pub hack_min_area_m2: f64,
    /// The number of basins of the fit: the land cells whose drainage area
    /// is at least `hack_min_area_m2` and whose main-stream length is
    /// positive.
    pub hack_basins: usize,
    /// [`HackFit::new`] over the basins; `None` with fewer than two.
    pub hack: Option<HackFit>,
    /// The largest drainage area of a land cell (m²); `None` without land.
    pub largest_drainage_area_m2: Option<f64>,
    /// [`exceedance`] of the land cells' drainage areas, weighted by their
    /// areas; empty without land.
    pub exceedance: Vec<Exceedance>,
}

impl DrainageMetrics {
    /// The metrics of the surface `height_m` (m, one height per cell of
    /// `grid`, on the datum of `sea_level_m`) of a planet of radius
    /// `radius_m` (m), with the basins of the fit at least
    /// [`HACK_MIN_CELLS`] mean cell areas:
    /// [`with_min_cells`](Self::with_min_cells).
    ///
    /// # Panics
    ///
    /// As [`with_min_cells`](Self::with_min_cells).
    pub fn new(grid: &SphereGrid, height_m: &[f64], sea_level_m: f64, radius_m: f64) -> Self {
        Self::with_min_cells(grid, height_m, sea_level_m, radius_m, HACK_MIN_CELLS)
    }

    /// The metrics of the surface `height_m` with the basins of the fit at
    /// least `min_cells` mean cell areas.
    ///
    /// The drainage graph is that of the surface routed across its closed
    /// depressions (`stages::erosion::routed_drainage`, the graph of the
    /// Erosion stage's rivers); the drainage area of a cell is the sum of
    /// `a_i R²` (m²) over it and every cell upstream, as in the field
    /// `drainage_area`; the main-stream lengths come from
    /// [`main_stream_lengths_m`] with the distances of
    /// `stages::erosion::receiver_distances_m`. Land is `h_i > z`, as in the
    /// drainage graph (the elevation metrics count a cell at the sea level
    /// as land; the graph makes it an outlet). The basins are the land cells
    /// with `A_i ≥ min_cells · 4πR²/n` and `L_i > 0` (a basin of one cell
    /// has no stream); the exceedance is that of the land cells, weighted by
    /// their areas in sr.
    ///
    /// # Panics
    ///
    /// If `height_m` does not hold one height per cell, the heights are
    /// outside [`heights_in_limits`](Self::heights_in_limits), the sea level
    /// is not finite, the radius is outside
    /// [`radius_in_limits`](Self::radius_in_limits), `min_cells` is not
    /// finite and positive, or the threshold `min_cells · 4πR²/n` overflows
    /// (above `f64::MAX · n / 4πR²` cells, about 1.8 × 10⁸ · n on the
    /// largest surface, 10³⁰⁰ m²; a basin holds at most `n` cells, so any
    /// `min_cells` above `n` gives no basin).
    pub fn with_min_cells(
        grid: &SphereGrid,
        height_m: &[f64],
        sea_level_m: f64,
        radius_m: f64,
        min_cells: f64,
    ) -> Self {
        // Written so that NaN fails.
        assert!(
            min_cells > 0.0 && min_cells.is_finite(),
            "DrainageMetrics needs a finite, positive basin size, got {min_cells} cells"
        );
        let (area_m2, length_m) = areas_and_lengths(grid, height_m, sea_level_m, radius_m);
        // Exact: fewer than 2⁵³ cells.
        let mean_cell_m2 = 4.0 * PI * radius_m * radius_m / grid.len() as f64;
        let hack_min_area_m2 = min_cells * mean_cell_m2;
        assert!(
            hack_min_area_m2.is_finite(),
            "DrainageMetrics needs a basin size whose area is finite, got {min_cells:e} cells \
             of {mean_cell_m2:e} m²"
        );
        let land: Vec<usize> = (0..grid.len())
            .filter(|&i| height_m[i] > sea_level_m)
            .collect();
        let (basin_area_m2, basin_length_m): (Vec<f64>, Vec<f64>) = land
            .iter()
            .filter(|&&i| area_m2[i] >= hack_min_area_m2 && length_m[i] > 0.0)
            .map(|&i| (area_m2[i], length_m[i]))
            .unzip();
        let largest_drainage_area_m2 = land.iter().map(|&i| area_m2[i]).fold(None, |m, a| {
            Some(match m {
                Some(m) if m >= a => m,
                _ => a,
            })
        });
        let land_sr: Vec<f64> = land.iter().map(|&i| grid.area_sr(i as u32)).collect();
        let land_area_m2: Vec<f64> = land.iter().map(|&i| area_m2[i]).collect();
        Self {
            hack_min_cells: min_cells,
            hack_min_area_m2,
            hack_basins: basin_area_m2.len(),
            hack: HackFit::new(&basin_area_m2, &basin_length_m),
            largest_drainage_area_m2,
            exceedance: exceedance(&land_sr, &land_area_m2),
        }
    }

    /// The metrics of `world` ([`new`](Self::new)): its field
    /// [`SURFACE_HEIGHT`] relative to `world.sea_level_m`, on its grid and
    /// the radius of its planet; `None` if the world has no such field or
    /// no field [`DRAINAGE_AREA`] (the Erosion stage has not run; the values
    /// of `drainage_area` are not read).
    ///
    /// # Panics
    ///
    /// As [`new`](Self::new).
    pub fn of_world(world: &World) -> Option<Self> {
        let height = world.field(SURFACE_HEIGHT)?;
        world.field(DRAINAGE_AREA)?;
        Some(Self::new(
            &world.grid,
            &height.values,
            world.sea_level_m,
            world.params.radius_m,
        ))
    }

    /// Whether the heights `height_m` (m) are within the limits of the
    /// metrics: finite and at most 10³⁰⁰ m, so that filling the depressions
    /// (`drainage::fill_depressions`, which raises a cell by at most one
    /// ULP per cell above its spill height) stays finite. Real surfaces are
    /// far inside (Earth: below 10⁴ m).
    pub fn heights_in_limits(height_m: &[f64]) -> bool {
        height_m.iter().all(|&h| height_in_limits(h))
    }

    /// Whether the radius `radius_m` (m) of a planet with the cells of
    /// `grid` is within the limits of the metrics: finite and positive,
    /// with a surface `4πR²` of at most 10³⁰⁰ m², every cell area `a_i R²`
    /// at least 10⁻³⁰⁰ m², and every distance between neighbouring centres,
    /// `d_ij R`, positive. Real planets are far inside (Earth: 5.1 × 10¹⁴ m²
    /// in all, 7.8 × 10⁸ m² per cell at 650,000 cells). Within them every
    /// drainage area, main-stream length, logarithm and value of the
    /// exceedance ladder is finite and normal, so no metric is NaN or
    /// infinite.
    pub fn radius_in_limits(grid: &SphereGrid, radius_m: f64) -> bool {
        // Written so that NaN fails.
        if !(radius_m > 0.0 && radius_m.is_finite()) {
            return false;
        }
        let surface_m2 = 4.0 * PI * radius_m * radius_m;
        surface_m2 <= MAX_SURFACE_M2
            && (0..grid.len() as u32).all(|i| {
                grid.area_sr(i) * radius_m * radius_m >= MIN_CELL_AREA_M2
                    && grid
                        .centre_distances_rad(i)
                        .iter()
                        .all(|&d| d * radius_m > 0.0)
            })
    }
}

/// Whether the height `height_m` (m) is finite and at most 10³⁰⁰ m
/// ([`DrainageMetrics::heights_in_limits`]).
fn height_in_limits(height_m: f64) -> bool {
    height_m.is_finite() && height_m <= MAX_HEIGHT_M
}

/// The drainage area (m²) and the main-stream length (m) of every cell of
/// the surface `height_m`, on the graph of [`DrainageMetrics::with_min_cells`].
///
/// # Panics
///
/// As [`DrainageMetrics::with_min_cells`], but for `min_cells`.
fn areas_and_lengths(
    grid: &SphereGrid,
    height_m: &[f64],
    sea_level_m: f64,
    radius_m: f64,
) -> (Vec<f64>, Vec<f64>) {
    assert!(
        DrainageMetrics::radius_in_limits(grid, radius_m),
        "DrainageMetrics needs a radius within its limits, got {radius_m}"
    );
    if let Some(i) = height_m.iter().position(|&h| !height_in_limits(h)) {
        panic!(
            "DrainageMetrics needs finite heights of at most 1e300 m: cell {i} has {:e}",
            height_m[i]
        );
    }
    let (drainage, _) = routed_drainage(grid, height_m, sea_level_m, radius_m);
    // The expression of `stages::erosion::Rivers`, so that the areas are
    // those of the field `drainage_area` bit for bit.
    let cell_m2: Vec<f64> = (0..grid.len() as u32)
        .map(|i| grid.area_sr(i) * radius_m * radius_m)
        .collect();
    let distance_m = receiver_distances_m(grid, &drainage, radius_m);
    (
        drainage.accumulate(&cell_m2),
        main_stream_lengths_m(drainage.receivers(), drainage.stack(), &distance_m),
    )
}

#[cfg(test)]
mod tests;

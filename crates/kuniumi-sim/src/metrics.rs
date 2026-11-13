//! Validation metrics: numbers that compare a world with Earth (design:
//! "検証方法"; ticket M2-06).
//!
//! The design asks that worlds be judged by statistics that can be measured
//! on Earth, not by eye. This module holds pure functions over per-cell
//! values; they read no files and draw no random numbers. So far: the land
//! fraction and the hypsometry of the solid surface ([`ElevationMetrics`]).
//! Drainage and climate metrics follow (tickets M3-05 and M5-06).
//!
//! # Elevation
//!
//! Every metric is weighted by the cells' areas `a_i` and uses the heights
//! `h_i` (m, on the datum of the reference radius) relative to the sea level
//! `z` (m): the elevation `x_i = h_i − z`, correctly rounded.
//!
//! - A cell is **ocean** where it lies below the sea level, `h_i < z`
//!   (water stands `z − h_i > 0` deep on it), and **land** everywhere else,
//!   a cell exactly at the sea level included; a planet without water, whose
//!   sea level is its lowest height (`stages::crust::sea_level_m`), is all
//!   land. `x_i < 0` exactly for the ocean cells: the difference of two
//!   finite doubles is zero only when they are equal.
//! - Sums run in cell order, sorts by [`f64::total_cmp`] with ties broken by
//!   the cell index (docs/conventions.md, determinism rules 5 and 6). Only
//!   `+ − × ÷`, comparisons, `floor` and conversions of integer-valued
//!   `f64` are used, so the results are the same on every platform.
//!
//! [`earth`] holds the values of Earth that the design compares with.
//!
//! # Limits
//!
//! Every metric takes the areas `a_i`, the heights `h_i` and the sea level
//! `z` and panics (naming the first offending cell where there is one)
//! unless (ticket M2-15):
//!
//! 1. there is one area per height, and at least one cell;
//! 2. every area is finite and positive, and their total, summed in cell
//!    order, lies below 2¹⁰²³ ≈ 8.99e307, half the largest double. Merely
//!    finite would not do: the percentiles sum the areas in the order of
//!    the elevations, which rounds differently (the areas `MAX/2, MAX/2,
//!    2⁹⁶⁹, 2⁹⁶⁹` sum to `f64::MAX` in that order, to infinity with the
//!    small ones first); with the margin, sums of the areas in any order
//!    stay finite (for fewer than 2⁵⁰ cells);
//! 3. every height and the sea level are finite, and so is every elevation
//!    `x_i = h_i − z` (heights and sea level below 2¹⁰²³ in magnitude are
//!    enough).
//!
//! [`Hypsometry::new`] also needs every `x_i / w` to be finite (true for
//! any bin width `w ≥ 1 m`) and at most [`MAX_HYPSOMETRY_BINS`] bins, 10⁶:
//! a relief of about 100,000 km in the world's bins of 100 m (ticket
//! M2-18), which [`Hypsometry::bins_needed`] counts beforehand;
//! [`percentiles_m`] fractions in (0, 1]. The grid of a world has areas in
//! steradians summing to 4π, so for [`ElevationMetrics::of_world`] only the
//! heights and the sea level matter.
//!
//! Within these limits no metric returns NaN or infinity. The land fraction
//! and the shares are ratios of finite sums of areas; the percentiles and
//! the extremes are elevations. The means `Σ a_i v_i / Σ a_i` of elevations
//! and depths `v_i ≥ 0` keep this form where it is finite; where it
//! overflows (two cells of area 1 at 1e308 m, or a sea level 1.7e308 m
//! above a grid of 4π sr) they are `Σ (a_i / W) v_i`, with `W = Σ a_i`,
//! capped at the largest `v_i`, since a weighted mean lies within its
//! values. Products `a_i v_i` in the subnormal range (below about 2e-308)
//! lose precision.

use kuniumi_world::World;

use crate::stages::crust::SURFACE_HEIGHT;

/// The width (m) of the bins of the elevation histogram ([`Hypsometry`]):
/// 100 m (ticket M2-06).
pub const HYPSOMETRY_BIN_M: f64 = 100.0;

/// The fractions of the surface at which [`ElevationMetrics`] reports the
/// percentiles of the elevation ([`percentiles_m`]): 1, 5, 10, 25, 50, 75,
/// 90, 95 and 99 %.
pub const PERCENTILES: [f64; 9] = [0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99];

/// The largest share of the surface, relative to the smaller of the two
/// modes, that the valley between them may hold for the histogram to have
/// two separate peaks ([`Hypsometry::has_separate_peaks`]): ½.
pub const VALLEY_RATIO: f64 = 0.5;

/// The most bins a histogram may hold ([`Hypsometry::new`]): 10⁶ (ticket
/// M2-18).
///
/// A bound on memory and on the length of `metrics.txt` (in
/// `kuniumi_io::world_dir`), which lists every bin, empty ones included: at
/// the bound the shares take 8 MB and `metrics.txt` about 29 MB, where 2³¹
/// bins, the bound before, took 17 GB. In the world's bins of 100 m
/// ([`HYPSOMETRY_BIN_M`]) it is a relief of 100,000 km, almost four times the
/// diameter of the largest planet that `PlanetParams` allows (a radius of
/// 13,000 km, 260,000 bins across), which no solid surface comes near; the
/// relief of the Crust stage needs about 150 bins.
///
/// Elevations spanning `D` metres need between `⌊D / w⌋ + 1` and
/// `⌈D / w⌉ + 1` bins of width `w`, depending on where the edges fall: in
/// 100 m bins a relief below 99,999.9 km always fits, one of 100,000 km or
/// more never does. [`Hypsometry::bins_needed`] gives the exact number.
pub const MAX_HYPSOMETRY_BINS: usize = 1_000_000;

/// The bound on the total area of the cells, 2¹⁰²³ (module documentation,
/// "Limits").
const MAX_TOTAL_AREA: f64 = 8.988_465_674_311_58e307;

/// Approximate values of Earth, which the metrics are compared with.
///
/// From the table of the design document ("検証方法"), where each is given
/// as "about" (約), and for the two modes of the hypsometry from ticket
/// M2-06; heights are elevations above the sea level (m).
pub mod earth {
    /// The share of the surface that is land: about 0.29 (design: 約29%).
    pub const LAND_FRACTION: f64 = 0.29;

    /// The mean elevation of the land (m): about 800 m (design: 約0.8 km).
    pub const MEAN_LAND_ELEVATION_M: f64 = 800.0;

    /// The mean depth of the ocean (m): about 3,700 m (design: 約3.7 km).
    pub const MEAN_OCEAN_DEPTH_M: f64 = 3_700.0;

    /// The highest elevation (m): about +8,800 m (design: 最高峰 約8.8 km).
    pub const HIGHEST_M: f64 = 8_800.0;

    /// The lowest elevation (m): about −11,000 m (design: 最深部 約11 km).
    pub const LOWEST_M: f64 = -11_000.0;

    /// The mode of the land's elevation (m): about +100 m (ticket M2-06).
    pub const LAND_MODE_M: f64 = 100.0;

    /// The mode of the ocean floor's elevation (m): about −4,500 m (ticket
    /// M2-06).
    pub const OCEAN_MODE_M: f64 = -4_500.0;
}

/// The share of the surface that is land (module documentation):
/// `Σ_land a_i / Σ a_i`, both sums in cell order.
///
/// `area` holds the cells' areas (any unit, for example sr), `height_m`
/// their heights (m) and `sea_level_m` the sea level (m) on the same datum.
///
/// # Panics
///
/// Beyond the limits of the module documentation: if the slices are empty
/// or differ in length, an area is not finite and positive, the areas sum
/// to 2¹⁰²³ or more, a height or the sea level is not finite, or an
/// elevation `h_i − z` overflows.
pub fn land_fraction(area: &[f64], height_m: &[f64], sea_level_m: f64) -> f64 {
    check(area, height_m, sea_level_m);
    let (mut land, mut total) = (0.0, 0.0);
    for (&a, &h) in area.iter().zip(height_m) {
        total += a;
        if !is_ocean(h, sea_level_m) {
            land += a;
        }
    }
    land / total
}

/// The mean elevation of the land (m above the sea level),
/// `Σ_land a_i x_i / Σ_land a_i`; `None` without land.
///
/// # Panics
///
/// As [`land_fraction`].
pub fn mean_land_elevation_m(area: &[f64], height_m: &[f64], sea_level_m: f64) -> Option<f64> {
    check(area, height_m, sea_level_m);
    weighted_mean(area, height_m, |h| {
        (!is_ocean(h, sea_level_m)).then_some(h - sea_level_m)
    })
}

/// The mean depth of the ocean (m below the sea level, positive),
/// `Σ_ocean a_i (z − h_i) / Σ_ocean a_i`; `None` without ocean.
///
/// # Panics
///
/// As [`land_fraction`].
pub fn mean_ocean_depth_m(area: &[f64], height_m: &[f64], sea_level_m: f64) -> Option<f64> {
    check(area, height_m, sea_level_m);
    weighted_mean(area, height_m, |h| {
        is_ocean(h, sea_level_m).then_some(sea_level_m - h)
    })
}

/// Whether a cell of height `height_m` lies below the sea level
/// `sea_level_m` (module documentation).
fn is_ocean(height_m: f64, sea_level_m: f64) -> bool {
    height_m < sea_level_m
}

/// `Σ a_i v_i / Σ a_i` over the cells for which `value` of the height gives
/// some `v_i`, sums in cell order; `None` if it gives none.
///
/// The values are finite and not negative, and `W = Σ a_i` is finite (the
/// limits of the module documentation). Should the quotient overflow, the
/// mean is `Σ (a_i / W) v_i` in cell order instead, whose weights are at
/// most 1, capped at the largest `v_i`.
fn weighted_mean(
    area: &[f64],
    height_m: &[f64],
    value: impl Fn(f64) -> Option<f64>,
) -> Option<f64> {
    let (mut sum, mut weight, mut any) = (0.0, 0.0, false);
    for (&a, &h) in area.iter().zip(height_m) {
        if let Some(v) = value(h) {
            sum += a * v;
            weight += a;
            any = true;
        }
    }
    if !any {
        return None;
    }
    let mean = sum / weight;
    if mean.is_finite() {
        return Some(mean);
    }
    let (mut mean, mut largest) = (0.0, 0.0);
    for (&a, &h) in area.iter().zip(height_m) {
        if let Some(v) = value(h) {
            mean += a / weight * v;
            if v > largest {
                largest = v;
            }
        }
    }
    Some(if mean > largest { largest } else { mean })
}

/// The area-weighted percentiles of the elevation (m) at `fractions` of the
/// surface, one per fraction, in their order.
///
/// The cells are sorted by elevation ([`f64::total_cmp`], ties by cell
/// index) and their areas summed in that order; the percentile at `p` is
/// the elevation of the first cell at which the running sum reaches `p`
/// times the final sum. It is the lower quantile: with two cells of equal
/// area the median is the lower one.
///
/// # Panics
///
/// As [`land_fraction`], and if a fraction is not in (0, 1].
pub fn percentiles_m(
    area: &[f64],
    height_m: &[f64],
    sea_level_m: f64,
    fractions: &[f64],
) -> Vec<f64> {
    check(area, height_m, sea_level_m);
    if let Some(p) = fractions.iter().find(|&&p| !(p > 0.0 && p <= 1.0)) {
        panic!("percentiles need a fraction in (0, 1], got {p}");
    }
    let elevation: Vec<f64> = height_m.iter().map(|&h| h - sea_level_m).collect();
    let mut order: Vec<usize> = (0..elevation.len()).collect();
    order.sort_unstable_by(|&i, &j| elevation[i].total_cmp(&elevation[j]).then(i.cmp(&j)));
    let mut running = Vec::with_capacity(order.len());
    let mut sum = 0.0;
    for &i in &order {
        sum += area[i];
        running.push(sum);
    }
    fractions
        .iter()
        .map(|&p| {
            let target = p * sum;
            // The running sums never decrease; the last one is `sum`, which
            // reaches every target since `p ≤ 1`.
            let k = running.partition_point(|&s| s < target);
            elevation[order[k]]
        })
        .collect()
}

/// The area-weighted histogram of the elevation: the share of the surface
/// in bins of equal width with edges at whole multiples of the width from
/// the sea level.
///
/// Cell `i` falls into the bin `k_i = floor(x_i / w)`, with the quotient
/// correctly rounded (so an elevation within a rounding error of an edge
/// may fall into the bin above it); an ocean cell whose quotient rounds to
/// −0, which only a subnormal `x_i` does, falls into bin −1. So the bins
/// with `k ≥ 0` hold land and those with `k < 0` ocean, and bin `k` covers
/// the elevations `[k w, (k + 1) w)`. The histogram runs from the lowest
/// `k_i` to the highest, bins without cells included; each holds the share
/// `Σ_{k_i = k} a_i / Σ a_i` (sums in cell order), and the shares sum to 1
/// within rounding.
///
/// The indices `k` are kept as integer-valued `f64`, so every finite input
/// works without overflow: on a water world whose sea level lies 1e293 m
/// above the surface every cell falls into one bin (that far out, more than
/// 2⁵³ widths from the sea level, the edges are rounded and `k + 1 = k`).
/// The number of bins is the span of the elevations over the width, plus
/// one ([`bins_needed`](Self::bins_needed)): about 150 for the relief of the
/// Crust stage (−8 km to +5.6 km) in 100 m bins, and at most
/// [`MAX_HYPSOMETRY_BINS`].
#[derive(Clone, Debug, PartialEq)]
pub struct Hypsometry {
    /// The width `w` of the bins (m).
    bin_width_m: f64,
    /// The index `k` of the first bin, an integer.
    first_bin: f64,
    /// The share of the surface in each bin, lowest first.
    shares: Vec<f64>,
}

impl Hypsometry {
    /// The histogram of the elevations of the cells with the areas `area`
    /// and the heights `height_m` (m) relative to the sea level
    /// `sea_level_m` (m), in bins `bin_width_m` (m) wide (type
    /// documentation; the world's histogram uses [`HYPSOMETRY_BIN_M`]).
    ///
    /// # Panics
    ///
    /// As [`land_fraction`]; if the bin width is not finite and positive,
    /// if an elevation over the width overflows, or if more than
    /// [`MAX_HYPSOMETRY_BINS`] (10⁶) bins would be needed
    /// ([`bins_needed`](Self::bins_needed)).
    pub fn new(area: &[f64], height_m: &[f64], sea_level_m: f64, bin_width_m: f64) -> Self {
        check(area, height_m, sea_level_m);
        let (first, count) = bin_span(height_m, sea_level_m, bin_width_m);
        assert!(
            count <= MAX_HYPSOMETRY_BINS as f64,
            "the elevations need {count:e} bins of {bin_width_m} m, more than \
             {MAX_HYPSOMETRY_BINS:e}"
        );
        let mut areas = vec![0.0; count as usize];
        let mut total = 0.0;
        for (&a, &h) in area.iter().zip(height_m) {
            // `k − first` lies in [0, count − 1]: both are integers, and
            // rounding is monotonic.
            let k = bin_index(h - sea_level_m, bin_width_m);
            areas[(k - first) as usize] += a;
            total += a;
        }
        Self {
            bin_width_m,
            first_bin: first,
            shares: areas.into_iter().map(|a| a / total).collect(),
        }
    }

    /// The number of bins that [`new`](Self::new) needs for the elevations of
    /// the heights `height_m` (m) above the sea level `sea_level_m` (m), in
    /// bins `bin_width_m` (m) wide: `k_max − k_min + 1` for the bin indices
    /// `k_i` of the type documentation, by the same arithmetic as `new`. The
    /// areas do not matter.
    ///
    /// `new` panics where the number exceeds [`MAX_HYPSOMETRY_BINS`]; a
    /// caller that rejects such inputs instead checks it first
    /// (`kuniumi_io::world_dir::write` does). It is an integer-valued `f64`,
    /// like the indices, since elevations far apart need more bins than any
    /// integer type holds: up to about 3.6e306 in bins of 100 m, and
    /// infinitely many in bins of about 2 m or less, where `k_max − k_min`
    /// can overflow (never NaN). Beyond 2⁵³ it is rounded.
    ///
    /// # Panics
    ///
    /// As `new` for the heights (limits 1 and 3 of the module
    /// documentation): if `height_m` is empty, a height or the sea level is
    /// not finite, or an elevation `h − z` overflows; if the bin width is not
    /// finite and positive, or an elevation over the width overflows.
    pub fn bins_needed(height_m: &[f64], sea_level_m: f64, bin_width_m: f64) -> f64 {
        assert!(
            !height_m.is_empty(),
            "the histogram needs at least one height"
        );
        check_heights(height_m, sea_level_m);
        bin_span(height_m, sea_level_m, bin_width_m).1
    }

    /// The width of the bins (m).
    pub fn bin_width_m(&self) -> f64 {
        self.bin_width_m
    }

    /// The share of the surface in each bin, from the lowest bin to the
    /// highest; at least one bin.
    pub fn shares(&self) -> &[f64] {
        &self.shares
    }

    /// The lower edge (m above the sea level) of the bin `bin` (counted from
    /// the lowest, 0): `k w`.
    pub fn lower_edge_m(&self, bin: usize) -> f64 {
        self.index(bin) * self.bin_width_m
    }

    /// The upper edge (m above the sea level) of the bin `bin`: `(k + 1) w`.
    pub fn upper_edge_m(&self, bin: usize) -> f64 {
        (self.index(bin) + 1.0) * self.bin_width_m
    }

    /// The centre (m above the sea level) of the bin `bin`: `(k + ½) w`.
    pub fn centre_m(&self, bin: usize) -> f64 {
        (self.index(bin) + 0.5) * self.bin_width_m
    }

    /// Whether the bin `bin` holds land (`k ≥ 0`) rather than ocean.
    pub fn is_land(&self, bin: usize) -> bool {
        self.index(bin) >= 0.0
    }

    /// The index `k` of the bin `bin`.
    fn index(&self, bin: usize) -> f64 {
        self.first_bin + bin as f64
    }

    /// The land mode: the land bin with the largest share (the lowest of
    /// equal ones); `None` without land.
    pub fn land_mode(&self) -> Option<usize> {
        self.extreme(0..self.shares.len(), |bin| self.is_land(bin), |a, b| a > b)
    }

    /// The ocean mode: the ocean bin with the largest share (the lowest of
    /// equal ones); `None` without ocean.
    pub fn ocean_mode(&self) -> Option<usize> {
        self.extreme(0..self.shares.len(), |bin| !self.is_land(bin), |a, b| a > b)
    }

    /// The valley: the bin strictly between the ocean mode and the land
    /// mode with the smallest share (the lowest of equal ones); `None`
    /// without one of the modes or when they are adjacent. The ocean bins
    /// all lie below the land bins, so the valley lies above the ocean mode.
    pub fn valley(&self) -> Option<usize> {
        let (ocean, land) = (self.ocean_mode()?, self.land_mode()?);
        self.extreme(ocean + 1..land, |_| true, |a, b| a < b)
    }

    /// Whether the histogram has two separate peaks: both modes exist and
    /// the [`valley`](Self::valley) between them holds at most
    /// [`VALLEY_RATIO`] of the smaller of their shares. Adjacent modes are
    /// not separate.
    pub fn has_separate_peaks(&self) -> bool {
        match (self.ocean_mode(), self.valley(), self.land_mode()) {
            (Some(ocean), Some(valley), Some(land)) => {
                let (o, l) = (self.shares[ocean], self.shares[land]);
                let smaller = if o < l { o } else { l };
                self.shares[valley] <= VALLEY_RATIO * smaller
            }
            _ => false,
        }
    }

    /// Of the bins in `bins` for which `keep` holds, the first whose share
    /// is `better` than that of every bin before it; `None` if `keep` holds
    /// for none.
    fn extreme(
        &self,
        bins: core::ops::Range<usize>,
        keep: impl Fn(usize) -> bool,
        better: impl Fn(f64, f64) -> bool,
    ) -> Option<usize> {
        bins.filter(|&bin| keep(bin))
            .fold(None, |best, bin| match best {
                Some(b) if !better(self.shares[bin], self.shares[b]) => Some(b),
                _ => Some(bin),
            })
    }
}

/// The land fraction and the hypsometry of a surface (module
/// documentation): what the design's validation table asks of the relief.
#[derive(Clone, Debug, PartialEq)]
pub struct ElevationMetrics {
    /// The sea level (m, on the datum of the heights).
    pub sea_level_m: f64,
    /// [`land_fraction`].
    pub land_fraction: f64,
    /// [`mean_land_elevation_m`] (m).
    pub mean_land_elevation_m: Option<f64>,
    /// [`mean_ocean_depth_m`] (m).
    pub mean_ocean_depth_m: Option<f64>,
    /// The highest elevation (m above the sea level).
    pub highest_m: f64,
    /// The lowest elevation (m above the sea level, negative below it).
    pub lowest_m: f64,
    /// [`percentiles_m`] at [`PERCENTILES`] (m).
    pub percentiles_m: Vec<f64>,
    /// The histogram in bins of [`HYPSOMETRY_BIN_M`].
    pub hypsometry: Hypsometry,
}

impl ElevationMetrics {
    /// The metrics of the cells with the areas `area` and the heights
    /// `height_m` (m) relative to the sea level `sea_level_m` (m).
    ///
    /// # Panics
    ///
    /// As [`land_fraction`], and as [`Hypsometry::new`] with
    /// [`HYPSOMETRY_BIN_M`]: if the elevations need more than
    /// [`MAX_HYPSOMETRY_BINS`] (10⁶) bins of 100 m, a relief of about
    /// 100,000 km ([`Hypsometry::bins_needed`]).
    pub fn new(area: &[f64], height_m: &[f64], sea_level_m: f64) -> Self {
        check(area, height_m, sea_level_m);
        let (lowest_m, highest_m) = height_m
            .iter()
            .map(|&h| h - sea_level_m)
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), x| {
                (if x < lo { x } else { lo }, if x > hi { x } else { hi })
            });
        Self {
            sea_level_m,
            land_fraction: land_fraction(area, height_m, sea_level_m),
            mean_land_elevation_m: mean_land_elevation_m(area, height_m, sea_level_m),
            mean_ocean_depth_m: mean_ocean_depth_m(area, height_m, sea_level_m),
            highest_m,
            lowest_m,
            percentiles_m: percentiles_m(area, height_m, sea_level_m, &PERCENTILES),
            hypsometry: Hypsometry::new(area, height_m, sea_level_m, HYPSOMETRY_BIN_M),
        }
    }

    /// The metrics of `world`: its field [`SURFACE_HEIGHT`] (m) on the areas
    /// of its grid (sr), relative to `world.sea_level_m`; `None` if the
    /// world has no such field (the Crust stage has not run).
    ///
    /// # Panics
    ///
    /// As [`new`](Self::new): if the field does not hold one value per
    /// cell, a height or the sea level is not finite, an elevation above the
    /// sea level overflows, or the elevations need more than
    /// [`MAX_HYPSOMETRY_BINS`] bins. (The areas of a grid sum to 4π, far
    /// below the limit.)
    pub fn of_world(world: &World) -> Option<Self> {
        let height = world.field(SURFACE_HEIGHT)?;
        let area: Vec<f64> = (0..world.grid.len() as u32)
            .map(|cell| world.grid.area_sr(cell))
            .collect();
        Some(Self::new(&area, &height.values, world.sea_level_m))
    }
}

/// Checks the inputs of the metrics (module documentation, "Limits").
///
/// # Panics
///
/// If `area` is empty or not as long as `height_m`, an area is not finite
/// and positive, the areas sum (in cell order) to 2¹⁰²³ or more, or as
/// [`check_heights`].
fn check(area: &[f64], height_m: &[f64], sea_level_m: f64) {
    assert!(
        !area.is_empty() && area.len() == height_m.len(),
        "one area per height, at least one: {} areas, {} heights",
        area.len(),
        height_m.len()
    );
    if let Some(cell) = area.iter().position(|a| !(a.is_finite() && *a > 0.0)) {
        panic!(
            "the area of cell {cell} is not finite and positive: {}",
            area[cell]
        );
    }
    let total = area.iter().fold(0.0, |sum, &a| sum + a);
    assert!(
        total < MAX_TOTAL_AREA,
        "the areas of the cells sum to {total:e}, not below 2^1023"
    );
    check_heights(height_m, sea_level_m);
}

/// Checks the heights and the sea level (module documentation, "Limits",
/// limit 3).
///
/// # Panics
///
/// If a height or `sea_level_m` is not finite, or an elevation
/// `h − sea_level_m` overflows.
fn check_heights(height_m: &[f64], sea_level_m: f64) {
    if let Some(cell) = height_m.iter().position(|h| !h.is_finite()) {
        panic!(
            "the height of cell {cell} is not finite: {}",
            height_m[cell]
        );
    }
    assert!(
        sea_level_m.is_finite(),
        "the sea level is not finite: {sea_level_m}"
    );
    if let Some(cell) = height_m
        .iter()
        .position(|&h| !(h - sea_level_m).is_finite())
    {
        panic!(
            "the elevation of cell {cell} above the sea level overflows: height {:e} m, \
             sea level {sea_level_m:e} m",
            height_m[cell]
        );
    }
}

/// The index `k` of the bin of [`Hypsometry`] that holds the elevation
/// `elevation_m` (m) in bins `bin_width_m` (m) wide: `floor(x / w)`, or −1
/// for an ocean cell whose quotient rounds to −0 (type documentation).
///
/// # Panics
///
/// If the quotient overflows.
fn bin_index(elevation_m: f64, bin_width_m: f64) -> f64 {
    let k = (elevation_m / bin_width_m).floor();
    assert!(
        k.is_finite(),
        "the elevation {elevation_m} m over the bin width {bin_width_m} m overflows"
    );
    if elevation_m < 0.0 && k >= 0.0 {
        -1.0
    } else {
        k
    }
}

/// The index of the lowest bin of the elevations of `height_m` above
/// `sea_level_m` in bins `bin_width_m` wide, and the number of bins up to
/// the highest, `k_max − k_min + 1` ([`bin_index`], in cell order). The
/// heights are checked already and not empty.
///
/// # Panics
///
/// If the bin width is not finite and positive, or as [`bin_index`].
fn bin_span(height_m: &[f64], sea_level_m: f64, bin_width_m: f64) -> (f64, f64) {
    assert!(
        bin_width_m.is_finite() && bin_width_m > 0.0,
        "the histogram needs a bin width that is finite and positive, got {bin_width_m}"
    );
    let (first, last) = height_m
        .iter()
        .map(|&h| bin_index(h - sea_level_m, bin_width_m))
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), k| {
            (if k < lo { k } else { lo }, if k > hi { k } else { hi })
        });
    (first, last - first + 1.0)
}

#[cfg(test)]
mod tests;

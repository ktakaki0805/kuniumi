//! Colour maps for per-cell fields: values to 8-bit sRGB colours.
//!
//! Four maps ([`Colormap`]):
//!
//! - [`Colormap::Viridis`], perceptually uniform and sequential: viridis by
//!   Nathaniel J. Smith, Stefan van der Walt and Eric Firing (2015), released
//!   under CC0. `VIRIDIS` holds every 5th of the 256 entries of
//!   `_viridis_data` in `colormaps.py` of the BIDS colormap project
//!   (<https://github.com/BIDS/colormap>, commit `bc549477`; the table
//!   matplotlib ships), so 52 control points at `t = k/51`. Interpolated
//!   linearly between them, the map deviates from the full table by at most
//!   0.62 of an 8-bit unit (`tools/colormap_tables.py`).
//! - [`Colormap::Diverging`]: ColorBrewer RdBu with 11 classes, reversed so
//!   that low values are blue and high values red, white at the centre.
//! - [`Colormap::Elevation`]: blues below sea level and greens and browns
//!   above it, with a hard break at sea level; the colours are kuniumi's own.
//! - [`Colormap::Categorical`]: ColorBrewer Paired with 12 classes, for class
//!   ids.
//!
//! The ColorBrewer colours are by Cynthia Brewer, Mark Harrower and The
//! Pennsylvania State University, under the Apache License 2.0 (taken from
//! `colorbrewer_schemes.js`, <https://github.com/axismaps/colorbrewer>,
//! commit `7d135fc4`): "This product includes color specifications and
//! designs developed by Cynthia Brewer (<http://colorbrewer.org/>)." See
//! `THIRD_PARTY.md`.
//!
//! # Arithmetic
//!
//! Interpolation is linear in sRGB, in `f64`; each channel, in 8-bit units,
//! is clamped to [0, 255] and rounded half away from zero
//! ([`f64::round`]). Besides `+ − × ÷` and comparisons, the only operations
//! are that rounding and conversions from `f64` to integers of values that
//! are already integers or are truncated on purpose (the segment index),
//! all exact, so the colours are the same on every platform.

/// An 8-bit sRGB colour: red, green, blue.
pub type Rgb8 = [u8; 3];

/// The colour of a value that has no colour: any non-finite value, and a
/// categorical value that is not a class id. Magenta, so that it stands out.
pub const NO_DATA: Rgb8 = [255, 0, 255];

/// Viridis control points, sRGB channels in [0, 1] at `t = k/51`:
/// entries 0, 5, …, 255 of the BIDS table (module documentation).
const VIRIDIS: [[f64; 3]; 52] = [
    [0.267004, 0.004874, 0.329415],
    [0.273809, 0.031497, 0.358853],
    [0.278791, 0.062145, 0.386592],
    [0.281924, 0.089666, 0.412415],
    [0.283197, 0.115680, 0.436115],
    [0.282623, 0.140926, 0.457517],
    [0.280255, 0.165693, 0.476498],
    [0.276194, 0.190074, 0.493001],
    [0.270595, 0.214069, 0.507052],
    [0.263663, 0.237631, 0.518762],
    [0.255645, 0.260703, 0.528312],
    [0.246811, 0.283237, 0.535941],
    [0.237441, 0.305202, 0.541921],
    [0.227802, 0.326594, 0.546532],
    [0.218130, 0.347432, 0.550038],
    [0.208623, 0.367752, 0.552675],
    [0.199430, 0.387607, 0.554642],
    [0.190631, 0.407061, 0.556089],
    [0.182256, 0.426184, 0.557120],
    [0.174274, 0.445044, 0.557792],
    [0.166617, 0.463708, 0.558119],
    [0.159194, 0.482237, 0.558073],
    [0.151918, 0.500685, 0.557587],
    [0.144759, 0.519093, 0.556572],
    [0.137770, 0.537492, 0.554906],
    [0.131172, 0.555899, 0.552459],
    [0.125394, 0.574318, 0.549086],
    [0.121148, 0.592739, 0.544641],
    [0.119423, 0.611141, 0.538982],
    [0.121380, 0.629492, 0.531973],
    [0.128087, 0.647749, 0.523491],
    [0.140210, 0.665859, 0.513427],
    [0.157851, 0.683765, 0.501686],
    [0.180653, 0.701402, 0.488189],
    [0.208030, 0.718701, 0.472873],
    [0.239374, 0.735588, 0.455688],
    [0.274149, 0.751988, 0.436601],
    [0.311925, 0.767822, 0.415586],
    [0.352360, 0.783011, 0.392636],
    [0.395174, 0.797475, 0.367757],
    [0.440137, 0.811138, 0.340967],
    [0.487026, 0.823929, 0.312321],
    [0.535621, 0.835785, 0.281908],
    [0.585678, 0.846661, 0.249897],
    [0.636902, 0.856542, 0.216620],
    [0.688944, 0.865448, 0.182725],
    [0.741388, 0.873449, 0.149561],
    [0.793760, 0.880678, 0.120005],
    [0.845561, 0.887322, 0.099702],
    [0.896320, 0.893616, 0.096335],
    [0.945636, 0.899815, 0.112838],
    [0.993248, 0.906157, 0.143936],
];

/// ColorBrewer RdBu, 11 classes, reversed: control points at `t = k/10`.
const DIVERGING: [Rgb8; 11] = [
    [5, 48, 97],
    [33, 102, 172],
    [67, 147, 195],
    [146, 197, 222],
    [209, 229, 240],
    [247, 247, 247],
    [253, 219, 199],
    [244, 165, 130],
    [214, 96, 77],
    [178, 24, 43],
    [103, 0, 31],
];

/// Below sea level, from the deepest value to just below sea level: deep
/// navy to light blue (kuniumi's own choice).
const ELEVATION_OCEAN: [Rgb8; 4] = [[8, 24, 68], [22, 70, 140], [58, 128, 196], [150, 200, 235]];

/// At and above sea level, from sea level to the highest value: green,
/// yellow-green, tan, brown and a pale grey-brown for the highest ground
/// (kuniumi's own choice).
const ELEVATION_LAND: [Rgb8; 5] = [
    [46, 120, 56],
    [138, 172, 82],
    [212, 190, 122],
    [150, 102, 62],
    [222, 212, 200],
];

/// The categorical palette: ColorBrewer Paired, 12 classes. Class `k` has
/// the colour `CATEGORICAL[k % 12]`.
pub const CATEGORICAL: [Rgb8; 12] = [
    [166, 206, 227],
    [31, 120, 180],
    [178, 223, 138],
    [51, 160, 44],
    [251, 154, 153],
    [227, 26, 28],
    [253, 191, 111],
    [255, 127, 0],
    [202, 178, 214],
    [106, 61, 154],
    [255, 255, 153],
    [177, 89, 40],
];

/// The largest class id a categorical value may hold: below 2⁵³, so that
/// every id is an exact `f64`.
const MAX_CLASS: f64 = 9_007_199_254_740_992.0;

/// How the values of a field become colours.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Colormap {
    /// Viridis, perceptually uniform and sequential: dark purple at the
    /// minimum, yellow at the maximum.
    Viridis,
    /// Blue below `centre`, white at it, red above it (ColorBrewer RdBu),
    /// with the same scale on both sides: the half-width of the map is
    /// `m = max(centre − min, max − centre)`, so the value `centre ± m` takes
    /// the darkest colour of its side.
    Diverging {
        /// The value drawn white, in the unit of the field.
        centre: f64,
    },
    /// Elevation, in metres, with a hard break at `sea_level`: blues from
    /// the minimum to just below sea level, greens and browns from sea level
    /// to the maximum. Each side has its own scale.
    Elevation {
        /// The height of the sea surface, on the same reference as the
        /// values, in metres.
        sea_level: f64,
    },
    /// Class ids: a value that is a non-negative integer `k` below 2⁵³ takes
    /// the colour `CATEGORICAL[k % 12]`; any other value is [`NO_DATA`]. The
    /// range is ignored.
    Categorical,
}

/// The range of values a colour map spans.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Range {
    /// From the smallest to the largest finite value of the field.
    Auto,
    /// From `min` to `max`, both finite, `min ≤ max`. Values outside are
    /// clamped by each colour map: [`Colormap::Viridis`] and
    /// [`Colormap::Elevation`] draw them with the colour of the nearer
    /// bound, [`Colormap::Diverging`] clamps at `centre ± m` (so beyond the
    /// nearer bound a value still darkens), and [`Colormap::Categorical`]
    /// ignores the range.
    Fixed {
        /// The lower bound, in the unit of the field.
        min: f64,
        /// The upper bound, in the unit of the field.
        max: f64,
    },
}

impl Range {
    /// The bounds `(min, max)` for the field `values`: for [`Range::Auto`]
    /// the smallest and largest finite value (`(0, 0)` if there is none),
    /// for [`Range::Fixed`] the given bounds.
    ///
    /// # Panics
    ///
    /// If a fixed range has a non-finite bound or `min > max`.
    pub fn resolve(self, values: &[f64]) -> (f64, f64) {
        match self {
            Range::Auto => {
                // Plain comparisons in index order: the first of equal
                // values (+0 or −0) is kept.
                let mut bounds: Option<(f64, f64)> = None;
                for &v in values.iter().filter(|v| v.is_finite()) {
                    bounds = Some(match bounds {
                        None => (v, v),
                        Some((min, max)) => {
                            (if v < min { v } else { min }, if v > max { v } else { max })
                        }
                    });
                }
                bounds.unwrap_or((0.0, 0.0))
            }
            Range::Fixed { min, max } => {
                assert!(
                    min.is_finite() && max.is_finite() && min <= max,
                    "a fixed range needs finite bounds with min <= max, got [{min}, {max}]"
                );
                (min, max)
            }
        }
    }
}

impl Colormap {
    /// The colour of `value` on a map spanning `min` to `max` (the result of
    /// [`Range::resolve`]): [`NO_DATA`] for a non-finite value, else as
    /// described for each map. For [`Colormap::Viridis`], `value` is placed
    /// at `t = (value − min)/(max − min)`, clamped to [0, 1], and at
    /// `t = 0.5` when `max = min`.
    ///
    /// # Panics
    ///
    /// If the `centre` of [`Colormap::Diverging`] or the `sea_level` of
    /// [`Colormap::Elevation`] is not finite.
    pub fn colour(self, value: f64, min: f64, max: f64) -> Rgb8 {
        if let Colormap::Diverging { centre: level } | Colormap::Elevation { sea_level: level } =
            self
        {
            assert!(level.is_finite(), "{self:?} needs a finite level");
        }
        if !value.is_finite() {
            return NO_DATA;
        }
        match self {
            Colormap::Viridis => viridis(fraction(value, min, max)),
            Colormap::Diverging { centre } => {
                let half_width = f64::max(centre - min, max - centre);
                let t = if half_width > 0.0 {
                    0.5 + (value - centre) / (2.0 * half_width)
                } else {
                    0.5
                };
                diverging(t)
            }
            Colormap::Elevation { sea_level } => {
                if value < sea_level {
                    let t = if sea_level > min {
                        (value - min) / (sea_level - min)
                    } else {
                        1.0
                    };
                    interpolate(&ELEVATION_OCEAN, 1.0, t)
                } else {
                    let t = if max > sea_level {
                        (value - sea_level) / (max - sea_level)
                    } else {
                        0.0
                    };
                    interpolate(&ELEVATION_LAND, 1.0, t)
                }
            }
            Colormap::Categorical => {
                if (0.0..MAX_CLASS).contains(&value) && value.fract() == 0.0 {
                    // An exact conversion: a non-negative integer below 2⁵³.
                    categorical(value as u64)
                } else {
                    NO_DATA
                }
            }
        }
    }

    /// The colour of each of `values`, on the range `range` resolved over
    /// `values` ([`Range::resolve`]).
    ///
    /// # Panics
    ///
    /// As [`Range::resolve`] and [`colour`](Self::colour).
    pub fn colours(self, values: &[f64], range: Range) -> Vec<Rgb8> {
        let (min, max) = range.resolve(values);
        values.iter().map(|&v| self.colour(v, min, max)).collect()
    }
}

/// Viridis at `t`, clamped to [0, 1]: dark purple `#440154` at 0, yellow
/// `#FDE725` at 1.
pub fn viridis(t: f64) -> Rgb8 {
    interpolate(&VIRIDIS, 255.0, t)
}

/// The diverging map at `t`, clamped to [0, 1]: dark blue at 0,
/// `[247, 247, 247]` at 0.5, dark red at 1.
pub fn diverging(t: f64) -> Rgb8 {
    interpolate(&DIVERGING, 1.0, t)
}

/// The colour of class `class`: `CATEGORICAL[class % 12]`.
pub fn categorical(class: u64) -> Rgb8 {
    CATEGORICAL[(class % CATEGORICAL.len() as u64) as usize]
}

/// `(value − min)/(max − min)`, or 0.5 when `max ≤ min` (the colour maps
/// clamp it to [0, 1]).
fn fraction(value: f64, min: f64, max: f64) -> f64 {
    if max > min {
        (value - min) / (max - min)
    } else {
        0.5
    }
}

/// The colour at `t` (clamped to [0, 1]) on the polyline through `points`,
/// which are equally spaced from `t = 0` to `t = 1`; a channel of `points`
/// times `scale` is in 8-bit units (`scale` is 255 for channels in [0, 1]
/// and 1 for 8-bit channels).
///
/// A NaN `t` gives [`NO_DATA`]. Only overflow produces one: `∞/∞` when both
/// the distance of a value from a bound and the width of the range exceed
/// the largest `f64` (about 1.8e308).
fn interpolate<T: Copy + Into<f64>>(points: &[[T; 3]], scale: f64, t: f64) -> Rgb8 {
    if t.is_nan() {
        return NO_DATA;
    }
    let segments = points.len() - 1;
    let s = t.clamp(0.0, 1.0) * segments as f64;
    // A truncating cast of a value in [0, segments]; the last point is the
    // end of the last segment.
    let k = (s as usize).min(segments - 1);
    let f = s - k as f64;
    let (a, b) = (points[k], points[k + 1]);
    core::array::from_fn(|c| {
        let (a, b) = (a[c].into() * scale, b[c].into() * scale);
        to_u8(a + f * (b - a))
    })
}

/// An 8-bit channel value: clamped to [0, 255] and rounded half away from
/// zero.
fn to_u8(channel: f64) -> u8 {
    // In [0, 255] after rounding, so the cast is exact.
    channel.clamp(0.0, 255.0).round() as u8
}

#[cfg(test)]
#[allow(clippy::disallowed_methods)]
mod tests;

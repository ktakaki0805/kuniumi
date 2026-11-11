//! Charts drawn into pixel buffers, for images that humans and Claude Code
//! look at: so far bar charts of histograms ([`histogram`]).
//!
//! A chart is a row-major buffer of [`Rgb8`] pixels, top row first, ready for
//! [`crate::png::encode_rgb8`]. It holds no text: the numbers it shows, and
//! its axes, go into a text report next to it (ticket M2-06).
//!
//! # Layout
//!
//! An image `W × H` pixels holds the plot, `P = W − 2M` columns by
//! `Q = H − 2M` rows, inside a margin of [`MARGIN_PX`] (`M`) on every side,
//! all on a white [`BACKGROUND`]. Plot column `c` (image column `M + c`)
//! covers the values `[x_min + c Δ, x_min + (c + 1) Δ)` of the x axis, with
//! `Δ = (x_max − x_min)/P`; plot row `r` (image row `M + r`) counts from the
//! top. In drawing order:
//!
//! 1. [`GRID`] lines across the plot at every whole multiple of `x_step`
//!    from `x_min` to `x_max` (vertical) and of `y_step` above 0 up to
//!    `y_max` (horizontal). The value `x` lies in the column
//!    `⌊(x − x_min)/(x_max − x_min) · P⌋`, `x_max` in the last one; the value
//!    `y` on the row `Q − round(y / y_max · Q)`, `y_max` on the first.
//! 2. [`AXIS`] lines just outside the plot, which nothing else touches:
//!    below it (image row `M + Q`, columns `M − 1` to `M + P − 1`) and left
//!    of it (image column `M − 1`, rows `M` to `M + Q`).
//! 3. The bars ([`Bar`]): each column shows the tallest bar that overlaps it
//!    (the first of equal ones), from the bottom of the plot up,
//!    `round(value / y_max · Q)` rows high, at least one row for a positive
//!    value and at most `Q`, in the bar's colour.
//! 4. The markers ([`Marker`]), vertical lines over the whole plot, in
//!    their colour, over the bars.
//!
//! Only `+ − × ÷`, comparisons, `floor`, `ceil`, `round` and conversions of
//! integer-valued `f64` are used, so the pixels are the same on every
//! platform. Any finite bounds work: subnormal ones, and ones whose
//! difference `x_max − x_min` would overflow (the differences are then
//! taken of their halves).

use crate::Rgb8;

/// The margin around the plot (pixels): 20.
pub const MARGIN_PX: u32 = 20;

/// The colour of the background: white.
pub const BACKGROUND: Rgb8 = [255, 255, 255];

/// The colour of the grid lines: a light grey.
pub const GRID: Rgb8 = [221, 221, 221];

/// The colour of the axis lines: black.
pub const AXIS: Rgb8 = [0, 0, 0];

/// A bar: the interval `[lower, upper)` of the x axis and its height
/// `value` on the y axis, in the units of the [`Axes`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Bar {
    /// The lower end of the interval.
    pub lower: f64,
    /// The upper end of the interval (excluded).
    pub upper: f64,
    /// The height; a bar of height 0 or less is not drawn.
    pub value: f64,
    /// The colour.
    pub colour: Rgb8,
}

/// The ranges of the axes and the spacing of their grid lines.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Axes {
    /// The value at the left edge of the plot.
    pub x_min: f64,
    /// The value at the right edge of the plot, above `x_min`.
    pub x_max: f64,
    /// The spacing of the vertical grid lines, positive.
    pub x_step: f64,
    /// The value at the top of the plot, positive; the bottom is 0.
    pub y_max: f64,
    /// The spacing of the horizontal grid lines, positive.
    pub y_step: f64,
}

/// A vertical line at `x`, drawn over the bars; not drawn outside the range
/// of the x axis.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Marker {
    /// The value on the x axis.
    pub x: f64,
    /// The colour.
    pub colour: Rgb8,
}

/// A bar chart of `bars` on `axes`, with `markers`, `width_px × height_px`
/// pixels (module documentation).
///
/// # Panics
///
/// If the image is not larger than the margins (at least 41 × 41 pixels),
/// if a bound or step of `axes` is not finite, `x_min ≥ x_max`, `y_max`
/// or a step is not positive, if there would be more grid lines than
/// columns or rows, or if a bar or marker holds a value that is not
/// finite.
pub fn histogram(
    bars: &[Bar],
    axes: &Axes,
    markers: &[Marker],
    width_px: u32,
    height_px: u32,
) -> Vec<Rgb8> {
    assert!(
        width_px > 2 * MARGIN_PX && height_px > 2 * MARGIN_PX,
        "a chart of {width_px} x {height_px} pixels is not larger than its margins"
    );
    let Axes {
        x_min,
        x_max,
        x_step,
        y_max,
        y_step,
    } = *axes;
    assert!(
        x_min.is_finite() && x_max.is_finite() && x_min < x_max,
        "a chart needs finite x_min < x_max, got {x_min} and {x_max}"
    );
    assert!(
        y_max.is_finite() && y_max > 0.0,
        "a chart needs a finite, positive y_max, got {y_max}"
    );
    assert!(
        x_step.is_finite() && x_step > 0.0 && y_step.is_finite() && y_step > 0.0,
        "a chart needs finite, positive steps, got {x_step} and {y_step}"
    );
    let (width, margin) = (width_px as usize, MARGIN_PX as usize);
    let (columns, rows) = (width_px - 2 * MARGIN_PX, height_px - 2 * MARGIN_PX);
    let first_line = (x_min / x_step).ceil();
    let x_lines = (x_max / x_step).floor() - first_line + 1.0;
    let y_lines = (y_max / y_step).floor();
    assert!(
        x_lines <= f64::from(columns) && y_lines <= f64::from(rows),
        "{x_lines} vertical and {y_lines} horizontal lines: more grid lines than \
         the {columns} columns or {rows} rows"
    );
    for bar in bars {
        assert!(
            bar.lower.is_finite() && bar.upper.is_finite() && bar.value.is_finite(),
            "a bar holds a value that is not finite: {bar:?}"
        );
    }
    for marker in markers {
        assert!(
            marker.x.is_finite(),
            "a marker holds a value that is not finite: {marker:?}"
        );
    }

    // The position of `x` in columns from the left edge of the plot. Only
    // when `x_max − x_min` overflows are the differences taken of halves,
    // which keeps them finite; those bounds are far above the subnormal
    // range, where halving is exact and commutes with rounding.
    let range = x_max - x_min;
    let at = |x: f64| {
        if range.is_finite() {
            (x - x_min) / range * f64::from(columns)
        } else {
            (0.5 * x - 0.5 * x_min) / (0.5 * x_max - 0.5 * x_min) * f64::from(columns)
        }
    };
    // The column of `x`, which lies in the range: `x_max` in the last one.
    let column = |x: f64| {
        let c = at(x).floor();
        if c < f64::from(columns) {
            c as usize
        } else {
            columns as usize - 1
        }
    };
    // The plot row of the height `y` in [0, y_max]: the top row of a bar of
    // that height.
    let row = |y: f64| {
        let r = (y / y_max * f64::from(rows)).round();
        if r < f64::from(rows) {
            rows as usize - r as usize
        } else {
            0
        }
    };
    let mut pixels = vec![BACKGROUND; width * height_px as usize];
    let mut set = |c: usize, r: usize, colour: Rgb8| {
        pixels[(margin + r) * width + margin + c] = colour;
    };

    // 1. Grid lines.
    let mut k = first_line;
    while k * x_step <= x_max {
        let c = column(k * x_step);
        for r in 0..rows as usize {
            set(c, r, GRID);
        }
        k += 1.0;
    }
    let mut k = 1.0;
    while k <= y_lines {
        let r = row(k * y_step);
        for c in 0..columns as usize {
            set(c, r, GRID);
        }
        k += 1.0;
    }

    // 3. Bars: the tallest in each column, the first of equal ones.
    let mut tallest: Vec<Option<(f64, Rgb8)>> = vec![None; columns as usize];
    for bar in bars.iter().filter(|b| b.value > 0.0 && b.lower < b.upper) {
        let (lo, hi) = (at(bar.lower).floor(), at(bar.upper).ceil());
        let first = if lo > 0.0 { lo } else { 0.0 };
        let end = if hi < f64::from(columns) {
            hi
        } else {
            f64::from(columns)
        };
        let mut c = first;
        while c < end {
            let slot = &mut tallest[c as usize];
            if slot.is_none_or(|(value, _)| bar.value > value) {
                *slot = Some((bar.value, bar.colour));
            }
            c += 1.0;
        }
    }
    for (c, slot) in tallest.into_iter().enumerate() {
        if let Some((value, colour)) = slot {
            let r = (value / y_max * f64::from(rows)).round();
            let height = if r < 1.0 {
                1
            } else if r < f64::from(rows) {
                r as usize
            } else {
                rows as usize
            };
            for r in rows as usize - height..rows as usize {
                set(c, r, colour);
            }
        }
    }

    // 4. Markers.
    for marker in markers.iter().filter(|m| m.x >= x_min && m.x <= x_max) {
        let c = column(marker.x);
        for r in 0..rows as usize {
            set(c, r, marker.colour);
        }
    }

    // 2. Axis lines, outside the plot: below it and left of it.
    let bottom = (margin + rows as usize) * width;
    for x in margin - 1..margin + columns as usize {
        pixels[bottom + x] = AXIS;
    }
    for y in margin..=margin + rows as usize {
        pixels[y * width + margin - 1] = AXIS;
    }
    pixels
}

#[cfg(test)]
mod tests;

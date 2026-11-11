//! Tests of the bar charts: the pixels of small charts built by hand from
//! the layout of the module documentation.

use super::*;

/// The small chart's size: 60 × 50 pixels, a plot of 20 × 10.
const W: usize = 60;
const H: usize = 50;
const P: usize = W - 2 * MARGIN_PX as usize;
const Q: usize = H - 2 * MARGIN_PX as usize;
const M: usize = MARGIN_PX as usize;

/// x from 0 to 20 with lines every 5, y from 0 to 1 with lines every 0.5.
const AXES: Axes = Axes {
    x_min: 0.0,
    x_max: 20.0,
    x_step: 5.0,
    y_max: 1.0,
    y_step: 0.5,
};

const RED: Rgb8 = [200, 0, 0];
const BLUE: Rgb8 = [0, 0, 200];
const GREEN: Rgb8 = [0, 160, 0];
const GREY: Rgb8 = [90, 90, 90];
const ORANGE: Rgb8 = [250, 140, 0];
const PURPLE: Rgb8 = [120, 0, 160];
const FIRST: Rgb8 = [10, 20, 30];
const SECOND: Rgb8 = [30, 20, 10];
const MARK: Rgb8 = [1, 2, 3];

/// The image under construction, addressed by plot column and row.
struct Image(Vec<Rgb8>);

impl Image {
    /// A blank image of `W × H` with the grid and the axes of [`AXES`]:
    /// vertical lines at 0, 5, 10, 15 and 20 (the last column), horizontal
    /// lines at 0.5 (row 10 − 5) and 1 (row 0), the axes below the plot
    /// (image row 30, columns 19 to 39) and left of it (image column 19,
    /// rows 20 to 30).
    fn blank() -> Self {
        let mut image = Self(vec![BACKGROUND; W * H]);
        for c in [0, 5, 10, 15, 19] {
            for r in 0..Q {
                image.set(c, r, GRID);
            }
        }
        for r in [5, 0] {
            for c in 0..P {
                image.set(c, r, GRID);
            }
        }
        for x in M - 1..M + P {
            image.0[(M + Q) * W + x] = AXIS;
        }
        for y in M..=M + Q {
            image.0[y * W + M - 1] = AXIS;
        }
        image
    }

    /// Sets the pixel of plot column `c` and plot row `r`.
    fn set(&mut self, c: usize, r: usize, colour: Rgb8) {
        self.0[(M + r) * W + M + c] = colour;
    }

    /// A bar in column `c`, `rows` high.
    fn bar(&mut self, c: usize, rows: usize, colour: Rgb8) {
        for r in Q - rows..Q {
            self.set(c, r, colour);
        }
    }
}

/// `pixels` against `expected`, naming the first differing pixel.
#[track_caller]
fn assert_pixels(pixels: &[Rgb8], expected: &[Rgb8], width: usize) {
    assert_eq!(pixels.len(), expected.len());
    if let Some(i) = (0..pixels.len()).find(|&i| pixels[i] != expected[i]) {
        panic!(
            "pixel ({}, {}): {:?} against {:?}",
            i % width,
            i / width,
            pixels[i],
            expected[i]
        );
    }
}

fn bar(lower: f64, upper: f64, value: f64, colour: Rgb8) -> Bar {
    Bar {
        lower,
        upper,
        value,
        colour,
    }
}

#[test]
fn an_empty_chart_has_the_grid_and_the_axes() {
    let pixels = histogram(&[], &AXES, &[], W as u32, H as u32);
    assert_pixels(&pixels, &Image::blank().0, W);
}

#[test]
fn bars_by_hand() {
    let bars = [
        // Columns 2 and 3, 0.3 · 10 rows.
        bar(2.0, 4.0, 0.3, RED),
        // Columns 4 and 5 (it starts inside column 4), round(2.5) = 3 rows;
        // in column 5 the taller green bar wins.
        bar(4.5, 6.0, 0.25, BLUE),
        bar(5.0, 7.0, 0.8, GREEN),
        // Too small for a row: one row.
        bar(10.0, 11.0, 0.001, GREY),
        // Taller than the axis: the whole plot.
        bar(12.0, 13.0, 2.0, ORANGE),
        // Not drawn: height 0 or below, outside the range.
        bar(14.0, 15.0, 0.0, RED),
        bar(15.0, 16.0, -1.0, RED),
        bar(25.0, 30.0, 0.5, RED),
        // Partly outside: column 0 only.
        bar(-3.0, 1.0, 0.5, PURPLE),
        // Equal heights: the first in column 16, the second alone in 17.
        bar(16.0, 17.0, 0.4, FIRST),
        bar(16.0, 18.0, 0.4, SECOND),
    ];
    let pixels = histogram(&bars, &AXES, &[], W as u32, H as u32);
    let mut expected = Image::blank();
    for (c, rows, colour) in [
        (2, 3, RED),
        (3, 3, RED),
        (4, 3, BLUE),
        (5, 8, GREEN),
        (6, 8, GREEN),
        (10, 1, GREY),
        (12, 10, ORANGE),
        (0, 5, PURPLE),
        (16, 4, FIRST),
        (17, 4, SECOND),
    ] {
        expected.bar(c, rows, colour);
    }
    assert_pixels(&pixels, &expected.0, W);
}

#[test]
fn markers_are_drawn_over_the_bars() {
    let bars = [bar(5.0, 7.0, 0.8, GREEN)];
    let markers = [
        // Columns 7 and 6 (over the bar), the last column for x_max; none
        // outside the range.
        Marker {
            x: 7.5,
            colour: MARK,
        },
        Marker {
            x: 6.0,
            colour: MARK,
        },
        Marker {
            x: 20.0,
            colour: MARK,
        },
        Marker {
            x: -1.0,
            colour: RED,
        },
        Marker {
            x: 21.0,
            colour: RED,
        },
    ];
    let pixels = histogram(&bars, &AXES, &markers, W as u32, H as u32);
    let mut expected = Image::blank();
    expected.bar(5, 8, GREEN);
    expected.bar(6, 8, GREEN);
    for c in [7, 6, 19] {
        expected.bar(c, Q, MARK);
    }
    assert_pixels(&pixels, &expected.0, W);
}

#[test]
fn bars_narrower_than_a_column_show_the_tallest() {
    // 1,000 bars of width 1 on 10 columns of 100: column c shows the bar
    // 100 c + 99, the tallest of its hundred, round((100 c + 99)/1000 · 10)
    // = c + 1 rows high.
    let colour = |k: usize| [(k % 256) as u8, (k / 256) as u8, 7];
    let bars: Vec<Bar> = (0..1_000)
        .map(|k| bar(k as f64, k as f64 + 1.0, k as f64 / 1_000.0, colour(k)))
        .collect();
    let axes = Axes {
        x_min: 0.0,
        x_max: 1_000.0,
        x_step: 500.0,
        y_max: 1.0,
        y_step: 1.0,
    };
    let (w, h) = (2 * M + 10, 2 * M + 10);
    let pixels = histogram(&bars, &axes, &[], w as u32, h as u32);
    for c in 0..10 {
        for r in 0..10 {
            let pixel = pixels[(M + r) * w + M + c];
            if r >= 10 - (c + 1) {
                assert_eq!(pixel, colour(100 * c + 99), "column {c}, row {r}");
            } else {
                let line = r == 0 || c == 0 || c == 5 || c == 9;
                assert_eq!(
                    pixel,
                    if line { GRID } else { BACKGROUND },
                    "column {c}, row {r}"
                );
            }
        }
    }
}

#[test]
fn bins_of_100_m_on_1000_columns() {
    // The layout of `hypsometry.png`: −11 km to +9 km on 1,000 columns of
    // 20 m; the bin from −7,400 m to −7,300 m fills columns 180 to 184.
    let axes = Axes {
        x_min: -11_000.0,
        x_max: 9_000.0,
        x_step: 1_000.0,
        y_max: 0.25,
        y_step: 0.05,
    };
    let (w, h) = (1_040, 540);
    let pixels = histogram(&[bar(-7_400.0, -7_300.0, 0.25, RED)], &axes, &[], w, h);
    let bottom = (M + 499) * w as usize;
    let columns: Vec<usize> = (0..1_000)
        .filter(|&c| pixels[bottom + M + c] == RED)
        .collect();
    assert_eq!(columns, (180..185).collect::<Vec<_>>());
    // The full height of 500 rows, up to the top line.
    assert_eq!(pixels[M * w as usize + M + 180], RED);
    // Grid lines every 1,000 m (every 50 columns), every 0.05 (100 rows).
    assert_eq!(pixels[(M + 250) * w as usize + M + 50], GRID);
    assert_eq!(pixels[(M + 400) * w as usize + M + 3], GRID);
    assert_eq!(pixels[(M + 250) * w as usize + M + 51], BACKGROUND);
}

#[test]
fn bounds_whose_difference_overflows() {
    // From −1.7e308 to +1.7e308 (a difference beyond the largest double):
    // a bar from −1e308 to 0 covers columns ⌊(0.35/1.7) · 20⌋ = 4 to 9, a
    // marker at 0 lies in column 10.
    let axes = Axes {
        x_min: -1.7e308,
        x_max: 1.7e308,
        x_step: 1e308,
        y_max: 1.0,
        y_step: 0.5,
    };
    let bars = [bar(-1e308, 0.0, 0.5, RED)];
    let markers = [Marker {
        x: 0.0,
        colour: MARK,
    }];
    let pixels = histogram(&bars, &axes, &markers, W as u32, H as u32);
    let at = |c: usize, r: usize| pixels[(M + r) * W + M + c];
    let red: Vec<usize> = (0..P).filter(|&c| at(c, Q - 1) == RED).collect();
    assert_eq!(red, (4..10).collect::<Vec<_>>());
    assert!((0..Q).all(|r| at(10, r) == MARK));
    // Five rows high; column 4 also holds the grid line at −1e308.
    assert_eq!((at(5, Q - 5), at(5, Q - 6)), (RED, BACKGROUND));
    assert_eq!(at(4, Q - 6), GRID);
}

#[test]
fn subnormal_bounds() {
    // From 0 to 4 · 5e-324 on 20 columns: a bar from 5e-324 to 3 · 5e-324
    // covers columns 5 to 14, a marker at 2 · 5e-324 lies in column 10.
    let tiny = 5e-324;
    let axes = Axes {
        x_min: 0.0,
        x_max: 4.0 * tiny,
        x_step: 4.0 * tiny,
        y_max: 1.0,
        y_step: 0.5,
    };
    let bars = [bar(tiny, 3.0 * tiny, 0.5, RED)];
    let markers = [Marker {
        x: 2.0 * tiny,
        colour: MARK,
    }];
    let pixels = histogram(&bars, &axes, &markers, W as u32, H as u32);
    let at = |c: usize, r: usize| pixels[(M + r) * W + M + c];
    let red: Vec<usize> = (0..P).filter(|&c| at(c, Q - 1) == RED).collect();
    assert_eq!(red, [5, 6, 7, 8, 9, 11, 12, 13, 14]);
    assert!((0..Q).all(|r| at(10, r) == MARK));
}

#[test]
#[should_panic(expected = "larger than its margins")]
fn needs_room_for_a_plot() {
    let _ = histogram(&[], &AXES, &[], 40, 50);
}

#[test]
#[should_panic(expected = "x_min < x_max")]
fn needs_an_x_range() {
    let axes = Axes { x_max: 0.0, ..AXES };
    let _ = histogram(&[], &axes, &[], W as u32, H as u32);
}

#[test]
#[should_panic(expected = "positive y_max")]
fn needs_a_y_range() {
    let axes = Axes { y_max: 0.0, ..AXES };
    let _ = histogram(&[], &axes, &[], W as u32, H as u32);
}

#[test]
#[should_panic(expected = "positive steps")]
fn needs_positive_steps() {
    let axes = Axes {
        x_step: f64::NAN,
        ..AXES
    };
    let _ = histogram(&[], &axes, &[], W as u32, H as u32);
}

#[test]
#[should_panic(expected = "more grid lines than")]
fn needs_few_grid_lines() {
    let axes = Axes {
        y_step: 0.01,
        ..AXES
    };
    let _ = histogram(&[], &axes, &[], W as u32, H as u32);
}

#[test]
#[should_panic(expected = "not finite")]
fn needs_finite_bars() {
    let _ = histogram(
        &[bar(1.0, 2.0, f64::INFINITY, RED)],
        &AXES,
        &[],
        W as u32,
        H as u32,
    );
}

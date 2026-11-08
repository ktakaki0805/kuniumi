//! Tests of the colour maps. Reference values come from the sources through
//! `tools/colormap_tables.py`; std's `powf` and `cbrt` serve as references
//! for the colour science (sRGB and CIELAB), which the module does not use.

use super::*;

/// Viridis at `t = j/20`, interpolated linearly through all 256 entries of
/// the BIDS table and rounded half up (`tools/colormap_tables.py`).
const VIRIDIS_REFERENCE: [(u32, Rgb8); 21] = [
    (0, [68, 1, 84]),
    (1, [72, 20, 102]),
    (2, [72, 37, 117]),
    (3, [70, 52, 128]),
    (4, [65, 68, 135]),
    (5, [59, 82, 139]),
    (6, [52, 95, 141]),
    (7, [47, 108, 142]),
    (8, [42, 120, 142]),
    (9, [37, 132, 142]),
    (10, [33, 144, 140]),
    (11, [30, 156, 137]),
    (12, [34, 168, 132]),
    (13, [47, 180, 124]),
    (14, [67, 191, 113]),
    (15, [93, 201, 99]),
    (16, [122, 209, 81]),
    (17, [154, 217, 61]),
    (18, [188, 223, 39]),
    (19, [221, 227, 24]),
    (20, [253, 231, 37]),
];

/// ColorBrewer Paired, 12 classes (`tools/colormap_tables.py`).
const PAIRED: [Rgb8; 12] = [
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

fn max_channel_difference(a: Rgb8, b: Rgb8) -> i32 {
    (0..3)
        .map(|c| (i32::from(a[c]) - i32::from(b[c])).abs())
        .max()
        .unwrap()
}

/// sRGB channel to linear light (IEC 61966-2-1).
fn linear(c: u8) -> f64 {
    let c = f64::from(c) / 255.0;
    if c <= 0.04045 {
        c / 12.92
    } else {
        ((c + 0.055) / 1.055).powf(2.4)
    }
}

/// Relative luminance `Y` of an sRGB colour (ITU-R BT.709 primaries).
fn luminance(c: Rgb8) -> f64 {
    0.2126 * linear(c[0]) + 0.7152 * linear(c[1]) + 0.0722 * linear(c[2])
}

/// CIELAB of an sRGB colour, D65 white (CIE 15:2004).
fn lab(c: Rgb8) -> [f64; 3] {
    let (r, g, b) = (linear(c[0]), linear(c[1]), linear(c[2]));
    let x = (0.4124 * r + 0.3576 * g + 0.1805 * b) / 0.95047;
    let y = 0.2126 * r + 0.7152 * g + 0.0722 * b;
    let z = (0.0193 * r + 0.1192 * g + 0.9505 * b) / 1.08883;
    let d = 6.0 / 29.0;
    let f = |t: f64| {
        if t > d * d * d {
            t.cbrt()
        } else {
            t / (3.0 * d * d) + 4.0 / 29.0
        }
    };
    let (fx, fy, fz) = (f(x), f(y), f(z));
    [116.0 * fy - 16.0, 500.0 * (fx - fy), 200.0 * (fy - fz)]
}

fn delta_e(a: Rgb8, b: Rgb8) -> f64 {
    let (a, b) = (lab(a), lab(b));
    ((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]))
        .sqrt()
}

#[test]
fn viridis_matches_the_full_table() {
    for (j, expected) in VIRIDIS_REFERENCE {
        let t = f64::from(j) / 20.0;
        let got = viridis(t);
        assert!(
            max_channel_difference(got, expected) <= 1,
            "viridis({t}) = {got:?}, the full table gives {expected:?}"
        );
    }
    // The published end points, #440154 and #FDE725, exactly.
    assert_eq!(viridis(0.0), [0x44, 0x01, 0x54]);
    assert_eq!(viridis(1.0), [0xFD, 0xE7, 0x25]);
}

#[test]
fn viridis_luminance_increases() {
    let mut previous = luminance(viridis(0.0));
    for k in 1..=51 {
        let y = luminance(viridis(f64::from(k) / 51.0));
        assert!(y > previous, "luminance falls at t = {k}/51");
        previous = y;
    }
}

#[test]
fn maps_clamp_their_argument() {
    assert_eq!(viridis(-0.5), viridis(0.0));
    assert_eq!(viridis(7.0), viridis(1.0));
    assert_eq!(diverging(-1.0), diverging(0.0));
    assert_eq!(diverging(2.0), diverging(1.0));
}

#[test]
fn diverging_has_its_sources_colours() {
    assert_eq!(diverging(0.0), [5, 48, 97]);
    assert_eq!(diverging(0.1), [33, 102, 172]);
    assert_eq!(diverging(0.5), [247, 247, 247]);
    assert_eq!(diverging(0.9), [178, 24, 43]);
    assert_eq!(diverging(1.0), [103, 0, 31]);
    // Between the 5th and 6th control points, 3/8 of the way:
    // 209 + 0.375·38, 229 + 0.375·18, 240 + 0.375·7, rounded.
    assert_eq!(diverging(0.4375), [223, 236, 243]);
}

#[test]
fn diverging_is_symmetric_about_its_centre() {
    let map = Colormap::Diverging { centre: 0.0 };
    // m = max(0 − (−1), 3 − 0) = 3, so v ↦ 0.5 + v/6.
    assert_eq!(map.colour(0.0, -1.0, 3.0), [247, 247, 247]);
    assert_eq!(map.colour(-3.0, -1.0, 3.0), diverging(0.0));
    assert_eq!(map.colour(3.0, -1.0, 3.0), diverging(1.0));
    for v in [0.25, 0.75, 1.5] {
        assert_eq!(map.colour(v, -1.0, 3.0), diverging(0.5 + v / 6.0));
        assert_eq!(map.colour(-v, -1.0, 3.0), diverging(0.5 - v / 6.0));
    }
    // Another centre, and a field that is the centre everywhere.
    let map = Colormap::Diverging { centre: 10.0 };
    assert_eq!(map.colour(10.0, 8.0, 11.0), [247, 247, 247]);
    assert_eq!(map.colour(8.0, 8.0, 11.0), diverging(0.5 - 2.0 / 4.0));
    assert_eq!(map.colour(10.0, 10.0, 10.0), [247, 247, 247]);
}

#[test]
fn elevation_breaks_at_sea_level() {
    for sea_level in [0.0, 120.0] {
        let map = Colormap::Elevation { sea_level };
        let (min, max) = (-6_000.0, 4_000.0);
        assert_eq!(map.colour(min, min, max), ELEVATION_OCEAN[0]);
        assert_eq!(map.colour(sea_level, min, max), ELEVATION_LAND[0]);
        assert_eq!(map.colour(max, min, max), *ELEVATION_LAND.last().unwrap());
        let shallow = map.colour(sea_level - 1e-6, min, max);
        assert_eq!(shallow, *ELEVATION_OCEAN.last().unwrap());
        assert!(max_channel_difference(shallow, ELEVATION_LAND[0]) > 100);
        for k in 0..=1_000 {
            let v = min + (max - min) * f64::from(k) / 1_000.0;
            let [r, g, b] = map.colour(v, min, max);
            if v < sea_level {
                assert!(b > r && b > g, "{v} m below sea level is not blue");
            } else {
                assert!(b < r.max(g), "{v} m above sea level is blue");
            }
        }
    }
}

#[test]
fn elevation_handles_one_sided_ranges() {
    let map = Colormap::Elevation { sea_level: 0.0 };
    // Only land: the land scale runs from sea level to the maximum.
    assert_eq!(map.colour(0.0, 0.0, 0.0), ELEVATION_LAND[0]);
    assert_eq!(
        map.colour(500.0, 100.0, 500.0),
        *ELEVATION_LAND.last().unwrap()
    );
    // Only ocean: the ocean scale runs from the minimum to sea level.
    assert_eq!(map.colour(-10.0, -10.0, -5.0), ELEVATION_OCEAN[0]);
    // Values beyond an explicit range are clamped.
    assert_eq!(map.colour(-50.0, -10.0, 10.0), ELEVATION_OCEAN[0]);
    assert_eq!(
        map.colour(50.0, -10.0, 10.0),
        *ELEVATION_LAND.last().unwrap()
    );
}

#[test]
fn categorical_uses_paired() {
    assert_eq!(CATEGORICAL, PAIRED);
    for class in 0..40 {
        assert_eq!(categorical(class), PAIRED[(class % 12) as usize]);
    }
    let map = Colormap::Categorical;
    assert_eq!(map.colour(13.0, 0.0, 0.0), PAIRED[1]);
    assert_eq!(map.colour(-0.0, 5.0, 9.0), PAIRED[0]);
    assert_eq!(map.colour(9_007_199_254_740_991.0, 0.0, 0.0), PAIRED[7]);
    for bad in [1.5, -1.0, 9_007_199_254_740_992.0, 1e300] {
        assert_eq!(map.colour(bad, 0.0, 0.0), NO_DATA, "class {bad}");
    }
}

#[test]
fn categorical_colours_are_distinct() {
    let mut closest = f64::INFINITY;
    for (i, &a) in CATEGORICAL.iter().enumerate() {
        assert_ne!(a, NO_DATA);
        for &b in &CATEGORICAL[i + 1..] {
            closest = closest.min(delta_e(a, b));
        }
    }
    // 23.64 by tools/colormap_tables.py.
    assert!(
        (closest - 23.64).abs() < 0.01,
        "closest pair ΔE = {closest}"
    );
    assert!(closest >= 20.0);
}

#[test]
fn non_finite_values_have_no_data_colour() {
    let maps = [
        Colormap::Viridis,
        Colormap::Diverging { centre: 0.0 },
        Colormap::Elevation { sea_level: 0.0 },
        Colormap::Categorical,
    ];
    for map in maps {
        for v in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert_eq!(map.colour(v, -1.0, 1.0), NO_DATA, "{map:?} of {v}");
        }
    }
}

#[test]
fn overflowing_ranges_give_no_data() {
    // (1e308 + 1e308)/(1e308 + 1e308) is ∞/∞.
    assert_eq!(Colormap::Viridis.colour(1e308, -1e308, 1e308), NO_DATA);
    // Without overflow, the largest values still have colours.
    assert_eq!(Colormap::Viridis.colour(1e308, 0.0, 1e308), viridis(1.0));
}

#[test]
fn sequential_normalisation() {
    let map = Colormap::Viridis;
    assert_eq!(map.colour(2.0, 2.0, 6.0), viridis(0.0));
    assert_eq!(map.colour(3.0, 2.0, 6.0), viridis(0.25));
    assert_eq!(map.colour(6.0, 2.0, 6.0), viridis(1.0));
    assert_eq!(map.colour(-9.0, 2.0, 6.0), viridis(0.0));
    assert_eq!(map.colour(9.0, 2.0, 6.0), viridis(1.0));
    // A constant field takes the middle colour.
    assert_eq!(map.colour(5.0, 5.0, 5.0), viridis(0.5));
}

#[test]
fn auto_range_uses_the_finite_values() {
    let values = [3.0, f64::NAN, -1.0, f64::INFINITY, 2.0, f64::NEG_INFINITY];
    assert_eq!(Range::Auto.resolve(&values), (-1.0, 3.0));
    assert_eq!(Range::Auto.resolve(&[]), (0.0, 0.0));
    assert_eq!(Range::Auto.resolve(&[f64::NAN]), (0.0, 0.0));
    let fixed = Range::Fixed {
        min: -2.0,
        max: 8.0,
    };
    assert_eq!(fixed.resolve(&values), (-2.0, 8.0));
}

#[test]
#[should_panic(expected = "range")]
fn reversed_fixed_range_panics() {
    let _ = Range::Fixed { min: 1.0, max: 0.0 }.resolve(&[]);
}

#[test]
#[should_panic(expected = "range")]
fn non_finite_fixed_range_panics() {
    let _ = Range::Fixed {
        min: f64::NAN,
        max: 0.0,
    }
    .resolve(&[]);
}

#[test]
fn colours_maps_every_value() {
    let values = [0.0, 1.0, f64::NAN, 2.0];
    let colours = Colormap::Viridis.colours(&values, Range::Auto);
    assert_eq!(
        colours,
        [viridis(0.0), viridis(0.5), NO_DATA, viridis(1.0)].to_vec()
    );
    let fixed = Colormap::Viridis.colours(&values, Range::Fixed { min: 0.0, max: 4.0 });
    assert_eq!(
        fixed,
        [viridis(0.0), viridis(0.25), NO_DATA, viridis(0.5)].to_vec()
    );
}

#[test]
#[should_panic(expected = "finite level")]
fn non_finite_sea_level_panics() {
    let _ = Colormap::Elevation {
        sea_level: f64::NAN,
    }
    .colour(0.0, -1.0, 1.0);
}

#[test]
#[should_panic(expected = "finite level")]
fn non_finite_centre_panics() {
    let _ = Colormap::Diverging {
        centre: f64::INFINITY,
    }
    .colours(&[0.0], Range::Auto);
}

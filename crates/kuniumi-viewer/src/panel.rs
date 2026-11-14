//! The panel over the globe (ticket MV-06): the text of the cell under the
//! cursor ([`inspect::cell_text`]) and the legend of the layer shown as a
//! colour bar with its range, laid out as quads of the bitmap font
//! ([`font`]) on a translucent background in the top-left corner. Nothing
//! here needs a GPU; `gpu.rs` draws the quads.
//!
//! Lengths are in font pixels, multiplied by the integer text scale `s`
//! ([`text_scale`]) to give physical pixels, so each font pixel covers
//! `s × s` whole pixels and the glyphs stay sharp.

use kuniumi_io::Rgb8;
use kuniumi_io::colormap::Colormap;
use kuniumi_world::World;

use crate::font::{self, GLYPH_PX};
use crate::inspect;
use crate::layers::{Layer, Layers};

/// The panel's distance from the viewport's top-left corner (font px).
pub const MARGIN_PX: u32 = 8;
/// The space between the background's edge and its content (font px).
pub const PADDING_PX: u32 = 6;
/// The pitch of text lines (font px): 8 px glyphs and 2 px between them.
pub const LINE_PX: u32 = 10;
/// The space from the last text row to the colour bar and from the bar to
/// its labels (font px).
pub const BAR_GAP_PX: u32 = 4;
/// The height of the colour bar (font px).
pub const BAR_PX: u32 = 12;
/// The fewest columns a line is wrapped at, however narrow the viewport.
pub const MIN_COLUMNS: usize = 16;
/// The most columns of a line before it is wrapped.
pub const MAX_COLUMNS: usize = 64;
/// The indent of the continuation of a wrapped line (columns).
pub const INDENT: usize = 4;
/// The largest text scale.
pub const MAX_TEXT_SCALE: u32 = 8;
/// The colour of the text (sRGB and alpha).
pub const TEXT_COLOUR: [u8; 4] = [230, 230, 230, 255];
/// The colour of the background (sRGB and alpha): black, translucent.
pub const BACKGROUND_COLOUR: [u8; 4] = [0, 0, 0, 176];

/// A rectangle of the panel in physical pixels from the viewport's
/// top-left corner, filled with a glyph of the atlas stretched over it and
/// drawn in a colour.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Quad {
    /// The left edge (px).
    pub x_px: u32,
    /// The top edge (px).
    pub y_px: u32,
    /// The width (px).
    pub width_px: u32,
    /// The height (px).
    pub height_px: u32,
    /// The glyph ([`font::glyph`], or [`font::SOLID`] for a filled
    /// rectangle).
    pub glyph: u32,
    /// The colour of the glyph's set pixels: sRGB and alpha.
    pub colour: [u8; 4],
}

/// The text scale for a window whose scale factor (physical pixels per
/// logical pixel) is `scale_factor`: rounded half away from zero, clamped
/// to [1, [`MAX_TEXT_SCALE`]]; 1 for NaN.
pub fn text_scale(scale_factor: f64) -> u32 {
    let rounded = scale_factor.round();
    if rounded >= f64::from(MAX_TEXT_SCALE) {
        MAX_TEXT_SCALE
    } else if rounded >= 1.0 {
        // An integer from 1 to MAX_TEXT_SCALE − 1.
        rounded as u32
    } else {
        // Below 1, or NaN.
        1
    }
}

/// The lines of the panel: the text of `picked` if a cell is under the
/// cursor, else [`inspect::NO_CELL`] and the line of the layer shown.
pub fn lines(world: &World, layers: &Layers, picked: Option<u32>) -> Vec<String> {
    match picked {
        Some(cell) => inspect::cell_text(world, layers, cell),
        None => {
            let mut lines = vec![inspect::NO_CELL.to_string()];
            lines.extend(inspect::layer_line(layers));
            lines
        }
    }
}

/// `line` broken into lines of at most `columns` characters (at least
/// [`INDENT`] + 1): a longer line breaks at its last space that leaves at
/// most `columns` characters and some that are not spaces before it, the
/// spaces at the break dropped, else after `columns` characters; the
/// continuation starts with [`INDENT`] spaces and breaks the same way.
pub fn wrap(line: &str, columns: usize) -> Vec<String> {
    let columns = columns.max(INDENT + 1);
    let chars: Vec<char> = line.chars().collect();
    let mut wrapped = Vec::new();
    let mut rest = &chars[..];
    let mut indent = 0;
    loop {
        let room = columns - indent;
        let pad = " ".repeat(indent);
        if rest.len() <= room {
            wrapped.push(pad + &String::from_iter(rest));
            return wrapped;
        }
        // The last space within reach that has text before it.
        let space = (1..=room)
            .rev()
            .find(|&k| rest[k] == ' ' && rest[..k].iter().any(|&c| c != ' '));
        let (part, next) = match space {
            Some(k) => {
                let end = rest[..k]
                    .iter()
                    .rposition(|&c| c != ' ')
                    .map_or(0, |e| e + 1);
                let skip = rest[k..].iter().take_while(|&&c| c == ' ').count();
                (&rest[..end], k + skip)
            }
            None => (&rest[..room], room),
        };
        wrapped.push(pad + &String::from_iter(part));
        rest = &rest[next..];
        if rest.is_empty() {
            // Only spaces followed the break.
            return wrapped;
        }
        indent = INDENT;
    }
}

/// The colours of a bar of `n` columns for `layer`, from its minimum on the
/// left to its maximum on the right: for a continuous map, column `i` has
/// the colour of the value `min + (i + ½)/n · (max − min)`; for the
/// categorical map, of the class `min + ⌊(i + ½)·k/n⌋` with `k = max − min
/// + 1` classes, so every class from `min` to `max` gets an equal share.
pub fn bar_colours(layer: &Layer, n: usize) -> Vec<Rgb8> {
    let (map, min, max) = (layer.colormap, layer.min, layer.max);
    let n_f = n as f64;
    (0..n)
        .map(|i| {
            let middle = i as f64 + 0.5;
            let value = match map {
                Colormap::Categorical => min + (middle * (max - min + 1.0) / n_f).floor(),
                _ => min + middle / n_f * (max - min),
            };
            map.colour(value, min, max)
        })
        .collect()
}

/// The quads of the panel showing `lines` and, with a layer, its colour bar
/// and range, at the text scale `scale` (0 counts as 1) in a viewport of
/// `viewport_px` (width, height) physical pixels, in drawing order: the
/// background, the bar's columns, then one quad per character that is not
/// a space.
///
/// Lines are wrapped ([`wrap`]) at the columns that fit the viewport's
/// width, at least [`MIN_COLUMNS`] and at most [`MAX_COLUMNS`]. The bar
/// spans the panel's text width, one column per font pixel; under it the
/// minimum is written at the left and the maximum with the layer's unit at
/// the right.
pub fn layout(
    lines: &[String],
    layer: Option<&Layer>,
    scale: u32,
    viewport_px: (u32, u32),
) -> Vec<Quad> {
    let scale = scale.max(1);
    let inset = MARGIN_PX + PADDING_PX;
    let fit = ((viewport_px.0 / scale).saturating_sub(2 * inset) / GLYPH_PX) as usize;
    let limit = fit.clamp(MIN_COLUMNS, MAX_COLUMNS);
    let mut text: Vec<String> = lines.iter().flat_map(|line| wrap(line, limit)).collect();
    let length = |line: &String| line.chars().count();
    let labels = layer.map(|layer| {
        let left = inspect::format_field_value(layer.kind, layer.min);
        let mut right = inspect::format_field_value(layer.kind, layer.max);
        if !layer.unit.is_empty() {
            right = format!("{right} {}", layer.unit);
        }
        (left, right)
    });
    let shortest_labels = labels
        .as_ref()
        .map_or(0, |(left, right)| length(left) + 1 + length(right));
    let columns = text
        .iter()
        .map(length)
        .max()
        .unwrap_or(0)
        .max(shortest_labels);
    let rows = text.len() as u32;
    let text_y = |row: u32| inset + LINE_PX * row;
    // The bottom of the last row's glyphs.
    let mut bottom = text_y(rows.saturating_sub(1)) + GLYPH_PX;
    let mut quads = Vec::new();
    let mut label_y = None;
    let mut bar = Vec::new();
    if let (Some(layer), Some((left, right))) = (layer, &labels) {
        let bar_y = bottom + BAR_GAP_PX;
        let n = GLYPH_PX as usize * columns;
        bar = (0..)
            .zip(bar_colours(layer, n))
            .map(|(i, [r, g, b])| Quad {
                x_px: (inset + i) * scale,
                y_px: bar_y * scale,
                width_px: scale,
                height_px: BAR_PX * scale,
                glyph: font::SOLID,
                colour: [r, g, b, 255],
            })
            .collect();
        let y = bar_y + BAR_PX + BAR_GAP_PX;
        label_y = Some(y);
        bottom = y + GLYPH_PX;
        let gap = columns - length(left) - length(right);
        text.push(format!("{left}{}{right}", " ".repeat(gap)));
    }
    quads.push(Quad {
        x_px: MARGIN_PX * scale,
        y_px: MARGIN_PX * scale,
        width_px: (2 * PADDING_PX + GLYPH_PX * columns as u32) * scale,
        height_px: (bottom + PADDING_PX - MARGIN_PX) * scale,
        glyph: font::SOLID,
        colour: BACKGROUND_COLOUR,
    });
    quads.extend(bar);
    for (row, line) in (0..).zip(&text) {
        let y = if row < rows {
            text_y(row)
        } else {
            label_y.expect("the label line follows a bar")
        };
        for (column, c) in (0..).zip(line.chars()) {
            if c != ' ' {
                quads.push(Quad {
                    x_px: (inset + GLYPH_PX * column) * scale,
                    y_px: y * scale,
                    width_px: GLYPH_PX * scale,
                    height_px: GLYPH_PX * scale,
                    glyph: font::glyph(c),
                    colour: TEXT_COLOUR,
                });
            }
        }
    }
    quads
}

/// The bytes of one vertex: position (px) and texel coordinate as two
/// `f32` each, then the colour as four bytes.
pub const VERTEX_BYTES: usize = 20;

/// The vertices of each quad: two triangles.
pub const VERTICES_PER_QUAD: usize = 6;

/// The vertex buffer of `quads`, in the machine's byte order: per quad the
/// corners `(x0, y0) (x1, y0) (x0, y1)` and `(x1, y0) (x1, y1) (x0, y1)`,
/// each with its position in pixels (`Float32x2`), the matching corner of
/// the quad's glyph in the atlas in texels (`Float32x2`) and the quad's
/// colour (`Unorm8x4`).
pub fn vertices(quads: &[Quad]) -> Vec<u8> {
    let mut bytes = Vec::with_capacity(quads.len() * VERTICES_PER_QUAD * VERTEX_BYTES);
    for quad in quads {
        let (x0, y0) = (quad.x_px, quad.y_px);
        let (x1, y1) = (x0 + quad.width_px, y0 + quad.height_px);
        let [u0, v0] = font::glyph_origin(quad.glyph);
        let (u1, v1) = (u0 + GLYPH_PX, v0 + GLYPH_PX);
        let corners = [
            (x0, y0, u0, v0),
            (x1, y0, u1, v0),
            (x0, y1, u0, v1),
            (x1, y0, u1, v0),
            (x1, y1, u1, v1),
            (x0, y1, u0, v1),
        ];
        for (x, y, u, v) in corners {
            // Pixel and texel positions are small integers, exact in f32.
            for value in [x, y, u, v] {
                bytes.extend((value as f32).to_ne_bytes());
            }
            bytes.extend(quad.colour);
        }
    }
    bytes
}

#[cfg(test)]
mod tests {
    use kuniumi_io::colormap::{CATEGORICAL, viridis};
    use kuniumi_world::{FieldKind, PlanetParams, Settings};

    use super::*;

    fn layer(colormap: Colormap, min: f64, max: f64, unit: &str) -> Layer {
        Layer {
            name: "f".into(),
            unit: unit.into(),
            kind: match colormap {
                Colormap::Categorical => FieldKind::Class,
                _ => FieldKind::Quantity,
            },
            colormap,
            min,
            max,
        }
    }

    fn strings(lines: &[&str]) -> Vec<String> {
        lines.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn text_scales_by_hand() {
        let cases = [
            (0.5, 1),
            (1.0, 1),
            (1.25, 1),
            (1.5, 2),
            (1.75, 2),
            (2.0, 2),
            (2.49, 2),
            (3.0, 3),
            (100.0, MAX_TEXT_SCALE),
            (f64::INFINITY, MAX_TEXT_SCALE),
            (-2.0, 1),
            (f64::NAN, 1),
        ];
        for (factor, scale) in cases {
            assert_eq!(text_scale(factor), scale, "{factor}");
        }
    }

    #[test]
    fn wrapping_by_hand() {
        let cases: [(&str, usize, &[&str]); 10] = [
            ("abc", 10, &["abc"]),
            ("abcdefghij", 10, &["abcdefghij"]),
            ("abcde fghij klm", 10, &["abcde", "    fghij", "    klm"]),
            ("abcdefghij klm", 10, &["abcdefghij", "    klm"]),
            (
                "abcdefghijklmnopqrstuvwxyz",
                10,
                &["abcdefghij", "    klmnop", "    qrstuv", "    wxyz"],
            ),
            // The spaces at a break are dropped; leading spaces stay.
            ("abcd     efgh", 6, &["abcd", "    ef", "    gh"]),
            ("  ab cdefghijkl", 10, &["  ab", "    cdefgh", "    ijkl"]),
            ("abcd      ", 6, &["abcd"]),
            // A break that would leave only spaces is not taken.
            ("  abcdefghijkl", 10, &["  abcdefgh", "    ijkl"]),
            // Characters, not bytes.
            ("ééééé ééééé", 9, &["ééééé", "    ééééé"]),
        ];
        for (line, columns, expected) in cases {
            assert_eq!(wrap(line, columns), strings(expected), "{line:?}");
            for part in wrap(line, columns) {
                assert!(part.chars().count() <= columns, "{part:?}");
            }
        }
        // At least one column after the indent.
        assert_eq!(wrap("abcdefg", 2), strings(&["abcde", "    f", "    g"]));
        assert_eq!(wrap("", 10), strings(&[""]));
    }

    #[test]
    fn the_lines_with_and_without_a_cell() {
        let settings = Settings {
            cells: 1_000,
            allow_partial: true,
            ..Settings::default()
        };
        let world =
            kuniumi_sim::generate(&PlanetParams::earth(), &settings, 3, &mut |_| {}).unwrap();
        let layers = Layers::new(&world);
        assert_eq!(
            lines(&world, &layers, Some(517)),
            inspect::cell_text(&world, &layers, 517)
        );
        let none = lines(&world, &layers, None);
        assert_eq!(none.len(), 2);
        assert_eq!(none[0], "no cell under the cursor");
        assert_eq!(
            Some(&none[1]),
            inspect::cell_text(&world, &layers, 0).last()
        );
        let mut empty = world;
        empty.fields.clear();
        let layers = Layers::new(&empty);
        assert_eq!(lines(&empty, &layers, None), ["no cell under the cursor"]);
    }

    #[test]
    fn layout_by_hand() {
        let text = strings(&["ab c", "x"]);
        let bar = layer(Colormap::Viridis, 0.0, 4.0, "m");
        for s in [1, 2] {
            let quads = layout(&text, Some(&bar), s, (1280, 800));
            // Five columns: the label line "0 4 m" is the widest.
            let columns = 5;
            let quad = |x, y, w, h, glyph, colour| Quad {
                x_px: x * s,
                y_px: y * s,
                width_px: w * s,
                height_px: h * s,
                glyph,
                colour,
            };
            // Text rows at 8 + 6 = 14 and 24; the bar 4 below the second
            // row's glyphs, at 36, 12 high; the labels at 52; the content
            // ends at 60 and the background 6 below.
            let background = quad(
                8,
                8,
                6 + 8 * columns + 6,
                58,
                font::SOLID,
                BACKGROUND_COLOUR,
            );
            assert_eq!(quads[0], background);
            let n = 8 * columns as usize;
            let colours = bar_colours(&bar, n);
            for (i, colour) in colours.iter().enumerate() {
                let [r, g, b] = *colour;
                assert_eq!(
                    quads[1 + i],
                    quad(14 + i as u32, 36, 1, 12, font::SOLID, [r, g, b, 255]),
                    "column {i}"
                );
            }
            let glyph = |c, x, y| quad(x, y, 8, 8, font::glyph(c), TEXT_COLOUR);
            assert_eq!(
                quads[1 + n..],
                [
                    glyph('a', 14, 14),
                    glyph('b', 22, 14),
                    glyph('c', 38, 14),
                    glyph('x', 14, 24),
                    glyph('0', 14, 52),
                    glyph('4', 30, 52),
                    glyph('m', 46, 52),
                ]
            );
        }
        // Text wider than the labels: the maximum ends where the bar ends
        // (review R1-1). Eight columns; "0" at column 0, "4 m" at 5 to 7.
        let wide = layout(&strings(&["abcdefgh"]), Some(&bar), 1, (1280, 800));
        assert_eq!(wide[0].width_px, 6 + 8 * 8 + 6);
        let labels: Vec<(u32, u32, u32)> = wide[wide.len() - 3..]
            .iter()
            .map(|q| (q.glyph, q.x_px, q.x_px + q.width_px))
            .collect();
        assert_eq!(
            labels,
            [
                (font::glyph('0'), 14, 22),
                (font::glyph('4'), 54, 62),
                (font::glyph('m'), 70, 78),
            ]
        );
        let bar_end = wide[1..1 + 64].iter().map(|q| q.x_px + q.width_px).max();
        assert_eq!(bar_end, Some(14 + 8 * 8));
        // A scale of 0 counts as 1.
        assert_eq!(
            layout(&text, Some(&bar), 0, (1280, 800)),
            layout(&text, Some(&bar), 1, (1280, 800))
        );
        // Without a layer: no bar, no labels; the content ends with the
        // last row's glyphs at 24 + 8.
        let quads = layout(&text, None, 1, (1280, 800));
        assert_eq!(quads.len(), 1 + 4);
        assert_eq!(
            (
                quads[0].x_px,
                quads[0].y_px,
                quads[0].width_px,
                quads[0].height_px
            ),
            (8, 8, 6 + 8 * 4 + 6, 32 + 6 - 8)
        );
        // A class layer: integer labels, no unit, and a long unit widens
        // the panel.
        let classes = layer(Colormap::Categorical, 0.0, 11.0, "");
        let quads = layout(&strings(&["a"]), Some(&classes), 1, (1280, 800));
        let labels: Vec<u32> = quads.iter().rev().take(3).map(|q| q.glyph).collect();
        assert_eq!(
            labels,
            [font::glyph('1'), font::glyph('1'), font::glyph('0')]
        );
        assert_eq!(quads[0].width_px, 6 + 8 * 4 + 6, "\"0 11\"");
        let rate = layer(Colormap::Viridis, -0.5, 2.25, "m/s");
        let quads = layout(&strings(&["a"]), Some(&rate), 1, (1280, 800));
        assert_eq!(quads[0].width_px, 6 + 8 * 13 + 6, "\"-0.5 2.25 m/s\"");
    }

    #[test]
    fn the_width_follows_the_viewport() {
        // 100 characters in words of four.
        let long = vec!["abcd ".repeat(20).trim_end().to_string()];
        // (scale, viewport width, the columns that fit): ⌊(⌊W/s⌋ − 28)/8⌋
        // clamped to [16, 64].
        let cases = [
            (1, 1280, 64),
            (1, 300, 34),
            (2, 1280, 64),
            (2, 600, 34),
            (1, 50, 16),
        ];
        for (s, width, limit) in cases {
            let quads = layout(&long, None, s, (width, 800));
            let wrapped = wrap(&long[0], limit);
            let columns = wrapped.iter().map(|l| l.chars().count()).max().unwrap();
            assert!(columns <= limit && columns + 5 > limit, "{s} {width}");
            assert_eq!(quads[0].width_px, (12 + 8 * columns as u32) * s);
            let glyphs = &quads[1..];
            assert_eq!(glyphs.len(), 80, "one quad per letter");
            let right = glyphs.iter().map(|q| q.x_px + q.width_px).max().unwrap();
            assert_eq!(right, (14 + 8 * columns as u32) * s);
            let rows = wrapped.len() as u32;
            let bottom = glyphs.iter().map(|q| q.y_px + q.height_px).max().unwrap();
            assert_eq!(bottom, (14 + 10 * (rows - 1) + 8) * s);
        }
    }

    #[test]
    fn the_colour_bar_by_hand() {
        // Viridis over [0, 4], four columns: the values 0.5, 1.5, 2.5, 3.5.
        let bar = layer(Colormap::Viridis, 0.0, 4.0, "m");
        assert_eq!(
            bar_colours(&bar, 4),
            [
                viridis(0.125),
                viridis(0.375),
                viridis(0.625),
                viridis(0.875)
            ]
        );
        // Classes 3 to 6 over eight columns: two columns each.
        let classes = layer(Colormap::Categorical, 3.0, 6.0, "");
        let c = |k: usize| CATEGORICAL[k];
        assert_eq!(
            bar_colours(&classes, 8),
            [c(3), c(3), c(4), c(4), c(5), c(5), c(6), c(6)]
        );
        // One class: every column.
        let one = layer(Colormap::Categorical, 7.0, 7.0, "");
        assert_eq!(bar_colours(&one, 3), [c(7); 3]);
        // Elevation over [−3, 1] with the sea level at 0, four columns: the
        // values −2.5, −1.5, −0.5 below it and 0.5 above. By hand from the
        // tables of kuniumi_io::colormap: below the sea level the four
        // ocean colours at t = (v + 3)/3, above it the five land colours at
        // t = v/1, interpolated linearly and rounded half away from zero.
        let elevation = layer(Colormap::Elevation { sea_level: 0.0 }, -3.0, 1.0, "m");
        assert_eq!(
            bar_colours(&elevation, 4),
            [
                [15, 47, 104],   // halfway from [8, 24, 68] to [22, 70, 140]
                [40, 99, 168],   // halfway from [22, 70, 140] to [58, 128, 196]
                [104, 164, 216], // halfway from [58, 128, 196] to [150, 200, 235]
                [212, 190, 122], // the third land colour
            ]
        );
        // A range of one value: the colour of that value everywhere.
        let flat = layer(Colormap::Viridis, 2.0, 2.0, "m");
        assert_eq!(bar_colours(&flat, 2), [viridis(0.5); 2]);
    }

    #[test]
    fn vertices_by_hand() {
        let quad = Quad {
            x_px: 10,
            y_px: 20,
            width_px: 16,
            height_px: 8,
            glyph: 17,
            colour: [1, 2, 3, 4],
        };
        let bytes = vertices(&[quad, quad]);
        assert_eq!(bytes.len(), 2 * VERTICES_PER_QUAD * VERTEX_BYTES);
        let float = |k: usize| f32::from_ne_bytes(bytes[k..k + 4].try_into().unwrap());
        // Glyph 17 has its top-left texel at (8, 8).
        let expected: [[f32; 4]; 6] = [
            [10.0, 20.0, 8.0, 8.0],
            [26.0, 20.0, 16.0, 8.0],
            [10.0, 28.0, 8.0, 16.0],
            [26.0, 20.0, 16.0, 8.0],
            [26.0, 28.0, 16.0, 16.0],
            [10.0, 28.0, 8.0, 16.0],
        ];
        for (v, corner) in expected.iter().enumerate() {
            let at = v * VERTEX_BYTES;
            let got = [float(at), float(at + 4), float(at + 8), float(at + 12)];
            assert_eq!(
                got.map(f32::to_bits),
                corner.map(f32::to_bits),
                "vertex {v}"
            );
            assert_eq!(bytes[at + 16..at + 20], [1, 2, 3, 4]);
        }
        assert_eq!(bytes[..120], bytes[120..]);
        assert!(vertices(&[]).is_empty());
    }
}

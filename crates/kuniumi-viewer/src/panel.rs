//! The panel over the globe (ticket MV-06): the text of the cell under the
//! cursor ([`inspect::cell_text`]) and the legend of the layer shown as a
//! colour bar with its range, laid out as quads of the bitmap font
//! ([`font`]) on a translucent background in the top-left corner. Its width
//! follows the viewport's width; its height is bounded by the viewport's
//! height (ticket MV-07): the field lines that do not fit are cut around the
//! line of the layer shown and replaced by lines saying how many are hidden.
//! Nothing here needs a GPU; `gpu.rs` draws the quads.
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

/// The text of the panel in three parts (ticket MV-07): the lines always
/// shown at the top, the field lines, which are cut in a short viewport
/// around the line of the layer shown, and the lines always shown under
/// them.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Lines {
    /// Always shown, at the top: the cell's line (and
    /// [`inspect::NO_FIELDS`] for a world without fields), or
    /// [`inspect::NO_CELL`].
    pub head: Vec<String>,
    /// The field lines, in the world's order; cut in a short viewport.
    pub fields: Vec<String>,
    /// The index in `fields` of the line kept visible: the field of the
    /// layer shown. Clamped to the last line.
    pub shown: usize,
    /// Always shown, under the field lines: the layer line.
    pub tail: Vec<String>,
}

/// The lines of the panel: the text of `picked` if a cell is under the
/// cursor ([`inspect::cell_text`] in its three parts), else
/// [`inspect::NO_CELL`] and the line of the layer shown.
pub fn lines(world: &World, layers: &Layers, picked: Option<u32>) -> Lines {
    let tail = inspect::layer_line(layers).into_iter().collect();
    match picked {
        Some(cell) => {
            let mut head = vec![inspect::cell_line(world, cell)];
            let fields = inspect::field_lines(world, cell);
            if fields.is_empty() {
                head.push(inspect::NO_FIELDS.into());
            }
            Lines {
                head,
                fields,
                // The layers are the fields in the world's order.
                shown: layers.current_index(),
                tail,
            }
        }
        None => Lines {
            head: vec![inspect::NO_CELL.into()],
            fields: Vec::new(),
            shown: 0,
            tail,
        },
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

/// The height under the last text row's glyphs taken by a layer's colour
/// bar and its labels (font px): a gap, the bar, a gap and a row of glyphs.
pub const BAR_ROWS_PX: u32 = BAR_GAP_PX + BAR_PX + BAR_GAP_PX + GLYPH_PX;

/// The text rows that fit a viewport `height_px` physical pixels tall at
/// the text scale `scale` (0 counts as 1), with a layer's colour bar and
/// labels under them or without (ticket MV-07).
///
/// With `r` rows the background's bottom edge is at `MARGIN + PADDING +
/// LINE·(r − 1) + GLYPH + bar + PADDING` font pixels, `bar` being
/// [`BAR_ROWS_PX`] or 0; it must be at most `⌊H/s⌋ − MARGIN`, the same
/// margin below the panel as above it. So `r ≤ ⌊(⌊H/s⌋ − 26 − bar)/10⌋`:
/// with a layer `⌊(⌊H/s⌋ − 54)/10⌋`, 14 rows for 200 px at scale 1.
pub fn rows_that_fit(height_px: u32, scale: u32, with_bar: bool) -> usize {
    let bar = if with_bar { BAR_ROWS_PX } else { 0 };
    let fixed = 2 * MARGIN_PX + 2 * PADDING_PX + GLYPH_PX + bar - LINE_PX;
    ((height_px / scale.max(1)).saturating_sub(fixed) / LINE_PX) as usize
}

/// The marker line that stands for `count` field lines hidden at `place`
/// (`above`, `below`, or `hidden` for all of them).
fn marker(count: usize, place: &str) -> String {
    let noun = if count == 1 { "field" } else { "fields" };
    format!("  ... {count} {noun} {place}")
}

/// The rows of the field lines shown in `budget` rows (ticket MV-07):
/// `fields` are the field lines, each wrapped at `columns` into its rows,
/// and `shown` (clamped to the last) the one kept visible.
///
/// If all rows fit, all of them. Otherwise a window `[a, b)` of whole lines
/// around `shown`, with a marker line ([`marker`], wrapped at `columns`)
/// before it if lines above are hidden and after it if lines below are; a
/// window costs its rows and its markers' rows. From `[shown, shown + 1)`,
/// rounds add the next line below, then the next line above, each if the
/// window then still fits, until a round adds none: neither neighbour of
/// the final window fits. If not even `[shown, shown + 1)` fits, no line
/// is shown, only the marker `... n fields hidden` if it fits.
fn visible_fields(
    fields: &[Vec<String>],
    shown: usize,
    budget: usize,
    columns: usize,
) -> Vec<String> {
    let n = fields.len();
    let rows = |a: usize, b: usize| -> usize { fields[a..b].iter().map(Vec::len).sum() };
    if rows(0, n) <= budget {
        return fields.concat();
    }
    // n ≥ 1 from here: no lines would have fitted.
    let above = |a: usize| -> Vec<String> {
        if a > 0 {
            wrap(&marker(a, "above"), columns)
        } else {
            Vec::new()
        }
    };
    let below = |b: usize| -> Vec<String> {
        if b < n {
            wrap(&marker(n - b, "below"), columns)
        } else {
            Vec::new()
        }
    };
    let cost = |a: usize, b: usize| rows(a, b) + above(a).len() + below(b).len();
    let shown = shown.min(n - 1);
    let (mut a, mut b) = (shown, shown + 1);
    if cost(a, b) > budget {
        let hidden = wrap(&marker(n, "hidden"), columns);
        return if hidden.len() <= budget {
            hidden
        } else {
            Vec::new()
        };
    }
    loop {
        let mut grew = false;
        if b < n && cost(a, b + 1) <= budget {
            b += 1;
            grew = true;
        }
        if a > 0 && cost(a - 1, b) <= budget {
            a -= 1;
            grew = true;
        }
        if !grew {
            break;
        }
    }
    let mut visible = above(a);
    visible.extend(fields[a..b].iter().flatten().cloned());
    visible.extend(below(b));
    visible
}

/// The quads of the panel showing `lines` and, with a layer, its colour bar
/// and range, at the text scale `scale` (0 counts as 1) in a viewport of
/// `viewport_px` (width, height) physical pixels, in drawing order: the
/// background, the bar's columns, then one quad per character that is not
/// a space.
///
/// Lines are wrapped ([`wrap`]) at the columns that fit the viewport's
/// width, at least [`MIN_COLUMNS`] and at most [`MAX_COLUMNS`]. The head
/// and the tail of `lines` are always shown; the field lines fill the rows
/// left of those that fit the viewport's height ([`rows_that_fit`]), cut
/// around the shown one when they do not all fit, with lines saying how
/// many are hidden (ticket MV-07). The panel is as wide as its longest
/// visible row. The bar spans the panel's text width, one column per font
/// pixel; under it the minimum is written at the left and the maximum with
/// the layer's unit at the right.
///
/// The panel ends at least [`MARGIN_PX`] above the viewport's bottom edge
/// whenever the head and the tail fit (with no field line and no marker at
/// the least height); in a shorter viewport they are still drawn and run
/// past it.
pub fn layout(
    lines: &Lines,
    layer: Option<&Layer>,
    scale: u32,
    viewport_px: (u32, u32),
) -> Vec<Quad> {
    let scale = scale.max(1);
    let inset = MARGIN_PX + PADDING_PX;
    let fit = ((viewport_px.0 / scale).saturating_sub(2 * inset) / GLYPH_PX) as usize;
    let limit = fit.clamp(MIN_COLUMNS, MAX_COLUMNS);
    let wrap_all = |lines: &[String]| -> Vec<String> {
        lines.iter().flat_map(|line| wrap(line, limit)).collect()
    };
    let head = wrap_all(&lines.head);
    let tail = wrap_all(&lines.tail);
    let fields: Vec<Vec<String>> = lines.fields.iter().map(|line| wrap(line, limit)).collect();
    let budget = rows_that_fit(viewport_px.1, scale, layer.is_some())
        .saturating_sub(head.len() + tail.len());
    let mut text = head;
    text.extend(visible_fields(&fields, lines.shown, budget, limit));
    text.extend(tail);
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
    use kuniumi_world::{Field, FieldKind, PlanetParams, Settings};

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

    /// `lines` all in the head: never cut.
    fn plain(lines: &[&str]) -> Lines {
        Lines {
            head: strings(lines),
            ..Lines::default()
        }
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
        let mut layers = Layers::new(&world);
        let text = inspect::cell_text(&world, &layers, 517);
        let last = text.len() - 1;
        // The three parts of the cell's text; the layer shown at start is
        // surface_height.
        let shown = world
            .fields
            .iter()
            .position(|f| f.name == "surface_height")
            .unwrap();
        assert_eq!(
            lines(&world, &layers, Some(517)),
            Lines {
                head: text[..1].to_vec(),
                fields: text[1..last].to_vec(),
                shown,
                tail: text[last..].to_vec(),
            }
        );
        // The line kept visible follows the layer.
        layers.select("plate_id").unwrap();
        let picked = lines(&world, &layers, Some(517));
        assert_eq!(picked.shown, 2);
        assert!(picked.fields[2].starts_with("  plate_id "), "{picked:?}");
        assert!(picked.tail[0].starts_with("layer 3/"), "{picked:?}");
        // No cell: the note and the layer line, no field lines.
        let none = lines(&world, &layers, None);
        assert_eq!(none.head, ["no cell under the cursor"]);
        assert!(none.fields.is_empty());
        assert_eq!(
            Some(&none.tail[0]),
            inspect::cell_text(&world, &layers, 0).last()
        );
        assert_eq!(none.tail.len(), 1);
        // A world without fields: the cell's line and "no fields", no layer.
        let mut empty = world;
        empty.fields.clear();
        let layers = Layers::new(&empty);
        assert_eq!(
            lines(&empty, &layers, None),
            Lines {
                head: strings(&["no cell under the cursor"]),
                ..Lines::default()
            }
        );
        assert_eq!(
            lines(&empty, &layers, Some(3)),
            Lines {
                head: inspect::cell_text(&empty, &layers, 3),
                ..Lines::default()
            }
        );
        assert_eq!(lines(&empty, &layers, Some(3)).head[1], "no fields");
    }

    #[test]
    fn layout_by_hand() {
        let text = plain(&["ab c", "x"]);
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
        let wide = layout(&plain(&["abcdefgh"]), Some(&bar), 1, (1280, 800));
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
        let quads = layout(&plain(&["a"]), Some(&classes), 1, (1280, 800));
        let labels: Vec<u32> = quads.iter().rev().take(3).map(|q| q.glyph).collect();
        assert_eq!(
            labels,
            [font::glyph('1'), font::glyph('1'), font::glyph('0')]
        );
        assert_eq!(quads[0].width_px, 6 + 8 * 4 + 6, "\"0 11\"");
        let rate = layer(Colormap::Viridis, -0.5, 2.25, "m/s");
        let quads = layout(&plain(&["a"]), Some(&rate), 1, (1280, 800));
        assert_eq!(quads[0].width_px, 6 + 8 * 13 + 6, "\"-0.5 2.25 m/s\"");
    }

    #[test]
    fn the_width_follows_the_viewport() {
        // 100 characters in words of four.
        let line = "abcd ".repeat(20).trim_end().to_string();
        let long = plain(&[&line]);
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
            let wrapped = wrap(&line, limit);
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

    /// The text of `quads` read back: per row of glyphs (by their top edge,
    /// in pixels, rows in the order drawn, which must be from the top) its
    /// characters, each at the column of its left edge with spaces between,
    /// the glyphs mapped back to characters by the inverse of `font::glyph`
    /// over U+0020 to U+007E and the degree sign.
    fn read_back(quads: &[Quad], scale: u32) -> Vec<(u32, String)> {
        let mut chars = ['\0'; 97];
        for c in (' '..='~').chain(['°']) {
            chars[font::glyph(c) as usize] = c;
        }
        let mut rows: Vec<(u32, String)> = Vec::new();
        for quad in quads.iter().filter(|q| q.glyph != font::SOLID) {
            assert_eq!((quad.width_px, quad.height_px), (8 * scale, 8 * scale));
            match rows.last() {
                Some((y, _)) if *y == quad.y_px => {}
                last => {
                    assert!(last.is_none_or(|(y, _)| *y < quad.y_px), "{quad:?}");
                    rows.push((quad.y_px, String::new()));
                }
            }
            let line = &mut rows.last_mut().expect("a row was pushed").1;
            let column = ((quad.x_px - 14 * scale) / (8 * scale)) as usize;
            let length = line.chars().count();
            assert!(length <= column, "{quad:?} after {line:?}");
            line.extend(std::iter::repeat_n(' ', column - length));
            line.push(chars[quad.glyph as usize]);
        }
        rows
    }

    /// The text rows of `quads` without their positions.
    fn text_of(quads: &[Quad], scale: u32) -> Vec<String> {
        read_back(quads, scale)
            .into_iter()
            .map(|(_, r)| r)
            .collect()
    }

    /// The marker line for `count` hidden field lines at `place`.
    fn marker(count: usize, place: &str) -> String {
        let noun = if count == 1 { "field" } else { "fields" };
        format!("  ... {count} {noun} {place}")
    }

    #[test]
    fn the_height_follows_the_viewport() {
        // By hand (ticket MV-07): one head line, six one-row field lines, a
        // one-row layer line and a viridis bar over [0, 4] m. The rows that
        // fit are ⌊(⌊H/s⌋ − 54)/10⌋ (the background's bottom edge, 10·rows
        // + 46, at most H − 8); the field lines get two rows fewer, F.
        let bar = layer(Colormap::Viridis, 0.0, 4.0, "m");
        let fields = ["f0", "f1", "f2", "f3", "f4", "f5"];
        let text = |shown| Lines {
            head: strings(&["head"]),
            fields: strings(&fields),
            shown,
            tail: strings(&["tail"]),
        };
        let (a2, b2, b3) = (
            "  ... 2 fields above",
            "  ... 2 fields below",
            "  ... 3 fields below",
        );
        let hidden = "  ... 6 fields hidden";
        // (H at scale 1, shown, the field rows).
        let cases: [(u32, usize, &[&str]); 15] = [
            // F = 6: all six.
            (134, 0, &fields),
            (134, 3, &fields),
            // F = 5: four lines and a marker, the window around the shown
            // line grown below first, then above.
            (124, 0, &["f0", "f1", "f2", "f3", b2]),
            (124, 5, &[a2, "f2", "f3", "f4", "f5"]),
            (124, 2, &["f0", "f1", "f2", "f3", b2]),
            (124, 3, &[a2, "f2", "f3", "f4", "f5"]),
            // F = 4: below first, then above (above first would give
            // `... 1 field above`, f1, f2, `... 3 fields below`).
            (114, 4, &["  ... 3 fields above", "f3", "f4", "f5"]),
            (114, 2, &[a2, "f2", "f3", "  ... 2 fields below"]),
            // F = 3: markers on both sides.
            (104, 2, &[a2, "f2", b3]),
            // F = 2: the shown line and its two markers do not fit.
            (94, 2, &[hidden]),
            (94, 0, &["f0", "  ... 5 fields below"]),
            // F = 1.
            (84, 2, &[hidden]),
            (84, 0, &[hidden]),
            // F = 0: no field row, no marker; the panel ends at H − 8.
            (74, 2, &[]),
            (74, 0, &[]),
        ];
        for s in [1, 2] {
            for (height, shown, expected) in cases {
                let quads = layout(&text(shown), Some(&bar), s, (1280, height * s));
                let rows = text_of(&quads, s);
                let (label, rows) = rows.split_last().unwrap();
                let mut want = vec!["head"];
                want.extend(expected);
                want.push("tail");
                assert_eq!(rows, want, "H = {height}, shown {shown}, scale {s}");
                assert!(label.starts_with('0') && label.ends_with("4 m"), "{label}");
                let bottom = quads[0].y_px + quads[0].height_px;
                assert!(bottom <= (height - 8) * s, "H = {height}: {bottom}");
                if height == 134 || height == 74 {
                    assert_eq!(bottom, (height - 8) * s, "H = {height}");
                }
                if height == 134 {
                    // As in a tall viewport.
                    assert_eq!(quads, layout(&text(shown), Some(&bar), s, (1280, 800)));
                }
            }
        }
        // A field line of two rows hidden alone: its marker takes one row,
        // so it says "1 field" (a one-row line is never hidden alone).
        let long = Lines {
            head: strings(&["head"]),
            fields: vec!["f0".into(), "f1".into(), format!("f2{}", "x".repeat(70))],
            shown: 0,
            tail: strings(&["tail"]),
        };
        assert_eq!(wrap(&long.fields[2], 64).len(), 2);
        let quads = layout(&long, Some(&bar), 1, (1280, 104));
        let rows = text_of(&quads, 1);
        assert_eq!(
            rows[..rows.len() - 1],
            ["head", "f0", "f1", "  ... 1 field below", "tail"]
        );
        // The hidden line does not widen the panel: 19 columns, those of
        // the marker, not 64.
        assert_eq!(quads[0].width_px, 6 + 8 * 19 + 6);
        // `shown` past the end counts as the last line.
        assert_eq!(
            layout(&text(17), Some(&bar), 1, (1280, 124)),
            layout(&text(5), Some(&bar), 1, (1280, 124))
        );
        // Without a layer there is no bar: the rows that fit are ⌊(⌊H/s⌋ −
        // 26)/10⌋. H = 76: five rows, three of them for fields.
        let rows = text_of(&layout(&text(0), None, 1, (1280, 76)), 1);
        assert_eq!(rows, ["head", "f0", "f1", "  ... 4 fields below", "tail"]);
        // Below the minimum the head, the layer line, the bar and the labels
        // are still drawn and run past the bottom margin: H = 64 has room
        // for one row, F < 0.
        let quads = layout(&text(2), Some(&bar), 1, (1280, 64));
        let rows = text_of(&quads, 1);
        assert_eq!(rows[..2], ["head", "tail"]);
        assert_eq!(rows.len(), 3);
        assert_eq!(quads[0].y_px + quads[0].height_px, 10 * 2 + 46);
    }

    #[test]
    fn the_panel_fits_a_short_window() {
        // A world with many fields (ticket MV-07): the 13 of a generated
        // world and 47 copies of cell_area under other names, 60 in all.
        let settings = Settings {
            cells: 1_000,
            allow_partial: true,
            ..Settings::default()
        };
        let mut world =
            kuniumi_sim::generate(&PlanetParams::earth(), &settings, 3, &mut |_| {}).unwrap();
        let area = world.fields[0].clone();
        assert_eq!(area.name, "cell_area");
        for k in 0..47 {
            world.fields.push(Field {
                name: format!("extra_{k:02}"),
                ..area.clone()
            });
        }
        let n = world.fields.len();
        assert_eq!(n, 60);
        let mut layers = Layers::new(&world);
        let mut checked = 0;
        // The layers 1, 3 (plate_id), 13 (surface_height) and 60 shown.
        for shown in [0, 2, 12, 59] {
            let name = world.fields[shown].name.clone();
            layers.select(&name).unwrap();
            let layer = layers.current().unwrap();
            let left = inspect::format_field_value(layer.kind, layer.min);
            let right = format!(
                "{} {}",
                inspect::format_field_value(layer.kind, layer.max),
                layer.unit
            );
            let right = right.trim_end().to_string();
            for picked in [Some(517), None] {
                let text = lines(&world, &layers, picked);
                if picked.is_some() {
                    assert_eq!((text.fields.len(), text.shown), (n, shown));
                } else {
                    assert!(text.fields.is_empty());
                }
                // Narrow viewports wrap the field lines: at 400 px and scale
                // 2 into 21 columns; at 180 px and scale 1 into 19, where
                // the markers wrap too (`  ... 60 fields hidden` has 22
                // characters, `  ... 47 fields above` 21).
                for (s, width) in [(1, 1280), (2, 1280), (1, 400), (2, 400), (1, 180)] {
                    // The columns: ⌊(⌊W/s⌋ − 28)/8⌋ clamped to [16, 64].
                    let limit = ((width / s - 28) / 8).clamp(16, 64) as usize;
                    // The panel's right edge: inside the viewport, unless
                    // the label line alone is longer than the column limit
                    // and widens the panel past it (MV-06: the columns are
                    // at least the shortest label line), as for cell_area's
                    // `0.00652904 0.0178412 sr`, 23 characters, at 21 and 19
                    // columns.
                    let label_columns = (left.len() + 1 + right.len()) as u32;
                    let right_px = if label_columns as usize > limit {
                        s * (MARGIN_PX + 2 * PADDING_PX + 8 * label_columns)
                    } else {
                        width
                    };
                    let wrapped = |lines: &[String]| -> Vec<String> {
                        lines.iter().flat_map(|l| wrap(l, limit)).collect()
                    };
                    let head = wrapped(&text.head);
                    let tail = wrapped(&text.tail);
                    let fields: Vec<Vec<String>> =
                        text.fields.iter().map(|l| wrap(l, limit)).collect();
                    let n = fields.len();
                    let rows_of =
                        |a: usize, b: usize| -> usize { fields[a..b].iter().map(Vec::len).sum() };
                    let above = |a: usize| -> Vec<String> {
                        if a > 0 {
                            wrap(&marker(a, "above"), limit)
                        } else {
                            Vec::new()
                        }
                    };
                    let below = |b: usize| -> Vec<String> {
                        if b < n {
                            wrap(&marker(n - b, "below"), limit)
                        } else {
                            Vec::new()
                        }
                    };
                    let cost = |a: usize, b: usize| rows_of(a, b) + above(a).len() + below(b).len();
                    let (h, l) = (head.len() as u32, tail.len() as u32);
                    let all = rows_of(0, n);
                    // The stated minimum: the head, the layer line and one
                    // marker row; up to 50 px past the height where all fit.
                    let min_px = s * (54 + 10 * (h + l + 1));
                    let all_px = s * (54 + 10 * (h + l + all as u32));
                    for height in min_px..=min_px.max(all_px) + 50 {
                        let context = format!(
                            "shown {shown}, cell {picked:?}, scale {s}, {width} x {height}"
                        );
                        let quads = layout(&text, Some(layer), s, (width, height));
                        // Inside the viewport, the margin left below.
                        for q in &quads {
                            assert!(
                                q.x_px + q.width_px <= right_px
                                    && q.y_px + q.height_px <= height - 8 * s,
                                "{context}: {q:?}"
                            );
                        }
                        let read = read_back(&quads, s);
                        let ((label_y, label), read) = read.split_last().unwrap();
                        let last_row_y = read.last().unwrap().0;
                        let rows: Vec<String> = read.iter().map(|(_, r)| r.clone()).collect();
                        // The labels last, the head first, the layer line
                        // above the labels.
                        assert!(label.starts_with(&left), "{context}: {label}");
                        assert!(label.ends_with(&right), "{context}: {label}");
                        assert_eq!(rows[..head.len()], head, "{context}");
                        assert_eq!(rows[rows.len() - tail.len()..], tail, "{context}");
                        // The bar: 8 columns per text column, opaque, between
                        // the last text row and the labels.
                        let columns = (quads[0].width_px / s - 12) / 8;
                        let bar = &quads[1..1 + 8 * columns as usize];
                        for q in bar {
                            assert_eq!(q.glyph, font::SOLID, "{context}");
                            assert_eq!((q.width_px, q.height_px), (s, 12 * s), "{context}");
                            assert_eq!(q.colour[3], 255, "{context}");
                            assert!(q.y_px >= last_row_y + 8 * s, "{context}");
                            assert!(q.y_px + q.height_px <= *label_y, "{context}");
                        }
                        assert_ne!(quads[1 + bar.len()].glyph, font::SOLID, "{context}");
                        // The field lines.
                        let middle = &rows[head.len()..rows.len() - tail.len()];
                        let fit = ((height / s).saturating_sub(54) / 10) as usize;
                        let budget = fit.saturating_sub(head.len() + tail.len());
                        if all <= budget {
                            // All of them, no marker: the cell's text.
                            assert_eq!(middle, fields.concat(), "{context}");
                            if let Some(cell) = picked {
                                let text = inspect::cell_text(&world, &layers, cell);
                                assert_eq!(rows, wrapped(&text), "{context}");
                            }
                        } else if cost(shown, shown + 1) > budget {
                            // Not even the shown line with its markers.
                            let hidden = wrap(&marker(n, "hidden"), limit);
                            if hidden.len() <= budget {
                                assert_eq!(middle, hidden, "{context}");
                            } else {
                                assert!(middle.is_empty(), "{context}");
                            }
                        } else {
                            // A window [a, b) around the shown line, with
                            // the markers of what it hides.
                            let a = (0..=shown)
                                .find(|&a| {
                                    let m = above(a);
                                    middle.starts_with(&m)
                                        && middle[m.len()..].starts_with(&fields[a])
                                })
                                .unwrap_or_else(|| panic!("{context}: {middle:?}"));
                            let mut at = above(a).len();
                            let mut b = a;
                            while b < n && middle[at..].starts_with(&fields[b]) {
                                at += fields[b].len();
                                b += 1;
                            }
                            assert!(a <= shown && shown < b, "{context}: {middle:?}");
                            assert_eq!(middle[at..], below(b), "{context}");
                            assert!(cost(a, b) <= budget, "{context}");
                            // Neither neighbour fits.
                            assert!(a == 0 || cost(a - 1, b) > budget, "{context}");
                            assert!(b == n || cost(a, b + 1) > budget, "{context}");
                        }
                        checked += 1;
                    }
                }
            }
        }
        // Each case at every height from its minimum.
        assert!(checked > 10_000, "{checked}");
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

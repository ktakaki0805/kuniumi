//! The text that describes one cell: its id, latitude and longitude, the
//! value of every field with its unit, and the legend of the layer shown
//! (ticket MV-04). Printed to the terminal on a click and with
//! `--screenshot --cursor`; MV-06 draws the same lines in the window.
//!
//! Angles are converted to degrees here, at the boundary
//! (docs/conventions.md, Units); values keep the field's SI unit.

use kuniumi_math::units::RADIANS_PER_DEGREE;
use kuniumi_sphere::coords;
use kuniumi_world::{Field, FieldKind, World};

use crate::layers::Layers;

/// The significant digits of a quantity in the text.
pub const SIGNIFICANT_DIGITS: usize = 6;

/// `x` with [`SIGNIFICANT_DIGITS`] significant digits and trailing zeros
/// trimmed: in fixed notation from 10⁻³ up to 10⁶ (after rounding), in
/// scientific notation (`1.93302e-5`) outside; `0` for either zero; a value
/// that is not finite as Rust prints it (`NaN`, `inf`, `-inf`).
pub fn format_value(x: f64) -> String {
    if !x.is_finite() {
        return format!("{x}");
    }
    if x == 0.0 {
        return "0".into();
    }
    // Rounded to the significant digits first, so that the exponent is that
    // of the rounded value (999,999.5 becomes 1e6, not 1000000).
    let digits = SIGNIFICANT_DIGITS - 1;
    let scientific = format!("{x:.digits$e}");
    let (mantissa, exponent) = scientific
        .split_once('e')
        .expect("`{:e}` writes an exponent");
    let exponent: i32 = exponent.parse().expect("`{:e}` writes an integer exponent");
    if (-3..6).contains(&exponent) {
        let decimals = (digits as i32 - exponent) as usize;
        trim_zeros(&format!("{x:.decimals$}")).into()
    } else {
        format!("{}e{exponent}", trim_zeros(mantissa))
    }
}

/// `number` without the trailing zeros of its fraction, and without the
/// point if nothing is left after it.
fn trim_zeros(number: &str) -> &str {
    if number.contains('.') {
        number.trim_end_matches('0').trim_end_matches('.')
    } else {
        number
    }
}

/// The lines that describe `cell` of `world`, with `layers` its layers:
///
/// ```text
/// cell 123 of 1000, lat -12.3456 deg, lon 98.7654 deg
///   cell_area             0.0125664 sr
///   plate_id              7
///   …
/// layer 13/13: surface_height (m): elevation, …
/// ```
///
/// The latitude and longitude are those of the cell's centre, with four
/// decimals; the field names are padded to the longest, quantities given by
/// [`format_value`] with their unit, class fields as their integer id. The
/// last line is the layer shown with its [`legend`](crate::layers::Layer::legend).
/// A world without fields has the line `no fields` after the first instead.
///
/// # Panics
///
/// If `cell` is not a cell of `world`.
pub fn cell_text(world: &World, layers: &Layers, cell: u32) -> Vec<String> {
    let (lat_rad, lon_rad) = coords::lat_lon(world.grid.position(cell));
    let mut lines = vec![format!(
        "cell {cell} of {}, lat {:.4} deg, lon {:.4} deg",
        world.grid.len(),
        lat_rad / RADIANS_PER_DEGREE,
        lon_rad / RADIANS_PER_DEGREE
    )];
    if world.fields.is_empty() {
        lines.push("no fields".into());
        return lines;
    }
    let width = world.fields.iter().map(|f| f.name.len()).max().unwrap_or(0);
    for field in &world.fields {
        let value = field.values[cell as usize];
        let shown = match field.kind {
            // Class ids are integers below 2⁵³, exact as u64.
            FieldKind::Class if Field::is_class_value(value) => format!("{}", value as u64),
            _ => format_value(value),
        };
        let unit = if field.unit.is_empty() {
            String::new()
        } else {
            format!(" {}", field.unit)
        };
        lines.push(format!("  {:<width$}  {shown}{unit}", field.name));
    }
    if let Some(layer) = layers.current() {
        lines.push(format!(
            "layer {}/{}: {}",
            layers.current_index() + 1,
            layers.len(),
            layer.legend()
        ));
    }
    lines
}

// The tests use std's `to_degrees` as an independent reference.
#[cfg(test)]
#[allow(clippy::disallowed_methods)]
mod tests {
    use kuniumi_world::{PlanetParams, Settings};

    use super::*;

    #[test]
    fn values_by_hand() {
        let cases = [
            (0.0, "0"),
            (-0.0, "0"),
            (-4321.5734, "-4321.57"),
            (1.9330183543e-5, "1.93302e-5"),
            (6.371e6, "6.371e6"),
            (123_456.0, "123456"),
            (999_999.5, "1e6"),
            (999_999.4, "999999"),
            (9_999.995_1, "10000"),
            (0.001, "0.001"),
            (0.000_999_999_9, "0.001"),
            (0.000_999_99, "9.9999e-4"),
            (-2.5e-12, "-2.5e-12"),
            (1.0e300, "1e300"),
            (0.5, "0.5"),
            (7.0, "7"),
            (f64::NAN, "NaN"),
            (f64::INFINITY, "inf"),
            (f64::NEG_INFINITY, "-inf"),
        ];
        for (x, text) in cases {
            assert_eq!(format_value(x), text, "{x:e}");
        }
    }

    #[test]
    fn the_text_of_a_cell() {
        let settings = Settings {
            cells: 1_000,
            allow_partial: true,
            ..Settings::default()
        };
        let world =
            kuniumi_sim::generate(&PlanetParams::earth(), &settings, 3, &mut |_| {}).unwrap();
        let mut layers = Layers::new(&world);
        layers.select("plate_id").unwrap();
        let width = world.fields.iter().map(|f| f.name.len()).max().unwrap();
        for cell in [0, 1, 517, 999] {
            let text = cell_text(&world, &layers, cell);
            assert_eq!(text.len(), world.fields.len() + 2);
            let (lat, lon) = coords::lat_lon(world.grid.position(cell));
            assert_eq!(
                text[0],
                format!(
                    "cell {cell} of 1000, lat {:.4} deg, lon {:.4} deg",
                    lat.to_degrees(),
                    lon.to_degrees()
                )
            );
            for (line, field) in text[1..].iter().zip(&world.fields) {
                let value = field.values[cell as usize];
                let shown = match field.kind {
                    FieldKind::Class => format!("{}", value as u64),
                    FieldKind::Quantity => format_value(value),
                };
                let unit = if field.unit.is_empty() {
                    String::new()
                } else {
                    format!(" {}", field.unit)
                };
                assert_eq!(*line, format!("  {:<width$}  {shown}{unit}", field.name));
            }
            let index = layers.current_index();
            assert_eq!(
                text[text.len() - 1],
                format!(
                    "layer {}/{}: {}",
                    index + 1,
                    layers.len(),
                    layers.current().unwrap().legend()
                )
            );
        }
        // Hand-checked pieces: a class id is an integer, a quantity has its unit.
        let text = cell_text(&world, &layers, 517);
        let plate = text.iter().find(|l| l.contains("plate_id ")).unwrap();
        assert!(
            plate.trim_end().chars().last().unwrap().is_ascii_digit(),
            "{plate}"
        );
        assert!(!plate.contains('.'), "{plate}");
        let height = text.iter().find(|l| l.contains("surface_height ")).unwrap();
        assert!(height.ends_with(" m"), "{height}");
        let last = text.last().unwrap();
        assert!(last.starts_with("layer "), "{last}");
        assert!(last.contains("plate_id (class ids 0 to "), "{last}");

        // A world without fields.
        let mut empty = world;
        empty.fields.clear();
        let text = cell_text(&empty, &Layers::new(&empty), 3);
        assert_eq!(text.len(), 2);
        assert!(text[0].starts_with("cell 3 of 1000, lat "));
        assert_eq!(text[1], "no fields");
    }
}

//! Renders test maps of a Fibonacci grid with every colour map, as
//! equirectangular maps and globes, for looking at (ticket M1-08).
//!
//! ```sh
//! cargo run --release -p kuniumi-io --example render_maps [-- OUT_DIR [CELLS [WIDTH]]]
//! ```
//!
//! Defaults: `OUT_DIR` = `out/M1-08`, `CELLS` = 10,000, `WIDTH` = 1024 (the
//! globes are `WIDTH / 2` square); grid seed 1. The fields are `z` (the sine
//! of latitude, viridis), `cell % 12` (categorical), `x` (diverging, centre
//! 0) and a synthetic elevation in metres with land and sea (sea level 0).
//! For each image the program prints the time to render it, the time to
//! encode it and the size of the PNG.

use std::path::{Path, PathBuf};
use std::time::Instant;

use kuniumi_io::colormap::{Colormap, Range};
use kuniumi_io::{Rgb8, png, render};
use kuniumi_math::det;
use kuniumi_math::units::deg_to_rad;
use kuniumi_sphere::SphereGrid;

/// A field and how to draw it.
struct Layer {
    name: &'static str,
    values: Vec<f64>,
    colormap: Colormap,
}

/// Renders with `draw`, writes `<dir>/<file>` and prints the timings.
fn write(dir: &Path, file: &str, width: u32, height: u32, draw: impl FnOnce() -> Vec<Rgb8>) {
    let start = Instant::now();
    let pixels = draw();
    let rendered = start.elapsed();
    let bytes = png::encode_rgb8(width, height, &pixels);
    let encoded = start.elapsed() - rendered;
    let path = dir.join(file);
    std::fs::write(&path, &bytes).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    println!(
        "{file:32} {width:5} x {height:<5} render {:7.3} s  encode {:6.3} s  {:9} bytes",
        rendered.as_secs_f64(),
        encoded.as_secs_f64(),
        bytes.len()
    );
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let dir = PathBuf::from(args.first().map_or("out/M1-08", String::as_str));
    let cells: u32 = args.get(1).map_or(10_000, |s| s.parse().expect("CELLS"));
    let width: u32 = args.get(2).map_or(1024, |s| s.parse().expect("WIDTH"));
    let size = width / 2;
    std::fs::create_dir_all(&dir).expect("output directory");

    let start = Instant::now();
    let grid = SphereGrid::fibonacci(cells, 1).expect("grid");
    println!(
        "grid of {cells} cells in {:.3} s",
        start.elapsed().as_secs_f64()
    );

    let positions = grid.positions();
    let layers = [
        Layer {
            name: "z",
            values: positions.iter().map(|p| p.z).collect(),
            colormap: Colormap::Viridis,
        },
        Layer {
            name: "classes",
            values: (0..grid.len()).map(|i| (i % 12) as f64).collect(),
            colormap: Colormap::Categorical,
        },
        Layer {
            name: "x_diverging",
            values: positions.iter().map(|p| p.x).collect(),
            colormap: Colormap::Diverging { centre: 0.0 },
        },
        Layer {
            name: "elevation",
            values: positions
                .iter()
                .map(|p| {
                    4_000.0 * (p.x * p.y + 0.6 * p.z * p.z * p.z - 0.4 * p.x) - 500.0
                        + 800.0 * det::sin(7.0 * p.x + 5.0 * p.y) * det::sin(6.0 * p.z)
                })
                .collect(),
            colormap: Colormap::Elevation { sea_level: 0.0 },
        },
    ];

    for layer in &layers {
        let (values, colormap) = (&layer.values, layer.colormap);
        write(
            &dir,
            &format!("{}_equirect.png", layer.name),
            width,
            width / 2,
            || render::equirect(&grid, values, colormap, Range::Auto, width),
        );
        let (lat, lon) = (deg_to_rad(30.0), deg_to_rad(20.0));
        write(
            &dir,
            &format!("{}_globe.png", layer.name),
            size,
            size,
            || render::orthographic(&grid, values, colormap, Range::Auto, size, lat, lon),
        );
    }
    // The poles and the antimeridian, where projections tend to break.
    let classes = &layers[1];
    for (file, lat, lon) in [
        ("classes_north_pole.png", 90.0, 0.0),
        ("classes_south_pole.png", -90.0, 0.0),
        ("classes_antimeridian.png", 0.0, 180.0),
    ] {
        write(&dir, file, size, size, || {
            render::orthographic(
                &grid,
                &classes.values,
                classes.colormap,
                Range::Auto,
                size,
                deg_to_rad(lat),
                deg_to_rad(lon),
            )
        });
    }
}

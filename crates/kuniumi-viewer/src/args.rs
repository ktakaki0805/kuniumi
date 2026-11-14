//! The command line of `kuniumi-viewer`, parsed by hand (a handful of flags,
//! so no parser crate, D-10).

use std::path::PathBuf;

use kuniumi_math::units::deg_to_rad;

use crate::view::RELIEF_STEPS;

/// The usage text printed by `--help`.
pub const USAGE: &str = "\
Usage: kuniumi-viewer <world dir> [options]

Shows a world directory written by `kuniumi generate` as a 3D globe.

Options:
  --layer <name>       field shown first (default: surface_height if present,
                       else the first field)
  --relief <factor>    relief exaggeration, from 0 to 100 (default: 20)
  --no-sea             start without the sea-level surface
  --lat-deg <deg>      latitude under the camera (default: 20)
  --lon-deg <deg>      longitude under the camera (default: 0)
  --distance <radii>   camera distance from the centre in planet radii,
                       above 1 (default: 3)
  --size <w>x<h>       window or screenshot size in pixels (default: 1280x800)
  --screenshot <png>   render one frame without a window into a PNG and exit
  --cursor <x>,<y>     with --screenshot: pick the cell under that pixel,
                       highlight it, mark the pixel and print the cell
  --benchmark <frames> spin the globe in the window for that many frames,
                       print the times between presented frames and exit
  --offscreen          with --benchmark: draw without a window instead and
                       time each frame until the GPU has finished it
  -h, --help           print this text

Keys: Left/Right layer, Up/Down relief, S sea, R reset view, Esc or Q quit.
Mouse: drag to rotate, wheel to zoom; the cell under the cursor is
highlighted, and a click prints it.";

/// What the viewer does after loading the world.
#[derive(Clone, Debug, PartialEq)]
pub enum Mode {
    /// An interactive window.
    Window,
    /// One frame rendered without a window into the PNG file at this path.
    Screenshot(PathBuf),
    /// Spins the globe for `frames` frames and prints the frame times.
    Benchmark {
        /// The number of frames to time, at least 1.
        frames: u32,
        /// Whether the frames are drawn without a window (timed until the
        /// GPU has finished each) instead of presented in the window (timed
        /// between presents, paced by the display).
        offscreen: bool,
    },
}

/// The parsed command line.
#[derive(Clone, Debug, PartialEq)]
pub struct Args {
    /// The world directory.
    pub world_dir: PathBuf,
    /// The field shown first, if given.
    pub layer: Option<String>,
    /// The relief exaggeration, from 0 to 100.
    pub relief: f64,
    /// Whether the sea-level surface is drawn.
    pub sea: bool,
    /// The latitude under the camera (rad).
    pub lat_rad: f64,
    /// The longitude under the camera (rad).
    pub lon_rad: f64,
    /// The camera's distance from the centre, in planet radii.
    pub distance: f64,
    /// The width of the window or screenshot (px).
    pub width_px: u32,
    /// The height of the window or screenshot (px).
    pub height_px: u32,
    /// With a screenshot: the pixel (column, row from the top left) whose
    /// cell is picked and printed.
    pub cursor: Option<(u32, u32)>,
    /// What to do.
    pub mode: Mode,
}

/// The result of parsing: run the viewer, or print the usage.
#[derive(Clone, Debug, PartialEq)]
pub enum Parsed {
    /// Run with these arguments.
    Run(Args),
    /// `--help` was given.
    Help,
}

/// The default relief exaggeration.
pub const DEFAULT_RELIEF: f64 = 20.0;
/// The largest relief exaggeration, the largest key step: at ×1000, Earth's
/// trenches (7 km deep on a radius of 6,371 km) would be drawn below the
/// centre of the planet, turned inside out.
const MAX_RELIEF: f64 = RELIEF_STEPS[RELIEF_STEPS.len() - 1];
/// The default latitude under the camera (degrees).
const DEFAULT_LAT_DEG: f64 = 20.0;
/// The default camera distance (planet radii from the centre).
pub const DEFAULT_DISTANCE: f64 = 3.0;
/// The largest window or screenshot side (px).
const MAX_SIDE_PX: u32 = 16_384;

/// Parses the arguments after the program name.
///
/// # Errors
///
/// A message naming the offending argument.
pub fn parse(args: impl IntoIterator<Item = String>) -> Result<Parsed, String> {
    let mut args = args.into_iter();
    let mut world_dir = None;
    let mut layer = None;
    let mut relief = DEFAULT_RELIEF;
    let mut sea = true;
    let mut lat_deg = DEFAULT_LAT_DEG;
    let mut lon_deg = 0.0;
    let mut distance = DEFAULT_DISTANCE;
    let (mut width_px, mut height_px) = (1280, 800);
    let mut screenshot = None;
    let mut benchmark = None;
    let mut offscreen = false;
    let mut cursor = None;
    while let Some(arg) = args.next() {
        let mut value = |name: &str| args.next().ok_or(format!("{name} needs a value"));
        match arg.as_str() {
            "-h" | "--help" => return Ok(Parsed::Help),
            "--layer" => layer = Some(value("--layer")?),
            "--relief" => {
                relief = number("--relief", &value("--relief")?)?;
                if !(0.0..=MAX_RELIEF).contains(&relief) {
                    return Err(format!(
                        "--relief must lie in [0, {MAX_RELIEF}], got {relief}"
                    ));
                }
            }
            "--no-sea" => sea = false,
            "--lat-deg" => {
                lat_deg = number("--lat-deg", &value("--lat-deg")?)?;
                if !(-90.0..=90.0).contains(&lat_deg) {
                    return Err(format!("--lat-deg must lie in [-90, 90], got {lat_deg}"));
                }
            }
            "--lon-deg" => lon_deg = number("--lon-deg", &value("--lon-deg")?)?,
            "--distance" => {
                distance = number("--distance", &value("--distance")?)?;
                if distance <= 1.0 {
                    return Err(format!("--distance must be above 1, got {distance}"));
                }
            }
            "--size" => (width_px, height_px) = size(&value("--size")?)?,
            "--screenshot" => screenshot = Some(PathBuf::from(value("--screenshot")?)),
            "--benchmark" => {
                let text = value("--benchmark")?;
                benchmark = Some(match text.parse::<u32>() {
                    Ok(frames) if frames > 0 => frames,
                    _ => {
                        return Err(format!(
                            "--benchmark needs a frame count of at least 1, got `{text}`"
                        ));
                    }
                });
            }
            "--offscreen" => offscreen = true,
            "--cursor" => cursor = Some(point(&value("--cursor")?)?),
            _ if arg.starts_with('-') => return Err(format!("unknown option `{arg}`")),
            _ if world_dir.is_none() => world_dir = Some(PathBuf::from(arg)),
            _ => {
                return Err(format!(
                    "unexpected argument `{arg}`: give one world directory"
                ));
            }
        }
    }
    let world_dir = world_dir.ok_or("no world directory given")?;
    let mode = match (screenshot, benchmark) {
        (Some(_), Some(_)) => return Err("--screenshot and --benchmark exclude each other".into()),
        (Some(path), None) => Mode::Screenshot(path),
        (None, Some(frames)) => Mode::Benchmark { frames, offscreen },
        (None, None) => Mode::Window,
    };
    if offscreen && !matches!(mode, Mode::Benchmark { .. }) {
        return Err("--offscreen needs --benchmark".into());
    }
    if let Some((x, y)) = cursor {
        if !matches!(mode, Mode::Screenshot(_)) {
            return Err("--cursor needs --screenshot".into());
        }
        if x >= width_px || y >= height_px {
            return Err(format!(
                "--cursor {x},{y} lies outside the image of {width_px}x{height_px} px"
            ));
        }
    }
    Ok(Parsed::Run(Args {
        world_dir,
        layer,
        relief,
        sea,
        lat_rad: deg_to_rad(lat_deg),
        lon_rad: deg_to_rad(lon_deg),
        distance,
        width_px,
        height_px,
        cursor,
        mode,
    }))
}

/// `text` as a finite number.
fn number(name: &str, text: &str) -> Result<f64, String> {
    match text.parse::<f64>() {
        Ok(x) if x.is_finite() => Ok(x),
        _ => Err(format!("{name} needs a finite number, got `{text}`")),
    }
}

/// `text` as `<x>,<y>`, a pixel's column and row.
fn point(text: &str) -> Result<(u32, u32), String> {
    text.split_once(',')
        .and_then(|(x, y)| Some((x.parse().ok()?, y.parse().ok()?)))
        .ok_or(format!("--cursor needs <x>,<y> in pixels, got `{text}`"))
}

/// `text` as `<w>x<h>`, both sides from 1 to [`MAX_SIDE_PX`].
fn size(text: &str) -> Result<(u32, u32), String> {
    let side = |s: &str| {
        s.parse::<u32>()
            .ok()
            .filter(|v| (1..=MAX_SIDE_PX).contains(v))
    };
    text.split_once('x')
        .and_then(|(w, h)| Some((side(w)?, side(h)?)))
        .ok_or(format!(
            "--size needs <width>x<height> with sides from 1 to {MAX_SIDE_PX}, got `{text}`"
        ))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run(args: &[&str]) -> Result<Parsed, String> {
        parse(args.iter().map(|s| s.to_string()))
    }

    fn args(list: &[&str]) -> Args {
        match run(list) {
            Ok(Parsed::Run(args)) => args,
            other => panic!("{list:?} gave {other:?}"),
        }
    }

    #[test]
    fn defaults() {
        let a = args(&["out/world"]);
        assert_eq!(a.world_dir, PathBuf::from("out/world"));
        assert_eq!(a.layer, None);
        assert_eq!(a.relief.to_bits(), 20.0_f64.to_bits());
        assert!(a.sea);
        assert_eq!(a.lat_rad.to_bits(), deg_to_rad(20.0).to_bits());
        assert_eq!(a.lon_rad.to_bits(), 0.0_f64.to_bits());
        assert_eq!(a.distance.to_bits(), 3.0_f64.to_bits());
        assert_eq!((a.width_px, a.height_px), (1280, 800));
        assert_eq!(a.cursor, None);
        assert_eq!(a.mode, Mode::Window);
    }

    #[test]
    fn every_flag() {
        let a = args(&[
            "--layer",
            "plate_id",
            "--relief",
            "0",
            "--no-sea",
            "--lat-deg",
            "-45.5",
            "--lon-deg",
            "190",
            "--distance",
            "1.5",
            "--size",
            "640x480",
            "w",
            "--screenshot",
            "out/a.png",
            "--cursor",
            "639,0",
        ]);
        assert_eq!(a.world_dir, PathBuf::from("w"));
        assert_eq!(a.layer.as_deref(), Some("plate_id"));
        assert_eq!(a.relief.to_bits(), 0.0_f64.to_bits());
        assert!(!a.sea);
        assert_eq!(a.lat_rad.to_bits(), deg_to_rad(-45.5).to_bits());
        assert_eq!(a.lon_rad.to_bits(), deg_to_rad(190.0).to_bits());
        assert_eq!(a.distance.to_bits(), 1.5_f64.to_bits());
        assert_eq!((a.width_px, a.height_px), (640, 480));
        assert_eq!(a.mode, Mode::Screenshot(PathBuf::from("out/a.png")));
        assert_eq!(a.cursor, Some((639, 0)));
        let c = args(&["w", "--cursor", "0,799", "--screenshot", "a.png"]);
        assert_eq!(c.cursor, Some((0, 799)));
        // The largest relief is accepted (review R2-1).
        assert_eq!(MAX_RELIEF.to_bits(), 100.0_f64.to_bits());
        let top = args(&["w", "--relief", "100"]);
        assert_eq!(top.relief.to_bits(), 100.0_f64.to_bits());
        let b = args(&["w", "--benchmark", "600", "--offscreen"]);
        assert_eq!(
            b.mode,
            Mode::Benchmark {
                frames: 600,
                offscreen: true
            }
        );
        let c = args(&["w", "--benchmark", "1"]);
        assert_eq!(
            c.mode,
            Mode::Benchmark {
                frames: 1,
                offscreen: false
            }
        );
        assert_eq!(run(&["w", "--help"]), Ok(Parsed::Help));
        assert_eq!(run(&["-h"]), Ok(Parsed::Help));
    }

    #[test]
    fn bad_input() {
        let cases: [(&[&str], &str); 23] = [
            (&[], "no world directory given"),
            (
                &["a", "b"],
                "unexpected argument `b`: give one world directory",
            ),
            (&["w", "--fly"], "unknown option `--fly`"),
            (&["w", "--layer"], "--layer needs a value"),
            (
                &["w", "--relief", "-1"],
                "--relief must lie in [0, 100], got -1",
            ),
            (
                &["w", "--relief", "100.5"],
                "--relief must lie in [0, 100], got 100.5",
            ),
            (
                &["w", "--relief", "inf"],
                "--relief needs a finite number, got `inf`",
            ),
            (
                &["w", "--relief", "x"],
                "--relief needs a finite number, got `x`",
            ),
            (
                &["w", "--lat-deg", "91"],
                "--lat-deg must lie in [-90, 90], got 91",
            ),
            (
                &["w", "--distance", "1"],
                "--distance must be above 1, got 1",
            ),
            (
                &["w", "--size", "0x10"],
                "--size needs <width>x<height> with sides from 1 to 16384, got `0x10`",
            ),
            (
                &["w", "--size", "640"],
                "--size needs <width>x<height> with sides from 1 to 16384, got `640`",
            ),
            (
                &["w", "--size", "16385x1"],
                "--size needs <width>x<height> with sides from 1 to 16384, got `16385x1`",
            ),
            (
                &["w", "--benchmark", "0"],
                "--benchmark needs a frame count of at least 1, got `0`",
            ),
            (
                &["w", "--benchmark", "1", "--screenshot", "a.png"],
                "--screenshot and --benchmark exclude each other",
            ),
            (&["w", "--offscreen"], "--offscreen needs --benchmark"),
            (&["w", "--cursor", "1,2"], "--cursor needs --screenshot"),
            (
                &["w", "--cursor", "1,2", "--benchmark", "3"],
                "--cursor needs --screenshot",
            ),
            (
                &["w", "--screenshot", "a.png", "--cursor", "1280,0"],
                "--cursor 1280,0 lies outside the image of 1280x800 px",
            ),
            (
                &[
                    "w",
                    "--screenshot",
                    "a.png",
                    "--size",
                    "64x48",
                    "--cursor",
                    "3,48",
                ],
                "--cursor 3,48 lies outside the image of 64x48 px",
            ),
            (
                &["w", "--cursor", "3x4"],
                "--cursor needs <x>,<y> in pixels, got `3x4`",
            ),
            (
                &["w", "--cursor", "-1,4"],
                "--cursor needs <x>,<y> in pixels, got `-1,4`",
            ),
            (
                &["w", "--screenshot", "a.png", "--offscreen"],
                "--offscreen needs --benchmark",
            ),
        ];
        for (list, message) in cases {
            assert_eq!(run(list), Err(message.to_string()), "{list:?}");
        }
    }
}

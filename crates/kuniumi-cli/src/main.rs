//! `kuniumi`: command-line interface for the planet generator.
//!
//! Arguments are parsed with `clap` (derive, one subcommand per command).
//! Every error, including a usage error, exits with code 2.

use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::Instant;

use clap::{Args, CommandFactory, Parser, Subcommand};
use kuniumi_io::colormap::{Colormap, Range};
use kuniumi_io::{Rgb8, png, render, world_dir};
use kuniumi_math::units::{self, deg_to_rad};
use kuniumi_sim::SimError;
use kuniumi_world::{GENERATOR_VERSION, PlanetParams, Settings};

mod grid_report;

use grid_report::{GridStats, RunInfo};

/// kuniumi — physically based planet generator
#[derive(Debug, Parser)]
#[command(name = "kuniumi", version = GENERATOR_VERSION)]
struct Cli {
    /// The command to run; without one, the help is printed.
    #[command(subcommand)]
    command: Option<Command>,
}

/// The commands.
#[derive(Debug, Subcommand)]
enum Command {
    /// Print the generator version
    Version,
    /// Print the Earth-like default parameters (SI)
    Params,
    /// Generate a world and write it to a directory (the stages that do not
    /// exist yet are missing: writes what exists, then exits with code 2)
    Generate(GenerateArgs),
    /// Build the cell grid and write maps of its cells and statistics
    Grid(GridArgs),
    /// Write a test image to check the image pipeline
    TestImage {
        /// The PNG file to write
        #[arg(value_name = "FILE.png")]
        file: PathBuf,
    },
}

/// The options of `kuniumi generate`.
#[derive(Debug, Args)]
struct GenerateArgs {
    /// Random seed
    #[arg(long, value_name = "N")]
    seed: u64,
    /// Output directory for world.json, a .kfield file and a .png map per
    /// field, and metrics.txt and hypsometry.png once the Crust stage has
    /// run (created if missing)
    #[arg(long, value_name = "DIR")]
    out: PathBuf,
    /// Number of grid cells, 1000 to 715827884
    #[arg(long, value_name = "N", default_value_t = Settings::DEFAULT_CELLS,
          value_parser = cells_parser())]
    cells: u32,
}

/// The options of `kuniumi grid`.
#[derive(Debug, Args)]
struct GridArgs {
    /// Random seed
    #[arg(long, value_name = "N")]
    seed: u64,
    /// Output directory for cell_area.png, cell_area_globe.png,
    /// neighbors.png and stats.txt (created if missing)
    #[arg(long, value_name = "DIR")]
    out: PathBuf,
    /// Number of grid cells, 1000 to 715827884
    #[arg(long, value_name = "N", default_value_t = Settings::DEFAULT_CELLS,
          value_parser = cells_parser())]
    cells: u32,
    /// Width of the equirectangular maps in pixels, 2 to 16384 (height
    /// width/2; the globe is width/2 square)
    #[arg(long, value_name = "PX", default_value_t = 2048,
          value_parser = clap::value_parser!(u32).range(2..=i64::from(MAX_WIDTH_PX)))]
    width: u32,
}

/// The largest `--width`: a 16384 × 8192 map takes about 0.5 GB for its
/// cell indices and 0.4 GB for its pixels.
const MAX_WIDTH_PX: u32 = 16_384;

/// `--cells`: [`Settings::MIN_CELLS`] to [`Settings::MAX_CELLS`].
fn cells_parser() -> clap::builder::RangedI64ValueParser<u32> {
    clap::value_parser!(u32).range(i64::from(Settings::MIN_CELLS)..=i64::from(Settings::MAX_CELLS))
}

fn main() -> ExitCode {
    // `parse` prints usage errors and exits with code 2, and prints the
    // help or the version and exits with code 0.
    let cli = Cli::parse();
    match run(cli.command) {
        Ok(()) => ExitCode::SUCCESS,
        Err(msg) => {
            eprintln!("error: {msg}");
            ExitCode::from(2)
        }
    }
}

fn run(command: Option<Command>) -> Result<(), String> {
    match command {
        None => {
            Cli::command()
                .print_help()
                .map_err(|e| format!("cannot print the help: {e}"))?;
            Ok(())
        }
        Some(Command::Version) => {
            println!("kuniumi {GENERATOR_VERSION}");
            Ok(())
        }
        Some(Command::Params) => {
            print_params(&PlanetParams::earth());
            Ok(())
        }
        Some(Command::Generate(args)) => generate(&args),
        Some(Command::Grid(args)) => grid(&args),
        Some(Command::TestImage { file }) => test_image(&file),
    }
}

fn print_params(p: &PlanetParams) {
    println!("radius_m                   = {}", p.radius_m);
    println!("surface_gravity_m_s2       = {}", p.surface_gravity_m_s2);
    println!("rotation_period_s          = {}", p.rotation_period_s);
    println!(
        "obliquity_rad              = {} ({} deg)",
        p.obliquity_rad,
        p.obliquity_rad / units::RADIANS_PER_DEGREE
    );
    println!("orbital_period_s           = {}", p.orbital_period_s);
    println!("orbital_eccentricity       = {}", p.orbital_eccentricity);
    println!("stellar_flux_w_m2          = {}", p.stellar_flux_w_m2);
    println!("ocean_volume_m3            = {:e}", p.ocean_volume_m3);
    println!("surface_pressure_pa        = {}", p.surface_pressure_pa);
    println!(
        "continental_crust_fraction = {}",
        p.continental_crust_fraction
    );
    println!("rotation                   = {:?}", p.rotation);
    println!("tectonics                  = {:?}", p.tectonics);
}

/// `kuniumi generate`: runs the pipeline as far as it exists, writes the
/// world with [`world_dir::write`], and fails naming the first missing stage
/// if there is one.
fn generate(args: &GenerateArgs) -> Result<(), String> {
    let settings = Settings {
        cells: args.cells,
        allow_partial: true,
        ..Settings::default()
    };
    let world = kuniumi_sim::generate(&PlanetParams::earth(), &settings, args.seed, &mut |p| {
        eprintln!("[{}] {:.0}%", p.stage.name(), p.fraction * 100.0);
    })
    .map_err(|e| e.to_string())?;
    for path in world_dir::write(&args.out, &world).map_err(|e| e.to_string())? {
        println!("wrote {}", path.display());
    }
    match kuniumi_sim::next_stage(&world) {
        Some(stage) => Err(SimError::NotImplemented { stage }.to_string()),
        None => Ok(()),
    }
}

/// The view centre of `cell_area_globe.png`: 45° N, 180° E, so that the
/// north pole and the antimeridian are both in view.
const GLOBE_CENTRE_DEG: (f64, f64) = (45.0, 180.0);

/// `kuniumi grid`: builds the grid of the Grid stage and writes maps of its
/// cells and `stats.txt`.
fn grid(args: &GridArgs) -> Result<(), String> {
    let settings = Settings {
        cells: args.cells,
        ..Settings::default()
    };
    settings.validate().map_err(|e| e.to_string())?;
    let dir = &args.out;
    std::fs::create_dir_all(dir).map_err(|e| format!("{}: {e}", dir.display()))?;

    let start = Instant::now();
    let grid = kuniumi_sim::stages::grid::run(&settings, args.seed)
        .map_err(|e| format!("cannot build the grid: {e}"))?;
    let build_time_s = start.elapsed().as_secs_f64();
    eprintln!("[grid] {} cells in {build_time_s:.3} s", grid.len());

    let stats = GridStats::of(&grid);
    let mean_sr = stats.area_mean_sr();
    let relative_area: Vec<f64> = (0..grid.len() as u32)
        .map(|i| grid.area_sr(i) / mean_sr)
        .collect();
    let ring_size: Vec<f64> = (0..grid.len() as u32)
        .map(|i| grid.neighbors(i).len() as f64)
        .collect();
    let area_map = Colormap::Diverging { centre: 1.0 };
    let (width, height, size) = (args.width, args.width / 2, args.width / 2);
    let (lat_rad, lon_rad) = (
        deg_to_rad(GLOBE_CENTRE_DEG.0),
        deg_to_rad(GLOBE_CENTRE_DEG.1),
    );

    write_png(
        &dir.join("cell_area.png"),
        width,
        height,
        &render::equirect(&grid, &relative_area, area_map, Range::Auto, width),
    )?;
    write_png(
        &dir.join("cell_area_globe.png"),
        size,
        size,
        &render::orthographic(
            &grid,
            &relative_area,
            area_map,
            Range::Auto,
            size,
            lat_rad,
            lon_rad,
        ),
    )?;
    write_png(
        &dir.join("neighbors.png"),
        width,
        height,
        &render::equirect(&grid, &ring_size, Colormap::Categorical, Range::Auto, width),
    )?;

    let validation = grid.validate();
    let text = grid_report::report(
        &stats,
        &RunInfo {
            seed: args.seed,
            build_time_s,
            validation: &validation,
            radius_m: PlanetParams::earth().radius_m,
        },
    );
    let path = dir.join("stats.txt");
    std::fs::write(&path, text).map_err(|e| format!("{}: {e}", path.display()))?;
    println!("wrote {}", path.display());
    Ok(())
}

/// Encodes `pixels` as a PNG of `width × height` and writes it to `path`.
fn write_png(path: &Path, width: u32, height: u32, pixels: &[Rgb8]) -> Result<(), String> {
    std::fs::write(path, png::encode_rgb8(width, height, pixels))
        .map_err(|e| format!("{}: {e}", path.display()))?;
    println!("wrote {} ({width}x{height})", path.display());
    Ok(())
}

/// Writes a 360×180 colour gradient (longitude → red, latitude → blue) so the
/// image pipeline can be checked end to end without a generated world.
fn test_image(path: &Path) -> Result<(), String> {
    let (w, h) = (360u32, 180u32);
    let mut pixels = Vec::with_capacity((w * h) as usize);
    for y in 0..h {
        for x in 0..w {
            let r = (x * 255 / (w - 1)) as u8;
            let b = (y * 255 / (h - 1)) as u8;
            pixels.push([r, 64, b]);
        }
    }
    std::fs::write(path, png::encode_rgb8(w, h, &pixels))
        .map_err(|e| format!("{}: {e}", path.display()))?;
    println!("wrote {} ({w}x{h})", path.display());
    Ok(())
}

#[cfg(test)]
mod tests {
    use clap::error::ErrorKind;

    use super::*;

    /// Every subcommand with the arguments its help must list, written by
    /// hand: long options as `--name <VALUE>`, positionals as `<VALUE>`.
    const EXPECTED: [(&str, &[&str]); 5] = [
        ("version", &[]),
        ("params", &[]),
        ("generate", &["--seed <N>", "--out <DIR>", "--cells <N>"]),
        (
            "grid",
            &["--seed <N>", "--out <DIR>", "--cells <N>", "--width <PX>"],
        ),
        ("test-image", &["<FILE.png>"]),
    ];

    /// The text clap prints for `args`, which must end in `--help`.
    fn help(args: &[&str]) -> String {
        let err = Cli::try_parse_from(args).expect_err("--help stops parsing");
        assert_eq!(err.kind(), ErrorKind::DisplayHelp, "{args:?}");
        err.to_string()
    }

    #[test]
    fn command_is_well_formed() {
        Cli::command().debug_assert();
    }

    #[test]
    fn root_help_lists_every_command() {
        let text = help(&["kuniumi", "--help"]);
        for (name, _) in EXPECTED {
            assert!(
                text.lines()
                    .any(|line| line.trim_start().starts_with(&format!("{name} "))),
                "`{name}` missing from\n{text}"
            );
        }
        for option in ["-h, --help", "-V, --version"] {
            assert!(text.contains(option), "`{option}` missing from\n{text}");
        }
        // The hand-written list covers every subcommand (clap's generated
        // `help` is not among them).
        let mut names: Vec<_> = Cli::command()
            .get_subcommands()
            .map(|c| c.get_name().to_owned())
            .collect();
        names.sort();
        let mut expected: Vec<_> = EXPECTED.iter().map(|(n, _)| n.to_string()).collect();
        expected.sort();
        assert_eq!(names, expected);
    }

    #[test]
    fn subcommand_help_lists_every_argument() {
        for (name, arguments) in EXPECTED {
            let text = help(&["kuniumi", name, "--help"]);
            for argument in arguments.iter().chain(&["-h, --help"]) {
                assert!(
                    text.contains(argument),
                    "`{argument}` missing from `{name} --help`:\n{text}"
                );
            }
            // The hand-written list covers every argument.
            let command = Cli::command();
            let sub = command.find_subcommand(name).expect("subcommand");
            let count = sub.get_arguments().filter(|a| a.get_id() != "help").count();
            assert_eq!(count, arguments.len(), "{name}");
        }
    }

    #[test]
    fn help_shows_the_defaults_and_bounds() {
        // The bounds are written out in the doc comments of `--cells`.
        let bounds = format!("{} to {}", Settings::MIN_CELLS, Settings::MAX_CELLS);
        let default = format!("[default: {}]", Settings::DEFAULT_CELLS);
        for command in ["generate", "grid"] {
            let text = help(&["kuniumi", command, "--help"]);
            assert!(text.contains(&bounds), "{command}: {bounds}: {text}");
            assert!(text.contains(&default), "{command}: {default}: {text}");
        }
        let text = help(&["kuniumi", "grid", "--help"]);
        assert!(text.contains("[default: 2048]"), "{text}");
        let width_bounds = format!("2 to {MAX_WIDTH_PX}");
        assert!(text.contains(&width_bounds), "{width_bounds}: {text}");
    }

    #[test]
    fn parses_every_command() {
        let parse = |args: &[&str]| Cli::try_parse_from(args).map(|c| c.command);
        assert!(matches!(parse(&["kuniumi"]), Ok(None)));
        assert!(matches!(
            parse(&["kuniumi", "version"]),
            Ok(Some(Command::Version))
        ));
        assert!(matches!(
            parse(&["kuniumi", "params"]),
            Ok(Some(Command::Params))
        ));
        match parse(&["kuniumi", "generate", "--seed", "3", "--out", "x"]) {
            Ok(Some(Command::Generate(a))) => {
                assert_eq!((a.seed, a.cells), (3, 650_000));
                assert_eq!(a.out, PathBuf::from("x"));
            }
            other => panic!("{other:?}"),
        }
        match parse(&[
            "kuniumi",
            "grid",
            "--cells",
            "1000",
            "--seed",
            "18446744073709551615",
            "--out",
            "d",
            "--width",
            "64",
        ]) {
            Ok(Some(Command::Grid(a))) => {
                assert_eq!((a.seed, a.cells, a.width), (u64::MAX, 1_000, 64));
                assert_eq!(a.out, PathBuf::from("d"));
            }
            other => panic!("{other:?}"),
        }
        match parse(&["kuniumi", "grid", "--seed", "1", "--out", "d"]) {
            Ok(Some(Command::Grid(a))) => assert_eq!((a.cells, a.width), (650_000, 2048)),
            other => panic!("{other:?}"),
        }
        match parse(&["kuniumi", "test-image", "a.png"]) {
            Ok(Some(Command::TestImage { file })) => assert_eq!(file, PathBuf::from("a.png")),
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn rejects_bad_arguments() {
        let kind = |args: &[&str]| Cli::try_parse_from(args).expect_err("rejected").kind();
        assert_eq!(
            kind(&["kuniumi", "frobnicate"]),
            ErrorKind::InvalidSubcommand
        );
        assert_eq!(
            kind(&["kuniumi", "generate", "--out", "x"]),
            ErrorKind::MissingRequiredArgument
        );
        assert_eq!(
            kind(&["kuniumi", "grid", "--seed", "1"]),
            ErrorKind::MissingRequiredArgument
        );
        assert_eq!(
            kind(&["kuniumi", "test-image"]),
            ErrorKind::MissingRequiredArgument
        );
        assert_eq!(
            kind(&["kuniumi", "generate", "--seed", "-1", "--out", "x"]),
            ErrorKind::UnknownArgument
        );
        for cells in ["999", "715827885", "0"] {
            for command in ["generate", "grid"] {
                assert_eq!(
                    kind(&[
                        "kuniumi", command, "--seed", "1", "--out", "x", "--cells", cells
                    ]),
                    ErrorKind::ValueValidation,
                    "{command} --cells {cells}"
                );
            }
        }
        for width in ["0", "1", "16385", "4294967295"] {
            assert_eq!(
                kind(&[
                    "kuniumi", "grid", "--seed", "1", "--out", "x", "--width", width
                ]),
                ErrorKind::ValueValidation,
                "--width {width}"
            );
        }
        for width in ["2", "16384"] {
            assert!(
                Cli::try_parse_from([
                    "kuniumi", "grid", "--seed", "1", "--out", "x", "--width", width
                ])
                .is_ok(),
                "--width {width}"
            );
        }
        for cells in ["1000", "715827884"] {
            assert!(
                Cli::try_parse_from([
                    "kuniumi", "grid", "--seed", "1", "--out", "x", "--cells", cells
                ])
                .is_ok(),
                "--cells {cells}"
            );
        }
    }

    #[test]
    fn version_is_the_generator_version() {
        let err = Cli::try_parse_from(["kuniumi", "--version"]).expect_err("--version stops");
        assert_eq!(err.kind(), ErrorKind::DisplayVersion);
        assert_eq!(err.to_string(), format!("kuniumi {GENERATOR_VERSION}\n"));
    }
}

//! `kuniumi`: command-line interface for the planet generator.
//!
//! Arguments are parsed by hand to keep the dependency count at zero.

use std::path::PathBuf;
use std::process::ExitCode;

use kuniumi_io::png;
use kuniumi_math::units;
use kuniumi_world::{GENERATOR_VERSION, PlanetParams};

const USAGE: &str = "\
kuniumi — physically based planet generator

USAGE:
    kuniumi <COMMAND> [OPTIONS]

COMMANDS:
    version                          Print the generator version
    params                           Print the Earth-like default parameters (SI)
    generate --seed <N> --out <DIR>  Generate a world (not implemented yet)
    test-image <FILE.png>            Write a test image to check the image pipeline
    help                             Print this message
";

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match run(&args) {
        Ok(()) => ExitCode::SUCCESS,
        Err(msg) => {
            eprintln!("error: {msg}");
            ExitCode::from(2)
        }
    }
}

fn run(args: &[String]) -> Result<(), String> {
    let Some((command, rest)) = args.split_first() else {
        print!("{USAGE}");
        return Ok(());
    };
    match command.as_str() {
        "version" | "--version" | "-V" => {
            println!("kuniumi {GENERATOR_VERSION}");
            Ok(())
        }
        "params" => {
            print_params(&PlanetParams::earth());
            Ok(())
        }
        "generate" => generate(rest),
        "test-image" => test_image(rest),
        "help" | "--help" | "-h" => {
            print!("{USAGE}");
            Ok(())
        }
        other => Err(format!("unknown command `{other}`\n\n{USAGE}")),
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

fn generate(args: &[String]) -> Result<(), String> {
    let mut seed: Option<u64> = None;
    let mut out: Option<PathBuf> = None;
    let mut it = args.iter();
    while let Some(arg) = it.next() {
        match arg.as_str() {
            "--seed" => {
                let v = it.next().ok_or("--seed needs a value")?;
                seed = Some(v.parse().map_err(|_| format!("invalid seed `{v}`"))?);
            }
            "--out" => out = Some(it.next().ok_or("--out needs a value")?.into()),
            other => return Err(format!("unknown option `{other}`")),
        }
    }
    let seed = seed.ok_or("--seed is required")?;
    let _out = out.ok_or("--out is required")?;
    let world = kuniumi_sim::generate(&PlanetParams::earth(), seed, &mut |p| {
        eprintln!("[{}] {:.0}%", p.stage.name(), p.fraction * 100.0);
    })
    .map_err(|e| e.to_string())?;
    let _ = world;
    Ok(())
}

/// Writes a 360×180 colour gradient (longitude → red, latitude → blue) so the
/// image pipeline can be checked end to end without a generated world.
fn test_image(args: &[String]) -> Result<(), String> {
    let [path] = args else {
        return Err("usage: kuniumi test-image <FILE.png>".into());
    };
    let (w, h) = (360u32, 180u32);
    let mut pixels = Vec::with_capacity((w * h) as usize);
    for y in 0..h {
        for x in 0..w {
            let r = (x * 255 / (w - 1)) as u8;
            let b = (y * 255 / (h - 1)) as u8;
            pixels.push([r, 64, b]);
        }
    }
    std::fs::write(path, png::encode_rgb8(w, h, &pixels)).map_err(|e| format!("{path}: {e}"))?;
    println!("wrote {path} ({w}x{h})");
    Ok(())
}

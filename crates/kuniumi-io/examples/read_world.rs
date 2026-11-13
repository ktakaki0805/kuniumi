//! Reads a world directory back and reports the time it took, for the load
//! time of ticket MV-01; optionally writes the world again and compares the
//! files.
//!
//! ```sh
//! cargo run --release -p kuniumi-io --example read_world -- DIR [COPY_DIR]
//! ```
//!
//! Prints the time of `world_dir::read` for `DIR` (three times, as the
//! first may read from the disk and the others from the page cache), the
//! time of reading `grid.kgrid` alone, the size of the grid file and of all
//! `.kfield` files, and the world's cells, stages and fields. With
//! `COPY_DIR`, writes the world read back there with `world_dir::write`
//! and compares every file it writes with the one of the same name in
//! `DIR`, byte for byte; exits with code 1 if one differs.

use std::io::BufReader;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::Instant;

use kuniumi_io::{grid_file, world_dir};

fn main() -> ExitCode {
    let mut args = std::env::args().skip(1);
    let Some(dir) = args.next().map(PathBuf::from) else {
        eprintln!("usage: read_world DIR [COPY_DIR]");
        return ExitCode::from(2);
    };
    let copy = args.next().map(PathBuf::from);

    let mut world = None;
    for run in 1..=3 {
        let start = Instant::now();
        let read = world_dir::read(&dir).unwrap_or_else(|e| panic!("{e}"));
        println!("read {}: {:.3} s", run, start.elapsed().as_secs_f64());
        world = Some(read);
    }
    let world = world.expect("read three times");

    let grid_path = dir.join(world_dir::GRID_FILE);
    let start = Instant::now();
    let file = std::fs::File::open(&grid_path).unwrap_or_else(|e| panic!("{e}"));
    let grid = grid_file::read_grid(BufReader::new(file)).unwrap_or_else(|e| panic!("{e}"));
    println!(
        "grid.kgrid alone: {:.3} s ({} cells)",
        start.elapsed().as_secs_f64(),
        grid.len()
    );
    let size = |path: &Path| std::fs::metadata(path).map_or(0, |m| m.len());
    let fields_bytes: u64 = world
        .fields
        .iter()
        .map(|f| size(&dir.join(format!("{}.kfield", f.name))))
        .sum();
    println!(
        "bytes: grid.kgrid {}, {} .kfield files {}",
        size(&grid_path),
        world.fields.len(),
        fields_bytes
    );
    println!(
        "world: {} cells, seed {}, stages {:?}, sea level {} m",
        world.grid.len(),
        world.seed,
        world.stages,
        world.sea_level_m
    );
    for field in &world.fields {
        println!("  {} [{}] {:?}", field.name, field.unit, field.kind);
    }

    let Some(copy) = copy else {
        return ExitCode::SUCCESS;
    };
    let start = Instant::now();
    let written = world_dir::write(&copy, &world).unwrap_or_else(|e| panic!("{e}"));
    println!(
        "wrote {} files into {} in {:.3} s",
        written.len(),
        copy.display(),
        start.elapsed().as_secs_f64()
    );
    let mut differ = 0;
    for path in &written {
        let name = path.file_name().expect("a file name");
        let same = std::fs::read(path).ok() == std::fs::read(dir.join(name)).ok();
        if !same {
            differ += 1;
            println!("  differs: {}", name.to_string_lossy());
        }
    }
    println!(
        "{} of {} files identical",
        written.len() - differ,
        written.len()
    );
    if differ == 0 {
        ExitCode::SUCCESS
    } else {
        ExitCode::FAILURE
    }
}

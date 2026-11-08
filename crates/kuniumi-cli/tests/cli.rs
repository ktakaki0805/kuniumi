//! Runs the `kuniumi` binary: exit codes, output files and messages.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

/// Runs `kuniumi` with `args`.
fn kuniumi(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_kuniumi"))
        .args(args)
        .output()
        .expect("the binary runs")
}

/// A fresh, empty directory for the test `name`.
fn scratch(name: &str) -> PathBuf {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join(name);
    if dir.exists() {
        std::fs::remove_dir_all(&dir).expect("remove old scratch directory");
    }
    dir
}

fn stdout(o: &Output) -> String {
    String::from_utf8_lossy(&o.stdout).into_owned()
}

fn stderr(o: &Output) -> String {
    String::from_utf8_lossy(&o.stderr).into_owned()
}

/// The width and height in the IHDR chunk of a PNG file.
fn png_size(path: &Path) -> (u32, u32) {
    let bytes = std::fs::read(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    assert_eq!(&bytes[..8], b"\x89PNG\r\n\x1a\n", "{}", path.display());
    assert_eq!(&bytes[12..16], b"IHDR", "{}", path.display());
    let be = |b: &[u8]| u32::from_be_bytes([b[0], b[1], b[2], b[3]]);
    (be(&bytes[16..20]), be(&bytes[20..24]))
}

#[test]
fn grid_writes_maps_and_statistics() {
    let dir = scratch("grid_writes_maps_and_statistics");
    let out = dir.to_str().expect("UTF-8 path");
    let o = kuniumi(&[
        "grid", "--cells", "1000", "--seed", "1", "--out", out, "--width", "64",
    ]);
    assert!(o.status.success(), "{o:?}");
    assert_eq!(png_size(&dir.join("cell_area.png")), (64, 32));
    assert_eq!(png_size(&dir.join("cell_area_globe.png")), (32, 32));
    assert_eq!(png_size(&dir.join("neighbors.png")), (64, 32));
    let stats = std::fs::read_to_string(dir.join("stats.txt")).expect("stats.txt");
    // Values recorded for `SphereGrid::fibonacci(1_000, 1)` by ticket M1-06.
    for line in [
        "seed = 1",
        "cells = 1000",
        "validate = ok",
        "area_ratio_max_min = 2.5600",
        "ring_size_4 = 12 (1.200 %, #fb9a99)",
        "ring_size_5 = 228 (22.800 %, #e31a1c)",
        "ring_size_6 = 535 (53.500 %, #fdbf6f)",
        "ring_size_7 = 210 (21.000 %, #ff7f00)",
        "ring_size_8 = 15 (1.500 %, #cab2d6)",
        "neighbour_entries = 5988",
        "neighbour_entries_expected = 5988",
    ] {
        assert!(
            stats.lines().any(|l| l == line),
            "`{line}` missing:\n{stats}"
        );
    }
    for file in [
        "cell_area.png",
        "cell_area_globe.png",
        "neighbors.png",
        "stats.txt",
    ] {
        assert!(stdout(&o).contains(file), "{file}: {}", stdout(&o));
    }
}

#[test]
fn generate_stops_after_the_grid_with_exit_code_2() {
    let dir = scratch("generate_stops_after_the_grid");
    let out = dir.to_str().expect("UTF-8 path");
    let o = kuniumi(&["generate", "--seed", "1", "--cells", "1000", "--out", out]);
    assert_eq!(o.status.code(), Some(2), "{o:?}");
    assert_eq!(
        stderr(&o),
        "[grid] 0%\n[grid] 100%\nerror: stage `tectonics` is not implemented yet\n"
    );
}

#[test]
fn bad_arguments_exit_with_code_2() {
    for args in [
        &["grid", "--cells", "999", "--seed", "1", "--out", "x"][..],
        &["generate", "--seed", "1"],
        &["frobnicate"],
        &["test-image"],
    ] {
        let o = kuniumi(args);
        assert_eq!(o.status.code(), Some(2), "{args:?}: {o:?}");
        assert!(stderr(&o).starts_with("error:"), "{args:?}: {}", stderr(&o));
    }
}

#[test]
fn help_exits_with_code_0() {
    for args in [&[][..], &["--help"], &["-h"], &["help"]] {
        let o = kuniumi(args);
        assert!(o.status.success(), "{args:?}: {o:?}");
        let text = stdout(&o);
        for command in ["version", "params", "generate", "grid", "test-image"] {
            assert!(text.contains(command), "{args:?}: {command}: {text}");
        }
    }
    for command in ["generate", "grid", "test-image"] {
        let o = kuniumi(&["help", command]);
        assert!(o.status.success(), "{command}: {o:?}");
        assert!(
            stdout(&o).contains("Usage: kuniumi "),
            "{command}: {}",
            stdout(&o)
        );
    }
}

#[test]
fn old_commands_still_work() {
    let version = format!("kuniumi {}\n", env!("CARGO_PKG_VERSION"));
    for args in [&["version"][..], &["--version"], &["-V"]] {
        let o = kuniumi(args);
        assert!(o.status.success(), "{args:?}: {o:?}");
        assert_eq!(stdout(&o), version, "{args:?}");
    }

    let o = kuniumi(&["params"]);
    assert!(o.status.success(), "{o:?}");
    assert!(stdout(&o).starts_with("radius_m                   = 6371000\n"));

    let dir = scratch("old_commands_still_work");
    std::fs::create_dir_all(&dir).expect("scratch directory");
    let file = dir.join("test.png");
    let o = kuniumi(&["test-image", file.to_str().expect("UTF-8 path")]);
    assert!(o.status.success(), "{o:?}");
    assert_eq!(png_size(&file), (360, 180));
}

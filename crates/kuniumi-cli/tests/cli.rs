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

/// Reads a `.kfield` file.
fn read_kfield(path: &Path) -> kuniumi_world::Field {
    let bytes = std::fs::read(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    kuniumi_io::field_file::read_field(bytes.as_slice())
        .unwrap_or_else(|e| panic!("{}: {e}", path.display()))
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

/// The fields of a world stopped after the Crust stage, in the order
/// `generate` writes them.
const WORLD_FIELDS: [&str; 13] = [
    "cell_area",
    "neighbor_count",
    "plate_id",
    "plate_velocity_east",
    "plate_velocity_north",
    "crust_type",
    "boundary_type",
    "convergence_rate",
    "crust_thickness",
    "crust_age",
    "uplift_rate",
    "trench",
    "surface_height",
];

/// The files of a world stopped after the Crust stage, in the order
/// `generate` writes them: `<name>.kfield` and `<name>.png` of each field,
/// the metrics (ticket M2-06), then `world.json`.
fn world_files() -> Vec<String> {
    let mut files: Vec<String> = WORLD_FIELDS
        .iter()
        .flat_map(|f| [format!("{f}.kfield"), format!("{f}.png")])
        .collect();
    files.extend(["metrics.txt".into(), "hypsometry.png".into()]);
    files.push("world.json".into());
    files
}

#[test]
fn generate_writes_the_world_so_far_then_exits_with_code_2() {
    let dir = scratch("generate_writes_the_world_so_far");
    let out = dir.to_str().expect("UTF-8 path");
    let o = kuniumi(&["generate", "--seed", "1", "--cells", "1000", "--out", out]);
    assert_eq!(o.status.code(), Some(2), "{o:?}");
    assert_eq!(
        stderr(&o),
        "[grid] 0%\n[grid] 100%\n[tectonics] 0%\n[tectonics] 100%\n\
         [crust] 0%\n[crust] 100%\n\
         error: stage `climate` is not implemented yet\n"
    );
    let expected: String = world_files()
        .iter()
        .map(|f| format!("wrote {}\n", dir.join(f).display()))
        .collect();
    assert_eq!(stdout(&o), expected);
    for field in WORLD_FIELDS {
        let map = format!("{field}.png");
        assert_eq!(png_size(&dir.join(&map)), (2048, 1024), "{map}");
    }
    assert_eq!(png_size(&dir.join("hypsometry.png")), (1040, 540));
    let metrics = std::fs::read_to_string(dir.join("metrics.txt")).expect("metrics.txt");
    for line in ["seed = 1", "cells = 1000", "bin_width_m = 100"] {
        assert!(
            metrics.lines().any(|l| l == line),
            "`{line}` missing:\n{metrics}"
        );
    }
    let json = std::fs::read_to_string(dir.join("world.json")).expect("world.json");
    for line in [
        "  \"seed\": 1,",
        "    \"cells\": 1000,",
        "    \"plate_count\": 12,",
        "    \"allow_partial\": true",
        "    \"radius_m\": 6371000.0,",
        "    \"grid\",",
        "    \"tectonics\",",
        "    \"crust\"",
        "      \"name\": \"cell_area\",",
        "      \"unit\": \"sr\",",
        "      \"kind\": \"class\"",
        "      \"name\": \"plate_id\",",
        "      \"name\": \"plate_velocity_east\",",
        "      \"unit\": \"m/s\",",
        "      \"name\": \"crust_type\",",
        "      \"name\": \"boundary_type\",",
        "      \"name\": \"convergence_rate\",",
        "      \"name\": \"crust_thickness\",",
        "      \"unit\": \"m\",",
        "      \"name\": \"crust_age\",",
        "      \"unit\": \"s\",",
        "      \"name\": \"uplift_rate\",",
        "      \"name\": \"trench\",",
        "      \"name\": \"surface_height\",",
    ] {
        assert!(json.lines().any(|l| l == line), "`{line}` missing:\n{json}");
    }
    // The cell areas of 1,000 cells sum to 4π, as the grid's do (M1-06).
    let areas = read_kfield(&dir.join("cell_area.kfield"));
    assert_eq!(
        (areas.name.as_str(), areas.unit.as_str()),
        ("cell_area", "sr")
    );
    assert_eq!(areas.values.len(), 1_000);
    let sum: f64 = areas.values.iter().sum();
    assert!((sum - 4.0 * core::f64::consts::PI).abs() < 1e-12, "{sum}");
    // Ring sizes 4 to 8, as `kuniumi grid` reports for this grid (M1-09).
    let counts = read_kfield(&dir.join("neighbor_count.kfield"));
    let mut histogram = [0usize; 10];
    for &c in &counts.values {
        histogram[c as usize] += 1;
    }
    assert_eq!(histogram, [0, 0, 0, 0, 12, 228, 535, 210, 15, 0]);
    // Twelve plates, ids 0 to 11, each with a cell; no velocity above
    // 10 cm/yr (3.17e-9 m/s).
    let plates = read_kfield(&dir.join("plate_id.kfield"));
    let mut sizes = [0usize; 12];
    for &p in &plates.values {
        sizes[p as usize] += 1;
    }
    assert!(sizes.iter().all(|&s| s > 0), "{sizes:?}");
    let east = read_kfield(&dir.join("plate_velocity_east.kfield"));
    let north = read_kfield(&dir.join("plate_velocity_north.kfield"));
    assert_eq!((east.unit.as_str(), north.unit.as_str()), ("m/s", "m/s"));
    for (e, n) in east.values.iter().zip(&north.values) {
        let speed = (e * e + n * n).sqrt();
        assert!(speed <= 0.1 / 31_557_600.0 * (1.0 + 1e-12), "{speed}");
    }
    // Continental crust on 40 % of the area, within one percentage point
    // (Earth's `continental_crust_fraction`).
    let crust = read_kfield(&dir.join("crust_type.kfield"));
    assert_eq!(crust.unit, "");
    let (ocean, land) = (0.0_f64.to_bits(), 1.0_f64.to_bits());
    assert!(
        crust
            .values
            .iter()
            .all(|c| [ocean, land].contains(&c.to_bits()))
    );
    let land_sr: f64 = crust
        .values
        .iter()
        .zip(&areas.values)
        .map(|(c, a)| c * a)
        .sum();
    let fraction = land_sr / sum;
    assert!((fraction - 0.4).abs() <= 0.01, "{fraction}");
    // Boundary types 0 to 3, each present; no rate away from boundaries.
    let kinds = read_kfield(&dir.join("boundary_type.kfield"));
    let rates = read_kfield(&dir.join("convergence_rate.kfield"));
    assert_eq!((kinds.unit.as_str(), rates.unit.as_str()), ("", "m/s"));
    let mut types = [0usize; 4];
    for (&k, &r) in kinds.values.iter().zip(&rates.values) {
        types[k as usize] += 1;
        assert!(r.is_finite());
        if k.to_bits() == 0 {
            assert_eq!(r.to_bits(), 0);
        }
    }
    assert!(types.iter().all(|&t| t > 0), "{types:?}");
    // Crust: oceanic 7 km, continental 15 to 70 km (thinned towards the
    // margins, ticket M2-11); ages up to 180 Myr;
    // uplift up to 10 mm/yr, only on continental crust; trenches only on
    // oceanic convergent cells, and some trench.
    let thickness = read_kfield(&dir.join("crust_thickness.kfield"));
    let age = read_kfield(&dir.join("crust_age.kfield"));
    let uplift = read_kfield(&dir.join("uplift_rate.kfield"));
    let trench = read_kfield(&dir.join("trench.kfield"));
    assert_eq!(
        [&thickness, &age, &uplift, &trench].map(|f| f.unit.as_str()),
        ["m", "s", "m/s", ""]
    );
    let max_age_s = 180.0e6 * 31_557_600.0;
    let max_uplift_m_s = 0.01 / 31_557_600.0;
    let mut trenches = 0;
    for i in 0..1_000 {
        let continental = crust.values[i].to_bits() == land;
        let (t, a, u) = (thickness.values[i], age.values[i], uplift.values[i]);
        if continental {
            assert!((15e3..=70e3).contains(&t), "cell {i}: {t}");
        } else {
            assert_eq!(t.to_bits(), 7e3_f64.to_bits(), "cell {i}");
            assert_eq!(u.to_bits(), 0, "cell {i}");
        }
        assert!((0.0..=max_age_s).contains(&a), "cell {i}: {a}");
        assert!((0.0..=max_uplift_m_s).contains(&u), "cell {i}: {u}");
        if trench.values[i].to_bits() == 1.0_f64.to_bits() {
            trenches += 1;
            assert!(!continental && kinds.values[i].to_bits() == 1.0_f64.to_bits());
        } else {
            assert_eq!(trench.values[i].to_bits(), 0, "cell {i}");
        }
    }
    assert!(trenches > 0);
    // Surface heights (ticket M2-05): continental crust 15–70 km thick
    // stands at −2,730 m to +5,603 m (Airy isostasy), the ocean floor at
    // 2,500 m to 6,400 m deep (Parsons & Sclater), trenches at −8 km; the sea
    // level lies between the lowest and the highest cell, with both land
    // and sea.
    let height = read_kfield(&dir.join("surface_height.kfield"));
    assert_eq!(height.unit, "m");
    for i in 0..1_000 {
        let h = height.values[i];
        if crust.values[i].to_bits() == land {
            assert!((-2_730.31..=5_603.04).contains(&h), "cell {i}: {h}");
        } else if trench.values[i].to_bits() == 1.0_f64.to_bits() {
            assert_eq!(h.to_bits(), (-8_000.0_f64).to_bits(), "cell {i}");
        } else {
            assert!((-6_400.0..=-2_500.0).contains(&h), "cell {i}: {h}");
        }
    }
    let sea_level = sea_level_m(&json);
    assert!(height.values.iter().any(|&h| h > sea_level), "{sea_level}");
    assert!(height.values.iter().any(|&h| h < sea_level), "{sea_level}");
}

/// The sea level (m) in the text of `world.json`: the number on the line of
/// `sea_level_m`, written in the shortest form that parses back to the same
/// bits.
fn sea_level_m(json: &str) -> f64 {
    let line = json
        .lines()
        .find_map(|l| l.strip_prefix("  \"sea_level_m\": "))
        .unwrap_or_else(|| panic!("no sea level in\n{json}"));
    line.trim_end_matches(',')
        .parse()
        .unwrap_or_else(|e| panic!("{line}: {e}"))
}

#[test]
fn generate_writes_the_same_bytes_twice() {
    let (a, b) = (scratch("generate_twice_a"), scratch("generate_twice_b"));
    for dir in [&a, &b] {
        let out = dir.to_str().expect("UTF-8 path");
        let o = kuniumi(&["generate", "--seed", "3", "--cells", "1500", "--out", out]);
        assert_eq!(o.status.code(), Some(2), "{o:?}");
    }
    for dir in [&a, &b] {
        let mut files: Vec<_> = std::fs::read_dir(dir)
            .expect("the output directory")
            .map(|e| {
                e.expect("an entry")
                    .file_name()
                    .into_string()
                    .expect("UTF-8")
            })
            .collect();
        files.sort();
        let mut expected = world_files();
        expected.sort();
        assert_eq!(files, expected, "{}", dir.display());
    }
    for file in world_files() {
        let bytes = |dir: &Path| std::fs::read(dir.join(&file)).expect("written");
        assert!(bytes(&a) == bytes(&b), "{file} differs");
    }
}

/// Fingerprints of the fields for 10,000 cells, seed 1, as `kuniumi
/// generate` writes them. The Grid stage's two were recorded by ticket M2-01
/// after the invariant tests passed and the maps of `out/M2-01` were looked
/// at; both are exact copies of the grid, whose fingerprints M1-06 recorded
/// (`stages::grid::tests::golden_fingerprints`). The Tectonics stage's first
/// three were recorded by ticket M2-02 after its invariant tests passed
/// (`stages::tectonics::tests`) and the maps of `out/M2-02` were looked at;
/// its last three (continents and boundaries) by ticket M2-03 after the
/// tests of `stages::tectonics::continents` and `::boundaries` passed and
/// the maps of `out/M2-03` were looked at; the uplift by ticket M2-04 after
/// the tests of `stages::tectonics::crust` passed and the maps of
/// `out/M2-04` were looked at; the thickness, age and trenches again by
/// ticket M2-11 (margins, Earth's area–age distribution and the
/// time-nearest ridge, world-breaking) after the tests of
/// `stages::tectonics::crust` (`tests.rs`, `age_tests.rs`,
/// `margin_tests.rs`) passed and the maps of `out/M2-11` were looked at;
/// the surface height by ticket M2-05 after the tests of `stages::crust`
/// passed and the maps of `out/M2-05` were looked at.
const GOLDEN_FIELDS: [(&str, u64); 13] = [
    ("cell_area", 0x3293_5063_C4BF_C82D),
    ("neighbor_count", 0xBCFA_B43E_DA1C_53F6),
    ("plate_id", 0x5E9F_DB0A_B6B7_D521),
    ("plate_velocity_east", 0x8FFF_DF16_93AF_9B25),
    ("plate_velocity_north", 0x1A26_5DF6_79D2_2BB9),
    ("crust_type", 0x791F_6B17_52CF_D2B8),
    ("boundary_type", 0x5163_2C85_7C3F_C400),
    ("convergence_rate", 0x04BB_B545_5B56_84CC),
    ("crust_thickness", 0x9118_E6AA_2F43_6389),
    ("crust_age", 0x4293_1D66_883A_3DF9),
    ("uplift_rate", 0x6725_720D_C775_FE38),
    ("trench", 0xA2D4_773C_1D13_3C93),
    ("surface_height", 0xE587_6747_8054_1DE7),
];

/// The sea level (m) of the same world, as `world.json` writes it (the
/// shortest decimal that parses back to the same bits), recorded by ticket
/// M2-05 with `surface_height` (see [`GOLDEN_FIELDS`]).
const GOLDEN_SEA_LEVEL_M: f64 = -781.763_767_817_906_6;

#[test]
fn generate_golden_fingerprints() {
    let dir = scratch("generate_golden_fingerprints");
    let out = dir.to_str().expect("UTF-8 path");
    let o = kuniumi(&["generate", "--seed", "1", "--cells", "10000", "--out", out]);
    assert_eq!(o.status.code(), Some(2), "{o:?}");
    for (name, golden) in GOLDEN_FIELDS {
        let field = read_kfield(&dir.join(format!("{name}.kfield")));
        let fingerprint = kuniumi_io::field_file::fingerprint(&field);
        assert_eq!(fingerprint, golden, "{name}: {fingerprint:#018X}");
    }
    let json = std::fs::read_to_string(dir.join("world.json")).expect("world.json");
    let sea_level = sea_level_m(&json);
    assert_eq!(
        sea_level.to_bits(),
        GOLDEN_SEA_LEVEL_M.to_bits(),
        "sea level {sea_level}"
    );
}

/// The value of `key` in `metrics.txt`.
fn metric(text: &str, key: &str) -> String {
    let prefix = format!("{key} = ");
    text.lines()
        .find_map(|l| l.strip_prefix(prefix.as_str()))
        .unwrap_or_else(|| panic!("no `{key}` in\n{text}"))
        .to_owned()
}

/// The numbers of `metrics.txt` of the golden world (10,000 cells, seed 1),
/// from `tools/hypsometry_reference.py out/M2-06/10k` (ticket M2-06): exact
/// rational sums and bins over the `.kfield` files of the world, which
/// [`GOLDEN_FIELDS`] fixes. Shares are printed with 6 decimals, metres with
/// 1, so each agrees within half the last digit (and a rounding error).
const REFERENCE_METRICS: [(&str, f64, f64); 21] = [
    ("land_fraction", 0.296_053_819_165_441, 5e-7),
    ("mean_land_elevation_m", 1_061.460_292_780_781, 0.05),
    ("mean_ocean_depth_m", 3_718.062_915_493_437_2, 0.05),
    ("highest_m", 6_384.794_070_848_209, 0.05),
    ("lowest_m", -7_218.236_232_182_094, 0.05),
    ("percentile_1_m", -7_218.236_232_182_094, 0.05),
    ("percentile_5_m", -5_322.235_678_982_984, 0.05),
    ("percentile_10_m", -5_122.144_186_468_786, 0.05),
    ("percentile_25_m", -4_524.030_165_973_157, 0.05),
    ("percentile_50_m", -3_133.760_124_026_571, 0.05),
    ("percentile_75_m", 1_077.949_929_428_346_9, 0.05),
    ("percentile_90_m", 1_081.763_767_817_906_7, 0.05),
    ("percentile_95_m", 1_081.763_767_817_906_7, 0.05),
    ("percentile_99_m", 1_081.763_767_817_906_7, 0.05),
    ("ocean_mode_m", -1_750.0, 0.0),
    ("ocean_mode_fraction", 0.028_606_005_230_792_698, 5e-7),
    ("land_mode_m", 1_050.0, 0.0),
    ("land_mode_fraction", 0.245_278_992_164_091_3, 5e-7),
    ("valley_m", -850.0, 0.0),
    ("valley_fraction", 0.000_314_532_086_759_371_95, 5e-7),
    ("sea_level_m", -781.763_767_817_906_6, 5e-4),
];

#[test]
fn generate_writes_the_metrics_of_the_reference() {
    let dir = scratch("generate_writes_the_metrics_of_the_reference");
    let out = dir.to_str().expect("UTF-8 path");
    let o = kuniumi(&["generate", "--seed", "1", "--cells", "10000", "--out", out]);
    assert_eq!(o.status.code(), Some(2), "{o:?}");
    let text = std::fs::read_to_string(dir.join("metrics.txt")).expect("metrics.txt");
    for (key, reference, tolerance) in REFERENCE_METRICS {
        let value: f64 = metric(&text, key).parse().expect("a number");
        assert!(
            (value - reference).abs() <= tolerance * (1.0 + 1e-9),
            "{key}: {value} against {reference}"
        );
    }
    assert_eq!(metric(&text, "separate_peaks"), "yes");
    assert_eq!(metric(&text, "cells"), "10000");
    // 137 bins of 100 m, from −7,300 m to +6,400 m.
    let bins: Vec<&str> = text
        .lines()
        .skip_while(|l| !l.starts_with("bin_width_m = "))
        .skip(1)
        .collect();
    assert_eq!(bins.len(), 137);
    assert!(bins[0].starts_with("-7300 -7200 "), "{}", bins[0]);
    assert!(bins[136].starts_with("6300 6400 "), "{}", bins[136]);
    assert_eq!(png_size(&dir.join("hypsometry.png")), (1040, 540));
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

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

/// The fields of a world stopped after the Erosion stage, in the order
/// `generate` writes them: the Climate stage's (ticket M4-07) after the
/// Crust stage's `surface_height`.
const WORLD_FIELDS: [&str; 46] = [
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
    "temperature_annual",
    "temperature_m01",
    "temperature_m02",
    "temperature_m03",
    "temperature_m04",
    "temperature_m05",
    "temperature_m06",
    "temperature_m07",
    "temperature_m08",
    "temperature_m09",
    "temperature_m10",
    "temperature_m11",
    "temperature_m12",
    "precipitation_annual",
    "precipitation_m01",
    "precipitation_m02",
    "precipitation_m03",
    "precipitation_m04",
    "precipitation_m05",
    "precipitation_m06",
    "precipitation_m07",
    "precipitation_m08",
    "precipitation_m09",
    "precipitation_m10",
    "precipitation_m11",
    "precipitation_m12",
    "wind_east",
    "wind_north",
    "current_east",
    "current_north",
    "drainage_area",
    "discharge",
    "depression_depth",
];

/// The files of a world stopped after the Erosion stage, in the order
/// `generate` writes them: the grid (ticket MV-01), `<name>.kfield` and
/// `<name>.png` of each field, the metrics (ticket M2-06), the river map
/// (ticket M3-05), then `world.json`.
fn world_files() -> Vec<String> {
    let mut files = vec!["grid.kgrid".to_owned()];
    files.extend(
        WORLD_FIELDS
            .iter()
            .flat_map(|f| [format!("{f}.kfield"), format!("{f}.png")]),
    );
    files.extend([
        "metrics.txt".into(),
        "hypsometry.png".into(),
        "rivers.png".into(),
    ]);
    files.push("world.json".into());
    files
}

#[test]
fn generate_writes_the_world_so_far_then_exits_with_code_2() {
    let dir = scratch("generate_writes_the_world_so_far");
    let out = dir.to_str().expect("UTF-8 path");
    let o = kuniumi(&["generate", "--seed", "1", "--cells", "1000", "--out", out]);
    assert_eq!(o.status.code(), Some(2), "{o:?}");
    // Each stage's progress in whole percent, each line once: the Climate
    // stage's after each of its models (tickets M4-07, M4-17), and the 1,000 steps
    // of the Erosion stage give every percent from 0 to 100.
    let mut expected = String::from(
        "[grid] 0%\n[grid] 100%\n[tectonics] 0%\n[tectonics] 100%\n\
         [crust] 0%\n[crust] 100%\n\
         [climate] 0%\n[climate] 16%\n[climate] 94%\n[climate] 100%\n",
    );
    for percent in 0..=100 {
        expected.push_str(&format!("[erosion] {percent}%\n"));
    }
    expected.push_str("error: stage `hydrology` is not implemented yet\n");
    assert_eq!(stderr(&o), expected);
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
    assert_eq!(png_size(&dir.join("rivers.png")), (2048, 1024));
    let metrics = std::fs::read_to_string(dir.join("metrics.txt")).expect("metrics.txt");
    for line in ["seed = 1", "cells = 1000", "bin_width_m = 100"] {
        assert!(
            metrics.lines().any(|l| l == line),
            "`{line}` missing:\n{metrics}"
        );
    }
    let json = std::fs::read_to_string(dir.join("world.json")).expect("world.json");
    // The default erosion settings and Earth's longitude of perihelion as
    // Python's repr prints them: repr(5e-3 / 31557600.0),
    // repr(1e4 * 31557600.0), repr(102.94 * (math.pi / 180)).
    for line in [
        "  \"format_version\": 5,",
        "  \"seed\": 1,",
        "    \"cells\": 1000,",
        "    \"plate_count\": 12,",
        "    \"allow_partial\": true,",
        "    \"erosion\": {",
        "      \"erodibility\": 7.5e-10,",
        "      \"discharge_exponent\": 0.45,",
        "      \"diffusivity_m2_s\": 1.5844043907014476e-10,",
        "      \"time_step_s\": 315576000000.0,",
        "      \"steps\": 1000",
        "    \"radius_m\": 6371000.0,",
        "    \"perihelion_longitude_rad\": 1.7966419320029627,",
        "    \"grid\",",
        "    \"tectonics\",",
        "    \"crust\",",
        "    \"climate\",",
        "    \"erosion\"",
        // Earth's climate settles: no warning (ticket M4-16).
        "  \"warnings\": [],",
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
        "      \"name\": \"temperature_annual\",",
        "      \"unit\": \"K\",",
        "      \"name\": \"precipitation_m12\",",
        "      \"name\": \"current_north\",",
        "      \"name\": \"drainage_area\",",
        "      \"unit\": \"m^2\",",
        "      \"name\": \"discharge\",",
        "      \"unit\": \"m^3/s\",",
        "      \"name\": \"depression_depth\",",
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
    // and sea. The Erosion stage (ticket M3-05) keeps the ocean as it is and
    // the land above the sea level.
    let height = read_kfield(&dir.join("surface_height.kfield"));
    assert_eq!(height.unit, "m");
    let sea_level = sea_level_m(&json);
    for i in 0..1_000 {
        let h = height.values[i];
        if h > sea_level {
            assert!(crust.values[i].to_bits() == land, "cell {i}: {h}");
        } else if crust.values[i].to_bits() == land {
            assert!((-2_730.31..=sea_level).contains(&h), "cell {i}: {h}");
        } else if trench.values[i].to_bits() == 1.0_f64.to_bits() {
            assert_eq!(h.to_bits(), (-8_000.0_f64).to_bits(), "cell {i}");
        } else {
            assert!((-6_400.0..=-2_500.0).contains(&h), "cell {i}: {h}");
        }
    }
    assert!(height.values.iter().any(|&h| h > sea_level), "{sea_level}");
    assert!(height.values.iter().any(|&h| h < sea_level), "{sea_level}");
    // The rivers (ticket M3-05): every cell drains at least its own area,
    // the discharge is 1 m/yr of rain on the drainage area, the depression
    // depths are not negative.
    let area = read_kfield(&dir.join("drainage_area.kfield"));
    let discharge = read_kfield(&dir.join("discharge.kfield"));
    let depth = read_kfield(&dir.join("depression_depth.kfield"));
    assert_eq!(
        [&area, &discharge, &depth].map(|f| f.unit.as_str()),
        ["m^2", "m^3/s", "m"]
    );
    let r2 = 6_371_000.0_f64 * 6_371_000.0;
    for i in 0..1_000 {
        assert!(
            area.values[i] >= areas.values[i] * r2 * (1.0 - 1e-15),
            "cell {i}"
        );
        let q = area.values[i] / 31_557_600.0;
        assert!((discharge.values[i] - q).abs() <= 1e-12 * q, "cell {i}");
        assert!(depth.values[i] >= 0.0, "cell {i}");
    }
    // The map of the surface height (ticket M2-14): the elevation map at
    // the sea level over the chart's x axis, −11 km to +9 km above the sea
    // level for this world, the bounds written here.
    assert_eq!(metric(&metrics, "chart_x_min_m"), "-11000");
    assert_eq!(metric(&metrics, "chart_x_max_m"), "9000");
    let grid = kuniumi_sphere::SphereGrid::fibonacci(1_000, 1).expect("the grid");
    let pixels = kuniumi_io::render::equirect(
        &grid,
        &height.values,
        kuniumi_io::colormap::Colormap::Elevation { sea_level },
        kuniumi_io::colormap::Range::Fixed {
            min: sea_level - 11_000.0,
            max: sea_level + 9_000.0,
        },
        2048,
    );
    let map = std::fs::read(dir.join("surface_height.png")).expect("the map");
    assert!(
        map == kuniumi_io::png::encode_rgb8(2048, 1024, &pixels),
        "surface_height.png"
    );
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

#[test]
fn generate_writes_a_world_that_reads_back() {
    // `world_dir::read` gives the world of `kuniumi_sim::generate` with the
    // same seed and settings (those of `generate`: Earth, the default
    // settings with `--cells` and `allow_partial`), bit for bit.
    let dir = scratch("generate_reads_back");
    let out = dir.to_str().expect("UTF-8 path");
    let o = kuniumi(&["generate", "--seed", "5", "--cells", "1200", "--out", out]);
    assert_eq!(o.status.code(), Some(2), "{o:?}");
    let back = kuniumi_io::world_dir::read(&dir).unwrap_or_else(|e| panic!("{e}"));
    let settings = kuniumi_world::Settings {
        cells: 1_200,
        allow_partial: true,
        ..kuniumi_world::Settings::default()
    };
    let world = kuniumi_sim::generate(
        &kuniumi_world::PlanetParams::earth(),
        &settings,
        5,
        &mut |_| {},
    )
    .expect("the partial world");
    assert_eq!(back, world);
    assert_eq!(back.sea_level_m.to_bits(), world.sea_level_m.to_bits());
    let names: Vec<_> = back.fields.iter().map(|f| f.name.as_str()).collect();
    assert_eq!(names, WORLD_FIELDS);
    for (a, b) in back.fields.iter().zip(&world.fields) {
        assert_eq!(a.kind, b.kind, "{}", a.name);
        let bits =
            |f: &kuniumi_world::Field| f.values.iter().map(|v| v.to_bits()).collect::<Vec<_>>();
        assert!(bits(a) == bits(b), "{}", a.name);
    }
    let grid_bits = |w: &kuniumi_world::World| {
        (0..w.grid.len() as u32)
            .flat_map(|i| {
                let p = w.grid.position(i);
                [p.x, p.y, p.z, w.grid.area_sr(i)]
                    .into_iter()
                    .chain(w.grid.edge_lengths_rad(i).iter().copied())
                    .chain(w.grid.centre_distances_rad(i).iter().copied())
                    .map(f64::to_bits)
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>()
    };
    assert!(grid_bits(&back) == grid_bits(&world));
}

/// Fingerprints of the fields for 10,000 cells, seed 1, as `kuniumi
/// generate` writes them with 10 erosion steps ([`erosion_steps_file`]).
/// The Grid stage's two were recorded by ticket M2-01
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
/// passed and the maps of `out/M2-05` were looked at; the crust type,
/// thickness, age, uplift and surface height again by ticket M2-10 (the
/// oceanic pockets filled and the coast trimmed, world-breaking) after the
/// tests of `stages::tectonics::continents` passed and the maps of
/// `out/M2-10` were looked at; the age and surface height again by ticket
/// M2-12 (the floor of plates without a ridge ranked as the oldest,
/// world-breaking) after the tests of `stages::tectonics::crust`
/// (`age_tests.rs`) passed and the maps of `out/M2-12` were looked at; and
/// again by ticket M2-13 (a floor of 1 cm/yr on the spreading rate of the
/// ridges' fronts, world-breaking) after the tests of
/// `stages::tectonics::crust` (`age_tests.rs`, `tests.rs`) passed and the
/// maps of `out/M2-13` were looked at. That surface height, of the Crust
/// stage, is [`GOLDEN_CRUST_SURFACE_HEIGHT`] since ticket M3-05, which
/// erodes it: the surface height after 10 steps of the Erosion stage and
/// its three rivers were recorded by ticket M3-05 after the tests of
/// `stages::erosion` (`stage_tests.rs`) passed and the maps of `out/M3-05`
/// were looked at (world-breaking), and again by ticket M3-06 (the land
/// roughened by up to 1 m before the first erosion step, world-breaking)
/// after the tests of `stages::erosion` (`roughness_tests.rs`,
/// `stage_tests.rs`) passed and the maps of `out/M3-06` were looked at. The
/// 30 fields of the Climate stage were recorded by ticket M4-07 after the
/// tests of `stages::climate` (`stage_tests.rs`) passed and the maps of
/// `out/M4-07` were looked at; the 16 others are unchanged by it. The 30
/// were recorded again by ticket M4-15 (the balance step after each year
/// of the energy balance and the sea surface, world-breaking) after the
/// tests of `stages::climate` passed and the maps of `out/M4-15` were
/// looked at; the 16 others are unchanged by it. And again by ticket M4-17
/// (the coupled state of the sea surface, winds and currents,
/// world-breaking) after the tests of `stages::climate` passed and the maps
/// of `out/M4-17` were looked at; the 16 others are unchanged by it.
const GOLDEN_FIELDS: [(&str, u64); 46] = [
    ("cell_area", 0x3293_5063_C4BF_C82D),
    ("neighbor_count", 0xBCFA_B43E_DA1C_53F6),
    ("plate_id", 0x5E9F_DB0A_B6B7_D521),
    ("plate_velocity_east", 0x8FFF_DF16_93AF_9B25),
    ("plate_velocity_north", 0x1A26_5DF6_79D2_2BB9),
    ("crust_type", 0x4C68_F0C4_A729_DBF8),
    ("boundary_type", 0x5163_2C85_7C3F_C400),
    ("convergence_rate", 0x04BB_B545_5B56_84CC),
    ("crust_thickness", 0x865C_5974_07DF_8E7B),
    ("crust_age", 0x8B16_21E7_DFBF_D241),
    ("uplift_rate", 0x05E7_9013_D625_1ECB),
    ("trench", 0xA2D4_773C_1D13_3C93),
    ("surface_height", 0xDB1E_0D68_718E_F38C),
    ("temperature_annual", 0x7CE0_1DDE_DBB3_729C),
    ("temperature_m01", 0xC5BB_47B3_AEB1_ECAB),
    ("temperature_m02", 0x9EE8_A8C5_F7A2_6690),
    ("temperature_m03", 0x7BBE_9D3C_005C_20EF),
    ("temperature_m04", 0xA0A4_9B15_0B02_C3E6),
    ("temperature_m05", 0x9B20_7725_95BA_D377),
    ("temperature_m06", 0xAE99_4EF5_8D04_18D7),
    ("temperature_m07", 0x0F5F_F224_6C0C_915C),
    ("temperature_m08", 0xD228_C598_6542_8BE0),
    ("temperature_m09", 0xDDBB_7B54_A0A7_49FA),
    ("temperature_m10", 0x717B_F52A_183F_7735),
    ("temperature_m11", 0x31C4_F28C_89C0_1129),
    ("temperature_m12", 0x4611_28CA_64DD_6DD0),
    ("precipitation_annual", 0x8E67_3DD9_3481_E067),
    ("precipitation_m01", 0xCD31_8F87_72EE_320C),
    ("precipitation_m02", 0x53FB_E448_4E66_B695),
    ("precipitation_m03", 0x4B35_BA7F_BAC2_A065),
    ("precipitation_m04", 0xB2E9_7C93_E507_BDA4),
    ("precipitation_m05", 0x2669_1406_DF47_8B2E),
    ("precipitation_m06", 0xC99E_C844_7BEA_7EB2),
    ("precipitation_m07", 0x23A8_8710_177E_913D),
    ("precipitation_m08", 0x4E6B_12C0_5441_B263),
    ("precipitation_m09", 0x0F23_D098_1233_0F74),
    ("precipitation_m10", 0x7E04_FA03_A5D7_0ABC),
    ("precipitation_m11", 0x026C_B47F_0DA8_2B61),
    ("precipitation_m12", 0x5448_8ABF_41A4_6FCF),
    ("wind_east", 0x370B_EDC8_7190_13D4),
    ("wind_north", 0x72C1_226A_0F3A_7F6C),
    ("current_east", 0xC800_0976_C85F_A52C),
    ("current_north", 0xEA32_99F2_B702_4EF2),
    ("drainage_area", 0x74B7_97DE_863F_9C2A),
    ("discharge", 0x714B_06C0_204A_48F6),
    ("depression_depth", 0x2128_8731_7CB9_00EB),
];

/// The fingerprint of `surface_height` of the same world without erosion
/// steps, as the Crust stage leaves it: recorded by ticket M2-05 and again by
/// tickets M2-10, M2-12 and M2-13 (see [`GOLDEN_FIELDS`]).
const GOLDEN_CRUST_SURFACE_HEIGHT: u64 = 0xA90E_8D61_5D3A_E286;

/// A settings file of Earth and the default settings with `steps` erosion
/// steps, in a scratch directory of its own.
fn erosion_steps_file(name: &str, steps: u32) -> PathBuf {
    settings_file(name, &format!("[generation.erosion]\nsteps = {steps}\n"))
}

/// The sea level (m) of the same world, as `world.json` writes it (the
/// shortest decimal that parses back to the same bits), recorded by ticket
/// M2-05 with `surface_height` and again by tickets M2-10, M2-12 and M2-13
/// (see [`GOLDEN_FIELDS`]).
const GOLDEN_SEA_LEVEL_M: f64 = -747.218_732_917_737_4;

/// Runs `kuniumi generate` for 10,000 cells, seed 1, with `steps` erosion
/// steps, into a scratch directory named `name`, and returns the directory.
fn golden_world(name: &str, steps: u32) -> PathBuf {
    let planet = erosion_steps_file(&format!("{name}_settings"), steps);
    let dir = scratch(name);
    let out = dir.to_str().expect("UTF-8 path");
    let planet = planet.to_str().expect("UTF-8 path");
    let args = [
        "generate", "--seed", "1", "--cells", "10000", "--planet", planet, "--out", out,
    ];
    let o = kuniumi(&args);
    assert_eq!(o.status.code(), Some(2), "{o:?}");
    dir
}

/// The fingerprint of the field `name` in the world directory `dir`.
fn fingerprint_of(dir: &Path, name: &str) -> u64 {
    kuniumi_io::field_file::fingerprint(&read_kfield(&dir.join(format!("{name}.kfield"))))
}

/// Checks the sea level of the world directory `dir` against
/// [`GOLDEN_SEA_LEVEL_M`].
fn assert_golden_sea_level(dir: &Path) {
    let json = std::fs::read_to_string(dir.join("world.json")).expect("world.json");
    let sea_level = sea_level_m(&json);
    assert_eq!(
        sea_level.to_bits(),
        GOLDEN_SEA_LEVEL_M.to_bits(),
        "sea level {sea_level}"
    );
}

#[test]
fn generate_golden_fingerprints() {
    let dir = golden_world("generate_golden_fingerprints", 10);
    let fingerprints: Vec<(&str, u64)> = GOLDEN_FIELDS
        .iter()
        .map(|&(name, _)| (name, fingerprint_of(&dir, name)))
        .collect();
    for (name, fingerprint) in &fingerprints {
        println!("{name}: {fingerprint:#018X}");
    }
    assert_eq!(fingerprints, GOLDEN_FIELDS);
    assert_golden_sea_level(&dir);
    // Its climate settles (ticket M4-16).
    let json = std::fs::read_to_string(dir.join("world.json")).expect("world.json");
    assert!(json.lines().any(|l| l == "  \"warnings\": [],"), "{json}");
}

#[test]
fn generate_without_erosion_steps_keeps_the_relief_of_the_crust_stage() {
    // With 0 steps the Erosion stage leaves the surface height of the Crust
    // stage bit for bit; the fields of the earlier stages, the climate's,
    // which is that of the Crust stage's relief, and the sea level never
    // depend on the erosion settings.
    let dir = golden_world("generate_without_erosion_steps", 0);
    for (name, golden) in GOLDEN_FIELDS[..12].iter().chain(&GOLDEN_FIELDS[13..43]) {
        assert_eq!(fingerprint_of(&dir, name), *golden, "{name}");
    }
    let height = fingerprint_of(&dir, "surface_height");
    assert_eq!(height, GOLDEN_CRUST_SURFACE_HEIGHT, "{height:#018X}");
    assert_golden_sea_level(&dir);
}

/// The value of `key` in `metrics.txt`.
fn metric(text: &str, key: &str) -> String {
    let prefix = format!("{key} = ");
    text.lines()
        .find_map(|l| l.strip_prefix(prefix.as_str()))
        .unwrap_or_else(|| panic!("no `{key}` in\n{text}"))
        .to_owned()
}

/// The numbers of `metrics.txt` of the golden world (10,000 cells, seed 1)
/// without erosion steps, whose surface is that of the Crust stage
/// ([`GOLDEN_CRUST_SURFACE_HEIGHT`]), from `tools/hypsometry_reference.py
/// out/M2-13/10k` (tickets M2-06, M2-10, M2-12 and M2-13, run again on the
/// world of the new golden values): exact
/// rational sums and bins over the `.kfield` files of the world, which
/// [`GOLDEN_FIELDS`] fixes. Shares are printed with 6 decimals, metres with
/// 1, so each agrees within half the last digit (and a rounding error).
const REFERENCE_METRICS: [(&str, f64, f64); 21] = [
    ("land_fraction", 0.291_091_576_212_747_34, 5e-7),
    ("mean_land_elevation_m", 1_032.724_038_002_781, 0.05),
    ("mean_ocean_depth_m", 3_692.037_083_550_971, 0.05),
    ("highest_m", 6_350.249_035_948_041, 0.05),
    ("lowest_m", -7_252.781_267_082_262, 0.05),
    ("percentile_1_m", -7_252.781_267_082_262, 0.05),
    ("percentile_5_m", -5_282.689_705_833_666, 0.05),
    ("percentile_10_m", -5_097.280_904_695_039, 0.05),
    ("percentile_25_m", -4_507.220_560_074_717, 0.05),
    ("percentile_50_m", -3_143.220_779_432_618, 0.05),
    ("percentile_75_m", 992.208_819_451_451_1, 0.05),
    ("percentile_90_m", 1_047.218_732_917_737_3, 0.05),
    ("percentile_95_m", 1_047.218_732_917_737_3, 0.05),
    ("percentile_99_m", 1_047.218_732_917_737_3, 0.05),
    ("ocean_mode_m", -1_750.0, 0.0),
    ("ocean_mode_fraction", 0.028_872_874_312_862_1, 5e-7),
    ("land_mode_m", 1_050.0, 0.0),
    ("land_mode_fraction", 0.241_139_973_356_619_46, 5e-7),
    ("valley_m", -950.0, 0.0),
    ("valley_fraction", 0.000_281_920_946_343_612_6, 5e-7),
    ("sea_level_m", -747.218_732_917_737_4, 5e-4),
];

#[test]
fn generate_writes_the_metrics_of_the_reference() {
    let dir = golden_world("generate_writes_the_metrics_of_the_reference", 0);
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
fn generate_writes_the_drainage_metrics() {
    // The golden world after 10 erosion steps (ticket M3-08): its largest
    // basin drains 66 cells, so no basin reaches the 100 cells of the fit;
    // the largest drainage area and the exceedance shares agree with those
    // computed here from the `.kfield` files and the sea level of
    // `world.json`, independently of the graph that `metrics.txt` rebuilds.
    let dir = golden_world("generate_writes_the_drainage_metrics", 10);
    let text = std::fs::read_to_string(dir.join("metrics.txt")).expect("metrics.txt");
    let json = std::fs::read_to_string(dir.join("world.json")).expect("world.json");
    let sea_level = sea_level_m(&json);
    let height = read_kfield(&dir.join("surface_height.kfield")).values;
    let area = read_kfield(&dir.join("drainage_area.kfield")).values;
    let cell = read_kfield(&dir.join("cell_area.kfield")).values;
    let land: Vec<usize> = (0..height.len())
        .filter(|&i| height[i] > sea_level)
        .collect();
    assert_eq!(metric(&text, "hack_min_cells"), "100");
    // 100 · 4π (6,371 km)² / 10,000 = 5.1006e12 m².
    assert_eq!(metric(&text, "hack_min_area_m2"), "5.1006e12");
    assert_eq!(metric(&text, "hack_basins"), "0");
    assert_eq!(metric(&text, "hack_exponent"), "none");
    assert_eq!(metric(&text, "hack_coefficient"), "none");
    let largest = land.iter().map(|&i| area[i]).fold(0.0_f64, f64::max);
    assert_eq!(
        metric(&text, "largest_drainage_area_m2"),
        format!("{largest:.4e}")
    );
    let land_sr: f64 = land.iter().map(|&i| cell[i]).sum();
    let lines: Vec<(f64, f64)> = text
        .lines()
        .filter_map(|l| l.strip_prefix("exceedance_"))
        .map(|l| {
            let (key, share) = l.split_once("_m2 = ").expect("an exceedance line");
            (
                key.parse().expect("an area"),
                share.parse().expect("a share"),
            )
        })
        .collect();
    assert!(lines.len() >= 5, "{text}");
    for &(a, share) in &lines {
        let above: f64 = land
            .iter()
            .filter(|&&i| area[i] >= a)
            .map(|&i| cell[i])
            .sum();
        assert!(
            (share - above / land_sr).abs() <= 5e-7 * (1.0 + 1e-9),
            "{a:e}: {share} against {}",
            above / land_sr
        );
    }
    // From the largest value not above the smallest drainage area on land
    // to the smallest above the largest: shares 1 and 0.
    let smallest = land.iter().map(|&i| area[i]).fold(f64::INFINITY, f64::min);
    assert!(lines[0].0 <= smallest && smallest < lines[1].0, "{lines:?}");
    let last = lines.len() - 1;
    assert!(
        lines[last - 1].0 <= largest && largest < lines[last].0,
        "{lines:?}"
    );
    assert_eq!((lines[0].1, lines[last].1), (1.0, 0.0));
    // The section comes before the histogram, which ends the file.
    let section = text.find("\nhack_min_cells = ").expect("the section");
    assert!(section < text.find("\nbin_width_m = ").expect("the histogram"));
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
    assert!(stdout(&o).starts_with("# kuniumi settings file"));
    assert!(stdout(&o).contains("\nradius_m = 6371000.0 "));

    let dir = scratch("old_commands_still_work");
    std::fs::create_dir_all(&dir).expect("scratch directory");
    let file = dir.join("test.png");
    let o = kuniumi(&["test-image", file.to_str().expect("UTF-8 path")]);
    assert!(o.status.success(), "{o:?}");
    assert_eq!(png_size(&file), (360, 180));
}

#[test]
fn params_prints_the_earth_settings_file() {
    let o = kuniumi(&["params"]);
    assert!(o.status.success(), "{o:?}");
    assert_eq!(stderr(&o), "");
    let (earth, default) = (
        kuniumi_world::PlanetParams::earth(),
        kuniumi_world::Settings::default(),
    );
    // The file of `kuniumi_config::write` (whose text its own tests compare
    // with a hand-written one), which loads back as Earth.
    assert_eq!(stdout(&o), kuniumi_config::write(&earth, &default));
    assert_eq!(kuniumi_config::load(&stdout(&o)), Ok((earth, default)));
}

/// Writes `text` to `name` in a fresh scratch directory and returns the
/// file's path.
fn settings_file(name: &str, text: &str) -> PathBuf {
    let dir = scratch(name);
    std::fs::create_dir_all(&dir).expect("scratch directory");
    let path = dir.join("planet.toml");
    std::fs::write(&path, text).expect("settings file");
    path
}

#[test]
fn generate_with_the_earth_file_writes_the_world_without_a_file() {
    let o = kuniumi(&["params"]);
    let planet = settings_file("generate_earth_file", &stdout(&o));
    let planet = planet.to_str().expect("UTF-8 path");
    let (a, b) = (
        scratch("generate_earth_file_a"),
        scratch("generate_earth_file_b"),
    );
    let base = ["generate", "--seed", "3", "--cells", "1500", "--out"];
    let with_file = [&base[..], &[a.to_str().expect("UTF-8"), "--planet", planet]].concat();
    let without = [&base[..], &[b.to_str().expect("UTF-8")]].concat();
    for args in [with_file, without] {
        let o = kuniumi(&args);
        assert_eq!(o.status.code(), Some(2), "{args:?}: {o:?}");
    }
    for file in world_files() {
        let bytes = |dir: &Path| std::fs::read(dir.join(&file)).expect("written");
        assert!(bytes(&a) == bytes(&b), "{file} differs");
    }
}

#[test]
fn generate_reads_the_planet_file_and_the_command_line_overrides_cells() {
    // A Mars-sized planet in human units, with five plates and erosion
    // settings of its own; the file asks for `allow_partial = false`, which
    // `generate` overrides.
    let planet = settings_file(
        "generate_planet_file",
        "[planet]\nradius_km = 3389.5\nobliquity_deg = 25.19\n\
         perihelion_longitude_deg = 336.04\n\
         [generation]\ncells = 1000\nplate_count = 5\nallow_partial = false\n\
         [generation.erosion]\ntime_step_yr = 20000\nsteps = 5\n",
    );
    let planet = planet.to_str().expect("UTF-8 path");
    for (cells, expected) in [(None, 1_000), (Some("1200"), 1_200)] {
        let dir = scratch("generate_planet_file_out");
        let out = dir.to_str().expect("UTF-8 path");
        let mut args = vec!["generate", "--seed", "1", "--planet", planet, "--out", out];
        if let Some(cells) = cells {
            args.extend(["--cells", cells]);
        }
        let o = kuniumi(&args);
        assert_eq!(o.status.code(), Some(2), "{args:?}: {o:?}");
        assert!(
            stderr(&o).ends_with("error: stage `hydrology` is not implemented yet\n"),
            "{}",
            stderr(&o)
        );
        let json = std::fs::read_to_string(dir.join("world.json")).expect("world.json");
        // 3,389.5 km = 3,389,500 m; 25.19° and 336.04° (Mars's longitude
        // of perihelion) from Python, repr(25.19 * (math.pi / 180)) and
        // repr(336.04 * (math.pi / 180)); 20,000 Julian years,
        // 631,152,000,000 s.
        for line in [
            format!("    \"cells\": {expected},"),
            "    \"plate_count\": 5,".into(),
            "    \"allow_partial\": true,".into(),
            "      \"erodibility\": 7.5e-10,".into(),
            "      \"time_step_s\": 631152000000.0,".into(),
            "      \"steps\": 5".into(),
            "    \"radius_m\": 3389500.0,".into(),
            "    \"obliquity_rad\": 0.4396484385773716,".into(),
            "    \"perihelion_longitude_rad\": 5.865004418401745,".into(),
            "    \"surface_gravity_m_s2\": 9.80665,".into(),
        ] {
            assert!(json.lines().any(|l| l == line), "`{line}` missing:\n{json}");
        }
        let areas = read_kfield(&dir.join("cell_area.kfield"));
        assert_eq!(areas.values.len(), expected);
        let plates = read_kfield(&dir.join("plate_id.kfield"));
        let mut sizes = [0usize; 5];
        for &p in &plates.values {
            sizes[p as usize] += 1;
        }
        assert!(sizes.iter().all(|&s| s > 0), "{sizes:?}");
    }
}

#[test]
fn generate_prints_the_warnings_of_the_stages() {
    // A dry planet whose moisture is not periodic after its 10 years on
    // 1,000 cells (seed 3): 13,000 km, 11 m/s², 240 h, 90 days,
    // 2,500 W/m², no water, 10 bar, as corner 85 of the bounds tests of
    // kuniumi-sim but for the gravity; its vapour changed by 54.7 kg/m² in
    // the last year, its other models settled (the scratch tool of ticket
    // M4-16). One erosion step: the climate does not depend on erosion.
    let planet = settings_file(
        "generate_warnings_planet",
        "[planet]\nradius_m = 13000000.0\nsurface_gravity_m_s2 = 11.0\n\
         rotation_period_s = 864000.0\norbital_period_s = 7776000.0\n\
         stellar_flux_w_m2 = 2500.0\nocean_volume_m3 = 0.0\n\
         surface_pressure_pa = 1000000.0\n\
         [generation.erosion]\nsteps = 1\n",
    );
    let dir = scratch("generate_warnings");
    let out = dir.to_str().expect("UTF-8 path");
    let planet = planet.to_str().expect("UTF-8 path");
    let o = kuniumi(&[
        "generate", "--seed", "3", "--cells", "1000", "--planet", planet, "--out", out,
    ]);
    assert_eq!(o.status.code(), Some(2), "{o:?}");
    // The progress lines, then the warning, then the missing stage.
    let text = stderr(&o);
    let lines: Vec<&str> = text.lines().collect();
    let progress = [
        "[grid] 0%",
        "[grid] 100%",
        "[tectonics] 0%",
        "[tectonics] 100%",
        "[crust] 0%",
        "[crust] 100%",
        "[climate] 0%",
        "[climate] 16%",
        "[climate] 94%",
        "[climate] 100%",
        "[erosion] 0%",
        "[erosion] 100%",
    ];
    assert_eq!(lines.len(), progress.len() + 2, "{text}");
    assert_eq!(lines[..progress.len()], progress, "{text}");
    let warning = lines[progress.len()];
    assert!(
        warning.starts_with(
            "warning: climate: the moisture did not reach its periodic state in 10 years: "
        ),
        "{text}"
    );
    assert_eq!(
        lines[progress.len() + 1],
        "error: stage `hydrology` is not implemented yet"
    );
    // The world directory keeps the same warning.
    let world = kuniumi_io::world_dir::read(&dir).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(world.warnings.len(), 1, "{:?}", world.warnings);
    let kept = &world.warnings[0];
    assert_eq!(
        format!("warning: {}: {}", kept.stage, kept.message),
        warning
    );

    // The warnings come before the world is written (the owner's choice in
    // ticket M4-16), so they show even if writing fails: here the output
    // directory would lie under a regular file.
    let blocked = scratch("generate_warnings_blocked");
    std::fs::create_dir_all(&blocked).expect("a scratch directory");
    let file = blocked.join("file");
    std::fs::write(&file, "not a directory").expect("a regular file");
    let out = file.join("world");
    let o = kuniumi(&[
        "generate",
        "--seed",
        "3",
        "--cells",
        "1000",
        "--planet",
        planet,
        "--out",
        out.to_str().expect("UTF-8 path"),
    ]);
    assert_eq!(o.status.code(), Some(2), "{o:?}");
    let text = stderr(&o);
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(lines.len(), progress.len() + 2, "{text}");
    assert_eq!(lines[..progress.len()], progress, "{text}");
    assert_eq!(lines[progress.len()], warning, "{text}");
    let error = lines[progress.len() + 1];
    assert!(
        error.starts_with(&format!("error: {}", out.display())),
        "{text}"
    );
    assert_eq!(stdout(&o), "", "nothing written");
}

#[test]
fn generate_reports_a_bad_planet_file_and_writes_nothing() {
    let missing = scratch("generate_missing_planet_file").join("missing.toml");
    let range = settings_file("generate_bad_planet_range", "[planet]\nradius_km = -1\n");
    let unknown = settings_file("generate_bad_planet_key", "[planet]\nradius = 1.0\n");
    let both = settings_file(
        "generate_bad_planet_both",
        "[planet]\nradius_m = 1.0\nradius_km = 1.0\n",
    );
    // Earth's gravity on a planet of 13,000 km: 2,698 kg/m³, below the
    // densities of rocky planets.
    let inconsistent = settings_file(
        "generate_bad_planet_density",
        "[planet]\nradius_km = 13000\n",
    );
    let erosion = settings_file(
        "generate_bad_erosion_range",
        "[generation.erosion]\ntime_step_yr = -1\n",
    );
    // 1e300 m²/s times the default time step overflows.
    let step = settings_file(
        "generate_bad_erosion_step",
        "[generation.erosion]\ndiffusivity_m2_s = 1e300\n",
    );
    // What follows `error: <path>: `: its start, and a part of it.
    for (path, start, part) in [
        // The message of the operating system (error 2 on Unix and Windows).
        (&missing, "", "(os error 2)"),
        (
            &range,
            "`planet.radius_km` = -1 is out of range (radius_m = -1000)\n",
            "",
        ),
        (
            &unknown,
            "TOML parse error at line 2, column 1\n",
            "unknown field `radius`",
        ),
        (
            &both,
            "`planet.radius_m` and `planet.radius_km` are the same quantity; give one of them\n",
            "",
        ),
        (
            &inconsistent,
            "`planet.radius_km` and `planet.surface_gravity_m_s2` do not fit together: \
             bulk_density_kg_m3 = 2698.",
            " is out of range\n",
        ),
        (
            &erosion,
            "`generation.erosion.time_step_yr` = -1 is out of range \
             (erosion.time_step_s = -31557600)\n",
            "",
        ),
        (
            &step,
            "`generation.erosion.diffusivity_m2_s` and `generation.erosion.time_step_s` do not \
             fit together: erosion.diffusion_step_m2 = inf is out of range\n",
            "",
        ),
    ] {
        let dir = scratch("generate_bad_planet_file_out");
        let path = path.to_str().expect("UTF-8 path");
        let out = dir.to_str().expect("UTF-8 path");
        let o = kuniumi(&["generate", "--seed", "1", "--planet", path, "--out", out]);
        assert_eq!(o.status.code(), Some(2), "{path}: {o:?}");
        assert_eq!(stdout(&o), "", "{path}");
        let text = stderr(&o);
        let message = text
            .strip_prefix(&format!("error: {path}: "))
            .unwrap_or_else(|| panic!("{path}: {text}"));
        assert!(message.starts_with(start), "{path}: {text}");
        assert!(message.contains(part), "{path}: {text}");
        assert!(!dir.exists(), "{path}: {} was created", dir.display());
    }
}

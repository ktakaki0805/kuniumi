//! Tests of world directories:
//!
//! - `world_json_parses_back_bit_exactly` (Earth and a world of awkward
//!   numbers);
//! - `serde_json_round_trips_every_finite_f64` (random bit patterns, random
//!   decimals and hard cases);
//! - `earth_world_json_has_the_documented_layout`;
//! - `awkward_world_json_has_the_shortest_forms`;
//! - `the_same_world_gives_the_same_bytes`;
//! - `rejects_non_finite_metadata`;
//! - `rejects_bad_fields_before_writing_anything`;
//! - `writes_every_field_and_world_json_last`;
//! - `file_system_errors_name_the_path`;
//! - `maps_use_the_categorical_palette_for_classes` (and viridis for
//!   quantities, also beside a surface height);
//! - `surface_height_takes_the_elevation_map`,
//!   `land_and_sea_split_at_the_sea_level`,
//!   `an_elevation_has_one_colour_in_every_earth_like_world` and
//!   `the_range_falls_back_to_the_field_s_own` (ticket M2-14);
//! - `writes_maps_of_2048_by_1024`;
//! - `metrics_txt_by_hand` and `metrics_txt_without_land_or_peaks`;
//! - `chart_axes` and `hypsometry_chart_draws_the_bins_and_the_markers`;
//! - `writes_metrics_and_the_chart_before_world_json`;
//! - `no_metrics_without_a_surface_height_and_none_left_over`;
//! - `rejects_reserved_field_names_before_writing_anything`,
//!   `reserved_names_by_hand`, `fixed_files_are_every_file_but_the_fields`
//!   and `fields_may_take_the_stems_of_the_other_fixed_files` (ticket
//!   M2-19);
//! - `rejects_non_finite_heights_before_writing_anything`;
//! - `rejects_elevations_that_overflow_before_writing_anything` and
//!   `huge_elevations_have_finite_metrics` (ticket M2-15);
//! - `rejects_elevations_beyond_the_histogram_before_writing_anything` and
//!   `accepts_elevations_at_the_bound` (ticket M2-18);
//! - `river_colours_by_hand`, `river_colours_without_rivers`,
//!   `river_colours_need_one_discharge_per_height`,
//!   `writes_the_river_map_after_the_metrics` and
//!   `rejects_the_name_rivers_before_writing_anything` (ticket M3-05);
//! - `metrics_txt_with_drainage_by_hand`, `drainage_lines_by_hand`,
//!   `writes_drainage_metrics_after_erosion`,
//!   `rejects_a_radius_beyond_the_drainage_limits_before_writing_anything`
//!   and `rejects_heights_beyond_the_drainage_limits_before_writing_anything`
//!   (ticket M3-08).
//!
//! Reading directories back (ticket MV-01) is tested in `read_tests.rs`,
//! with the worlds of this file.
#![allow(clippy::disallowed_methods)]

use std::path::PathBuf;

use kuniumi_math::rng::Rng;
use kuniumi_sphere::SphereGrid;
use kuniumi_world::{
    ErosionSettings, Field, FieldKind, PlanetParams, RotationRegime, Settings, TectonicRegime,
    World,
};

use kuniumi_sim::metrics::{DrainageMetrics, ElevationMetrics, Exceedance, HackFit, earth};

use super::*;
use crate::chart;
use crate::colormap::{Colormap, Range};
use crate::render;

/// A world on a 1,000-cell grid with a quantity field and a class field,
/// and erosion settings that are valid but not the defaults.
pub(super) fn world() -> World {
    let grid = SphereGrid::fibonacci(1_000, 1).expect("a valid grid");
    let n = grid.len();
    let mut world = World {
        params: PlanetParams::earth(),
        settings: Settings {
            cells: 1_000,
            plate_count: 12,
            allow_partial: true,
            erosion: ErosionSettings {
                erodibility: 2.5e-10,
                discharge_exponent: 0.5,
                diffusivity_m2_s: 3.2e-10,
                time_step_s: 6.311_52e11,
                steps: 250,
            },
        },
        seed: 1,
        grid,
        sea_level_m: 0.0,
        stages: vec!["grid".into()],
        fields: Vec::new(),
    };
    let z: Vec<f64> = world.grid.positions().iter().map(|p| p.z).collect();
    world.add_field("height", "m", z).unwrap();
    let classes = (0..n).map(|i| (i % 5) as f64).collect();
    world.add_class_field("zone", classes).unwrap();
    world
}

/// A world whose metadata holds awkward numbers: subnormals, the extremes,
/// −0, values that need 17 digits, the v2 regimes, settings that are not
/// valid.
pub(super) fn awkward_world() -> World {
    World {
        params: PlanetParams {
            radius_m: f64::MIN_POSITIVE,
            surface_gravity_m_s2: 5e-324,
            rotation_period_s: f64::MAX,
            obliquity_rad: 0.1,
            orbital_period_s: 1e23,
            orbital_eccentricity: 1.0 / 3.0,
            perihelion_longitude_rad: -1e300,
            stellar_flux_w_m2: 2.225_073_858_507_201e-308,
            ocean_volume_m3: 9_007_199_254_740_991.0,
            surface_pressure_pa: -0.0,
            continental_crust_fraction: 0.1 + 0.2,
            rotation: RotationRegime::TidallyLocked,
            tectonics: TectonicRegime::StagnantLid,
        },
        settings: Settings {
            cells: u32::MAX,
            plate_count: u32::MAX - 1,
            allow_partial: false,
            erosion: ErosionSettings {
                erodibility: 5e-324,
                discharge_exponent: -0.0,
                diffusivity_m2_s: f64::MAX,
                time_step_s: 1.0 / 3.0,
                steps: u32::MAX,
            },
        },
        seed: u64::MAX,
        sea_level_m: -1.5e-320,
        stages: Vec::new(),
        ..world()
    }
}

/// The `f64` metadata of `world`, as bits, listed by hand from the `World`.
fn world_floats(world: &World) -> Vec<u64> {
    let e = &world.settings.erosion;
    let p = &world.params;
    [
        e.erodibility,
        e.discharge_exponent,
        e.diffusivity_m2_s,
        e.time_step_s,
        p.radius_m,
        p.surface_gravity_m_s2,
        p.rotation_period_s,
        p.obliquity_rad,
        p.orbital_period_s,
        p.orbital_eccentricity,
        p.perihelion_longitude_rad,
        p.stellar_flux_w_m2,
        p.ocean_volume_m3,
        p.surface_pressure_pa,
        p.continental_crust_fraction,
        world.sea_level_m,
    ]
    .map(f64::to_bits)
    .to_vec()
}

/// The `f64` values of parsed metadata, as bits, listed by hand.
fn json_floats(json: &WorldJson) -> Vec<u64> {
    let e = &json.settings.erosion;
    let p = &json.planet;
    [
        e.erodibility,
        e.discharge_exponent,
        e.diffusivity_m2_s,
        e.time_step_s,
        p.radius_m,
        p.surface_gravity_m_s2,
        p.rotation_period_s,
        p.obliquity_rad,
        p.orbital_period_s,
        p.orbital_eccentricity,
        p.perihelion_longitude_rad,
        p.stellar_flux_w_m2,
        p.ocean_volume_m3,
        p.surface_pressure_pa,
        p.continental_crust_fraction,
        json.sea_level_m,
    ]
    .map(f64::to_bits)
    .to_vec()
}

#[test]
fn world_json_parses_back_bit_exactly() {
    for world in [world(), awkward_world()] {
        let bytes = world_json(&world).unwrap();
        let back: WorldJson = serde_json::from_slice(&bytes).unwrap();
        assert_eq!(back, WorldJson::of(&world));
        assert_eq!(json_floats(&back), world_floats(&world));
        assert_eq!(back.format_version, 4);
        assert_eq!(back.generator_version, kuniumi_world::GENERATOR_VERSION);
        assert_eq!(back.seed, world.seed);
        let s = &back.settings;
        assert_eq!(
            (s.cells, s.plate_count, s.allow_partial, s.erosion.steps),
            (
                world.settings.cells,
                world.settings.plate_count,
                world.settings.allow_partial,
                world.settings.erosion.steps
            )
        );
        assert_eq!(s.to_settings(), world.settings);
        assert_eq!(back.stages, world.stages);
        let fields: Vec<_> = back
            .fields
            .iter()
            .map(|f| (f.name.as_str(), f.unit.as_str(), f.kind))
            .collect();
        assert_eq!(
            fields,
            [
                ("height", "m", KindJson::Quantity),
                ("zone", "", KindJson::Class)
            ]
        );
    }
    let back: WorldJson = serde_json::from_slice(&world_json(&awkward_world()).unwrap()).unwrap();
    assert_eq!(back.planet.rotation, RotationJson::TidallyLocked);
    assert_eq!(back.planet.tectonics, TectonicsJson::StagnantLid);
    assert_eq!((back.seed, back.settings.cells), (u64::MAX, u32::MAX));
    assert_eq!(back.settings.erosion.steps, u32::MAX);
}

/// Serialises `values` and parses them back with serde_json, and checks
/// that every value keeps its bits.
fn assert_round_trip(values: &[f64]) {
    let text = serde_json::to_string(values).unwrap();
    let back: Vec<f64> = serde_json::from_str(&text).unwrap();
    assert_eq!(back.len(), values.len());
    for (a, b) in values.iter().zip(&back) {
        assert_eq!(a.to_bits(), b.to_bits(), "{a:e} came back as {b:e}");
    }
}

#[test]
fn serde_json_round_trips_every_finite_f64() {
    let mut rng = Rng::from_seed(20_261_109);
    // Random bit patterns: every exponent equally likely, subnormals
    // included.
    let bits: Vec<f64> = (0..100_000)
        .map(|_| f64::from_bits(rng.next_u64()))
        .filter(|v| v.is_finite())
        .collect();
    assert!(bits.len() > 99_000);
    assert_round_trip(&bits);
    // Random numbers of the sizes a world holds, 1e-20 to 1e20.
    let mut decimals = Vec::new();
    let mut scale = 1e-20;
    for _ in 0..41 {
        decimals.extend((0..1_000).map(|_| rng.next_f64() * scale));
        scale *= 10.0;
    }
    assert_round_trip(&decimals);
    // Hard cases of decimal conversion.
    assert_round_trip(&[
        0.0,
        -0.0,
        5e-324,
        -5e-324,
        2.225_073_858_507_201e-308,
        2.225_073_858_507_201_4e-308,
        f64::MIN_POSITIVE,
        f64::MAX,
        f64::MIN,
        f64::EPSILON,
        0.1,
        0.3,
        0.1 + 0.2,
        1.0 / 3.0,
        1e23,
        8.988_465_674_311_58e307,
        9_007_199_254_740_992.0,
        9_007_199_254_740_994.0,
        123_456_789_012_345_680.0,
        1e-7,
        1e21,
        1e22,
    ]);
}

/// `world.json` for Earth after the Grid stage, written by hand: the decimal
/// literals of `PlanetParams::earth()` and `ErosionSettings::default()`,
/// the obliquity `23.44 · (π/180)` as `repr(23.44 * (math.pi / 180))`
/// prints it in Python 3 (0.40910517666747087), the longitude of
/// perihelion as `repr(102.94 * (math.pi / 180))`, and the default
/// diffusivity and time step as `repr(5e-3 / 31557600.0)` and
/// `repr(1e4 * 31557600.0)` print them. It is the example in the module
/// documentation.
const EARTH_JSON: &str = r#"{
  "format_version": 4,
  "generator_version": "0.0.1",
  "seed": 1,
  "settings": {
    "cells": 20000,
    "plate_count": 12,
    "allow_partial": true,
    "erosion": {
      "erodibility": 7.5e-10,
      "discharge_exponent": 0.45,
      "diffusivity_m2_s": 1.5844043907014476e-10,
      "time_step_s": 315576000000.0,
      "steps": 1000
    }
  },
  "planet": {
    "radius_m": 6371000.0,
    "surface_gravity_m_s2": 9.80665,
    "rotation_period_s": 86164.0905,
    "obliquity_rad": 0.40910517666747087,
    "orbital_period_s": 31558149.8,
    "orbital_eccentricity": 0.0167,
    "perihelion_longitude_rad": 1.7966419320029627,
    "stellar_flux_w_m2": 1361.0,
    "ocean_volume_m3": 1.335e+18,
    "surface_pressure_pa": 101325.0,
    "continental_crust_fraction": 0.4,
    "rotation": "free",
    "tectonics": "plate_tectonics"
  },
  "sea_level_m": 0.0,
  "stages": [
    "grid"
  ],
  "fields": [
    {
      "name": "cell_area",
      "unit": "sr",
      "kind": "quantity"
    },
    {
      "name": "neighbor_count",
      "unit": "",
      "kind": "class"
    }
  ]
}
"#;

#[test]
fn earth_world_json_has_the_documented_layout() {
    // The metadata does not depend on the grid's size, so a small grid
    // stands in for the 20,000 cells of the settings.
    let mut world = World {
        settings: Settings {
            cells: 20_000,
            plate_count: 12,
            allow_partial: true,
            erosion: ErosionSettings::default(),
        },
        fields: Vec::new(),
        ..world()
    };
    let n = world.grid.len();
    world.add_field("cell_area", "sr", vec![0.0; n]).unwrap();
    world
        .add_class_field("neighbor_count", vec![6.0; n])
        .unwrap();
    assert_eq!(kuniumi_world::GENERATOR_VERSION, "0.0.1");
    let text = String::from_utf8(world_json(&world).unwrap()).unwrap();
    assert_eq!(text, EARTH_JSON);
    // The module documentation shows the same text.
    let doc = include_str!("../world_dir.rs");
    let documented: String = doc
        .lines()
        .skip_while(|l| *l != "//! ```json")
        .skip(1)
        .take_while(|l| *l != "//! ```")
        .map(|l| format!("{}\n", l.strip_prefix("//! ").unwrap_or("").trim_end()))
        .collect();
    assert_eq!(documented, EARTH_JSON);
}

/// `world.json` of [`awkward_world`], written by hand: each number as
/// Python 3's `repr` prints it, the shortest decimal that parses back to the
/// same bits (`repr(5e-324)`, `repr(1.7976931348623157e308)`, `repr(1e23)`,
/// `repr(1 / 3)`, `repr(0.1 + 0.2)`, `repr(-1.5e-320)`, …). It checks
/// serde_json's float formatter, whose digit writers differ between
/// targets, on subnormals, 17-digit values, the extremes and −0.
const AWKWARD_JSON: &str = r#"{
  "format_version": 4,
  "generator_version": "0.0.1",
  "seed": 18446744073709551615,
  "settings": {
    "cells": 4294967295,
    "plate_count": 4294967294,
    "allow_partial": false,
    "erosion": {
      "erodibility": 5e-324,
      "discharge_exponent": -0.0,
      "diffusivity_m2_s": 1.7976931348623157e+308,
      "time_step_s": 0.3333333333333333,
      "steps": 4294967295
    }
  },
  "planet": {
    "radius_m": 2.2250738585072014e-308,
    "surface_gravity_m_s2": 5e-324,
    "rotation_period_s": 1.7976931348623157e+308,
    "obliquity_rad": 0.1,
    "orbital_period_s": 1e+23,
    "orbital_eccentricity": 0.3333333333333333,
    "perihelion_longitude_rad": -1e+300,
    "stellar_flux_w_m2": 2.225073858507201e-308,
    "ocean_volume_m3": 9007199254740991.0,
    "surface_pressure_pa": -0.0,
    "continental_crust_fraction": 0.30000000000000004,
    "rotation": "tidally_locked",
    "tectonics": "stagnant_lid"
  },
  "sea_level_m": -1.5e-320,
  "stages": [],
  "fields": [
    {
      "name": "height",
      "unit": "m",
      "kind": "quantity"
    },
    {
      "name": "zone",
      "unit": "",
      "kind": "class"
    }
  ]
}
"#;

#[test]
fn awkward_world_json_has_the_shortest_forms() {
    let text = String::from_utf8(world_json(&awkward_world()).unwrap()).unwrap();
    assert_eq!(text, AWKWARD_JSON);
}

#[test]
fn the_same_world_gives_the_same_bytes() {
    for world in [world(), awkward_world()] {
        assert_eq!(
            world_json(&world).unwrap(),
            world_json(&world.clone()).unwrap()
        );
    }
    let (a, b) = (scratch("same_bytes_a"), scratch("same_bytes_b"));
    let world = world();
    let written_a = write_maps_of_width(&a, &world, 64).unwrap();
    let written_b = write_maps_of_width(&b, &world.clone(), 64).unwrap();
    assert_eq!(written_a.len(), 6);
    for (pa, pb) in written_a.iter().zip(&written_b) {
        assert_eq!(pa.file_name(), pb.file_name());
        assert_eq!(read(pa), read(pb), "{}", pa.display());
    }
    remove(&a);
    remove(&b);
}

/// A change that makes a world unwritable.
type Spoil = fn(&mut World);

#[test]
fn rejects_non_finite_metadata() {
    let cases: [(Spoil, &str); 9] = [
        (|w| w.sea_level_m = f64::NAN, "`sea_level_m`"),
        (
            |w| w.settings.erosion.erodibility = f64::NAN,
            "`settings.erosion.erodibility`",
        ),
        (
            |w| w.settings.erosion.discharge_exponent = f64::NEG_INFINITY,
            "-inf of `settings.erosion.discharge_exponent`",
        ),
        (
            |w| w.settings.erosion.diffusivity_m2_s = f64::INFINITY,
            "inf of `settings.erosion.diffusivity_m2_s`",
        ),
        (
            |w| w.settings.erosion.time_step_s = f64::NAN,
            "`settings.erosion.time_step_s`",
        ),
        (
            |w| w.sea_level_m = f64::NEG_INFINITY,
            "-inf of `sea_level_m`",
        ),
        (
            |w| w.params.radius_m = f64::INFINITY,
            "inf of `planet.radius_m`",
        ),
        (
            |w| w.params.continental_crust_fraction = f64::NAN,
            "`planet.continental_crust_fraction`",
        ),
        (
            |w| w.params.perihelion_longitude_rad = f64::NEG_INFINITY,
            "-inf of `planet.perihelion_longitude_rad`",
        ),
    ];
    for (spoil, key) in cases {
        let mut world = world();
        spoil(&mut world);
        let err = world_json(&world).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidInput);
        let message = err.to_string();
        assert!(message.starts_with("world.json cannot hold the non-finite value "));
        assert!(message.ends_with(key), "{message}");

        let dir = scratch("non_finite");
        let err = write(&dir, &world).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidInput);
        assert!(!dir.exists(), "nothing is written");
    }
    // The first in the order written: the settings before the planet.
    let mut world = world();
    world.params.radius_m = f64::NAN;
    world.settings.erosion.time_step_s = f64::INFINITY;
    let message = world_json(&world).unwrap_err().to_string();
    assert!(
        message.ends_with("inf of `settings.erosion.time_step_s`"),
        "{message}"
    );
}

#[test]
fn rejects_bad_fields_before_writing_anything() {
    let cases: [(Spoil, &str); 7] = [
        (
            |w| w.fields[1].name = "Zone".into(),
            "`Zone` is not a valid field name",
        ),
        (
            // Reserved in a world without a `surface_height` too (ticket
            // M2-19).
            |w| w.fields[1].name = "hypsometry".into(),
            HYPSOMETRY_RESERVED,
        ),
        (
            |w| w.fields[0].name = "../height".into(),
            "`../height` is not a valid field name",
        ),
        (
            |w| w.fields[1].name = "height".into(),
            "the field `height` appears twice",
        ),
        (
            |w| {
                w.fields[1].values.pop();
            },
            "field `zone` has 999 values for 1000 cells",
        ),
        (
            |w| w.fields[1].values[7] = 0.5,
            "class field `zone` has the value 0.5 at cell 7, not a class id",
        ),
        (
            |w| {
                // A quantity turned into a class field after it was added.
                let field = &mut w.fields[0];
                field.kind = FieldKind::Class;
                field.values.fill(2.0);
                field.values[4] = -1.0;
            },
            "class field `height` has the value -1 at cell 4, not a class id",
        ),
    ];
    for (spoil, message) in cases {
        let mut world = world();
        spoil(&mut world);
        let dir = scratch("bad_fields");
        let err = write(&dir, &world).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidInput);
        assert_eq!(err.to_string(), message);
        assert!(!dir.exists(), "nothing is written");
    }
}

#[test]
fn writes_every_field_and_world_json_last() {
    let dir = scratch("every_field");
    // A file of another program in the directory is left alone.
    std::fs::create_dir_all(&dir).unwrap();
    std::fs::write(dir.join("notes.txt"), "kept").unwrap();
    let world = world();
    let written = write_maps_of_width(&dir, &world, 64).unwrap();
    let names: Vec<_> = written
        .iter()
        .map(|p| p.file_name().unwrap().to_str().unwrap())
        .collect();
    assert_eq!(
        names,
        [
            "grid.kgrid",
            "height.kfield",
            "height.png",
            "zone.kfield",
            "zone.png",
            "world.json"
        ]
    );
    assert!(written.iter().all(|p| p.parent() == Some(dir.as_path())));
    let grid = grid_file::read_grid(read(&dir.join(GRID_FILE)).as_slice()).unwrap();
    assert_eq!(grid, world.grid);
    for field in &world.fields {
        let back =
            field_file::read_field(read(&dir.join(format!("{}.kfield", field.name))).as_slice())
                .unwrap();
        assert_eq!((&back.name, &back.unit), (&field.name, &field.unit));
        let bits = |f: &Field| f.values.iter().map(|v| v.to_bits()).collect::<Vec<_>>();
        assert_eq!(bits(&back), bits(field));
    }
    assert_eq!(read(&dir.join("world.json")), world_json(&world).unwrap());
    assert_eq!(read(&dir.join("notes.txt")), b"kept");
    // Writing again replaces the files.
    let mut changed = world.clone();
    changed.seed = 2;
    write_maps_of_width(&dir, &changed, 64).unwrap();
    assert_eq!(read(&dir.join("world.json")), world_json(&changed).unwrap());
    remove(&dir);
}

#[test]
fn file_system_errors_name_the_path() {
    // The directory is a file.
    let file = scratch("dir_is_a_file");
    std::fs::write(&file, "not a directory").unwrap();
    let err = write_maps_of_width(&file, &world(), 64).unwrap_err();
    assert!(
        err.to_string()
            .starts_with(&format!("{}: ", file.display())),
        "{err}"
    );
    std::fs::remove_file(&file).unwrap();
    // A field's file is a directory.
    let dir = scratch("kfield_is_a_dir");
    std::fs::create_dir_all(dir.join("zone.kfield")).unwrap();
    let err = write_maps_of_width(&dir, &world(), 64).unwrap_err();
    let path = dir.join("zone.kfield");
    assert!(
        err.to_string()
            .starts_with(&format!("{}: ", path.display())),
        "{err}"
    );
    // The files before it are written, world.json is not.
    assert!(dir.join("height.png").exists());
    assert!(!dir.join("world.json").exists());
    remove(&dir);
    // The grid's file is a directory: nothing after it is written.
    let dir = scratch("kgrid_is_a_dir");
    std::fs::create_dir_all(dir.join(GRID_FILE)).unwrap();
    let err = write_maps_of_width(&dir, &world(), 64).unwrap_err();
    let path = dir.join(GRID_FILE);
    assert!(
        err.to_string()
            .starts_with(&format!("{}: ", path.display())),
        "{err}"
    );
    assert!(!dir.join("height.kfield").exists());
    remove(&dir);
    // A rewrite that fails halfway leaves no world.json of the earlier
    // world beside the new world's files.
    let dir = scratch("failed_rewrite");
    write_maps_of_width(&dir, &world(), 64).unwrap();
    assert!(dir.join("world.json").exists());
    std::fs::remove_file(dir.join("zone.kfield")).unwrap();
    std::fs::create_dir(dir.join("zone.kfield")).unwrap();
    let mut other = world();
    other.seed = 2;
    other.fields[0].values[0] = 0.25;
    assert!(write_maps_of_width(&dir, &other, 64).is_err());
    assert!(!dir.join("world.json").exists());
    let height = field_file::read_field(read(&dir.join("height.kfield")).as_slice()).unwrap();
    assert_eq!(height.values[0].to_bits(), 0.25f64.to_bits());
    remove(&dir);
}

#[test]
fn maps_use_the_categorical_palette_for_classes() {
    // Other fields keep their maps, with or without a surface height.
    for world in [world(), world_with_relief()] {
        let metrics = ElevationMetrics::of_world(&world);
        let expected = [
            (&world.fields[0], Colormap::Viridis),
            (&world.fields[1], Colormap::Categorical),
        ];
        for (field, colormap) in expected {
            assert_eq!(
                map_colours(field, &world, metrics.as_ref()),
                (colormap, Range::Auto),
                "{}",
                field.name
            );
        }
        let dir = scratch("palettes");
        write_maps_of_width(&dir, &world, 64).unwrap();
        for (field, colormap) in expected {
            let (width, height, pixels) = decode_png(&dir.join(format!("{}.png", field.name)));
            assert_eq!((width, height), (64, 32));
            let reference = render::equirect(&world.grid, &field.values, colormap, Range::Auto, 64);
            assert_eq!(pixels, reference, "{}", field.name);
        }
        // Class ids 0 to 4 take the first five colours of the palette.
        let (_, _, pixels) = decode_png(&dir.join("zone.png"));
        let mut seen: Vec<Rgb8> = pixels.clone();
        seen.sort_unstable();
        seen.dedup();
        let mut first_five = crate::colormap::CATEGORICAL[..5].to_vec();
        first_five.sort_unstable();
        assert_eq!(seen, first_five);
        remove(&dir);
    }
    // A class field is categorical whatever its name.
    let mut world = world();
    world.fields.clear();
    let classes = (0..world.grid.len()).map(|i| (i % 3) as f64).collect();
    world.add_class_field(SURFACE_HEIGHT, classes).unwrap();
    let metrics = ElevationMetrics::of_world(&world);
    assert_eq!(
        map_colours(&world.fields[0], &world, metrics.as_ref()),
        (Colormap::Categorical, Range::Auto)
    );
}

/// The colours of the six cells of [`octahedron_world`] (sea level
/// −500 m) on the elevation map over −11.5 km to +8.5 km (the chart's
/// axis, −11 km to +9 km above the sea level), by hand from the palettes
/// of `colormap` (each channel `a + f (b − a)`, rounded half away from
/// zero):
///
/// - +700 m (cells 0 and 5): land, `t = 700/9,000`, `4t = 0.3111` on the
///   first segment, green `[46, 120, 56]` to `[138, 172, 82]`:
///   74.62, 136.18, 64.09;
/// - −4,500 m (cell 1): ocean, `t = 6,500/11,000`, `3t = 1.7727`, `f =
///   0.7727` from `[22, 70, 140]` to `[58, 128, 196]`: 49.82, 114.82,
///   183.27;
/// - +1,500 m (cell 2): `4t = 2/3` on the first segment: 107.33, 154.67,
///   73.33;
/// - −3,500 m (cell 3): `t = 7,500/11,000`, `3t = 2.0455`, `f = 0.0455`
///   from `[58, 128, 196]` to `[150, 200, 235]`: 62.18, 131.27, 197.77;
/// - +3,000 m (cell 4): `4t = 4/3`, `f = 1/3` from `[138, 172, 82]` to
///   `[212, 190, 122]`: 162.67, 178, 95.33.
const OCTAHEDRON_COLOURS: [Rgb8; 6] = [
    [75, 136, 64],
    [50, 115, 183],
    [107, 155, 73],
    [62, 131, 198],
    [163, 178, 95],
    [75, 136, 64],
];

#[test]
fn surface_height_takes_the_elevation_map() {
    let world = octahedron_world();
    let metrics = ElevationMetrics::of_world(&world).unwrap();
    // The chart's axis, −11 km to +9 km, shifted by the sea level, −500 m.
    let expected = (
        Colormap::Elevation { sea_level: -500.0 },
        Range::Fixed {
            min: -11_500.0,
            max: 8_500.0,
        },
    );
    assert_eq!(
        map_colours(&world.fields[0], &world, Some(&metrics)),
        expected
    );
    let dir = scratch("elevation_map");
    write_maps_of_width(&dir, &world, 64).unwrap();
    let (width, height, pixels) = decode_png(&dir.join("surface_height.png"));
    assert_eq!((width, height), (64, 32));
    let reference = render::equirect(
        &world.grid,
        &world.fields[0].values,
        expected.0,
        expected.1,
        64,
    );
    assert_eq!(pixels, reference);
    // The five colours of the six elevations, and each cell's colour at its
    // vertex: +X at (32, 16), −X at (63, 16), +Y at (48, 16), −Y at
    // (16, 16), +Z in the top row and −Z in the bottom row (pixel centres
    // 2.8° from the vertices; `render`'s layout).
    let mut seen = pixels.clone();
    seen.sort_unstable();
    seen.dedup();
    let mut colours = OCTAHEDRON_COLOURS.to_vec();
    colours.sort_unstable();
    colours.dedup();
    assert_eq!(seen, colours);
    let at = |x: usize, y: usize| pixels[y * 64 + x];
    let vertices = [(32, 16), (63, 16), (48, 16), (16, 16), (0, 0), (0, 31)];
    for (cell, (x, y)) in vertices.into_iter().enumerate() {
        assert_eq!(at(x, y), OCTAHEDRON_COLOURS[cell], "cell {cell}");
    }
    remove(&dir);
}

/// Whether `colour` is one of the ocean colours of the elevation map: all
/// have blue above green (by at least 35 at the control points), all land
/// colours blue below green (by at least 12).
fn is_ocean_colour(colour: Rgb8) -> bool {
    colour[2] > colour[1]
}

#[test]
fn land_and_sea_split_at_the_sea_level() {
    let relief = world_with_relief();
    // No water: the sea level at the lowest height, which every cell below
    // −2 km takes, so many cells lie exactly at the sea level.
    let mut dry = world_with_relief();
    for h in &mut dry.fields[2].values {
        if *h < -2_000.0 {
            *h = -2_000.0;
        }
    }
    dry.sea_level_m = -2_000.0;
    // A water world: the sea level 1 km above the highest height.
    let mut wet = world_with_relief();
    wet.sea_level_m = 5_400.0;
    let cells = render::equirect_cells(&relief.grid, 64);
    let mut counts = Vec::new();
    for world in [relief, dry, wet] {
        let dir = scratch("land_and_sea");
        write_maps_of_width(&dir, &world, 64).unwrap();
        let (_, _, pixels) = decode_png(&dir.join("surface_height.png"));
        let heights = &world.fields[2].values;
        let z = world.sea_level_m;
        for (p, (&cell, &colour)) in cells.iter().zip(&pixels).enumerate() {
            let h = heights[cell as usize];
            assert_eq!(
                is_ocean_colour(colour),
                h < z,
                "pixel {p}: cell {cell} at {h} m, sea level {z} m, {colour:?}"
            );
        }
        let at_sea_level = cells
            .iter()
            .filter(|&&c| heights[c as usize].to_bits() == z.to_bits())
            .count();
        let ocean = pixels.iter().filter(|&&c| is_ocean_colour(c)).count();
        counts.push((ocean, at_sea_level));
        remove(&dir);
    }
    // Land and sea in the first world, cells at the sea level drawn as land
    // in the second, only sea in the third.
    assert!(counts[0].0 > 0 && counts[0].0 < cells.len(), "{counts:?}");
    assert!(counts[1].0 == 0 && counts[1].1 > 100, "{counts:?}");
    assert_eq!(counts[2].0, cells.len(), "{counts:?}");
}

#[test]
fn an_elevation_has_one_colour_in_every_earth_like_world() {
    let world = world_with_relief();
    let z = world.sea_level_m;
    // The same world with one cell 10 km below the sea level and one 8 km
    // above it: new extremes, still inside Earth's range.
    let mut other = world.clone();
    other.fields[2].values[10] = z - 10_000.0;
    other.fields[2].values[20] = z + 8_000.0;
    let colours = |world: &World| {
        let metrics = ElevationMetrics::of_world(world).unwrap();
        let (colormap, range) = map_colours(&world.fields[2], world, Some(&metrics));
        (
            colormap.colours(&world.fields[2].values, range),
            colormap.colours(&world.fields[2].values, Range::Auto),
        )
    };
    let ((fixed, auto), (other_fixed, other_auto)) = (colours(&world), colours(&other));
    let others = |c: &[Rgb8]| -> Vec<Rgb8> {
        c.iter()
            .enumerate()
            .filter(|&(cell, _)| cell != 10 && cell != 20)
            .map(|(_, &colour)| colour)
            .collect()
    };
    assert_eq!(others(&fixed), others(&other_fixed));
    // On each world's own range most colours would change.
    let changed = others(&auto)
        .iter()
        .zip(&others(&other_auto))
        .filter(|(a, b)| a != b)
        .count();
    assert!(changed > 900, "{changed}");
    // And so in the written maps.
    let (a, b) = (scratch("one_colour_a"), scratch("one_colour_b"));
    write_maps_of_width(&a, &world, 64).unwrap();
    write_maps_of_width(&b, &other, 64).unwrap();
    let (_, _, pixels) = decode_png(&a.join("surface_height.png"));
    let (_, _, other_pixels) = decode_png(&b.join("surface_height.png"));
    let cells = render::equirect_cells(&world.grid, 64);
    for (p, &cell) in cells.iter().enumerate() {
        if cell != 10 && cell != 20 {
            assert_eq!(pixels[p], other_pixels[p], "pixel {p}, cell {cell}");
        }
    }
    remove(&a);
    remove(&b);
}

#[test]
fn the_range_falls_back_to_the_field_s_own() {
    let world = world_with_relief();
    let elevation = Colormap::Elevation { sea_level: -700.0 };
    assert_eq!(
        map_colours(&world.fields[2], &world, None),
        (elevation, Range::Auto)
    );
    // A sea level of 1.7e308 m: the chart's axis ends at +1e307 m above
    // it, beyond the largest double.
    let mut far = world_with_relief();
    far.sea_level_m = 1.7e308;
    let metrics = ElevationMetrics::of_world(&far).unwrap();
    let axes = hypsometry_axes(&metrics);
    assert!((far.sea_level_m + axes.x_max).is_infinite(), "{axes:?}");
    assert_eq!(
        map_colours(&far.fields[2], &far, Some(&metrics)),
        (Colormap::Elevation { sea_level: 1.7e308 }, Range::Auto)
    );
    let dir = scratch("far_sea_level");
    write_maps_of_width(&dir, &far, 64).unwrap();
    let (_, _, pixels) = decode_png(&dir.join("surface_height.png"));
    let reference = render::equirect(
        &far.grid,
        &far.fields[2].values,
        Colormap::Elevation { sea_level: 1.7e308 },
        Range::Auto,
        64,
    );
    assert_eq!(pixels, reference);
    assert!(pixels.iter().all(|&c| is_ocean_colour(c)));
    remove(&dir);
    // A sea level of −1.7e308 m: the chart's axis starts at −1e307 m below
    // it, beyond the most negative double.
    let mut deep = world_with_relief();
    deep.sea_level_m = -1.7e308;
    let metrics = ElevationMetrics::of_world(&deep).unwrap();
    let axes = hypsometry_axes(&metrics);
    assert!((deep.sea_level_m + axes.x_min).is_infinite(), "{axes:?}");
    assert!((deep.sea_level_m + axes.x_max).is_finite(), "{axes:?}");
    assert_eq!(
        map_colours(&deep.fields[2], &deep, Some(&metrics)),
        (
            Colormap::Elevation {
                sea_level: -1.7e308
            },
            Range::Auto
        )
    );
    let dir = scratch("deep_sea_level");
    write_maps_of_width(&dir, &deep, 64).unwrap();
    let (_, _, pixels) = decode_png(&dir.join("surface_height.png"));
    let reference = render::equirect(
        &deep.grid,
        &deep.fields[2].values,
        Colormap::Elevation {
            sea_level: -1.7e308,
        },
        Range::Auto,
        64,
    );
    assert_eq!(pixels, reference);
    assert!(pixels.iter().all(|&c| !is_ocean_colour(c)));
    remove(&dir);
}

#[test]
fn writes_maps_of_2048_by_1024() {
    let dir = scratch("full_width");
    let world = World {
        fields: world().fields[..1].to_vec(),
        ..world()
    };
    write(&dir, &world).unwrap();
    let (width, height, pixels) = decode_png(&dir.join("height.png"));
    assert_eq!((width, height, MAP_WIDTH_PX), (2048, 1024, 2048));
    let reference = render::equirect(
        &world.grid,
        &world.fields[0].values,
        Colormap::Viridis,
        Range::Auto,
        2048,
    );
    assert_eq!(pixels, reference);
    remove(&dir);
}

/// A world on the octahedron grid whose six cells stand at −4,500, −3,500,
/// +700 (twice), +1,500 and +3,000 m above a sea level of −500 m.
pub(super) fn octahedron_world() -> World {
    let mut world = World {
        grid: SphereGrid::octahedron(),
        sea_level_m: -500.0,
        stages: vec!["grid".into(), "tectonics".into(), "crust".into()],
        fields: Vec::new(),
        ..world()
    };
    world
        .add_field(
            "surface_height",
            "m",
            vec![200.0, -5_000.0, 1_000.0, -4_000.0, 2_500.0, 200.0],
        )
        .unwrap();
    world
}

/// The text of `metrics.txt` for [`octahedron_world`], written by hand from
/// the definitions of `kuniumi_sim::metrics`: land 4/6, its mean
/// (700 + 700 + 1,500 + 3,000)/4 m, the ocean (4,500 + 3,500)/2 m deep;
/// the percentiles from the sorted elevations −4,500, −3,500, 700, 700,
/// 1,500, 3,000 (1/6 each); bins −45 to 30, 1/6 at −4,500, −3,500, +1,500
/// and +3,000 m and 2/6 at +700 m; the ocean mode the lower of the two
/// ocean bins, the valley the lowest empty bin above it; the chart from
/// −11 km to +9 km and up to 7 steps of 0.05.
fn octahedron_metrics_txt() -> String {
    let mut text = String::from(
        "# kuniumi metrics: land fraction and hypsometry (ticket M2-06)
# Area-weighted over the cells. Elevations are heights above the sea level
# (m), negative below it; cells below the sea level are ocean, all others
# land. Earth's values are approximate (docs/design-ja.md, 検証方法).
generator_version = 0.0.1
seed = 1
cells = 6
sea_level_m = -500.000

# Earth: land fraction about 0.29, mean land elevation about 800 m, mean
# ocean depth about 3700 m, highest about 8800 m, lowest about -11000 m.
land_fraction = 0.666667
mean_land_elevation_m = 1475.0
mean_ocean_depth_m = 4000.0
highest_m = 3000.0
lowest_m = -4500.0

# Percentiles: the elevation of the first cell, from the lowest, at which
# the share of the surface at or below it reaches the given percentage.
percentile_1_m = -4500.0
percentile_5_m = -4500.0
percentile_10_m = -4500.0
percentile_25_m = -3500.0
percentile_50_m = 700.0
percentile_75_m = 1500.0
percentile_90_m = 3000.0
percentile_95_m = 3000.0
percentile_99_m = 3000.0

# Peaks: the centres of the bins with the largest share below and above
# the sea level (Earth: about -4500 m and +100 m), and the valley, the bin
# between them with the smallest share. The peaks are separate when the
# valley holds at most 0.5 of the smaller one.
ocean_mode_m = -4450.0
ocean_mode_fraction = 0.166667
land_mode_m = 750.0
land_mode_fraction = 0.333333
valley_m = -4350.0
valley_fraction = 0.000000
separate_peaks = yes

# hypsometry.png: the histogram as bars, blue below the sea level and green
# to brown above it; a black line at the sea level and red lines at Earth's
# modes; grey lines every chart_x_step_m and every chart_y_step_fraction.
chart_x_min_m = -11000
chart_x_max_m = 9000
chart_x_step_m = 1000
chart_y_max_fraction = 0.350
chart_y_step_fraction = 0.050

# Histogram: bins of bin_width_m from the lowest elevation to the highest:
# lower edge (m), upper edge (m), share of the surface.
bin_width_m = 100
",
    );
    for k in -45..=30 {
        let share = match k {
            -45 | -35 | 15 | 30 => "0.16666667",
            7 => "0.33333333",
            _ => "0.00000000",
        };
        text.push_str(&format!("{} {} {share}\n", 100 * k, 100 * (k + 1)));
    }
    text
}

#[test]
fn metrics_txt_by_hand() {
    let world = octahedron_world();
    let metrics = ElevationMetrics::of_world(&world).unwrap();
    assert_eq!(kuniumi_world::GENERATOR_VERSION, "0.0.1");
    let text = metrics_txt(&world, &metrics, None);
    let expected = octahedron_metrics_txt();
    for (k, (a, b)) in text.lines().zip(expected.lines()).enumerate() {
        assert_eq!(a, b, "line {}", k + 1);
    }
    assert_eq!(text, expected);
}

#[test]
fn metrics_txt_without_land_or_peaks() {
    // A water world: no land, so no mean elevation, land mode or valley.
    let mut world = octahedron_world();
    world.sea_level_m = 10_000.0;
    let text = metrics_txt(&world, &ElevationMetrics::of_world(&world).unwrap(), None);
    for line in [
        "land_fraction = 0.000000",
        "mean_land_elevation_m = none",
        // (9,800 · 2 + 15,000 + 9,000 + 14,000 + 7,500)/6 m.
        "mean_ocean_depth_m = 10850.0",
        "land_mode_m = none",
        "land_mode_fraction = none",
        "valley_m = none",
        "valley_fraction = none",
        "separate_peaks = no",
    ] {
        assert!(text.lines().any(|l| l == line), "`{line}` missing:\n{text}");
    }
    // No water: no ocean.
    world.sea_level_m = -5_000.0;
    let text = metrics_txt(&world, &ElevationMetrics::of_world(&world).unwrap(), None);
    for line in [
        "land_fraction = 1.000000",
        "mean_ocean_depth_m = none",
        "ocean_mode_m = none",
        "ocean_mode_fraction = none",
        "separate_peaks = no",
    ] {
        assert!(text.lines().any(|l| l == line), "`{line}` missing:\n{text}");
    }
}

#[test]
fn chart_axes() {
    // Earth's range for the octahedron world, steps of 1 km and 0.05.
    let metrics = ElevationMetrics::of_world(&octahedron_world()).unwrap();
    let axes = hypsometry_axes(&metrics);
    assert_eq!(
        (axes.x_min, axes.x_max, axes.x_step, axes.y_step),
        (-11_000.0, 9_000.0, 1_000.0, 0.05)
    );
    assert!((axes.y_max - 0.35).abs() < 1e-15, "{}", axes.y_max);
    // From −30 km: 41 lines of 1 km would be too many, 21 of 2 km do.
    let deep = ElevationMetrics::new(&[1.0, 3.0], &[-30_000.0, 1_500.0], 0.0);
    let axes = hypsometry_axes(&deep);
    assert_eq!(
        (axes.x_min, axes.x_max, axes.x_step, axes.y_max, axes.y_step),
        (-30_000.0, 10_000.0, 2_000.0, 0.8, 0.1)
    );
    // A sea level 1,000 km up: steps of 50 km; one bin with everything.
    let wet = ElevationMetrics::new(&[1.0], &[-8_000.0], 1.0e6);
    let axes = hypsometry_axes(&wet);
    assert_eq!(
        (axes.x_min, axes.x_max, axes.x_step, axes.y_max, axes.y_step),
        (-1_050_000.0, 50_000.0, 50_000.0, 1.0, 0.1)
    );
    // A sea level near the largest double: finite axes (rounded, so far
    // out), and a chart.
    let far = ElevationMetrics::new(&[1.0], &[-8_000.0], 1.7e308);
    let axes = hypsometry_axes(&far);
    assert!(axes.x_min.is_finite() && axes.x_max.is_finite(), "{axes:?}");
    assert!(
        (axes.x_min / -1.7e308 - 1.0).abs() < 1e-15 && axes.x_max >= earth::HIGHEST_M,
        "{axes:?}"
    );
    let pixels = hypsometry_chart(&far);
    assert_eq!(pixels.len(), (CHART_WIDTH_PX * CHART_HEIGHT_PX) as usize);
    // The markers at the sea level and at Earth's modes, 1e307 m from the
    // right edge, fall into one plot column near it, drawn top to bottom.
    let m = chart::MARGIN_PX as usize;
    let row = |r: usize| &pixels[(m + r) * CHART_WIDTH_PX as usize..][m..m + 1_000];
    let marked: Vec<usize> = (0..1_000)
        .filter(|&c| (0..500).all(|r| row(r)[c] != chart::BACKGROUND && row(r)[c] != chart::GRID))
        .collect();
    assert_eq!(marked.len(), 1, "{marked:?}");
    assert!(marked[0] > 900, "{marked:?}");
    // Little area in the fullest bin: the smallest step.
    let flat_area = vec![1.0; 2_000];
    let flat_height: Vec<f64> = (0..2_000).map(|k| f64::from(k) * 100.0 + 50.0).collect();
    let flat = ElevationMetrics::new(&flat_area, &flat_height, 0.0);
    let axes = hypsometry_axes(&flat);
    assert_eq!((axes.y_max, axes.y_step), (0.001, 0.001));
}

#[test]
fn hypsometry_chart_draws_the_bins_and_the_markers() {
    let metrics = ElevationMetrics::of_world(&octahedron_world()).unwrap();
    let pixels = hypsometry_chart(&metrics);
    assert_eq!(pixels.len(), (CHART_WIDTH_PX * CHART_HEIGHT_PX) as usize);
    let axes = hypsometry_axes(&metrics);
    let bars: Vec<chart::Bar> = (0..metrics.hypsometry.shares().len())
        .map(|bin| {
            let hyps = &metrics.hypsometry;
            chart::Bar {
                lower: hyps.lower_edge_m(bin),
                upper: hyps.upper_edge_m(bin),
                value: hyps.shares()[bin],
                colour: Colormap::Elevation { sea_level: 0.0 }.colour(
                    hyps.centre_m(bin),
                    axes.x_min,
                    axes.x_max,
                ),
            }
        })
        .collect();
    let markers = [
        chart::Marker {
            x: 0.0,
            colour: SEA_LEVEL_MARKER,
        },
        chart::Marker {
            x: earth::OCEAN_MODE_M,
            colour: EARTH_MARKER,
        },
        chart::Marker {
            x: earth::LAND_MODE_M,
            colour: EARTH_MARKER,
        },
    ];
    let reference = chart::histogram(&bars, &axes, &markers, CHART_WIDTH_PX, CHART_HEIGHT_PX);
    assert_eq!(pixels, reference);
    // By hand: 20 m per column from −11 km, so the sea level in plot column
    // 550, Earth's modes in 325 and 555; the bin at +700 m in columns 585 to
    // 589, 2/6 / 0.35 · 500 = 476 rows high, in the land colour at +750 m.
    let m = chart::MARGIN_PX as usize;
    let at = |c: usize, r: usize| pixels[(m + r) * CHART_WIDTH_PX as usize + m + c];
    for r in 0..500 {
        assert_eq!(at(550, r), SEA_LEVEL_MARKER, "row {r}");
        assert_eq!(at(325, r), EARTH_MARKER, "row {r}");
        assert_eq!(at(555, r), EARTH_MARKER, "row {r}");
    }
    let land = Colormap::Elevation { sea_level: 0.0 }.colour(750.0, -11_000.0, 9_000.0);
    for c in 585..590 {
        assert_eq!(at(c, 499), land, "column {c}");
        assert_eq!(at(c, 500 - 476), land, "column {c}");
        assert_ne!(at(c, 500 - 477), land, "column {c}");
    }
    assert_ne!(at(584, 499), land);
    assert_ne!(at(590, 499), land);
}

/// [`world`] with a surface 5 km · z − 600 m high and the sea level at
/// −700 m.
pub(super) fn world_with_relief() -> World {
    let mut world = world();
    let heights = world
        .grid
        .positions()
        .iter()
        .map(|p| 5_000.0 * p.z - 600.0)
        .collect();
    world.add_field("surface_height", "m", heights).unwrap();
    world.sea_level_m = -700.0;
    world
}

#[test]
fn writes_metrics_and_the_chart_before_world_json() {
    let dir = scratch("metrics");
    let world = world_with_relief();
    let written = write_maps_of_width(&dir, &world, 64).unwrap();
    let names: Vec<_> = written
        .iter()
        .map(|p| p.file_name().unwrap().to_str().unwrap())
        .collect();
    assert_eq!(
        names,
        [
            "grid.kgrid",
            "height.kfield",
            "height.png",
            "zone.kfield",
            "zone.png",
            "surface_height.kfield",
            "surface_height.png",
            "metrics.txt",
            "hypsometry.png",
            "world.json"
        ]
    );
    let metrics = ElevationMetrics::of_world(&world).unwrap();
    assert_eq!(
        read(&dir.join(METRICS_FILE)),
        metrics_txt(&world, &metrics, None).as_bytes()
    );
    let (width, height, pixels) = decode_png(&dir.join(HYPSOMETRY_FILE));
    assert_eq!((width, height), (CHART_WIDTH_PX, CHART_HEIGHT_PX));
    assert_eq!(pixels, hypsometry_chart(&metrics));
    // Land above z = −0.02, (1 + 0.02)/2 of the surface (Archimedes),
    // within the discretisation of 1,000 cells.
    let text = String::from_utf8(read(&dir.join(METRICS_FILE))).unwrap();
    let land: f64 = text
        .lines()
        .find_map(|l| l.strip_prefix("land_fraction = "))
        .unwrap()
        .parse()
        .unwrap();
    assert!((land - 0.51).abs() < 0.01, "{land}");
    remove(&dir);
}

#[test]
fn no_metrics_without_a_surface_height_and_none_left_over() {
    let dir = scratch("no_metrics");
    write_maps_of_width(&dir, &world_with_relief(), 64).unwrap();
    assert!(dir.join(METRICS_FILE).exists() && dir.join(HYPSOMETRY_FILE).exists());
    // A world without a surface height removes the earlier world's metrics.
    let written = write_maps_of_width(&dir, &world(), 64).unwrap();
    assert_eq!(written.len(), 6);
    assert!(!dir.join(METRICS_FILE).exists());
    assert!(!dir.join(HYPSOMETRY_FILE).exists());
    assert!(
        dir.join("surface_height.kfield").exists(),
        "other files stay"
    );
    // So does a rewrite that fails halfway.
    write_maps_of_width(&dir, &world_with_relief(), 64).unwrap();
    std::fs::remove_file(dir.join("zone.kfield")).unwrap();
    std::fs::create_dir(dir.join("zone.kfield")).unwrap();
    assert!(write_maps_of_width(&dir, &world_with_relief(), 64).is_err());
    for file in [METRICS_FILE, HYPSOMETRY_FILE, METADATA_FILE] {
        assert!(!dir.join(file).exists(), "{file}");
    }
    remove(&dir);
}

/// The error of [`write()`] and [`read()`] for a field named `hypsometry`,
/// written out by hand (ticket M2-19).
pub(super) const HYPSOMETRY_RESERVED: &str = concat!(
    "the field name `hypsometry` is reserved: ",
    "hypsometry.png is a file of the world directory itself"
);

/// The names and bytes of the files in `dir`, sorted by name.
fn files_of(dir: &Path) -> Vec<(String, Vec<u8>)> {
    let mut files: Vec<_> = std::fs::read_dir(dir)
        .unwrap()
        .map(|entry| {
            let entry = entry.unwrap();
            (
                entry.file_name().into_string().unwrap(),
                read(&entry.path()),
            )
        })
        .collect();
    files.sort();
    files
}

#[test]
fn rejects_reserved_field_names_before_writing_anything() {
    // The map of a field `hypsometry` would be overwritten by the chart of
    // the metrics.
    let mut world = world_with_relief();
    world.fields[1].name = "hypsometry".into();
    let dir = scratch("reserved_name");
    let err = write(&dir, &world).unwrap_err();
    assert_eq!(err.kind(), io::ErrorKind::InvalidInput);
    assert_eq!(err.to_string(), HYPSOMETRY_RESERVED);
    assert!(!dir.exists(), "nothing is written");

    // Into a directory holding a world: every file stays as it was.
    write_maps_of_width(&dir, &world_with_relief(), 64).unwrap();
    let before = files_of(&dir);
    assert_eq!(before.len(), 10);
    let err = write_maps_of_width(&dir, &world, 64).unwrap_err();
    assert_eq!(err.kind(), io::ErrorKind::InvalidInput);
    assert_eq!(err.to_string(), HYPSOMETRY_RESERVED);
    assert!(files_of(&dir) == before, "a file changed");
    assert!(!dir.join("hypsometry.kfield").exists());
    remove(&dir);
}

#[test]
fn reserved_names_by_hand() {
    // Field names hold no dot, so a fixed file is a field's file only if
    // it ends in `.kfield` or `.png`: of world.json, grid.kgrid,
    // metrics.txt, hypsometry.png and rivers.png, only the last two.
    assert_eq!(fixed_file_of_field("hypsometry"), Some("hypsometry.png"));
    assert_eq!(fixed_file_of_field("rivers"), Some("rivers.png"));
    for name in [
        "river",
        "rivers_map",
        "rivers1",
        "world",
        "grid",
        "metrics",
        "hypsometry_map",
        "hypsometr",
        "hypsometry1",
        "kfield",
        "png",
        "surface_height",
    ] {
        assert_eq!(fixed_file_of_field(name), None, "{name}");
    }
}

#[test]
fn fixed_files_are_every_file_but_the_fields() {
    // With the metrics, every fixed file is written; without, some. A
    // fixed file that only some worlds write is seen only if one of these
    // worlds writes it: a ticket that adds one adds such a world here.
    for (name, world, every) in [
        ("fixed_rivers", world_with_rivers(), true),
        ("fixed_relief", world_with_relief(), false),
        ("fixed_plain", world(), false),
    ] {
        let dir = scratch(name);
        let written = write_maps_of_width(&dir, &world, 64).unwrap();
        let of_fields: Vec<String> = world
            .fields
            .iter()
            .flat_map(|f| [format!("{}.kfield", f.name), format!("{}.png", f.name)])
            .collect();
        let mut rest: Vec<&str> = written
            .iter()
            .map(|p| p.file_name().unwrap().to_str().unwrap())
            .filter(|file| !of_fields.iter().any(|f| f == file))
            .collect();
        rest.sort_unstable();
        let mut fixed = FIXED_FILES.to_vec();
        fixed.sort_unstable();
        if every {
            assert_eq!(rest, fixed);
        } else {
            assert!(rest.iter().all(|file| fixed.contains(file)), "{rest:?}");
        }
        remove(&dir);
    }
}

#[test]
fn fields_may_take_the_stems_of_the_other_fixed_files() {
    // world.json, grid.kgrid and metrics.txt do not collide with
    // world.kfield, grid.png, metrics.png and the others.
    let mut world = world_with_relief();
    world.fields[0].name = "world".into();
    world.fields[1].name = "grid".into();
    let x: Vec<f64> = world.grid.positions().iter().map(|p| p.x).collect();
    world.add_field("metrics", "", x).unwrap();
    let n = world.grid.len();
    world
        .add_class_field("hypsometry_map", (0..n).map(|i| (i % 3) as f64).collect())
        .unwrap();
    let dir = scratch("other_stems");
    let written = write_maps_of_width(&dir, &world, 64).unwrap();
    let names: Vec<_> = written
        .iter()
        .map(|p| p.file_name().unwrap().to_str().unwrap())
        .collect();
    assert_eq!(
        names,
        [
            "grid.kgrid",
            "world.kfield",
            "world.png",
            "grid.kfield",
            "grid.png",
            "surface_height.kfield",
            "surface_height.png",
            "metrics.kfield",
            "metrics.png",
            "hypsometry_map.kfield",
            "hypsometry_map.png",
            "metrics.txt",
            "hypsometry.png",
            "world.json"
        ]
    );
    assert_eq!(read(&dir.join(METADATA_FILE)), world_json(&world).unwrap());
    let grid = grid_file::read_grid(read(&dir.join(GRID_FILE)).as_slice()).unwrap();
    assert_eq!(grid, world.grid);
    let metrics = ElevationMetrics::of_world(&world).unwrap();
    assert_eq!(
        read(&dir.join(METRICS_FILE)),
        metrics_txt(&world, &metrics, None).as_bytes()
    );
    let (width, height, pixels) = decode_png(&dir.join(HYPSOMETRY_FILE));
    assert_eq!((width, height), (CHART_WIDTH_PX, CHART_HEIGHT_PX));
    assert_eq!(pixels, hypsometry_chart(&metrics));
    let bits = |f: &Field| f.values.iter().map(|v| v.to_bits()).collect::<Vec<_>>();
    for field in &world.fields {
        let back =
            field_file::read_field(read(&dir.join(format!("{}.kfield", field.name))).as_slice())
                .unwrap();
        assert_eq!((&back.name, &back.unit), (&field.name, &field.unit));
        assert_eq!(bits(&back), bits(field));
        let (width, height, _) = decode_png(&dir.join(format!("{}.png", field.name)));
        assert_eq!((width, height), (64, 32), "{}", field.name);
    }
    assert_eq!(super::read(&dir).unwrap(), world);
    remove(&dir);
}

#[test]
fn rejects_non_finite_heights_before_writing_anything() {
    for (value, shown) in [(f64::NAN, "NaN"), (f64::NEG_INFINITY, "-inf")] {
        let mut world = world_with_relief();
        world.fields[2].values[3] = value;
        let dir = scratch("non_finite_height");
        let err = write(&dir, &world).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidInput);
        assert_eq!(
            err.to_string(),
            format!("field `surface_height` has the non-finite value {shown} at cell 3")
        );
        assert!(!dir.exists(), "nothing is written");
    }
}

#[test]
fn rejects_elevations_that_overflow_before_writing_anything() {
    // 1e308 m above a sea level of −1e308 m: the elevation overflows.
    let mut world = world_with_relief();
    world.fields[2].values[3] = 1e308;
    world.sea_level_m = -1e308;
    let dir = scratch("overflowing_elevation");
    let err = write(&dir, &world).unwrap_err();
    assert_eq!(err.kind(), io::ErrorKind::InvalidInput);
    assert_eq!(
        err.to_string(),
        "field `surface_height` has the value 1e308 at cell 3, whose elevation above \
         the sea level of -1e308 m overflows"
    );
    assert!(!dir.exists(), "nothing is written");
}

/// [`world_with_relief`] with every height at `base_m` except cell 3 at
/// `peak_m`, and the sea level at `sea_level_m`.
fn world_with_one_peak(base_m: f64, peak_m: f64, sea_level_m: f64) -> World {
    let mut world = world_with_relief();
    world.fields[2].values.fill(base_m);
    world.fields[2].values[3] = peak_m;
    world.sea_level_m = sea_level_m;
    world
}

#[test]
fn rejects_elevations_beyond_the_histogram_before_writing_anything() {
    // Each case: the base, the peak and the sea level (m), and the number
    // of bins of 100 m by hand.
    let cases = [
        // Bins 0 and 10⁶: 1,000,001, one more than the bound (the old code
        // wrote it).
        (0.0, 1e8, 0.0, "1.000001e6"),
        // Bins 0 and 1e10 (the old code panicked, more than 2³¹ bins).
        (0.0, 1e12, 0.0, "1.0000000001e10"),
        // Counted above the sea level: the elevations 60 m and
        // 100,000,010 m are in bins 0 and 10⁶, while the heights alone
        // (0 m and 99,999,950 m, bins 0 and 999,999) would fit.
        (0.0, 99_999_950.0, -60.0, "1.000001e6"),
    ];
    for (base, peak, sea_level, bins) in cases {
        let world = world_with_one_peak(base, peak, sea_level);
        let dir = scratch("beyond_the_histogram");
        let err = write(&dir, &world).unwrap_err();
        assert_eq!(
            err.kind(),
            io::ErrorKind::InvalidInput,
            "{peak} {sea_level}"
        );
        assert_eq!(
            err.to_string(),
            format!(
                "the elevations of field `surface_height` need {bins} bins of 100 m for the \
                 histogram of metrics.txt, more than 1e6"
            )
        );
        assert!(!dir.exists(), "nothing is written");
    }
}

#[test]
fn accepts_elevations_at_the_bound() {
    // Each case needs exactly 10⁶ bins of 100 m. Checked through
    // `elevation_metrics`, which `write` calls first, to keep the test
    // from writing a `metrics.txt` of 29 MB.
    let cases = [
        // Bins 0 and 999,999.
        (0.0, 99_999_900.0, 0.0),
        // The elevations −60 m and 99,999,890 m: bins −1 and 999,998.
        (0.0, 99_999_950.0, 60.0),
        // The elevations 0 m and 99,999,950 m: bins 0 and 999,999, while
        // the heights alone (50 m and 1e8 m, bins 0 and 10⁶) would not fit.
        (50.0, 1e8, 50.0),
    ];
    for (base, peak, sea_level) in cases {
        let world = world_with_one_peak(base, peak, sea_level);
        let metrics = elevation_metrics(&world).unwrap().unwrap();
        let shares = metrics.hypsometry.shares();
        assert_eq!(shares.len(), 1_000_000, "{peak} {sea_level}");
        assert!(shares[0] > 0.99 && shares[999_999] > 0.0);
        assert_eq!(metrics.highest_m.to_bits(), (peak - sea_level).to_bits());
    }
}

#[test]
fn huge_elevations_have_finite_metrics() {
    // Every cell 1e308 m above the sea level: the area-weighted sum of the
    // elevations overflows (4π sr · 1e308 m), their mean does not.
    let mut world = world_with_relief();
    world.fields[2].values.fill(1e308);
    world.sea_level_m = 0.0;
    let dir = scratch("huge_elevations");
    write_maps_of_width(&dir, &world, 64).unwrap();
    let text = String::from_utf8(read(&dir.join(METRICS_FILE))).unwrap();
    assert!(!text.contains("inf") && !text.contains("NaN"), "{text}");
    let mean = text
        .lines()
        .find_map(|line| line.strip_prefix("mean_land_elevation_m = "))
        .unwrap();
    let mean: f64 = mean.parse().unwrap();
    assert!((mean / 1e308 - 1.0).abs() < 1e-15, "{mean}");
    remove(&dir);
}

/// A fresh path for the test `name` under the system's temporary directory
/// (unique to this process); the directory itself does not exist.
fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("kuniumi-world_dir-{}-{name}", std::process::id()));
    remove(&dir);
    dir
}

fn remove(dir: &Path) {
    if dir.exists() {
        std::fs::remove_dir_all(dir).unwrap();
    }
}

fn read(path: &Path) -> Vec<u8> {
    std::fs::read(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

/// Decodes an 8-bit RGB PNG: width, height and pixels.
fn decode_png(path: &Path) -> (u32, u32, Vec<Rgb8>) {
    let decoder = ::png::Decoder::new(std::io::Cursor::new(read(path)));
    let mut reader = decoder.read_info().unwrap();
    let mut buf = vec![0; reader.output_buffer_size().unwrap()];
    let info = reader.next_frame(&mut buf).unwrap();
    assert_eq!(info.color_type, ::png::ColorType::Rgb);
    assert_eq!(info.bit_depth, ::png::BitDepth::Eight);
    buf.truncate(info.buffer_size());
    let pixels = buf.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
    (info.width, info.height, pixels)
}

// ---------------------------------------------------------------------------
// The river map (ticket M3-05).

/// The error of [`write()`] and [`read()`] for a field named `rivers`,
/// written out by hand.
pub(super) const RIVERS_RESERVED: &str = concat!(
    "the field name `rivers` is reserved: ",
    "rivers.png is a file of the world directory itself"
);

/// [`octahedron_world`] after an Erosion stage: the land cells 0, 2, 4 and
/// 5 (+700, +1,500, +3,000 and +700 m above the sea level) carry 1, 1,000,
/// 0.1 and 100 m³/s, the ocean cells 1 and 3 carry 5,000 and 7 m³/s.
fn octahedron_world_with_rivers() -> World {
    let mut world = octahedron_world();
    world.stages.push("erosion".into());
    world
        .add_field(
            "discharge",
            "m^3/s",
            vec![1.0, 5_000.0, 1_000.0, 7.0, 0.1, 100.0],
        )
        .unwrap();
    world
}

/// [`world_with_relief`] with a discharge that grows from 1 m³/s at the
/// north pole to 10¹² m³/s at the south pole, `10^(6 (1 − z))`: on the
/// land, `z > −0.02`, up to about 10^6.12 m³/s at the coast.
pub(super) fn world_with_rivers() -> World {
    let mut world = world_with_relief();
    let discharge = world
        .grid
        .positions()
        .iter()
        .map(|p| 10.0_f64.powf(6.0 * (1.0 - p.z)))
        .collect();
    world.add_field("discharge", "m^3/s", discharge).unwrap();
    world
}

/// The colour `base` with the river colour mixed in with the weight `w`,
/// rounded half away from zero.
fn mixed(base: Rgb8, w: f64) -> Rgb8 {
    core::array::from_fn(|k| {
        (f64::from(base[k]) * (1.0 - w) + f64::from(RIVER_COLOUR[k]) * w).round() as u8
    })
}

/// The colours of the map of `surface_height` in `world`.
fn elevation_colours(world: &World) -> Vec<Rgb8> {
    let metrics = ElevationMetrics::of_world(world).unwrap();
    let height = world.field("surface_height").unwrap();
    let (colormap, range) = map_colours(height, world, Some(&metrics));
    colormap.colours(&height.values, range)
}

#[test]
fn river_colours_by_hand() {
    let world = octahedron_world_with_rivers();
    let metrics = ElevationMetrics::of_world(&world).unwrap();
    let base = elevation_colours(&world);
    // The largest discharge on land is 1,000 m³/s (cell 2; the ocean's
    // 5,000 m³/s does not count), the smallest 0.1 m³/s (cell 4), so rivers
    // are drawn down to max(1,000 / 1,000, 5 · 0.1) = 1 m³/s, with the
    // weight 0.35 + 0.65 log₁₀₀₀ Q: cell 0 0.35, cell 2 exactly 1 (the
    // river colour), cell 5 0.35 + 0.65 · 2/3; cell 4 and the ocean keep
    // their colours.
    let expected = vec![
        mixed(base[0], 0.35),
        base[1],
        RIVER_COLOUR,
        base[3],
        base[4],
        mixed(base[5], 0.35 + 0.65 * 2.0 / 3.0),
    ];
    assert_eq!(river_colours(&world, Some(&metrics)), Some(expected));
    // Cell 0 by hand: +700 m on the land scale of 0 to +9,000 m (the
    // chart's axis), t = 7/90, 4t = 0.3111 of the way from [46, 120, 56] to
    // [138, 172, 82]: [74.62, 136.18, 64.09] → [75, 136, 64]; then 0.65 of
    // it and 0.35 of [10, 40, 190]: [52.25, 102.4, 108.1] → [52, 102, 108].
    assert_eq!(base[0], [75, 136, 64]);
    assert_eq!(
        river_colours(&world, Some(&metrics)).unwrap()[0],
        [52, 102, 108]
    );
    // Without the metrics the elevation map spans the heights' own range.
    let auto = {
        let height = world.field("surface_height").unwrap();
        let (colormap, range) = map_colours(height, &world, None);
        colormap.colours(&height.values, range)
    };
    assert_eq!(river_colours(&world, None).unwrap()[1], auto[1]);
    assert_eq!(river_colours(&world, None).unwrap()[2], RIVER_COLOUR);
}

#[test]
fn river_colours_without_rivers() {
    // No discharge, or no surface height: no river map.
    assert_eq!(river_colours(&octahedron_world(), None), None);
    let mut world = octahedron_world_with_rivers();
    world.fields.retain(|f| f.name != "surface_height");
    assert_eq!(river_colours(&world, None), None);

    // A cell exactly at the sea level is ocean to the drainage graph: its
    // 5,000 m³/s are neither drawn nor counted, as in `river_colours_by_hand`.
    let mut world = octahedron_world_with_rivers();
    world.fields[0].values[1] = world.sea_level_m;
    let metrics = ElevationMetrics::of_world(&world).unwrap();
    let base = elevation_colours(&world);
    let expected = vec![
        mixed(base[0], 0.35),
        base[1],
        RIVER_COLOUR,
        base[3],
        base[4],
        mixed(base[5], 0.35 + 0.65 * 2.0 / 3.0),
    ];
    assert_eq!(river_colours(&world, Some(&metrics)), Some(expected));

    // A coarse world: the largest discharge on land is 10 m³/s and the
    // smallest 0.5 m³/s, so rivers are drawn down to max(10 / 1,000,
    // 5 · 0.5) = 2.5 m³/s, not to 0.01 m³/s, with the weight
    // 0.35 + 0.65 ln(Q / 2.5) / ln 4: cell 0 (3 m³/s) 0.35 + 0.65 · 0.1315,
    // cell 2 the river colour, cell 5 (2 m³/s) not drawn.
    let mut world = octahedron_world_with_rivers();
    world.fields[1].values = vec![3.0, 5_000.0, 10.0, 7.0, 0.5, 2.0];
    let metrics = ElevationMetrics::of_world(&world).unwrap();
    let base = elevation_colours(&world);
    let weight = |q: f64| 0.35 + 0.65 * (q / 2.5).ln() / 4.0_f64.ln();
    assert!((weight(3.0) - (0.35 + 0.65 * 0.1315)).abs() < 1e-4);
    let expected = vec![
        mixed(base[0], weight(3.0)),
        base[1],
        RIVER_COLOUR,
        base[3],
        base[4],
        base[5],
    ];
    assert_eq!(river_colours(&world, Some(&metrics)), Some(expected));

    // Discharges that are not finite are not drawn and do not count: the
    // largest is then 100 m³/s (cell 5) and the smallest 1 m³/s (cell 0),
    // so rivers are drawn down to max(0.1, 5) = 5 m³/s: cell 5 alone, in
    // the river colour.
    let mut world = octahedron_world_with_rivers();
    world.fields[1].values[2] = f64::INFINITY;
    world.fields[1].values[4] = f64::NAN;
    let metrics = ElevationMetrics::of_world(&world).unwrap();
    let base = elevation_colours(&world);
    let expected = vec![base[0], base[1], base[2], base[3], base[4], RIVER_COLOUR];
    assert_eq!(river_colours(&world, Some(&metrics)), Some(expected));

    // The threshold at the largest discharge, 5 · 2 = 10 m³/s: cell 2 alone
    // is drawn, in the river colour (t = 1 where Q_min = Q_max).
    let mut world = octahedron_world_with_rivers();
    world.fields[1].values = vec![2.0, 5_000.0, 10.0, 7.0, 4.0, 3.0];
    let metrics = ElevationMetrics::of_world(&world).unwrap();
    let base = elevation_colours(&world);
    let expected = vec![base[0], base[1], RIVER_COLOUR, base[3], base[4], base[5]];
    assert_eq!(river_colours(&world, Some(&metrics)), Some(expected));

    // No land, no positive discharge on land, or a single cell with a
    // discharge (Q_min = 5 · 10 m³/s, above it): the elevation map alone.
    for (heights, discharge) in [
        (vec![-600.0; 6], vec![1.0; 6]),
        (vec![200.0; 6], vec![0.0; 6]),
        (vec![200.0; 6], vec![-1.0, 0.0, -0.0, -5.0, 0.0, -1e300]),
        (vec![200.0; 6], vec![0.0, 0.0, 10.0, 0.0, -1.0, f64::NAN]),
    ] {
        let mut world = octahedron_world_with_rivers();
        world.fields[0].values = heights;
        world.fields[1].values = discharge;
        let metrics = ElevationMetrics::of_world(&world).unwrap();
        assert_eq!(
            river_colours(&world, Some(&metrics)),
            Some(elevation_colours(&world))
        );
    }
}

#[test]
#[should_panic(
    expected = "river_colours needs one discharge per height: 5 discharges for 6 heights"
)]
fn river_colours_need_one_discharge_per_height() {
    // `write()` refuses such a world first (`check_fields`).
    let mut world = octahedron_world_with_rivers();
    world.fields[1].values.pop();
    let _ = river_colours(&world, None);
}

#[test]
fn writes_the_river_map_after_the_metrics() {
    let dir = scratch("rivers");
    let world = world_with_rivers();
    let written = write_maps_of_width(&dir, &world, 64).unwrap();
    let names: Vec<_> = written
        .iter()
        .map(|p| p.file_name().unwrap().to_str().unwrap())
        .collect();
    assert_eq!(
        names,
        [
            "grid.kgrid",
            "height.kfield",
            "height.png",
            "zone.kfield",
            "zone.png",
            "surface_height.kfield",
            "surface_height.png",
            "discharge.kfield",
            "discharge.png",
            "metrics.txt",
            "hypsometry.png",
            "rivers.png",
            "world.json"
        ]
    );
    let metrics = ElevationMetrics::of_world(&world).unwrap();
    let colours = river_colours(&world, Some(&metrics)).unwrap();
    let cells = render::equirect_cells(&world.grid, 64);
    let (width, height, pixels) = decode_png(&dir.join(RIVERS_FILE));
    assert_eq!((width, height), (64, 32));
    let expected: Vec<Rgb8> = cells.iter().map(|&c| colours[c as usize]).collect();
    assert_eq!(pixels, expected);
    // The rivers are drawn over the elevation map where the discharge is
    // at least 10⁻³ of the largest on land (at the coast, z ≈ −0.02), so
    // on land south of z ≈ 1 − 3.12/6 = 0.48: there the colour differs from
    // that of surface_height.png; the ocean and the land north of it look
    // the same.
    let (_, _, map) = decode_png(&dir.join("surface_height.png"));
    let positions = world.grid.positions();
    let heights = &world.field("surface_height").unwrap().values;
    let mut drawn = 0;
    for (k, &c) in cells.iter().enumerate() {
        let (z, land) = (positions[c as usize].z, heights[c as usize] > -700.0);
        if land && z < 0.45 {
            assert_ne!(pixels[k], map[k], "pixel {k}");
            drawn += 1;
        } else if !land || z > 0.51 {
            assert_eq!(pixels[k], map[k], "pixel {k}");
        }
    }
    assert!(drawn > 100, "{drawn}");

    // A world without a discharge removes the earlier world's river map,
    // and so does a rewrite that fails halfway.
    let written = write_maps_of_width(&dir, &world_with_relief(), 64).unwrap();
    assert!(!written.iter().any(|p| p.ends_with(RIVERS_FILE)));
    assert!(!dir.join(RIVERS_FILE).exists());
    write_maps_of_width(&dir, &world, 64).unwrap();
    std::fs::remove_file(dir.join("zone.kfield")).unwrap();
    std::fs::create_dir(dir.join("zone.kfield")).unwrap();
    assert!(write_maps_of_width(&dir, &world, 64).is_err());
    assert!(!dir.join(RIVERS_FILE).exists());
    remove(&dir);
}

#[test]
fn rejects_the_name_rivers_before_writing_anything() {
    // The map of a field `rivers` would be overwritten by the river map,
    // in any world.
    for mut world in [world(), world_with_rivers()] {
        world.fields[1].name = "rivers".into();
        let dir = scratch("rivers_reserved");
        let err = write(&dir, &world).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidInput);
        assert_eq!(err.to_string(), RIVERS_RESERVED);
        assert!(!dir.exists(), "nothing is written");
    }
}

// ---------------------------------------------------------------------------
// The drainage metrics (ticket M3-08).

/// [`octahedron_world`] after the Erosion stage: with a `drainage_area`
/// (whose values the metrics do not read).
fn octahedron_world_with_drainage() -> World {
    let mut world = octahedron_world();
    world.stages.push("erosion".into());
    world
        .add_field("drainage_area", "m^2", vec![1.0; 6])
        .unwrap();
    world
}

/// The drainage section of `metrics.txt` for
/// [`octahedron_world_with_drainage`], written by hand: on Earth's radius
/// the four land cells (+X, +Y, +Z and −Z) each drain straight into the
/// ocean (−X or −Y, lower than every other neighbour), so each drains its
/// own area, 4πR²/6 = 8.501 074 5 × 10¹³ m² (Python, `4 * math.pi *
/// 6371000.0**2 / 6`), and no main stream is longer than 0: no basin. The
/// threshold is 100 such areas; the exceedance runs from 5e13 m², the
/// largest value of the ladder not above the land's drainage areas, to
/// 1e14 m², the smallest above.
fn octahedron_drainage_lines() -> &'static str {
    "
# Drainage (ticket M3-08): the rivers of the final surface, routed across
# its closed depressions. Land is above the sea level here, as in the
# drainage graph. Hack's law L = c A^h: least squares of ln L on ln A over
# the land cells whose drainage area A (m²) is at least hack_min_area_m2
# (hack_min_cells mean cell areas) and whose main-stream length L (m), the
# longest flow path up to a source, is positive; c is in m^(1-2h). Earth:
# h about 0.6 (Hack 1957).
hack_min_cells = 100
hack_min_area_m2 = 8.5011e15
hack_basins = 0
hack_exponent = none
hack_coefficient = none
largest_drainage_area_m2 = 8.5011e13

# Exceedance: the share of the land whose drainage area is at least the
# area in the key (m²), for 1, 2 and 5 times the powers of ten from the
# smallest drainage area on land to the largest. Earth: the share falls
# as a^-0.43 (Rodríguez-Iturbe & Rinaldo 1997).
exceedance_5e13_m2 = 1.000000
exceedance_1e14_m2 = 0.000000
"
}

/// `text` with `section` inserted before the blank line that precedes the
/// histogram's comment.
fn with_section_before_histogram(text: &str, section: &str) -> String {
    let at = text
        .find("\n\n# Histogram")
        .expect("the histogram's comment");
    format!("{}{section}{}", &text[..at + 1], &text[at + 1..])
}

#[test]
fn metrics_txt_with_drainage_by_hand() {
    let world = octahedron_world_with_drainage();
    let metrics = ElevationMetrics::of_world(&world).unwrap();
    let drainage = DrainageMetrics::of_world(&world).unwrap();
    let text = metrics_txt(&world, &metrics, Some(&drainage));
    let expected =
        with_section_before_histogram(&octahedron_metrics_txt(), octahedron_drainage_lines());
    for (k, (a, b)) in text.lines().zip(expected.lines()).enumerate() {
        assert_eq!(a, b, "line {}", k + 1);
    }
    assert_eq!(text, expected);
    // The histogram stays last.
    assert!(text.ends_with("3000 3100 0.16666667\n"), "{text}");
    // Without drainage metrics the file is that of ticket M2-06.
    assert_eq!(
        metrics_txt(&world, &metrics, None),
        octahedron_metrics_txt()
    );
}

#[test]
fn drainage_lines_by_hand() {
    // A fit, a ladder below and above 1 m², shares in 6 decimals.
    let world = octahedron_world_with_drainage();
    let metrics = ElevationMetrics::of_world(&world).unwrap();
    let ladder = |m: u32, k: i32, a: f64, share: f64| Exceedance {
        mantissa: m,
        power_of_ten: k,
        drainage_area_m2: a,
        share,
    };
    let drainage = DrainageMetrics {
        hack_min_cells: 12.5,
        hack_min_area_m2: 9.812_345e10,
        hack_basins: 9_339,
        hack: Some(HackFit {
            exponent: 0.524_702_864_681_634_2,
            coefficient: 0.907_098_300_575_099,
        }),
        largest_drainage_area_m2: Some(4.933_9e12),
        exceedance: vec![
            ladder(5, -1, 0.5, 1.0),
            ladder(1, 0, 1.0, 0.533_562_4),
            ladder(2, 0, 2.0, 0.000_000_4),
            ladder(5, 0, 5.0, 0.0),
        ],
    };
    let text = metrics_txt(&world, &metrics, Some(&drainage));
    for line in [
        "hack_min_cells = 12.5",
        "hack_min_area_m2 = 9.8123e10",
        "hack_basins = 9339",
        "hack_exponent = 0.5247",
        "hack_coefficient = 9.0710e-1",
        "largest_drainage_area_m2 = 4.9339e12",
        "exceedance_5e-1_m2 = 1.000000",
        "exceedance_1e0_m2 = 0.533562",
        "exceedance_2e0_m2 = 0.000000",
        "exceedance_5e0_m2 = 0.000000",
    ] {
        assert!(text.lines().any(|l| l == line), "`{line}` missing:\n{text}");
    }
    // Without land: no fit, no largest area, no exceedance lines.
    let drainage = DrainageMetrics {
        hack_basins: 0,
        hack: None,
        largest_drainage_area_m2: None,
        exceedance: Vec::new(),
        ..drainage
    };
    let text = metrics_txt(&world, &metrics, Some(&drainage));
    for line in [
        "hack_basins = 0",
        "hack_exponent = none",
        "hack_coefficient = none",
        "largest_drainage_area_m2 = none",
    ] {
        assert!(text.lines().any(|l| l == line), "`{line}` missing:\n{text}");
    }
    assert!(!text.contains("exceedance_"), "{text}");
    assert!(text.contains("(Rodríguez-Iturbe & Rinaldo 1997).\n\n# Histogram"));
}

/// [`world_with_relief`] after the Erosion stage: with a `drainage_area`.
fn world_with_drainage() -> World {
    let mut world = world_with_relief();
    world.stages.push("erosion".into());
    let n = world.grid.len();
    world
        .add_field("drainage_area", "m^2", vec![1.0; n])
        .unwrap();
    world
}

#[test]
fn writes_drainage_metrics_after_erosion() {
    let dir = scratch("drainage");
    let world = world_with_drainage();
    let written = write_maps_of_width(&dir, &world, 64).unwrap();
    assert!(written.iter().any(|p| p.ends_with(METRICS_FILE)));
    let metrics = ElevationMetrics::of_world(&world).unwrap();
    let drainage = DrainageMetrics::of_world(&world).unwrap();
    assert_eq!(
        read(&dir.join(METRICS_FILE)),
        metrics_txt(&world, &metrics, Some(&drainage)).as_bytes()
    );
    // The cap of 1,000 cells south of z ≈ −0.02 is land, all draining
    // towards the coast: rivers, and a section with lines.
    let text = String::from_utf8(read(&dir.join(METRICS_FILE))).unwrap();
    assert!(text.contains("\nhack_min_cells = 100\n"), "{text}");
    assert!(text.contains("\nexceedance_"), "{text}");
    assert!(drainage.largest_drainage_area_m2.is_some());
    // A world without `drainage_area` writes no section.
    write_maps_of_width(&dir, &world_with_relief(), 64).unwrap();
    let text = String::from_utf8(read(&dir.join(METRICS_FILE))).unwrap();
    assert!(!text.contains("hack_"), "{text}");
    remove(&dir);
}

#[test]
fn rejects_a_radius_beyond_the_drainage_limits_before_writing_anything() {
    for radius_m in [0.0, -1.0, 1e160, 1e-160] {
        let mut world = world_with_drainage();
        world.params.radius_m = radius_m;
        let dir = scratch("drainage_radius");
        let err = write(&dir, &world).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidInput);
        assert_eq!(
            err.to_string(),
            format!(
                "the planet's radius of {radius_m:e} m is beyond the limits of the drainage \
                 metrics of metrics.txt (DrainageMetrics::radius_in_limits)"
            )
        );
        assert!(!dir.exists(), "nothing is written");
        // Without a surface height there are no metrics, and no limit.
        world.fields.retain(|f| f.name != "surface_height");
        write_maps_of_width(&dir, &world, 64).unwrap();
        assert!(!dir.join(METRICS_FILE).exists());
        remove(&dir);
    }
}

#[test]
fn rejects_heights_beyond_the_drainage_limits_before_writing_anything() {
    // Heights above 1e300 m pass the checks of the elevation metrics (one
    // bin of the histogram when all are equal) but would overflow the fill
    // of the depressions that the drainage metrics rebuild (review round 1,
    // R1-1).
    for h in [f64::MAX, 1.5e300] {
        let mut world = octahedron_world_with_drainage();
        world.fields[0].values = vec![h; 6];
        let dir = scratch("drainage_heights");
        let err = write(&dir, &world).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidInput);
        assert_eq!(
            err.to_string(),
            format!(
                "field `surface_height` has the value {h:e} at cell 0, above the 1e300 m that \
                 the drainage metrics of metrics.txt take (DrainageMetrics::heights_in_limits)"
            )
        );
        assert!(!dir.exists(), "nothing is written");
        // Without `drainage_area` the elevation metrics are written.
        world.fields.retain(|f| f.name != "drainage_area");
        write_maps_of_width(&dir, &world, 64).unwrap();
        assert!(dir.join(METRICS_FILE).exists());
        remove(&dir);
    }
    // At 1e300 m the metrics are written.
    let mut world = octahedron_world_with_drainage();
    world.fields[0].values = vec![1e300; 6];
    let dir = scratch("drainage_heights_at_the_limit");
    write_maps_of_width(&dir, &world, 64).unwrap();
    let text = String::from_utf8(read(&dir.join(METRICS_FILE))).unwrap();
    assert!(text.contains("\nhack_basins = 0\n"), "{text}");
    remove(&dir);
}

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
//! - `maps_use_the_categorical_palette_for_classes`;
//! - `writes_maps_of_2048_by_1024`.

use std::path::PathBuf;

use kuniumi_math::rng::Rng;
use kuniumi_sphere::SphereGrid;
use kuniumi_world::{
    Field, FieldKind, PlanetParams, RotationRegime, Settings, TectonicRegime, World,
};

use super::*;
use crate::colormap::{Colormap, Range};

/// A world on a 1,000-cell grid with a quantity field and a class field.
fn world() -> World {
    let grid = SphereGrid::fibonacci(1_000, 1).expect("a valid grid");
    let n = grid.len();
    let mut world = World {
        params: PlanetParams::earth(),
        settings: Settings {
            cells: 1_000,
            allow_partial: true,
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
/// −0, values that need 17 digits, the v2 regimes.
fn awkward_world() -> World {
    World {
        params: PlanetParams {
            radius_m: f64::MIN_POSITIVE,
            surface_gravity_m_s2: 5e-324,
            rotation_period_s: f64::MAX,
            obliquity_rad: 0.1,
            orbital_period_s: 1e23,
            orbital_eccentricity: 1.0 / 3.0,
            stellar_flux_w_m2: 2.225_073_858_507_201e-308,
            ocean_volume_m3: 9_007_199_254_740_991.0,
            surface_pressure_pa: -0.0,
            continental_crust_fraction: 0.1 + 0.2,
            rotation: RotationRegime::TidallyLocked,
            tectonics: TectonicRegime::StagnantLid,
        },
        settings: Settings {
            cells: u32::MAX,
            allow_partial: false,
        },
        seed: u64::MAX,
        sea_level_m: -1.5e-320,
        stages: Vec::new(),
        ..world()
    }
}

/// The `f64` metadata of `world`, as bits, listed by hand from the `World`.
fn world_floats(world: &World) -> Vec<u64> {
    let p = &world.params;
    [
        p.radius_m,
        p.surface_gravity_m_s2,
        p.rotation_period_s,
        p.obliquity_rad,
        p.orbital_period_s,
        p.orbital_eccentricity,
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
    let p = &json.planet;
    [
        p.radius_m,
        p.surface_gravity_m_s2,
        p.rotation_period_s,
        p.obliquity_rad,
        p.orbital_period_s,
        p.orbital_eccentricity,
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
        assert_eq!(back.format_version, 1);
        assert_eq!(back.generator_version, kuniumi_world::GENERATOR_VERSION);
        assert_eq!(back.seed, world.seed);
        assert_eq!(
            (back.settings.cells, back.settings.allow_partial),
            (world.settings.cells, world.settings.allow_partial)
        );
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
/// literals of `PlanetParams::earth()`, and its obliquity
/// `23.44 · (π/180)` as `repr(23.44 * (math.pi / 180))` prints it in
/// Python 3 (0.40910517666747087). It is the example in the module
/// documentation.
const EARTH_JSON: &str = r#"{
  "format_version": 1,
  "generator_version": "0.0.1",
  "seed": 1,
  "settings": {
    "cells": 20000,
    "allow_partial": true
  },
  "planet": {
    "radius_m": 6371000.0,
    "surface_gravity_m_s2": 9.80665,
    "rotation_period_s": 86164.0905,
    "obliquity_rad": 0.40910517666747087,
    "orbital_period_s": 31558149.8,
    "orbital_eccentricity": 0.0167,
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
            allow_partial: true,
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
  "format_version": 1,
  "generator_version": "0.0.1",
  "seed": 18446744073709551615,
  "settings": {
    "cells": 4294967295,
    "allow_partial": false
  },
  "planet": {
    "radius_m": 2.2250738585072014e-308,
    "surface_gravity_m_s2": 5e-324,
    "rotation_period_s": 1.7976931348623157e+308,
    "obliquity_rad": 0.1,
    "orbital_period_s": 1e+23,
    "orbital_eccentricity": 0.3333333333333333,
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
    assert_eq!(written_a.len(), 5);
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
    let cases: [(Spoil, &str); 4] = [
        (|w| w.sea_level_m = f64::NAN, "`sea_level_m`"),
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
}

#[test]
fn rejects_bad_fields_before_writing_anything() {
    let cases: [(Spoil, &str); 6] = [
        (
            |w| w.fields[1].name = "Zone".into(),
            "`Zone` is not a valid field name",
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
            "height.kfield",
            "height.png",
            "zone.kfield",
            "zone.png",
            "world.json"
        ]
    );
    assert!(written.iter().all(|p| p.parent() == Some(dir.as_path())));
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
    let world = world();
    assert_eq!(colormap(&world.fields[0]), Colormap::Viridis);
    assert_eq!(colormap(&world.fields[1]), Colormap::Categorical);
    let dir = scratch("palettes");
    write_maps_of_width(&dir, &world, 64).unwrap();
    for (field, expected) in [
        (&world.fields[0], Colormap::Viridis),
        (&world.fields[1], Colormap::Categorical),
    ] {
        let (width, height, pixels) = decode_png(&dir.join(format!("{}.png", field.name)));
        assert_eq!((width, height), (64, 32));
        let reference = render::equirect(&world.grid, &field.values, expected, Range::Auto, 64);
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

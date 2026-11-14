//! Tests of reading world directories back (ticket MV-01):
//!
//! - `reads_back_a_generated_world_bit_exactly` (20,000 cells, from
//!   `kuniumi_sim::generate`) and `reads_back_worlds_of_every_shape`;
//! - `rejects_other_format_versions` and `rejects_malformed_world_json`
//!   (`world.json`);
//! - `rejects_a_listed_reserved_name_even_with_its_file` (ticket M2-19);
//! - `rejects_bad_kfield_files` and `rejects_bad_grid_files`;
//! - `a_missing_directory_is_an_io_error` and `errors_name_their_file`.

use std::path::{Path, PathBuf};

use kuniumi_world::{Field, PlanetParams, Settings, World};

use super::tests::{
    HYPSOMETRY_RESERVED, awkward_world, octahedron_world, world, world_with_relief,
};
use super::*;

/// A fresh path for the test `name` under the system's temporary directory
/// (unique to this process); the directory itself does not exist.
fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "kuniumi-world_dir-read-{}-{name}",
        std::process::id()
    ));
    remove(&dir);
    dir
}

fn remove(dir: &Path) {
    if dir.exists() {
        std::fs::remove_dir_all(dir).unwrap();
    }
}

fn bytes_of(path: &Path) -> Vec<u8> {
    std::fs::read(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

/// The bits of every number of `world` and its other contents, listed by
/// hand from the public members and accessors: the metadata, every array
/// of the grid and every field, so that two worlds with the same list are
/// the same bit for bit (`==` on `f64` cannot tell −0 from 0).
fn contents(world: &World) -> Vec<String> {
    let p = &world.params;
    let mut out = vec![
        format!(
            "params {:?} {:?}",
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
            ]
            .map(f64::to_bits),
            (p.rotation, p.tectonics)
        ),
        format!("settings {:?}", world.settings),
        format!("seed {}", world.seed),
        format!("sea level {:#x}", world.sea_level_m.to_bits()),
        format!("stages {:?}", world.stages),
    ];
    let grid = &world.grid;
    let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
    for i in 0..grid.len() as u32 {
        let c = grid.position(i);
        out.push(format!(
            "cell {i}: {:?} {:#x} {:?} {:?} {:?}",
            [c.x, c.y, c.z].map(f64::to_bits),
            grid.area_sr(i).to_bits(),
            grid.neighbors(i),
            bits(grid.edge_lengths_rad(i)),
            bits(grid.centre_distances_rad(i)),
        ));
    }
    for f in &world.fields {
        out.push(format!(
            "field {} {} {:?} {:?}",
            f.name,
            f.unit,
            f.kind,
            bits(&f.values)
        ));
    }
    out
}

/// Writes `world` into `dir` with small maps, reads it back and checks that
/// the two are the same, bit for bit; then writes the world read back into
/// `copy` and checks that every file is the same.
fn assert_round_trip(world: &World, dir: &Path, copy: &Path) {
    let written = write_maps_of_width(dir, world, 64).unwrap();
    let back = read(dir).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(&back, world);
    assert!(contents(&back) == contents(world), "not the same bits");
    let again = write_maps_of_width(copy, &back, 64).unwrap();
    assert_eq!(written.len(), again.len());
    for (a, b) in written.iter().zip(&again) {
        assert_eq!(a.file_name(), b.file_name());
        assert!(bytes_of(a) == bytes_of(b), "{} differs", b.display());
    }
}

#[test]
fn reads_back_a_generated_world_bit_exactly() {
    let settings = Settings {
        cells: 20_000,
        allow_partial: true,
        ..Settings::default()
    };
    let world = kuniumi_sim::generate(&PlanetParams::earth(), &settings, 1, &mut |_| {}).unwrap();
    // Every stage that exists so far, with its fields.
    assert_eq!(world.stages, ["grid", "tectonics", "crust"]);
    assert_eq!(world.grid.len(), 20_000);
    assert_eq!(world.fields.len(), 13);
    assert!(world.fields.iter().any(|f| f.kind == FieldKind::Class));
    let (dir, copy) = (scratch("generated"), scratch("generated_copy"));
    assert_round_trip(&world, &dir, &copy);
    remove(&dir);
    remove(&copy);
}

#[test]
fn reads_back_worlds_of_every_shape() {
    let mut no_fields = world();
    no_fields.fields.clear();
    no_fields.stages.clear();
    let mut negative_zero = world();
    negative_zero.sea_level_m = -0.0;
    negative_zero.fields[0].values[3] = -0.0;
    negative_zero.fields[0].values[4] = f64::MIN_POSITIVE / 4.0;
    let cases = [
        ("small", world()),
        // Metadata at the limits of f64, the v2 regimes, settings that are
        // not valid: read as written, not validated.
        ("awkward", awkward_world()),
        // A grid that is not the Fibonacci grid of the settings and seed.
        ("octahedron", octahedron_world()),
        ("no_fields", no_fields),
        ("negative_zero", negative_zero),
        ("relief", world_with_relief()),
    ];
    for (name, world) in cases {
        let (dir, copy) = (scratch(name), scratch(&format!("{name}_copy")));
        assert_round_trip(&world, &dir, &copy);
        remove(&dir);
        remove(&copy);
    }
}

/// Writes [`world`] into a fresh directory for the test `name`, applies
/// `spoil` to it and reads it back.
fn read_spoiled(name: &str, spoil: impl FnOnce(&Path)) -> (PathBuf, ReadError) {
    let dir = scratch(name);
    write_maps_of_width(&dir, &world(), 64).unwrap();
    spoil(&dir);
    let err = read(&dir).expect_err("the spoiled directory is rejected");
    (dir, err)
}

/// `world.json` of [`world`] with `edit` applied to its text.
fn edit_world_json(dir: &Path, edit: impl FnOnce(String) -> String) {
    let path = dir.join(METADATA_FILE);
    let text = String::from_utf8(bytes_of(&path)).unwrap();
    std::fs::write(&path, edit(text)).unwrap();
}

/// Checks that `err`, of the case `what`, is about the file `file` of
/// `dir`, is the variant `kind` and has the message `message` after the
/// path.
fn assert_error(what: &str, err: &ReadError, dir: &Path, file: &str, kind: &str, message: &str) {
    let path = dir.join(file);
    assert_eq!(err.path(), path, "{what}: {err}");
    let variant = match err {
        ReadError::Io { .. } => "io",
        ReadError::Version { .. } => "version",
        ReadError::Invalid { .. } => "invalid",
    };
    assert_eq!(variant, kind, "{what}: {err}");
    assert_eq!(
        err.to_string(),
        format!("{}: {message}", path.display()),
        "{what}"
    );
}

#[test]
fn rejects_other_format_versions() {
    for version in ["0", "1", "3", "18446744073709551615"] {
        let (dir, err) = read_spoiled("version", |dir| {
            edit_world_json(dir, |text| {
                text.replace(
                    "\"format_version\": 2,",
                    &format!("\"format_version\": {version},"),
                )
            })
        });
        assert_error(
            version,
            &err,
            &dir,
            "world.json",
            "version",
            &format!(
                "world directory format version {version} is not supported: this build reads \
                 only version 2 (generate the world again)"
            ),
        );
        remove(&dir);
    }
    // The version is checked before the rest of the file: a version 1
    // directory of ticket M2-01, with a key this build does not know.
    let (dir, err) = read_spoiled("version_first", |dir| {
        edit_world_json(dir, |text| {
            text.replace(
                "\"format_version\": 2,",
                "\"format_version\": 1, \"old\": 0,",
            )
        })
    });
    assert_error(
        "version 1",
        &err,
        &dir,
        "world.json",
        "version",
        "world directory format version 1 is not supported: this build reads only version 2 \
         (generate the world again)",
    );
    remove(&dir);
    // So is a later version whose layout changed: structs as arrays, enums
    // with data, keys missing. Only the top level must be an object.
    let (dir, err) = read_spoiled("version_layout", |dir| {
        edit_world_json(dir, |text| {
            let text = text
                .replace("\"format_version\": 2,", "\"format_version\": 3,")
                .replace("\"seed\": 1,", "")
                .replace(
                    "\"rotation\": \"free\"",
                    "\"rotation\": {\"free\": {\"period_s\": 1.0}}",
                )
                .replace(
                    "\"kind\": \"class\"",
                    "\"kind\": {\"class\": {\"palette\": 2}}",
                );
            let start = text.find("\"settings\": {").unwrap();
            let end = start + text[start..].find('}').unwrap() + 1;
            format!(
                "{}\"settings\": [1000, 12, true]{}",
                &text[..start],
                &text[end..]
            )
        })
    });
    assert_error(
        "version 3, another layout",
        &err,
        &dir,
        "world.json",
        "version",
        "world directory format version 3 is not supported: this build reads only version 2 \
         (generate the world again)",
    );
    remove(&dir);
    // A top level that is not an object has no version to read.
    let (dir, err) = read_spoiled("version_array", |dir| {
        edit_world_json(dir, |_| "[3]".into())
    });
    assert_error(
        "an array",
        &err,
        &dir,
        "world.json",
        "invalid",
        "not valid world metadata: the top level is not a JSON object",
    );
    remove(&dir);
}

#[test]
fn rejects_malformed_world_json() {
    // Every prefix that cuts the closing brace.
    let dir = scratch("prefixes");
    write_maps_of_width(&dir, &world(), 64).unwrap();
    let path = dir.join(METADATA_FILE);
    let full = bytes_of(&path);
    let end = full.iter().rposition(|&b| b == b'}').unwrap();
    for len in 0..=end {
        std::fs::write(&path, &full[..len]).unwrap();
        let err = read(&dir).expect_err("a truncated world.json is rejected");
        assert!(
            matches!(&err, ReadError::Invalid { path: p, message }
                if *p == path && message.starts_with("not valid world metadata: ")),
            "{len}: {err}"
        );
    }
    std::fs::write(&path, &full[..len_without_newline(&full)]).unwrap();
    read(&dir).expect("the final newline is not needed");
    remove(&dir);

    let replace = |from: &'static str, to: &'static str| {
        move |text: String| {
            assert!(text.contains(from), "{from}");
            text.replacen(from, to, 1)
        }
    };
    type Edit = Box<dyn FnOnce(String) -> String>;
    let cases: Vec<(&str, Edit, &str)> = vec![
        (
            "not JSON",
            Box::new(|_| "not json".into()),
            "expected ident at line 1 column 2",
        ),
        (
            "empty",
            Box::new(|_| String::new()),
            "EOF while parsing a value at line 1 column 0",
        ),
        // serde's derive would read these arrays as the structs' values.
        (
            "an array",
            Box::new(|_| "[2]".into()),
            "the top level is not a JSON object",
        ),
        (
            "settings as an array",
            Box::new(replace(
                "\"settings\": {\n    \"cells\": 1000,\n    \"plate_count\": 12,\n    \"allow_partial\": true\n  }",
                "\"settings\": [1000, 12, true]",
            )),
            "`settings` is not a JSON object",
        ),
        (
            "planet as an array",
            Box::new(replace(
                "\"planet\": {",
                "\"planet\": [], \"old_planet\": {",
            )),
            "`planet` is not a JSON object",
        ),
        (
            "a field as an array",
            Box::new(replace(
                "{\n      \"name\": \"zone\",\n      \"unit\": \"\",\n      \"kind\": \"class\"\n    }",
                "[\"zone\", \"\", \"class\"]",
            )),
            "entry 1 of `fields` is not a JSON object",
        ),
        // And these objects as the enums' unit variants.
        (
            "a rotation as an object",
            Box::new(replace(
                "\"rotation\": \"free\"",
                "\"rotation\": {\"free\": null}",
            )),
            "`planet.rotation` is not a JSON string",
        ),
        (
            "a tectonic regime as an object",
            Box::new(replace(
                "\"tectonics\": \"plate_tectonics\"",
                "\"tectonics\": {\"plate_tectonics\": null}",
            )),
            "`planet.tectonics` is not a JSON string",
        ),
        (
            "a kind as an object",
            Box::new(replace(
                "\"kind\": \"class\"",
                "\"kind\": {\"class\": null}",
            )),
            "`kind` of entry 1 of `fields` is not a JSON string",
        ),
        (
            "no version",
            Box::new(replace("\"format_version\": 2,", "")),
            "missing field `format_version`",
        ),
        (
            "a version that is not a number",
            Box::new(replace(
                "\"format_version\": 2,",
                "\"format_version\": \"2\",",
            )),
            "`format_version` is \"2\", not an unsigned integer",
        ),
        (
            "a fractional version",
            Box::new(replace(
                "\"format_version\": 2,",
                "\"format_version\": 2.0,",
            )),
            "`format_version` is 2.0, not an unsigned integer",
        ),
        (
            "a missing key",
            Box::new(replace("\"seed\": 1,", "")),
            "missing field `seed`",
        ),
        (
            "an unknown key",
            Box::new(replace("\"seed\": 1,", "\"seed\": 1, \"colour\": 1,")),
            "unknown field `colour`",
        ),
        (
            "an unknown key of the planet",
            Box::new(replace("\"radius_m\":", "\"mass_kg\": 1.0, \"radius_m\":")),
            "unknown field `mass_kg`",
        ),
        (
            "an unknown key of the settings",
            Box::new(replace(
                "\"allow_partial\": true",
                "\"allow_partial\": true, \"extra\": 1",
            )),
            "unknown field `extra`",
        ),
        (
            "an unknown key of a field",
            Box::new(replace(
                "\"kind\": \"class\"",
                "\"kind\": \"class\", \"extra\": 1",
            )),
            "unknown field `extra`",
        ),
        (
            "a string for a number",
            Box::new(replace("\"seed\": 1,", "\"seed\": \"1\",")),
            "invalid type: string \"1\", expected u64",
        ),
        (
            "a negative seed",
            Box::new(replace("\"seed\": 1,", "\"seed\": -1,")),
            "invalid value: integer `-1`, expected u64",
        ),
        (
            "cells beyond u32",
            Box::new(replace("\"cells\": 1000,", "\"cells\": 4294967296,")),
            "invalid value: integer `4294967296`, expected u32",
        ),
        (
            "a number beyond f64",
            Box::new(replace("\"sea_level_m\": 0.0,", "\"sea_level_m\": 1e400,")),
            "number out of range",
        ),
        (
            "null for a number",
            Box::new(replace("\"sea_level_m\": 0.0,", "\"sea_level_m\": null,")),
            "invalid type: null, expected f64",
        ),
        (
            "an unknown kind",
            Box::new(replace("\"kind\": \"class\"", "\"kind\": \"vector\"")),
            "unknown variant `vector`, expected `quantity` or `class`",
        ),
        (
            "an unknown regime",
            Box::new(replace(
                "\"rotation\": \"free\"",
                "\"rotation\": \"retrograde\"",
            )),
            "unknown variant `retrograde`, expected `free` or `tidally_locked`",
        ),
        (
            "text after the object",
            Box::new(|text| text + "{}"),
            "trailing characters",
        ),
    ];
    for (what, edit, message) in cases {
        let (dir, err) = read_spoiled("malformed", |dir| edit_world_json(dir, edit));
        assert!(
            matches!(&err, ReadError::Invalid { path, message: m }
                if *path == dir.join(METADATA_FILE)
                    && m.starts_with(&format!("not valid world metadata: {message}"))),
            "{what}: {err}"
        );
        remove(&dir);
    }
    // The field list: names that are not valid (and could leave the
    // directory), a reserved name (ticket M2-19), twice the same name.
    let cases: [(&str, &str, &str); 5] = [
        (
            "\"name\": \"zone\"",
            "\"name\": \"../zone\"",
            "`../zone` is not a valid field name",
        ),
        (
            "\"name\": \"zone\"",
            "\"name\": \"Zone\"",
            "`Zone` is not a valid field name",
        ),
        (
            "\"name\": \"zone\"",
            "\"name\": \"con\"",
            "`con` is not a valid field name",
        ),
        (
            "\"name\": \"zone\"",
            "\"name\": \"hypsometry\"",
            HYPSOMETRY_RESERVED,
        ),
        (
            "\"name\": \"zone\"",
            "\"name\": \"height\"",
            "the field `height` is listed twice",
        ),
    ];
    for (from, to, message) in cases {
        let (dir, err) = read_spoiled("field_list", |dir| {
            edit_world_json(dir, |text| text.replacen(from, to, 1))
        });
        assert_error(to, &err, &dir, "world.json", "invalid", message);
        remove(&dir);
    }
}

#[test]
fn rejects_a_listed_reserved_name_even_with_its_file() {
    // A directory as a build before ticket M2-19 wrote a world with a field
    // `hypsometry` and a `surface_height`: its .kfield intact, its map
    // overwritten by the chart of the metrics.
    let dir = scratch("reserved_name");
    let mut world = world_with_relief();
    write_maps_of_width(&dir, &world, 64).unwrap();
    world.fields[1].name = "hypsometry".into();
    put(&dir, "hypsometry.kfield", &kfield(&world.fields[1]));
    std::fs::remove_file(dir.join("zone.kfield")).unwrap();
    edit_world_json(&dir, |text| {
        text.replacen("\"name\": \"zone\"", "\"name\": \"hypsometry\"", 1)
    });
    let err = read(&dir).expect_err("a reserved name is rejected");
    assert_error(
        "hypsometry",
        &err,
        &dir,
        "world.json",
        "invalid",
        HYPSOMETRY_RESERVED,
    );
    remove(&dir);
}

/// The length of `bytes` without its final newline.
fn len_without_newline(bytes: &[u8]) -> usize {
    assert_eq!(bytes.last(), Some(&b'\n'));
    bytes.len() - 1
}

/// Replaces the file `file` of `dir` with `bytes`.
fn put(dir: &Path, file: &str, bytes: &[u8]) {
    std::fs::write(dir.join(file), bytes).unwrap();
}

/// The bytes of `field` as a `.kfield` file.
fn kfield(field: &Field) -> Vec<u8> {
    let mut bytes = Vec::new();
    field_file::write_field(&mut bytes, field).unwrap();
    bytes
}

#[test]
fn rejects_bad_kfield_files() {
    // zone.kfield: 32 bytes of header (8 magic, 4 version, 4 + 4 name,
    // 4 + 0 unit, 8 count), then 1,000 values.
    assert_eq!(kfield(&world().fields[1]).len(), 32 + 8_000);
    type Spoil = Box<dyn FnOnce(&Path)>;
    let cases: Vec<(&str, Spoil, &str, &str, &str)> = vec![
        (
            "truncated",
            Box::new(|dir| put(dir, "zone.kfield", &kfield(&world().fields[1])[..50])),
            "zone.kfield",
            "invalid",
            "truncated .kfield file: it ends inside the values (after 2 of 1000)",
        ),
        (
            "empty",
            Box::new(|dir| put(dir, "zone.kfield", b"")),
            "zone.kfield",
            "invalid",
            "truncated .kfield file: it ends inside the magic number",
        ),
        (
            "trailing bytes",
            Box::new(|dir| {
                let mut bytes = kfield(&world().fields[1]);
                bytes.extend(b"abc");
                put(dir, "zone.kfield", &bytes);
            }),
            "zone.kfield",
            "invalid",
            "3 bytes follow the end of the .kfield data",
        ),
        (
            "another format",
            Box::new(|dir| put(dir, "zone.kfield", b"\x89PNG\r\n\x1a\n and so on")),
            "zone.kfield",
            "invalid",
            "not a .kfield file: it does not start with the magic number KUNIFLD\\x00",
        ),
        (
            "another version",
            Box::new(|dir| {
                let mut bytes = kfield(&world().fields[0]);
                bytes[8] = 2;
                put(dir, "height.kfield", &bytes);
            }),
            "height.kfield",
            "version",
            ".kfield format version 2 is not supported: this build reads only version 1",
        ),
        (
            "another field's file",
            Box::new(|dir| put(dir, "zone.kfield", &kfield(&world().fields[0]))),
            "zone.kfield",
            "invalid",
            "holds the field `height`, but world.json lists `zone`",
        ),
        (
            "another unit",
            Box::new(|dir| {
                let mut field = world().fields[0].clone();
                field.unit = "km".into();
                put(dir, "height.kfield", &kfield(&field));
            }),
            "height.kfield",
            "invalid",
            "has the unit `km`, but world.json gives `m`",
        ),
        (
            "too few values",
            Box::new(|dir| {
                let mut field = world().fields[0].clone();
                field.values.pop();
                put(dir, "height.kfield", &kfield(&field));
            }),
            "height.kfield",
            "invalid",
            "has 999 values for the 1000 cells of the grid",
        ),
        (
            "not a class",
            Box::new(|dir| {
                let mut field = world().fields[1].clone();
                field.values[3] = 0.5;
                put(dir, "zone.kfield", &kfield(&field));
            }),
            "zone.kfield",
            "invalid",
            "class field `zone` has the value 0.5 at cell 3, not a class id",
        ),
        (
            "negative zero in a class field",
            Box::new(|dir| {
                let mut field = world().fields[1].clone();
                field.values[0] = -0.0;
                put(dir, "zone.kfield", &kfield(&field));
            }),
            "zone.kfield",
            "invalid",
            "class field `zone` has the value -0 at cell 0, not a class id",
        ),
    ];
    for (what, spoil, file, kind, message) in cases {
        let (dir, err) = read_spoiled("kfield", spoil);
        assert_error(what, &err, &dir, file, kind, message);
        remove(&dir);
    }
    // A missing file is an error of the file system.
    let (dir, err) = read_spoiled("kfield_missing", |dir| {
        std::fs::remove_file(dir.join("zone.kfield")).unwrap();
    });
    assert!(
        matches!(&err, ReadError::Io { path, source }
            if *path == dir.join("zone.kfield") && source.kind() == io::ErrorKind::NotFound),
        "{err}"
    );
    remove(&dir);
    // So is a directory in its place (it can be opened on some systems,
    // but not read).
    let (dir, err) = read_spoiled("kfield_dir", |dir| {
        std::fs::remove_file(dir.join("zone.kfield")).unwrap();
        std::fs::create_dir(dir.join("zone.kfield")).unwrap();
    });
    assert!(
        matches!(&err, ReadError::Io { path, .. } if *path == dir.join("zone.kfield")),
        "{err}"
    );
    remove(&dir);
}

#[test]
fn rejects_bad_grid_files() {
    let grid = || {
        let mut bytes = Vec::new();
        grid_file::write_grid(&mut bytes, &world().grid).unwrap();
        bytes
    };
    // 28 bytes of header, 1,000 centres of 24 bytes, then the offsets.
    type Spoil = Box<dyn FnOnce(&Path)>;
    let cases: Vec<(&str, Spoil, &str, &str, &str)> = vec![
        (
            "truncated",
            Box::new(move |dir| put(dir, GRID_FILE, &grid()[..28 + 24_000 + 10])),
            GRID_FILE,
            "invalid",
            "truncated .kgrid file: it ends inside the offsets (after 2 of 1001)",
        ),
        (
            "trailing bytes",
            Box::new(move |dir| {
                let mut bytes = grid();
                bytes.push(0);
                put(dir, GRID_FILE, &bytes);
            }),
            GRID_FILE,
            "invalid",
            "1 byte follows the end of the .kgrid data",
        ),
        (
            "another version",
            Box::new(move |dir| {
                let mut bytes = grid();
                bytes[8] = 9;
                put(dir, GRID_FILE, &bytes);
            }),
            GRID_FILE,
            "version",
            ".kgrid format version 9 is not supported: this build reads only version 1",
        ),
        (
            "a field file in its place",
            Box::new(|dir| put(dir, GRID_FILE, &kfield(&world().fields[0]))),
            GRID_FILE,
            "invalid",
            "not a .kgrid file: it does not start with the magic number KUNIGRD\\x00",
        ),
        (
            "not a valid grid",
            Box::new(move |dir| {
                // The first area (after the centres, the 1,001 offsets and
                // the 5,988 neighbours) doubled.
                let mut bytes = grid();
                let at = 28 + 24_000 + 4 * 1_001 + 4 * 5_988;
                let area = f64::from_le_bytes(bytes[at..at + 8].try_into().unwrap());
                bytes[at..at + 8].copy_from_slice(&(2.0 * area).to_le_bytes());
                put(dir, GRID_FILE, &bytes);
            }),
            GRID_FILE,
            "invalid",
            "the .kgrid file does not hold a valid grid: cell areas are not positive and finite \
             or do not sum to 4π",
        ),
        (
            "the grid of another world",
            Box::new(|dir| {
                let mut bytes = Vec::new();
                let other = kuniumi_sphere::SphereGrid::fibonacci(1_500, 1).unwrap();
                grid_file::write_grid(&mut bytes, &other).unwrap();
                put(dir, GRID_FILE, &bytes);
            }),
            "height.kfield",
            "invalid",
            "has 1000 values for the 1500 cells of the grid",
        ),
    ];
    for (what, spoil, file, kind, message) in cases {
        let (dir, err) = read_spoiled("grid", spoil);
        assert_error(what, &err, &dir, file, kind, message);
        remove(&dir);
    }
    let (dir, err) = read_spoiled("grid_missing", |dir| {
        std::fs::remove_file(dir.join(GRID_FILE)).unwrap();
    });
    assert!(
        matches!(&err, ReadError::Io { path, source }
            if *path == dir.join(GRID_FILE) && source.kind() == io::ErrorKind::NotFound),
        "{err}"
    );
    remove(&dir);
}

#[test]
fn a_missing_directory_is_an_io_error() {
    let dir = scratch("missing");
    let err = read(&dir).unwrap_err();
    let path = dir.join(METADATA_FILE);
    assert!(
        matches!(&err, ReadError::Io { path: p, source }
            if *p == path && source.kind() == io::ErrorKind::NotFound),
        "{err}"
    );
    let source = std::error::Error::source(&err).expect("the error of the file system");
    assert_eq!(err.to_string(), format!("{}: {source}", path.display()));
}

#[test]
fn errors_name_their_file() {
    let path = PathBuf::from("some/dir/zone.kfield");
    let err = ReadError::Invalid {
        path: path.clone(),
        message: "what is wrong".into(),
    };
    assert_eq!(err.path(), path);
    assert_eq!(err.to_string(), "some/dir/zone.kfield: what is wrong");
    assert!(std::error::Error::source(&err).is_none());
    let err = ReadError::Version {
        path: path.clone(),
        message: "a version".into(),
    };
    assert_eq!(err.to_string(), "some/dir/zone.kfield: a version");
    assert!(std::error::Error::source(&err).is_none());
    let err = ReadError::Io {
        path: path.clone(),
        source: io::Error::other("no disk"),
    };
    assert_eq!(err.path(), path);
    assert_eq!(err.to_string(), "some/dir/zone.kfield: no disk");
    assert!(std::error::Error::source(&err).is_some());
}

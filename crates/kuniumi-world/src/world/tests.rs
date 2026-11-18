//! Tests of `World`'s field helpers, of `Field`'s name and class rules, and
//! of the warnings of the stages (ticket M4-16).

use super::*;

/// A world on the octahedron (6 cells) with no fields.
fn world() -> World {
    World {
        params: PlanetParams::earth(),
        settings: Settings::default(),
        seed: 1,
        grid: SphereGrid::octahedron(),
        sea_level_m: 0.0,
        stages: Vec::new(),
        warnings: Vec::new(),
        fields: Vec::new(),
    }
}

#[test]
fn adds_fields_in_order() {
    let mut w = world();
    w.add_field("surface_height", "m", vec![1.0; 6]).unwrap();
    w.add_class_field("plate_id", vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
        .unwrap();
    w.add_field("mask", "", vec![0.5; 6]).unwrap();
    let names: Vec<_> = w.fields.iter().map(|f| f.name.as_str()).collect();
    assert_eq!(names, ["surface_height", "plate_id", "mask"]);
    let h = w.field("surface_height").unwrap();
    assert_eq!((h.unit.as_str(), h.kind), ("m", FieldKind::Quantity));
    let p = w.field("plate_id").unwrap();
    assert_eq!((p.unit.as_str(), p.kind), ("", FieldKind::Class));
    assert_eq!(p.values, [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]);
    assert!(w.field("missing").is_none());
}

#[test]
fn field_mut_updates_in_place() {
    let mut w = world();
    w.add_field("surface_height", "m", vec![1.0; 6]).unwrap();
    w.field_mut("surface_height").unwrap().values[2] = -4_000.0;
    assert_eq!(
        w.field("surface_height").unwrap().values,
        [1.0, 1.0, -4_000.0, 1.0, 1.0, 1.0]
    );
    assert!(w.field_mut("missing").is_none());
}

#[test]
fn rejects_duplicates_wrong_lengths_and_bad_names() {
    let mut w = world();
    w.add_field("a", "m", vec![0.0; 6]).unwrap();
    let before = w.clone();
    assert_eq!(
        w.add_field("a", "K", vec![0.0; 6]),
        Err(FieldError::Duplicate { name: "a".into() })
    );
    assert_eq!(
        w.add_class_field("a", vec![0.0; 6]),
        Err(FieldError::Duplicate { name: "a".into() })
    );
    for got in [0, 5, 7] {
        assert_eq!(
            w.add_field("b", "m", vec![0.0; got]),
            Err(FieldError::Length {
                name: "b".into(),
                expected: 6,
                got
            })
        );
    }
    for name in [
        "", "A", "1a", "_a", "a-b", "a b", "a.kfield", "../a", "é", "aB",
    ] {
        assert_eq!(
            w.add_field(name, "", vec![0.0; 6]),
            Err(FieldError::InvalidName { name: name.into() }),
            "{name:?}"
        );
    }
    assert_eq!(w, before);
}

#[test]
fn valid_names() {
    for name in ["a", "z", "plate_id", "temperature_m01", "a_", "a__b", "x9"] {
        assert!(Field::is_valid_name(name), "{name:?}");
    }
    // Near the reserved device names of Windows, but not one of them.
    for name in [
        "co", "cons", "con_", "nul1", "com", "com10", "lpt", "lpt_1", "comx",
    ] {
        assert!(Field::is_valid_name(name), "{name:?}");
    }
}

#[test]
fn rejects_the_device_names_of_windows() {
    let mut reserved = vec!["con".to_owned(), "prn".into(), "aux".into(), "nul".into()];
    for d in 0..10 {
        reserved.push(format!("com{d}"));
        reserved.push(format!("lpt{d}"));
    }
    for name in &reserved {
        assert!(!Field::is_valid_name(name), "{name:?}");
    }
    let mut w = world();
    assert_eq!(
        w.add_field("nul", "", vec![0.0; 6]),
        Err(FieldError::InvalidName { name: "nul".into() })
    );
}

#[test]
fn class_values_are_non_negative_integers_below_2_pow_53() {
    for v in [
        0.0,
        1.0,
        11.0,
        12.0,
        4_294_967_295.0,
        9_007_199_254_740_991.0,
    ] {
        assert!(Field::is_class_value(v), "{v}");
    }
    for v in [
        -0.0,
        -1.0,
        0.5,
        1.0 + f64::EPSILON,
        9_007_199_254_740_992.0,
        f64::MAX,
        f64::INFINITY,
        f64::NAN,
        f64::MIN_POSITIVE,
    ] {
        assert!(!Field::is_class_value(v), "{v}");
    }
}

#[test]
fn rejects_the_first_value_that_is_not_a_class() {
    let mut w = world();
    let err = w
        .add_class_field("plate_id", vec![0.0, 1.0, 2.5, -1.0, 0.0, 0.0])
        .unwrap_err();
    assert_eq!(
        err,
        FieldError::NotAClass {
            name: "plate_id".into(),
            cell: 2,
            value: 2.5
        }
    );
    let err = w
        .add_class_field("plate_id", vec![0.0, 1.0, 2.0, 3.0, 4.0, -0.0])
        .unwrap_err();
    assert!(
        matches!(err, FieldError::NotAClass { cell: 5, .. }),
        "{err:?}"
    );
    // A quantity field takes any value.
    w.add_field("plate_id", "", vec![0.0, 1.0, 2.5, -1.0, f64::NAN, -0.0])
        .unwrap();
    assert!(w.fields.len() == 1 && w.fields[0].kind == FieldKind::Quantity);
}

#[test]
fn errors_display() {
    let cases = [
        (
            FieldError::InvalidName { name: "A".into() },
            "`A` is not a valid field name (lowercase letters, digits and `_`, \
             starting with a letter, not a Windows device name)",
        ),
        (
            FieldError::Duplicate { name: "a".into() },
            "the world already has a field `a`",
        ),
        (
            FieldError::Length {
                name: "a".into(),
                expected: 6,
                got: 5,
            },
            "field `a` has 5 values for 6 cells",
        ),
        (
            FieldError::NotAClass {
                name: "plate_id".into(),
                cell: 2,
                value: 2.5,
            },
            "class field `plate_id` has the value 2.5 at cell 2, not a class id",
        ),
    ];
    for (err, text) in cases {
        assert_eq!(err.to_string(), text);
    }
}

#[test]
fn warnings_are_part_of_the_world() {
    let w = world();
    assert!(w.warnings.is_empty());
    let warning = |stage: &str, message: &str| Warning {
        stage: stage.into(),
        message: message.into(),
    };
    let mut a = world();
    a.warnings = vec![warning("climate", "one"), warning("climate", "two")];
    assert_eq!(a.clone(), a);
    assert_ne!(a, w);
    // A world differs from another in a warning's stage, its message and
    // the order of the warnings.
    let mut b = a.clone();
    b.warnings[0].stage = "erosion".into();
    assert_ne!(a, b);
    let mut b = a.clone();
    b.warnings[1].message = "two ".into();
    assert_ne!(a, b);
    let mut b = a.clone();
    b.warnings.reverse();
    assert_ne!(a, b);
}

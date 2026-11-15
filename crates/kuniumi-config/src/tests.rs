use kuniumi_math::Rng;
use kuniumi_world::{ErosionSettings, RotationRegime, TectonicRegime};

use super::*;

/// The settings file of Earth with the default settings, written by hand:
/// the decimal literals of `PlanetParams::earth()`, `Settings::default()`
/// and `ErosionSettings::default()`, the obliquity from Python's
/// `repr(23.44 * (math.pi / 180))`, and the default diffusivity and time
/// step from `repr(5e-3 / 31557600.0)` and `repr(1e4 * 31557600.0)`.
const EARTH: &str = "\
# kuniumi settings file: a planet and how to generate it (TOML).
# Values are SI. A quantity whose comment names another key may be given
# under that key instead, in its unit, but not under both. Missing keys take
# Earth's values and the default settings. The seed is not part of the file.

[planet]
radius_m = 6371000.0                # mean radius; or radius_km
surface_gravity_m_s2 = 9.80665      # at the surface
rotation_period_s = 86164.0905      # sidereal day; or rotation_period_h
obliquity_rad = 0.40910517666747087 # axial tilt; or obliquity_deg
orbital_period_s = 31558149.8       # year; or orbital_period_days
orbital_eccentricity = 0.0167       # 0 for a circular orbit
stellar_flux_w_m2 = 1361.0          # mean flux at the orbit
ocean_volume_m3 = 1.335e18          # surface water; or ocean_volume_km3
surface_pressure_pa = 101325.0      # mean air pressure; or surface_pressure_bar
continental_crust_fraction = 0.4    # share of the surface
rotation = \"free\"                   # \"free\" or \"tidally_locked\" (v2)
tectonics = \"plate_tectonics\"       # \"plate_tectonics\" or \"stagnant_lid\" (v2)

[generation]
cells = 650000                      # grid cells, 1000 to 715827884
plate_count = 12                    # tectonic plates, 2 to 64
allow_partial = false               # return the world so far at a missing stage

[generation.erosion]
erodibility = 7.5e-10               # stream power K, m^(1-3m) s^(m-1)
discharge_exponent = 0.45           # m of the discharge in K Q^m S
diffusivity_m2_s = 1.5844043907014476e-10 # hillslope; or diffusivity_m2_yr
time_step_s = 315576000000.0        # time step; or time_step_yr
steps = 1000                        # number of time steps
";

/// The bits of every float of `p`, field by field. The destructuring lists
/// every field, so a new field fails to compile here until it is added.
fn bits(p: &PlanetParams) -> [u64; 10] {
    let PlanetParams {
        radius_m,
        surface_gravity_m_s2,
        rotation_period_s,
        obliquity_rad,
        orbital_period_s,
        orbital_eccentricity,
        stellar_flux_w_m2,
        ocean_volume_m3,
        surface_pressure_pa,
        continental_crust_fraction,
        rotation: _,
        tectonics: _,
    } = p;
    [
        radius_m,
        surface_gravity_m_s2,
        rotation_period_s,
        obliquity_rad,
        orbital_period_s,
        orbital_eccentricity,
        stellar_flux_w_m2,
        ocean_volume_m3,
        surface_pressure_pa,
        continental_crust_fraction,
    ]
    .map(|v| v.to_bits())
}

/// The bits of every float of `s`, field by field, listed as [`bits`]
/// lists those of a planet.
fn settings_bits(s: &Settings) -> [u64; 4] {
    let Settings {
        cells: _,
        plate_count: _,
        allow_partial: _,
        erosion,
    } = s;
    let ErosionSettings {
        erodibility,
        discharge_exponent,
        diffusivity_m2_s,
        time_step_s,
        steps: _,
    } = erosion;
    [
        erodibility,
        discharge_exponent,
        diffusivity_m2_s,
        time_step_s,
    ]
    .map(|v| v.to_bits())
}

/// Asserts that `loaded` is `(p, s)` bit for bit: every float by its bits
/// (so that −0 is not +0), and everything by `==`.
fn assert_same(loaded: &(PlanetParams, Settings), p: &PlanetParams, s: &Settings, text: &str) {
    assert_eq!(bits(&loaded.0), bits(p), "{text}");
    assert_eq!(loaded.0, *p, "{text}");
    assert_eq!(settings_bits(&loaded.1), settings_bits(s), "{text}");
    assert_eq!(loaded.1, *s, "{text}");
}

/// Loads `text`, which must be valid.
fn ok(text: &str) -> (PlanetParams, Settings) {
    load(text).unwrap_or_else(|e| panic!("{e}\n{text}"))
}

/// Loads `text`, which must fail.
fn err(text: &str) -> ConfigError {
    match load(text) {
        Ok(v) => panic!("loaded {v:?} from\n{text}"),
        Err(e) => e,
    }
}

/// The report of the toml crate for `text`, which must fail there.
fn toml_message(text: &str) -> String {
    match err(text) {
        ConfigError::Toml { message } => message,
        e => panic!("{e:?} from\n{text}"),
    }
}

#[test]
fn earth_file_has_the_documented_layout() {
    let text = write(&PlanetParams::earth(), &Settings::default());
    assert_eq!(text, EARTH);
    // Every line, comments included, fits 80 columns.
    for line in text.lines() {
        assert!(line.chars().count() <= 80, "{line}");
    }
}

#[test]
fn earth_round_trips_bit_exactly() {
    let (p, s) = (PlanetParams::earth(), Settings::default());
    assert_same(&ok(&write(&p, &s)), &p, &s, EARTH);
    assert_same(&ok(EARTH), &p, &s, EARTH);
}

/// The next double above the positive finite `v`.
fn next_up(v: f64) -> f64 {
    f64::from_bits(v.to_bits() + 1)
}

/// The next double below the positive finite `v`.
fn next_down(v: f64) -> f64 {
    f64::from_bits(v.to_bits() - 1)
}

/// Planets of awkward numbers, each valid: the extremes of every range,
/// subnormals and −0 where zero is allowed (the ocean volume, the
/// obliquity, the eccentricity and the continental fraction; the discharge
/// exponent and the diffusivity), 17-digit values and numbers whose
/// shortest decimal needs an exponent.
fn awkward_planets() -> Vec<(PlanetParams, Settings)> {
    let earth = PlanetParams::earth();
    let tiny = f64::from_bits(1); // 5e-324
    let largest_subnormal = f64::from_bits((1 << 52) - 1);
    vec![
        (
            PlanetParams {
                radius_m: PlanetParams::MIN_RADIUS_M,
                // 3.0000000000000004: 7,154 kg/m³ at 1,500 km.
                surface_gravity_m_s2: (0.1 + 0.2) * 10.0,
                rotation_period_s: PlanetParams::MIN_ROTATION_PERIOD_S,
                obliquity_rad: -0.0,
                orbital_period_s: PlanetParams::MIN_ORBITAL_PERIOD_S,
                orbital_eccentricity: -0.0,
                stellar_flux_w_m2: 4_000.0 / 3.0,
                ocean_volume_m3: -0.0,
                surface_pressure_pa: PlanetParams::MIN_SURFACE_PRESSURE_PA,
                continental_crust_fraction: -0.0,
                ..earth.clone()
            },
            Settings {
                cells: Settings::MIN_CELLS,
                plate_count: Settings::MIN_PLATE_COUNT,
                allow_partial: true,
                erosion: ErosionSettings {
                    erodibility: f64::MIN_POSITIVE,
                    discharge_exponent: -0.0,
                    diffusivity_m2_s: -0.0,
                    time_step_s: tiny,
                    steps: 0,
                },
            },
        ),
        (
            PlanetParams {
                radius_m: PlanetParams::MAX_RADIUS_M,
                // 11,006 kg/m³ at 13,000 km.
                surface_gravity_m_s2: PlanetParams::MAX_SURFACE_GRAVITY_M_S2,
                rotation_period_s: PlanetParams::MAX_ROTATION_PERIOD_S,
                obliquity_rad: core::f64::consts::PI,
                orbital_period_s: PlanetParams::MAX_ORBITAL_PERIOD_S,
                orbital_eccentricity: 0.99,
                stellar_flux_w_m2: PlanetParams::MAX_STELLAR_FLUX_W_M2,
                ocean_volume_m3: PlanetParams::MAX_OCEAN_VOLUME_M3, // 1e20
                surface_pressure_pa: PlanetParams::MAX_SURFACE_PRESSURE_PA,
                continental_crust_fraction: 1.0,
                ..earth.clone()
            },
            Settings {
                cells: Settings::MAX_CELLS,
                plate_count: Settings::MAX_PLATE_COUNT,
                allow_partial: false,
                // The largest diffusion step, f64::MAX m².
                erosion: ErosionSettings {
                    erodibility: f64::MAX,
                    discharge_exponent: f64::MAX,
                    diffusivity_m2_s: 1.0,
                    time_step_s: f64::MAX,
                    steps: u32::MAX,
                },
            },
        ),
        (
            PlanetParams {
                radius_m: next_up(earth.radius_m),
                surface_gravity_m_s2: next_up(earth.surface_gravity_m_s2),
                rotation_period_s: next_down(earth.rotation_period_s),
                obliquity_rad: tiny,
                orbital_period_s: 123_456_789.012_345_68,
                orbital_eccentricity: tiny,
                stellar_flux_w_m2: next_down(PlanetParams::MAX_STELLAR_FLUX_W_M2),
                ocean_volume_m3: tiny,
                surface_pressure_pa: next_up(PlanetParams::MIN_SURFACE_PRESSURE_PA),
                continental_crust_fraction: largest_subnormal,
                ..earth
            },
            Settings {
                cells: 1_001,
                plate_count: 3,
                allow_partial: true,
                erosion: ErosionSettings {
                    erodibility: tiny,
                    discharge_exponent: 1.0 / 3.0,
                    diffusivity_m2_s: largest_subnormal,
                    time_step_s: 123_456_789.012_345_68,
                    steps: 1,
                },
            },
        ),
    ]
}

/// A random finite float in `[lo, hi]`: half of them uniform in the range,
/// half from random bit patterns (any exponent and mantissa, sign clear).
fn within(rng: &mut Rng, lo: f64, hi: f64) -> f64 {
    loop {
        let v = if rng.below(2) == 0 {
            rng.range_f64(lo, hi)
        } else {
            f64::from_bits(rng.next_u64() >> 1)
        };
        if v.is_finite() && v >= lo && v <= hi {
            return v;
        }
    }
}

/// A random valid planet with random settings: every float within its
/// bounds, a radius and gravity whose bulk density is in range, and a
/// diffusion step `Δt K_d` that is finite.
fn random_planet(rng: &mut Rng) -> (PlanetParams, Settings) {
    let p = loop {
        let p = PlanetParams {
            radius_m: within(rng, PlanetParams::MIN_RADIUS_M, PlanetParams::MAX_RADIUS_M),
            surface_gravity_m_s2: within(
                rng,
                PlanetParams::MIN_SURFACE_GRAVITY_M_S2,
                PlanetParams::MAX_SURFACE_GRAVITY_M_S2,
            ),
            rotation_period_s: within(
                rng,
                PlanetParams::MIN_ROTATION_PERIOD_S,
                PlanetParams::MAX_ROTATION_PERIOD_S,
            ),
            obliquity_rad: within(rng, 0.0, core::f64::consts::PI),
            orbital_period_s: within(
                rng,
                PlanetParams::MIN_ORBITAL_PERIOD_S,
                PlanetParams::MAX_ORBITAL_PERIOD_S,
            ),
            orbital_eccentricity: within(rng, 0.0, 0.99),
            stellar_flux_w_m2: within(
                rng,
                PlanetParams::MIN_STELLAR_FLUX_W_M2,
                PlanetParams::MAX_STELLAR_FLUX_W_M2,
            ),
            ocean_volume_m3: within(rng, 0.0, PlanetParams::MAX_OCEAN_VOLUME_M3),
            surface_pressure_pa: within(
                rng,
                PlanetParams::MIN_SURFACE_PRESSURE_PA,
                PlanetParams::MAX_SURFACE_PRESSURE_PA,
            ),
            continental_crust_fraction: within(rng, 0.0, 1.0),
            rotation: RotationRegime::Free,
            tectonics: TectonicRegime::PlateTectonics,
        };
        let rho = p.bulk_density_kg_m3();
        if (PlanetParams::MIN_BULK_DENSITY_KG_M3..=PlanetParams::MAX_BULK_DENSITY_KG_M3)
            .contains(&rho)
        {
            break p;
        }
    };
    let span = |lo: u32, hi: u32| u64::from(hi - lo) + 1;
    let s = Settings {
        cells: Settings::MIN_CELLS
            + rng.below(span(Settings::MIN_CELLS, Settings::MAX_CELLS)) as u32,
        plate_count: Settings::MIN_PLATE_COUNT
            + rng.below(span(Settings::MIN_PLATE_COUNT, Settings::MAX_PLATE_COUNT)) as u32,
        allow_partial: rng.below(2) == 1,
        erosion: loop {
            let tiny = f64::from_bits(1);
            let e = ErosionSettings {
                erodibility: within(rng, tiny, f64::MAX),
                discharge_exponent: within(rng, 0.0, f64::MAX),
                diffusivity_m2_s: within(rng, 0.0, f64::MAX),
                time_step_s: within(rng, tiny, f64::MAX),
                steps: rng.below(1 << 32) as u32,
            };
            if (e.time_step_s * e.diffusivity_m2_s).is_finite() {
                break e;
            }
        },
    };
    (p, s)
}

#[test]
fn random_and_awkward_planets_round_trip_bit_exactly() {
    let mut rng = Rng::from_seed(7);
    let mut planets = awkward_planets();
    planets.extend((0..10_000).map(|_| random_planet(&mut rng)));
    for (p, s) in &planets {
        assert_eq!(p.validate(), Ok(()), "{p:?}");
        assert_eq!(s.validate(), Ok(()), "{s:?}");
        let text = write(p, s);
        assert_same(&ok(&text), p, s, &text);
    }
}

#[test]
fn every_finite_float_is_written_as_a_toml_float_that_parses_back() {
    // `write` formats floats itself; each must be a TOML float (not an
    // integer, so that values beyond i64 parse) with the same bits.
    let mut rng = Rng::from_seed(11);
    let mut values = vec![0.0, -0.0, f64::MAX, -f64::MAX, f64::MIN_POSITIVE];
    values.extend((0..100_000).map(|_| f64::from_bits(rng.next_u64())));
    for v in values.into_iter().filter(|v| v.is_finite()) {
        let written = float(v);
        let table: toml::Table =
            toml::from_str(&format!("x = {written}")).unwrap_or_else(|e| panic!("{written}: {e}"));
        match table.get("x") {
            Some(toml::Value::Float(back)) => {
                assert_eq!(back.to_bits(), v.to_bits(), "{written}");
            }
            other => panic!("{written} read as {other:?}"),
        }
    }
}

#[test]
fn invalid_parameters_are_written_and_rejected_by_key() {
    for (p, s, key) in [
        (
            PlanetParams {
                radius_m: f64::NAN,
                ..PlanetParams::earth()
            },
            Settings::default(),
            "planet.radius_m",
        ),
        (
            PlanetParams {
                stellar_flux_w_m2: f64::INFINITY,
                ..PlanetParams::earth()
            },
            Settings::default(),
            "planet.stellar_flux_w_m2",
        ),
        (
            PlanetParams {
                ocean_volume_m3: f64::NEG_INFINITY,
                ..PlanetParams::earth()
            },
            Settings::default(),
            "planet.ocean_volume_m3",
        ),
        (
            PlanetParams {
                tectonics: TectonicRegime::StagnantLid,
                ..PlanetParams::earth()
            },
            Settings::default(),
            "planet.tectonics",
        ),
        (
            PlanetParams::earth(),
            Settings {
                plate_count: 0,
                ..Settings::default()
            },
            "generation.plate_count",
        ),
        (
            PlanetParams::earth(),
            Settings {
                erosion: ErosionSettings {
                    erodibility: f64::NAN,
                    ..ErosionSettings::default()
                },
                ..Settings::default()
            },
            "generation.erosion.erodibility",
        ),
        (
            PlanetParams::earth(),
            Settings {
                erosion: ErosionSettings {
                    time_step_s: f64::NEG_INFINITY,
                    ..ErosionSettings::default()
                },
                ..Settings::default()
            },
            "generation.erosion.time_step_s",
        ),
    ] {
        let text = write(&p, &s);
        match err(&text) {
            ConfigError::Invalid { key: k, .. } => assert_eq!(k, key, "{text}"),
            e => panic!("{e:?} from\n{text}"),
        }
    }
    let text = write(
        &PlanetParams {
            radius_m: f64::NAN,
            stellar_flux_w_m2: f64::INFINITY,
            ocean_volume_m3: f64::NEG_INFINITY,
            ..PlanetParams::earth()
        },
        &Settings::default(),
    );
    for line in [
        "radius_m = nan ",
        "stellar_flux_w_m2 = inf ",
        "ocean_volume_m3 = -inf ",
    ] {
        assert!(text.contains(&format!("\n{line}")), "{line}\n{text}");
    }
    let text = write(
        &PlanetParams::earth(),
        &Settings {
            erosion: ErosionSettings {
                discharge_exponent: f64::NAN,
                diffusivity_m2_s: f64::INFINITY,
                ..ErosionSettings::default()
            },
            ..Settings::default()
        },
    );
    for line in ["discharge_exponent = nan ", "diffusivity_m2_s = inf "] {
        assert!(text.contains(&format!("\n{line}")), "{line}\n{text}");
    }
    // A diffusion step that overflows is written as it is, and rejected.
    let s = Settings {
        erosion: ErosionSettings {
            diffusivity_m2_s: f64::MAX,
            ..ErosionSettings::default()
        },
        ..Settings::default()
    };
    let text = write(&PlanetParams::earth(), &s);
    assert!(
        text.contains("\ndiffusivity_m2_s = 1.7976931348623157e308 "),
        "{text}"
    );
    assert!(
        matches!(err(&text), ConfigError::Inconsistent { .. }),
        "{text}"
    );
}

#[test]
fn human_unit_keys_convert_to_si() {
    // Earth in human units. Hand-computed SI values: 6,371 km = 6,371,000 m
    // and 1.335e9 km³ = 1.335e18 m³ exactly; the rest are the IEEE 754
    // products from Python, repr(23.9344696 * 3600.0),
    // repr(23.44 * (math.pi / 180)), repr(365.256363004 * 86400.0),
    // repr(1.01325 * 1e5).
    let (p, s) = ok("\
[planet]
radius_km = 6371
rotation_period_h = 23.9344696
obliquity_deg = 23.44
orbital_period_days = 365.256363004
ocean_volume_km3 = 1.335e9
surface_pressure_bar = 1.01325
");
    let expected = PlanetParams {
        radius_m: 6_371_000.0,
        rotation_period_s: 86_164.090_56,
        obliquity_rad: 0.409_105_176_667_470_87,
        orbital_period_s: 31_558_149.763_545_6,
        ocean_volume_m3: 1.335e18,
        surface_pressure_pa: 101_325.0,
        ..PlanetParams::earth()
    };
    assert_eq!(bits(&p), bits(&expected));
    assert_eq!(p, expected);
    assert_eq!(s, Settings::default());
    // 23.44° is Earth's obliquity bit for bit.
    assert_eq!(
        p.obliquity_rad.to_bits(),
        PlanetParams::earth().obliquity_rad.to_bits()
    );
    // The bounds of a range under a human-unit key: 180° is π exactly
    // (Python: 180 * (math.pi / 180) == math.pi), and 0 km³ is 0 m³.
    let (p, _) = ok("[planet]\nobliquity_deg = 180\nocean_volume_km3 = 0\n");
    assert_eq!(p.obliquity_rad.to_bits(), core::f64::consts::PI.to_bits());
    assert_eq!(p.ocean_volume_m3.to_bits(), 0.0_f64.to_bits());

    // Mars-like values, with a thicker air than Mars's 0.006 bar, which
    // is below the bounds: 3,389.5 km = 3,389,500 m exactly; from Python,
    // repr(24.6229 * 3600.0), repr(25.19 * (math.pi / 180)),
    // repr(686.98 * 86400.0), repr(0.7 * 1e5).
    let (p, _) = ok("\
[planet]
radius_km = 3389.5
surface_gravity_m_s2 = 3.721
rotation_period_h = 24.6229
obliquity_deg = 25.19
orbital_period_days = 686.98
surface_pressure_bar = 0.7
");
    let expected = PlanetParams {
        radius_m: 3_389_500.0,
        surface_gravity_m_s2: 3.721,
        rotation_period_s: 88_642.44,
        obliquity_rad: 0.439_648_438_577_371_6,
        orbital_period_s: 59_355_072.0,
        surface_pressure_pa: 70_000.0,
        ..PlanetParams::earth()
    };
    assert_eq!(bits(&p), bits(&expected));
}

#[test]
fn the_bounds_in_human_units_load_to_the_si_bounds() {
    // Each bound, written as a round number under its human-unit key,
    // loads to the constant bit for bit. The radius and gravity keep a
    // bulk density in range: 2.5 m/s² at 1,500 km is 5,962 kg/m³, 20 m/s²
    // at 13,000 km 5,503 kg/m³ (hand calculation, 3g / (4πGR)).
    let get = |p: &PlanetParams, key: &str| match key {
        "radius_km" => p.radius_m,
        "rotation_period_h" => p.rotation_period_s,
        "orbital_period_days" => p.orbital_period_s,
        "ocean_volume_km3" => p.ocean_volume_m3,
        "surface_pressure_bar" => p.surface_pressure_pa,
        "surface_gravity_m_s2" => p.surface_gravity_m_s2,
        "stellar_flux_w_m2" => p.stellar_flux_w_m2,
        _ => unreachable!("{key}"),
    };
    for (line, key, bound) in [
        (
            "radius_km = 1500\nsurface_gravity_m_s2 = 2.5",
            "radius_km",
            PlanetParams::MIN_RADIUS_M,
        ),
        (
            "radius_km = 13000\nsurface_gravity_m_s2 = 20",
            "radius_km",
            PlanetParams::MAX_RADIUS_M,
        ),
        (
            "radius_km = 13000\nsurface_gravity_m_s2 = 40",
            "surface_gravity_m_s2",
            PlanetParams::MAX_SURFACE_GRAVITY_M_S2,
        ),
        (
            "rotation_period_h = 8",
            "rotation_period_h",
            PlanetParams::MIN_ROTATION_PERIOD_S,
        ),
        (
            "rotation_period_h = 240",
            "rotation_period_h",
            PlanetParams::MAX_ROTATION_PERIOD_S,
        ),
        (
            "orbital_period_days = 90",
            "orbital_period_days",
            PlanetParams::MIN_ORBITAL_PERIOD_S,
        ),
        (
            "orbital_period_days = 3650",
            "orbital_period_days",
            PlanetParams::MAX_ORBITAL_PERIOD_S,
        ),
        (
            "stellar_flux_w_m2 = 400",
            "stellar_flux_w_m2",
            PlanetParams::MIN_STELLAR_FLUX_W_M2,
        ),
        (
            "stellar_flux_w_m2 = 2500",
            "stellar_flux_w_m2",
            PlanetParams::MAX_STELLAR_FLUX_W_M2,
        ),
        (
            "ocean_volume_km3 = 1e11",
            "ocean_volume_km3",
            PlanetParams::MAX_OCEAN_VOLUME_M3,
        ),
        (
            "surface_pressure_bar = 0.1",
            "surface_pressure_bar",
            PlanetParams::MIN_SURFACE_PRESSURE_PA,
        ),
        (
            "surface_pressure_bar = 10",
            "surface_pressure_bar",
            PlanetParams::MAX_SURFACE_PRESSURE_PA,
        ),
    ] {
        let text = format!("[planet]\n{line}\n");
        let (p, _) = ok(&text);
        assert_eq!(get(&p, key).to_bits(), bound.to_bits(), "{text}");
    }
    // Just beyond a bound, the error names the key as written; the SI
    // value is the written number times the unit (`kuniumi_math::units`).
    for (line, key, written, name, factor) in [
        (
            "radius_km = 1499.999\nsurface_gravity_m_s2 = 2.5",
            "planet.radius_km",
            1_499.999,
            "radius_m",
            1e3,
        ),
        (
            "radius_km = 13000.001\nsurface_gravity_m_s2 = 20",
            "planet.radius_km",
            13_000.001,
            "radius_m",
            1e3,
        ),
        (
            "surface_gravity_m_s2 = 0.999",
            "planet.surface_gravity_m_s2",
            0.999,
            "surface_gravity_m_s2",
            1.0,
        ),
        (
            "surface_gravity_m_s2 = 40.001",
            "planet.surface_gravity_m_s2",
            40.001,
            "surface_gravity_m_s2",
            1.0,
        ),
        (
            "rotation_period_h = 7.999",
            "planet.rotation_period_h",
            7.999,
            "rotation_period_s",
            3_600.0,
        ),
        (
            "rotation_period_h = 240.001",
            "planet.rotation_period_h",
            240.001,
            "rotation_period_s",
            3_600.0,
        ),
        (
            "orbital_period_days = 89.999",
            "planet.orbital_period_days",
            89.999,
            "orbital_period_s",
            86_400.0,
        ),
        (
            "orbital_period_days = 3650.001",
            "planet.orbital_period_days",
            3_650.001,
            "orbital_period_s",
            86_400.0,
        ),
        (
            "stellar_flux_w_m2 = 399.999",
            "planet.stellar_flux_w_m2",
            399.999,
            "stellar_flux_w_m2",
            1.0,
        ),
        (
            "stellar_flux_w_m2 = 2500.001",
            "planet.stellar_flux_w_m2",
            2_500.001,
            "stellar_flux_w_m2",
            1.0,
        ),
        (
            "ocean_volume_km3 = 1.000001e11",
            "planet.ocean_volume_km3",
            1.000_001e11,
            "ocean_volume_m3",
            1e9,
        ),
        (
            "surface_pressure_bar = 0.0999",
            "planet.surface_pressure_bar",
            0.0999,
            "surface_pressure_pa",
            1e5,
        ),
        (
            "surface_pressure_bar = 10.001",
            "planet.surface_pressure_bar",
            10.001,
            "surface_pressure_pa",
            1e5,
        ),
    ] {
        let text = format!("[planet]\n{line}\n");
        assert_eq!(
            err(&text),
            ConfigError::Invalid {
                key,
                written: Some(written),
                error: ParamError::OutOfRange {
                    name,
                    value: written * factor,
                },
            },
            "{text}"
        );
    }
}

#[test]
fn erosion_human_unit_keys_convert_to_si() {
    // The defaults in human units, bit for bit: 10,000 Julian years and
    // 5e-3 m²/yr are converted by the same operations as the constants.
    let (_, s) = ok("[generation.erosion]\ntime_step_yr = 10000\ndiffusivity_m2_yr = 0.005\n");
    assert_same(
        &(PlanetParams::earth(), s),
        &PlanetParams::earth(),
        &Settings::default(),
        "the defaults in years",
    );
    // Other values against Python's IEEE 754 arithmetic:
    // repr(2.5 * 31557600.0), repr(1.0 / 31557600.0),
    // repr(2e4 * 31557600.0), repr(0.01 / 31557600.0).
    for (text, time_step_s, diffusivity_m2_s) in [
        (
            "time_step_yr = 2.5\ndiffusivity_m2_yr = 1",
            78_894_000.0,
            3.168_808_781_402_895e-8,
        ),
        (
            "time_step_yr = 2e4\ndiffusivity_m2_yr = 0.01",
            631_152_000_000.0,
            3.168_808_781_402_895e-10,
        ),
        // Zero is a valid diffusivity under either key, −0 too.
        ("time_step_s = 1\ndiffusivity_m2_yr = -0.0", 1.0, -0.0),
    ] {
        let text = format!("[generation.erosion]\n{text}\n");
        let (_, s) = ok(&text);
        let expected = Settings {
            erosion: ErosionSettings {
                time_step_s,
                diffusivity_m2_s,
                ..ErosionSettings::default()
            },
            ..Settings::default()
        };
        assert_eq!(settings_bits(&s), settings_bits(&expected), "{text}");
        assert_eq!(s, expected, "{text}");
    }
}

#[test]
fn a_diffusion_step_that_overflows_names_both_keys() {
    // The keys as written, in field order, the SI key for a value left at
    // its default. 1e300 m²/yr is 3.2e292 m²/s and 1e300 yr 3.2e307 s
    // (Python: repr(1e300 / 31557600.0), repr(1e300 * 31557600.0)), each
    // finite, their product not; 1e300 m²/s times the default 3.2e11 s
    // overflows too.
    for (text, expected) in [
        (
            "diffusivity_m2_s = 1e200\ntime_step_s = 1e200",
            [
                "generation.erosion.diffusivity_m2_s",
                "generation.erosion.time_step_s",
            ],
        ),
        (
            "time_step_yr = 1e300\ndiffusivity_m2_yr = 1e300",
            [
                "generation.erosion.diffusivity_m2_yr",
                "generation.erosion.time_step_yr",
            ],
        ),
        (
            "diffusivity_m2_s = 1e300",
            [
                "generation.erosion.diffusivity_m2_s",
                "generation.erosion.time_step_s",
            ],
        ),
    ] {
        let text = format!("[generation.erosion]\n{text}\n");
        let e = err(&text);
        assert_eq!(
            e,
            ConfigError::Inconsistent {
                keys: expected,
                error: ParamError::Inconsistent {
                    name: "erosion.diffusion_step_m2",
                    value: f64::INFINITY,
                    fields: ["erosion.diffusivity_m2_s", "erosion.time_step_s"],
                },
            },
            "{text}"
        );
        assert_eq!(
            e.to_string(),
            format!(
                "`{}` and `{}` do not fit together: erosion.diffusion_step_m2 = inf is out of \
                 range",
                expected[0], expected[1]
            ),
            "{text}"
        );
    }
    // A value out of its own domain is reported first.
    assert!(matches!(
        err("[generation.erosion]\ndiffusivity_m2_s = 1e300\nerodibility = 0\n"),
        ConfigError::Invalid {
            key: "generation.erosion.erodibility",
            ..
        }
    ));
}

#[test]
fn a_radius_and_gravity_that_do_not_fit_name_both_keys() {
    // The keys as written, in field order, with Earth's gravity when the
    // file leaves it out. Densities from Python,
    // repr(3 * g / (4 * math.pi * 6.6743e-11 * r)).
    for (text, keys, density, message) in [
        (
            "[planet]\nradius_km = 13000\n",
            ["planet.radius_km", "planet.surface_gravity_m_s2"],
            2_698.254_998_805_064,
            "`planet.radius_km` and `planet.surface_gravity_m_s2` do not fit together: \
             bulk_density_kg_m3 = 2698.254998805064 is out of range",
        ),
        (
            "[planet]\nsurface_gravity_m_s2 = 1.0\nradius_m = 1.5e6\n",
            ["planet.radius_m", "planet.surface_gravity_m_s2"],
            2_384.593_786_492_896,
            "`planet.radius_m` and `planet.surface_gravity_m_s2` do not fit together: \
             bulk_density_kg_m3 = 2384.593786492896 is out of range",
        ),
        (
            "[planet]\nradius_km = 1500\nsurface_gravity_m_s2 = 7\n",
            ["planet.radius_km", "planet.surface_gravity_m_s2"],
            16_692.156_505_450_27,
            "`planet.radius_km` and `planet.surface_gravity_m_s2` do not fit together: \
             bulk_density_kg_m3 = 16692.15650545027 is out of range",
        ),
    ] {
        let e = err(text);
        match &e {
            ConfigError::Inconsistent {
                keys: k,
                error:
                    ParamError::Inconsistent {
                        name: "bulk_density_kg_m3",
                        value,
                        fields: ["radius_m", "surface_gravity_m_s2"],
                    },
            } => {
                assert_eq!(*k, keys, "{text}");
                assert!(
                    ((value - density) / density).abs() < 1e-12,
                    "{text}: {value}"
                );
            }
            e => panic!("{e:?} from\n{text}"),
        }
        assert_eq!(e.to_string(), message, "{text}");
    }
    // A field out of its own bounds is reported first.
    assert!(matches!(
        err("[planet]\nradius_km = 13000\nsurface_gravity_m_s2 = 0.5\n"),
        ConfigError::Invalid {
            key: "planet.surface_gravity_m_s2",
            ..
        }
    ));
}

#[test]
fn float_keys_accept_integers() {
    let (p, _) = ok("[planet]\nradius_m = 6371000\nocean_volume_m3 = 0\n");
    assert_eq!(p.radius_m.to_bits(), 6_371_000.0_f64.to_bits());
    assert_eq!(p.ocean_volume_m3.to_bits(), 0.0_f64.to_bits());
    // The largest TOML integer becomes the nearest float, 2^63.
    let (p, _) = ok("[planet]\nocean_volume_m3 = 9223372036854775807\n");
    assert_eq!(
        p.ocean_volume_m3.to_bits(),
        9_223_372_036_854_775_808.0_f64.to_bits()
    );
}

#[test]
fn missing_keys_take_earth_and_the_default_settings() {
    let (earth, default) = (PlanetParams::earth(), Settings::default());
    for text in [
        "",
        "\n",
        "# nothing\n",
        "[planet]\n",
        "[generation]\n",
        "[planet]\n[generation]\n",
    ] {
        assert_same(&ok(text), &earth, &default, text);
    }
    let text = "[planet]\nradius_km = 3389.5\n[generation]\ncells = 20000\n";
    assert_same(
        &ok(text),
        &PlanetParams {
            radius_m: 3_389_500.0,
            ..earth.clone()
        },
        &Settings {
            cells: 20_000,
            ..default.clone()
        },
        text,
    );
    let text = "[generation]\nplate_count = 5\nallow_partial = true\n";
    assert_same(
        &ok(text),
        &earth,
        &Settings {
            plate_count: 5,
            allow_partial: true,
            ..default
        },
        text,
    );
    // Dotted keys and inline tables are the same tables in TOML.
    let text = "planet.radius_km = 3389.5\ngeneration = { cells = 20000 }\n";
    assert_eq!(ok(text).0.radius_m.to_bits(), 3_389_500.0_f64.to_bits());
    assert_eq!(ok(text).1.cells, 20_000);
    // The erosion settings, each missing key at its default.
    for text in [
        "[generation.erosion]\n",
        "[generation]\n[generation.erosion]\n",
        "generation.erosion = {}\n",
    ] {
        assert_same(&ok(text), &earth, &Settings::default(), text);
    }
    let some = Settings {
        cells: 2_000,
        erosion: ErosionSettings {
            discharge_exponent: 0.5,
            steps: 5,
            ..ErosionSettings::default()
        },
        ..Settings::default()
    };
    for text in [
        "[generation]\ncells = 2000\n[generation.erosion]\ndischarge_exponent = 0.5\nsteps = 5\n",
        "[generation.erosion]\nsteps = 5\ndischarge_exponent = 0.5\n[generation]\ncells = 2000\n",
        "generation.cells = 2000\ngeneration.erosion.steps = 5\n\
         generation.erosion.discharge_exponent = 0.5\n",
        "[generation]\ncells = 2000\nerosion = { steps = 5, discharge_exponent = 0.5 }\n",
    ] {
        assert_same(&ok(text), &earth, &some, text);
    }
}

#[test]
fn regimes_are_read_by_name() {
    let (p, _) = ok("[planet]\nrotation = \"free\"\ntectonics = \"plate_tectonics\"\n");
    assert_eq!(
        (p.rotation, p.tectonics),
        (RotationRegime::Free, TectonicRegime::PlateTectonics)
    );
}

#[test]
fn both_keys_of_a_quantity_are_an_error() {
    for (table, si, human) in [
        ("planet", "radius_m", "radius_km"),
        ("planet", "rotation_period_s", "rotation_period_h"),
        ("planet", "obliquity_rad", "obliquity_deg"),
        ("planet", "orbital_period_s", "orbital_period_days"),
        ("planet", "ocean_volume_m3", "ocean_volume_km3"),
        ("planet", "surface_pressure_pa", "surface_pressure_bar"),
        (
            "generation.erosion",
            "diffusivity_m2_s",
            "diffusivity_m2_yr",
        ),
        ("generation.erosion", "time_step_s", "time_step_yr"),
    ] {
        // Even with the same value under both, and in either order.
        for text in [
            format!("[{table}]\n{si} = 1.0\n{human} = 1.0\n"),
            format!("[{table}]\n{human} = 0.001\n{si} = 1.0\n"),
        ] {
            let e = err(&text);
            match &e {
                ConfigError::BothKeys { si_key, human_key } => {
                    assert_eq!(*si_key, format!("{table}.{si}"), "{text}");
                    assert_eq!(*human_key, format!("{table}.{human}"), "{text}");
                }
                e => panic!("{e:?} from\n{text}"),
            }
            assert_eq!(
                e.to_string(),
                format!(
                    "`{table}.{si}` and `{table}.{human}` are the same quantity; give one of them"
                )
            );
        }
    }
}

#[test]
fn unknown_tables_and_keys_are_errors_with_their_position() {
    for (text, name, position) in [
        (
            "[planet]\nradius = 6371000.0\n",
            "radius",
            "line 2, column 1",
        ),
        (
            "[planet]\nobliquity_degrees = 23.44\n",
            "obliquity_degrees",
            "line 2, column 1",
        ),
        (
            "[planet]\nsurface_gravity_g = 1.0\n",
            "surface_gravity_g",
            "line 2, column 1",
        ),
        ("[generation]\nseed = 1\n", "seed", "line 2, column 1"),
        ("seed = 1\n", "seed", "line 1, column 1"),
        ("[orbit]\nperiod_s = 1.0\n", "orbit", "line 1, column 2"),
        (
            "[planet]\nradius_m = 1.0\n[planets]\n",
            "planets",
            "line 3, column 2",
        ),
        ("[generation.erosion]\nk = 1.0\n", "k", "line 2, column 1"),
        (
            "[generation.erosion]\ntime_step_kyr = 10\n",
            "time_step_kyr",
            "line 2, column 1",
        ),
        // The erosion settings are a table of `[generation]` only.
        ("[erosion]\nsteps = 1\n", "erosion", "line 1, column 2"),
        (
            "[planet.erosion]\nsteps = 1\n",
            "erosion",
            "line 1, column 9",
        ),
    ] {
        let message = toml_message(text);
        assert!(
            message.contains(&format!("unknown field `{name}`")),
            "{text}\n{message}"
        );
        assert!(message.contains(position), "{text}\n{message}");
    }
}

#[test]
fn wrong_types_are_errors_with_their_position() {
    for (text, expected, position) in [
        (
            "[planet]\nradius_m = \"big\"\n",
            "invalid type: string \"big\", expected f64",
            "line 2, column 12",
        ),
        (
            "[planet]\nradius_km = true\n",
            "invalid type: boolean `true`, expected f64",
            "line 2, column 13",
        ),
        (
            "[generation]\ncells = 1.5e4\n",
            "invalid type: floating point `15000.0`, expected u32",
            "line 2, column 9",
        ),
        (
            "[generation]\ncells = -1\n",
            "invalid value: integer `-1`, expected u32",
            "line 2, column 9",
        ),
        (
            "[generation]\nplate_count = 4294967296\n",
            "invalid value: integer `4294967296`, expected u32",
            "line 2, column 15",
        ),
        (
            "[generation]\nallow_partial = 1\n",
            "invalid type: integer `1`, expected a boolean",
            "line 2, column 17",
        ),
        (
            "[planet]\nrotation = \"spinning\"\n",
            "unknown variant `spinning`, expected `free` or `tidally_locked`",
            "line 2, column 12",
        ),
        (
            "[planet]\ntectonics = \"plates\"\n",
            "unknown variant `plates`, expected `plate_tectonics` or `stagnant_lid`",
            "line 2, column 13",
        ),
        (
            "[planet]\nrotation = 1\n",
            "wanted string or table",
            "line 2, column 12",
        ),
        (
            "planet = 1\n",
            "invalid type: integer `1`, expected the table [planet]",
            "line 1, column 10",
        ),
        (
            "generation = \"fast\"\n",
            "expected the table [generation]",
            "line 1, column 14",
        ),
        (
            "[generation.erosion]\nsteps = -1\n",
            "invalid value: integer `-1`, expected u32",
            "line 2, column 9",
        ),
        (
            "[generation.erosion]\nsteps = 1.5\n",
            "invalid type: floating point `1.5`, expected u32",
            "line 2, column 9",
        ),
        (
            "[generation.erosion]\nerodibility = \"small\"\n",
            "invalid type: string \"small\", expected f64",
            "line 2, column 15",
        ),
        (
            "[generation]\nerosion = 1\n",
            "invalid type: integer `1`, expected the table [generation.erosion]",
            "line 2, column 11",
        ),
    ] {
        let message = toml_message(text);
        assert!(message.contains(expected), "{text}\n{message}");
        assert!(message.contains(position), "{text}\n{message}");
    }
}

#[test]
fn malformed_toml_is_an_error_with_its_position() {
    for (text, position) in [
        ("[planet\nradius_m = 1.0\n", "line 1, column 8"),
        ("[planet]\nradius_m = \n", "line 2, column 12"),
        ("[planet]\nradius_m = 6,371\n", "line 2, column 13"),
        ("[planet]\nradius_m 6371000.0\n", "line 2, column 10"),
        (
            "[planet]\nradius_m = 1.0\nradius_m = 2.0\n",
            "line 3, column 1",
        ),
        // The duplicate key is the table's name.
        ("[planet]\n[planet]\n", "line 2, column 2"),
        // At the end of the line, where the closing quote is missing.
        ("[planet]\nrotation = \"free\n", "line 2, column 17"),
    ] {
        let message = toml_message(text);
        assert!(
            message.starts_with("TOML parse error at "),
            "{text}\n{message}"
        );
        assert!(message.contains(position), "{text}\n{message}");
    }
}

#[test]
fn out_of_range_values_name_the_key_as_written() {
    // The text, the key, the number written under it, and the SI name and
    // value that `validate` rejects.
    for (text, key, written, name, value) in [
        (
            "[planet]\nradius_m = 0.0\n",
            "planet.radius_m",
            0.0,
            "radius_m",
            0.0,
        ),
        (
            "[planet]\nradius_km = -1\n",
            "planet.radius_km",
            -1.0,
            "radius_m",
            -1_000.0,
        ),
        // 1e306 km is beyond the largest float in metres.
        (
            "[planet]\nradius_km = 1e306\n",
            "planet.radius_km",
            1e306,
            "radius_m",
            f64::INFINITY,
        ),
        (
            "[planet]\nsurface_gravity_m_s2 = -inf\n",
            "planet.surface_gravity_m_s2",
            f64::NEG_INFINITY,
            "surface_gravity_m_s2",
            f64::NEG_INFINITY,
        ),
        (
            "[planet]\nrotation_period_h = 0\n",
            "planet.rotation_period_h",
            0.0,
            "rotation_period_s",
            0.0,
        ),
        (
            "[planet]\nobliquity_deg = 180.1\n",
            "planet.obliquity_deg",
            180.1,
            "obliquity_rad",
            180.1 * (core::f64::consts::PI / 180.0),
        ),
        (
            "[planet]\nobliquity_rad = -1e-300\n",
            "planet.obliquity_rad",
            -1e-300,
            "obliquity_rad",
            -1e-300,
        ),
        (
            "[planet]\norbital_period_days = -365\n",
            "planet.orbital_period_days",
            -365.0,
            "orbital_period_s",
            -365.0 * 86_400.0,
        ),
        (
            "[planet]\norbital_eccentricity = 0.995\n",
            "planet.orbital_eccentricity",
            0.995,
            "orbital_eccentricity",
            0.995,
        ),
        (
            "[planet]\nstellar_flux_w_m2 = inf\n",
            "planet.stellar_flux_w_m2",
            f64::INFINITY,
            "stellar_flux_w_m2",
            f64::INFINITY,
        ),
        (
            "[planet]\nocean_volume_km3 = -1\n",
            "planet.ocean_volume_km3",
            -1.0,
            "ocean_volume_m3",
            -1e9,
        ),
        (
            "[planet]\nsurface_pressure_bar = 0\n",
            "planet.surface_pressure_bar",
            0.0,
            "surface_pressure_pa",
            0.0,
        ),
        (
            "[planet]\ncontinental_crust_fraction = 1.5\n",
            "planet.continental_crust_fraction",
            1.5,
            "continental_crust_fraction",
            1.5,
        ),
        (
            "[generation]\ncells = 999\n",
            "generation.cells",
            999.0,
            "cells",
            999.0,
        ),
        (
            "[generation]\ncells = 715827885\n",
            "generation.cells",
            715_827_885.0,
            "cells",
            715_827_885.0,
        ),
        (
            "[generation]\nplate_count = 1\n",
            "generation.plate_count",
            1.0,
            "plate_count",
            1.0,
        ),
        (
            "[generation]\nplate_count = 65\n",
            "generation.plate_count",
            65.0,
            "plate_count",
            65.0,
        ),
        // The erosion settings, named by their path from `Settings`; the
        // SI values of the human-unit keys from Python,
        // repr(-1.0 / 31557600.0) and repr(1e301 * 31557600.0).
        (
            "[generation.erosion]\nerodibility = 0\n",
            "generation.erosion.erodibility",
            0.0,
            "erosion.erodibility",
            0.0,
        ),
        (
            "[generation.erosion]\nerodibility = -inf\n",
            "generation.erosion.erodibility",
            f64::NEG_INFINITY,
            "erosion.erodibility",
            f64::NEG_INFINITY,
        ),
        (
            "[generation.erosion]\ndischarge_exponent = -0.1\n",
            "generation.erosion.discharge_exponent",
            -0.1,
            "erosion.discharge_exponent",
            -0.1,
        ),
        (
            "[generation.erosion]\ndischarge_exponent = inf\n",
            "generation.erosion.discharge_exponent",
            f64::INFINITY,
            "erosion.discharge_exponent",
            f64::INFINITY,
        ),
        (
            "[generation.erosion]\ndiffusivity_m2_s = -1e-300\n",
            "generation.erosion.diffusivity_m2_s",
            -1e-300,
            "erosion.diffusivity_m2_s",
            -1e-300,
        ),
        (
            "[generation.erosion]\ndiffusivity_m2_yr = -1\n",
            "generation.erosion.diffusivity_m2_yr",
            -1.0,
            "erosion.diffusivity_m2_s",
            -3.168_808_781_402_895e-8,
        ),
        (
            "[generation.erosion]\ndiffusivity_m2_yr = inf\n",
            "generation.erosion.diffusivity_m2_yr",
            f64::INFINITY,
            "erosion.diffusivity_m2_s",
            f64::INFINITY,
        ),
        (
            "[generation.erosion]\ntime_step_s = 0\n",
            "generation.erosion.time_step_s",
            0.0,
            "erosion.time_step_s",
            0.0,
        ),
        (
            "[generation.erosion]\ntime_step_yr = -1\n",
            "generation.erosion.time_step_yr",
            -1.0,
            "erosion.time_step_s",
            -31_557_600.0,
        ),
        // 1e301 years is beyond the largest float in seconds.
        (
            "[generation.erosion]\ntime_step_yr = 1e301\n",
            "generation.erosion.time_step_yr",
            1e301,
            "erosion.time_step_s",
            f64::INFINITY,
        ),
    ] {
        assert_eq!(
            err(text),
            ConfigError::Invalid {
                key,
                written: Some(written),
                error: ParamError::OutOfRange { name, value },
            },
            "{text}"
        );
    }
    // NaN, which is never equal to itself.
    for (text, key) in [
        ("[planet]\nradius_m = nan\n", "planet.radius_m"),
        ("[planet]\nobliquity_deg = -nan\n", "planet.obliquity_deg"),
        (
            "[generation.erosion]\nerodibility = nan\n",
            "generation.erosion.erodibility",
        ),
        (
            "[generation.erosion]\ntime_step_yr = nan\n",
            "generation.erosion.time_step_yr",
        ),
    ] {
        match err(text) {
            ConfigError::Invalid {
                key: k,
                written: Some(written),
                error: ParamError::OutOfRange { value, .. },
            } => {
                assert_eq!(k, key, "{text}");
                assert!(written.is_nan() && value.is_nan(), "{text}");
            }
            e => panic!("{e:?} from\n{text}"),
        }
    }
    // The planet is checked before the settings, each in field order, the
    // erosion settings after the cells and the plate count.
    for (text, key) in [
        (
            "[generation]\ncells = 1\n[planet]\nradius_km = 0\nobliquity_deg = 200\n",
            "planet.radius_km",
        ),
        (
            "[generation.erosion]\nerodibility = 0\n[planet]\nradius_km = 0\n",
            "planet.radius_km",
        ),
        (
            "[generation.erosion]\nerodibility = 0\n[generation]\nplate_count = 1\n",
            "generation.plate_count",
        ),
        (
            "[generation.erosion]\ntime_step_yr = -1\ndiffusivity_m2_yr = -1\n\
             discharge_exponent = -1\n",
            "generation.erosion.discharge_exponent",
        ),
        (
            "[generation.erosion]\ntime_step_yr = -1\ndiffusivity_m2_yr = -1\n",
            "generation.erosion.diffusivity_m2_yr",
        ),
    ] {
        let e = err(text);
        assert!(
            matches!(e, ConfigError::Invalid { key: k, .. } if k == key),
            "{text}: {e:?}"
        );
    }
}

#[test]
fn out_of_range_messages_give_the_number_as_written() {
    // A human-unit key adds the SI value that was checked:
    // repr(200 * (math.pi / 180)) in Python.
    for (text, message) in [
        (
            "[planet]\nobliquity_deg = 200\n",
            "`planet.obliquity_deg` = 200 is out of range (obliquity_rad = 3.490658503988659)",
        ),
        (
            "[planet]\nradius_km = -1\n",
            "`planet.radius_km` = -1 is out of range (radius_m = -1000)",
        ),
        (
            "[planet]\nradius_m = 0\n",
            "`planet.radius_m` = 0 is out of range",
        ),
        (
            "[planet]\nradius_m = nan\n",
            "`planet.radius_m` = NaN is out of range",
        ),
        (
            "[generation]\ncells = 999\n",
            "`generation.cells` = 999 is out of range",
        ),
        // Exponents outside 1e-5 to 1e16, so that no message runs to
        // hundreds of digits; from Python, repr(-5e-6 * 1e5) and
        // repr(-1e8 * 1e9).
        (
            "[planet]\nradius_km = 1e306\n",
            "`planet.radius_km` = 1e306 is out of range (radius_m = inf)",
        ),
        (
            "[planet]\nobliquity_rad = -1e-300\n",
            "`planet.obliquity_rad` = -1e-300 is out of range",
        ),
        (
            "[planet]\nsurface_pressure_bar = -5e-6\n",
            "`planet.surface_pressure_bar` = -5e-6 is out of range (surface_pressure_pa = -0.5)",
        ),
        (
            "[planet]\nocean_volume_km3 = -1e8\n",
            "`planet.ocean_volume_km3` = -100000000 is out of range (ocean_volume_m3 = -1e17)",
        ),
        // The erosion settings: an SI key gives the number alone, a
        // human-unit key adds the SI value under its path from `Settings`
        // (repr(-1.0 / 31557600.0) in Python).
        (
            "[generation.erosion]\nerodibility = 0\n",
            "`generation.erosion.erodibility` = 0 is out of range",
        ),
        (
            "[generation.erosion]\ntime_step_s = -1\n",
            "`generation.erosion.time_step_s` = -1 is out of range",
        ),
        (
            "[generation.erosion]\ntime_step_yr = -1\n",
            "`generation.erosion.time_step_yr` = -1 is out of range \
             (erosion.time_step_s = -31557600)",
        ),
        (
            "[generation.erosion]\ndiffusivity_m2_yr = -1\n",
            "`generation.erosion.diffusivity_m2_yr` = -1 is out of range \
             (erosion.diffusivity_m2_s = -3.168808781402895e-8)",
        ),
        (
            "[generation.erosion]\ntime_step_yr = 1e301\n",
            "`generation.erosion.time_step_yr` = 1e301 is out of range \
             (erosion.time_step_s = inf)",
        ),
    ] {
        assert_eq!(err(text).to_string(), message, "{text}");
    }
}

#[test]
fn v2_regimes_are_unsupported_errors_naming_the_key() {
    for (text, key) in [
        (
            "[planet]\nrotation = \"tidally_locked\"\n",
            "planet.rotation",
        ),
        (
            "[planet]\ntectonics = \"stagnant_lid\"\n",
            "planet.tectonics",
        ),
        // Rotation is checked first.
        (
            "[planet]\ntectonics = \"stagnant_lid\"\nrotation = \"tidally_locked\"\n",
            "planet.rotation",
        ),
    ] {
        let e = err(text);
        match &e {
            ConfigError::Invalid {
                key: k,
                written: None,
                error: ParamError::Unsupported { planned: "v2", .. },
            } => assert_eq!(*k, key, "{text}"),
            e => panic!("{e:?} from\n{text}"),
        }
        assert!(e.to_string().starts_with(&format!("`{key}`: ")), "{e}");
        assert!(
            e.to_string()
                .ends_with("is not supported yet (planned for v2)"),
            "{e}"
        );
    }
}

#[test]
fn toml_errors_display_the_report_of_the_toml_crate() {
    let text = "[planet]\nradius = 1.0\n";
    let e = err(text);
    let message = toml_message(text);
    assert_eq!(e.to_string(), message);
    // The report shows the line, with the key marked.
    assert!(message.contains("2 | radius = 1.0"), "{message}");
    assert!(!message.ends_with('\n'), "{message:?}");
}

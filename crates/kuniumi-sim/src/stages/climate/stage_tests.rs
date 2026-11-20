//! Tests of the Climate stage (tickets M4-07 and M4-17):
//!
//! - `the_field_names_are_the_ticket_s`;
//! - `the_stage_composes_the_models`: every field bit for bit against the
//!   models of tickets M4-02 to M4-06 and M4-17 composed by hand, the annual means
//!   recomputed from the months, the names, units, kinds and order;
//! - `the_climate_is_earth_like`: finite, positive temperatures,
//!   non-negative precipitation, still water on land, Earth's global means
//!   and seasons;
//! - `climate_reports_its_progress`, `existing_outputs_leave_the_world_unchanged`,
//!   `the_stage_needs_the_heights`;
//! - `a_climate_that_settles_has_no_warnings` and
//!   `warnings_name_each_model_that_did_not_settle` (ticket M4-16; the
//!   warnings of the planets that do not settle are checked by the bounds
//!   tests of `lib.rs`);
//! - `golden_fingerprints_of_the_climate_stage`.
//!
//! The stage runs once per test process on the Earth-parameter world of
//! seed 1 on 10,000 cells through the Crust stage ([`earth_10k`]), in debug
//! and release builds alike; the composition by hand reuses the runs of
//! [`earth_10k`] and [`earth_coupled_10k`] that the models' own tests share.
#![allow(clippy::disallowed_methods)]

use core::f64::consts::PI;
use std::sync::OnceLock;

use kuniumi_math::rng::fnv1a64;
use kuniumi_math::units::SECONDS_PER_YEAR;
use kuniumi_world::{FieldKind, PlanetParams, Settings, World};

use super::moisture::monthly_moisture;
use super::test_worlds::{earth_10k, earth_coupled_10k};
use super::*;

/// The world of [`earth_10k`] after the Climate stage, and the stage's
/// progress reports, computed once.
fn stage_world() -> &'static (World, Vec<f64>) {
    static WORLD: OnceLock<(World, Vec<f64>)> = OnceLock::new();
    WORLD.get_or_init(|| {
        let mut world = earth_10k().world.clone();
        let mut reports = Vec::new();
        run(&mut world, &mut |f| reports.push(f)).expect("the stage runs");
        (world, reports)
    })
}

/// The values of the field `name` of `world`.
fn values<'a>(world: &'a World, name: &str) -> &'a [f64] {
    &world
        .field(name)
        .unwrap_or_else(|| panic!("no field `{name}`"))
        .values
}

/// The bits of `values`, to compare them exactly.
fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|v| v.to_bits()).collect()
}

/// The mean of twelve monthly values of one cell, as the stage defines it:
/// the sum in month order, divided by 12.
fn mean_of_months(months: &[Vec<f64>; MONTHS], i: usize) -> f64 {
    let mut sum = 0.0;
    for month in months {
        sum += month[i];
    }
    sum / 12.0
}

/// The 30 fields of the stage in their order, with their units.
fn expected_names() -> Vec<(&'static str, &'static str)> {
    let mut names = vec![(TEMPERATURE_ANNUAL, "K")];
    names.extend(TEMPERATURE_MONTHS.iter().map(|&n| (n, "K")));
    names.push((PRECIPITATION_ANNUAL, "m/s"));
    names.extend(PRECIPITATION_MONTHS.iter().map(|&n| (n, "m/s")));
    names.extend([
        (WIND_EAST, "m/s"),
        (WIND_NORTH, "m/s"),
        (CURRENT_EAST, "m/s"),
        (CURRENT_NORTH, "m/s"),
    ]);
    names
}

#[test]
fn the_field_names_are_the_ticket_s() {
    // The monthly names of the ticket, spelt out.
    for k in 0..MONTHS {
        assert_eq!(TEMPERATURE_MONTHS[k], format!("temperature_m{:02}", k + 1));
        assert_eq!(
            PRECIPITATION_MONTHS[k],
            format!("precipitation_m{:02}", k + 1)
        );
    }
    assert_eq!(expected_names().len(), 30);
}

#[test]
fn the_stage_composes_the_models() {
    let (world, _) = stage_world();
    let before = &earth_10k().world;
    let params = PlanetParams::earth();
    let grid = &before.grid;
    let heights = values(before, SURFACE_HEIGHT);
    let z = before.sea_level_m;

    // The fields of the Crust stage and before, the sea level and the
    // stages: unchanged, with `climate` added.
    assert_eq!(&world.fields[..before.fields.len()], &before.fields[..]);
    assert_eq!(world.sea_level_m.to_bits(), z.to_bits());
    assert_eq!(world.stages, before.stages, "the runner records the stage");
    let added: Vec<_> = world.fields[before.fields.len()..]
        .iter()
        .map(|f| (f.name.as_str(), f.unit.as_str(), f.kind))
        .collect();
    let expected: Vec<_> = expected_names()
        .into_iter()
        .map(|(n, u)| (n, u, FieldKind::Quantity))
        .collect();
    assert_eq!(added, expected);

    // The models composed by hand, as the tests of the coupling compose
    // them: the energy balance's periodic state, the coupled state of the
    // sea surface, winds and currents from it, the moisture of its winds
    // over its temperatures.
    let coupled = earth_coupled_10k();
    let (circulation, currents) = (&coupled.circulation, &coupled.currents);
    let sea = &coupled.sea_surface.temperature;
    let moisture = monthly_moisture(&params, grid, heights, z, circulation, sea);

    let n = grid.len();
    for k in 0..MONTHS {
        assert!(
            bits(values(world, TEMPERATURE_MONTHS[k])) == bits(&sea.surface_k[k]),
            "{}",
            TEMPERATURE_MONTHS[k]
        );
        assert!(
            bits(values(world, PRECIPITATION_MONTHS[k])) == bits(&moisture.precipitation_m_s[k]),
            "{}",
            PRECIPITATION_MONTHS[k]
        );
    }
    let annual = |months: &[Vec<f64>; MONTHS]| -> Vec<f64> {
        (0..n).map(|i| mean_of_months(months, i)).collect()
    };
    assert!(bits(values(world, TEMPERATURE_ANNUAL)) == bits(&annual(&sea.surface_k)));
    assert!(
        bits(values(world, PRECIPITATION_ANNUAL)) == bits(&annual(&moisture.precipitation_m_s))
    );
    assert!(bits(values(world, WIND_EAST)) == bits(&annual(&circulation.wind_east_m_s)));
    assert!(bits(values(world, WIND_NORTH)) == bits(&annual(&circulation.wind_north_m_s)));
    assert!(bits(values(world, CURRENT_EAST)) == bits(&currents.current_east_m_s));
    assert!(bits(values(world, CURRENT_NORTH)) == bits(&currents.current_north_m_s));
    // Over the ocean the surface temperature is the sea-surface
    // temperature, bit for bit (ticket M4-05): no field of its own.
    for (surface, sea_level) in sea.surface_k.iter().zip(&sea.sea_level_k) {
        for ((s, l), h) in surface.iter().zip(sea_level).zip(heights) {
            if *h <= z {
                assert_eq!(s.to_bits(), l.to_bits());
            }
        }
    }
}

/// The area-weighted mean of `values` over the cells whose sine of latitude
/// passes `keep`.
fn area_mean(world: &World, values: &[f64], keep: impl Fn(f64) -> bool) -> f64 {
    let (mut sum, mut area) = (0.0, 0.0);
    for (i, v) in values.iter().enumerate() {
        let a = world.grid.area_sr(i as u32);
        if keep(world.grid.position(i as u32).z) {
            sum += a * v;
            area += a;
        }
    }
    sum / area
}

#[test]
fn the_climate_is_earth_like() {
    let (world, _) = stage_world();
    let heights = values(world, SURFACE_HEIGHT);
    let z = world.sea_level_m;
    for (name, _) in expected_names() {
        let v = values(world, name);
        assert!(v.iter().all(|x| x.is_finite()), "{name}");
        if name.starts_with("temperature") {
            assert!(v.iter().all(|&t| t > 0.0), "{name}");
        }
        if name.starts_with("precipitation") {
            assert!(v.iter().all(|&p| p >= 0.0), "{name}");
        }
    }
    for name in [CURRENT_EAST, CURRENT_NORTH] {
        let v = values(world, name);
        for (i, (c, h)) in v.iter().zip(heights).enumerate() {
            if *h > z {
                assert_eq!(c.to_bits(), 0, "{name} on land, cell {i}");
            }
        }
    }
    // Earth's global means: the calibration's 288 K (M4-05 and M4-17,
    // 40,000 cells; the coupled state of this 10,000-cell world 289.02 K in
    // the prototype of ticket M4-17, the chain of M4-07 288.98 K) and about
    // 1 m/yr (M4-06; the prototype 0.958 m/yr on this world).
    let mean_k = area_mean(world, values(world, TEMPERATURE_ANNUAL), |_| true);
    let rain_m_yr =
        area_mean(world, values(world, PRECIPITATION_ANNUAL), |_| true) * SECONDS_PER_YEAR;
    println!("mean surface temperature {mean_k:.3} K, precipitation {rain_m_yr:.4} m/yr");
    assert!((mean_k - 288.0).abs() <= 1.5, "{mean_k}");
    assert!((0.8..=1.2).contains(&rain_m_yr), "{rain_m_yr}");
    // The seasons: the southern summer is January, the northern July.
    let s30 = -(PI / 6.0).sin();
    let n30 = (PI / 6.0).sin();
    let south = |m: &str| area_mean(world, values(world, m), |sin_lat| sin_lat < s30);
    let north = |m: &str| area_mean(world, values(world, m), |sin_lat| sin_lat > n30);
    let (jan, jul) = (TEMPERATURE_MONTHS[0], TEMPERATURE_MONTHS[6]);
    println!(
        "south of 30° S: January {:.2} K, July {:.2} K; north of 30° N: January {:.2} K, July {:.2} K",
        south(jan),
        south(jul),
        north(jan),
        north(jul)
    );
    assert!(south(jan) > south(jul));
    assert!(north(jan) < north(jul));
}

/// A 1,000-cell Earth world through the Crust stage.
fn small_world() -> World {
    let settings = Settings {
        cells: 1_000,
        allow_partial: true,
        ..Settings::default()
    };
    crate::generate_through_crust(&PlanetParams::earth(), &settings, 1).expect("a partial world")
}

/// The climate of [`small_world`] and the progress that [`climate`]
/// reported, computed once.
fn small_climate() -> &'static (Climate, Vec<f64>) {
    static CLIMATE: OnceLock<(Climate, Vec<f64>)> = OnceLock::new();
    CLIMATE.get_or_init(|| {
        let world = small_world();
        let mut reports = Vec::new();
        let result = climate(
            &world.params,
            &world.grid,
            values(&world, SURFACE_HEIGHT),
            world.sea_level_m,
            &mut |f| reports.push(f),
        );
        (result, reports)
    })
}

#[test]
fn climate_reports_its_progress() {
    let (result, reports) = small_climate();
    assert_eq!(reports, &PROGRESS);
    assert!(result.energy_balance.converged && result.sea_surface.temperature.converged);
    // The stage passes them on.
    let (_, stage_reports) = stage_world();
    assert_eq!(stage_reports, &PROGRESS);
}

#[test]
fn a_climate_that_settles_has_no_warnings() {
    // Earth on 1,000 cells (seed 1): every model periodic, every solve
    // converged (the scratch tool of ticket M4-16).
    let (result, _) = small_climate();
    assert!(result.energy_balance.converged);
    assert!(result.sea_surface.temperature.converged);
    assert!(result.moisture.converged);
    assert_eq!(
        (
            result.moisture.unsettled_steps,
            result.moisture.failed_solves
        ),
        (0, 0)
    );
    assert!(result.currents.solves.iter().all(|s| s.converged()));
    assert_eq!(result.warnings(), Vec::<String>::new());
    // Nor does the stage add any to the 10,000-cell Earth world.
    let (world, _) = stage_world();
    assert_eq!(world.warnings, []);
}

/// The six messages of [`Climate::warnings`] for the flags that
/// `warnings_name_each_model_that_did_not_settle` sets, written out by hand
/// from the ticket (M4-16, "The Climate stage's warnings"; the energy
/// balance's and the sea surface's as ticket M4-17 renamed them, the
/// moisture's solves with the balance steps of ticket M4-22), in their
/// order.
const UNSETTLED_MESSAGES: [&str; 6] = [
    "the energy balance without currents did not reach its periodic state in 100 years: a \
     monthly-mean sea-level temperature changed by 11.500 K in the last year, whose \
     temperatures are kept as the start of the coupled sea surface; they may be those of a \
     climate in transition (an ice front still advancing on a planet with little heat \
     transport, or a planet near an ice-albedo bifurcation)",
    "the currents of 2 of the 3 ocean basins did not reach the solver's tolerance (the \
     largest relative residual 3.200e-9): those currents, and the sea surface and moisture \
     made with them, may be inaccurate",
    "the sea surface coupled with its winds and currents did not reach its periodic state in \
     100 years: a monthly-mean sea-level temperature changed by 7.405 K in the last year, \
     whose temperatures, winds and currents are kept as the climate's, with the moisture made \
     from them; they may be those of a climate in transition (an ice front still advancing on \
     a planet with little heat transport, a planet near an ice-albedo bifurcation, or winds, \
     upwelling and currents that keep changing each other)",
    "the moisture did not reach its periodic state in 10 years: a monthly-mean column water \
     vapour changed by 2.006e-6 kg/m^2 in the last year, whose precipitation is kept as the \
     climate's; it may be that of an atmosphere still filling or drying",
    "3 of the moisture's 480 sub-steps did not settle the active set within 50 solves: each \
     kept the set of its last solve, in which a cell below the threshold rains a negative \
     excess",
    "2 of the moisture's monthly potentials, sub-steps and balance steps had a solve that did \
     not reach its tolerance: its precipitation may be inaccurate",
];

#[test]
fn warnings_name_each_model_that_did_not_settle() {
    use crate::operators::Solve;

    let (settled, _) = small_climate();
    // Each condition of the ticket, set by hand on the settled climate:
    // the energy balance after 100 years (Earth at 1,296 W/m², ticket
    // M4-02: 11.5 K), two of three basins at the solver's limit (the
    // larger failing residual second, the converged basin's larger than
    // both: the message takes the largest of the failing ones), the sea
    // surface after 100 years (corner 57 of the bounds tests: 7.4052 K),
    // the moisture after 10 years (corner 17: 2.006e-6 kg/m²), 3 sub-steps
    // unsettled and 2 solves failed.
    let edits: [fn(&mut Climate); 6] = [
        |c| {
            c.energy_balance.years = 100;
            c.energy_balance.change_k = 11.5;
            c.energy_balance.converged = false;
        },
        |c| {
            c.currents.solves = vec![
                Solve::new(12, 9e-9, true),
                Solve::new(20_000, 1.5e-9, false),
                Solve::new(20_000, 3.2e-9, false),
            ];
        },
        |c| {
            let t = &mut c.sea_surface.temperature;
            t.years = 100;
            t.change_k = 7.4052;
            t.converged = false;
        },
        |c| {
            c.moisture.years = 10;
            c.moisture.change_kg_m2 = 2.006e-6;
            c.moisture.converged = false;
        },
        |c| {
            // Of 10 years, alone too: 480 sub-steps.
            c.moisture.years = 10;
            c.moisture.unsettled_steps = 3;
        },
        |c| c.moisture.failed_solves = 2,
    ];
    // Each alone gives its message alone.
    for (edit, message) in edits.iter().zip(UNSETTLED_MESSAGES) {
        let mut climate = settled.clone();
        edit(&mut climate);
        assert_eq!(climate.warnings(), [message]);
    }
    // All together, the six in the order the models ran.
    let mut climate = settled.clone();
    for edit in edits {
        edit(&mut climate);
    }
    assert_eq!(climate.warnings(), UNSETTLED_MESSAGES);
}

#[test]
fn existing_outputs_leave_the_world_unchanged() {
    for (name, _) in expected_names() {
        let mut world = small_world();
        world
            .add_field(name, "", vec![0.0; world.grid.len()])
            .expect("a new field");
        let before = world.clone();
        let mut reports = 0;
        let result = run(&mut world, &mut |_| reports += 1);
        assert_eq!(
            result,
            Err(SimError::Field(FieldError::Duplicate { name: name.into() }))
        );
        // The fields, the stages and the warnings (none added).
        assert_eq!(world, before, "{name}");
        assert_eq!(reports, 0, "{name}");
    }
}

#[test]
#[should_panic(expected = "the Climate stage needs the field `surface_height` of the Crust stage")]
fn the_stage_needs_the_heights() {
    let mut world = small_world();
    world.fields.retain(|f| f.name != SURFACE_HEIGHT);
    let _ = run(&mut world, &mut |_| {});
}

#[test]
fn golden_fingerprints_of_the_climate_stage() {
    // The Earth world of seed 1 at 10,000 cells after the Climate stage:
    // FNV-1a over the bits of each field's values (f64, little-endian). A
    // change means that the same world has another climate. Recorded by
    // ticket M4-07 after the other tests of this file passed and the maps
    // of out/M4-07 were viewed, and again by ticket M4-15 (the balance step
    // after each year: Earth's temperatures within the periodic tolerance
    // of their periodic state, world-breaking) after the other tests passed
    // and the maps of out/M4-15 were viewed, and by ticket M4-17 (the
    // coupled state of the sea surface, winds and currents, world-breaking)
    // after the other tests passed and the maps of out/M4-17 were viewed,
    // and by ticket M4-21 (the thermal equator as the warm centre,
    // world-breaking) after the other tests passed and the maps of
    // out/M4-21 were viewed; `kuniumi generate` writes the same values (the
    // CLI's golden test). The
    // temperatures, winds and currents are those of the coupled state of
    // `coupling::tests` on this world.
    let (world, _) = stage_world();
    let fingerprints: Vec<(&str, u64)> = GOLDEN_FINGERPRINTS
        .iter()
        .map(|&(name, _)| {
            let bytes: Vec<u8> = values(world, name)
                .iter()
                .flat_map(|v| v.to_bits().to_le_bytes())
                .collect();
            (name, fnv1a64(&bytes))
        })
        .collect();
    for (name, fingerprint) in &fingerprints {
        println!("    (\"{name}\", {fingerprint:#018X}),");
    }
    let names: Vec<&str> = expected_names().into_iter().map(|(n, _)| n).collect();
    let golden: Vec<&str> = GOLDEN_FINGERPRINTS.iter().map(|&(n, _)| n).collect();
    assert_eq!(golden, names);
    assert_eq!(fingerprints, GOLDEN_FINGERPRINTS);
}

/// The fingerprints of `golden_fingerprints_of_the_climate_stage`.
/// Recorded again in ticket M4-10 (the zenith angle law of the albedo,
/// world-breaking) after the other tests passed and the images of out/M4-10
/// were viewed.
const GOLDEN_FINGERPRINTS: [(&str, u64); 30] = [
    ("temperature_annual", 0x5797_419B_B591_6559),
    ("temperature_m01", 0x36F6_6608_2320_0DD2),
    ("temperature_m02", 0x388E_DED5_9021_BCC9),
    ("temperature_m03", 0xE0BA_4185_3E6F_C364),
    ("temperature_m04", 0x849F_C360_2CE6_F138),
    ("temperature_m05", 0x7447_05FC_20AD_1C24),
    ("temperature_m06", 0x9CFB_4B29_07FC_DB8E),
    ("temperature_m07", 0xA4C4_6494_27AC_816D),
    ("temperature_m08", 0x852E_2B2E_0EC4_A5C4),
    ("temperature_m09", 0xEFF8_8C86_A169_79A1),
    ("temperature_m10", 0xC750_1249_BD5D_3BEF),
    ("temperature_m11", 0xDFD8_251F_FA37_4785),
    ("temperature_m12", 0x2693_3436_A3F9_A11E),
    ("precipitation_annual", 0xB7FE_D044_E927_1C7E),
    ("precipitation_m01", 0xB08B_A6DE_F32E_4CEB),
    ("precipitation_m02", 0xA9F3_E9B5_476D_6EA9),
    ("precipitation_m03", 0xB9C5_E54F_2084_5811),
    ("precipitation_m04", 0xC4FB_CE56_FAD6_741C),
    ("precipitation_m05", 0x0F0E_DD80_45AD_41DD),
    ("precipitation_m06", 0xA8F9_720C_5346_EA31),
    ("precipitation_m07", 0x670F_A14B_EDCA_7F26),
    ("precipitation_m08", 0x38B8_7043_63E9_DE4E),
    ("precipitation_m09", 0x2838_F629_513E_3207),
    ("precipitation_m10", 0x2019_EB83_BDAB_5410),
    ("precipitation_m11", 0x0696_B7C4_8CD0_69B3),
    ("precipitation_m12", 0x668B_2930_AFF9_C486),
    ("wind_east", 0x2650_ADCC_22AF_D87B),
    ("wind_north", 0xD19E_856D_517B_6511),
    ("current_east", 0x2572_528E_5F45_35C8),
    ("current_north", 0x5F07_3B45_307D_2BA5),
];

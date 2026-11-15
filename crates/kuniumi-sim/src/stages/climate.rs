//! Climate (design: "気候モデル"; milestone 4).
//!
//! Intermediate-complexity, monthly (12 steps per year) model run to a
//! periodic steady state: analytic insolation from obliquity and orbit; a 2D
//! energy-balance model with land/ocean heat capacity and ice albedo; pressure
//! belts whose number scales with rotation rate (not fixed at three cells);
//! surface winds from geostrophy with friction; wind-driven ocean gyres from a
//! Stommel model per basin; SST with upwelling; moisture advection,
//! evaporation, precipitation, orographic lift and rain shadows; ice sheets.
//!
//! Calibrated by running on Earth's real topography and comparing with
//! observed climatology (see docs/roadmap.md, milestone 4).
//!
//! Outputs (planned): `temperature_m01`..`temperature_m12` (K),
//! `precipitation_m01`..`precipitation_m12` (m/s), wind and SST fields.
//!
//! So far (ticket M4-01): the monthly-mean insolation at the top of the
//! atmosphere of every cell, [`insolation::monthly_insolation`], with the
//! calendar of the model (twelve equal months, the northern spring equinox
//! at 79.5/365 of the year). The stage itself is wired by ticket M4-07;
//! until then the runner skips it.

pub mod insolation;

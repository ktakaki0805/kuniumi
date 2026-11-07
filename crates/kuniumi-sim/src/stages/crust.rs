//! Crust to elevation (design: "地殻から標高へ"; milestone 2).
//!
//! - Ocean floor depth from crust age (Parsons & Sclater 1977).
//! - Continental elevation from crust thickness by Airy isostasy
//!   (crust 2800 kg/m³, mantle 3300 kg/m³), using `surface_gravity_m_s2`.
//! - Sea level solved by bisection so that the volume below it equals
//!   `ocean_volume_m3`. Land fraction is a result, never an input.
//!
//! Outputs (planned): `surface_height` (m, relative to the reference radius)
//! and `World::sea_level_m`.

//! Plate tectonics (design: "プレートテクトニクス"; milestones 2 and 6).
//!
//! Milestone 2: one-shot plates. Voronoi plates with Euler-pole rotations;
//! crust thickness from boundary convergence/divergence. No time evolution.
//!
//! Milestone 6: kinematic simulation after Cortial et al. (2019). Rigid plates
//! of crust samples rotate about Euler poles in 1–2 Myr steps for 300–500 Myr,
//! starting from a supercontinent. Handles ocean–continent and ocean–ocean
//! subduction, continental collision and suturing, ridge spreading, rifting
//! along old sutures, and fixed hotspots. Samples are resampled onto the grid
//! about every 10 steps. Records the uplift rate of the final 20–50 Myr for
//! the erosion stage.
//!
//! Outputs (planned): `crust_type`, `crust_thickness` (m), `crust_age` (s),
//! `orogeny_age` (s), `plate_id`, `uplift_rate` (m/s).

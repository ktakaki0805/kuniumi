//! Erosion and landscape evolution (design: "侵食と地形発達"; milestone 3).
//!
//! dh/dt = U − K Q^m S^n + K_d ∇²h + G, with discharge Q from precipitation.
//! Per step: steepest-descent receivers ([`crate::drainage`]); depression
//! routing (Cordonnier, Bovy & Braun 2019, or Priority-Flood); stack order
//! and discharge on the drainage graph of [`crate::drainage`]; implicit
//! stream-power solve (Braun & Willett 2013, n = 1); hillslope diffusion;
//! deposition (Yuan et al. 2019); marine sediment; isostatic rebound.
//! Later: lithology-dependent K and glacial erosion.
//!
//! The implicit solve is sequential along the drainage stack; the future GPU
//! path will need a different algorithm and is only statistically equivalent.

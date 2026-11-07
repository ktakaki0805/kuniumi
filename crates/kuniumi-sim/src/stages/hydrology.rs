//! Rivers and lakes (design: "水文（川・湖）"; milestone 5).
//!
//! Runoff = precipitation − evapotranspiration (Budyko curve); discharge
//! accumulated along the final drainage; depressions filled by water balance
//! into open (fresh) or closed (salt) lakes, processed bottom-up through the
//! depression hierarchy (Fill–Spill–Merge, Barnes et al. 2021); rivers where
//! discharge exceeds a threshold; seasonal rivers from monthly runoff.

//! Erosion and landscape evolution (design: "侵食と地形発達"; milestone 3).
//!
//! dh/dt = U − K Q^m S^n + K_d ∇²h + G, with discharge Q from precipitation.
//! Per step: steepest-descent receivers ([`crate::drainage`]); depression
//! routing on the surface filled by Priority-Flood
//! ([`crate::drainage::fill_depressions`]; the basin graph of Cordonnier,
//! Bovy & Braun 2019 if routing comes to dominate); stack order and
//! discharge on the drainage graph of [`crate::drainage`]; implicit
//! stream-power solve (Braun & Willett 2013, n = 1; [`incise`]); hillslope
//! diffusion; deposition (Yuan et al. 2019); marine sediment; isostatic
//! rebound. Later: lithology-dependent K and glacial erosion.
//!
//! So far (ticket M3-03): river incision against uplift,
//! `∂h/∂t = U − K Qᵐ S` with n = 1. [`erode`] runs
//! [`ErosionSettings::steps`] steps of [`ErosionSettings::time_step_s`] on a
//! surface, each one the first four steps of the design (侵食と地形発達 /
//! 1ステップの処理): it fills the closed depressions
//! ([`fill_depressions`]), builds the drainage graph on the filled surface
//! ([`Drainage::build`]), accumulates the discharge from the precipitation
//! down that graph ([`Drainage::accumulate`]), and erodes the original
//! surface along it ([`incise`]). Until the Climate stage exists the
//! precipitation is [`PLACEHOLDER_PRECIPITATION_M_S`] everywhere. Ticket
//! M3-05 wires the stage into the pipeline, M3-07 makes
//! [`ErosionSettings`] part of the world's settings.
//!
//! Everything here depends only on its inputs, through integer
//! bookkeeping, IEEE `+ − × ÷`, comparisons, [`det::pow`] and exact
//! operations on bit patterns, in the fixed order of the drainage stack
//! (docs/conventions.md, "Determinism", rule 2), so it gives the same
//! result on every platform; it draws no random numbers.
//!
//! The implicit solve is sequential along the drainage stack; the future GPU
//! path will need a different algorithm and is only statistically equivalent.

use kuniumi_math::det;
use kuniumi_math::units::SECONDS_PER_YEAR;
use kuniumi_sphere::SphereGrid;
use kuniumi_world::ErosionSettings;

use crate::drainage::{Drainage, fill_depressions, next_up};

/// The precipitation (m/s) that drives erosion until the Climate stage
/// exists: 1 m/yr everywhere (3.17 × 10⁻⁸ m/s), about Earth's global mean
/// (ticket M3-03).
pub const PLACEHOLDER_PRECIPITATION_M_S: f64 = 1.0 / SECONDS_PER_YEAR;

/// The drainage network that one step of [`incise`] follows, one entry per
/// cell.
#[derive(Clone, Copy, Debug)]
pub struct Network<'a> {
    /// The receiver of each cell, the cell itself for a base level
    /// ([`Drainage::receivers`]).
    pub receivers: &'a [u32],
    /// Every cell exactly once, each after its receiver
    /// ([`Drainage::stack`]).
    pub stack: &'a [u32],
    /// The distance from each cell to its receiver (m), `L_i`: finite and
    /// positive except at the base levels, where it is not read
    /// ([`receiver_distances_m`]).
    pub receiver_distance_m: &'a [f64],
    /// The discharge through each cell (m³/s), `Q_i`: finite and not
    /// negative ([`Drainage::accumulate`] of the precipitation).
    pub discharge_m3_s: &'a [f64],
}

/// One implicit step of river incision against uplift on `height_m` (m, one
/// height per cell, on the datum of `sea_level_m`, in place):
/// `∂h/∂t = U − K Qᵐ S` with n = 1, `U` the rock uplift rate
/// `uplift_rate_m_s` (m/s), `K`, `m` and `Δt` from `settings`, `Q` the
/// discharge (m³/s) and `S` the slope to the receiver (m/m).
///
/// The method is the implicit one of Braun & Willett (2013, "A very
/// efficient O(n), implicit and parallel method to solve the stream power
/// equation governing fluvial incision and landscape evolution",
/// *Geomorphology* 180–181, 170–179) for n = 1: along the stack, each cell
/// after its receiver, whose new height is then known, the backward-Euler
/// update of the cell is explicit. With `a = h_i + Δt U_i` and `b` the
/// water surface the river reaches at the receiver `r`:
///
/// - An ocean cell (`h_i ≤ sea_level_m` at the start of the step) is a
///   fixed base level: unchanged, uplift included. Its water surface is the
///   sea level.
/// - A land base level (a pit, which the routing on the filled surface
///   leaves only on land without ocean) has no receiver and is only
///   uplifted: `h_i = a`. Its water surface is `a`, or the sea level if
///   negative uplift takes it to the sea level or below (it becomes
///   ocean), so that its donors never grade below the sea.
/// - Rivers grade to the water surface: `b` is the sea level where the
///   receiver is an ocean cell, not its floor (which lies kilometres
///   deeper: grading to it would turn the coast into deep sea), and the
///   surface of the lake where the receiver lies under one. The water
///   surface of a land cell is its new height, or that of its receiver if
///   higher.
/// - A cell with `a ≤ b` lies under the water of its receiver, in a lake of
///   a closed depression across which the drainage graph routes the water
///   ([`fill_depressions`]): it is not eroded, only uplifted, `h_i = a`.
///   The update below would raise it towards `b` (deposition, which this
///   law does not describe) or cut the lake floor. A lake drains as its
///   outlet is cut down.
/// - Otherwise `h_i = b + (a − b) / (1 + F)` with `F = Δt K Q_iᵐ / L_i`.
///   This is Braun & Willett's `(a + F b) / (1 + F)` written relative to
///   the base, which keeps the result above `b` in floating point; where
///   rounding still gives `h_i ≤ b` (for `F` beyond about `(a − b)` over a
///   ULP of `b`), `h_i` is the least double above `b`, and it is never
///   above `a`. `F` is 0, and the cell not eroded, where `Δt K` or `Q_iᵐ`
///   is 0, so that no `0 · ∞` occurs where the other overflows; where `F`
///   overflows to `∞`, the cell takes the limit of the update, the least
///   double above `b`, so that no `∞ / ∞` occurs either.
///
/// So every eroded cell ends strictly above the water surface of its
/// receiver and not above `a`, for any time step: erosion never takes a
/// cell below its receiver and never raises one. At steady state
/// (`h_i` unchanged) the slope is `(h_i − b) / L_i = U_i / (K Q_iᵐ)`.
/// Without uplift heights never rise, land stays land, no land cell loses
/// its last strictly lower neighbour and no closed depression appears.
/// With uplift equal on every land cell, no closed depression appears
/// either; that the rounding of `h + Δt U` could merge two heights one ULP
/// apart into a flat where the sum crosses a power of two is not excluded
/// (the tests find none). Uplift that differs between neighbours can dam a
/// valley into a depression: that is tectonics, not the scheme. Uplift on
/// ocean cells is ignored, and the sea level is not changed, so the ocean
/// below it keeps its volume.
///
/// Time O(n): the inputs are checked and every cell is visited once.
/// Finite inputs give finite heights unless `h_i + Δt U_i` exceeds
/// `f64::MAX`.
///
/// # Panics
///
/// If a slice of `network` or `uplift_rate_m_s` does not hold one value per
/// height; if the stack does not hold every cell once, or a cell comes
/// before its receiver in it; if a height, an uplift rate or the sea level
/// is not finite; if a discharge is negative or not finite; if the
/// distance of a cell that is not a base level to its receiver is not
/// finite and positive; if `settings` is not valid
/// ([`ErosionSettings::validate`]).
pub fn incise(
    network: &Network<'_>,
    uplift_rate_m_s: &[f64],
    sea_level_m: f64,
    settings: &ErosionSettings,
    height_m: &mut [f64],
) {
    check_network(network, uplift_rate_m_s, sea_level_m, settings, height_m);
    let Network {
        receivers,
        stack,
        receiver_distance_m,
        discharge_m3_s,
    } = *network;
    let dt_s = settings.time_step_s;
    let dt_k = dt_s * settings.erodibility;
    let mut water_m = vec![0.0; height_m.len()];
    for &i in stack {
        let i = i as usize;
        let h = height_m[i];
        if h <= sea_level_m {
            water_m[i] = sea_level_m;
            continue;
        }
        let a = h + dt_s * uplift_rate_m_s[i];
        let r = receivers[i] as usize;
        if r == i {
            height_m[i] = a;
            water_m[i] = if a > sea_level_m { a } else { sea_level_m };
            continue;
        }
        let b = water_m[r];
        if a <= b {
            height_m[i] = a;
            water_m[i] = b;
            continue;
        }
        let q_m = det::pow(discharge_m3_s[i], settings.discharge_exponent);
        let f = if dt_k == 0.0 || q_m == 0.0 {
            0.0
        } else {
            dt_k * q_m / receiver_distance_m[i]
        };
        let mut eroded = if f == 0.0 {
            a
        } else if f.is_infinite() {
            b
        } else {
            b + (a - b) / (1.0 + f)
        };
        if eroded <= b {
            eroded = next_up(b);
        }
        if eroded > a {
            eroded = a;
        }
        height_m[i] = eroded;
        water_m[i] = eroded;
    }
}

/// The checks of [`incise`] (its "Panics").
fn check_network(
    network: &Network<'_>,
    uplift_rate_m_s: &[f64],
    sea_level_m: f64,
    settings: &ErosionSettings,
    height_m: &[f64],
) {
    let n = height_m.len();
    for (what, len) in [
        ("receiver", network.receivers.len()),
        ("receiver distance", network.receiver_distance_m.len()),
        ("discharge", network.discharge_m3_s.len()),
        ("uplift rate", uplift_rate_m_s.len()),
    ] {
        assert_eq!(
            len, n,
            "incise needs one {what} per cell: {len} values for {n} heights"
        );
    }
    if let Err(e) = settings.validate() {
        panic!("incise needs valid settings: {e}");
    }
    assert!(
        sea_level_m.is_finite(),
        "incise needs a finite sea level, got {sea_level_m}"
    );
    if let Some(i) = height_m.iter().position(|h| !h.is_finite()) {
        panic!("incise needs finite heights: cell {i} has {}", height_m[i]);
    }
    if let Some(i) = uplift_rate_m_s.iter().position(|u| !u.is_finite()) {
        panic!(
            "incise needs finite uplift rates: cell {i} has {}",
            uplift_rate_m_s[i]
        );
    }
    // Written so that NaN fails.
    if let Some(i) = network
        .discharge_m3_s
        .iter()
        .position(|&q| !(q >= 0.0 && q.is_finite()))
    {
        panic!(
            "incise needs finite, non-negative discharges: cell {i} has {}",
            network.discharge_m3_s[i]
        );
    }
    let mut position = vec![usize::MAX; n];
    assert_eq!(
        network.stack.len(),
        n,
        "incise needs a stack of every cell once: {} entries for {n} cells",
        network.stack.len()
    );
    for (k, &i) in network.stack.iter().enumerate() {
        let i = i as usize;
        assert!(
            i < n && position[i] == usize::MAX,
            "incise needs a stack of every cell once: entry {k} is {i}"
        );
        position[i] = k;
    }
    for (i, &r) in network.receivers.iter().enumerate() {
        let r = r as usize;
        if r == i {
            continue;
        }
        assert!(
            r < n && position[r] < position[i],
            "incise needs each cell after its receiver in the stack: cell {i} drains to {r}"
        );
        let l = network.receiver_distance_m[i];
        assert!(
            l > 0.0 && l.is_finite(),
            "incise needs a finite, positive distance to the receiver: cell {i} has {l}"
        );
    }
}

/// The distance (m) from each cell of `grid` to its receiver in
/// `drainage`, `L_i = d_{i r(i)} R`, with `d` the centre distance of the
/// two cells ([`SphereGrid::centre_distances_rad`]) and `R = radius_m`, the
/// product that [`Drainage::build`] divides the drop by; 0 at the base
/// levels.
///
/// # Panics
///
/// If `drainage` was not built on a grid of as many cells as `grid`, if a
/// receiver is not a neighbour of its cell in `grid`, or if `radius_m` is
/// not finite and positive.
pub fn receiver_distances_m(grid: &SphereGrid, drainage: &Drainage, radius_m: f64) -> Vec<f64> {
    assert_eq!(
        drainage.len(),
        grid.len(),
        "receiver_distances_m needs the drainage graph of the grid: {} cells for {}",
        drainage.len(),
        grid.len()
    );
    // Written so that NaN fails.
    assert!(
        radius_m > 0.0 && radius_m.is_finite(),
        "receiver_distances_m needs a finite positive radius, got {radius_m}"
    );
    (0..grid.len() as u32)
        .map(|i| {
            let r = drainage.receiver(i);
            if r == i {
                return 0.0;
            }
            let k = grid
                .neighbors(i)
                .iter()
                .position(|&j| j == r)
                .unwrap_or_else(|| panic!("the receiver {r} of cell {i} is not its neighbour"));
            grid.centre_distances_rad(i)[k] * radius_m
        })
        .collect()
}

/// What the erosion of a surface depends on besides its heights: the grid,
/// the planet and the forcing, the rates one value per cell of the grid.
#[derive(Clone, Copy, Debug)]
pub struct ErosionInputs<'a> {
    /// The grid of the surface.
    pub grid: &'a SphereGrid,
    /// The planet's radius (m), finite and positive.
    pub radius_m: f64,
    /// The sea level (m), on the datum of the heights, finite.
    pub sea_level_m: f64,
    /// The rock uplift rate of each cell (m/s), finite; the Tectonics
    /// stage's `uplift_rate`.
    pub uplift_rate_m_s: &'a [f64],
    /// The precipitation of each cell (m/s), finite and not negative;
    /// [`PLACEHOLDER_PRECIPITATION_M_S`] until the Climate stage exists.
    pub precipitation_m_s: &'a [f64],
}

/// Runs `settings.steps` steps of river erosion on `height_m` (m, one
/// height per cell of `inputs.grid`, in place), and reports
/// `(k + 1) / steps` to `progress` after step `k`; with 0 steps the heights
/// are unchanged and nothing is reported.
///
/// Each step (module documentation): the closed depressions of the
/// current surface are filled ([`fill_depressions`]); the drainage graph is
/// built on the filled surface ([`Drainage::build`]); the discharge
/// `Q_i` (m³/s) is accumulated down it from `P_i a_i R²`, the
/// precipitation `P_i` (m/s) times the area of the cell (m²)
/// ([`Drainage::accumulate`]); and the original surface is eroded along the
/// graph by [`incise`], with the distances of [`receiver_distances_m`].
/// The routing follows the filled surface, the erosion acts on the
/// original one (design, 侵食と地形発達 / 1ステップの処理).
///
/// The sea level is fixed, and so are the ocean cells; see [`incise`] for
/// what the step guarantees. A step costs a fill (O(n log n)), the graph,
/// one accumulation and the solve (O(n)).
///
/// # Panics
///
/// If `height_m`, the uplift or the precipitation does not hold one value
/// per cell; if the radius is not finite and positive; if a precipitation
/// rate is negative or not finite; if a height, an uplift rate or the sea
/// level is not finite; if `settings` is not valid. All are checked before
/// the first step, so with 0 steps as well.
pub fn erode(
    inputs: &ErosionInputs<'_>,
    settings: &ErosionSettings,
    height_m: &mut [f64],
    progress: &mut dyn FnMut(f64),
) {
    let ErosionInputs {
        grid,
        radius_m,
        sea_level_m,
        uplift_rate_m_s,
        precipitation_m_s,
    } = *inputs;
    let n = grid.len();
    for (what, len) in [
        ("height", height_m.len()),
        ("uplift rate", uplift_rate_m_s.len()),
        ("precipitation rate", precipitation_m_s.len()),
    ] {
        assert_eq!(
            len, n,
            "erode needs one {what} per cell: {len} values for {n} cells"
        );
    }
    // Written so that NaN fails.
    assert!(
        radius_m > 0.0 && radius_m.is_finite(),
        "erode needs a finite positive radius, got {radius_m}"
    );
    if let Some(i) = precipitation_m_s
        .iter()
        .position(|&p| !(p >= 0.0 && p.is_finite()))
    {
        panic!(
            "erode needs finite, non-negative precipitation: cell {i} has {}",
            precipitation_m_s[i]
        );
    }
    assert!(
        sea_level_m.is_finite(),
        "erode needs a finite sea level, got {sea_level_m}"
    );
    if let Some(i) = height_m.iter().position(|h| !h.is_finite()) {
        panic!("erode needs finite heights: cell {i} has {}", height_m[i]);
    }
    if let Some(i) = uplift_rate_m_s.iter().position(|u| !u.is_finite()) {
        panic!(
            "erode needs finite uplift rates: cell {i} has {}",
            uplift_rate_m_s[i]
        );
    }
    if let Err(e) = settings.validate() {
        panic!("erode needs valid settings: {e}");
    }

    let inflow_m3_s: Vec<f64> = (0..n as u32)
        .map(|i| precipitation_m_s[i as usize] * grid.area_sr(i) * radius_m * radius_m)
        .collect();
    for k in 0..settings.steps {
        let filled_m = fill_depressions(grid, height_m, sea_level_m);
        let drainage = Drainage::build(grid, &filled_m, sea_level_m, radius_m);
        let discharge_m3_s = drainage.accumulate(&inflow_m3_s);
        let distance_m = receiver_distances_m(grid, &drainage, radius_m);
        let network = Network {
            receivers: drainage.receivers(),
            stack: drainage.stack(),
            receiver_distance_m: &distance_m,
            discharge_m3_s: &discharge_m3_s,
        };
        incise(&network, uplift_rate_m_s, sea_level_m, settings, height_m);
        progress(f64::from(k + 1) / f64::from(settings.steps));
    }
    debug_assert!(
        height_m.iter().all(|h| h.is_finite()),
        "erode produced a height that is not finite"
    );
}

#[cfg(test)]
mod tests;

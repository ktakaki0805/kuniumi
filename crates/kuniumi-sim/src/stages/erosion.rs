//! Erosion and landscape evolution (design: "侵食と地形発達"; milestone 3).
//!
//! dh/dt = U − K Q^m S^n + K_d ∇²h + G, with discharge Q from precipitation.
//! Per step: steepest-descent receivers ([`crate::drainage`]); depression
//! routing on the surface filled by Priority-Flood
//! ([`crate::drainage::fill_depressions`]; the basin graph of Cordonnier,
//! Bovy & Braun 2019 if routing comes to dominate); stack order and
//! discharge on the drainage graph of [`crate::drainage`]; implicit
//! stream-power solve (Braun & Willett 2013, n = 1; [`incise`]); hillslope
//! diffusion ([`diffuse`]); deposition (Yuan et al. 2019); marine sediment;
//! isostatic rebound. Later: lithology-dependent K and glacial erosion.
//!
//! So far (tickets M3-03 and M3-04): river incision against uplift and
//! hillslope diffusion, `∂h/∂t = U − K Qᵐ S + K_d ∇²h` with n = 1.
//! [`erode`] runs [`ErosionSettings::steps`] steps of
//! [`ErosionSettings::time_step_s`] on a surface, each one the first five
//! steps of the design (侵食と地形発達 / 1ステップの処理): it fills the closed
//! depressions ([`fill_depressions`]), builds the drainage graph on the
//! filled surface ([`Drainage::build`]), accumulates the discharge from the
//! precipitation down that graph ([`Drainage::accumulate`]), erodes the
//! original surface along it ([`incise`]), and diffuses the land
//! ([`diffuse`]). Until the Climate stage exists the
//! precipitation is [`PLACEHOLDER_PRECIPITATION_M_S`] everywhere.
//!
//! The stage ([`run`], ticket M3-05) runs [`erode`] on the world's
//! `surface_height` with its `uplift_rate` and its erosion settings,
//! `Settings::erosion` (ticket M3-07), after the Crust stage (the runner
//! skips the Climate stage until it exists), and adds the rivers of the
//! final surface ([`Rivers`]): [`DRAINAGE_AREA`] (m²), [`DISCHARGE`] (m³/s)
//! and [`DEPRESSION_DEPTH`] (m). With the default settings, 1,000 steps of
//! 10⁴ yr, the continental interiors, which have no uplift, are planed down
//! to the sea level and the highest peaks depend on the resolution (5.2 to
//! 15.0 km at 650,000 cells, seeds 1–5); calibrating them waits for
//! isostatic rebound and the calibration of milestone 7 (tickets M7-03 and
//! M7-05).
//!
//! Everything here depends only on its inputs, through integer
//! bookkeeping, IEEE `+ − × ÷`, `sqrt`, comparisons, [`det::pow`] and exact
//! operations on bit patterns, in the fixed order of the drainage stack and
//! of the cells (docs/conventions.md, "Determinism", rule 2), so it gives
//! the same result on every platform; it draws no random numbers.
//!
//! The implicit solves are sequential, along the drainage stack and over
//! the cells; the future GPU path will need different algorithms and is
//! only statistically equivalent.

use kuniumi_math::det;
use kuniumi_math::units::SECONDS_PER_YEAR;
use kuniumi_sphere::SphereGrid;
use kuniumi_world::{ErosionSettings, FieldError, World};

use crate::SimError;
use crate::drainage::{Drainage, depression_depth_m, fill_depressions, next_up};
use crate::operators::{Laplacian, Solve};
use crate::stages::crust::SURFACE_HEIGHT;
use crate::stages::tectonics::UPLIFT_RATE;

/// The field of the drainage area (m²): the area of the cell and of every
/// cell upstream of it on the final surface ([`Rivers::drainage_area_m2`]).
pub const DRAINAGE_AREA: &str = "drainage_area";

/// The field of the discharge (m³/s): the precipitation on the cell and on
/// every cell upstream of it ([`Rivers::discharge_m3_s`]).
pub const DISCHARGE: &str = "discharge";

/// The field of the depression depth (m): how deep the water of a closed
/// depression, filled to its spill level, stands on the cell, and a few
/// ULPs on drained flats ([`Rivers::depression_depth_m`]).
pub const DEPRESSION_DEPTH: &str = "depression_depth";

/// The rivers of a surface: what the Erosion stage adds besides the eroded
/// heights, on the drainage graph of the surface routed across its closed
/// depressions ([`routed_drainage`]), one value per cell.
#[derive(Clone, Debug, PartialEq)]
pub struct Rivers {
    /// The drainage area (m²): `a_i R²` summed over the cell and every cell
    /// upstream of it ([`Drainage::accumulate`]). An ocean cell counts the
    /// rivers that enter it; the base levels sum to `4πR²` within rounding.
    pub drainage_area_m2: Vec<f64>,
    /// The discharge (m³/s): `P_i a_i R²` summed likewise, `P_i` the
    /// precipitation (m/s); with the uniform placeholder `P` it is `P`
    /// times the drainage area, within rounding.
    pub discharge_m3_s: Vec<f64>,
    /// The depression depth (m): the filled surface minus the surface
    /// ([`depression_depth_m`]), positive where a closed depression would
    /// hold water, the ε gradient of the fill alone (at most `2n` ULPs of
    /// the largest height) on drained flats, and +0 elsewhere and on the
    /// ocean. The lakes of the Hydrology stage start from it (ticket M5-02).
    pub depression_depth_m: Vec<f64>,
}

impl Rivers {
    /// The rivers of the surface `height_m` (m, one height per cell of
    /// `inputs.grid`) under `inputs.precipitation_m_s`, on the planet and
    /// sea level of `inputs` (the uplift is not used).
    ///
    /// # Panics
    ///
    /// As [`routed_drainage`], and if the precipitation does not hold one
    /// value per cell or a rate is negative or not finite.
    pub fn new(inputs: &ErosionInputs<'_>, height_m: &[f64]) -> Self {
        let ErosionInputs {
            grid,
            radius_m,
            sea_level_m,
            ..
        } = *inputs;
        let (drainage, filled_m) = routed_drainage(grid, height_m, sea_level_m, radius_m);
        let area_m2: Vec<f64> = (0..grid.len() as u32)
            .map(|i| grid.area_sr(i) * radius_m * radius_m)
            .collect();
        Self {
            drainage_area_m2: drainage.accumulate(&area_m2),
            discharge_m3_s: drainage.accumulate(&inflow_m3_s(inputs)),
            depression_depth_m: depression_depth_m(height_m, &filled_m),
        }
    }
}

/// The drainage graph of the surface `height_m` (m, one height per cell of
/// `grid`, on the datum of `sea_level_m`) routed across its closed
/// depressions, on a planet of radius `radius_m` (m): the graph that
/// [`Drainage::build`] builds on the surface filled by [`fill_depressions`],
/// and that filled surface. Each step of [`erode`] follows it, and
/// [`Rivers`] accumulates along it.
///
/// # Panics
///
/// As [`fill_depressions`] and [`Drainage::build`]: if `height_m` does not
/// hold one height per cell, a height or the sea level is not finite, or
/// the radius is not finite and positive.
pub fn routed_drainage(
    grid: &SphereGrid,
    height_m: &[f64],
    sea_level_m: f64,
    radius_m: f64,
) -> (Drainage, Vec<f64>) {
    let filled_m = fill_depressions(grid, height_m, sea_level_m);
    let drainage = Drainage::build(grid, &filled_m, sea_level_m, radius_m);
    (drainage, filled_m)
}

/// `P_i a_i R²` (m³/s) for every cell of `inputs.grid`: the water that
/// falls on it.
///
/// # Panics
///
/// If the precipitation does not hold one value per cell, or a rate is
/// negative or not finite.
fn inflow_m3_s(inputs: &ErosionInputs<'_>) -> Vec<f64> {
    let ErosionInputs {
        grid,
        radius_m,
        precipitation_m_s,
        ..
    } = *inputs;
    assert_eq!(
        precipitation_m_s.len(),
        grid.len(),
        "the rivers need one precipitation rate per cell: {} values for {} cells",
        precipitation_m_s.len(),
        grid.len()
    );
    // Written so that NaN fails.
    if let Some(i) = precipitation_m_s
        .iter()
        .position(|&p| !(p >= 0.0 && p.is_finite()))
    {
        panic!(
            "the rivers need finite, non-negative precipitation: cell {i} has {}",
            precipitation_m_s[i]
        );
    }
    (0..grid.len() as u32)
        .map(|i| precipitation_m_s[i as usize] * grid.area_sr(i) * radius_m * radius_m)
        .collect()
}

/// Runs the Erosion stage on `world` (a
/// [`StageFn`](crate::pipeline::StageFn); module documentation): erodes
/// the field [`SURFACE_HEIGHT`] of the Crust stage by [`erode`], with the
/// field [`UPLIFT_RATE`] of the Tectonics stage, the settings
/// `world.settings.erosion`, the planet's radius, the world's sea level and
/// the placeholder precipitation [`PLACEHOLDER_PRECIPITATION_M_S`] on every
/// cell (until the Climate stage exists), replacing the heights in place;
/// then adds the [`Rivers`] of the final surface as the fields
/// [`DRAINAGE_AREA`] (`m^2`), [`DISCHARGE`] (`m^3/s`) and
/// [`DEPRESSION_DEPTH`] (`m`), in that order. The progress of [`erode`],
/// `(k + 1) / steps` after step `k`, is passed on. The sea level and the
/// other fields are not changed.
///
/// One step costs about 0.1 s at 650,000 cells, so the default 1,000 steps
/// take about two minutes there (107 s of the 112 s that `kuniumi
/// generate` takes for seed 1, ticket M3-05).
///
/// # Errors
///
/// [`SimError::Field`] if the world already has one of the three fields;
/// the world is then unchanged and no progress is reported.
///
/// # Panics
///
/// If the world lacks [`SURFACE_HEIGHT`] or [`UPLIFT_RATE`], which the
/// runner's earlier stages always add; as [`erode`] for settings, heights
/// or uplift that are not valid (the runner checks the settings first).
pub fn run(world: &mut World, progress: &mut dyn FnMut(f64)) -> Result<(), SimError> {
    for name in [DRAINAGE_AREA, DISCHARGE, DEPRESSION_DEPTH] {
        if world.field(name).is_some() {
            return Err(FieldError::Duplicate { name: name.into() }.into());
        }
    }
    let precipitation_m_s = vec![PLACEHOLDER_PRECIPITATION_M_S; world.grid.len()];
    let (height_m, rivers) = {
        let inputs = ErosionInputs {
            grid: &world.grid,
            radius_m: world.params.radius_m,
            sea_level_m: world.sea_level_m,
            uplift_rate_m_s: input(world, UPLIFT_RATE, "Tectonics"),
            precipitation_m_s: &precipitation_m_s,
        };
        let mut height_m = input(world, SURFACE_HEIGHT, "Crust").to_vec();
        erode(&inputs, &world.settings.erosion, &mut height_m, progress);
        let rivers = Rivers::new(&inputs, &height_m);
        (height_m, rivers)
    };
    if let Some(field) = world.field_mut(SURFACE_HEIGHT) {
        field.values = height_m;
    }
    world.add_field(DRAINAGE_AREA, "m^2", rivers.drainage_area_m2)?;
    world.add_field(DISCHARGE, "m^3/s", rivers.discharge_m3_s)?;
    world.add_field(DEPRESSION_DEPTH, "m", rivers.depression_depth_m)?;
    Ok(())
}

/// The values of the field `name` of `world`, an output of the `stage`
/// stage.
///
/// # Panics
///
/// If the world has no such field.
fn input<'a>(world: &'a World, name: &str, stage: &str) -> &'a [f64] {
    match world.field(name) {
        Some(field) => &field.values,
        None => panic!("the Erosion stage needs the field `{name}` of the {stage} stage"),
    }
}

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

/// Runs `settings.steps` steps of erosion on `height_m` (m, one height per
/// cell of `inputs.grid`, in place), and reports `(k + 1) / steps` to
/// `progress` after step `k`; with 0 steps the heights are unchanged and
/// nothing is reported.
///
/// Each step (module documentation): the closed depressions of the
/// current surface are filled ([`fill_depressions`]); the drainage graph is
/// built on the filled surface ([`Drainage::build`]); the discharge
/// `Q_i` (m³/s) is accumulated down it from `P_i a_i R²`, the
/// precipitation `P_i` (m/s) times the area of the cell (m²)
/// ([`Drainage::accumulate`]); the original surface is eroded along the
/// graph by [`incise`], with the distances of [`receiver_distances_m`];
/// and the land is diffused by [`diffuse`], with the [`Laplacian`] of the
/// grid, built once per call (neither where
/// [`ErosionSettings::diffusivity_m2_s`] is 0, when `diffuse` would change
/// nothing). The routing follows the filled surface, the
/// erosion acts on the original one (design, 侵食と地形発達 /
/// 1ステップの処理).
///
/// The sea level is fixed, and so are the ocean cells. See [`incise`] for
/// what the incision guarantees (all of it holds for the whole step when
/// [`ErosionSettings::diffusivity_m2_s`] is 0) and [`diffuse`] for the
/// diffusion, which keeps the land above the sea level but may raise a
/// valley floor or, rarely, leave a cell without a lower neighbour (the
/// next step's fill routes across it). A step costs a fill
/// (O(n log n)), the graph, one accumulation, the incision (O(n)) and the
/// diffusion (O(n) per iteration of its solver; one or two with the
/// default settings).
///
/// # Panics
///
/// If `height_m`, the uplift or the precipitation does not hold one value
/// per cell; if the radius is not finite and positive; if a precipitation
/// rate is negative or not finite; if a height, an uplift rate or the sea
/// level is not finite; if `settings` is not valid; if `Δt K_d / R²` is
/// not finite. All are checked before the first step, so
/// with 0 steps as well.
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
    let diffusion = step_diffusion(radius_m, settings);
    assert!(
        diffusion.is_finite(),
        "erode needs a finite step Δt K_d / R², got {diffusion}"
    );

    // Without diffusivity `diffuse` changes nothing: skip it, and the
    // Laplacian (about 30 ms at 650,000 cells).
    let laplacian = (settings.diffusivity_m2_s > 0.0).then(|| Laplacian::new(grid));
    let inflow_m3_s = inflow_m3_s(inputs);
    for k in 0..settings.steps {
        let (drainage, _) = routed_drainage(grid, height_m, sea_level_m, radius_m);
        let discharge_m3_s = drainage.accumulate(&inflow_m3_s);
        let distance_m = receiver_distances_m(grid, &drainage, radius_m);
        let network = Network {
            receivers: drainage.receivers(),
            stack: drainage.stack(),
            receiver_distance_m: &distance_m,
            discharge_m3_s: &discharge_m3_s,
        };
        incise(&network, uplift_rate_m_s, sea_level_m, settings, height_m);
        if let Some(laplacian) = &laplacian {
            diffuse(laplacian, radius_m, sea_level_m, settings, height_m);
        }
        progress(f64::from(k + 1) / f64::from(settings.steps));
    }
    debug_assert!(
        height_m.iter().all(|h| h.is_finite()),
        "erode produced a height that is not finite"
    );
}

/// `c = Δt K_d / R²` (dimensionless): the diffusion of one step on the
/// unit sphere, as [`Laplacian::implicit_step`] computes it.
fn step_diffusion(radius_m: f64, settings: &ErosionSettings) -> f64 {
    settings.time_step_s * settings.diffusivity_m2_s / (radius_m * radius_m)
}

/// One step of hillslope diffusion, `∂h/∂t = K_d ∇²h`, on `height_m` (m,
/// one height per cell of the grid of `laplacian`, on the datum of
/// `sea_level_m`, in place), with `K_d` and `Δt` from `settings`
/// ([`ErosionSettings::diffusivity_m2_s`], [`ErosionSettings::time_step_s`])
/// on a planet of radius `radius_m` (design document, 侵食と地形発達 /
/// 1ステップの処理, step 5: 「斜面拡散。拡散項を陰的に解きます。」).
///
/// The land cells (`h_i > sea_level_m`) diffuse by one backward-Euler step
/// of the finite-volume Laplacian ([`Laplacian::implicit_step`]), solved by
/// the conjugate gradient method; the ocean cells are a fixed boundary and
/// keep their heights. A land cell sees an ocean neighbour at the sea level,
/// not at its floor: hillslope creep carries soil to the shore, as rivers
/// grade to the water surface in [`incise`], and a floor kilometres below
/// would drag the coast under the sea, which is the transport of marine
/// sediment (design, step 7; milestone 7). With the sea level as the
/// boundary value the exact step keeps every land cell above the sea level
/// and not above the highest land cell (a discrete maximum principle), for
/// every time step; where the solver's tolerance would put a land cell at
/// or below the sea level, it is set to the least double above it, so land
/// stays land. On a surface without ocean the step conserves
/// `Σ h_i a_i` up to rounding.
///
/// With [`ErosionSettings::diffusivity_m2_s`] 0 the heights are unchanged
/// bit for bit. Returns what the solver did.
///
/// # Panics
///
/// If `height_m` does not hold one height per cell; if `radius_m` is not
/// finite and positive; if a height or the sea level is not finite; if
/// `settings` is not valid ([`ErosionSettings::validate`]); if
/// `Δt K_d / R²` is not finite.
pub fn diffuse(
    laplacian: &Laplacian,
    radius_m: f64,
    sea_level_m: f64,
    settings: &ErosionSettings,
    height_m: &mut [f64],
) -> Solve {
    let n = laplacian.len();
    assert_eq!(
        height_m.len(),
        n,
        "diffuse needs one height per cell: {} values for {n} cells",
        height_m.len()
    );
    // Written so that NaN fails.
    assert!(
        radius_m > 0.0 && radius_m.is_finite(),
        "diffuse needs a finite positive radius, got {radius_m}"
    );
    assert!(
        sea_level_m.is_finite(),
        "diffuse needs a finite sea level, got {sea_level_m}"
    );
    if let Some(i) = height_m.iter().position(|h| !h.is_finite()) {
        panic!("diffuse needs finite heights: cell {i} has {}", height_m[i]);
    }
    if let Err(e) = settings.validate() {
        panic!("diffuse needs valid settings: {e}");
    }
    let c = step_diffusion(radius_m, settings);
    assert!(
        c.is_finite(),
        "diffuse needs a finite step Δt K_d / R², got {c}"
    );

    let land: Vec<bool> = height_m.iter().map(|&h| h > sea_level_m).collect();
    let mut values: Vec<f64> = height_m
        .iter()
        .zip(&land)
        .map(|(&h, &is_land)| if is_land { h } else { sea_level_m })
        .collect();
    let solve = laplacian.implicit_step(
        radius_m,
        settings.diffusivity_m2_s,
        settings.time_step_s,
        &land,
        &mut values,
    );
    for i in 0..n {
        if land[i] {
            let v = values[i];
            height_m[i] = if v > sea_level_m {
                v
            } else {
                next_up(sea_level_m)
            };
        }
    }
    debug_assert!(
        height_m.iter().all(|h| h.is_finite()),
        "diffuse produced a height that is not finite"
    );
    solve
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod diffusion_tests;

#[cfg(test)]
mod stage_tests;

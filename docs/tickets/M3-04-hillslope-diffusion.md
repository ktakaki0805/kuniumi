# M3-04 Hillslope diffusion

- Status: done
- Milestone: M3
- Depends on: M3-03
- Crates: kuniumi-world (`ErosionSettings`), kuniumi-sim (stages/erosion, a reusable `operators` module)
- Design: 侵食と地形発達 / 基本方程式; 侵食と地形発達 / 1ステップの処理
- Decisions: none

## Goal

Hillslopes soften over time (`K_d ∇²h`), using a finite-volume Laplacian on the
Voronoi grid that later stages (temperature diffusion, sediment) reuse.

## Scope

In:
- Finite-volume operator: flux across edge (i, j) =
  `K_d · (h_j − h_i) / (d_ij R) · (ℓ_ij R)`; the change of cell i is the sum of
  fluxes divided by `a_i R²`, using the edge lengths and centre distances of
  M1-06.
- Backward-Euler step solved by conjugate gradient on the area-weighted
  symmetric system, with dot products computed sequentially in a fixed order;
  tolerance and iteration cap documented.
- Applied to land cells; ocean cells act as a fixed boundary for now.
- `K_d` in m²/s with an Earth-like default (value and source documented).

Out:
- Nonlinear (critical-slope) diffusion.

## Implementation notes

Gaps found while refining (the ticket was `ready`): where `K_d` lives and
whether `erode` (M3-03) runs the diffusion; the value the ocean boundary
holds; the solver's formulation, preconditioner, stopping rule and cap;
how mass conservation is made to hold within 1e-12 for every surface; the
default `K_d` and its source; how the criteria are measured.

### API

`kuniumi-sim`, new module `operators` (tests in `operators/tests.rs`, so
`cargo test -p kuniumi-sim operators` selects them):

```rust
/// Relative residual at which the implicit step stops.
pub const SOLVER_TOLERANCE: f64 = 1.0e-12;
/// Iterations after which it stops regardless.
pub const SOLVER_MAX_ITERATIONS: u32 = 10_000;

/// The finite-volume Laplacian of a grid: per neighbour pair the
/// coupling w_ij = ℓ_ij / d_ij (dimensionless, bit-for-bit symmetric), per
/// cell its area a_i (sr) and Σ_j w_ij, in the CSR order of the grid.
pub struct Laplacian { /* private */ }
impl Laplacian {
    pub fn new(grid: &SphereGrid) -> Self;
    pub fn len(&self) -> usize;
    pub fn is_empty(&self) -> bool;
    /// w_ij of the neighbours of i, in the order of grid.neighbors(i).
    pub fn weights(&self, i: u32) -> &[f64];
    /// (∇²u)_i = Σ_j w_ij (u_j − u_i) / (a_i R²), in [u]/m².
    pub fn apply(&self, radius_m: f64, values: &[f64]) -> Vec<f64>;
    /// One backward-Euler step of ∂u/∂t = κ ∇²u on the cells with
    /// free[i]; the others are fixed (Dirichlet) at their value in `values`.
    pub fn implicit_step(&self, radius_m: f64, diffusivity_m2_s: f64,
                         time_step_s: f64, free: &[bool],
                         values: &mut [f64]) -> Solve;
}
/// What a solve did: iterations, final relative residual.
pub struct Solve { pub iterations: u32, pub relative_residual: f64 }
impl Solve { pub fn converged(&self) -> bool; }
```

`kuniumi-world`: `ErosionSettings` gains `diffusivity_m2_s` (`K_d`, m²/s,
after `discharge_exponent`), `DEFAULT_DIFFUSIVITY_M2_S` and its check in
`validate` (finite and not negative, by name in field order; 0 turns the
diffusion off). M3-07 then records five settings instead of four (its text
is updated).

`kuniumi-sim`, `stages::erosion`:

```rust
/// One backward-Euler step of hillslope diffusion, ∂h/∂t = K_d ∇²h, on
/// the land cells (h > sea level); the ocean is a fixed boundary.
pub fn diffuse(laplacian: &Laplacian, radius_m: f64, sea_level_m: f64,
               settings: &ErosionSettings, height_m: &mut [f64]) -> Solve;
```

`erode` builds the `Laplacian` once and runs `diffuse` after `incise` in
every step: step 5 of the design (侵食と地形発達 / 1ステップの処理,
「斜面拡散。拡散項を陰的に解きます。」). Without it M3-05, whose scope does
not mention diffusion, would wire a stage without step 5.

### Discretisation

The ticket's finite-volume operator: across the Voronoi edge (i, j) the
flux is `K_d (h_j − h_i)/(d_ij R) · ℓ_ij R` (m³/s); cell i changes by the sum
over its edges divided by `a_i R²`. So
`dh_i/dt = (K_d / R²) Σ_j w_ij (h_j − h_i) / a_i`, `w_ij = ℓ_ij / d_ij`,
with ℓ and d the edge lengths and centre distances of M1-06 (radians, the
same bits in both directions, so `w_ij = w_ji` bit for bit). This is the
two-point flux of a Voronoi (orthogonal) mesh; for `h = z` its error
falls to first order in the spacing on the jittered grid (measured before
refining: the largest error, relative to `2/R²`, is 1.2e-2, 4.0e-3 and
1.3e-3 at 10⁴, 10⁵ and 10⁶ cells). (Corrected in review round 1: this
holds for `z` because its Hessian is isotropic; for general fields the
error of the operator cell by cell does not shrink, while the solutions of
the implicit step converge. See the Result.)

Backward Euler with `c = Δt κ / R²` (dimensionless): for each free cell
`a_i h'_i + c Σ_j w_ij (h'_i − h'_j) = a_i h_i`, where `h'_j` of a fixed
neighbour is its given value. That is `A h' = b` with
`A = M + c L` (M the diagonal of areas, L the weighted graph Laplacian of
the free cells, its diagonal including the couplings to fixed cells),
symmetric positive definite, and
`b_i = a_i h_i + c Σ_{j fixed} w_ij h_j`. A is an M-matrix: in exact
arithmetic the new values lie within the range of the old values and the
boundary values (discrete maximum principle), for every Δt.

### Solver

- Conjugate gradient preconditioned by the diagonal of A (Jacobi; Shewchuk
  1994, "An Introduction to the Conjugate Gradient Method Without the
  Agonizing Pain", B3), started from the old values, unknowns the free
  cells; every vector operation and dot product sequential in index order.
- The system is divided by `1 + c` (entries `a_i/(1 + c)` and
  `c/(1 + c) · w_ij`), so that its entries stay of the order of the areas
  and weights for any finite `c`; the stopping rule is relative and does
  not change.
- Stops when `‖b − A h'‖₂ ≤ SOLVER_TOLERANCE · ‖b‖₂` (the recursively
  updated residual; squared norms compared), 1e-12, or after
  `SOLVER_MAX_ITERATIONS` = 10,000 iterations; `Solve` reports which.
  Measured before refining (10⁴ and 10⁵ cells): about 15 iterations at the
  explicit limit, a few hundred at 10⁴ times it, under 2,500 at 10¹²
  times it. In the Erosion stage `c Σ_j w_ij / a_i` is about 10⁻⁶
  (default settings, 650,000 cells): one or two iterations.
- Mass: in exact arithmetic `Σ a_i h_i` does not change where no free cell
  has a fixed neighbour (on the connected grid: every cell free, the
  closed domain of the criterion). The solver's tolerance alone does not
  keep it within 1e-12: for a lone spike its relative error reaches 1e-11
  (measured). On such a closed domain the solved values are therefore
  shifted by the constant `(Σ a_i h_i − Σ a_i h'_i) / Σ a_i` (sums in index
  order). The constant vector is the eigenvector of `M⁻¹A` (of the pencil
  `(A, M)`) with the smallest eigenvalue, 1, where the solver's error
  collects (wording corrected in review round 1, R1-3); the shift removes
  that error's mass and changes no difference between cells. (Rejected:
  re-evaluating the step in flux form, `h_i + c Σ_j w_ij (h'_j − h'_i) /
  a_i`, which conserves mass exactly but amplifies the rounding of h' by
  `c Σ w / a`: an undershoot of 3 % of a spike at 10¹⁴ times the explicit
  limit, measured.)
- With `c = 0` (no diffusivity) the values are unchanged bit for bit.

### Erosion: the boundary and land

- Land cells (`h > sea level` when `diffuse` is called) are free; ocean
  cells are fixed and keep their heights. A land cell sees an ocean
  neighbour at the sea level, not at its floor: hillslope creep carries soil
  to the shore, as rivers grade to the water surface in M3-03, and the deep
  floor kilometres below would drag the coast under the sea, which is the
  marine transport of step 7 (M7). By the maximum principle land then stays
  above the sea level for every Δt; where the solver's tolerance would put
  a land cell at or below it, it is set to the least double above it.
- `diffuse` panics on the conditions of `incise` that apply (lengths,
  finite heights and sea level, valid settings, radius finite and
  positive).

### Default diffusivity

`DEFAULT_DIFFUSIVITY_M2_S` = 5.0 × 10⁻³ m²/yr = 1.584 × 10⁻¹⁰ m²/s: the
mean of 49 ± 37 cm²/yr that McKean et al. (1993) computed for 34 sites in
California, Oregon and Washington from the rates of colluvium accumulation
in hollows reported by Reneau (1988), as summarised by Dietrich et al.
(2003, "Geomorphic transport laws for predicting landscape form and
dynamics", AGU Geophysical Monograph 135, 103–132), who also cite 0.0032
m²/yr for the Oregon Coast Range (Roering et al. 1999) as the coefficient
of the nonlinear, critical-slope law, the linear diffusivity at low
gradients; within the range of field estimates 4.4 × 10⁻⁴ to
3.6 × 10⁻² m²/yr (Fernandes & Dietrich 1997, as quoted by Ruetenik et al.
2023, Earth Surf. Dynam. 11, 865–880). At the
default grid (28 km between centres) this is sub-grid: over the default
10 Myr the diffusion length `√(K_d t)` is 220 m, so the step changes the
heights by little; the operator matters at finer grids and for its other
uses (temperature, sediment). Noted for the calibration of M3-05.

### Tests

`operators/tests.rs`:
- `laplacian_by_hand`: on the octahedron every weight is
  `acos(1/3) / (π/2)`; the Laplacian of given values by hand; a constant
  gives exactly 0; weights symmetric bit for bit on the 10,000-cell grid.
- `laplacian_of_z_converges` (criterion 1): `h = z` on grids of 10⁴, 10⁵
  and 10⁶ cells (10⁶ in release only: 33 s to build in debug), R = 6,371
  km; the largest and the area-weighted RMS error against `−2z/R²`,
  relative to `2/R²`, printed; each at most half the one of the coarser
  grid; the largest below 2e-2 at 10⁴ cells.
- `implicit_step_by_hand`: the octahedron with a spike, all free and with
  fixed cells, against the 3×3 system its symmetry reduces to, solved by
  hand (relative 1e-12).
- `implicit_step_solves_the_system`: on the 10,000-cell grid, the residual
  `a_i h'_i + c Σ w_ij (h'_i − h'_j) − b_i` recomputed in the test from the
  grid's ℓ, d and a, for several c, below the tolerance's bound.
- `mass_is_conserved_on_a_closed_domain` (criterion 2): every cell free, a
  spike, uniform random heights and the Earth world of seed 1 (10,000
  cells), 10 steps each at `N = c max_i(Σ_j w_ij / a_i)` of 10⁻⁶ to 10⁸;
  `|Σ a h' − Σ a h| / |Σ a h|` (compensated sums in the test) at most 1e-12
  per step; the largest reported.
- `large_steps_are_stable` (criterion 3): N from 10⁻² to 10⁸, where
  forward Euler at the same step grows (checked in the test): a spike of
  1,000 m on a closed domain and on a cap of free cells whose rim is fixed
  at 0, and random heights; 10 steps each: no NaN; every value within the
  range of the old values and boundary values to within 1e-9 of that range
  (the tolerance's effect); the spike cell stays the highest and its value
  falls strictly at every step until it is within that margin of its limit
  (the mean on the closed domain, 0 on the cap), and stays there.
- `a_degree_one_harmonic_decays_at_its_rate`: `h = z`, one step, ratio
  `h'/z` within 2e-2 of `1/(1 + 2c)` (the discrete eigenvalue is 2 up to
  the operator's error).
- `zero_diffusivity_changes_nothing` (bit for bit), `fixed_cells_stay`,
  the panics, and a golden fingerprint of one step on fixed inputs
  (recorded after the tests above pass).

`stages/erosion/diffusion_tests.rs`:
- `diffuse_by_hand`: the octahedron with one land cell and ocean around it
  at −4,000 m: the land cell relaxes towards the sea level,
  `h' = a h / (a + 4 c w)` (sea level 0), the ocean unchanged.
- `land_stays_land`: Earth worlds of seeds 1–3 and N up to 10⁸: ocean
  unchanged bit for bit, every land cell above the sea level.
- `erode_diffuses_after_incising`: one step of `erode` equals `incise` on
  the routing of that step followed by `diffuse`, bit for bit; with
  `diffusivity_m2_s = 0` the incision of M3-03 alone (its golden value
  unchanged: the M3-03 tests now set the diffusivity to 0, since they
  test the incision's own guarantees, such as "no height rises", which
  diffusion does not keep).
- A golden fingerprint of 10 default steps of `erode` (with diffusion) on
  the Earth world of seed 1 at 10,000 cells.
- `kuniumi-world`: the default, `validate` of the new field.

### Images and timing

A scratch tool outside the repository (a Cargo package in the scratchpad
with path dependencies on the workspace crates, as in M3-03) draws into
`out/M3-04/`: the error of the Laplacian of z at 10⁵ cells (equirectangular,
diverging colours); a spike before and after steps at several N
(orthographic close-ups); the Earth world of seed 1 at 100,000 cells after
diffusion alone with a large diffusivity (smoothing visible) and the
difference that the default diffusion makes to 100 steps of `erode`. It
times one `diffuse` and one `erode` step in release at 650,000 cells
(best of 5).

## Acceptance criteria

- [x] For `h = z` (a degree-1 spherical harmonic) the discrete Laplacian
      approximates `−2·z / R²`, with the error shrinking under refinement:
      largest and RMS errors relative to `2/R²` reported for 10⁴, 10⁵ and
      10⁶ cells, each at most half that of the coarser grid
      (`laplacian_of_z_converges`).
- [x] Mass conservation on a closed domain (all cells land): Σ h·a stays
      constant within 1e-12 relative per step, for steps from 10⁻⁶ to 10⁸
      times the explicit limit (`mass_is_conserved_on_a_closed_domain`).
- [x] Unconditionally stable for large Δt (up to 10⁸ times the explicit
      limit, where forward Euler grows): no NaN, values within the range
      of the old and boundary values; a single spike decays monotonically
      (`large_steps_are_stable`).
- [x] `erode` diffuses after incising (bit for bit `incise` then
      `diffuse`); ocean cells unchanged and land above the sea level; with
      zero diffusivity the M3-03 golden value is unchanged.
- [x] `mise run ci` passes.

## Verification

```sh
mise exec -- cargo test -p kuniumi-sim operators
mise exec -- cargo test -p kuniumi-sim --release operators
mise exec -- cargo test -p kuniumi-sim erosion
mise exec -- cargo test -p kuniumi-world
```

Images in `out/M3-04/` from the scratch tool (see Images and timing).

## Result

What was done:
- `kuniumi-sim`: the module `operators` with `Laplacian` (`new`, `len`,
  `is_empty`, `weights`, `apply`, `implicit_step`), `Solve`,
  `SOLVER_TOLERANCE` (1e-12) and `SOLVER_MAX_ITERATIONS` (10,000), as
  planned; 26 tests in `operators/tests.rs`.
- `kuniumi-sim`, `stages::erosion`: `diffuse` (land free, ocean fixed and
  seen at the sea level, land kept above it), run by `erode` after
  `incise` in every step with the `Laplacian` built once per call; `erode`
  checks that `Δt K_d / R²` is finite before the first step, and skips
  the diffusion and the Laplacian where `diffusivity_m2_s` is 0; 14 tests
  in `stages/erosion/diffusion_tests.rs`. The tests of M3-03 run the step
  without diffusion (`settings` sets `diffusivity_m2_s` to 0) and pass
  unchanged, their golden value included; four of their helpers became
  `pub(super)` for the new file.
- `kuniumi-world`: `ErosionSettings::diffusivity_m2_s`,
  `DEFAULT_DIFFUSIVITY_M2_S` = 5 × 10⁻³ m²/yr (1.584 × 10⁻¹⁰ m²/s, sources
  in its documentation) and its check in `validate`; the tests extended.
- `CHANGELOG.md` (Added, Changed); notes from this ticket in M3-05 (the
  sub-grid diffusivity and the cost), M3-06 (the sensitivity of the
  routing on flats, below) and M3-07 (five settings instead of four).

Deviations from the plan: none in the method. In the tests:
`weights_are_symmetric` and `apply_scales_with_the_radius` are separate
tests; `erode_steps_compose_with_diffusion` was added, and
`degree_two_fields_converge_in_the_step_not_in_apply` after review round 1
(R1-1); `a_degree_one_harmonic_decays_at_its_rate` bounds the deviation by
3e-2, not the planned 2e-2 (measured 6.3e-3 at most); the margin of the
stability test is 1e-8 m (the largest excursion measured is 9.6e-14 m);
`diffuse_by_hand` compares to within 1e-12 of the starting height, not of
the result, because the step adds the solved change to the height (at
c = 10⁶ the cell falls from 100 m to 6.7e-5 m and keeps the absolute
precision of 100 m). A planned hand test of the land guard (one cell one
ULP above the sea level) was dropped: there the change is below the
solver's tolerance and nothing is solved; the guard is exercised by
`land_stays_land_and_the_ocean_stays` instead (without it, a coastal cell
of seed 1 ends 2e-13 m below the sea level at N = 10⁶).

Evidence (N is the step in units of the explicit limit,
`Δt κ / R² · max_i Σ_j w_ij / a_i`, computed in the tests from the grid):

1. Criterion 1 (`laplacian_of_z_converges`, release; 10⁴ and 10⁵ cells in
   debug): largest error against `−2z/R²`, relative to `2/R²`: 1.196e-2,
   4.038e-3 and 1.299e-3 at 10⁴, 10⁵ and 10⁶ cells; area-weighted RMS
   2.863e-3, 9.108e-4 and 2.901e-4. Each refinement divides them by 2.96,
   3.11 and 3.14, 3.14: first order in the spacing (√10 = 3.16).
   This holds for `z` because its Hessian on the sphere is isotropic (`−z`
   times the identity), which cancels the error the two-point flux makes
   where a centre is not its cell's centroid. For other fields the
   operator is not consistent cell by cell (review round 1, R1-1, measured
   independently by the reviewer and here): for the degree-2 harmonics
   `3z² − 1` and `xy` the error of `apply`, relative to the largest
   `|∇²u|`, stays at about 3.2 % RMS and 17–24 % at most at 10⁴, 10⁵ and
   10⁶ cells. The solutions of the implicit step converge nonetheless, at
   second order in the spacing as measured (the bounds proven by Eymard,
   Gallouët & Herbin 2000 for the two-point flux on orthogonal meshes are
   of first order; review round 2, R2-1): one step of
   c = 0.1 lies within 2.3e-4, 2.4e-5 and 2.6e-6 of `u / (1 + 0.6)` and
   one of c = 1 within 8.9e-5, 9.1e-6 and 9.8e-7 of `u / 7` (relative to
   the largest `|u|`; `degree_two_fields_converge_in_the_step_not_in_apply`
   checks 10⁴ and 10⁵ cells). The module documentation and `apply` say so:
   a later stage that needs `∇²u` cell by cell needs a consistent operator
   or a centroidal grid.
   `laplacian_error_z_100000.png`: per-cell noise, no large-scale bias,
   strongest at the equator and fading towards the poles (the first-order
   term follows ∇z, zero at the poles); the close-up shows it cell by cell
   with no pattern. `laplacian_by_hand` (octahedron),
   `weights_are_symmetric`, `apply_scales_with_the_radius`.
2. Criterion 2 (`mass_is_conserved_on_a_closed_domain`): every cell free,
   10 steps at each N of 10⁻⁶, 10⁻³, 10⁻², 1, 10², 10⁴, 10⁶ and 10⁸ on a
   spike, uniform random heights and the Earth world of seed 1: the
   largest relative change of `Σ a h` in a step (compensated sums in the
   test) is 1.7e-16, 1.5e-16 and 1.9e-16 (bound 1e-12). Without the shift
   (a mutant) the spike loses up to 1e-11.
3. Criterion 3 (`large_steps_are_stable`): forward Euler at the same step
   takes the spike cell to 1,000 (1 − N) m (asserted to grow for N ≥ 10²).
   Backward Euler at every N above, 10 steps each, on a spike of 1,000 m
   on the closed domain and on a cap of free cells with a rim fixed at 0,
   and on random heights: all values finite; the largest step outside the
   range of the old values and the boundary 9.6e-14 m; the spike cell
   falls strictly at every step while more than 1e-8 m above its limit and
   stays the highest cell, then stays within 1e-8 m of it (closed: the
   mean, 5.28e-2 m, within 2.3e-7 m after 10 steps at N = 10⁴ and reached
   at 10⁶; cap: 0). Iterations: 1 at N = 10⁻⁶, about 15 at N = 1, up to
   131 at 10², up to 675 at 10⁸.
   `implicit_step_by_hand` (octahedron, closed and with a fixed cell,
   c from 10⁻³ to 10⁶, against the 3×3 symmetric system solved by hand);
   `implicit_step_solves_the_system` (residual recomputed from the grid:
   2.0e-14, 5.4e-13, 7.6e-13 and 7.0e-13 at N = 10⁻⁴, 1, 10², 10⁴);
   `a_degree_one_harmonic_decays_at_its_rate` (z multiplied by
   `1/(1 + 2c)` within 3.9e-5, 2.0e-4 and 6.3e-3 at c = 10⁻³, 1, 10³).
   Images `spike_N1_after1.png`, `spike_N1_after10.png`,
   `spike_N100_after1.png`, `spike_N100_after10.png` (100,000 cells,
   log scale from 1 mm to 1 km): a round blob that widens with steps and
   with N, highest at the spike, falling monotonically outwards, no ring
   or checkerboard; the lowest value of any run is positive (1e-14 to
   2e-13 m).
4. Criterion 4 (`erode_diffuses_after_incising`): one step of `erode`
   equals `incise` on that step's routing followed by `diffuse`, bit for
   bit, at the default diffusivity, at 10⁶ times it and at 0; the
   diffusion changes the step's heights by at most 1.4e-5 m (default) and
   14 m (10⁶ ×) on the Earth world of seed 1 at 10,000 cells.
   `land_stays_land_and_the_ocean_stays` (seeds 1–3 and the dry world, N
   from 10⁻² to 10⁸, 3 steps): ocean cells unchanged bit for bit, every
   land cell above the sea level and not above the highest land cell
   before. `erode_steps_compose_with_diffusion`. The M3-03 golden value
   `0x73A6_193E_F4E0_5D72` is unchanged with `diffusivity_m2_s` 0, and all
   M3-03 tests pass.
5. Images of the Earth world of seed 1 at 100,000 cells (`out/M3-04/`,
   the scratch tool `<scratchpad>/m3_04_maps`; elevation colours with hill
   shading), viewed:
   - `diffusion_alone_N10x10_100000.png` (10 steps of N = 10 without
     incision, `K_d` = 5.3 × 10⁵ m²/yr): the belts are blurred into broad
     swells (highest 6,320 → 4,819 m), the land falls smoothly towards the
     coasts, the coastline (29,623 land cells) and the ocean are unchanged.
   - `after_100_default_100000.png` (100 steps, default settings): the
     dendritic valleys of M3-03. `diffusion_effect_100_default_100000.png`
     (against the same run without diffusion): zero over most of the
     globe; 4,310 cells differ by more than 1 m, in a few patches by the
     belts, by up to 938 m (22.6° S, 178.8° E,
     `after_100_{nodiff,default}_100000_maxdiff.png`): there a knickpoint
     took another way into a flat plateau. The diffusion changes a step's
     heights by some 10⁻⁵ m, but on the flat interiors such changes decide
     the routing (M3-02's ε order); noted in M3-06. Land statistics agree
     (mean 421 m, median 243 m, highest 8,564 m in both).
   - `after_100_n01_100000.png` and its close-ups (`K_d` = 5,263 m²/yr,
     N = 0.1 per step, an effective large-scale diffusivity): the same
     network with wider, smoother valleys and rounded ridges (mean land
     elevation 528 m instead of 421 m, highest 8,342 m); the interior
     close-up (3.0° N, 89.3° E) shows the valleys of the run without
     diffusion softened. The default `K_d` is sub-grid at this resolution
     (diffusion length 220 m in 10 Myr); a coarse-grid value is left to
     M3-05's calibration.
6. Golden values, recorded after 1–5, the same in debug and release:
   `golden_fingerprint_of_an_implicit_step` `0x85B3_7FE4_3AB1_E898` (31
   iterations); `golden_fingerprint_of_an_eroded_and_diffused_earth_world`
   `0xB715_B888_F60A_A864`. The golden values of M3-01, M3-02 and M3-03
   are unchanged.
7. Timing (release, this machine under a load of about 3.7, Earth world
   of seed 1 at 650,000 cells, best of 5): `Laplacian::new` 31 ms; one
   `diffuse` 13 ms with the default settings (N = 6.6 × 10⁻⁷, one
   iteration), 74 ms at N = 1 (13 iterations), 596 ms at N = 100 (116);
   one `erode` step 127 ms, the Laplacian included (built once per call,
   so about 96 ms per step in a long run), against 81 ms without diffusion
   (as in M3-03: since review round 1, R1-5, `erode` skips the diffusion
   and the Laplacian where the diffusivity is 0).
8. Test-first: against the stubs 21 of the 25 `operators` tests failed
   (the four that passed hold for the identity: the unbroken call, zero
   diffusion, mass conservation, the radius scaling of zeros), 10 of the
   14 diffusion tests (the four that passed: land stays land, zero
   diffusivity, composition, the unbroken call), and the two
   `kuniumi-world` tests of the new field. Mutations, each caught by the
   final tests (release): weights `d/ℓ` (9 tests), no mass shift (2), a
   shift on open domains too (7), the tolerance not squared (9), the
   coupling `c` instead of `c/(1 + c)` (5), `apply` without `R²` (3), the
   flux sign reversed (10), the ocean floor as the boundary (2), no land
   guard (1), `erode` without `diffuse` (2), land at `h ≥` sea level (4).
9. `mise run ci` passes: fmt, clippy `-D warnings`, rustdoc, the browser
   checks, the tests in debug and release, MSRV 1.85.0, cargo-deny.

Review:
- Round 1 (CHANGES_REQUESTED): R1-1 major, the operator was called
  consistent to first order on the jittered grid, which holds only for
  `z` (isotropic Hessian); for degree-2 harmonics the error of `apply`
  cell by cell stays at about 3 % RMS while the implicit step's solutions
  converge at second order. Fixed: module documentation, a warning on
  `apply`, CHANGELOG and Result corrected, the test
  `degree_two_fields_converge_in_the_step_not_in_apply` added. R1-2 minor,
  the attribution of 49 ± 37 cm²/yr to McKean et al. (1993) doubted:
  checked against the text of Dietrich et al. (2003), which gives it so;
  the citation now says it was computed from Reneau's (1988) hollow
  fills (rebuttal accepted in round 2). Nits fixed: R1-3 the eigenvector
  wording (pencil `(A, M)`), R1-4 `converged()` now set by the loop, R1-5
  `erode` skips the Laplacian and the diffusion at zero diffusivity (one
  step back to 81 ms), R1-6 test counts, a listed deviation, a wrapped
  line. The reviewer re-ran CI and every criterion and probed degree-2
  fields with its own tool.
- Round 2 (APPROVE): R1-1 to R1-6 resolved; the reviewer reproduced the
  convergence orders (degree 2 to 4, with a Dirichlet rim too) and ran
  `diffuse` at 650,000 cells up to N = 10⁸ (at most 3,181 iterations).
  R2-1 minor, the second order was credited to Eymard, Gallouët & Herbin
  (2000), whose bounds are first order: now stated as measured. R2-2 nit,
  the eigenvector wording in the Implementation notes, and R2-3 nit,
  Roering et al.'s 3.2 × 10⁻³ m²/yr is the coefficient of the nonlinear
  law: both fixed. `mise run ci` passes after the fixes.

Follow-ups: none.

# M4-09 Multigrid for the climate model at 650,000 cells

- Status: done
- Milestone: M4
- Depends on: M4-02
- Crates: kuniumi-sim (multigrid, operators, sparse, stages/climate)
- Design: 気候モデル / ExoPlaSim との使い分け; 実装ロードマップ / 計算コストの目安
- Decisions: D-08

## Goal

The energy balance of M4-02, the sea surface with currents of M4-05 and the
potential of the divergent wind of M4-06 run at the default 650,000 cells in
minutes rather than hours, at the grid's full resolution, before M4-07 wires
the Climate stage.

## Scope

In:
- Measured before this ticket (release, per model year): the energy balance
  0.37 s at 10,000 cells, 3.2 s at 40,000, 33 s at 160,000, 277 s at
  650,000; a periodic state needs 12 to 20 years, so about an hour at
  650,000 cells. The Jacobi-preconditioned conjugate gradient needs more
  iterations the finer the grid (78, 153, 333, 662 per solve), because the
  diffusion `D L` is stiff at the grid's scale (10⁵ W m⁻² K⁻¹ against
  `C/Δt + B` ≈ 18 for land). The sea surface with currents (M4-05) solves a
  non-symmetric system every sub-step by BiCGSTAB with the Jacobi
  preconditioner (5.07 s per model year at 40,000 cells against 3.73 s
  alone). The moisture (M4-06) solves once a month a Poisson equation for
  the divergent wind (426, 832 and 1,660 iterations at 10,000, 40,000 and
  160,000 cells).
- A deterministic smoothed-aggregation algebraic multigrid (new module
  `kuniumi_sim::multigrid`), as the preconditioner of the conjugate gradient
  method for the symmetric systems and of BiCGSTAB for the non-symmetric
  system of the sea surface, keeping the full resolution and the design's
  table ("解像度: 全球格子と同じ").
- Used by: the energy balance's sub-steps, warm start and annual correction
  (M4-02); the sea surface's sub-steps, return vectors, warm start and
  correction (M4-05); the moisture's potential (M4-06).
- Timing at 10,000 to 650,000 cells recorded; results compared with the
  Jacobi solvers at 160,000 cells.

Out (the split approved by the owner; "Split" below):
- The moisture's upwind advection systems at 650,000 cells (new draft
  M4-14): the multigrid does not suit them.
- A solver of the periodic state itself for planets with little heat
  transport or near an ice–albedo bifurcation (new draft M4-15; was
  criterion 3 of this ticket).
- The Stommel solve of M4-04 (once per climate, 37 s at 650,000 cells):
  unchanged.
- The hillslope diffusion of the Erosion stage (`Laplacian::implicit_step`,
  one or two iterations a step): unchanged, bit for bit.
- The physics of the energy balance (M4-02, M4-10, M4-11).

## Implementation notes

The ticket was `draft` and left the method to the owner: a coarse climate
grid (about 40,000 cells, the climate interpolated back to the planet grid)
or a multigrid preconditioner at the full resolution. A prototype outside
the repository (`<scratchpad>/m4_09_proto`, a Cargo package with path
dependencies, as for M4-01 to M4-06) measured both sides on the
Earth-parameter world of seed 1 (relief of the Crust stage), release, this
machine.

### Measurements (prototype)

Smoothed aggregation (below) against the Jacobi solvers of the code, one
solve each, from a representative state (`288 K − 28 K P₂(sin φ)`, 5 K
cooler on land); the currents, upwelling and winds of the sea surface and
the moisture come from a synthetic seasonal temperature of the same form,
as the prototype does not run the periodic states at 650,000 cells. The
multigrid's setup is separate and done once per matrix (0.32 s at 650,000
cells, 0.47 s with the advection).

| System | Cells | Jacobi: iterations, s | Multigrid: iterations, s |
| --- | --- | --- | --- |
| Energy balance sub-step, Earth `D` | 40,000 | 182, 0.080 | 12, 0.027 |
| | 160,000 | 387, 0.836 | 13, 0.127 |
| | 650,000 | 762, 6.69 | 15, 0.574 |
| Sub-step, `990 D⊕` (240 h, 10 bar) | 650,000 | 3,300, 28.9 | 26, 0.99 |
| Sub-step, `D⊕/91` (8 h, 0.1 bar) | 650,000 | 97, 0.85 | 10, 0.41 |
| Warm start (`B a + D L`) | 650,000 | 2,144, 18.7 | 23, 0.88 |
| Sea surface sub-step (BiCGSTAB) | 40,000 | 122, 0.108 | 8, 0.037 |
| | 160,000 | 280, 1.07 | 8, 0.162 |
| | 650,000 | 771, 11.9 | 10, 0.740 |
| Moisture potential (Poisson) | 40,000 | 1,259, 0.548 | 24, 0.052 |
| | 160,000 | 2,597, 5.51 | 25, 0.245 |
| | 650,000 | 5,186, 46.1 (a synthetic divergence; January's winds 2,867, 27.2) | 30, 1.14 |
| Moisture advection (BiCGSTAB) | 160,000 | 157, 0.460 | 84, 1.50 |
| | 650,000 | 305, 3.98 | diverges |

The multigrid's iterations hardly grow with the grid (12, 13, 15), so its
speed-up grows with it: 12 times at 650,000 cells for the energy balance,
16 times for the sea surface, 40 times for the potential. Its solutions
differ from the Jacobi ones by at most 4.2 × 10⁻⁶ K (both stop at the same
relative residual). Estimated per model year at 650,000 cells: the energy
balance about 30 s (277 s), Earth's periodic state (12 years) about 6–7
minutes (about an hour); the sea surface with currents about 40 s (about
570 s), its periodic state about 8 minutes (about 2 hours); the twelve
potentials about 14 s (5 to 9 minutes).

The moisture's upwind advection at 650,000 cells has Courant numbers up to
8,450 per sub-step and takes 305 Jacobi iterations (4.0 s) per solve, a few
solves per sub-step (the active set, about 4.6 on average in M4-06): about
15 minutes per model year, 45 for its periodic state of 3 years. The
multigrid built from the symmetric part does not suit it (84 iterations and
slower at 160,000 cells; no convergence in 2,000 at 650,000); symmetric
Gauss–Seidel and Gauss–Seidel in a downwind order as preconditioners save a
third to a half (650,000 cells: 68 and 41 iterations, 2.16 and 2.75 s;
160,000: 38 and 28, 0.30 and 0.33 s). That is a different problem
(advection, not diffusion): proposed as M4-14.

A coarse climate grid of 40,000 cells, with the solvers of today: the energy
balance's periodic state 45 s, the sea surface's 61 s, the moisture's 29 s
(tickets M4-02, M4-05, M4-06), about 2–3 minutes in all, faster still with
the multigrid; but the western boundary currents, the coastal upwelling
bands (50–150 km on Earth) and the slopes that make rain shadows would be
resolved at about 110 km instead of 30 km, the coasts of the coarse grid
differ from the planet's (the sea surface needs extrapolating to the
planet's coastal ocean cells), the orographic precipitation would need the
planet grid's relief in a downscaling step, and the design's table changes.

### Decision (chosen with the owner)

The multigrid at the full resolution (this plan), with the split below and
the time targets of criterion 1, chosen with the owner when the ticket was
refined. Not chosen: the coarse climate grid (about 40,000 cells,
aggregation of heat capacity and albedo, interpolation back with the lapse
rate on each cell's own height, a fine-grid orographic term), which would
have changed the design's table.

### The multigrid (`multigrid.rs`, tests in `multigrid/tests.rs`)

Smoothed aggregation (P. Vaněk, J. Mandel & M. Brezina 1996, "Algebraic
multigrid by smoothed aggregation for second and fourth order elliptic
problems", Computing 56, 179–196), for the M-matrices of the climate model
(positive diagonal, non-positive entries off it, weakly diagonally
dominant rows):

- **Strength.** `j` is strongly connected to `i` if
  `a_ij² > θ² a_ii a_jj`, `θ` = 0.08 on the finest level, halved on each
  coarser one (Vaněk et al. 1996, §5). Computed on the symmetric matrix
  that defines the hierarchy (below).
- **Aggregation**, in index order (deterministic): pass 1, each cell that
  is not yet aggregated and whose strong neighbours are all free starts an
  aggregate of itself and them; pass 2, each remaining cell joins the
  aggregate of pass 1 of its strongest neighbour (ties to the lowest index);
  pass 3, each remaining cell with strong neighbours starts an aggregate of
  itself and its free strong neighbours. A cell without strong neighbours
  (a mass that dominates its row, on coarse levels) belongs to no
  aggregate: the smoother alone solves it.
- **Prolongator.** The tentative `P₀` (1 where a cell belongs to an
  aggregate, the constant near-null space of the diffusion), smoothed by one
  damped Jacobi step of the filtered matrix (weak entries added to the
  diagonal): `P = (I − ω D_F⁻¹ A_F) P₀`, `ω = 4/(3ρ)` with `ρ` the
  Gershgorin bound of `D_F⁻¹ A_F` (no eigenvalue iteration).
- **Coarse operators** `A_c = Pᵀ A P` (Galerkin), each row's columns sorted,
  each entry summed in a fixed order. Coarsening stops at 400 cells or
  less, or when a level keeps more than 80 % of the cells; the coarsest
  level is factorised densely (Cholesky for a symmetric hierarchy, LU
  without pivoting for an M-matrix otherwise).
- **V-cycle**, one forward Gauss–Seidel sweep before the coarse correction
  and one backward sweep after it (symmetric, so a valid preconditioner of
  the conjugate gradient method for a symmetric matrix), in index order.
- **Non-symmetric systems** (the sea surface): the aggregates and
  prolongators come from the symmetric matrix of the energy balance
  (`a(C/Δt + B) + D L`), the coarse operators are `Pᵀ A P` of the
  non-symmetric matrix, the coarsest level is solved by LU.

```rust
/// The hierarchy of a sparse matrix.
pub struct Multigrid { /* levels: A, P, Pᵀ; the coarsest factor */ }
impl Multigrid {
    /// The hierarchy of the symmetric M-matrix `matrix`.
    pub fn new(matrix: &SparseMatrix) -> Self;
    /// The hierarchy of `matrix` with the aggregates and prolongators of
    /// the symmetric M-matrix `symmetric` (the same size and pattern).
    pub fn with_prolongators_of(matrix: &SparseMatrix, symmetric: &SparseMatrix) -> Self;
    pub fn levels(&self) -> usize;
    /// Σ nnz of every level over nnz of the finest (1.35 at 650,000 cells).
    pub fn operator_complexity(&self) -> f64;
    /// One V-cycle for `A z = r` from `z = 0`.
    pub fn apply(&self, residual: &[f64], correction: &mut [f64]);
}
/// The conjugate gradient method preconditioned by `multigrid` (of the
/// symmetric `matrix`), as `Laplacian::solve` otherwise: the same stops,
/// the same `Solve`.
pub fn conjugate_gradient(matrix: &SparseMatrix, multigrid: &Multigrid, rhs: &[f64],
                          tolerance: f64, max_iterations: u32, values: &mut [f64]) -> Solve;
```

`sparse.rs` gains `bicgstab_multigrid` (the same method as `bicgstab`, with
the V-cycle in place of the diagonal, right-preconditioned, the same stops,
breakdowns and restarts) and, for the multigrid, a matrix product and
transpose. `operators::Laplacian` gains `matrix(mass, coupling) ->
SparseMatrix` (the matrix of `solve`, entries in neighbour order sorted by
column). `Laplacian::solve` and `implicit_step` stay as they are, bit for
bit (the Erosion stage's golden values do not change).

### Use in the climate model

- `EnergyBalance::new` builds the hierarchy of the sub-step matrix
  `a(C/Δt + B) + D L` once (the same in every sub-step: the albedo enters
  only the right-hand side); `warm_start` builds that of `B a + D L` once;
  `correct` that of `B' a + D L` each year (0.3 s at 650,000 cells). The
  solves keep their tolerances (10⁻¹⁰, the correction 10⁻⁶).
- `SeaSurfaceModel::new` builds one hierarchy for the twelve months'
  sub-step matrices from the annual-mean upwelling (the months differ only
  in the upwelling on the diagonal, so one preconditioner serves them all
  and the memory stays at one hierarchy, about 80 MB at 650,000 cells), and
  one each for the warm start and the correction.
- `MoistureModel::new` builds the hierarchy of the potential's matrix
  `10⁻⁹ a + L` once (the same every month).
- The `Solve` reports keep their meaning (iterations are the preconditioned
  method's).

### Comparison with the full-resolution Jacobi solvers

The model is the same; only the solver changes, so the comparison at
160,000 cells is against the Jacobi solvers of today on the same world:
zonal means (10° bands) and seasonal ranges of the periodic state of the
energy balance and of the sea surface with currents within 2 × 10⁻³ K
(twice the periodic tolerance, as the two runs may stop a year apart;
prototype single solves agree to 4 × 10⁻⁶ K). At 650,000 cells (a
different world, as the world depends on the cell count) the zonal means of
the energy balance are compared with North & Coakley's fit as in M4-02, and
the maps are viewed.

### Tests (`cargo test -p kuniumi-sim multigrid`, and the climate modules)

Reference values from Gaussian elimination in the tests (independent of the
code under test) and the Jacobi solvers of the code.

- `multigrid`: `aggregates_cover_the_connected_cells` (every cell with a
  strong neighbour in exactly one aggregate, aggregates connected by strong
  edges, isolated cells only where the mass dominates);
  `the_prolongator_keeps_constants_away_from_the_boundary`
  (`P 1 = 1` where no row is filtered, for the graph Laplacian);
  `coarse_operators_are_galerkin` (`Pᵀ A P` against a dense product, and
  symmetric bit for bit for a symmetric `A`); `solves_match_a_dense_solution`
  (random SPD and non-symmetric M-matrices of 5 to 400 rows, the grid
  systems on 2,000 cells: within 10⁻¹² of the largest value);
  `iterations_do_not_grow_with_the_grid` (the Earth sub-step on 2,500,
  10,000 and 40,000 cells: at most 20, within 3 of each other);
  `robust_over_the_planet_range` (`D⊕/91` to `990 D⊕`: at most 30);
  `the_v_cycle_is_symmetric` (`xᵀ B y = yᵀ B x` to 10⁻¹² for a symmetric
  hierarchy); `repeatable_bit_for_bit`; panics (sizes, non-finite entries,
  a non-positive diagonal).
- `sparse`: `bicgstab_multigrid` as `bicgstab`'s tests (stops, exact start,
  zero right-hand side, the breakdown and restart).
- Climate modules: every existing test, unchanged in its tolerances; new
  `multigrid_matches_jacobi` per module (energy balance, sea surface,
  moisture potential, 10,000 cells in debug and 40,000 in release): the
  periodic state within 2 × 10⁻³ K, the potential within 10⁻¹⁰ of its
  largest value.
- Golden values: the fingerprints of M4-02
  (`golden_fingerprint_of_earth_temperature`), M4-03 (circulation, from
  those temperatures), M4-04 (currents), M4-05 (sea surface) and M4-06
  (moisture) change, as the solutions change within the solvers'
  tolerances; recorded again after the tests above and the images. The
  pipeline's worlds do not change (the Climate stage is not wired yet);
  stated in the commit's `World-breaking:` line and in `CHANGELOG.md`.

### Timing and images

A scratch tool (`<scratchpad>/m4_09_tool`) times, per model year and for
the periodic state, the energy balance and the sea surface with currents at
10,000, 40,000, 160,000 and 650,000 cells, and the moisture's months
(potentials) at the same sizes, in release; at 160,000 cells also with the
Jacobi solvers for the comparison above. It renders into `out/M4-09/`: the
annual-mean surface temperature, the seasonal range and the annual SST at
650,000 cells, and the difference multigrid − Jacobi at 160,000 cells.

### Split (approved by the owner)

The draft had three independent parts; this ticket keeps the first:

1. M4-09 (this ticket): the multigrid and its use in the diffusion systems.
2. M4-14 (new, draft): the moisture's upwind advection at 650,000 cells
   (a downwind-ordered Gauss–Seidel or ILU preconditioner, or another
   formulation), with the measurements above; depends on M4-09.
3. M4-15 (new, draft): a solver of the periodic state itself (Anderson
   acceleration of the year map, Newton–Krylov shooting or similar,
   deterministic) for the thin, fast planets and the near-bifurcation
   cases of M4-02, with the criterion that was this ticket's third;
   depends on M4-09.

M4-07 (the Climate stage) then depends on M4-06, M4-09 and M4-14; M4-15
can follow M4-07, whose stage reports an unconverged climate as a warning.

### Size

About 1,400 lines: the multigrid with its documentation about 550, its
tests about 450, the preconditioned BiCGSTAB and its tests about 200, the
use in the three models about 150, documents the rest. Over the guideline
of 800; kept as one ticket because the solver has no use without the models
and the models' tests are what prove it.

### Design document

`実装ロードマップ / 計算コストの目安` gains the measured cost of the climate
at 650,000 cells (the energy balance, the sea surface; the moisture after
M4-14). The table of `ExoPlaSim との使い分け` stays ("全球格子と同じ").

## Acceptance criteria

- [ ] Earth's energy balance (`monthly_temperature`, the Earth-parameter
      world of seed 1) reaches its periodic state at 650,000 cells within
      10 minutes, and the sea surface with currents
      (`monthly_sea_surface_temperature`) within 15 minutes (release, this
      machine), recorded with the time per model year at 10,000 to 650,000
      cells; the twelve potentials of the moisture within 1 minute.
- [ ] The multigrid's iterations per solve stay within 20 for the Earth
      sub-step and within 30 over the v1 range of `D` at every size
      (`iterations_do_not_grow_with_the_grid`,
      `robust_over_the_planet_range`, and the timing tool at 160,000 and
      650,000 cells); its solutions match dense solutions
      (`solves_match_a_dense_solution`).
- [ ] At 160,000 cells the zonal means (10° bands) and seasonal ranges of
      the periodic states of the energy balance and of the sea surface with
      currents agree with those of the Jacobi solvers within 2 × 10⁻³ K;
      every existing climate test passes with its tolerances
      (`multigrid_matches_jacobi` per module).
- [ ] Images of `out/M4-09/` viewed and described; golden values recorded
      again and the change stated as world-breaking for the climate
      modules' outputs (no pipeline world changes); the Erosion stage's
      golden values unchanged.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test --release -p kuniumi-sim multigrid
cargo test --release -p kuniumi-sim climate
```

Timing and images from the scratch tool, in `out/M4-09/`.

## Result

What was done:
- `kuniumi-sim::multigrid` (new): smoothed-aggregation algebraic
  multigrid (Vaněk, Mandel & Brezina 1996) for the M-matrices of the
  climate model: `STRENGTH_THRESHOLD` (0.08, halved per level),
  `COARSEST_CELLS` (400), `COARSENING_LIMIT` (0.8), `COARSEST_SWEEPS` (4);
  `Multigrid` (`new` for a symmetric matrix, checked bit for bit;
  `with_prolongators_of` for a non-symmetric matrix with the aggregates
  and prolongators of a symmetric hierarchy; `matrix`, `len`, `is_empty`,
  `levels`, `level_sizes`, `operator_complexity`, `apply`);
  `conjugate_gradient` and `bicgstab` preconditioned by one V-cycle. The
  aggregation (three passes in index order), the Gershgorin-damped
  prolongator smoothing, the Galerkin products (sparse accumulation, rows
  sorted as integers), the Gauss–Seidel V-cycle and the dense LU with
  partial pivoting of the coarsest level are crate-internal; everything
  sequential in index order. The module documentation gives the method,
  its source, the non-symmetric use and the determinism. 25 tests in
  `multigrid/tests.rs`.
- `operators::Laplacian::matrix(mass, coupling)`: the matrix of `solve`
  as a `SparseMatrix` (symmetric bit for bit). `Laplacian::solve` and
  `implicit_step` are unchanged (the Erosion stage's golden values pass).
- `sparse.rs`: the BiCGSTAB loop moved into the crate-internal
  `bicgstab_with`, which takes the preconditioner; `bicgstab` calls it with
  the Jacobi preconditioner, the same operations in the same order (the
  golden values of M4-04 to M4-06 passed unchanged after this step, before
  the climate models changed solver); `check_system` names its solver.
- `energy_balance.rs`: `EnergyBalance::new` builds the multigrid of the
  sub-step matrix `a (C/Δt + B) + D L` once (the albedo enters only the
  right-hand side); `warm_start` builds that of `B a + D L` once for all
  its solves, `correct` that of `B' a + D L` each year
  (`diffusion_multigrid`); the documentation of the time stepping, the
  limitations (with the measured cost) and of `monthly_temperature`, and
  M4-15 for the solver of the periodic state.
- `sea_surface.rs`: one hierarchy for the twelve months' sub-steps and
  return vectors, of the matrix with the annual-mean upwelling with the
  prolongators of the energy balance's sub-step hierarchy; the warm start
  and the correction build theirs from `B a + D L` and `B' a + D L`; the
  V-cycle smooths each month's own matrix on the finest level (below).
- `moisture.rs`: the hierarchy of the potential's matrix `10⁻⁹ a + L` is
  built once in `MoistureModel::new`; the upwind advection keeps
  `sparse::bicgstab` (M4-14).
- Documentation of `operators`, `sparse`, `stages/climate.rs`,
  `test_worlds.rs` and the test modules; `docs/design-ja.md`
  (実装ロードマップ / 計算コストの目安: the measured cost of the climate at
  650,000 cells); `CHANGELOG.md` (Added, Changed, World-breaking). The
  split approved by the owner when the ticket was refined: M4-14 and M4-15
  (draft) in `docs/tickets/README.md`, M4-07 depends on M4-14 too.

Deviations from the plan:
- The V-cycle smooths the finest level with the matrix being solved, not
  the hierarchy's own finest matrix (`cycle(work, finest, …)`): with one
  hierarchy for twelve months the sea surface's preconditioner then
  differs from each month's matrix only on the coarser levels. Found by
  `upwelling_cools_its_cell_and_returns_the_heat_by_area`, whose system is
  diagonal (no diffusion, no currents) and was solved exactly by the
  Jacobi preconditioner: with the annual-mean diagonal in the smoother the
  solve stopped at its tolerance 7 × 10⁻⁷ K from the closed form, against
  the test's 10⁻⁹ K; with the month's own matrix it is exact again and the
  test passes unchanged. Covered on several levels by
  `the_cycle_smooths_the_finest_level_with_the_matrix_being_solved` (the
  cycle bit for bit against its parts composed in the test; added in
  review round 1) and on one level by
  `a_hierarchy_serves_a_matrix_close_to_its_own`. A hierarchy of a single
  level of at most 400 cells is solved exactly, so its cycle is the
  inverse of its own matrix whatever matrix is being solved (documented
  with the solvers).
- The coarsest level is solved by LU with partial pivoting for symmetric
  and non-symmetric hierarchies alike (the plan had Cholesky for the
  symmetric one): one code path; the cycle stays symmetric to rounding. A
  level that cannot coarsen and keeps more than 400 cells (masses
  dominating, as with a coupling of 10⁻⁹) is solved by 4 pairs of
  Gauss–Seidel sweeps (`COARSEST_SWEEPS`).
- The symmetry of the Galerkin product is checked to 10⁻¹⁴ of the largest
  entry, not bit for bit as planned: `Pᵀ (A P)` sums the two mirrored
  entries in different orders.
- `multigrid_matches_jacobi` became one test per module that compares a
  single solve with the Jacobi solver of before on the same system
  (`a_substep_matches_the_jacobi_solver` in the energy balance and the sea
  surface, `the_potential_matches_the_jacobi_solver` in the moisture);
  the periodic states are compared with main's Jacobi build by the scratch
  tool at 10,000 to 160,000 cells (criterion 3), which an in-repository
  test could not do without keeping the old path.
- `solves_match_a_dense_solution` runs on random matrices of up to 900
  rows and the grid's system on 1,000 cells, within 10⁻⁹ of the largest
  value, not on 2,000 cells within 10⁻¹² as planned: the dense
  elimination of the test on 2,000 rows is slow in debug builds, and
  10⁻¹² is below what a relative residual of 10⁻¹³ guarantees for these
  condition numbers.
- The planned tests of the multigrid BiCGSTAB's breakdowns and restarts
  are those of `sparse::bicgstab`: both run the same code
  (`sparse::bicgstab_with`) with a different preconditioner. `the_stops`
  checks its zero right-hand side, exact start and iteration limit (the
  first two added in review round 1).
- The solvers check that the matrix's diagonal is finite and positive, as
  the Gauss–Seidel sweeps divide by it, and the moisture builds the
  potential's hierarchy only where `δ` < 1 (review round 1).
- Tests beyond the plan: `laplacian_matrix_is_the_system_of_solve`,
  `the_hierarchy_is_reported`, `a_level_that_cannot_coarsen_is_smoothed`,
  `a_hierarchy_serves_a_matrix_close_to_its_own`, `the_stops`,
  `the_cycle_smooths_the_finest_level_with_the_matrix_being_solved`,
  `conjugate_gradient_needs_a_positive_diagonal`.
- The API differs from the plan's sketch: `Multigrid::new` and
  `with_prolongators_of` take the matrix by value (the hierarchy keeps it
  as its finest level, so the energy balance stores it once), and
  `with_prolongators_of` takes the symmetric `&Multigrid` rather than its
  matrix (its prolongators are reused, not built again); the
  preconditioned BiCGSTAB is `multigrid::bicgstab` over the crate-internal
  `sparse::bicgstab_with`, not `sparse::bicgstab_multigrid`; the
  rectangular product and transpose are the crate-internal `Rows` of
  `multigrid.rs`, not functions of `sparse.rs`.
- Size: about 3,100 lines added in all (the multigrid with its
  documentation 1,030, its tests 880, the climate models and their tests
  about 500, documents the rest), against the estimated 1,400: rustfmt's
  layout, the module documentation and the tests added in review.

Evidence:

1. Criterion 1 (release, this machine, the scratch tool
   `<scratchpad>/m4_09_tool`, a Cargo package with path dependencies, the
   real API; Earth's parameters, seed 1, relief of the Crust stage; each
   size run alone; the "Jacobi" column is the same tool built against
   main):

   | Cells | Energy balance, Jacobi | multigrid | Sea surface, Jacobi | multigrid | Twelve months of moisture, Jacobi | multigrid |
   | --- | --- | --- | --- | --- | --- | --- |
   | 10,000 | 4.6 s (10 years) | 2.4 s (10) | 5.4 s (9) | 2.8 s (9) | 0.52 s | 0.12 s |
   | 40,000 | 45.5 s (12) | 11.5 s (12) | 60.5 s (12) | 14.0 s (12) | 4.22 s | 0.54 s |
   | 160,000 | 312.8 s (8) | 38.6 s (8) | 433.0 s (8) | 47.7 s (8) | 39.9 s | 2.47 s |
   | 650,000 | about 55 min (277 s per year, M4-02) | 257.3 s (11) | about 2 h (estimated) | 318.7 s (12) | about 9 min (estimated) | 12.0 s |

   Per model year with the multigrid: the energy balance 0.24, 0.96, 4.82
   and 23.4 s, the sea surface 0.31, 1.17, 5.97 and 26.6 s. At 650,000
   cells Earth's energy balance is periodic in 4.3 minutes (target 10),
   the sea surface with currents in 5.3 minutes (target 15), the twelve
   months of the moisture (with their potentials, at most 29 iterations)
   in 12.0 s (target 1 minute); the currents of M4-04 take 39 s more. The
   650,000-cell Earth: mean surface temperature 288.00 K, `P₂` −26.83 K,
   every 10° band of the annual sea-level temperature within 3.30 K of
   North & Coakley's fit; with currents 288.23 K and −26.08 K.
2. Criterion 2. `iterations_do_not_grow_with_the_grid`: 12, 14 and 15
   iterations on 2,500, 10,000 and 40,000 cells; the scratch tool's
   `iterations` mode (sub-steps from the warm start, with `mise run ci`
   running alongside): 160,000 cells 11–12 (Earth's `D`), 13–16
   (990 `D⊕`), 8 (`D⊕/91`); 650,000 cells 13–14, 16–20 and 9.
   `robust_over_the_planet_range` (10,000 cells): 7, 14 and 17 for
   `D⊕/91`, `D⊕` and 990 `D⊕`, 15 for the warm start, 18 for the
   potential. `solves_match_a_dense_solution`: random symmetric M-matrices
   of 5 to 900 rows and the Earth-like sub-step on 1,000 cells, with and
   without advection, within 10⁻⁹ of the dense solutions.
   `the_prolongator_is_the_smoothed_aggregation`,
   `coarse_operators_are_galerkin`,
   `aggregates_cover_the_strongly_connected_cells`,
   `the_v_cycle_is_symmetric` (to 10⁻¹²), `repeatable_bit_for_bit`,
   `the_cycle_smooths_the_finest_level_with_the_matrix_being_solved`,
   `small_systems_are_solved_by_the_coarsest_level` (one iteration),
   `a_level_that_cannot_coarsen_is_smoothed`,
   `a_hierarchy_serves_a_matrix_close_to_its_own` (24 and 16 iterations
   with a diagonal up to 10 % off), `the_hierarchy_is_reported` (10,000
   cells: levels of 10,000, 1,651 and 384 cells, operator complexity
   1.41) and the panics pass.
3. Criterion 3 (the scratch tool's `compare`, this branch against main's
   Jacobi build on the same worlds; the same number of years at every
   size): on 160,000 cells the energy balance's monthly means differ by at
   most 3.7 × 10⁻⁷ K, the annual means of the 10° bands by 6.1 × 10⁻⁸ K,
   their mean seasonal ranges by 4.7 × 10⁻⁸ K; the sea surface with
   currents 3.9 × 10⁻⁶, 3.0 × 10⁻⁷ and 1.2 × 10⁻⁶ K (40,000 cells:
   7.5 × 10⁻⁷, 5.2 × 10⁻⁸, 4.3 × 10⁻⁸ and 5.0 × 10⁻⁶, 4.8 × 10⁻⁷,
   1.1 × 10⁻⁶ K; 10,000: 6.9 × 10⁻⁷, 1.6 × 10⁻⁷, 2.2 × 10⁻⁷ and
   3.5 × 10⁻⁶, 2.6 × 10⁻⁷, 9.1 × 10⁻⁷ K): all far within 2 × 10⁻³ K.
   `a_substep_matches_the_jacobi_solver` (10,000 cells): the energy
   balance 8 iterations against 78, 3.3 × 10⁻⁷ K apart, and the warm
   start's first solve 16 against 297, 3.6 × 10⁻⁹ K (added in review
   round 1); the sea surface 5 against 51, 1.4 × 10⁻⁶ K. `the_potential_matches_the_jacobi_solver`:
   20 iterations against 424, every flow within 1.6 × 10⁻¹³ of the
   largest. Every existing climate test passes with its tolerances, in
   debug and release.
4. Criterion 4. Images in `out/M4-09/` (scratch tool; 650,000 cells,
   Earth, seed 1), viewed:
   - `surface_annual.png`, `surface_m01.png`, `surface_m07.png` (RdBu
     centred at 273.15 K, 213 to 333 K): zonal bands, warmest along the
     equator; the northern continents at 40–60° N pale in January and
     warm in July, the large southern continent at 30–60° S blue in July
     and pale in January; the mountain belts along the coasts are narrow
     cold ribbons a few cells wide, the lapse rate on each 30-km cell's own
     height; polar land cold all year.
   - `seasonal_range.png` (0 to 40 K): smooth, largest (above 20 K) over
     the interior of the large southern continent at about 40° S, then the
     northern continents' interiors, smallest along the equator over the
     ocean; the continents' signal spreads about 1,000 km over the sea;
     no cell-scale noise.
   - `sst_annual.png` (ocean only, RdBu centred at 287 K, 268 to 306 K):
     zonal bands with the warm pool across the tropics, warmest in the
     western parts of the basins, cooler water along some eastern margins,
     below 268 K poleward of about 60° (no sea ice yet, M4-11).
   - `difference_annual_160k.png` (multigrid − Jacobi, annual sea-level
     temperature on 160,000 cells, ±8.2 × 10⁻⁸ K): a smooth pattern on
     the scale of the oceans, no structure at the coasts or at the scale of
     the cells: the smooth error that both solvers leave at their
     tolerance.
   Golden values recorded again after the tests above passed and the
   images were viewed (10,000 cells, seed 1; the same in debug and
   release): `golden_fingerprint_of_earth_temperature`
   `0xB43F_20A2_5237_116E` (was `0x9FAC_D213_E991_ED31`),
   `golden_fingerprint_of_earth_circulation` `0xD3B8_67CA_7631_CB7E`
   (was `0xE75E_0007_E628_9097`), `golden_fingerprint_of_earth_currents`
   `0x9195_CFEF_7FA2_AF8C` (was `0x043D_1F09_12AD_6477`),
   `golden_fingerprint_of_earth_sea_surface` `0x8DAD_2930_A815_9B7D` (was
   `0xD958_B31C_77E8_4D41`), `golden_fingerprint_of_earth_moisture`
   `0xFE04_B92F_C2B5_4AAA` (was `0xBB8E_6949_41D8_448D`); the insolation's
   and every other stage's golden values (the Erosion stage's among them)
   unchanged. Stated as world-breaking for the climate modules' outputs in
   `CHANGELOG.md`; no world of the pipeline changes.
5. `mise run ci` passes: fmt, clippy `-D warnings`, rustdoc, the browser
   build's lints, 1,957 tests summed over debug and release, MSRV 1.85.0,
   cargo-deny (advisories, bans, licences and sources ok).

Review:
- Round 1 (APPROVE; the reviewer re-ran `mise run ci`, its own copy of
  the scratch tool at 10,000 to 650,000 cells against main's Jacobi build,
  bit-identical outputs at 650,000 cells, the images, the golden values and
  mutants of the code). R1-1 minor: no test checked that
  `with_prolongators_of` forms its coarse operators from the non-symmetric
  matrix (a mutant using the symmetric one passed):
  `solves_match_a_dense_solution` now compares its first coarse level
  with `Pᵀ A P` of the non-symmetric matrix and checks that it is not
  symmetric (the mutant fails, checked). R1-2 minor: only the one-level
  path of the finest level's smoothing with the matrix being solved was
  tested: `the_cycle_smooths_the_finest_level_with_the_matrix_being_solved`
  added (the mutant fails, checked), and the documentation of the solvers
  says that a dense single level ignores the matrix being solved. R1-3
  minor: the warm start was claimed but not compared with the Jacobi
  solver: its first solve now is. R1-4 minor: the multigrid BiCGSTAB's
  zero right-hand side and exact start added to `the_stops`; the
  deviations from the planned tests recorded above. R1-5 nit: the
  changelog's iterations at 650,000 cells are the implementation's
  (13–14), its Changed entry appended. R1-6 nit: the Crates line names
  `operators`. R1-7 nit: the solvers check the matrix's diagonal, with
  `conjugate_gradient_needs_a_positive_diagonal`. R1-8 nit: the potential's
  hierarchy is built only where `δ` < 1. `mise run ci` ran again after the
  fixes (evidence 5).
- Round 2 (APPROVE; the reviewer re-ran `mise run ci` and the MSRV check
  from scratch, the scratch tool at 650,000 cells on an idle machine
  (255.8 s, 316.6 s and 11.85 s; outputs bit-identical to round 1's, so the
  fixes changed no number), the comparison at 160,000 cells, the images,
  the golden values, every round-1 mutant and four more, all caught; every
  round-1 finding resolved). R2-1 nit: the size above updated after the
  tests added in round 1. R2-2 nit: the API's departures from the plan's
  sketch recorded above.

Follow-ups: M4-14 (the moisture's advection at 650,000 cells), M4-15 (the
periodic state of slowly converging planets), both from the split approved
by the owner.

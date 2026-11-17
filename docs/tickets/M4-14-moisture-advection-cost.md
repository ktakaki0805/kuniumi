# M4-14 The moisture's advection at 650,000 cells

- Status: done
- Milestone: M4
- Depends on: M4-09
- Crates: kuniumi-sim (sparse, stages/climate)
- Design: 実装ロードマップ / 計算コストの目安
- Decisions: D-08

## Goal

The moisture of M4-06 reaches its periodic state at the default 650,000
cells in minutes rather than about 45 minutes, before M4-07 wires the
Climate stage.

## Scope

In:
- Measured in M4-09's prototype (release, the Earth-parameter world of
  seed 1, one solve of a sub-step's upwind system by BiCGSTAB with the
  Jacobi preconditioner): 157 iterations and 0.46 s at 160,000 cells, 305
  and 3.98 s at 650,000, with Courant numbers up to 2,100 and 8,450 per
  sub-step of 7.6 days; with about 4.6 solves per sub-step (the active set
  of the precipitation threshold, M4-06) about 15 minutes per model year at
  650,000 cells and 45 for the periodic state of 3 years. The twelve
  potentials of the divergent wind are M4-09's.
- Preconditioners tried there: M4-09's multigrid from the symmetric part
  does not suit the advection (84 iterations and slower at 160,000 cells,
  no convergence at 650,000); symmetric Gauss–Seidel in index order 68
  iterations and 2.16 s at 650,000 cells; one Gauss–Seidel sweep in a
  downwind order (Kahn's algorithm on the outflow graph, cycles broken at
  the cell of least remaining inflow) 41 and 2.75 s.
- The choice below ("Choice", chosen with the owner): a downwind
  order of each month's cells, nonlinear Gauss–Seidel sweeps of each
  sub-step in that order, then today's active set method with BiCGSTAB
  preconditioned by a Gauss–Seidel sweep in that order; the same discrete
  equations and the same stops, so the results within the solver's
  tolerance of today's.
- Timing at 10,000 to 650,000 cells recorded; results compared with the
  Jacobi solver at 160,000 cells.

Out:
- The physics of the moisture (M4-06, M4-13), its time stepping and its
  periodic-state test (3 years, the third confirming the second).
- The potential of the divergent wind (M4-09), the Stommel solve of M4-04
  (`sparse::bicgstab` with the Jacobi preconditioner stays for it).
- Parallel loops (`rayon`, D-10): the sweeps are sequential.

## Implementation notes

The ticket was `draft` and left the method and its time target to the
owner. A prototype outside the repository (`<scratchpad>/m4_14_proto`, a
Cargo package with path dependencies on a snapshot of `main`, as for M4-01
to M4-09) reimplements `MoistureModel::substep` and `periodic_state` over
the public accessors of `MoistureMonth`, with pluggable solvers; with the
Jacobi solver it reproduces the API bit for bit (10,000 cells). Inputs as
the tests build them: the Earth-parameter world of seed 1 (relief of the
Crust stage), the winds (M4-03) of the energy balance's periodic state
(M4-02) and the surface temperatures of the sea surface with currents
(M4-05), cached once per size. Release, this machine.

### Measurements (prototype)

The periodic state (3 years, 144 sub-steps), without the twelve months'
set-up (the potentials of M4-09: 0.13, 0.6, 2.8 and 12.3 s; the downwind
orders below add 0.05, 0.26, 1.9 and 9.4 s in the prototype). Each cell
gives the time, then the linear solves per sub-step and the iterations per
solve:

| Method | 10,000 | 40,000 | 160,000 | 650,000 |
| --- | --- | --- | --- | --- |
| Today: BiCGSTAB, Jacobi (the API) | 3.0 s; 3.3 solves a sub-step, 35 iterations a solve | 29.6 s | 331 s (other runs alongside; M4-06: 303 s) | about 45 min (M4-09) |
| BiCGSTAB, one Gauss–Seidel sweep in the downwind order | 0.8 s; 3.3, 6.8 | 5.9 s; 4.0, 9.4 | 50.8 s; 5.1, 12.9 | 310 s; 6.0, 17.5 |
| the same, ILU(0) in the downwind order | 1.0 s; 3.3, 5.4 | 6.3 s; 4.0, 7.4 | | |
| the same, symmetric Gauss–Seidel in the downwind order | 0.8 s; 3.3, 5.9 | 6.4 s; 4.0, 8.2 | | |
| the same, symmetric Gauss–Seidel in index order | 1.4 s; 3.3, 10.4 | | | |
| downwind Gauss–Seidel, the solves until the set repeats at 10⁻⁶ | | 3.1 s; 5.0, 3.6 | 22.9 s; 6.2, 4.4 | 147 s; 7.2, 6.2 |
| **40 nonlinear sweeps, then downwind Gauss–Seidel** | 0.4 s; 1.0, 0 | 2.0 s; 1.0, 0 | 10.2 s; 1.0, 1.3 | **72.9 s**; 1.2, 9.0 |

- The downwind order (below) leaves an upper part (entries towards cells
  later in the order: the closed streamlines) of 6 % of a row's diagonal
  on average and at most 40 to 46 % (10,000 to 650,000 cells, January),
  against 36 % and 108 % in index order (10,000 cells).
- The nonlinear sweeps converge by a factor of about 0.6 (160,000 cells) to
  0.8 (650,000) per sweep, and the active set settles within 10 to 30
  sweeps. Sweeps at 650,000 cells: 30 give 88.9 s (1.7 solves a sub-step),
  40 give 72.9 s, 50 give 67.9 s, 60 give 68.1 s; at 160,000 cells 20, 30,
  40 and 60 give 11.5, 10.1, 10.2 and 13.3 s; at 10,000 and 40,000 cells
  40 sweeps reach the tolerance themselves (the solve that follows stops at
  once). A sweep costs about 7 ms at 650,000 cells, a BiCGSTAB iteration
  about 21 ms.
- Every method reached its tolerance and every active set settled; each
  periodic state took 3 years. The differences from the Jacobi solver's
  periodic state are those of the solver's tolerance: at 160,000 cells
  (40 sweeps) every monthly-mean vapour within 5.0 × 10⁻⁹ kg/m², every
  monthly precipitation within 8.6 × 10⁻⁹ mm/day and evaporation within
  5.1 × 10⁻⁹ mm/day (the Jacobi run's own last year-to-year change
  6.4 × 10⁻⁹ kg/m²); at 10,000 and 40,000 cells 1.1 and 1.9 × 10⁻⁹ kg/m²;
  at 650,000 cells the sweeps against the plain downwind Gauss–Seidel
  4.6 × 10⁻⁹ kg/m².

### Choice (chosen with the owner)

The last row, chosen with the owner when the ticket was refined, with the
target: the moisture's periodic state with its twelve months at 650,000
cells **within 3 minutes** (prototype about 95 s: the potentials 12 s, the
orders 9 s, the periodic state 73 s). Not chosen:

- Downwind Gauss–Seidel alone (about 5.5 minutes at 650,000 cells) or with
  loose solves until the set repeats (about 2.8 minutes): each sub-step
  still solves the linear system 6 or 7 times.
- ILU(0) in the downwind order: fewer iterations than one Gauss–Seidel
  sweep (7.4 against 9.4 at 40,000 cells), slower, as its factors change
  with the active set every solve.
- Another time stepping of the vapour (semi-Lagrangian, explicit
  sub-cycling): the results would change beyond the solver's tolerance.

### The downwind order (`sparse::downwind_order`)

For an M-matrix whose entries off the diagonal are the upwind inflows
(`a_ij < 0`: cell `i` receives from cell `j`), Kahn's topological sort of
the inflow graph (A. B. Kahn 1962, "Topological sorting of large
networks", Commun. ACM 5, 558–562), each cell after the cells that flow
into it: the "downwind numbering" of W. Hackbusch & T. Probst (1997,
"Downwind Gauß–Seidel smoothing for convection dominated problems", Numer.
Linear Algebra Appl. 4, 85–102) and J. Bey & G. Wittum (1997, "Downwind
numbering: robust multigrid for convection-diffusion problems", Appl.
Numer. Math. 23, 177–192), in which a Gauss–Seidel sweep solves an
acyclic upwind system exactly. Where no cell is ready (closed
streamlines), the next is the unordered cell whose inflow from unordered
cells is the smallest share of its diagonal, `Σ_j (−a_ij) / a_ii`; ties
by index. A binary heap keyed by `(ready, the share's bits, index)` with
lazy deletion; the shares are non-negative, so their bits order as their
values (docs/conventions.md, rule 5); everything sequential.

### The solve (`sparse::bicgstab_gauss_seidel`)

BiCGSTAB as `sparse::bicgstab` (the same `bicgstab_with`, the same stops,
breakdowns and restarts), right-preconditioned by one forward Gauss–Seidel
sweep `(D + L)⁻¹` in row order, for a matrix renumbered by the order
(`SparseMatrix::permuted`, each row's columns sorted, so its lower part is
a prefix of the row).

### The moisture (`moisture.rs`)

- `MoistureModel::month` builds the order of the month's sub-step matrix
  without the saturation excess (the advection's outflows plus
  `a_i/Δt`, the lift and the ocean's exchange: the same in every sub-step
  of the month) and keeps the advection renumbered in it (its inflows
  only; the zero entries of the outflow edges, which today's matrix
  carries, are dropped). Same memory as today's matrix.
- `MoistureModel::substep`, from the previous vapour: first
  `NONLINEAR_SWEEPS` = 40 nonlinear Gauss–Seidel sweeps of the sub-step's
  equations in the order: each cell solves its own piecewise-linear
  equation with its neighbours' latest values,
  `s = b⁰_i − Σ_j a_ij W'_j`, `W'_i = s/d⁰_i` where that is at most the
  threshold `r_c W_s`, else `(s + k_i r_c W_s)/(d⁰_i + k_i)`
  (`d⁰`, `b⁰` the row without the excess, `k_i = a_i keep/τ_c`; both
  branches agree at the threshold, so the update is continuous and
  monotone). The matrix is an M-matrix whose columns are strictly
  diagonally dominant (each column's entries off the diagonal are the
  outflows that its diagonal holds, which also holds `a_i/Δt` and the
  sinks), and the excess an isotone function of each cell's own vapour:
  the system is a continuous, surjective M-function, for which the
  nonlinear Gauss–Seidel iteration converges in any order and from any
  start to the unique solution (W. C. Rheinboldt 1970, "On M-functions
  and their application to nonlinear Gauss–Seidel iterations and to
  network flows", J. Math. Anal. Appl. 32, 274–307). Then the active set method of today,
  started from the sweeps' set and vapour, each linear solve by
  `bicgstab_gauss_seidel` on the renumbered system to `SOLVER_TOLERANCE`
  (10⁻¹²), until the set repeats, at most `MAX_ACTIVE_SET_ITERATIONS`.
  The sweeps only choose the start: the result is the solution of the same
  system with the same stops as today.
- `MoistureModel::substep` no longer takes the active set of the previous
  sub-step (the sweeps choose it): `substep(month, vapour)`, with the set
  of its last solve in `MoistureStep::active`.
- `MoistureStep::active_set_iterations` keeps counting the linear solves;
  the documentation of the time stepping, of the limitations (cost) and of
  `stages/climate.rs` and `sparse.rs` say what changed.

### Tests

Reference values from Gaussian elimination in the tests, the Jacobi solver
of the code (`sparse::bicgstab`) and closed forms.

- `sparse`: `the_downwind_order_is_a_permutation` (random matrices);
  `an_acyclic_flow_is_ordered_downwind` (a random acyclic flow: every cell
  after the cells that flow into it, and `bicgstab_gauss_seidel` solves it
  in one iteration); `a_cycle_is_broken_at_its_weakest_inflow` (a ring of
  cells with one weak edge); `the_permuted_matrix_is_the_same_operator`;
  `bicgstab_gauss_seidel_matches_a_dense_solution` (random non-symmetric
  M-matrices of 5 to 400 rows, within 10⁻¹⁰ of the largest value); its
  stops (zero right-hand side, exact start, iteration limit) and panics;
  `bicgstab`'s tests unchanged.
- `moisture`: `the_sweeps_solve_a_cell_exactly` (no transport: one sweep
  gives the closed form on both sides of the threshold);
  `a_substep_matches_the_jacobi_solver` (10,000 cells, Earth: the sub-step
  against the active set method with `sparse::bicgstab` on the same
  system, built in the test: the same set, the vapour within 10⁻⁹ of its
  largest value); every existing test unchanged in its tolerances (the
  budget, the transport's conservation and monotonicity, the Earth
  criteria, the extreme planets, every solve converged and every set
  settled).
- Golden values: `golden_fingerprint_of_earth_moisture` changes, as the
  solution changes within the solver's tolerance; recorded again after the
  tests above and the images. The pipeline's worlds do not change (the
  Climate stage is not wired yet), nor any other golden value; stated in
  the commit's `World-breaking:` line and in `CHANGELOG.md`.

### Timing and images

A scratch tool (`<scratchpad>/m4_14_tool`, the real API) times the twelve
months and the periodic state of the moisture at 10,000, 40,000, 160,000
and 650,000 cells, each size alone, release; at 160,000 cells also with
the build of `main` (Jacobi), and compares the periodic states (largest
differences of the monthly vapour, precipitation and evaporation); it
measures main's 650,000 cells once (about 45 minutes). It renders into
`out/M4-14/`: the annual precipitation at 650,000 cells (√(mm/day), as
M4-06's maps), January's downwind order at 160,000 cells (a cell's position
in the order over `n`), and the difference of the annual precipitation
from the Jacobi solver's at 160,000 cells.

### Size

About 900 lines: `sparse.rs` (order, renumbering, the preconditioned
solve) with its documentation about 250 and its tests about 300; the
moisture with its documentation about 200 and its tests about 100;
documents the rest. Slightly over the guideline of 800; one ticket, as the
order serves only these solves and the moisture's tests are what prove
them.

### Design document

`実装ロードマップ / 計算コストの目安`: the measured cost of the moisture's
periodic state at 650,000 cells replaces "水蒸気の風上移流（約45分）は M4-14
で扱います".

## Acceptance criteria

- [x] At 650,000 cells (Earth, seed 1, release, this machine) the
      moisture's twelve months (their potentials and downwind orders) and
      its periodic state take at most 3 minutes together (prototype about
      95 s), recorded with the times at 10,000, 40,000 and 160,000 cells
      and the Jacobi solver's at 160,000 and 650,000 cells.
- [x] At 160,000 cells the periodic state takes the same years as the
      Jacobi solver's and every monthly-mean vapour lies within 10⁻⁶ kg/m²
      (`PERIODIC_TOLERANCE_KG_M2`), every monthly precipitation and
      evaporation within 10⁻⁶ mm/day of it (prototype: 5.0 × 10⁻⁹ kg/m²,
      8.6 × 10⁻⁹ and 5.1 × 10⁻⁹ mm/day), documented;
      `a_substep_matches_the_jacobi_solver`, the solver's tests
      (`bicgstab_gauss_seidel_matches_a_dense_solution`,
      `an_acyclic_flow_is_ordered_downwind`, …) and every existing moisture
      test pass with their tolerances.
- [x] Images of `out/M4-14/` viewed and described; the moisture's golden
      value recorded again and stated as world-breaking for the moisture
      module's outputs (no pipeline world changes); every other golden
      value unchanged.
- [x] `mise run ci` passes.

## Verification

```sh
cargo test --release -p kuniumi-sim sparse
cargo test --release -p kuniumi-sim moisture
cargo test -p kuniumi-sim moisture
```

Timing, the comparison and the images from the scratch tool, in
`out/M4-14/`.

## Result

What was done:
- `kuniumi-sim::sparse`: `downwind_order` (Kahn's topological sort of the
  inflow graph, the ready row of the lowest index first, closed
  streamlines cut at the unordered row whose inflow from unordered rows is
  the smallest share of its diagonal, ties by index; a binary heap of
  integers and the bits of non-negative shares with lazy deletion),
  `SparseMatrix::permuted` (row and column `q` are row and column
  `order[q]`, each row's entries bit for bit with sorted columns) and
  `bicgstab_gauss_seidel` (`bicgstab_with` with one forward Gauss–Seidel
  sweep `(D + L)⁻¹` in row order as the right preconditioner). The module
  documentation gives the downwind numbering with its sources. 12 tests
  added to `sparse/tests.rs` (seven and five panics); `bicgstab`'s
  unchanged.
- `moisture.rs`: `MoistureModel::month` builds the advection of its
  inflows only (the zero entries of the outflow edges dropped), the
  downwind order of the sub-step matrix without the excess and the
  advection renumbered in it (`MoistureMonth::downwind_order`);
  `MoistureModel::substep(month, vapour)` runs `NONLINEAR_SWEEPS` (40)
  nonlinear Gauss–Seidel sweeps in that order, then the active set method
  from their result with every linear solve by `bicgstab_gauss_seidel` on
  the renumbered system to `SOLVER_TOLERANCE`, and returns the set of its
  last solve in `MoistureStep::active` (the parameter `active` is gone:
  the sweeps choose the set; Implementation notes); `periodic_state`
  follows. The rows of a sub-step in the downwind order and one nonlinear
  sweep are the private `MoistureModel::rows` and `nonlinear_sweep`; the
  coefficient of a cell's own vapour besides the outflows and the excess
  is `MoistureModel::own_coefficient`, shared by the month's order and the
  sub-steps (review round 1). The
  module documentation (time stepping with the M-function argument and
  its source, the comparison with the Jacobi solver, a new "Cost"
  section), `SOLVER_TOLERANCE` and `SOLVER_MAX_ITERATIONS` updated. Two
  tests added to `moisture/tests.rs`, four call sites follow the new
  signature, the golden value recorded again.
- `stages/climate.rs` and `multigrid.rs` say how the moisture's advection
  is solved; `docs/design-ja.md` (実装ロードマップ / 計算コストの目安: the
  measured cost replaces "約45分 … M4-14 で扱います"); `CHANGELOG.md`
  (Added, Changed, World-breaking).

Deviations from the plan:
- `bicgstab_gauss_seidel_matches_a_dense_solution` holds the strictly
  dominant random matrices to 10⁻¹⁰ of the largest value as planned and
  the weakly dominant, nearly singular ones to 10⁻⁸, as
  `bicgstab_matches_a_dense_solution` does (the same matrices); it also
  checks that it needs no more iterations than the Jacobi preconditioner,
  in index order and in the downwind order.
- `a_substep_matches_the_jacobi_solver` compares the four sub-steps of a
  month, not one, and the sets of raining cells exactly.
- The downwind order is rendered at 650,000 cells as well as at 160,000,
  and the annual precipitation at 160,000 cells as well as at 650,000.
- Size: about 1,100 lines without the ticket (`sparse.rs` and its tests
  560, the moisture and its tests 480, documents the rest) against the
  estimated 900: the tests of the order and the solver.

Evidence:

1. Criterion 1 (release, this machine, the scratch tool
   `<scratchpad>/m4_14_tool`, a Cargo package with path dependencies, the
   real API, built against the branch and against a snapshot of `main`;
   Earth's parameters, seed 1, relief of the Crust stage; the inputs of
   the tests, the winds of the energy balance's periodic state and the
   surface temperatures of the sea surface with currents, computed once
   per size and cached, the branch's bit-identical to them at 10,000 and
   40,000 cells; each size run alone):

   | Cells | Twelve months (potentials, orders) | Periodic state, this branch | Periodic state, Jacobi (main) |
   | --- | --- | --- | --- |
   | 10,000 | 0.15 s | 0.3 s (3 years) | 3.0 s (3 years) |
   | 40,000 | 0.70 s | 1.5 s (3) | 26.6 s (3) |
   | 160,000 | 3.4 s | 8.2 s (3) | 303.1 s (3) |
   | 650,000 | 18.7 s | 62.2 s (3) | 2,395 s, 40 minutes (3) |

   At 650,000 cells the model, its twelve months and the periodic state
   take 81.3 s together (target 3 minutes; main's 2,408 s, its twelve
   months 12.4 s); per model year 0.11, 0.50, 2.73 and 20.7 s against
   0.99, 8.9, 101 and 798 s. Every solve converged and every active set
   settled at every size. (Main's 650,000 cells ran once, `mise run ci`
   alongside for 3.5 of its 40 minutes; M4-09 had estimated 45.)
2. Criterion 2 (the scratch tool's `compare`, the branch against main's
   Jacobi build on the same inputs; both 3 years at every size): on
   160,000 cells every monthly-mean vapour within 5.0 × 10⁻⁹ kg/m²,
   precipitation within 8.6 × 10⁻⁹ mm/day and evaporation within
   5.1 × 10⁻⁹ mm/day (the Jacobi run's own last change 6.4 × 10⁻⁹ kg/m²);
   40,000 cells 1.9 × 10⁻⁹ kg/m², 1.5 × 10⁻⁹ and 8.9 × 10⁻¹⁰ mm/day;
   10,000 cells 1.1 × 10⁻⁹ kg/m², 2.6 × 10⁻⁹ and 1.6 × 10⁻⁹ mm/day;
   650,000 cells 4.9 × 10⁻⁹ kg/m², 1.4 × 10⁻⁸ and 8.4 × 10⁻⁹ mm/day: far
   within 10⁻⁶. The global mean precipitation agrees to the printed digits
   (2.5360, 2.8789, 2.6369 and 2.5835 mm/day). `a_substep_matches_the_jacobi_solver`
   (10,000 cells, the four sub-steps of February): the same sets, every
   vapour within 1.8 × 10⁻¹¹ of the largest (88.6 kg/m²), 4 solves against
   the Jacobi method's 14. `the_sweeps_solve_a_cell_exactly`: after one
   nonlinear sweep (`nonlinear_sweep`, private, tested directly) and after
   the sub-step every cell within 10⁻¹² of the closed form, 1,686 dry and
   1,648 raining land cells, 2,943 dry and 3,723 raining ocean cells, one
   solve; a mutant whose sweep never takes the branch above the threshold
   fails it (checked).
   `an_acyclic_flow_is_ordered_downwind`: every cell after its upwind
   cells, the renumbered matrix lower triangular, one iteration to the
   dense solution (10, 100 and 400 rows; 4, 7 and 8 in index order).
   `a_cycle_is_broken_at_its_weakest_inflow` and
   `a_cut_weighs_only_the_inflow_left_behind` (the orders by hand; the
   second added in review round 1, a mutant that weighs the whole inflow
   fails it, checked),
   `the_downwind_order_is_a_permutation`,
   `the_permuted_matrix_is_the_same_operator`,
   `bicgstab_gauss_seidel_matches_a_dense_solution`,
   `bicgstab_gauss_seidel_stops` and five panics pass. Every existing
   moisture test passes with its tolerances, in debug and release:
   `the_budget_closes_in_every_substep` (largest imbalance
   1.0 × 10⁻¹⁴ of the precipitation, one solve per sub-step),
   `transport_conserves_water_and_is_monotone` (drift 5.5 × 10⁻¹⁵ and
   1.3 × 10⁻¹², range 0.4952–0.5147 without the divergence),
   `earth_converges_to_a_periodic_state` (3 years, last change
   7.1 × 10⁻¹⁵ kg/m² on 40,000 cells), `extreme_planets_stay_finite`
   (every planet converged, settled and finite; the same precipitation as
   before to the printed digits), the Earth criteria of M4-06.
3. Criterion 3. Images in `out/M4-14/` (scratch tool), viewed:
   - `precipitation_annual_650000.png` (√(mm/day), 0–16 mm/day, as
     M4-06's maps): the equatorial rain band, dry subtropical belts at
     15–35° over the oceans and the continents' interiors, the storm
     tracks at 40–60° wettest on the west coasts and windward slopes,
     narrow coastal ribbons of orographic rain at the 30-km cells' scale,
     dry polar caps; streaks along the wind over the equatorial
     continents (the relief's lift carried downwind). The same patterns
     at 160,000 cells (`precipitation_annual_160000.png`, another world).
   - `downwind_order_m01_160000.png`, `downwind_order_m01_650000.png`
     (January, a cell's place in the order over the number of cells, 0
     to 1): smooth regions that grow along the wind, earliest under the
     divergent subtropical highs, latest in the convergent storm tracks
     and polar caps, with steps where closed streamlines are cut and
     spirals around the subpolar lows.
   - `difference_annual_160000.png` (branch − Jacobi, annual
     precipitation, ±6.9 × 10⁻¹⁰ mm/day): zero nearly everywhere, thin
     streaks along the southern westerlies and the northern polar flow,
     where the closed streamlines carry what the solvers leave at their
     tolerance; nothing at the coasts or at the cells' scale.
   `golden_fingerprint_of_earth_moisture` recorded again after the tests
   above passed and the images were viewed: `0x5B63_B12B_058C_535F` (was
   `0xFE04_B92F_C2B5_4AAA`; 10,000 cells, seed 1, 3 years, the same in
   debug and release); every other golden value unchanged (`mise run ci`).
   Stated as world-breaking for the moisture's outputs in `CHANGELOG.md`;
   no world of the pipeline changes.
4. `mise run ci` passes: fmt, clippy `-D warnings`, rustdoc, the browser
   build's lints, 1,985 tests summed over debug and release, MSRV 1.85.0,
   cargo-deny (advisories, bans, licences and sources ok).

Review:
- Round 1 (APPROVE; the reviewer re-ran `mise run ci`, its own copy of the
  scratch tool against the branch and against main at 10,000 to 160,000
  cells and the branch at 650,000 (80.2 s in all; outputs bit-identical
  to the ones above), the comparisons, the images (byte-identical), and
  mutants of the sweeps, the inflow's sign and the sorting, each caught).
  R1-1 minor: the module's determinism paragraph said that everything
  runs in index order: it now names the downwind order of the sweeps and
  solves and why it is deterministic. R1-2 minor: no test checked that a
  cut weighs the inflow from unordered rows only (a mutant that never
  lowered it passed every sparse test):
  `a_cut_weighs_only_the_inflow_left_behind` added (the mutant fails it,
  checked). R1-3 nit: the dry diagonal was written twice (the month's
  order and the sub-steps): one helper, `own_coefficient`. The fixes
  change no output (the scratch tool's outputs at 10,000 to 160,000 cells
  bit-identical, the golden value unchanged); `mise run ci` ran again
  (evidence 4).
- Round 2 (APPROVE, no findings; the reviewer re-ran `mise run ci`, its own
  copy of the scratch tool against the fixed branch at every size (650,000
  cells 81.8 s in all, every output bit-identical to round 1's) and main's
  at 160,000 cells (310.8 s), the comparisons, the images
  (byte-identical), the round-1 mutant of the cut (now caught by
  `a_cut_weighs_only_the_inflow_left_behind`) and four more mutants of
  `own_coefficient`, `keep`, the preconditioner and the order, each
  caught; every round-1 finding resolved).

Follow-ups: none.

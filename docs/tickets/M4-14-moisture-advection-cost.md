# M4-14 The moisture's advection at 650,000 cells

- Status: draft
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
- Choose, with the owner, a deterministic preconditioner or formulation
  (downwind-ordered Gauss–Seidel or ILU(0), fewer active-set solves, a
  better start of each solve, or another time stepping of the vapour) and
  its time target; the results within the solver's tolerance of today's.

Out:
- The physics of the moisture (M4-06, M4-13).

## Implementation notes

To be refined: the choice above (owner), and then its accuracy targets.

## Acceptance criteria

- [ ] The moisture's periodic state of Earth at 650,000 cells within a time
      agreed with the owner, recorded.
- [ ] Results within the solver's tolerance of the Jacobi solver's at
      160,000 cells, documented.
- [ ] `mise run ci` passes.

## Verification

Timing runs at 160,000 and 650,000 cells.

## Result

_(filled in when the ticket is implemented)_

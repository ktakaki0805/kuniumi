# M4-25 A robust solve of the currents of an ocean over the whole sphere

- Status: draft
- Milestone: M4
- Depends on: M4-10
- Crates: kuniumi-sim (stages/climate: ocean; the bounds tests of lib.rs)
- Design: 気候モデル / 構成要素
- Decisions: D-08, D-12

## Goal

The currents of every basin of the planets of the bounds tests reach the
solver's tolerance, so that corner 105 no longer warns of its currents.

## Scope

In:
- Found in M4-10 (1,000 cells, seed 3). Corner 105 (13,000 km, 8-hour
  days, a 3,650-day year, 400 W/m², all ocean, 10 bar) is one basin over
  the whole sphere. Its system is solved by unpreconditioned BiCGSTAB from
  zero (M4-04), which needed 2,872 iterations with the climate of the
  `P₂` form of the albedo; with the winds of the zenith angle law it
  stagnates at a relative residual of 2.1 × 10⁻³ after its 20,000
  iterations, while with 400,000 iterations every solve of the coupled
  loop converges (the last in 2,784): BiCGSTAB's erratic convergence, not
  an ill-posed system.
- A preconditioner for the basins' systems (the multigrid of M4-09, or an
  incomplete factorisation), or a restart of BiCGSTAB on stagnation, so
  that every basin of the bounds tests converges well within the cap;
  Earth's currents unchanged within the solver's tolerance.

Out:
- The ocean model itself (M4-04).

## Implementation notes

To be refined: which preconditioner suits the non-symmetric Stommel
operator with its held row and the free land masses; the cost on 650,000
cells (Earth's world ocean: 3,532 iterations, 37 s, M4-04).

## Acceptance criteria

- [ ] Corner 105 of the bounds tests no longer warns of its currents and
      leaves `UNCONVERGED_CURRENTS`; every basin of the bounds tests
      converges; Earth's currents within the solver's tolerance of today's.
- [ ] `mise run ci` passes.

## Verification

The bounds tests in release; the iterations of each basin printed by a
scratch tool.

## Result

_(filled in when the ticket is implemented)_

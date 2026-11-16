# M4-15 The periodic state of slowly converging planets

- Status: draft
- Milestone: M4
- Depends on: M4-09
- Crates: kuniumi-sim (stages/climate)
- Design: 気候モデル / 計算の流れ
- Decisions: D-08

## Goal

Planets with little heat transport or close to an ice–albedo bifurcation
reach their periodic state within the 100 years of M4-02, by a solver of
the periodic state itself rather than by running the years one after
another.

## Scope

In:
- The spin-up of planets with little heat transport (decided with the owner
  in M4-02): with a thin atmosphere and a fast spin, `D` falls to `D⊕ / 91`
  and sea ice spreads as a front of about a cell a year, so the periodic
  state took 50 to 270 years on 10,000 cells and up to 988 on 40,000
  (8 h, 0.1 bar, obliquity 60°), and 12 h at 0.3 bar and 60° up to 153
  years. M4-02 stops after 100 years and reports such planets unconverged.
- Planets close to an ice–albedo bifurcation are slow whatever their `D`
  (decided with the owner in M4-02): Earth at 1,296 W/m² (2,000 cells,
  seed 2) is mid-way into a snowball after 100 years (11.5 K change,
  periodic in year 113), Earth at 1,282 W/m² (seed 1, 10,000 cells) takes
  70 years, and 24 h at 0.5 atm with obliquity 90° (seed 2, 10,000 cells)
  64 and 59 years at 1,370 and 1,375 W/m².
- A deterministic solver of the periodic state (Anderson acceleration of
  the year map, Newton–Krylov shooting, or similar) for the energy balance
  (M4-02) and the sea surface with currents (M4-05), verified on these
  cases; split from M4-09 (approved by the owner in M4-09).

Out:
- The cost of each model year (M4-09).
- The physics of the energy balance (M4-02, M4-10, M4-11).

## Implementation notes

To be refined: the method and its source, how a periodic state that the
year map approaches through a transition (the snowball) is told from a
spurious fixed point, and the convergence criterion.

## Acceptance criteria

- [ ] The thin, fast planets of M4-02 (0.1 bar, 8 h, obliquity 45° to 90°;
      0.3 bar, 12 h, 60°) and its near-bifurcation cases (Earth at 1,296
      and 1,282 W/m²; 24 h at 0.5 atm with obliquity 90° at 1,370 and
      1,375 W/m², seed 2, 10,000 cells) reach a periodic state within 100
      years, or the remaining cases are documented with the owner.
- [ ] Earth's periodic state unchanged within the periodic tolerance.
- [ ] `mise run ci` passes.

## Verification

The cases above on 2,000 to 40,000 cells, with and without the solver.

## Result

_(filled in when the ticket is implemented)_

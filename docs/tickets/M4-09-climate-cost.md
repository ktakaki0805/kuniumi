# M4-09 The climate model at 650,000 cells

- Status: draft
- Milestone: M4
- Depends on: M4-02
- Crates: kuniumi-sim (operators, stages/climate)
- Design: 気候モデル / ExoPlaSim との使い分け; 実装ロードマップ / 計算コストの目安
- Decisions: D-08

## Goal

The energy balance of M4-02 (and the climate model's other solves) runs at
the default 650,000 cells in minutes rather than an hour, before M4-07 wires
the Climate stage.

## Scope

In:
- Measured in M4-02 (release, per model year): 0.37 s at 10,000 cells,
  3.2 s at 40,000, 33 s at 160,000, 277 s at 650,000; a periodic state needs
  12 to 20 years, so about 80 minutes at 650,000 cells. The
  Jacobi-preconditioned conjugate gradient needs more iterations the finer
  the grid (78, 153, 333, 662 per solve), because the diffusion `D L` is
  stiff at the grid's scale (10⁵ W m⁻² K⁻¹ against `C/Δt + B` ≈ 18 for
  land).
- Choose, with the owner, between:
  - a coarse climate grid (about 40,000 cells, whose angular resolution
    resolves the diffusion length `√(D/B)` ≈ 0.5 rad), aggregating heat
    capacity and albedo from the planet grid and interpolating the
    temperature back, with the lapse rate on each cell's own height; this
    changes the design's table ("解像度: 全球格子と同じ"); or
  - a multigrid (for example smoothed-aggregation) preconditioner for
    `Laplacian::solve`, deterministic, keeping the full resolution; also
    useful for the Stommel solve of M4-04.
- Timing at 650,000 cells recorded; results compared with the full-resolution
  model on 160,000 cells.
- The spin-up of planets with little heat transport (decided with the owner
  in M4-02): with a thin atmosphere and a fast spin, `D` falls to `D⊕ / 91`
  and sea ice spreads as a front of about a cell a year, so the periodic
  state took 50 to 270 years on 10,000 cells and up to 988 on 40,000
  (8 h, 0.1 bar, obliquity 60°), and 12 h at 0.3 bar and 60° up to 153
  years. M4-02 stops after 100 years and reports such planets unconverged.
  Planets close to an ice–albedo bifurcation are slow whatever their `D`
  (decided with the owner in M4-02): Earth at 1,296 W/m² (2,000 cells,
  seed 2) is mid-way into a snowball after 100 years (11.5 K change,
  periodic in year 113), Earth at 1,282 W/m² (seed 1, 10,000 cells) takes
  70 years, and 24 h at 0.5 atm with obliquity 90° (seed 2, 10,000 cells)
  64 and 59 years at 1,370 and 1,375 W/m². Look for a solver
  of the periodic state itself (Anderson acceleration of the year map,
  Newton–Krylov shooting, or similar), deterministic, that reaches it
  within the cap, and verify it on these near-bifurcation cases too.

Out:
- The physics of the energy balance (M4-02, M4-10, M4-11).

## Implementation notes

To be refined: the choice above (owner), and then its accuracy targets.

## Acceptance criteria

- [ ] A periodic state of Earth's energy balance at 650,000 cells within a
      time agreed with the owner, recorded.
- [ ] Zonal means and seasonal ranges agree with the full-resolution model
      at 160,000 cells within a documented tolerance.
- [ ] The thin, fast planets of M4-02 (0.1 bar, 8 h, obliquity 45° to 90°;
      0.3 bar, 12 h, 60°) and its near-bifurcation cases (Earth at 1,296
      and 1,282 W/m²; 24 h at 0.5 atm with obliquity 90° at 1,370 and
      1,375 W/m², seed 2, 10,000 cells) reach a periodic state within 100
      years, or the remaining cases are documented with the owner.
- [ ] `mise run ci` passes.

## Verification

Timing runs at 160,000 and 650,000 cells; maps compared side by side.

## Result

_(filled in when the ticket is implemented)_

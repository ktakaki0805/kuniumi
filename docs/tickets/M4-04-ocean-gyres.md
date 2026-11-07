# M4-04 Ocean basins and wind-driven gyres

- Status: draft
- Milestone: M4
- Depends on: M4-03
- Crates: kuniumi-sim (stages/climate)
- Design: 気候モデル / 構成要素
- Decisions: none

## Goal

Surface ocean currents from the wind field: subtropical and subpolar gyres with
western intensification, and circumpolar flow where a latitude band is open.

## Scope

In:
- Ocean basin labelling (connected ocean cells).
- Stommel (1948) barotropic model per basin,
  `r ∇²ψ + β ∂ψ/∂x = curl τ / (ρ H)`, with ψ = 0 on coasts, solved on the grid.
- Annual-mean wind stress as forcing; outputs current east and north (m/s).
- Circumpolar channels (no coast boundary around the planet).

## Implementation notes

To be refined: discretisation of `∂/∂x` on the Voronoi grid, solver choice
(reuse conjugate gradient or another iterative solver), parameter values.

## Acceptance criteria

- [ ] An Earth-like basin shows a clockwise subtropical gyre in the northern
      hemisphere with the strongest flow on its western side, reversed in the
      south.
- [ ] Currents are tangent to coasts (no flow into land).
- [ ] Current maps viewed and described.
- [ ] `mise run ci` passes.

## Verification

Current maps rendered to `out/M4-04/` and viewed.

## Result

_(filled in when the ticket is implemented)_

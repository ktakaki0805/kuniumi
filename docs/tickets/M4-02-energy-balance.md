# M4-02 Energy-balance temperature model

- Status: draft
- Milestone: M4
- Depends on: M4-01
- Crates: kuniumi-sim (stages/climate)
- Design: 気候モデル / 構成要素; 気候モデル / 計算の流れ
- Decisions: D-12

## Goal

Monthly surface temperature (K) for every cell from a seasonal 2D energy-balance
model with land and ocean heat capacity, heat diffusion, ice–albedo feedback and
the lapse rate.

## Scope

In:
- `C ∂T/∂t = Q(1 − α) − (A + B·T) + ∇·(D ∇T)` on the grid (finite-volume
  operator of M3-04), integrated over months to a periodic steady state.
- Heat capacity of an ocean mixed layer versus land; albedo from surface type
  and snow or ice (temperature threshold); greenhouse strength through `A` (the
  calibration knob).
- Sea-level temperature, then a lapse rate of 6.5 K/km applied to land height.
- The diffusion coefficient's dependence on rotation rate documented (needed
  for non-Earth planets).

Out:
- Ocean heat transport by currents (added as a source term by M4-05).

## Implementation notes

To be refined: constants A, B and D (North et al. 1981), time step, reuse of the
implicit solver of M3-04, convergence criterion.

## Acceptance criteria

- [ ] Earth: zonal-mean annual temperature within a few K of observations
      (about 300 K at the equator, about 250 K at the poles); global mean about
      288 K after calibrating `A`.
- [ ] Continental interiors have a larger seasonal amplitude than oceans at the
      same latitude.
- [ ] Converges to a periodic state (year-to-year change below a documented
      threshold).
- [ ] `mise run ci` passes.

## Verification

Temperature maps for January and July rendered to `out/M4-02/` and viewed.

## Result

_(filled in when the ticket is implemented)_

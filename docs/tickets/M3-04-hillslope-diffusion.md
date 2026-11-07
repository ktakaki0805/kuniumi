# M3-04 Hillslope diffusion

- Status: ready
- Milestone: M3
- Depends on: M3-03
- Crates: kuniumi-sim (stages/erosion, a reusable `operators` module)
- Design: 侵食と地形発達 / 基本方程式
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

## Acceptance criteria

- [ ] For `h = z` (a degree-1 spherical harmonic) the discrete Laplacian
      approximates `−2·z / R²`, with the error shrinking under refinement
      (errors reported for 10⁴, 10⁵ and 10⁶ cells).
- [ ] Mass conservation on a closed domain (all cells land): Σ h·a stays
      constant within 1e-12 relative per step.
- [ ] Unconditionally stable for large Δt; a single spike decays monotonically.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-sim operators
```

## Result

_(filled in when the ticket is implemented)_

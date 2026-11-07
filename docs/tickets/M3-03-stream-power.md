# M3-03 Implicit stream-power erosion with uplift

- Status: ready
- Milestone: M3
- Depends on: M3-02
- Crates: kuniumi-sim (stages/erosion)
- Design: 侵食と地形発達 / 基本方程式
- Decisions: none

## Goal

River incision balanced against tectonic uplift, solved implicitly so large time
steps stay stable: the core of realistic relief.

## Scope

In:
- `∂h/∂t = U − K Q^m S` (n = 1), with Q the discharge (m³/s) from the drainage
  graph, solved per step in stack order (Braun & Willett 2013):
  `h_i ← (h_i + Δt·U_i + F·h_r) / (1 + F)`, `F = Δt·K·Q_iᵐ / L_ir`.
- Units: with E in m/s and Q in m³/s, K is in m^(1−3m)·s^(m−1). Document the
  default K, m = 0.45, Δt and how they were chosen.
- Ocean cells are a fixed base level (not eroded here).
- `ErosionSettings` in `Settings` (K, m, Δt in s, step count) and a uniform
  precipitation placeholder of 1 m/yr (in m/s) until the climate stage exists.
- A step driver: rebuild drainage, route depressions, solve, repeat; progress
  reported per step.

Out:
- Diffusion (M3-04); deposition and isostatic rebound (M7).

## Acceptance criteria

- [ ] Steady state on a synthetic 1D channel (a chain of cells with known Q and
      uniform U): the slope converges to `S = U / (K Qᵐ)` within 1e-6 relative.
- [ ] Unconditional stability: time steps of 10⁴, 10⁵ and 10⁶ years (in s)
      never create new pits or NaN, and no cell drops below its receiver.
- [ ] Without uplift, total land volume decreases monotonically.
- [ ] Before and after images of `surface_height` on a 100,000-cell planet
      viewed and described (dendritic valleys appear).
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-sim erosion
```

## Result

_(filled in when the ticket is implemented)_

# M2-08 Physical bounds for the planet parameters

- Status: draft
- Milestone: M2
- Depends on: none
- Crates: kuniumi-world
- Design: プロジェクト方針
- Decisions: D-12

## Goal

`PlanetParams::validate` rejects values outside the range of Earth-like
planets that v1 generates, so that no stage meets a radius, gravity or period
for which its arithmetic overflows or underflows.

## Scope

In:
- Lower and upper bounds for `radius_m`, `surface_gravity_m_s2`,
  `rotation_period_s`, `orbital_period_s`, `stellar_flux_w_m2`,
  `ocean_volume_m3` and `surface_pressure_pa` (now any positive finite
  value, or any non-negative one for the ocean volume), chosen with the owner
  for the Earth-like range of D-12; `ParamError::OutOfRange` names the
  parameter as now.

Out:
- Settings files (M2-07), which check the same bounds through `validate`.

## Implementation notes

Found by the review of M2-02 (R2-1): with `radius_m = 5e-324` the Tectonics
stage's `radius_m * lever` underflows to 0, so the angular speed and every
velocity are infinite (caught only by the runner's debug assertion); with
`radius_m = f64::MAX` the angular speed is subnormal and the fastest speed is
off by about 1e-7 relative. Radii from 1e-310 m to 1e300 m were fine. The
bounds themselves (for example 0.1 to 10 Earth radii) are a decision for the
owner when this ticket is refined.

## Acceptance criteria

- [ ] Each bounded parameter: its bounds accepted, values just outside and
      non-finite values rejected by name.
- [ ] The Earth parameters and the bounds' extremes generate finite worlds
      through the Tectonics stage.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-world
```

## Result

_(filled in when the ticket is implemented)_

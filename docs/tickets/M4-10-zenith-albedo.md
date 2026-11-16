# M4-10 Albedo by the solar zenith angle

- Status: draft
- Milestone: M4
- Depends on: M4-02
- Crates: kuniumi-sim (stages/climate)
- Design: 気候モデル / 構成要素
- Decisions: D-12

## Goal

The ice-free albedo of the energy balance rises with the solar zenith angle
of each latitude and month rather than with `P₂(sin φ)`, so that planets of
high obliquity, whose poles see a high sun in summer, get the right
absorbed sunlight.

## Scope

In:
- M4-02 uses `α_free = α_type + α₂ P₂(sin φ)` (`α₂` = 0.08, after North &
  Coakley 1979), which stands for the larger zenith angles and cloud cover
  of Earth's high latitudes. At obliquity 60° to 90° the poles' annual
  insolation exceeds the equator's (M4-01), and this term makes their albedo
  too high.
- The insolation-weighted mean cosine of the zenith angle per latitude and
  month (from the daily-mean geometry of M4-01, `∫ cos² z / ∫ cos z` over
  the day), and an albedo law in it with a cited source, calibrated so that
  Earth keeps the zonal albedo and temperatures of M4-02.

Out:
- Clouds as a separate field.

## Implementation notes

To be refined: the source of the zenith-angle law (Fresnel reflectance of
water, cloud albedo with zenith angle as in Cess 1976 or Williams & Kasting
1997), how it enters the monthly means.

## Acceptance criteria

- [ ] Earth's zonal-mean temperatures and global albedo stay within the
      tolerances of M4-02.
- [ ] At obliquity 90° the poles' annual mean absorbed sunlight follows
      their insolation (documented comparison with the `P₂` form).
- [ ] `mise run ci` passes.

## Verification

Zonal profiles for obliquities 23.44°, 60° and 90° rendered to `out/M4-10/`.

## Result

_(filled in when the ticket is implemented)_

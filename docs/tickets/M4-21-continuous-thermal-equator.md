# M4-21 A continuous thermal equator

- Status: done
- Milestone: M4
- Depends on: M4-20
- Crates: kuniumi-sim (stages/climate: circulation, coupling; the bounds tests of lib.rs), kuniumi-cli (tests)
- Design: 気候モデル / 構成要素; 気候モデル / 計算の流れ
- Decisions: D-08, D-22

## Goal

The thermal equator of the circulation moves continuously with the
temperatures, so that the coupled loop can settle on the planets where
two latitudes are almost equally warm.

## Scope

In:
- The thermal equator of M4-03 is the warmest latitude of the zonal fit
  of a month's sea-level temperatures, capped at ±60°: the argument of a
  maximum, a discontinuous function of the temperatures. Where two
  latitudes are almost equally warm it jumps from one to the other, and
  the coupled loop of M4-17 may have no fixed point (M4-20, "Diagnosis"):
  at 1,000 cells, seed 3, corners 68, 84, 86, 92, 94 and orbit planet 1
  (240-hour days with 10 bar: the summer pole is the warmest, and in a
  month near an equinox the poles take turns, the thermal equator
  jumping by 120°) and corners 73, 75 and 91 (jumps of 1–40° near the
  equator) do not settle even with M4-20's relaxation; at 10,000 cells
  11 of the 13 planets of the bounds tests that do not settle are such
  planets (corners 48, 50, 84, 86, 92, 94, 96, 98, 112, 114, orbit 1), so
  the finer the grid, the more this matters.
- Measured in M4-20 (diagnosis only): freezing each month's thermal
  equator once the loop relaxes settles all but corner 68 and orbit 3,
  but leaves the frozen thermal equators up to 120° from the warmest
  latitude (0.17–3.3 K below it); a thermal equator that follows its own
  local maximum from the year before (by ascent) settles 19 of the 26
  planets of M4-20 instead of 16, but leaves some on a much lesser
  maximum (616 K below the warmest for orbit 0).
- A rule that is continuous where two maxima compete and keeps Earth's
  thermal equator near the warmest latitude where one dominates (below,
  "The plan"), its effect on Earth (world-breaking; M4-05's and M4-06's
  calibrations checked), and the bounds tests' lists measured again.

Out:
- The circulation's belts themselves (M4-03), Earth's calibration (M4-08).
- Earth at e = 0.99 (orbit planet 3), whose monthly winds reach 700 m/s:
  beyond the circulation's range of validity (M4-03, "Limitations"; M4-20).
- The corners whose loop does not settle for reasons other than the
  thermal equator (at seed 1, corners 64 and 66, whose thermal equator
  stays at −1°: the equatorial upwelling of M4-20); reported.
- Two zones of convergence at once (a double ITCZ): the belts have one
  equatorial low (M4-03).

## Implementation notes

The ticket was `draft`. A scratch prototype (not committed): the thermal
equator of `CirculationModel::month` chosen by a rule, and the public
`coupling::coupled_sea_surface` (M4-20's loop with its relaxation) run on
the planets of the bounds tests through the Crust stage (Earth, the 128
corners, the 6 orbit planets; the energy balance without currents as the
start, as `stages::climate::climate` does); release, this machine. With
the rule of M4-03 it reproduces M4-20's numbers exactly (seed 3: 125 of
135 settle, 2,872 coupled years).

### Rules considered

Notation: `f(φ)` the Legendre fit of degree 8 of a month's sea-level
temperatures (M4-03, `ZonalFit`), evaluated on the 1,201 latitudes
`φ_k = −60° + 0.1° k` of M4-03's search; `f* = max_k f(φ_k)`.

- **The energy-flux equator** (T. Bischoff & T. Schneider 2014,
  "Energetic constraints on the position of the intertropical convergence
  zone", J. Climate 27, 4937–4951): the ITCZ near the latitude where the
  atmosphere's meridional energy flux vanishes, to first order
  `δ ≈ −F₀ / (a NEI₀)` (the flux `F₀` across the equator and the net
  energy input `NEI₀` there). The energy balance's transport is
  diffusive (M4-02), so its flux vanishes exactly where `∂f/∂φ = 0`: the
  energy-flux equator is the warmest latitude again (the same
  discontinuity), and its linearisation is a Newton step `−f'(0)/f''(0)`
  from the equator, which diverges where `NEI₀` (`∝ −f''(0)`) changes sign
  and points away from a warm pole (the 240-hour planets near an
  equinox, whose equator is the coldest latitude). Not proposed.
- **A weighted centre of the warm latitudes**: a weighted mean latitude,
  as the ITCZ is located in GCM studies by the centroid of the zonal-mean
  precipitation (the precipitation centroid of D. M. W. Frierson &
  Y.-T. Hwang 2012, "Extratropical influence on ITCZ shifts in slab ocean
  simulations of global warming", J. Climate 25, 720–733; the weighted
  mean of O. Adam, T. Bischoff & T. Schneider 2016, "Seasonal and
  interannual variations of the energy flux equator and ITCZ. Part I",
  J. Climate 29, 3219–3230, the precipitation to the 10th power as the
  weight). The weight here is the warmth: a latitude counts if its fit is
  within a band `ΔT` of `f*`, by its excess `f(φ) − (f* − ΔT)`, and by its
  area (`cos φ`). Measured with `ΔT` = 0.25 to 4 K; also with the
  Boltzmann weight `exp((f − f*)/τ)` (`τ` = 0.25, 0.5, 1 K), with a band
  relative to the month's range of the fit (1 %, 3 %, 10 %) and with the
  centre taken over the whole sphere and then capped at ±60°.

The centre is continuous because every weight is: `f` is linear in the
temperatures, `f*` and the positive part are 1-Lipschitz, so a change of
the fit by at most `ε` (K) changes each weight by at most `2ε`, and the
largest value's own weight, `ΔT`, keeps the denominator at least
`½ cos 60° ΔT` times its quadrature weight. Where two separated warm
regions trade places as the warmer, any continuous rule must cross the
latitudes between them (by the intermediate value theorem), so at a tie
the thermal equator lies between them: a compromise that the loop may
settle on, as below.

### Measured

Settled (the criterion of M4-20: the largest residual and the change of
a year both at most 10⁻³ K within 100 years) of the 135 planets, and the
coupled years summed over them:

| Rule | seed 3 | seed 1 | seed 2 | Lost (settled with M4-20, not with the rule) |
| --- | --- | --- | --- | --- |
| M4-20 (the warmest latitude) | 125 (2,872 years) | 123 (3,103) | 126 (2,266) | |
| centre, `ΔT` = 0.25 / 0.5 K | 131 / 130 | | | c122 / c82 (seed 3) |
| centre, `ΔT` = 1 K | 133 (2,245) | 127 | 133 | c74, c91 (seed 1) |
| centre, `ΔT` = 1.5 K | 132 | 129 | 133 | c80, c123 (seed 3); c74, c91 (seed 1) |
| **centre, `ΔT` = 2 K (the plan)** | **134 (1,928)** | **129 (2,148)** | **133 (1,427)** | **c74, c91 (seed 1)** |
| centre, `ΔT` = 2.5 / 3 / 4 K | 134 / 134 / 134 | 129 / 128 / 129 | 133 | c74, c91 (and c88 at 3 K) (seed 1) |
| centre over the sphere, capped, 2 K | 132 | 129 | 133 | c80, c82 (seed 3); c74, c91 |
| Boltzmann, `τ` = 0.25 / 0.5 / 1 K | 132 / 133 / 134 | – / – / 129 | – / – / 133 | c80 (0.5 K); c74, c91 (1 K, seed 1) |
| band relative to the range, 1 / 3 / 10 % | 125 / 126 / 127 | | | c122 (1 %), c19 (10 %) |

With `ΔT` = 2 K at seed 3 every planet settles but orbit 3 (Earth at
e = 0.99, out of scope): the 10 of M4-20's lists but orbit 3 (corners 68,
73, 75, 84, 86, 91, 92, 94, orbit 1), none lost; 32 of the 135 relax
instead of 38. At seed 1 the 6 left are corners 64 and 66 (thermal
equator at −1° all year: the equatorial upwelling, not the thermal
equator), 74 and 91 (they settled with M4-20 in 98 and 90 years; with the
rule their residual after 100 years is 2.0 and 2.4 × 10⁻³ K, still
falling), orbits 0 (1.4 × 10⁻³ K) and 3. At seed 2 orbits 0 and 3. Every
variant of the table loses corners 74 and 91 at seed 1; narrower bands
settle fewer at seed 3, and the relative band is too narrow on the hot
planets (a range of thousands of kelvin).

At 10,000 cells (seed 3, every planet of the bounds tests, about 20
minutes): 133 of the 135 settle (M4-20: 122), the 2 left orbits 0 and 3;
the coupled years sum to 1,310 instead of 2,524.

Earth, `ΔT` = 2 K. From the periodic temperatures of the energy balance
(the circulation's own input in M4-03's tests), the twelve monthly
thermal equators against the warmest latitude of the same fit: within
0.61°, 0.15° and 0.35° at 10,000 cells (seeds 1, 2, 3), 0.085°, 0.13° and
0.27° at 40,000 (`ΔT` = 1 K: 0.31° at most; 4 K: 1.14°). In the coupled
state, against M4-20's: within 1.22°, 0.23°, 0.35° (10,000 cells) and
0.47°, 0.61°, 0.91° (40,000), the seasonal range 0.3–1.2° narrower (40,000
cells, seed 1: from 8.47° to 7.63°); the loop settles in the same 11–12
years, never relaxed. The calibration of M4-05 and M4-06 on the 40,000-cell
world (`earth_is_calibrated_with_its_coupled_climate`) stays: a mean
surface temperature of 288.011 K (288.016 K with M4-20), `P₂` −27.954 K
(−27.939), the worst 10° band 4.66 K from North & Coakley's fit (4.68),
precipitation 1.0896 m/yr (1.0900).

What the rule does elsewhere (seed 3, the returned state against the
warmest latitude of its own temperatures): where the fit is nearly flat
over a wide range of latitudes (corners 12–15, 28–31, 60–63: within 0.5–
0.6 K over tens of degrees) the thermal equator is the middle of the
plateau, 10–40° from the warmest latitude; where the summer pole is the
warmest the centre lies a few degrees inside the cap (a fit
`300 + 20 sin φ`: 56.45° instead of 60°), and a broad peak near the cap
is pulled equatorward, as the band beyond 60° is cut (`300 − 40 (sin φ −
sin 59°)²`: 51.39°; the centre taken over the whole sphere and then
capped avoids this but lost corners 80 and 82 at seed 3, above). Within
one warm region the thermal equator stays within 2 K of the warmest
temperature (the centre of an interval lies in it); and where two separated warm
regions are within 2 K of each other the loop settles on a thermal
equator between them, below the warmest by 9.3 K (corner 82, 22.5° from
the warmest latitude), 6.1 K (88), 1.6 K (91) and 0.9 K (90) — M4-20
settled 82 and 88 at their warmest latitude, and the frozen thermal
equators of M4-20's diagnosis were up to 3.3 K below.

### The plan (for the owner's approval)

#### The rule (`stages/climate/circulation.rs`)

```text
φ_T = Σ_k c_k cos φ_k e_k φ_k / Σ_k c_k cos φ_k e_k,
e_k = max(0, f(φ_k) − f* + ΔT),   f* = max_k f(φ_k),
```

over M4-03's 1,201 latitudes `φ_k` within ±60°, `c_k` the trapezoidal
weights (½ at the two ends, 1 otherwise), `ΔT` = 2 K: the area-weighted
mean latitude of the latitudes within 2 K of the warmest, each weighted
by how much warmer it is than `f* − ΔT`. It lies within ±60° (a convex
combination of the `φ_k`). `ΔT` is a choice of the model, like the
1,000 km of M4-17's smoothing, not a cited parameterisation: about the
offset of the threshold of deep convection over the tropical oceans
below the warmest sea surface (convection sets in near 27.5 °C, N. E.
Graham & T. P. Barnett 1987, "Sea surface temperature, surface wind
divergence, and convection over tropical oceans", Science 238, 657–659,
against 29–30 °C in the warm pool; the threshold moves with the tropical
mean, N. C. Johnson & S.-P. Xie 2010, "Changes in the sea surface
temperature threshold for tropical convection", Nature Geoscience 3,
842–845), and the band that settled the most planets without losing one
at seed 3 (above).

```rust
/// `ΔT`, the band below the warmest zonal-mean temperature whose
/// latitudes make the thermal equator (module documentation, "The thermal
/// equator"): 2 K.
pub const THERMAL_EQUATOR_BAND_K: f64 = 2.0;

impl ZonalFit {
    /// The warm centre (rad) within `|φ| ≤ cap_rad`: the rule above on
    /// the steps of [`maximum_rad`](Self::maximum_rad).
    ///
    /// # Panics
    ///
    /// If `cap_rad` is not finite and positive, or `band_k` not finite
    /// and positive.
    pub fn warm_centre_rad(&self, cap_rad: f64, band_k: f64) -> f64;
}
```

`CirculationModel::month` takes `fit.warm_centre_rad(THERMAL_EQUATOR_CAP_RAD,
THERMAL_EQUATOR_BAND_K)` instead of `fit.maximum_rad(…)`; `maximum_rad`
stays (the warmest latitude, for the tests and diagnostics). Cost: one
cosine per latitude and month, well under a millisecond a month. The
module documentation's "The thermal equator" gets the rule, its sources,
its continuity, Earth's numbers and what it does on the extreme planets
(above); "Limitations" the compromise between two warm regions. Sums in
index order with `det::sin` and `det::cos`; no randomness.

#### Tests (`circulation/tests.rs`)

New, with reference values from a new script
`tools/thermal_equator_reference.py` (mpmath: the centre of the
continuous profile, `∫ asin(z) e(z) dz / ∫ e(z) dz` over the band's
intervals in `z = sin φ`, the extremum and the band's ends as polynomial
roots; independent of the code's 0.1° sum on a fitted grid):
- `the_thermal_equator_is_continuous_between_two_competing_maxima`
  (criterion 1): on 10,000 cells the zonally symmetric
  `T = 300 − 40 (z² − ¼)² + δ z` (two peaks near ±30°, 2.5 K above the
  valley at the equator, the northern warmer by about `δ`), `δ` from −4 to
  +4 K in steps of 0.02 K: the thermal equator rises monotonically, by at
  most 0.5° a step (the script's centre by at most 0.376° a step, at
  `δ` ≈ 1.4 K), is within 0.05° of the equator at `δ` = 0 and within
  0.05° of the script's centre at `δ` = ±0.5, ±1, ±2, ±2.5, ±4 K (e.g.
  13.41° at 1 K, 30.35° at 2.5 K, where the southern peak has left the
  band); the warmest latitude (`maximum_rad`) jumps by more than 50° in
  one step at `δ` = 0.
- `the_thermal_equator_is_the_centre_of_the_warm_latitudes`: single
  peaks (parabolas in `z` at 12.3°, 20°, −7° and 59°, the warmest at the
  cap `300 + 20 z`; the script: 12.367°, 20.162°, 51.389°, 56.452°, …)
  within 0.05° of the script's centre and within 2 K of the warmest
  value of their fit; a uniform field gives the area-weighted centre of
  the cap, 0 within 10⁻¹²; `|φ_T| ≤ 60°`.
- `earth_s_thermal_equator_is_near_its_warmest_latitude` (criterion 1,
  Earth): the shared Earth run (40,000 cells in release, 10,000 in debug
  builds), every month within 1° of `maximum_rad` of the same fit
  (measured 0.085° and 0.61°).
- `warm_centre_needs_a_positive_band` (`#[should_panic]`).

Changed: `belt_winds_match_the_closed_form` (its thermal equators from
the script: 0°, 20.162° and 56.452° instead of 0°, 20° and 60° exactly);
`a_zero_length_is_the_circulation_of_m4_03` renamed
`a_zero_length_skips_the_smoothing` (the circulation of M4-03 with the
new thermal equator; its fingerprint recorded again);
`the_itcz_follows_the_seasons` unchanged (prints the new values).

#### Golden values (world-breaking)

Recorded again, after the invariant tests pass and the maps are viewed:
the circulation's (`golden_fingerprint_of_earth_circulation`, the
unsmoothed one), the ocean's, the sea surface's, the moisture's, the
relaxed loop's (`golden_fingerprint_of_a_relaxed_loop`, corner 72), the
climate stage's and the CLI world's climate fields. Unchanged: the
insolation's and the energy balance's, the erosion's and the drainage's
(erosion runs on the uniform precipitation, M4-19), the CLI's fields
before the climate. `an_oscillating_loop_settles_relaxed` (corner 72:
measured 35 years, relaxed from year 5, instead of 34) re-measured; if
corner 72 no longer oscillates without relaxation, another corner that
does (from the prototype's list) takes its place.

#### The bounds tests (`lib.rs`)

`UNSETTLED_SEA_SURFACE` and `UNSETTLED_ORBITS` measured again with the
implementation by the flags tool of M4-20 (each model's own flags
through the public `stages::climate::climate`; expected from the
prototype `[]` and `[3]`), and `UNSETTLED_MOISTURE` likewise, with their
documentation; the comment on the debug builds' corners.

#### Documents

- `coupling.rs`: "Relaxation" (the paragraph on the thermal equator that
  jumps), "Convergence" and "Cost" (the numbers measured again: 1,000
  and 10,000 cells, seeds 1–3); `climate.rs`: the count of unsettled
  planets; Earth at 650,000 cells (seed 1) generated once to check its
  years and cost.
- `docs/design-ja.md`: 気候モデル / 構成要素, the row 気圧帯
  (熱赤道 = 最も暖かい緯度から 2 K 以内の緯度の、暖かさと面積で重み付けした
  平均緯度), and 計算の流れ (the planets whose thermal equator jumped now
  settle; of the bounds tests only Earth at e = 0.99 does not, M4-21).
- `CHANGELOG.md`: Changed (the thermal equator) and World-breaking (the
  climate of every world, Earth's thermal equator by up to about 1° and
  its calibration kept); `tools/README.md`: the new script.

#### World-breaking

Every world's climate: the thermal equator of every month changes (on
Earth by up to 0.6° from the energy balance's temperatures and 1.2° in
the coupled state at 10,000 cells), and with it the pressure, the winds,
the currents, the sea surface, the moisture and the precipitation. The
golden values above change; the calibration stays within its targets.

#### Size

About 600 changed lines without the ticket: the circulation 100 (with its
documentation), its tests 200, the script 120, the bounds tests and the
coupling's tests 60, golden values 60, documents 60.

### Questions for the owner

1. The rule and its band: the centre of the latitudes within `ΔT` = 2 K
   of the warmest (world-breaking for every world, Earth's thermal
   equator within 0.61° of its warmest latitude, the calibration kept),
   which settles every planet of the bounds tests at seed 3 but Earth at
   e = 0.99 and loses none there (at seed 1, 2 marginal ones). The band
   is a choice of the model, not given by the design. Alternatives:
   `ΔT` = 1 K (Earth within 0.31°, corner 75 still unsettled at seed 3,
   127 at seed 1); the Boltzmann weight with `τ` = 1 K (as many settle,
   Earth within 0.73°, every latitude weighted); no change (the 10
   planets stay unsettled and reported).

## Acceptance criteria

- [x] The thermal equator is a continuous function of the temperatures:
      `the_thermal_equator_is_continuous_between_two_competing_maxima`
      passes; Earth's within 1° of its warmest latitude (M4-03's thermal
      equator): `earth_s_thermal_equator_is_near_its_warmest_latitude`.
- [x] More of the planets of the bounds tests settle than with M4-20
      (1,000 cells, seed 3: at least 130 of the 135, M4-20 125, none lost;
      the lists of the bounds tests measured again); at seeds 1 and 2 and
      at 10,000 cells (seed 3) more settle than with M4-20 (scratch
      tool).
- [x] Earth's calibration kept: `earth_is_calibrated_with_its_coupled_climate`
      passes (288.0 ± 0.3 K, `P₂` −28.0 ± 0.5 K, 1.0 ± 0.1 m/yr), Earth's
      coupled loop settles unrelaxed (`earth_s_loop_is_never_relaxed`).
- [x] The golden values that change are recorded again after the other
      tests pass and the maps of `out/M4-21/` are viewed (thermal equator,
      pressure and winds of Earth and of a corner that newly settles);
      the others unchanged.
- [x] `mise run ci` passes.

## Verification

```sh
cargo test --release -p kuniumi-sim circulation
cargo test --release -p kuniumi-sim coupling
cargo test --release -p kuniumi-sim bounds
tools/.venv/bin/python tools/thermal_equator_reference.py
```

The prototype and the flags tool of M4-20 on the planets of the bounds
tests at seeds 3, 1 and 2 and at 10,000 cells.

## Result

Approved by the owner as planned, with `ΔT` = 2 K (question 1).

What was done:
- `stages/climate/circulation.rs`: `ZonalFit::warm_centre_rad` (the rule of
  the plan, summed from the south over the steps of `maximum_rad`, with
  `det::sin` and `det::cos`) and `THERMAL_EQUATOR_BAND_K` (2 K);
  `CirculationModel::month` takes the warm centre as the thermal equator;
  `maximum_rad` stays, bit for bit, the two sharing a private `search`.
  The module documentation's "The thermal equator" (the rule, its sources,
  its continuity, the choice of `ΔT`, what it gives on Earth and on the
  extreme planets) and "Limitations" (one equatorial low between two warm
  regions).
- Tests (`circulation/tests.rs`), with reference values from
  `tools/thermal_equator_reference.py` (committed with the ticket's work in
  progress): `the_thermal_equator_is_continuous_between_two_competing_maxima`,
  `the_thermal_equator_is_the_centre_of_the_warm_latitudes`,
  `earth_s_thermal_equator_is_near_its_warmest_latitude`,
  `warm_centre_needs_a_positive_band`,
  `warm_centre_needs_a_cap_within_the_poles`; seen to fail against a stub
  that returned the warmest latitude (all but the Earth test, which the stub
  passes trivially). `belt_winds_match_the_closed_form` takes the script's
  thermal equators (0°, 20.162°, 56.452°); `a_zero_length_is_the_circulation_of_m4_03`
  renamed `a_zero_length_skips_the_smoothing`.
- `coupling/tests.rs`: `an_oscillating_loop_settles_relaxed` measured again
  (corner 72 still oscillates without relaxation; relaxed from year 5, 35
  years instead of 34; the exact last change, the residual; the residuals
  of years 3 and 4 in its comment).
- The bounds tests (`lib.rs`): `UNSETTLED_SEA_SURFACE` [], `UNSETTLED_ORBITS`
  [3], `UNSETTLED_MOISTURE` 16 corners (12 before), their documentation and
  that of the debug builds' corners.
- Golden values recorded again (below). Documents: `coupling.rs`
  ("Relaxation", "Convergence", "Cost"), `climate.rs`, `docs/design-ja.md`
  (構成要素, the row 気圧帯; 計算の流れ), `CHANGELOG.md` (Added, Changed,
  World-breaking), `tools/README.md`; follow-up M4-22 (draft) in
  `docs/tickets/README.md`.

Deviations from the plan:
- `warm_centre_rad` also panics for a cap beyond π/2 (where `cos φ` would
  weigh latitudes negatively, so the centre would not be a convex
  combination), with the test `warm_centre_needs_a_cap_within_the_poles`.
- `UNSETTLED_MOISTURE` grows: the prototype measured only the coupled loop,
  and the moisture of four dry corners at the largest radius, the longest
  orbit, the 240-hour day and 10 bar (77, 79, 93, 95) no longer reaches its
  periodic state in 10 years (77, 79: 3.6 and 3.5 × 10⁻⁶ kg/m², just above
  the tolerance of 10⁻⁶; 93, 95: 61 and 57 kg/m², as 85 and 87 since M4-16).
  A scratch test ran corner 93's moisture on for 30 years: a monotone drain
  of the initial vapour from a dry polar cell (−82°, from 609 kg/m² in
  year 7 by about 14 % a year), not an oscillation. At seeds 1 and 2 the
  list changes too (seed 1: 5 and 7 settle, 93 and 95 do not; seed 2: 5, 7,
  77, 79 do not). Recorded in the documentation and as follow-up M4-22,
  outside this ticket's scope (the moisture model).
- The plan described the thermal equator on nearly flat fits as "10–40°
  from the warmest latitude" (corners 12–15, 28–31, 60–63); measured over
  every planet it is up to 60° from it on the 33 planets with a 240-hour
  day and 10 bar, whose fit is within 2 K over most of the cap (within
  0.79 K of the warmest temperature): documented so.

Evidence:

1. Criterion 1.
   `the_thermal_equator_is_continuous_between_two_competing_maxima`
   (10,000 cells): over the 401 values of `δ` the centre rises in every
   step, by at most 0.3756° (script 0.3755°); at `δ` = 0 it is 0.0012°;
   at ±0.5, ±1, ±2, ±2.5 and ±4 K within 0.0045° of the script (13.4153°
   against 13.4109° at 1 K, 30.3554° against 30.3544° at 2.5 K); the
   warmest latitude jumps by 60.02° between `δ` = −0.02 and 0 K.
   `the_thermal_equator_is_the_centre_of_the_warm_latitudes`: the six
   profiles within 0.0097° of the script (51.3992° against 51.3895° for the
   peak at 59°, 56.4541° against 56.4522° for `300 + 20 z`; since review
   round 1 within 0.015° and 0.01° for these two, so that dropping the
   trapezoidal half weights of the ends fails, checked by that mutation:
   51.4229°), at most
   0.652 K below the warmest of their fit; a uniform field 1.9 × 10⁻¹⁷ rad.
   `earth_s_thermal_equator_is_near_its_warmest_latitude`: within 0.085°
   (release, 40,000 cells) and 0.609° (debug, 10,000 cells), as the
   prototype; `belt_winds_match_the_closed_form`: 0.000°, 20.163°, 56.454°.
2. Criterion 2 (scratch tool `<scratchpad>/m4_21/flags`, after M4-20's
   `m4_20_flags`: each model's own flags through the public
   `stages::climate::climate`; its coupled-only mode through
   `coupled_sea_surface`). 1,000 cells, seed 3: 134 of the 135 planets
   settle (M4-20: 125), none lost; corners 68, 84, 86, 92, 94 and orbit 1
   in 4 to 8 years unrelaxed, 73, 75 and 91 relaxed in 38 to 55; only
   Earth at e = 0.99 does not (201 K; its monthly winds 681 m/s); 32 relax
   (38); 1,928 years summed (2,872). Seeds 1 and 2: 129 and 133 (123 and
   126); lost at seed 1 corners 74 and 91 (2.0 and 2.4 × 10⁻³ K after 100
   years), as the prototype. 10,000 cells, seed 3: 133 (122), none lost,
   orbits 0 and 3 left; 1,310 years summed (2,524), 4.4 minutes on 16
   threads. The bounds tests pass with the new lists (release: Earth, the
   128 corners, the 6 orbit planets).
3. Criterion 3. `earth_is_calibrated_with_its_coupled_climate` (release,
   40,000 cells): 288.011 K, `P₂` −27.954 K, the worst band 4.66 K,
   1.0896 m/yr, 24.90 kg/m² (M4-20: 288.016 K, −27.939 K, 4.68 K,
   1.0900 m/yr). `earth_s_loop_is_never_relaxed`: 11 and 12 years, not
   relaxed. The scratch tool's Earth check: not relaxed on 1,000 cells
   (seeds 1, 2, 3: 11, 11, 11 years), 1,200 (seed 5: 10), 1,500 (seed 3:
   10), 10,000 (11, 12, 11) and 40,000 (12, 12, 12), as with M4-20; and on
   650,000 cells (seed 1, generated once through the Crust stage, its
   energy balance and its coupled state, release): 11 years, not relaxed,
   the loop 939 s (M4-17: 926 s), the energy balance 190 s (184 s), its
   thermal equators −7.56° to +3.91°, within 0.20° of the warmest latitude
   on the energy balance's temperatures and 0.85° in the coupled state;
   2.8 GB at most.
4. Criterion 4. Maps in `out/M4-21/` (scratch tool `<scratchpad>/m4_21/maps`,
   built against the branch and against `main`: the coupled state's
   sea-level pressure less the surface pressure on RdBu, wind arrows every
   6°, coasts, latitude lines every 30°, the thermal equator in green and
   the warmest latitude of the fit of the returned temperatures dashed in
   magenta):
   - `earth/branch_m01.png`, `branch_m07.png`, `main_m01.png`,
     `main_m07.png` (40,000 cells, seed 1, ±2,500 Pa, 2.4 px per m/s):
     three belts a hemisphere, the trades converging on the thermal
     equator (−3.15° and +2.35° on the branch, −3.46° and +2.74° on
     `main`), westerlies at 35–55°, the winter highs over the northern
     continents in January and the southern ones in July; branch and
     `main` alike to the eye, the magenta line one or two pixels from the
     green on the branch. The branch's monthly winds differ from `main`'s
     by at most 0.41 m/s (RMS 0.09 m/s), its pressure by at most 24 Pa;
     the seasonal range of the thermal equator 7.62° instead of 8.47°.
   - `c84/main_m07.png`, `main_m10.png`, `main_m11.png` and `branch_…`
     (corner 84, 1,000 cells, seed 3: 1,500 km, 240 h, 10 bar, no water;
     ±6,000 Pa, 6 px per m/s): one cell to each pole. On `main` (100
     years, not settled, last change 6.6 K) the equatorial low and the
     convergence of the meridional winds sit at +28.6° in October and at
     −29.0° in November, far from the warmest latitude (+30.2°, then
     +60°); on the branch (5 years, settled, not relaxed) at +1.9° and
     +0.8°, the flow from both poles meeting there, symmetric, with the
     warmest latitude at +60°, whose fit is within 0.3 K of the thermal
     equator's.
   Golden values recorded again after these and the tests above passed:
   the circulation (`0x8F79_4961_AB46_90A9`; unsmoothed
   `0x69B6_A1CD_AC60_773A`), the currents (`0xAE96_5844_0C53_DC15`), the
   sea surface (`0xAB61_A6F9_838F_1609`), the moisture
   (`0xC171_0A87_B57D_2B3D`), the relaxed loop (`0xB1B8_3871_864B_8328`),
   the Climate stage's 30 fields and the CLI's 30 climate fields.
   Unchanged and passing: the insolation, the energy balance, the erosion,
   the drainage, the CLI's other 16 fields and the sea level.
5. Criterion 5. `mise run ci` passes: fmt, clippy `-D warnings`, rustdoc,
   the browser build's lints, 2,077 tests summed over debug and release,
   MSRV 1.85.0, cargo-deny (advisories, bans, licences, sources ok); 8 min
   35 s. The debug tests of `kuniumi-sim` take 229 s (275 s in M4-20; the
   debug corners 96 and 112 settle in 6 and 8 years instead of 17 and 25),
   the release ones 98 s.

Review:
- Round 1 (APPROVE; the reviewer re-ran `mise run ci`, built its own flags
  tool against the branch and an extracted `main` at seeds 1–3 and 10,000
  cells (main reproducing M4-20's records and the branch these, byte for
  byte), re-ran the reference script and checked four of its values with
  an independent sum, regenerated and viewed the maps, and modelled the
  rule in Python to test mutations). R1-1 minor: the comment of
  `earth_is_calibrated_with_its_coupled_climate` still gave M4-20's
  measurements: the new ones written in (debug 289.03 K, −24.51 K,
  0.963 m/yr, measured again), M4-20's kept in parentheses. R1-2 nit: long
  comment lines and a broken reflow in `climate.rs`: rewrapped. R1-3 nit:
  the `Crates:` line now names `kuniumi-cli (tests)`. R1-4 nit: the comment
  on corner 72's change of its last year reworded. R1-5 nit: the coupled
  state's 1.2° is that of 10,000 cells; on 40,000 within 0.8° (seeds 1–3:
  0.44°, 0.55°, 0.79°), now said so in `circulation.rs`. R1-6 nit: only
  the golden values pinned the trapezoidal end weights: the two profiles at
  the cap now have tolerances of 0.015° and 0.01° (above, evidence 1).

Follow-ups: M4-22 (draft): the periodic state of the moisture of the dry
planets.

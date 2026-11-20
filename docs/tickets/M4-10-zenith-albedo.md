# M4-10 Albedo by the solar zenith angle

- Status: done
- Milestone: M4
- Depends on: M4-02
- Crates: kuniumi-sim (stages/climate: insolation, energy_balance, sea_surface, coupling and their tests; the bounds tests of lib.rs), kuniumi-cli (tests: golden fingerprints), tools (zenith_albedo_reference.py)
- Design: 気候モデル / 構成要素; 気候モデル / 較正：地球を正解として使う
- Decisions: D-08, D-11, D-12, D-22

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
- How the monthly albedo enters the sub-steps, the warm starts and the
  balance step of the energy balance and of the sea surface (M4-05).
- `A` and `D⊕` of the energy balance and of the sea surface fitted again on
  the world of M4-02's calibration (as M4-05 did), so that Earth's
  temperature targets hold (below, "Earth's temperatures").
- World-breaking: every climate changes; the golden fingerprints of the
  climate models, of the Climate stage and of the CLI are recorded again.

Out:
- Clouds as a separate field.
- A zenith-angle dependence of the snow and ice albedo (`α_ice` = 0.6
  stays constant) and of land types.
- Earth's calibration on its real topography (M4-08).

## Implementation notes

The ticket was `draft`. A scratch prototype (not committed: the plan below
on the branch, plus a probe that runs the Earth-parameter world of seed 1 on
40,000 cells through the Crust stage, release, this machine) measured the
numbers quoted here; the same probe on the code of `main` gave the numbers
of the `P₂` form.

### The zenith angle (insolation, M4-01)

At latitude `φ` and declination `δ`, with `A = sin φ sin δ`,
`B = cos φ cos δ` and the hour angle of sunset `h₀` (M4-01),
`cos z = A + B cos h`. Over the day the insolation is `∝ F = ∫₀^h₀ cos z dh
= h₀ A + B sin h₀` (Berger 1978, M4-01), and

```text
G = ∫₀^h₀ cos² z dh = A² h₀ + 2AB sin h₀ + (B²/2)(h₀ + sin h₀ cos h₀)
μ̄_day = G / F
```

is the insolation-weighted mean cosine of the zenith angle of the day
(T. W. Cronin 2014, "On the choice of average solar zenith angle", J.
Atmos. Sci. 71, 2994–3003: the insolation-weighted zenith angle reproduces
the absorbed sunlight of the diurnal cycle better than the daytime mean
`(1/h₀) ∫ cos z dh`, which Vladilo et al. 2013 use). Polar day `π A² +
π B²/2`, polar night 0. Hand values: `π/4` at the equator on an equinox,
`sin δ` at a pole in polar day.

The monthly value weights the days by their insolation: with the
sub-steps of M4-01 (the star's longitude in steps of at most 0.5°, where
Kepler's second law cancels `(a/r)²` from the time integral) it is
`μ̄_k = Σ G / Σ F` over the month's sub-steps, `+0.0` for a month of
polar night (no light, so the albedo does not matter), in `[0, 1]`. The
insolation's own values and golden fingerprint do not change: `F` keeps
its function, `G` shares its `h₀`.

### The albedo law

The ice-free albedo of a cell in month `k`:

```text
α_free,k = α_type + b (μ_r − μ̄_k)        (α_type: ocean 0.263, land 0.363, as in M4-02)
```

Linear in the cosine because the two-stream (Eddington) reflectance of a
non-absorbing cloud layer of optical depth `τ` and asymmetry `g`,

```text
R(μ₀) = [(1 − g) τ + (2/3 − μ₀)(1 − e^{−τ/μ₀})] / [4/3 + (1 − g) τ]
```

(E. P. Shettle & J. A. Weinman 1970, J. Atmos. Sci. 27, 1048–1055; the
conservative case of W. E. Meador & W. R. Weaver 1980, J. Atmos. Sci. 37,
630–643; derived again by hand for this ticket), is linear in `μ₀` for a
thick layer (`e^{−τ/μ₀}` → 0) and rises as the sun sinks; R. D. Cess
(1976, J. Atmos. Sci. 33, 1831–1843) found Earth's zonal cloud albedo to
rise with the zenith angle (used by Williams & Kasting 1997 and Vladilo et
al. 2013). A law linear in `μ` has a second property: the monthly mean of
the instantaneous absorbed sunlight `Q (1 − α(μ))` is exactly
`Q̄_k (1 − α(μ̄_k))` with the insolation-weighted `μ̄_k`, so the monthly
means lose nothing.

The alternative, linear in the zenith angle itself as Cess's clouds
(`α = α_type + β (Z̄ − Z_r)` with `Z̄ = acos μ̄`; Williams & Kasting 1997;
Vladilo et al. 2013, eq. A15, `−0.07 + 0.008 Z` with a floor of 0.19), is
not exact under the monthly means and falls steeply where the sun stands
high (an ocean at 0.035 under an overhead sun, as at the summer pole of a
planet of obliquity 90°), where the albedo of a cosine law is flat
(`dμ/dZ = −sin Z`); it would need a floor, another parameter. The plan
takes the cosine law (a question to the owner while refining).

**Calibration of the law on Earth's albedo.** `b` and `μ_r` make the
annual ice-free absorbed sunlight of Earth's orbit (obliquity 23.44°; the
annual integrals do not depend on the eccentricity, the perihelion or the
calendar, as `Σ_k Q_k μ̄_k ∝ ∫ G dλ` and `Σ_k Q_k ∝ ∫ F dλ`) have the same
global mean and the same `P₂` moment as the `P₂` form:

```text
∫ w(x) [b (μ_r I_F(x) − I_G(x)) − α₂ P₂(x) I_F(x)] dx = 0,   w = 1 and P₂(x),  x = sin φ
I_F(x) = ∫₀^2π F dλ,  I_G(x) = ∫₀^2π G dλ
```

The global and `P₂` moments are what the energy balance's mean and `P₂`
temperatures answer to (NCC81's two-mode solution), so the calibration of
`A` and `D` keeps its meaning. A new script, `tools/zenith_albedo_reference.py`
(mpmath, the day integrals checked against the hour-angle quadrature):
`b` = 0.33802, `μ_r` = 0.64410 (`Z_r` = 49.9°); the code takes `b` = 0.338,
`μ_r` = 0.644. The insolation-weighted annual ice-free albedo of each 10°
band against the `P₂` form (both relative to `α_type`; float sketch, the
script will give the test's values):

| Band | 0–10° | 10–20° | 20–30° | 30–40° | 40–50° | 50–60° | 60–70° | 70–80° | 80–90° |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| zenith law | −0.0362 | −0.0306 | −0.0196 | −0.0041 | +0.0150 | +0.0361 | +0.0585 | +0.0858 | +0.1062 |
| `P₂` form | −0.0388 | −0.0318 | −0.0186 | −0.0008 | +0.0194 | +0.0397 | +0.0576 | +0.0710 | +0.0782 |

Within 0.005 equatorward of 60°; the polar caps, whose sun never rises
above 23.44° plus their colatitude, are brighter by up to 0.03 (they are
under snow and ice most of the year). The law ranges from `α_type − 0.12`
under an overhead sun (`μ̄` = 1; ocean 0.143) to `α_type + 0.218` at a
grazing sun (ocean 0.481, land 0.581), below `α_ice`.

### How it enters the model

- `EnergyBalance::new` takes the monthly zenith cosines of the cells
  (`insolation::monthly_zenith_cosine(params, grid)`) and keeps the
  ice-free albedo of each cell and month (`free_albedos(month)`); the
  albedo of month `k` is `(1 − f) α_free,k + f α_ice` with the ice
  fraction of the surface temperature as before. The four sub-steps of a
  month use its albedo (`substep`, `month` and `albedos` take the month;
  the sea surface's `substep` already does).
- The warm starts (energy balance and sea surface): the annual absorbed
  sunlight is the mean over the months, `S(T) = mean_k Q_k (1 − α_free,k)
  − f(T) mean_k Q_k (α_ice − α_free,k)`, the same expression as before
  when `α_free` does not depend on the month.
- The balance step: `S_i(T + δ)` with each month's albedo.
- `EnergyBalanceParams`: `albedo_p2` is replaced by `albedo_zenith_slope`
  (`b`) and `reference_zenith_cosine` (`μ_r`); `free_albedo(ocean, μ̄)`
  and `albedo(ocean, μ̄, T_s)` take the cosine instead of `sin φ`.

### Earth's temperatures: `A` and `D⊕` fitted again

With the law and the `A`, `D⊕` of M4-02 and M4-05 (40,000 cells):

| | mean surface (K) | `P₂` (K) | worst band (K) | planetary albedo |
| --- | --- | --- | --- | --- |
| energy balance, `P₂` form (main) | 287.959 | −28.030 | 3.80 | 0.3006 |
| energy balance, zenith law | 288.220 | −27.675 | 4.21 | 0.2990 |
| coupled, `P₂` form (main) | 288.011 | −27.954 | 4.66 | 0.3006 |
| coupled, zenith law | 288.284 | −27.562 | **5.10** | 0.2989 |

The coupled state's worst band leaves M4-02's 5 K. Under snow and ice the
brighter caps of the new law do not count, while the mid-latitudes are
0.003–0.004 darker, so the planet warms by 0.26 K and its `P₂` weakens by
0.4 K. So `A` and `D⊕` are fitted again to the targets of M4-02 and M4-05
(a mean surface temperature of 288.0 K and a `P₂` coefficient of −28.0 K)
by Newton steps with differences, rounded to 0.1 W/m² and
0.001 W m⁻² K⁻¹ as before (prototype):

| | `A` (W/m²) | `D⊕` (W m⁻² K⁻¹) | mean | `P₂` | worst band |
| --- | --- | --- | --- | --- | --- |
| energy balance alone | 202.7 → **202.8** | 0.583 → **0.579** | 288.043 | −27.966 | 3.86 |
| sea surface, coupled | 202.6 → **202.8** | 0.508 → **0.502** | 287.991 | −27.973 | 4.44 |

With the refit the coupled climate keeps M4-06's precipitation targets
(1.0 ± 0.1 m/yr, 22–28 kg/m²): 1.096 m/yr and 24.9 kg/m² (1.090 m/yr and
24.9 kg/m² before); the sea surface with the winds of the energy balance
alone (M4-05's test, ±1 K) 287.72 K, −28.75 K, worst band 3.45 K.

### The other tests and the bounds tests (prototype, with the refit)

- Tests whose planets sit at an ice–albedo bifurcation of the `P₂` form
  still pass: Earth at 1,296 W/m² (seed 2, 2,000 cells) freezes over in
  20 years (236.3 K); Earth at 1,282 W/m² settles in 9 years, 24 h at
  0.5 atm and 90° at 1,370 and 1,375 W/m² in 10 and 25 (the test's
  limit: 25). Without the refit the first did not freeze (275 K) and
  Earth at 1,282 W/m² took 27 years.
- `an_oscillating_loop_settles_relaxed` (corner 72): the same 35 years
  and first relaxed year 5; its recorded last change (bits) is recorded
  again.
- The bounds tests (1,000 cells, seed 3) warn for three planets that did
  not warn before:
  - corner 16 (1,500 km, 8 h, 90-day year, 2,500 W/m², no ocean,
    0.1 bar): the energy balance without currents is not periodic after
    100 years (change 1.38 K), the slow ice front of a planet with little
    heat transport that M4-15 accepted as a limit;
  - corner 105 (13,000 km, 8 h, 3,650-day year, 400 W/m², all ocean,
    10 bar): the currents of its one basin miss the solver's tolerance
    (relative residual 2.1 × 10⁻³);
  - orbit planet 0 (obliquity 90°, e = 0.99, 2,500 W/m², 8 h, 0.1 bar,
    monthly means beyond 10,000 K, outside the range of validity of
    `A + BT` and of the circulation): the coupled loop does not settle
    (change 96.6 K after 100 years; it settled relaxed since M4-20).

  The plan: diagnose each within this ticket; where the cause is the
  known limit (an ice front, a climate beyond the models' range), record
  the planet in the lists of `lib.rs` with its diagnosis (as M4-15, M4-17
  and M4-20 did) and open a follow-up `draft` ticket for its
  convergence.

### High obliquity

On a circular orbit at obliquity 90°, by hand: the pole's annual
insolation-weighted `μ̄` is `∫ sin² λ / ∫ sin λ = π/4`, the equator's
`(π/4)²` = π²/16; so the ice-free ocean's annual co-albedo is 0.7848 at
the pole and 0.7278 at the equator, where the `P₂` form gives 0.657 and
0.777. The pole then absorbs 1.694 times the equator's sunlight, against
an insolation ratio of π/2 = 1.571 (`P₂` form 1.328).

Prototype, the coupled climate (Earth's `A` and `D`, 40,000 cells, seed 1):

| Obliquity | | `P₂` form | zenith law |
| --- | --- | --- | --- |
| 90° | 80–90° S: insolation, absorbed, surface | 431.6 W/m², 284.4 W/m², 300.7 K | 431.6, 338.3 W/m², 311.4 K |
| 90° | 5° S–5° N | 276.7, 201.9 W/m², 281.9 K | 276.7, 188.2 W/m², 280.0 K |
| 90° | global mean surface, planetary albedo | 288.22 K, 0.2990 | 290.35 K, 0.2859 |
| 60° | 80–90° S | 374.3, 246.6 W/m², 292.5 K | 374.3, 280.4 W/m², 299.8 K |
| 60° | 5° S–5° N | 334.0, 246.0 W/m², 289.0 K | 334.0, 231.9 W/m², 287.3 K |

### Tests

Reference values from `tools/zenith_albedo_reference.py` (mpmath; the
daily cosine by Gauss–Legendre quadrature over the hour angle, the monthly
one by time integrals with Kepler's equation, as
`insolation_reference.py`) and hand values.

- `insolation`: `the_daily_zenith_cosine_is_the_insolation_weighted_mean`
  (hand values exact: `π/4`, `sin δ`, 0 in polar night; the script's
  values within 10⁻¹²); `monthly_zenith_cosines_match_the_time_integral`
  (Earth and an eccentric orbit, within 2 × 10⁻⁵, the sub-steps' midpoint
  rule); `zenith_cosines_lie_in_the_unit_interval` (the corners of the v1
  orbits: e = 0.99, obliquity 0° and 90°, every latitude, 0 exactly in
  polar night); the insolation's golden fingerprint unchanged.
- `energy_balance`: `the_zenith_law_keeps_earth_s_albedo` (criterion 1:
  on 40,000 cells the global and `P₂` moments of the annual ice-free
  absorbed sunlight within 0.2 % of the `P₂` form's, which the test
  evaluates in closed form; every 10° band within 0.005 equatorward of 60°
  and within 0.03 poleward); `the_poles_absorb_with_their_insolation_at_high_obliquity`
  (criterion 2: obliquity 90° on 10,000 cells, the ice-free co-albedo of
  the cells poleward of 80° at least the equator's, and the hand values
  within 0.005); `the_free_albedo_follows_the_zenith_cosine` (exact values
  at `μ_r`, 0 and 1); `earth_keeps_its_planetary_albedo` (within 0.005 of
  M4-02's 0.3006); the existing tests with the month (the uniform and
  two-mode tests keep `b` = 0); the Earth criteria of M4-02 with the refit.
- `sea_surface`, `coupling`: the calibration tests with the refit (M4-05's
  and M4-17's tolerances, M4-06's precipitation); a test whose planet was
  chosen at the snowball edge of the `P₂` form
  (`a_planet_near_the_snowball_edge_freezes_over`, 1,296 W/m²) is moved to
  the edge of the new law (measured, documented).
- The bounds tests (lib.rs): the lists of unsettled planets measured
  again (above).
- Golden fingerprints of the energy balance, circulation, currents, sea
  surface, coupling, moisture, the Climate stage and the CLI recorded
  again after the other tests pass and the images are viewed.

### Images

`out/M4-10/`, by a scratch tool from the probe's tables: for obliquities
23.44°, 60° and 90° (40,000 cells, coupled), the zonal annual surface
temperature, absorbed sunlight and ice-free albedo of the `P₂` form (main)
and of the zenith law; and January and July maps of the ice-free albedo.

### Size

About 1,000 changed lines (insolation and its tests 200, the energy
balance, sea surface and coupling 250, their tests 300, the script 150,
documents and goldens the rest), over the guideline of 800; one ticket,
because the parts depend on each other (the law needs the cosines, the
calibration needs the law).

### Design document

`気候モデル / 構成要素`, row 気温: the albedo rises with the
insolation-weighted solar zenith angle of each month. `気候モデル / 較正`:
`A` and `D` fitted again with the zenith law (M4-10).

## Acceptance criteria

- [ ] Earth (the Earth-parameter world of seed 1, 40,000 cells, release):
      with the refitted `A` and `D⊕`, the energy balance alone and the
      coupled climate each have every 10° band within 5 K of North &
      Coakley's fit, a `P₂` coefficient of −28.0 ± 0.5 K and a mean surface
      temperature of 288.0 ± 0.3 K; the coupled climate keeps M4-06's
      1.0 ± 0.1 m/yr and 22–28 kg/m²; the planetary albedo is within 0.005
      of M4-02's 0.3006; Earth's ice-free albedo keeps the global and `P₂`
      moments of the `P₂` form within 0.2 % and every 10° band within
      0.005 equatorward of 60° (0.03 poleward).
- [ ] At obliquity 90° the poles' annual mean absorbed sunlight follows
      their insolation: the ice-free co-albedo of the polar cells is at
      least the equator's (the `P₂` form: 0.657 against 0.777), so their
      absorbed sunlight exceeds the equator's by at least the ratio of
      their insolation; documented comparison with the `P₂` form at 60° and
      90° (absorbed sunlight and temperatures) in the module and in the
      Result; images of `out/M4-10/` viewed and described.
- [ ] The zenith cosines match the reference within the stated
      tolerances and lie in `[0, 1]`; the insolation's golden is
      unchanged; the energy budget of every sub-step still closes to
      10⁻⁵ W/m²; the bounds tests pass, each planet that no longer
      settles diagnosed and recorded in their lists with a follow-up
      ticket; the golden fingerprints are recorded again after the other
      tests and the images, and `CHANGELOG.md` lists the change under
      World-breaking.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test --release -p kuniumi-sim climate
tools/.venv/bin/python tools/zenith_albedo_reference.py
```

Zonal profiles for obliquities 23.44°, 60° and 90° rendered to `out/M4-10/`.

## Result

The plan was approved by the owner as written ("計画どおり承認"): the law
linear in `μ̄` with `b` = 0.338 and `μ_r` = 0.644, `A` and `D⊕` fitted
again (the design's calibration section updated), and the planets of the
bounds tests that no longer settle diagnosed, recorded with their
diagnosis and given follow-up tickets.

What was done:
- `insolation`: `daily_zenith_cosine`, `MonthlyInsolation::zenith_cosines`
  and `monthly_zenith_cosine` (`G` from the private `day_integrals`, which
  computes `F` by the operations of `day_integral`, so the insolation and
  its golden value are unchanged); the module documentation's "The zenith
  angle"; three tests.
- `energy_balance`: `EnergyBalanceParams::albedo_zenith_slope` and
  `reference_zenith_cosine` (`ALBEDO_ZENITH_SLOPE` = 0.338,
  `REFERENCE_ZENITH_COSINE` = 0.644) instead of `albedo_p2`;
  `free_albedo` and `albedo` take the cosine; `EnergyBalance::new` takes
  the monthly cosines and keeps each month's ice-free albedo
  (`free_albedos`); `albedos`, `substep` and `month` take the month; the
  warm start and the balance step use each month's albedo
  (`annual_sunlight_w_m2`); `EARTH_OLR_AT_FREEZING_W_M2` = 202.8 and
  `EARTH_DIFFUSION_W_M2_K` = 0.579; the module documentation ("Albedo",
  "The zenith angle law", "Earth calibration", "Limitations"); six new
  tests, the existing ones with the month.
- `sea_surface`: the month's albedo in the sub-steps, the warm start's
  annual sunlight, `EARTH_OLR_AT_FREEZING_W_M2` = 202.8 and
  `EARTH_DIFFUSION_W_M2_K` = 0.502, "Calibration"; `coupling` and the
  test worlds pass the cosines.
- `lib.rs` (bounds tests): `UNSETTLED_ENERGY_BALANCE` = [16],
  `UNCONVERGED_CURRENTS` = [105] (new; `assert_warned` knows the
  currents' warning), `UNSETTLED_ORBITS` = [0, 3], each with its
  diagnosis.
- `tools/zenith_albedo_reference.py` (and `tools/README.md`).
- Golden values recorded again (below), `CHANGELOG.md` (Added, Changed,
  World-breaking), `docs/design-ja.md` (構成要素, 較正), follow-up tickets
  M4-24 and M4-25 in `docs/tickets/README.md`.

Deviations from the plan:
- Tolerances measured: the monthly cosines are within 10⁻⁴ of the
  script, and their error times the month's insolation within 0.01 W/m²
  (planned: 2 × 10⁻⁵): the midpoint rule's error is largest in a month
  whose light starts or ends within it (7.5 × 10⁻⁵ at 90° N, month 6, of
  the eccentric orbit, under 5.6 W/m²; Earth 1.4 × 10⁻⁵). The moments of
  the ice-free absorbed sunlight are checked within 10⁻³ on the grid
  (planned 0.2 %; measured 4.7 × 10⁻⁵ and 3.2 × 10⁻⁵).
- Test names: the exact values of the law are in
  `snow_and_ice_ramp_the_albedo` and the per-cell values in
  `the_model_takes_each_month_s_free_albedo` (planned together as
  `the_free_albedo_follows_the_zenith_cosine`); two panics added.
- `a_planet_leaves_an_unstable_periodic_state`: its closed-form balances
  computed again for `A` = 202.8 W/m² (296.011, 241.188 and 268.050 K
  instead of 296.059, 241.236 and 268.039 K).
- The bounds planets: corner 16 is not the slow ice front that the plan
  expected but a cycle of the balance step (its own years settle in 207
  years, after `MAX_YEARS`); recorded with follow-up M4-24. Corner 105 is
  BiCGSTAB stagnating on the all-ocean basin (without a cap its solves
  over the coupled loop take 284 to 31,181 iterations, against the cap of
  20,000); recorded with follow-up M4-25. Review round 1 found that
  neither is a known limit of the plan's kind (an ice front, a climate
  beyond the models' range) but a numerical regression on a valid
  planet; asked after round 1, the owner chose to record both as
  warnings in the lists of `lib.rs` and leave their fixes to M4-24 and
  M4-25 ("記録して後続へ"). Orbit 0 (14,000 K, 12,000 m/s winds) is beyond the models'
  range, as orbit 3; recorded with follow-up M4-26 (review round 1).
- Size: about 1,600 lines changed besides this ticket (the tests 700, the
  models and their documentation 500, the script 230, documents and
  goldens the rest), against the estimated 1,000.

Evidence:

1. Criterion 1 (release, the 40,000-cell world of seed 1, with the
   refitted `A` and `D⊕`):
   `earth_zonal_means_match_north_and_coakley` (the energy balance alone):
   every 10° band within 3.86 K of North & Coakley's fit (worst 60–50° S),
   `P₂` −27.966 K, mean surface 288.043 K (sea level 290.078 K);
   `earth_is_calibrated_with_its_coupled_climate`: 287.991 K, −27.973 K,
   worst band 4.44 K, 1.0957 m/yr, 24.93 kg/m² (with the `A` and `D` of
   `main` 288.284 K, −27.562 K, 5.10 K); `earth_stays_calibrated` (the sea
   surface with the winds of the energy balance alone, ±1 K): 287.718 K,
   −28.751 K, 3.45 K; `earth_keeps_its_planetary_albedo`: 0.2998 (the
   coupled state 0.3001; M4-02's 0.3006). Debug builds (10,000 cells,
   another world): 5.19 K, −25.828 K, 288.698 K; coupled 288.919 K and
   −24.663 K, 0.956 m/yr, 23.9 kg/m²; planetary albedo 0.2963.
   `the_zenith_law_keeps_earth_s_albedo` (40,000 cells): the global mean
   of the ice-free absorbed sunlight 253.4064 against the `P₂` form's
   253.3945 W/m² (4.7 × 10⁻⁵), the `P₂` moment −142.7172 against
   −142.7126 W/m² (3.2 × 10⁻⁵); each 10° band's insolation-weighted
   annual ice-free albedo within 0.0045 of the form's equatorward of 60°
   (largest 0.0044 at 40–50°), +0.028 at 80–90°, and within 10⁻⁴ of the
   script's exact bands (both hemispheres).
2. Criterion 2: `the_poles_absorb_with_their_insolation_at_high_obliquity`
   (obliquity 90°, Earth's orbit otherwise): the poles' annual weighted
   cosine 0.785 399 (π/4), the equator's 0.616 849 (π²/16 = 0.616 850);
   the ice-free ocean's co-albedo 0.78479 and 0.72782 (script 0.784 793
   and 0.727 823; `P₂` form 0.657 and 0.777); the pole absorbs 1.6937
   times the equator's sunlight against an insolation ratio of 1.5708
   (`P₂` form 1.3282); on 10,000 cells the model's co-albedo poleward of
   80° 0.78381 against 0.72681 within 10° of the equator. The coupled
   climate (40,000 cells, seed 1; scratch probe, `main` for the `P₂` form
   with its `A` and `D`), annual means:

   | Obliquity | | `P₂` form | zenith law |
   | --- | --- | --- | --- |
   | 90° | 80–90° S: absorbed, surface | 284.4 W/m², 300.7 K | 338.3 W/m², 311.4 K |
   | 90° | 5° S–5° N | 201.9 W/m², 281.9 K | 188.1 W/m², 279.8 K |
   | 90° | global surface, planetary albedo | 288.22 K, 0.2990 | 290.23 K, 0.2861 |
   | 60° | 80–90° S | 246.6 W/m², 292.5 K | 280.4 W/m², 299.8 K |
   | 60° | 5° S–5° N | 246.0 W/m², 289.0 K | 231.9 W/m², 287.2 K |
   | 60° | global surface, planetary albedo | 289.04 K, 0.2943 | 290.46 K, 0.2850 |
   | 23.44° | zonal surface (2° bands) | | −0.45 to +0.15 K |

   Documented in the module ("The zenith angle law"). Images in
   `out/M4-10/` (scratch tool `<scratchpad>/m4_10_tool`, a Cargo package
   with path dependencies, from the probe's zonal tables of 2° bands):
   - `zonal_23.44.png`, `zonal_60.png`, `zonal_90.png`: three panels
     against latitude (−90° to 90°, lines every 30°): the annual surface
     temperature (180–340 K, lines every 20 K), the insolation (blue) and
     the absorbed sunlight (0–600 W/m², every 50), and the
     insolation-weighted ice-free albedo (0.1–0.5, every 0.05); black the
     `P₂` form (`main`), red the zenith law. At 23.44° the two lie on top
     of each other but for the caps, whose ice-free albedo is about 0.03
     higher with the law (under snow most of the year, so the absorbed
     sunlight hardly moves). At 90° the law's ice-free albedo is lowest at
     the poles (about 0.22 against 0.34) and highest at the equator (0.30
     against 0.27), the absorbed sunlight follows the shape of the
     insolation, and the poles are about 11 K warmer, the equator 2 K
     cooler. At 60° the same, smaller: the caps 7 K warmer.
   - `free_albedo_23.44_m01.png`, `_m07.png` (viridis 0.1–0.6): the
     darkest band where the sun is highest, at 10–40° of the summer
     hemisphere (about 0.22 over the ocean), brighter towards the winter
     pole, the polar night uniform at 0.48 (ocean) and 0.58 (land); land
     0.1 brighter than the sea everywhere. `free_albedo_90_m07.png`: the
     darkest near the north pole under the midsummer sun, brighter towards
     the south, the southern hemisphere in polar night;
     `zenith_cosine_90_m01.png` (0–1): zero north of about 40° N (polar
     night through month 1, as the declination stays below −50°), rising
     to about 0.9 at the south pole.
3. Criterion 3: `the_daily_zenith_cosine_is_the_insolation_weighted_mean`
   (nine script values within 10⁻¹², π/4, sin δ and polar night's 0
   exact); `monthly_zenith_cosines_match_the_time_integral` (largest
   difference 1.35 × 10⁻⁵ for Earth and 7.5 × 10⁻⁵ for the eccentric
   orbit, times the insolation 9.1 × 10⁻⁴ and 3.6 × 10⁻³ W/m²; the zeros
   of polar night exact); `zenith_cosines_lie_in_the_unit_interval` (six
   orbits at the corners, 721 latitudes, and every cell of 2,000);
   `golden_fingerprint_of_earth_insolation` unchanged
   (`0xFAD5_AA6F_8873_F9B5`); `energy_is_conserved_in_every_substep`
   largest imbalance 9.97 × 10⁻⁸ W/m². The bounds tests pass with the
   lists measured again (the scratch probe ran `stages::climate::climate`
   on the worlds of the Crust stage and read each model's flags, on this
   branch and on `main`): corner 16 (energy balance), corner 105
   (currents), orbit 0 (coupled sea surface), diagnosed as in "Deviations"
   and in `lib.rs`. The new tests fail on mutations (`G` without its
   `2AB sin h₀` term; `b` = 0). Golden fingerprints recorded again after
   the tests above passed and the images were viewed: the energy balance
   `0xFC6C_0021_4D38_7341`, the circulation `0xD7D2_9D99_F54A_59B3` (without
   the smoothing `0x0812_F5E2_4A6D_4F1C`), the currents
   `0x1AC6_1A80_A9A0_A56A`, the sea surface `0x4DD4_C2E1_3E10_F277`, the
   moisture `0xDA05_CEE5_D1F1_A586`, the relaxed loop
   `0x5D11_FE79_8A27_CC94` (and its last change, 8.144889403922662e-4 K,
   35 years, relaxed from year 5 as before), the Climate stage's 30
   fields and the CLI's 30 climate fields (its 16 others and the sea level
   unchanged). `CHANGELOG.md`: Added, Changed, World-breaking.
4. `mise run ci` passes: fmt, clippy `-D warnings`, rustdoc, the browser
   build's lints, 2,109 tests summed over debug and release (10 ignored,
   as before), MSRV 1.85.0, cargo-deny (advisories, bans, licences and
   sources ok).

Timing (release, this machine): the zenith cosines of 650,000 cells take
2.70 s (the insolation 2.10 s), computed twice in a run of the Climate
stage; Earth's energy balance on 40,000 cells is periodic in 7 years in
7.9 s and its coupled state in 12 years in 33 s, as before.

Review:
- Round 1 (CHANGES_REQUESTED; the reviewer re-ran `mise run ci`, the
  criteria's tests, the reference script, an own quadrature of `b` and
  `μ_r`, and probes of the branch and of `main` that reproduced the
  comparison tables and images). R1-1 blocker: orbit planet 0 was recorded
  without a follow-up ticket, against criterion 3 (fixed: M4-26). R1-2
  major: corners 16 and 105 are numerical regressions, not known limits
  of the approved plan (asked to the owner, who chose to record them with
  the follow-ups M4-24 and M4-25 and land). R1-3 minor: the planetary
  albedo after the fit is 0.2998, not 0.2990 (fixed in the module and the
  test). R1-4 minor: the design document's and the CLI test's statements
  that only Earth at e = 0.99 warns (fixed, with the coupling module's
  count). R1-5 nit: a long documentation line rewrapped. R1-6 nit: the
  feedback of the ramp in the diagnosis of corner 16 stated correctly
  (850 W/m² × 0.27 over 10 K ≈ 23 W m⁻² K⁻¹, eleven times `B`). R1-7 nit:
  the script's running time given as one to three minutes in both
  places.
- Round 2 (APPROVE; the reviewer re-ran `mise run ci`, the criteria's
  tests, the reference script, an own computation of `b` and `μ_r`, the
  coupled climate on snapshots of the branch and of `main`, and looked at
  the images; R1-1 and R1-3 to R1-7 resolved, R1-2 settled by the owner's
  decision). R2-1 nit: an assertion that the `P₂` form's ratio falls short
  of the insolation ratio was always true (removed; the ratio stays in the
  printout and the comment). R2-2 nit: the coupling module's sentence on
  orbit planet 0 names its 1,000 cells.

Follow-ups: M4-24 (the balance step that cycles on corner 16), M4-25 (the
solve of the currents of an ocean over the whole sphere), M4-26 (the
coupled climate of the hottest orbit of the bounds tests).

# M4-20 The coupled climate of the extreme planets

- Status: done
- Milestone: M4
- Depends on: M4-17
- Crates: kuniumi-sim (stages/climate: coupling, climate; the bounds tests of lib.rs)
- Design: 気候モデル / 計算の流れ
- Decisions: D-08, D-22

## Goal

The coupled loop of temperature, winds and currents (M4-17) settles on the
planets at the corners of the v1 range as it does on Earth, or the planets
where it cannot are understood and reported.

## Scope

In:
- Measured in M4-17 (the design's yearly loop with the thermal anomaly
  smoothed over 1,000 km; 1,000 cells, seed 3): the coupled sea surface of
  23 of the 128 corners of the bounds tests does not settle within its
  100 years (corners 17, 19, 25, 27, 48, 50, 64, 66, 68, 72–75, 83, 84,
  86, 89–92, 94, 96 and 112; `UNSETTLED_SEA_SURFACE` in
  `kuniumi-sim/src/lib.rs`), their last change from 0.03 K to 65 K:
  mostly the 10-bar planets, the 8-hour ones and the worlds without
  water, whose wind stress and coastal upwelling are many times Earth's.
  Nor does that of 3 of the 6 planets at the bounds of the orbit
  (`UNSETTLED_ORBITS`): obliquity 90° with e = 0.99 at 2,500 W/m², 8 h,
  0.1 bar (last change 403 K) and at 400 W/m², 240 h, 10 bar (0.51 K),
  and Earth at e = 0.99 with the perihelion at 270° (743 K), which
  settled before M4-17. The set depends on the seed (the review of M4-17:
  22 planets at seed 1, 17 at seed 2, the 3 orbit planets at every
  seed).
- Why they do not settle (below, "Diagnosis"), and a method that settles
  those whose loop only oscillates, without changing Earth: the drive of
  the circulation relaxed once the loop oscillates (below, "The plan").
- The cost of the method, and the cost of the planets that settle with it.
- The bounds tests' lists of the unsettled planets measured again.

Out:
- Earth's coupled state and its calibration (M4-17, M4-08): unchanged, bit
  for bit.
- The planets whose thermal equator jumps between two latitudes from one
  year to the next (below, "Diagnosis"; follow-up ticket M4-21, draft): a
  continuous thermal equator is a change of the circulation (M4-03).
- Earth at e = 0.99 (orbit planet 3), whose monthly winds reach 700 m/s:
  beyond the circulation's range of validity (M4-03, "Limitations");
  reported.

## Implementation notes

The ticket was `draft`; it left open whether the loops that do not settle
are oscillations of the iteration (which damping removes) or
instabilities of the coupled model itself, and the cost of a method at
650,000 cells. A scratch test module inside `stages/climate` (not
committed; a copy of `coupling::coupled_run` with options, checked equal
to it bit for bit on Earth and on corners 0, 17 and 48) measured them:
release, this machine, the planets of the bounds tests through the Crust
stage on 1,000 cells (seed 3 unless said otherwise), the energy balance
without currents as the start, as `stages::climate::climate` does.

Notation: `D_y` the monthly sea-level temperatures that drive the
circulation of year `y` (the winds, the upwelling and the currents of
that year), `T_y` the monthly mean sea-level temperatures that year
gives; M4-17's loop takes `D_{y+1} = T_y` (and `D_1` the energy balance
without currents). The residual `R_y = T_y − D_y`, its largest value over
the months and cells `r_y`, and `λ_y = ⟨R_y, R_{y−1}⟩ / ⟨R_{y−1}, R_{y−1}⟩`
with the area-weighted sum over every month and cell. In M4-17's loop
`r_y` is the change from the year before, its criterion.

### Diagnosis

The 26 planets (23 corners and 3 orbit planets), M4-17's loop, the last
20 of their 100 years:

| Class | Planets | What the loop does |
| --- | --- | --- |
| a two-year cycle | 72, 92, 94, orbit 1 (`λ` = −1.000, `r` constant; `T_y = T_{y−2}` within 0.07 K for 72, 10⁻⁷ K or better for the others); 73, 74, 75, 84, 86 (`r` nearly constant, `λ` between −0.5 and −1.4) | the temperatures alternate between two years |
| longer cycles | 90 and 112 (period 3: `r` repeats to 3 digits), 91 (period 4), 48, 50, 89 (`r` alternating between two values, `λ` alternating in sign) | |
| irregular | 17, 19, 25, 27, 64, 66, 83, 96, orbits 0 and 3 | `λ` between −1.4 and +3.4 from year to year |
| a constant change | 68 | `r` = 2.11 K every year, `λ` ≈ +0.6 |

Where the largest residual sits: mostly a coastal ocean cell within 1°
of the equator (cells 476, 480, 506, 514 of the 1,000-cell grid; some
years a cell within 5°) for 17, 19, 25, 27, 64, 66, 72–75, 83, 89–91:
the upwelling–wind feedback of M4-17, on the 8-hour, 13,000-km worlds
without water under 2,500 W/m² and on the 10-bar ones; an open-ocean
cell at 30–34° for 48, 50, 96, 112 (the planets with an ocean of
10²⁰ m³); a polar coastal cell (74.8–81.8°) for 68, 84, 86, 92, 94 and
orbit 1; cells at 2–47° for orbits 0 and 3.

The thermal equator `φ_T` (the warmest latitude of the zonal fit, capped
at ±60°, M4-03): its largest change from one year to the next, over the
months, in the last 20 years, is 120° (from one cap to the other) in
nearly every year for 68, 84, 86, 92, 94, orbits 1 and 3 (the first six
have 240-hour days and 10 bar: the summer pole is the warmest, and in a
month near an equinox the two poles take turns), 22–41° for 73, 75, 90,
91, 10–13° for 96 and 112, 2–8° for 64, 66, 72, 74, 83, 89, and at most
1° for 17, 19, 25, 27, 48, 50 and orbit 0.

So most of these loops are oscillations of the iteration, not of the
coupled model: the year's lag between the temperatures and the winds
made from them, with a gain below −1 (`λ` < −1) at the equatorial
coasts, makes the loop swing from one year to the next around a state
that exists. Relaxing the drive (below) settles them to it, a fixed
point of the same loop.
The thermal equator, though, is a discontinuous function of the
temperatures (the argument of a maximum): where two latitudes are almost
equally warm the loop may have no fixed point at all, the circulation of
either making the other warmer, and no relaxation settles it.

### Methods measured

Every planet of the bounds tests (Earth, 128 corners, 6 orbit planets),
seeds 3, 1 and 2; "settle" is the criterion of the plan (below): `r_y` and
the change from the year before both at most 10⁻³ K within 100 years.

| Method | seed 3 | seed 1 | seed 2 | Earth |
| --- | --- | --- | --- | --- |
| M4-17's loop | 109 of 135 (median 6 years, 90 % within 23) | 113 | 118 | 11 years |
| relaxed from year 2, `β` = 0.5 / 0.25 / 0.1 (26 planets) | 17 / 17 / 1 of the 26 (median 34 / 52 / 77 years) | | | 21 / 43 years instead of 11 (`β` = 0.5 / 0.25); 18 instead of 11 at 10,000 cells, seed 1 (`β` = 0.5) |
| `β` halved whenever `r` grows (26 planets) | 7 of the 26 | | | |
| Aitken's dynamic relaxation (Irons & Tuck 1969) from year 2 / after a growth (26 planets) | 13 / 13 of the 26 | | | |
| Anderson acceleration (Walker & Ni 2011; depth 3, mixing 0.5) after the first growth of `r` | 125, none lost | 123, none lost | 127, none lost | unchanged |
| relaxed after the first growth of `r`, `β` = 0.35 | 126, none lost | 121, 2 lost | 126 | unchanged |
| **relaxed after an oscillation, `β` = 1/3** (the plan) | **125, none lost** | **123, none lost** | **126, none lost** | **unchanged** |
| the same, `β` = 0.25 / 0.3 / 0.4 / 0.5 | 125 / 126 / 125 / 124 | 120 / 123 / 122 / 121 | 127 / 126 / 126 / 126 | unchanged |
| the plan with `φ_T` frozen once relaxed (diagnosis only) | 133 | 130 | 133 | unchanged |

"Lost": a planet that settles in M4-17's loop and not with the method.
Relaxing from the start slows Earth; relaxing after the first growth of
`r` slows the planets whose loop converges monotonically but slowly
(`λ` ≈ +0.9; seed 1: corner 122 from 53 years to more than 100, 92 from
56), because relaxation slows a monotone mode (below). The plan's trigger
needs a reversal of the residual as well as a growth, and loses no planet
at any of the three seeds.

With the plan at seed 3: the 109 planets that settled before still do
(median 6 years, 90 % within 25, the slowest 74; 13 of them relax, 10
of which take other years, from corner 0's 30 instead of 41 to corner
82's 35 instead of 19); 16 more settle (corners 17, 19, 25, 27, 48,
50, 64, 66, 72, 74, 83, 89, 90, 96, 112 in 17 to 71 years and orbit 0 in
80); 10 do not: corners 68, 73, 75, 84, 86, 91, 92, 94 and orbit 1,
whose thermal equator jumps (last change 0.12 to 32 K), and orbit 3
(140 K). The coupled years of the 135 planets sum to 2,872 instead of
3,693. Seed 1: 10 more settle, 12 do not (corners 8, 64, 66, 68, 70,
73, 75, 76, 78, orbits 0, 1, 3); seed 2: 8 more, 9 do not (84, 86, 92,
94, 96, 112, orbits 0, 1, 3).

Freezing the thermal equator of each month once the loop relaxes settles
all but corner 68 (which never oscillates) and orbit 3 at seed 3, which
shows that the thermal equator is what keeps most of the others from
settling; but their frozen thermal equators end up as far as 120° from
the warmest latitude of their own temperatures (corners 84, 86, 92, 94,
orbit 1; the fit 0.17 to 3.3 K below its maximum there; 24–27° and
0.6–0.75 K for 73 and 75), so the state is not that of the model: not proposed. Nor is a
thermal equator that follows its own local maximum from year to year
(found by ascent from the year before's; Earth unchanged): with the plan
it settles 19 of the 26 instead of 16, but leaves some on a lesser
maximum, up to 616 K below the warmest (orbit 0) and 4.1 K (corner 90).
Both belong to follow-up M4-21.

Earth never triggers the relaxation, so its loop is M4-17's bit for bit
(the prototype's result equal to `coupled_sea_surface`): 1,000 cells
seeds 1, 2, 3; 1,200 cells seed 5; 1,500 cells seed 3; 10,000 cells
seeds 1, 2, 3; 40,000 cells seeds 1, 2, 3 (10 to 12 years; `r` falls
every year); at 650,000 cells (seed 1) M4-17's sequence of `r` falls
every year (8.0, 2.1, 0.95, 0.20, 0.13, 0.036, 0.011, 0.0039, 0.0015,
8.4 × 10⁻⁴ K), so it cannot trigger either.

At 10,000 cells (seed 3; every planet of the bounds tests, about 20
minutes a run), M4-17's loop settles 120 of the 135 planets (median 7
years), the plan 122 (corners 73 and 122 added, none lost, 10 relaxed
that settled before; the years sum to 2,524 instead of 2,600). Of the 13
that do not settle, 11 are planets whose thermal equator jumps: with it
frozen once relaxed, all but orbits 0 and 3 settle (corners 48, 50, 84,
86, 92, 94, 96, 98, 112, 114, orbit 1), and with the ascent 8. So on finer
grids the thermal equator, not the iteration, is what keeps most of the
extreme planets from settling: the relaxation settles the swings of the
iteration at every resolution measured (16 more planets at 1,000 cells,
2 at 10,000), and follow-up M4-21 is the larger step there.

### The plan (for the owner's approval)

#### The relaxation (`stages/climate/coupling.rs`)

The loop of M4-17 is a fixed-point iteration `D_{y+1} = G(D_y)` (with the
sea surface's own state carried from year to year). Its relaxation,

```text
D_{y+1} = D_y + β (T_y − D_y)
```

has the same fixed points (`T = D`) and multiplies a linear mode of
factor `λ` by `1 − β (1 − λ)` instead: with `β` = 1/3 an oscillation
(`λ` < 0) decays for any `λ` > −5 and one of `λ` = −2 vanishes, but a
monotone mode slows (`λ` = 0.9 becomes 0.967): the relaxation of a
fixed-point iteration (W. R. Mann 1953, "Mean value methods in
iteration", Proc. Amer. Math. Soc. 4, 506–510), the under-relaxation of
partitioned coupled solvers (U. Küttler & W. A. Wall 2008, "Fixed-point
fluid–structure interaction solvers with dynamic relaxation", Comput.
Mech. 43, 61–72). So the loop runs as M4-17's until it oscillates, and
relaxed from then on:

- year `y` (≥ 3) *oscillates* when `r_y > r_{y−1}` (the residual grows)
  and `⟨R_y, R_{y−1}⟩ < 0` or `⟨R_{y−1}, R_{y−2}⟩ < 0` (it reversed this
  year or the year before; the second catches cycles of about four
  years, corners 48 and 50), the inner product area-weighted, summed
  over the months and cells in order (`R_1 = T_1 − D_1`);
- from the year after the first year that oscillates, every drive is
  `D_y + β (T_y − D_y)` with `RELAXATION` = 1/3; before it, `T_y` itself,
  bit for bit;
- the loop stops when `r_y` and the change from the year before
  (`max |T_y − T_{y−1}|`) are both at most `PERIODIC_TOLERANCE_K`, or
  after `MAX_YEARS`; before the relaxation they are the same number.
  The returned circulation, currents and upwelling are those of `D_y`,
  within the tolerance of `T_y`; `Temperature::change_k` is the larger of
  the two.

```rust
/// The factor β of the relaxation of the drive (module documentation,
/// "Relaxation"): 1/3.
pub const RELAXATION: f64 = 1.0 / 3.0;

pub struct Coupled {
    pub circulation: Circulation,
    pub currents: OceanCurrents,
    pub sea_surface: SeaSurface,
    /// The first year whose winds were made from relaxed temperatures
    /// (module documentation, "Relaxation"); none where the loop did not
    /// oscillate, as on Earth.
    pub relaxed_from_year: Option<u32>,
}
```

`coupled_run` (crate-internal) takes the factor as a parameter (1 for no
relaxation, the loop of M4-17, for the tests). Cost: per year two passes
over the twelve months (the inner product and the relaxed drive), well
under a second at 650,000 cells against the year's 80 s, and one more
set of twelve monthly fields held (62 MB at 650,000 cells). A planet
that settles after relaxing runs its years, about 80 s each at 650,000
cells (corners at 1,000 cells: 17 to 80 years instead of 100); one that
does not still runs 100 (about 2 hours).

The module documentation gets a section "Relaxation" (the trigger, the
factor and its effect on a mode, the sources) and its "Convergence" the
measured numbers; `climate.rs`'s module documentation the new count of
unsettled corners. `climate::climate` ignores `relaxed_from_year`
(`Climate` keeps its five fields); the warnings are unchanged.

#### Tests

New, `coupling/tests.rs`:
- `earth_s_loop_is_never_relaxed`: the shared Earth runs (10,000 cells
  in debug builds, 40,000 in release builds) have `relaxed_from_year`
  none; with the golden fingerprints of the stage and of the CLI's world,
  which do not change, Earth's loop is M4-17's bit for bit.
- `an_oscillating_loop_settles_relaxed`: corner 72 of the bounds tests
  (1,000 cells, seed 3; a two-year cycle of 57.5 K in M4-17's loop):
  with `RELAXATION` the loop relaxes, settles (measured 34 years) and
  its state is within the tolerance of a fixed point: the winds of the
  returned temperatures within a bound set from the measurement of the
  returned winds, and the loop run on to 10⁻⁵ K within 10⁻³ K of the
  returned state (as M4-17's criterion 2); in release builds also: with
  the factor 1 the loop does not settle in 100 years and its last two
  years' residuals agree (the cycle).
- `the_trigger_needs_a_growth_and_a_reversal`: the crate-internal
  function that decides, on hand-made residuals: a falling residual, a
  growth without a reversal (a monotone drift), a reversal without a
  growth, and a growth after a reversal the year before.

Changed: the bounds tests' `UNSETTLED_SEA_SURFACE` (corners 68, 73, 75,
84, 86, 91, 92, 94) and `UNSETTLED_ORBITS` (1, 3), measured again with
the implementation by the scratch tool of M4-17 (`m4_17_flags`, each
model's own flags through the public `stages::climate::climate`), and
`UNSETTLED_MOISTURE` likewise; their documentation (why the remaining
ones do not settle), and the comment on the debug builds' corners (none
of 17, 48, 64, 96, 112 warns any more).

#### Documents

`docs/design-ja.md`, 気候モデル / 計算の流れ: the paragraph of the
coupled loop says that the loop is relaxed once it oscillates, and that
the planets whose thermal equator jumps between two latitudes (240-hour
days with dense atmospheres) and Earth at e = 0.99 still do not settle
(M4-20, M4-21). `CHANGELOG.md`: Changed (the relaxation) and
World-breaking (the climate of the planets whose loop oscillates; no
Earth world changes). Follow-up M4-21 (draft) in
`docs/tickets/README.md`.

#### World-breaking

The climate of a planet whose coupled loop oscillates (in the bounds
tests 38 of the 135 planets at seed 3: the 16 that now settle, 9 of the
10 that still do not, and 13 of those that settled before); every
planet whose loop does not oscillate, Earth among them, keeps its world
bit for bit. No golden value changes (they
are all Earth's).

#### Size

About 350 changed lines without the ticket: the coupling module 120 (with
its documentation), tests 150, the bounds tests 40, documents 40.

### Questions for the owner

1. The method: the loop relaxed with `β` = 1/3 once it oscillates (a
   growth of the residual with a reversal in that year or the one
   before), Earth unchanged bit for bit; at 1,000 cells, seed 3, 16 of
   the 26 planets settle and none is lost (seeds 1 and 2: 10 of 22, 8 of
   17), at 10,000 cells 2 of 15. The planets whose thermal equator jumps
   (most of those left, and 11 of the 13 left at 10,000 cells) and orbit
   3 stay unsettled and reported; the thermal equator's continuity goes
   to follow-up M4-21 (draft, added with this ticket), a change of the
   circulation (M4-03) with its own choice of rule. Alternatives: the
   same and a continuous thermal equator in this ticket (larger; the two
   rules measured here are not proposed, above); relaxation after any
   growth (one more corner at seed 3, two planets lost at seed 1);
   Anderson acceleration (as many settle; more machinery and memory).

## Acceptance criteria

- [x] The planets of the bounds tests that settle (1,000 cells, seed 3)
      are those of M4-17 and at least 10 more of its 23 corners and 3
      orbit planets (measured with the plan: 16), the lists of the bounds
      tests measured again; at seeds 1 and 2 none is lost (scratch
      tool).
- [x] Earth's years and state unchanged: `earth_s_loop_is_never_relaxed`,
      and the golden fingerprints of the climate's modules, the stage and
      the CLI's world unchanged.
- [x] `an_oscillating_loop_settles_relaxed` and
      `the_trigger_needs_a_growth_and_a_reversal` pass.
- [x] The planets that still do not settle are described (thermal
      equator, orbit 3) in the module documentation, the bounds tests and
      follow-up M4-21.
- [x] `mise run ci` passes.

## Verification

```sh
cargo test --release -p kuniumi-sim coupling
cargo test --release -p kuniumi-sim bounds
```

The scratch tool of M4-17 (`m4_17_flags`) on the planets of the bounds
tests at seeds 3, 1 and 2.

## Result

What was done (as planned and approved by the owner):
- `stages/climate/coupling.rs`: the loop relaxed once it oscillates, the
  module documentation's new section "Relaxation" (the iteration, the
  trigger, the factor and its effect on a mode, the sources, why a
  thermal equator that jumps is not settled by it) and the measured
  numbers in "Convergence" and "Cost"; `RELAXATION` (1/3),
  `Coupled::relaxed_from_year`, the crate-internal `Control` (tolerance,
  most years, factor) of `coupled_run` and `oscillates`. The residual of
  the year before is kept and replaced in place, so the loop holds one
  more set of twelve monthly fields until it relaxes (the relaxed drive
  after).
- Tests (`coupling/tests.rs`): `earth_s_loop_is_never_relaxed`,
  `an_oscillating_loop_settles_relaxed` (exact years and first relaxed
  year), `the_trigger_needs_a_growth_and_a_reversal`,
  `golden_fingerprint_of_a_relaxed_loop` (corner 72's coupled state,
  `0xC793_26CF_E725_5548`; added in review round 1),
  `the_reversal_is_weighted_by_the_cells_areas` (review round 2); the two existing
  calls of `coupled_run` take a `Control`, and the fixed-point test of
  M4-17 checks that Earth's loop run on to 10⁻⁵ K is not relaxed either.
- The bounds tests (`lib.rs`): `UNSETTLED_SEA_SURFACE` 8 corners,
  `UNSETTLED_ORBITS` [1, 3], their documentation, the debug corners'
  comment; `UNSETTLED_MOISTURE` measured again, unchanged.
- `climate.rs`'s module documentation (step 2 of the stage, the count of
  unsettled corners); `climate::climate` ignores `relaxed_from_year`.
- Documents: `docs/design-ja.md` (気候モデル / 計算の流れ), `CHANGELOG.md`
  (Added, Changed, World-breaking), follow-up M4-21 (draft) in
  `docs/tickets/README.md`.

Deviations from the plan:
- `coupled_run` takes the tolerance, the years and the factor as one
  crate-internal `Control` rather than three more arguments (clippy's
  limit of seven arguments; no lint is allowed in library code).
- `an_oscillating_loop_settles_relaxed`'s release-only check of the loop
  without relaxation asserts that it runs its 100 years unconverged with
  a last change above 10 K (measured 57.54 K, the two-year cycle of the
  prototype) instead of comparing the residuals of its last two years,
  which the loop does not return.
- A golden fingerprint of a relaxed loop (corner 72), not in the plan:
  every other golden value is Earth's, whose loop never relaxes (review
  round 1).

Evidence:

1. Bounds tests (release: Earth, the 128 corners and the 6 orbit
   planets; debug: Earth, 16 corners and the orbit planets) pass with the
   lists measured again by a scratch tool (`<scratchpad>/m4_20_flags`,
   after M4-17's `m4_17_flags`: each model's own flags through the public
   `stages::climate::climate`, and `relaxed_from_year` of
   `coupled_sea_surface`, equal to the climate's sea surface), seed 3:
   the coupled sea surface of corners 68, 73, 75, 84, 86, 91, 92, 94 and
   orbits 1 and 3 does not settle (last change 2.11, 6.72, 31.9, 6.64,
   3.03, 10.1, 7.06, 6.47, 0.119 and 141 K); every other planet's does,
   among them the 109 that settled in M4-17 and 16 more (corners 17, 19,
   25, 27, 48, 50, 64, 66, 72, 74, 83, 89, 90, 96, 112, orbit 0); 38 of
   the 135 relax; the relaxed planets that settle take 15 to 80 years.
   Of the planets that settled before, 13 relax at seed 3 (10 take other
   years, the most corner 6, 8 to 19 years, and 82, 19 to 35), 15 at seed
   1 (corner 88 from 19 to 54 years, 90 from 25 to 69, 84 from 15 to 34)
   and 12 at seed 2 (corner 6 from 10 to 29, 88 from 21 to 58): up to
   about three times their years. Some newly settled planets at seed 1
   use most of the 100 years (corners 74: 98, 91: 90, 72: 81).
   The energy balances, the currents' solves and the moisture's flags are
   as in M4-17 (`UNSETTLED_MOISTURE` unchanged). Seeds 1 and 2 (the same
   tool): the coupled sea surface of 12 and 9 planets does not settle,
   each of them unsettled in M4-17's loop too (22 and 17 there): none
   lost.
2. `earth_s_loop_is_never_relaxed` (10,000 and 40,000 cells: 11 and 12
   years, `relaxed_from_year` none). No golden value changed and every
   golden test passes (the climate modules', the stage's and the CLI's
   world, all Earth's); `the_coupled_state_is_within_the_tolerance_of_its_fixed_point`
   prints M4-17's numbers (40,000 cells: 17 years to 10⁻⁵ K,
   1.627 × 10⁻⁴ K away, the winds 1.853 × 10⁻⁴ m/s).
3. `an_oscillating_loop_settles_relaxed` (corner 72, 1,000 cells, seed
   3): relaxed from year 5, settles in 34 years (last change
   8.51 × 10⁻⁴ K), as in the prototype; run on to 10⁻⁵ K (48 years) it
   stays 2.36 × 10⁻⁴ K away; the winds of the returned temperatures
   9.5 × 10⁻⁵ m/s from the returned winds (bound 2 × 10⁻³ m/s); without
   relaxation 100 years, last change 57.54 K.
   The reported last change is exactly the residual, 8.507919344822312 ×
   10⁻⁴ K (the year's change 8.7 × 10⁻⁵ K).
   `the_trigger_needs_a_growth_and_a_reversal` and
   `the_reversal_is_weighted_by_the_cells_areas` pass.
   `golden_fingerprint_of_a_relaxed_loop`: `0xC793_26CF_E725_5548` in
   debug and release builds, recorded after the tests above passed and
   the maps of evidence 5 were viewed.
4. Described: `coupling`'s "Relaxation" and "Convergence", the bounds
   tests' `UNSETTLED_SEA_SURFACE` and `UNSETTLED_ORBITS`, `climate.rs`,
   `docs/design-ja.md` and M4-21.
5. Generated worlds (`kuniumi generate --seed 3 --planet
   out/M4-20/corner<N>.toml`, main's binary against the branch's;
   `out/M4-20/`), compared and drawn by the scratch tool's `compare` and
   `maps`:
   - corner 64 on 10,000 cells, which settles both ways (39 years
     without relaxation, 32 relaxed): the two states agree within
     6.6 × 10⁻⁴ K (monthly temperatures, RMS 3 × 10⁻⁵ K), the winds within
     2.3 × 10⁻⁴ m/s: the same fixed point. `maps_c64/diff_temperature_m01.png`
     (±10⁻³ K): white but one equatorial coastal spot of −6.6 × 10⁻⁴ K.
   - corner 72 on 1,000 cells (a two-year cycle on main, warned;
     settled on the branch, no warning): `maps_c72/before_temperature_m09.png`
     and `after_…` (one range, 99–258 K): the same pattern, but the warm
     spot of one cell on the equator (cell 476, 0.7°, longitude −66.9°) of
     main's last year, yellow, is gone; `diff_temperature_m09.png` (±30 K):
     −29 K at that cell with a halo of about 1,000 km, a lesser one at the
     second equatorial site (about +107°), and up to about +1 K elsewhere, mostly
     at high latitudes; the global annual mean 164.735 K against
     164.737 K. `before_wind_east.png`, `after_…` (±8 m/s): the five
     bands of an 8-hour day, the same (largest difference 1.5 m/s).
6. `mise run ci` passes (again after the fixes of review rounds 1 and 2): fmt,
   clippy `-D warnings`, rustdoc, the browser build's lints, 2,067 tests
   summed over debug and release, MSRV 1.85.0, cargo-deny (advisories,
   bans, licences, sources ok); 9 min 28 s. The debug tests of
   `kuniumi-sim` take 275 s (303 s in M4-17: the debug corners 17, 48,
   64, 96 and 112 no longer run their 100 years), the release ones 95 s.

Review:
- Round 1 (CHANGES_REQUESTED; the reviewer re-ran `mise run ci`, wrote
  its own flags tools against `main` and the branch at seeds 1–3, hashed
  every planet's coupled state (bit for bit `main`'s wherever the loop
  does not relax), regenerated and viewed the maps, and mutated a copy).
  R1-1 major: no golden or exact value covered the relaxed loop (every
  golden value is Earth's), so mutations of the relaxed drive passed:
  fixed with exact years and first relaxed year in
  `an_oscillating_loop_settles_relaxed` and the new
  `golden_fingerprint_of_a_relaxed_loop`. R1-2 minor: corner 68 was said
  to be relaxed; it never is (its residual never reverses): fixed in
  `lib.rs` and `coupling`. R1-3 minor: the stage's step 2,
  `Climate::circulation` and the warnings did not describe the relaxed
  loop's criterion and last change: documented (the warning text kept, as
  planned). R1-4 minor: the slowdown at seeds 1 and 2 was not stated:
  added to "Cost" and evidence 1. R1-5 nit: the fixed-point test now
  asserts that Earth's loop is not relaxed. R1-6, R1-7 nits: the map's
  longitude and criterion 1's wording corrected.
- Round 2 (APPROVE; mutations of round 1 now fail, the fingerprint
  reproduced by the reviewer's own tool). R2-1 minor: nothing checked
  that the reported change is the larger of the change and the residual:
  `an_oscillating_loop_settles_relaxed` asserts it exactly. R2-2 nit: the
  area weighting of the reversal was untested: new
  `the_reversal_is_weighted_by_the_cells_areas`. R2-3 nit: the design
  document's stopping rule now names the residual. R2-4 nit: "has no
  fixed point" softened to "may have".

Follow-ups: M4-21 (draft): a continuous thermal equator.

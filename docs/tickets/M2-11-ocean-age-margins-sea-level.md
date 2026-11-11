# M2-11 Ocean-floor age, passive margins and sea level

- Status: done
- Milestone: M2
- Depends on: M2-04
- Crates: kuniumi-sim (stages/tectonics), kuniumi-cli (tests)
- Design: 地殻から標高へ; プレートテクトニクス / 境界ごとの処理; プレートテクトニクス / プレート速度と分裂; 検証方法
- Decisions: none

## Goal

Decide, before M2-05, how the one-shot planet reaches Earth's land fraction
of about 29 % from 40 % continental crust: the ocean floor of M2-04 is much
older than Earth's, and no continental margin is submerged.

## Scope

In:
- The finding of M2-04 (its Implementation notes, "Prototype
  measurements"): with the age of M2-04 (distance to the nearest ridge of
  the plate ÷ its half spreading rate, capped at 180 Myr) the area-weighted
  mean age of the ocean floor is 79–144 Myr (Earth: about 60 Myr) and
  23–68 % of it sits at the cap; with M2-05's formulas (Parsons & Sclater
  depth, Airy isostasy, 1.335e18 m³ of water) the sea level falls 0.7–1.5 km
  below the datum and the land fraction is 40.0 % for every seed (seeds 1–5
  at 100,000 cells, 1–10 at 10,000).
- The age is discontinuous inside plates: cells take their age from the
  nearest ridge cell, whose rate may differ greatly from its neighbours'
  (a ridge's rate falls to zero at its relative pole). Measured in the
  review of M2-04 at 100,000 cells, seed 1: 1,955 pairs of neighbouring
  oceanic cells of the same plate differ by more than 20 Myr, seams
  thousands of kilometres long reach 127–159 Myr, and next to a ridge cell
  spreading at 0.10 cm/yr the age jumps from 0 to 180 Myr in one step, a
  step of about 3.7 km in M2-05's sea floor. Taking for each cell the ridge
  it is reached from first by growth at each ridge cell's own half rate
  (the "time-nearest" variant of M2-04's prototype) would largely remove these
  seams; smoothing the rates along the ridges is another option.
- Even Earth's triangular area–age distribution would give a mean depth of
  4.9 km against the 4.4 km that the water fills on 60 % of the surface, so
  the design's thinned continental margins (地殻から標高へ / 大陸: 25–30 km,
  submerged as continental shelves) are needed as well.
- Options to decide with the owner: thin passive margins (continental cells
  near oceanic crust of the same plate) and continental rifts; faster
  ridges, for example plate speeds that grow with the length of their
  subduction boundaries (design, プレート速度と分裂), which changes M2-02
  and is world-breaking; another age model; or a combination.

Out:
- The time-evolving tectonics of M6, which produces margins and ages by
  itself.

## Implementation notes

Refined from `draft`. The options of the Scope were measured with a
throwaway prototype in the scratchpad (not committed) on the real grids,
plates, continents, boundaries and crust of M2-02 to M2-04 (Earth's
parameters, 12 plates), with M2-05's formulas written into the prototype:
Parsons & Sclater depth switched where the two curves meet (26.27 Myr,
4,294 m), Airy isostasy (2,800 and 3,300 kg/m³, 35 km → +300 m), trench
cells at −8 km, and the sea level by bisection on 1.335e18 m³.

### Prototype measurements

Land fraction, sea level `z` (relative to the datum), area-weighted mean
age of the ocean floor, share of the ocean floor at the 180 Myr cap,
"seams" (pairs of neighbouring oceanic cells of one plate whose ages differ
by more than 20 Myr, seed 1) and the mode of land elevation above sea level
(100 m bins; Earth: about +0.1 km, mean +0.84 km). 100,000 cells, seeds
1–5, unless noted.

| Option | Land | `z` (km) | Mean age (Myr) | At cap | Seams | Land mode (km) |
| --- | --- | --- | --- | --- | --- | --- |
| M2-04 as is | 40.0 % | −0.73 … −1.46 | 82–144 | 23–68 % | 1,955 | +1.05 … +1.75 |
| a. time-nearest age | 40.0 % | −0.69 … −1.41 | 78–138 | 19–59 % | 396 | +0.95 … +1.75 |
| b. a, ranked onto Earth's area–age distribution | 40.0 % | −0.53 … −0.55 | 60.0 | 0 % | 82 | +0.85 |
| c. a with plate speeds × 2 | 40.0 % | −0.17 … −1.10 | 49–110 | 1–39 % | 83 | +0.45 … +1.45 |
| c'. a with plate speeds × 3 | 40.0 % | +0.14 … −0.88 | 33–95 | 0–31 % | 45 | +0.15 … +1.15 |
| c''. a, speeds 10 − 9 × (continental share of the plate) cm/yr | 40.0 % | −0.61 … −1.28 | 69–122 | 2–38 % | 94 | — |
| margins only (M2-04 age), 10 km → 35 km over 600 km | 29.9–31.6 % | −0.96 … −1.62 | 82–144 | 23–68 % | 1,955 | +1.25 … +1.95 |
| a + margins 10 km / 600 km | 29.8–31.3 % | −0.92 … −1.58 | 78–138 | 19–59 % | 396 | +1.25 … +1.85 |
| **b + margins 15 km / 600 km** | **28.2–30.0 %** | **−0.71 … −0.76** | **60.0** | **0 %** | **82** | **+1.05** |
| b + margins 15 km / 300 km | 33.1–34.2 % | −0.63 … −0.68 | 60.0 | 0 % | 82 | +0.95 |
| c + margins 15 km / 600 km | 28.6–30.6 % | −0.40 … −1.23 | 49–110 | 1–39 % | 83 | +0.75 … +1.55 |

The prototype's variants:

- *Time-nearest age*: growth from every ridge cell at its own half spreading
  rate; a cell takes the first front that reaches it, fronts pass only
  through cells they reached first (as the partition's growth). At 650,000
  cells, seed 1, the seams fall from 2,034 to 0: the first-arrival time
  changes by at most one cell's travel time between neighbours. It changes
  the mean age by less than 10 %, as M2-04 found.
- *Earth's area–age distribution*: the oceanic cells reached by a ridge
  sorted by their time-nearest age, then given the age at the middle of
  their area in Parsons' (1982) triangular distribution
  `dA/dt ∝ 1 − t/t_max`, `t_max` = 180 Myr (mean 60 Myr); cells of oceanic
  plates without a divergent cell stay at 180 Myr (up to 20 % of the ocean
  at 10,000 cells, none at 100,000 cells for seeds 1–5). The spatial pattern
  still comes from the ridges; only the scale is Earth's.
- *Faster plates* (M2-02 change): the speed of every plate multiplied, or
  set from its continental share (design, プレート速度と分裂: 大陸の多いプレ
  ートは遅くします). Speeds from the trench share of each plate's boundary
  (after Forsyth & Uyeda 1975, "On the relative importance of the driving
  forces of plate motion", GJRAS 43, 163–200; one pass) gave 54–110 Myr at 10,000 cells
  and 82–126 Myr at 100,000 cells. Earth-like mean ages need about
  twice to three times today's speeds (fastest plates 20–30 cm/yr, outside
  the design's 1–10 cm/yr), and the sea level then varies by about 1 km
  between seeds.
- *Margins*: every continental cell's thickness
  `T = T_c + (35 km − T_c) · min(1, d/W)` plus M2-04's orogen thickening,
  with `d` the distance (m) from the coastline through continental cells of
  the plate. Sources (variant "nc", the one in the table): continental
  cells that are not convergent boundary cells and either touch an oceanic
  cell of an ocean region (oceanic cells connected through any plate) of at
  least 0.5 % of the surface, or are divergent boundary cells (continental
  rifts); each starts at half the centre distance to that neighbour (the
  coastline is the Voronoi edge), so the profile does not depend on the cell
  size. Only coasts facing ocean of the same plate (passive coasts in the
  narrow sense) give 1–2 points more land; all coasts including convergent
  fronts give 0.4 points less, but thin the arcs of M2-04.

Findings:

1. Every option needs thinned margins. All continental crust stands at
   +300 m or more and the sea level stays below +0.2 km, so without them the
   land fraction is the continental fraction, 40.0 %, for every seed and
   every age model.
2. The margins set the land fraction; the age model sets the sea level and
   so the hypsometry. With the same margins the land fraction differs by
   1–2 points between age models, but the sea level by up to 1 km, and with
   it the height of the continental interiors above the sea.
3. Margins must be wide. With the sea level 0.5–1.6 km below the datum,
   crust floods only where it is thinner than 22–30 km, so the design's 25–30 km shelves are
   partly dry, and the 11 % of the surface that must flood to reach 29 %
   needs thinning over about 600 km from the coastline (of which about
   400 km flood). Earth's submerged continental crust (shelf, slope and part
   of the rise, about 11 % of the surface, 5.6e13 m²) lies along margins of
   which the passive ones alone are about 105,000 km long (Bradley 2008,
   "Passive margins through earth history", Earth-Science Reviews 91, 1–26),
   so it is a few hundred kilometres wide on average; the coastlines here
   are 230,000–295,000 km long at 100,000 cells, rough at the cell scale,
   and the bands of neighbouring coasts overlap.
4. Resolution: for b + margins 15 km / 600 km the land fraction is
   28.2–30.0 % at 10,000 cells (seeds 1–10), 28.2–30.0 % at 100,000 cells
   (seeds 1–5) and 30.2–30.8 % at 650,000 cells (seeds 1–2).
5. Even Earth's own area–age distribution leaves the sea level 0.5 km
   below the datum (mean depth 4,877 m on the oceanic 60 %, against the
   4,362 m that the water fills): Parsons & Sclater's depths are for
   sediment-free crust. The depth–age curve GDH1 (Stein & Stein 1992, "A
   model for the global variation in oceanic depth and heat flow with
   lithospheric age", Nature 359, 123–129) gives nearly the same (4,845 m),
   so changing the depth curve is not a lever.

### Plan (recommended option, b + margins)

Chosen by the owner (Decision below). No change to the plates, continents
or boundaries; `crust_age`, `crust_thickness` and possibly `trench` change
(world-breaking).

Ocean-floor age (in `stages::tectonics::crust`):

1. `ridge_times(grid, plates, radius_m, bounds) -> Vec<f64>` (s): the
   first-arrival time from the divergent cells of each plate. Multi-source
   Dijkstra: every divergent cell `s` enters a `BinaryHeap` at time 0 and
   path length 0; the entry of least time is taken (ordered by time with
   `total_cmp`, then cell, then source); an unsettled cell is settled with
   the entry's time, path length and source, and pushes each unsettled
   neighbour `j` of the same plate at path length `ℓ' = ℓ + d_ij R` and
   time `2ℓ' / (−c_s)`. Paths through any cell of the plate, as M2-04.
   `+∞` where no ridge of the plate reaches.
2. `earth_area_age(grid, continental, times_s) -> Vec<f64>` (s): the
   oceanic cells with a finite time, sorted by (time, cell); with `A` their
   area (summed in that order) and `G_i = (A_<i + a_i/2) / A`,
   `age_i = t_max (1 − √(1 − G_i))`, `t_max` = `MAX_CRUST_AGE_S` (the inverse
   of the triangular cumulative distribution `F(t) = 1 − (1 − t/t_max)²`,
   Parsons 1982, "Causes and consequences of the relation between area and
   age of the ocean floor", JGR 87(B1), 289–302). Ridge cells (time 0) get
   +0; unreached oceanic cells and continental cells `t_max`.
3. The trench rule of M2-04 compares the uncapped times of step 1 (the
   ordering is the same as of the ages, except that it keeps `+∞`).

Continental margins:

4. `margin_distances(grid, plates, continental, bounds, radius_m) -> Vec<f64>`
   (m): sources and starting distances as in the prototype's variant "nc"
   above (`MARGIN_MIN_OCEAN_FRACTION` = 0.005 of the surface for an ocean
   region, so the pockets of M2-10 make no margin); Dijkstra through
   continental cells of the plate; `+∞` where none reaches.
5. Thickness of continental cells:
   `T = T_c + (T₀ − T_c) · min(1, d / W) + max(ΔT_coll f_coll, ΔT_arc f_arc)`,
   `T₀` = 35 km, `T_c` = `MARGIN_COAST_THICKNESS_M` (15 km, crust thinned
   towards the continent–ocean transition), `W` = `MARGIN_WIDTH_M`
   (600 km). The two values are to be fixed during implementation so that
   the mean land fraction of the ten seeds is 28–31 %, and stated with
   Earth's margins as their order-of-magnitude source. Uplift unchanged.
6. `crust` returns the new age and thickness; its signature and `Crust`
   stay. The design document (地殻から標高へ / 海洋底 and 大陸) is updated
   with the age model and the margin profile, with the owner's agreement.

Measuring the land fraction: in the tests, with a reference of M2-05's
formulas as in the prototype (std functions allowed in tests; no trench
deepening, which M2-05 decides). M2-05 then implements the stage for real
and can reuse the reference as its independent check. The order of the
tickets stays.

### Question for the owner

How should the one-shot planet reach Earth's land fraction?

- Recommended: **b + margins** as planned above.
- a + margins: ages from the ridges only (time-nearest, no Earth
  calibration), margins 10 km / 600 km; land 29.8–31.3 %, sea level 0.9–1.6
  km below the datum, interiors 1.2–1.9 km above the sea.
- c + margins: as a with plate speeds doubled (changes M2-02, fastest
  plates up to 20 cm/yr); land 28.6–30.6 %, sea level 0.4–1.2 km below.
- Margins only, M2-04's ages kept (seams stay); land 29.9–31.6 %.

### Decision (owner, 2026-11-11)

The owner chose the recommended option, **b + margins** with 15 km at the
coastline returning linearly to 35 km 600 km inland, from coasts other than
convergent fronts and from continental rifts; ages from the time-nearest
arrival ranked onto Earth's triangular area–age distribution; the plates
unchanged. The choice approves the plan above, including the update of the
design document (地殻から標高へ / 大陸 and 海洋底).

## Acceptance criteria

- [x] Time-nearest ages: on M2-04's 32 two-plate cases with every
      divergent cell given one spreading rate, `ridge_times` equals
      `2d/(−c)` of `nearest_sources` within 1e-12 relative; on the real
      worlds (seeds 1–3 at 10,000 cells, seed 1 at 100,000) every reached
      cell's time is `2ℓ/(−c_s)` of its own path and no same-plate
      neighbour could be reached earlier through it
      (`t_j ≤ 2(ℓ_i + d_ij R)/(−c_{s_i})`, 1e-12 relative); ridge cells
      +0; unreached cells exactly those of plates without a divergent cell.
- [x] Earth's area–age distribution: on the same worlds the area-weighted
      cumulative distribution of the reached oceanic cells' ages differs from
      `1 − (1 − t/t_max)²` by at most the largest cell area over `A` at
      every cell's age; the ages keep the order of the times; mean
      60 ± 0.5 Myr where every oceanic plate has a ridge.
- [x] Seams at 100,000 cells, seed 1: fewer than 200 pairs (1,955 with
      M2-04).
- [x] Margins: a synthetic setup (one continental cap per plate on two
      plates, a straight coast inside each plate, a rift and a convergent
      front on the boundary) matches the profile worked out by hand with
      reference distances from a Bellman–Ford search written in the test
      (1e-9 relative); convergent fronts and ocean regions below 0.5 % of
      the surface make no margin; continental thickness on the real worlds
      within `[T_c, 70 km]`, oceanic exactly 7 km.
- [x] Land fraction with Earth's parameters, seeds 1–10 at 10,000 cells,
      measured with the test reference of M2-05's formulas: every seed in
      20–40 %, mean of the ten in 28–31 %; also reported for 100,000 cells
      (seeds 1–5) and 650,000 cells (seeds 1–10).
- [x] Maps of `crust_age` and `crust_thickness` and a land–sea map from the
      reference formulas (100,000 and 650,000 cells, seed 1) viewed and
      described.
- [x] Golden fingerprints of the changed fields re-recorded after the tests
      above pass and the maps are viewed; world-breaking line and
      `CHANGELOG.md` entry.
- [x] `mise run ci` passes.

## Verification

```sh
cargo run --release -p kuniumi-cli -- generate --seed 1 --cells 100000 --out out/M2-11/100k
cargo run --release -p kuniumi-cli -- generate --seed 1 --out out/M2-11/650k
```

## Result

What was done:
- `kuniumi-sim`, `stages::tectonics::crust`: `ridge_times` and
  `RidgeTimes` (first arrival of fronts from every divergent cell at its
  own half spreading rate, heap entries ordered by time, cell, ridge and
  path), `earth_area_age` (ranking onto Parsons' triangular distribution),
  `margin_distances` (open-ocean coasts except convergent fronts, and
  continental rifts, starting half a centre distance from the coastline),
  the private `margin_profile`, and the constants
  `MARGIN_COAST_THICKNESS_M` (15 km), `MARGIN_WIDTH_M` (600 km) and
  `MARGIN_MIN_OCEAN_FRACTION` (0.005), all re-exported from `tectonics`.
  `nearest_sources` now calls a private `grow` that takes starting
  distances (its results are unchanged: the golden `uplift_rate`, which
  depends on it, is the same). `crust` uses the times for the trench rule,
  the ranked ages and the margins; its signature and `Crust` stay.
- Tests: `crust/age_tests.rs` (7), `crust/margin_tests.rs` (6) and the
  M2-04 tests of `crust/tests.rs` adapted to the new ages and margins, with
  the shared references `reference_from` (start distances),
  `reference_times`/`best_offer` and `reference_margin_distances`;
  `kuniumi-cli` tests: the 1,000-cell range (15–70 km) and the golden
  values.
- `docs/design-ja.md` (地殻から標高へ / 海洋底 and 大陸, Parsons 1982 in the
  references), `CHANGELOG.md` (Changed, World-breaking), the module
  documentation, and the follow-up draft M2-12 in the README.

Deviations from the plan:
- The first criterion compared the single-rate case with `nearest_sources`;
  the test uses the independent `reference_nearest` instead, and the
  real-world check is stronger than planned: every reached cell holds
  exactly the least offer of its neighbours' own final labels (the
  equations whose solution is unique, since offers are later than the
  labels they come from), besides matching `reference_times` bit for bit.
  A first version of `reference_times` kept a label until a better offer
  came (Bellman–Ford), which leaves offers of fronts that later lost their
  cell; it disagreed with `ridge_times`, and it now recomputes each label
  from its neighbours at every visit.
- The prototype set trench cells to −8 km; the tests' reference of
  M2-05's formulas leaves trenches undeepened, as planned (M2-05 decides).
- `T_c` = 15 km and `W` = 600 km stayed as in the prototype: the ten-seed
  mean is 28.9 % at 10,000 cells, inside the 28–31 % of the plan.

Evidence:

1. Time-nearest ages: `ridge_times_on_the_octahedron` (by hand: the faster
   of two ridges wins cells equidistant from both, a path of two steps,
   ties to the lower ridge, a plate without a ridge is not reached);
   `ridge_times_match_a_reference_search` (the 32 two-plate cases of M2-04,
   seeds 1–3 at 10,000 cells and seed 1 at 100,000 cells, the last added
   in review round 1: times, paths and ridges bit for bit
   against `reference_times`; `t = 2ℓ/(−c_s)` bit for bit; ridge cells
   `(0, 0, itself)`; every other reached cell equals `best_offer` of the
   result's own labels, so `t_j ≤ 2(ℓ_i + d_ij R)/(−c_{s_i})` for every
   same-plate pair; 756 unreached cells, all on plates without a divergent
   cell); `one_spreading_rate_gives_the_nearest_ridge` (every divergent
   cell at 4 cm/yr: `2d/(4 cm/yr)` of `reference_nearest` within 1e-12
   relative on all 32 cases).
2. Earth's area–age distribution: `earth_area_age_by_hand` (octahedron,
   `G` = 3/8, 5/8, 7/8, ridge +0, continental and unreached 180 Myr);
   `ages_follow_the_area_age_distribution` (seeds 1–3 at 10,000 cells,
   seed 1 at 100,000, the 32 two-plate cases: continental and unreached
   cells exactly 180 Myr, ridge cells +0, ages in the order of the times,
   the ridge cells an atom at +0 holding 3.7–5.0 % of the reached floor at
   10,000 cells, 2.0 % at 100,000 and 0.4–0.8 % at 650,000 cells (seeds
   1–10, scratch tool), from which on `F(age)` lies within each cell's
   share of the area, the atom's area included; mean 59.88–60.00 Myr of
   the reached floor, asserted 60 ± 0.5 Myr on the 34 worlds without
   ridgeless plates); `crust_ages_are_the_ranked_ridge_times`
   (seeds 1–10, bit for bit); `ages_grow_away_from_a_ridge` and
   `ages_cross_continents_on_their_plate` (the M2-04 setups, now against
   the ranked reference times; ages rise strictly along the reference
   paths; the floor beyond a continental band is reached through it).
3. Seams (`seams_almost_vanish`, 100,000 cells, seed 1): 82 pairs of
   neighbouring oceanic cells of one plate more than 20 Myr apart (1,955
   with M2-04), 20 of them next to a ridge cell, the largest step 106 Myr.
   At 650,000 cells (seeds 1–10, scratch tool) there is no such pair in six
   worlds and 3–110 in the others (seed 1: 0, against 2,034 with M2-04),
   with largest steps of 25–75 Myr. The remaining steps lie next to ridge
   cells so slow that their own front is the first to reach their
   neighbours (clusters near a relative pole); found in review round 1 and
   left to the follow-up M2-13.
4. Margins: `margin_profile_by_hand` (15, 20, 25, 30, 35 km at 0, 150,
   300, 450, 600 km; exactly 35 km from 600 km and for `+∞`);
   `margin_distances_match_a_reference_search` (5 crust layouts on each of
   the 32 two-plate cases and seeds 1–3: 131 comparisons bit for bit);
   `convergent_fronts_make_no_margin` (one plate continental against an
   oceanic one: cells more than 600 km from every non-convergent boundary
   cell are at least 600 km from every margin, convergent fronts among
   them; it fails when convergent cells are made sources);
   `small_ocean_pockets_make_no_margin` (a cap of 0.30 % of the surface
   thins nothing, one of 1.0 % does, its coast starting at most half a
   centre distance away); `two_continents_thin_towards_their_coasts_and_rift`
   (continents on both sides of the boundary with a coast inside each
   plate, a rift and a collision front: every thickness
   `m(d_m) + 35 km · min(1, c/5 cm/yr) · f(d)` within 1e-9, uplift
   unchanged, thinnest crust below 17 km); the M2-04 collision and arc
   tests with the margin term. Real worlds
   (`thickness_and_uplift_stay_in_range_on_real_worlds`, seeds 1–20 at
   10,000 cells and seed 1 at 100,000): oceanic 7 km, continental within
   [15 km, 70 km] (thinnest 15.4 km, thickest 70.0 km), arcs 50.0 km
   before the margin's thinning, exactly 35 km away from fronts and
   margins.
5. Land fraction with Earth's parameters and M2-05's formulas
   (`land_fraction_of_earth_like_worlds`, 10,000 cells, seeds 1–10):
   27.9–29.7 % (29.1, 28.5, 29.1, 29.7, 28.6, 28.7, 29.6, 29.6, 28.0,
   27.9), mean 28.89 %, sea level −669 to −934 m. With the same reference
   in a scratch tool: 100,000 cells, seeds 1–5, 28.1–29.9 % (mean 29.37 %);
   650,000 cells, seeds 1–10, 28.7–31.4 % (mean 29.95 %), sea level −680 to
   −885 m. At 650,000 cells the mean ocean depth below sea level is
   3.67–3.81 km (Earth 3.7 km), the mean land elevation 0.96–1.20 km
   (Earth 0.84 km), the mode of land elevation +0.95 to +1.15 km (Earth
   about +0.1 km: Parsons & Sclater's depths are for sediment-free crust,
   finding 5), and 8.6–11.3 % of the surface is submerged continental
   crust (Earth about 11 %).
6. Maps (`kuniumi generate --seed 1 --cells N --out out/M2-11/N`, N =
   10,000, 100,000 and 650,000; land–sea maps from the reference formulas
   by the scratch tool), viewed:
   - `100000/crust_age.png`: young (dark) bands along the ridges and ages
     rising smoothly away from them inside each plate, without the long
     radial seams of M2-04 (no step reaches 180 Myr any more);
     continents uniformly 180 Myr (yellow); the oldest floor in
     the south-west and the south-east, far from ridges, light green
     rather than at the cap. One step remains where a ridge begins near
     13° N 80° W: ridge cell 39194 spreads at 0.096 cm/yr and its
     neighbour is 105.9 Myr old, a pale tongue beside the dark ridge (about
     3.2 km of sea floor in M2-05; review round 1, M2-13).
   - `100000/crust_thickness.png`: every continent framed by a soft band
     of thinner crust a few hundred kilometres wide (purple-blue fading to
     teal inland); no band along the convergent coasts (the west coast of
     the oval continent near 30° W, whose arc runs next to the coast) and
     the collision belt near 160° E–175° W, 15° S; a dark speck inside the
     north-eastern collision zone is an enclosed oceanic pocket (M2-10),
     with no margin around it.
   - `100000/land_sea_reference.png` and `land_sea_classes.png`:
     continents above water inside rings of submerged continental crust
     (shelves, light blue; the class map shows the ring, about 400 km wide,
     on every coast except the active ones); the deep floor darkens away
     from the ridges, which show as pale lines (the meridional ridge near
     75° W); collision plateaus white and brown, arcs tan; 29.6 % land.
   - `650000/crust_age.png` and `land_sea_reference.png`: the same at
     28 km cells, smooth age fields with dark ridges, light ridge lines in
     the bathymetry, shelves round every passive coast, a collision
     plateau near 35° S 130° W; 30.1 % land.
   - `10000/land_sea_reference.png` (the golden case): the same features
     at 2° cells; two light-blue patches inside continents (near 12° N
     2° W and 38° S 80° E) are continental rifts thinned to 18–21 km and
     submerged, like the Red Sea (checked cell by cell); `golden/trench.png`
     lines one or two cells wide along the convergent boundaries.
7. Golden fingerprints (`tests/cli.rs::generate_golden_fingerprints`,
   10,000 cells, seed 1), recorded after the tests above passed and the
   maps were looked at: `crust_thickness` `0x9118_E6AA_2F43_6389`,
   `crust_age` `0x4293_1D66_883A_3DF9`, `trench` `0xA2D4_773C_1D13_3C93`
   (world-breaking); the other nine values unchanged, `uplift_rate`
   included.
8. `mise run ci` passes: fmt, clippy `-D warnings`, the tests in debug and
   release, MSRV 1.85.0, cargo-deny (advisories, bans, licences, sources
   ok). `cargo doc -p kuniumi-sim` gives no warnings. Test-first: against
   stubs (no fronts, all ages 180 Myr, no margins) 18 of the 26 crust
   tests failed.

Timing (release, this machine): `kuniumi generate --seed 1` (650,000
cells) takes 5.0 s and 203 MB (4.8 s and 199 MB after M2-04); at 100,000
cells 1.0 s and 50 MB.

Review:
- Round 1 (CHANGES_REQUESTED): the reviewer re-ran `mise run ci`, the
  stub commit (18 of 26 failing), the land fractions at 100,000 and
  650,000 cells with its own reference, the maps and the timing, and
  checked the local properties of the times at 100,000 and 650,000 cells
  (no violation). R1-1 major, the 100,000-cell world of criterion 1 was
  missing from `ridge_times_match_a_reference_search` (fixed: added, with
  the reference search). R1-2 minor, the ridge cells form an atom at age
  +0 that the distribution check skips without saying so (fixed: stated
  in the test and in evidence 2 with its measured share). R1-3 minor, "the
  jump next to a slow ridge cell is gone" overclaimed: steps of up to
  106 Myr remain next to very slow ridge cells (fixed: the Result, the
  design sentence and the `ridge_times` doc corrected and quantified at
  650,000 cells; follow-up M2-13). R1-4 nit, the margin doc said fronts
  are never thinned (fixed); R1-5 nit, `roles`' parameter renamed
  `time_s` (fixed).
- Round 2 (APPROVE): all five resolved. R2-1 minor, the changelog still
  said "no more seams" (fixed); R2-2 nit, the test module docs still
  called `reference_times` Bellman–Ford and left out the 100,000-cell
  world (fixed).

Follow-ups: M2-12 (draft, ocean floor of plates without a ridge: at
650,000 cells such plates hold 10.3–14.7 % of the ocean floor at 180 Myr
in four of ten worlds); M2-13 (draft, age steps around very slow ridge
cells).

# M2-12 Ocean floor of plates without a ridge

- Status: done
- Milestone: M2
- Depends on: M2-11
- Crates: kuniumi-sim (stages/tectonics; stages/crust docs), kuniumi-cli (tests)
- Design: 地殻から標高へ / 海洋底：年代から水深; プレートテクトニクス / 境界ごとの処理; プレートテクトニクス / プレート速度と分裂
- Decisions: none

## Goal

Decide what age the ocean floor of a plate without a divergent boundary
gets, so that no large uniform patch of the oldest floor remains in the
one-shot planet.

## Scope

In:
- The finding of M2-11: the ages come from fronts that grow from the
  divergent cells of each plate, ranked onto Earth's area–age
  distribution; a plate whose boundaries are all convergent or transform
  has no front, and its ocean floor is 180 Myr old everywhere, a flat
  abyss about 6 km deep with a step at the plate boundary. Measured with
  Earth's parameters (12 plates): at 650,000 cells, seeds 1–10, such plates
  hold 10.3–14.7 % of the ocean floor in four worlds (seeds 3, 6, 9, 10),
  whose mean ocean-floor age becomes 72–78 Myr instead of 60 Myr; at
  10,000 cells up to 20.5 % (seed 4).
- Options to decide with the owner: rank the floor of such plates too, by
  the time from another source (for example their transform boundaries,
  or the boundaries of any type at the plate's mean rate); give them a
  ridge (a divergent boundary where the relative motion is slowest);
  or keep them at 180 Myr.

Out:
- The time-evolving tectonics of M6, in which every ocean floor has a
  ridge of origin.
- The age steps next to very slow ridge cells (M2-13).

## Implementation notes

Refined from `draft`. The options were measured with a throwaway
prototype in the scratchpad (not committed) on the real grids, plates,
continents, boundaries and crust of the current `main` (Earth's
parameters, 12 plates), with the real Crust stage of M2-05
(`surface_heights_m` with trenches at −8 km, `sea_level_m`) and the
metrics of M2-06.

"Ridgeless floor": the oceanic cells of a plate without a divergent
cell, which `ridge_times` does not reach (time `+∞`). Every variant
below only changes their order in the ranking of
`earth_area_age`; the floor reached by a ridge keeps its time-nearest
order.

### Prototype measurements

650,000 cells, seeds 1–10. Ridgeless plates exist in seeds 3 (plate 5:
11.9 % of the ocean floor, 74 % of the plate oceanic, mean speed
4.1 cm/yr, 443 convergent and 712 transform boundary cells), 6 (plate 1:
10.3 %, 76 % oceanic, 3.9 cm/yr), 9 (plates 6 and 7: 6.5 % and 4.8 %, 2.1
and 0.8 cm/yr) and 10 (plates 1 and 5: 11.3 % and 3.4 %; plate 1 is
entirely oceanic). In the other six worlds no variant changes anything.
Columns for seeds 3 / 6 / 9 / 10: share of the ocean floor at exactly
180 Myr, area-weighted mean age of all ocean floor, mean age of the
ridgeless floor and of the rest, land fraction, sea level, and the median
and largest age step between neighbouring oceanic cells across the
boundary of a ridgeless plate.

| Variant | At 180 Myr | Mean (Myr) | Ridgeless floor (Myr) | Other floor (Myr) | Land | Sea level (m) | Step at its boundary, median / max (Myr) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| keep (as now) | 11.9 / 10.3 / 11.2 / 14.7 % | 74.3 / 72.4 / 73.5 / 77.7 | 180 | 60.0 | 29.6 / 29.7 / 31.4 / 29.4 % | −857 / −834 / −819 / −892 | 63–123 / 100–180 |
| **O. oldest, by distance from its boundaries** | **0** | **60.0** | **134–141 (from 111–122)** | **47–51** | **29.0 / 29.1 / 30.8 / 28.5 %** | **−723 / −717 / −689 / −726** | **24–68 / 52–121** |
| V. virtual ridge, slowest neighbour | 0 | 60.0 | 43–97 | 55–62 | 29.0 / 29.1 / 30.8 / 28.6 % | −723 / −716 / −689 / −728 | 32–70 / 123–179 |
| B. fronts from all its boundaries | 0 | 60.0 | 14–46 | 62–67 | 29.0 / 29.1 / 30.9 / 28.6 % | −725 / −716 / −690 / −728 | 65–109 / 147–179 |

The variants:

- *Keep*: the ridgeless floor at `t_max` = 180 Myr. It is one atom of the
  hypsometry: 7.1–8.9 % of the surface in the 100 m bin at −5,350 m (the
  ocean mode), against 2.9 % in the mode of the worlds without such
  plates, and the sea level is 100–200 m lower than in those worlds.
- *O, oldest*: the ridgeless floor ranked after all floor that a ridge
  reaches, ordered by its distance (m) from its own plate's boundary
  (Dijkstra through the cells of the plate from its boundary cells, each
  starting at half the centre distance to its nearest neighbour on
  another plate, the Voronoi edge, as the margins of M2-11). So it takes
  the oldest part of Earth's distribution, from 111–122 Myr at its
  boundary to nearly 180 Myr in its interior, and the rest of the floor is
  ranked onto the younger part (its mean falls to 47–51 Myr in these four
  worlds). Ordering by time at the plate's mean speed instead of by
  distance gives the same in seeds 3 and 6 and within 1 Myr of the step
  statistics in 9 and 10. The trench rule of M2-04 (the older floor
  sinks at an oceanic convergent boundary) is unchanged and consistent
  with it: the ridgeless floor has time `+∞` there and is the oldest.
- *V, virtual ridge*: the boundary cells of the ridgeless plate facing
  the neighbour with the smallest mean relative speed `|Δv|` become ridge
  cells for the ages only, spreading at that mean speed; their times are
  ranked with all others. The plate then looks like one with a ridge
  along one side. But the boundary there stays transform or convergent
  (no divergence in the kinematics), the trench rule would have to use
  these times (changing `trench`), and at 10,000 cells the slow virtual
  ridges give steep ramps: 237 pairs of neighbouring cells of a
  ridgeless plate more than 20 Myr apart in seed 4 (none with O); with
  the per-cell `|Δv|` instead of the mean, 227.
- *B, all boundaries*: fronts from every boundary cell of the ridgeless
  plate at the plate's area-weighted mean speed, ranked with all others
  (the ticket's example). The floor is young along trenches and
  transforms alike, a ridge-like swell around the whole plate; fronts
  from its transform cells only give much the same (ridgeless floor
  14–47 Myr on average).

Common to O, V and B: the mean age is Earth's 60 Myr in every world, the
cap atom is gone, and the sea level of the ten worlds spans −689 to
−728 m (O: −726 m) instead of −689 to −892 m. Land fraction at 650,000
cells, seeds 1–10: keep 28.7–31.4 % (mean 29.96 %); O 28.5–30.8 % (mean 29.69 %).
At 10,000 cells, seeds 1–10 (ridgeless plates in seeds 1, 3, 4 and 9:
3.3, 6.8, 20.6 and 0.1 % of the ocean floor), with the real Crust stage:
keep 28.1–30.0 % (mean 29.14 %), O 28.1–30.0 % (mean 28.89 %; seed 4
from 30.0 to 28.2 %).

Maps viewed (equirectangular, 2,048 px, ages on a fixed 0–180 Myr
viridis scale):
- seed 3, keep: the ridgeless plate, which wraps the south-west of the
  central oval continent (about 20° S–40° N, 20° W–75° E), is uniformly
  yellow, indistinguishable from the continents; O: a
  smooth dome, green (about 120 Myr) along its boundary to pale
  yellow-green in its interior, the neighbouring floor across its
  boundary teal to green, so the boundary shows as a soft change of tone;
  V: young (dark) along its eastern boundary, ageing towards the
  north-west like a plate with a ridge; B: dark along every boundary,
  purple inside, a strong contrast with the teal neighbours. In the
  elevation map of O the plate is the darkest, deepest basin (about 6 km)
  with no visible step at its boundary.
- seed 10, keep: two yellow patches (the oceanic plate across the
  antimeridian at about 10° S–55° N, and a crescent at 35–70° S,
  150–60° W); O: both domes, teal at their edges and light yellow-green
  in their interiors; V: the antimeridian plate young on its eastern side
  (near 160° W), the crescent young at its eastern tip.

### Plan (recommended option, O)

To be confirmed by the owner (Question below). No change to the plates,
continents, boundaries, thickness, uplift or trenches; `crust_age`, and
with it `surface_height` and the sea level, change in worlds with a
ridgeless plate (world-breaking: the golden world, 10,000 cells, seed 1,
has one, 3.3 % of its ocean floor).

In `stages::tectonics::crust`:

1. `ridgeless_distances(grid, plates, bounds, radius_m) -> Vec<f64>` (m):
   for every cell of a plate with no divergent cell, the length of the
   shortest path through cells of its plate from the plate's boundary:
   sources are its boundary cells (`bounds.kind` not `None`), each
   starting at `½ d_ij R` for its nearest neighbour `j` on another plate
   (the plate boundary lies on the Voronoi edges); multi-source Dijkstra
   with the private `grow` of M2-11, passable: every cell of the plate.
   `+∞` for the cells of plates with a divergent cell. Paths through
   continental cells too, as `ridge_times`.
2. `earth_area_age(grid, continental, time_s, ridgeless_m) -> Vec<f64>`
   (s): the oceanic cells in one order, first those with a finite time
   sorted by (time, index), then those with an infinite time and a finite
   ridgeless distance sorted by (distance, index); with `A` their area
   (summed in that order) the age is `t_max (1 − √(1 − G_i))`,
   `G_i = (A_<i + a_i/2) / A`, as in M2-11. Ridge cells (time 0) +0;
   continental cells, and oceanic cells with neither (none, unless a
   spreading rate underflows), `t_max`.
3. `crust` passes the distances; the trench rule keeps the ridge times
   (`+∞` on ridgeless floor, which is the oldest). Its signature and
   `Crust` stay.
4. Design (地殻から標高へ / 海洋底): the sentence "海嶺のないプレートの海洋底
   は180 Myrです" replaced by the new rule, with the owner's agreement.

Tests (test-first, `crust/age_tests.rs`):

- `ridgeless_distances` against the shared Bellman–Ford reference with
  start distances of M2-11 (`reference_from`) on the 10,000-cell worlds
  with ridgeless plates (seeds 1, 3, 4, 9) and on a two-plate case whose
  boundary is made all convergent: bit for bit; `+∞` exactly on the plates
  with a divergent cell.
- `earth_area_age` by hand on the octahedron with a ridgeless cell: the
  ridgeless cells take the last shares.
- On the real worlds (seeds 1–10 at 10,000 cells, seed 1 at 100,000) and
  the 32 two-plate cases: no oceanic cell at `t_max`; ridgeless floor
  older than all other floor and in the order of its distances; the
  area-weighted cumulative distribution within each cell's share of `F`
  as in M2-11 (ridge atom included); mean 60 ± 0.5 Myr in every world
  (M2-11 asserted it only on worlds without ridgeless plates).
- The existing land-fraction test (`margin_tests.rs`, M2-05's reference
  formulas): every seed in 20–40 %, mean of the ten in 28–31 %.
- Golden fingerprints of `crust_age` and `surface_height`, the golden sea
  level and the reference metrics of `metrics.txt`
  (`tools/hypsometry_reference.py`, as M2-06 and M2-10) re-recorded after
  these pass and the maps are viewed; `trench` and the other fields
  unchanged.

### Question for the owner

What age should the ocean floor of a plate without a divergent boundary
get?

- Recommended: **O, the oldest part of Earth's distribution**, ordered
  by distance from the plate's own boundary, as planned above: the floor
  whose ridge is gone is older than all floor with an active ridge; the
  trench rule stays consistent; only `crust_age` changes.
- V: a virtual ridge on the boundary with the slowest relative motion,
  ranked with all other floor: looks like a plate with a ridge, but the
  ridge is not in the boundary types, `trench` changes, and slow virtual
  ridges ramp steeply at 10,000 cells.
- B: fronts from all the plate's boundaries at its mean speed, ranked
  with all other floor: young floor along its trenches too.
- Keep 180 Myr and record the limitation in the design document until M6.

### Decision (owner, 2026-11-12)

The owner chose the recommended option, **O**: the ridgeless floor takes
the oldest part of Earth's area–age distribution, older the farther it
lies from its plate's boundary (about 111–122 Myr at the boundary, close
to 180 Myr inside); only `crust_age` changes and the trench rule stays
consistent; world-breaking. The choice approves the plan above
(`ridgeless_distances`, the extended ranking of `earth_area_age`, the
update of the design document's sentence on the ocean floor).

## Acceptance criteria

- [x] With Earth's parameters, seeds 1–10 at 650,000 cells (scratch tool
      on the committed code): no ocean floor at exactly 180 Myr (0 %) and
      the mean age 60 ± 0.5 Myr in every world; the ridgeless share,
      land fraction and sea level of each world reported.
- [x] `ridgeless_distances` matches the reference search bit for bit on
      the worlds and the two-plate case above; `+∞` exactly on plates
      with a divergent cell.
- [x] `earth_area_age` by hand with ridgeless cells; on the real worlds
      and the 32 two-plate cases: no oceanic cell at `t_max`, ridgeless
      floor older than all other floor and in the order of its distances,
      the distribution within each cell's share of `F`, mean 60 ± 0.5 Myr.
- [x] Land fraction at 10,000 cells, seeds 1–10 (M2-11's test): every seed
      in 20–40 %, mean in 28–31 %; reported also at 650,000 cells.
- [x] Maps of `crust_age` and `surface_height` of seeds 3, 6, 9 and 10 at
      650,000 cells and of the golden world viewed and described.
- [x] Golden fingerprints of `crust_age` and `surface_height`, the golden
      sea level and the reference metrics re-recorded after the tests
      above pass and the maps are viewed; `trench` unchanged;
      world-breaking line and `CHANGELOG.md` entry; the design document
      updated.
- [x] `mise run ci` passes.

## Verification

```sh
cargo run --release -p kuniumi-cli -- generate --seed 3 --out out/M2-12/seed3
cargo run --release -p kuniumi-cli -- generate --seed 1 --cells 10000 --out out/M2-12/10k
```

## Result

What was done:
- `kuniumi-sim`, `stages::tectonics::crust`: `ridgeless_distances` (the
  distance of every cell of a plate without a divergent cell from its
  plate's boundary, Dijkstra through the plate from its boundary cells,
  each starting half a centre distance from its nearest neighbour on
  another plate, with the private `grow` of M2-11; `+∞` on plates with a
  divergent cell), re-exported from `tectonics`; `earth_area_age` takes
  these distances as a fourth argument and ranks the ridgeless floor after
  all floor that a ridge reaches, by (distance, index); `crust` passes them
  and keeps the trench rule on the ridge times. Documentation of the
  module, `crust`, `Crust::age_s` and `ridge_times`; the sea level and
  trench figures of `stages::crust` brought up to date.
- Tests (`crust/age_tests.rs`): `ridgeless_distances_on_the_octahedron`,
  `ridgeless_distances_match_a_reference_search`, the extended
  `earth_area_age_by_hand`, `ages_follow_the_area_age_distribution` and
  `crust_ages_are_the_ranked_ridge_times`; the two callers of
  `earth_area_age` in `crust/tests.rs` pass `+∞` (both plates have a
  ridge there). `kuniumi-cli` tests: the golden values and the reference
  metrics.
- `docs/design-ja.md` (地殻から標高へ / 海洋底: the rule for plates without
  a ridge), `CHANGELOG.md` (Changed, World-breaking).

Deviations from the plan:
- The reference search of `ridgeless_distances` runs on seeds 1–10 at
  10,000 cells (all of them, not only those with ridgeless plates) and on
  64 two-plate cases (the 32 of M2-04 with the divergent cells of plate 1,
  or of both plates, made convergent). Seed 8 at 10,000 cells also has a
  plate without a divergent cell, but it is all continental, so its ages
  do not change.
- The trench rule is checked directly as well
  (`crust_ages_are_the_ranked_ridge_times` compares `trench` with the
  roles on the ridge times alone).

Evidence:

1. 650,000 cells, Earth's parameters, seeds 1–10, measured with the
   scratch tool on the committed code (`crust` and the real Crust stage):
   no oceanic cell at exactly 180 Myr in any world (0 %, against
   10.3–14.7 % of the ocean floor in seeds 3, 6, 9 and 10 before), mean
   age 59.997–59.999 Myr in all ten (72.4–77.7 Myr in those four before).
   Per seed, share of the ocean floor on plates without a ridge, land
   fraction and sea level: 1: 0 %, 30.1 %, −699 m; 2: 0 %, 30.8 %, −691 m;
   3: 11.9 %, 29.0 %, −723 m; 4: 0 %, 28.7 %, −720 m; 5: 0 %, 30.0 %,
   −701 m; 6: 10.3 %, 29.1 %, −717 m; 7: 0 %, 29.2 %, −716 m; 8: 0 %,
   30.7 %, −689 m; 9: 11.2 %, 30.8 %, −689 m; 10: 14.7 %, 28.5 %, −726 m.
   The ridgeless floor is 118.0, 122.1, 119.6 and 110.9 Myr old at its
   youngest (seeds 3, 6, 9, 10), 134–141 Myr on average; the rest of the
   floor 47–51 Myr on average. The committed ages equal those of the
   prototype's variant O (written separately in the scratchpad) bit for
   bit in all ten worlds. Ocean mode: 2.9 % of the surface in every world
   (7.1–8.9 % at −5,350 m before).
2. `ridgeless_distances_on_the_octahedron` (by hand: boundary cells at
   `q/2`, the interior cell at `3q/2`, a plate with a divergent cell
   `+∞`, a plate of one transform cell at `q/2`);
   `ridgeless_distances_match_a_reference_search` (seeds 1–10 at 10,000
   cells and the 64 two-plate cases without ridges: bit for bit against
   `reference_from` with starts computed in the test; `+∞` exactly on the
   plates with a divergent cell; ridgeless cells in seeds 1, 3, 4, 8 and 9
   and in every two-plate case).
3. `earth_area_age_by_hand` (octahedron: the ridgeless cells take the
   shares 5/10, 7/10, 9/10 after the floor reached by a ridge, also when
   the latter's time is larger than every distance; a cell with a finite
   time keeps its share when it has a finite distance too, added in
   review round 1; a cell with neither key `t_max`); `ages_follow_the_area_age_distribution` (seeds 1–10 at
   10,000 cells, seed 1 at 100,000, the 32 two-plate cases and the 64
   without ridges): continental cells exactly 180 Myr, no oceanic cell at
   180 Myr, ridge cells +0, ages in the order of (time, then ridgeless
   distance), `F(age)` within each cell's share of the area, mean
   59.80–60.00 Myr, asserted 60 ± 0.5 Myr on every world (ridgeless floor
   3.3, 6.8, 20.6 and 0.1 % of the ocean floor in seeds 1, 3, 4 and 9,
   youngest at 147.1, 133.2, 98.3 and 175.1 Myr);
   `crust_ages_are_the_ranked_ridge_times` (seeds 1–10: `crust`'s ages
   bit for bit, trenches those of the ridge times alone).
4. Land fraction at 10,000 cells, seeds 1–10: M2-11's test with the
   reference formulas of M2-05 without deepened trenches
   (`margin_tests::land_fraction_of_earth_like_worlds`) 27.8–29.6 %, mean
   28.58 % (28.89 % before); the Crust stage's own test
   (`stages::crust::tests::land_fraction_of_earth_like_worlds`)
   28.14–30.04 %, mean 28.90 % (29.14 % before). At 650,000 cells
   28.5–30.8 %, mean 29.69 % (29.96 % before).
5. Maps (`kuniumi generate --seed S --out out/M2-12/seedS`, S = 3, 6, 9,
   10; `--seed 1 --cells 10000 --out out/M2-12/10k`), viewed:
   - `seed3/crust_age.png`: the ridgeless plate round the south-west of the
     central oval continent is a smooth dome, green (about 120 Myr) along
     its boundary to pale yellow-green inside, no longer the continents'
     yellow; the rest as before but younger overall. `surface_height.png`:
     that plate is the deepest, darkest basin, with only a change of tone
     at its boundary; `hypsometry.png`: the ocean side falls smoothly from
     its peak near −5.2 km, with no spike (7.6 % of the surface in one bin
     at −5,350 m before).
   - `seed6/crust_age.png`: the plate at about 25° N–45° S, 100–170° E is
     a dome from green at its edges to yellow along a narrow streak, the
     line farthest from its boundaries (its medial axis); in
     `surface_height` it is uniform deep floor, the streak invisible
     (170–180 Myr floor lies at 6.19–6.22 km).
   - `seed9/crust_age.png`: two domes, one in the north-west (a pale
     streak along its long axis near 65° N) and one north of the western
     continent (about 20–50° N, 95–140° W), both green to yellow-green.
   - `seed10/crust_age.png` and `surface_height.png`: the oceanic plate
     across the antimeridian (about 10° S–55° N) and the crescent at
     35–70° S, 150–60° W are domes, teal at their edges and light
     yellow-green inside; in the elevation map both are the deepest
     basins, without steps at their boundaries.
   - `10k/crust_age.png` (the golden world): the ridgeless plate is the arc
     of floor north of the central continent (about 7–57° N, 40° W–40° E),
     147–178 Myr, light green against the blue-green floor around;
     `surface_height.png`: a deep basin there, lined on its north by
     trenches (it sinks, as the oldest floor).
6. Golden values (`tests/cli.rs`, 10,000 cells, seed 1), recorded after the
   tests above passed and the maps were looked at: `crust_age`
   `0x6228_96F5_9A74_CB65`, `surface_height` `0x506C_FB1B_6C86_7FD6`, sea
   level −747.142 141 911 436 2 m (was −782.657 m), the reference metrics
   from `tools/hypsometry_reference.py out/M2-12/10k` (land fraction
   0.291 092, was 0.293 999); the other eleven fingerprints unchanged,
   `trench` `0xA2D4_773C_1D13_3C93` among them. World-breaking line and
   `CHANGELOG.md` entry; design document updated.
7. `mise run ci` passes: fmt, clippy `-D warnings`, rustdoc `-D
   warnings`, the tests in debug and release, MSRV 1.85.0, cargo-deny
   (advisories, bans, licences, sources ok). Test-first: against stubs
   (`ridgeless_distances` all `+∞`, the new argument of `earth_area_age`
   ignored) 4 of the 9 age tests failed (commit "tests of the ridgeless
   floor against stubs").

Timing (release, this machine): `kuniumi generate` at 650,000 cells takes
5.3–5.5 s and 204–209 MB (seed 1: 5.4–5.5 s, 209 MB; seed 3: 5.3 s,
204 MB); M2-11 measured 5.0 s and 203 MB for seed 1 before M2-05 and
M2-10.

Review:
- Round 1 (APPROVE): the reviewer re-ran `mise run ci`, the stub commit
  (4 of 9 age tests failing), the 650,000-cell statistics of seeds 1–10
  with its own script (0 cells at the cap, means 59.997–59.999 Myr, the
  per-seed values of evidence 1), the "before" figures on a build of
  `main`, the reference metrics, the maps and the timing, and compared
  the world files of seeds 3 and 10 with `main` (only `crust_age` and
  `surface_height` differ). R1-1 minor: no test gave a cell both a finite
  time and a finite ridgeless distance, so ranking such a cell twice went
  unnoticed (fixed: a case in `earth_area_age_by_hand`, which fails under
  that mutation). R1-2 nit: the `Crates:` line left out the doc changes of
  `stages/crust` (fixed).

Follow-ups: none.

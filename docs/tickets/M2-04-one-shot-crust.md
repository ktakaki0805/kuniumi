# M2-04 One-shot crust thickness, age and uplift

- Status: done
- Milestone: M2
- Depends on: M2-03
- Crates: kuniumi-sim (stages/tectonics), kuniumi-cli (tests)
- Design: プレートテクトニクス / 境界ごとの処理; 地殻から標高へ
- Decisions: none

## Goal

Without time evolution, give each cell a plausible crust thickness, ocean-floor
age and uplift rate from its distance to the boundaries, so later stages have
realistic input until M6 replaces this.

## Scope

In:
- Distances (m) from each boundary type by multi-source Dijkstra on the cell
  graph, with centre distances × radius as edge weights (deterministic heap,
  ties by index).
- Ocean-floor age = distance to the nearest divergent boundary ÷ that
  boundary's half spreading rate, capped at 180 Myr (stored in seconds).
- Continental thickness: 35 km base; continent–continent convergence thickens
  the crust up to about 70 km over a few hundred kilometres (plateau);
  ocean–continent convergence builds an arc (about 45–50 km) inland of the
  trench. Profiles are documented with their parameters in SI.
- A trench flag on the oceanic side of subduction boundaries.
- `uplift_rate` (m/s): non-zero near convergent boundaries, scaled by the
  convergence rate; it feeds the erosion stage until M6.
- Outputs: `crust_thickness` (m), `crust_age` (s), `uplift_rate` (m/s),
  `trench` (class).

Out:
- Conversion to elevation (M2-05).

## Implementation notes

Gaps found while refining (the ticket was `ready`): the ticket leaves open
which cells a distance reaches, the oceanic thickness, the age of
continental crust, how the crust across a boundary and the subducting side
are found, the exact profiles and the uplift factor; and a prototype shows
that the ages of this method put the sea level of M2-05 far below the
continents (below). Measured with a throwaway prototype in the scratchpad
on the real grids, plates, continents and boundaries of M2-02 and M2-03
(Earth's parameters, 12 plates).

### Distances (`nearest_sources`)

Multi-source Dijkstra as in the partition: every source enters a
`BinaryHeap` at distance 0; entries are ordered by distance
(`f64::total_cmp`), then cell, then source (the `Entry` of the partition);
the first entry popped for a cell settles it. A settled cell pushes each
neighbour `j` that is unsettled, **on the same plate** and passable, at
`d + centre_distance_rad · radius_m`. So each cell gets the distance (m) to
the nearest source of its own plate, along a path inside the plate, and
that source; cells that no source reaches get none. Distances never cross
a plate boundary: the age and the orogens of a plate come from its own
boundaries.

```rust
pub struct Nearest { pub distance_m: Vec<f64>, pub source: Vec<Option<u32>> }
pub fn nearest_sources(grid: &SphereGrid, plates: &[u32], radius_m: f64,
                       sources: &[u32], passable: &dyn Fn(u32) -> bool) -> Nearest;
```

### Roles of convergent cells

For a convergent cell `i` (plate `a`, plate across `b` from
`Boundaries::other_plate`), the *cell across* is its neighbour on `b` with
the longest shared Voronoi edge (the first in neighbour order on a tie).

| Cell `i` | Cell across | Role of `i` |
| --- | --- | --- |
| continental | continental | collision front |
| continental | oceanic | arc front (Andean, overriding) |
| oceanic | continental | trench (subducting) |
| oceanic | oceanic, younger than `i` | trench (the older, denser floor sinks; design 境界ごとの処理) |
| oceanic | oceanic, not younger | island-arc side: no thickening, no uplift |

The ages compared are the uncapped `d / h` of both cells (`+∞` on a plate
without a divergent boundary); on an exact tie the lower plate id subducts.
Island arcs are left to M6-03: the criterion below fixes oceanic crust at
6–8 km.

### Fields

All four on every cell, added after `convergence_rate` in the order of
the Scope (thickness, age, uplift, trench):

- `crust_age` (s): for oceanic cells `min(d / h, 180 Myr)`, with `d` the
  distance to the nearest divergent cell of the plate (any crust; paths
  through any cell of the plate) and `h = −c / 2` that cell's half spreading
  rate (`c < 0` on divergent cells, so `h > 0`); 180 Myr on a plate without
  a divergent cell. Continental cells: 180 Myr, "older than any ocean
  floor" (their real age needs M6).
- `crust_thickness` (m): oceanic cells 7 km, the mean thickness of oceanic
  crust (White, McKenzie & O'Nions 1992, JGR 97(B13), 19683: 7.1 ± 0.8 km).
  Continental cells `35 km + max(ΔT_coll f_coll(d_coll), ΔT_arc f_arc(d_arc))`,
  with `d_coll`, `d_arc` the distances to the nearest collision and arc
  fronts of the plate along paths through continental cells only:
  - collision: `ΔT_coll = 35 km · min(1, c / 5 cm/yr)` (up to 70 km; India
    and Eurasia converge at about 4–5 cm/yr and built Tibet's 70 km crust),
    `f_coll = 1` up to 200 km from the front, then a cosine taper
    `½(1 + cos(π (d − 200 km) / 200 km))` to 0 at 400 km, on both plates;
  - arc: `ΔT_arc = 15 km` (to 50 km), `f_arc = ½(1 + cos(π (d − 200 km) / 200 km))`
    for `d ≤ 400 km` (0 at the front and at 400 km, 1 at 200 km inland,
    where arc–trench gaps of 100–300 km put Andean arcs), not rate-scaled.
  `c` is the convergence rate of the front cell. The maximum, not the sum,
  so overlapping belts stay at most 70 km.
- `uplift_rate` (m/s): continental cells
  `min(10 mm/yr, max(0.1 c_coll f_coll(d_coll), 0.03 c_arc f_arc(d_arc)))`,
  0 on oceanic cells and beyond the belts. The ratios are order-of-
  magnitude Earth values (Himalaya: rock uplift of several mm/yr at about
  4–5 cm/yr of convergence; central Andes: about 1–2 mm/yr at about
  7 cm/yr); M3 calibrates erodibility against them. The cap of 10 mm/yr is
  the criterion's bound: closing rates reach 19 cm/yr in the prototype.
- `trench` (class): 1 on trench cells, 0 elsewhere.

### Constants and API

```rust
pub const CRUST_AGE: &str = "crust_age";               // s
pub const CRUST_THICKNESS: &str = "crust_thickness";   // m
pub const UPLIFT_RATE: &str = "uplift_rate";           // m/s
pub const TRENCH: &str = "trench";                     // class: 0, 1
pub const MAX_CRUST_AGE_S: f64;                // 180 Myr
pub const OCEANIC_THICKNESS_M: f64;            // 7 km
pub const CONTINENTAL_THICKNESS_M: f64;        // 35 km
pub const COLLISION_THICKNESS_M: f64;          // 70 km
pub const COLLISION_REFERENCE_RATE_M_S: f64;   // 5 cm/yr
pub const COLLISION_PLATEAU_M: f64;            // 200 km
pub const COLLISION_TAPER_M: f64;              // 200 km
pub const ARC_THICKNESS_M: f64;                // 50 km
pub const ARC_DISTANCE_M: f64;                 // 200 km
pub const ARC_HALF_WIDTH_M: f64;               // 200 km
pub const COLLISION_UPLIFT_RATIO: f64;         // 0.1
pub const ARC_UPLIFT_RATIO: f64;               // 0.03
pub const MAX_UPLIFT_RATE_M_S: f64;            // 10 mm/yr
pub struct Crust { pub age_s: Vec<f64>, pub thickness_m: Vec<f64>,
                   pub uplift_rate_m_s: Vec<f64>, pub trench: Vec<bool> }
pub fn crust(grid: &SphereGrid, plates: &[u32], continental: &[bool],
             bounds: &Boundaries, radius_m: f64) -> Crust;
```

In a new submodule `stages::tectonics::crust`, re-exported from
`tectonics`; `run` adds the four fields. Cosines through `det::cos`.

### Prototype measurements

At 100,000 cells (seeds 1–5) and 10,000 cells (seeds 1–10): oceanic
crust 7 km; continental 35–70 km (collision belts reach 70 km in 10 of
the 15 worlds); largest uplift per world 2.6–10 mm/yr (the cap is reached once, seed 2
at 100,000 cells, where a front converges at 19 cm/yr); convergent cells per
role at 100,000 cells, seed 1: collision 176, arc 310, island-arc side
294, trench 601. Maps viewed (`out/M2-04/proto100k`): thin arc belts
inland of active margins, a few plateaus of 60–70 km, ages rising from
dark ridges, with radial streaks where neighbouring ridge cells have
different rates.

**Finding: the ocean floor is old.** The area-weighted mean age of the
ocean floor is 79–144 Myr (Earth: about 60 Myr, the triangular area–age
distribution of Parsons 1982), and 23–68 % of the ocean floor is at the
180 Myr cap at 100,000 cells: half spreading rates are 1.1–3.4 cm/yr
(median per seed) and many ocean cells lie thousands of kilometres from
their nearest ridge, often from a slow part of it near the relative pole.
Taking the time-nearest ridge instead (growth from each ridge cell at its
own half rate) or the plate's mean rate changes the mean age by less
than 10 %. With M2-05's formulas (Parsons & Sclater depth, Airy isostasy,
1.335e18 m³ of water), the prototype puts the sea level 0.7–1.5 km below
the datum: every continental cell stays dry and the land fraction is
40.0 % for every seed, at the upper bound of M2-05's criterion (20–40 %).
Even Earth's area–age distribution would give a mean depth of 4.9 km
against the 4.4 km that the water fills on 60 % of the surface, so the
land fraction of 29 % needs the submerged continental margins of the
design (地殻から標高へ / 大陸: 分裂で薄くなった縁辺部（25〜30 km）は海面下に沈み、
大陸棚になります). Thinning the passive margins in this ticket (25 km at a
coast inside a plate, back to 35 km over 200 km) gives only 36.8–40.0 %
with today's ages, because the sea level is so low.

### Question for the owner

Recommended: the plan above as it stands (the age method of the ticket,
no passive margins in this ticket), and a new draft ticket M2-11 "Ocean
floor age, passive margins and sea level of the one-shot planet", to be
done before M2-05, that decides how to reach Earth's land fraction
(passive margins, spreading rates, or both). Alternatives: the same plan
with the passive-margin thinning above added to this ticket; or the same
plan with collision thickness not scaled by the convergence rate (70 km
at every collision).

### Decision (owner, 2026-11-11)

The owner approved the plan above as written and kept the age method of the
ticket. The continental margins are not thinned in this ticket; the draft
ticket M2-11 (ocean-floor age, passive margins and sea level of the
one-shot planet), placed before M2-05, decides how the land fraction is
reached.

## Acceptance criteria

- [x] On synthetic two-plate divergent setups (all oceanic; two hemispheres
      split by four great circles, relative pole on the circle, 10,000
      cells with grid seeds 1–3 and 100,000 cells with seed 1), every
      cell's age equals `d / h` within 1e-12 relative, with `d` and the
      ridge cell from an independent reference search written in the test
      (or exactly 180 Myr where `d / h` exceeds it); age is 0 on ridge
      cells and increases strictly from each cell to the next along the
      reference shortest paths.
- [x] Thickness ranges for seeds 1–20 at 10,000 cells and seed 1 at
      100,000 cells: oceanic cells within 6–8 km (exactly 7 km),
      continental cells within 25–75 km (35–70 km by construction), with
      collision belts above 60 km and arcs of 50 km present among
      these worlds.
- [x] Synthetic convergent setups match the profiles worked out by hand:
      collision (both continental), arc and trench (continental against
      oceanic), and ocean–ocean with the older side flagged as trench.
- [x] Uplift rates between 0 and 10 mm/yr (in m/s) on the same worlds and on
      a synthetic collision converging at 20 cm/yr (capped); 0 on oceanic
      cells and beyond 400 km of every front.
- [x] Images of the four fields (100,000 cells, seed 1, and the 10,000-cell
      golden case) viewed and described.
- [x] Golden fingerprints for 10,000 cells, seed 1, recorded after the
      tests above pass and the images are viewed.
- [x] `mise run ci` passes.

## Verification

```sh
cargo run --release -p kuniumi-cli -- generate --seed 1 --cells 100000 --out out/M2-04/100k
cargo run --release -p kuniumi-cli -- generate --seed 1 --cells 10000 --out out/M2-04/10k
```

## Result

What was done:
- `kuniumi-sim`, `stages::tectonics`: the submodule `crust`
  (`nearest_sources`, `Nearest`, `crust`, `Crust`, the four field names and
  the constants of the plan), re-exported from `tectonics`, with its tests
  (`crust/tests.rs`); `run` adds `crust_thickness`, `crust_age`,
  `uplift_rate` and `trench` after `convergence_rate`; the module, `run`
  and `generate` documentation name them. The roles of convergent cells
  are a private `Role` (collision, arc, island-arc side, trench).
- `kuniumi-cli`: tests only (the twelve fields written, the ranges of the
  new fields on the 1,000-cell world, the golden values).
- `docs/tickets/M2-11-ocean-age-margins-sea-level.md` (draft, listed before
  M2-05 in the README), as the owner decided; `CHANGELOG.md` (Added).

Deviations from the plan:
- The age is computed as `2d / (−c)`, equal to `d / (−c/2)` within rounding,
  so that no NaN can arise even if a half rate underflowed to zero.
- Larger and smaller of two values by explicit comparisons (`larger`,
  `smaller`) instead of `f64::max`/`min`, so no zero sign depends on the
  platform (determinism rule 2); the raised cosine is clamped at `+0`.

Evidence:

1. Ages (`crust::tests::ages_grow_away_from_a_ridge`): 32 two-plate cases
   (10,000 cells with grid seeds 1–3 and 100,000 cells with seed 1; the
   equator and three tilted circles; both signs of the relative rotation,
   whose pole lies on the circle, plus a common rotation), both plates
   oceanic. For every cell the age equals `d / h` within 1e-12 relative,
   with `d` and the ridge cell from `reference_nearest` (Bellman–Ford
   sweeps written in the test, ties to the lower source), or exactly
   180 Myr where `d / h` is larger (such cells occur); ridge cells have age
   `+0`; along the reference shortest paths the age increases strictly from
   each cell to the next below the cap; `d` is at least the great-circle
   distance to the ridge cell; both plates have ridges; thickness 7 km and
   no uplift everywhere. `ages_cross_continents_on_their_plate` (the same
   cases with a continental band inside plate 0, 1,300–2,600 km from the
   boundary): oceanic ages beyond the band are `d / h` along paths through
   the band, continental ages exactly 180 Myr; it fails when the ridge
   search is restricted to oceanic cells.
   `nearest_sources_matches_a_reference_search`
   (seeds 1–3 at 10,000 cells, every 97th cell a source, through all or
   only continental cells: distances bit for bit, same sources) and
   `nearest_sources_on_the_octahedron` (by hand: 0, `πR/2`, `πR`; an
   impassable cell; a plate of one cell; a tie to the lower source).
2. Thickness ranges (`thickness_and_uplift_stay_in_range_on_real_worlds`,
   Earth's parameters, seeds 1–20 at 10,000 cells and seed 1 at 100,000
   cells): oceanic cells exactly 7 km, continental cells within 35–70 km
   and 180 Myr old; thickest crust 70.0 km (a collision plateau); thickest
   continental crust farther than 400 km from every collision front
   50.0 km (an arc crest).
3. Profiles (`profiles_worked_out_by_hand`: collision 1 at 0–200 km,
   `½ + √2/4` at 250 km, ½ at 300 km, 0 from 400 km; arc 0 at 0, ½ at
   100 km, 1 at 200 km, ½ at 300 km, 0 from 400 km) and the 32 two-plate
   cases: `a_collision_thickens_and_uplifts_both_plates` (both
   continental: every cell `35 km + 35 km · min(1, c/5 cm/yr) · f(d)` and
   `min(10 mm/yr, 0.1 c f(d))` within 1e-9 against the formula evaluated
   with std `cos` and the reference distances; exactly 35 km and no uplift
   beyond 400 km; both the saturated and the rate-scaled thickening and
   the uplift cap occur); `an_arc_rises_inland_of_a_trench` (each plate in
   turn continental: `35 km + 15 km · f_arc(d)` and `0.03 c f_arc(d)` on
   the continental plate, 7 km and no uplift on the oceanic one, trenches
   exactly its convergent cells); `the_older_ocean_floor_subducts` (both
   oceanic: a convergent cell is a trench exactly when its reference age
   exceeds that of the cell across; both outcomes occur);
   `roles_of_hand_chosen_ages` (plate 0 older, plate 1 older, equal ages
   and both `+∞`: the lower plate id subducts on ties; collision and arc
   roles).
4. Uplift: on the real worlds of 2 every rate lies in `[0, 10 mm/yr]` with
   no `−0`, is 0 on oceanic cells and on continental cells farther than
   400 km from every front of their plate; fastest 9.92 mm/yr. In the
   two-plate collision, fronts converge at up to 20 cm/yr and the cap of
   exactly 10 mm/yr is reached. The CLI test checks the same ranges on the
   1,000-cell world.
5. Maps (`kuniumi generate --seed 1 --cells 100000 --out out/M2-04/100k`,
   release, exit code 2 naming `crust`), viewed:
   - `crust_thickness.png`: oceans uniform dark purple (7 km), continents
     uniform teal (35 km); narrow green belts a few hundred kilometres
     inland of the convergent boundaries of `boundary_type.png` where a
     continent faces ocean (the north-west coast of the oval continent
     near 30° W, the north coast of the small continent near 140° W
     10° S, the west coast of the north-eastern continent near 135° E,
     the south polar coasts); a yellow plateau (70 km) where continents
     meet across a convergent boundary, from about 150° E to 175° W at
     13°–17° S (across the antimeridian), and a lighter, rate-scaled one
     (about 60 km) at 15°–30° N, 135°–160° E.
   - `crust_age.png`: dark (young) bands along the divergent boundaries,
     brightening away from them inside each plate, and sharp contrasts at
     plate boundaries; continents and much of the ocean of slow or distant
     ridges at the 180 Myr cap (yellow). Inside plates, seams run away
     from the ridges where cells take their age from ridge cells of very
     different rates (a ridge's rate falls to zero at its relative pole):
     measured in the review at 100,000 cells, seed 1, 1,955 pairs of
     neighbouring oceanic cells of the same plate differ by more than
     20 Myr; the diagonal seam of plate 3 (52° S to 17° S, 38° W to 48° E)
     has 251 pairs above 50 Myr, up to 127 Myr; jumps reach 159 Myr on
     plate 5 (ridge cells spreading at 4.73 and 0.05 cm/yr), and next to a
     ridge cell spreading at 0.10 cm/yr the age jumps from 0 to 180 Myr in
     one step; 127 Myr at 10,000 cells. With M2-05's depths a 180 Myr step
     is a step of about 3.7 km in the sea floor. The method is the
     owner's (Decision above); M2-11 lists the seams as an input.
   - `uplift_rate.png`: zero (dark) everywhere except the arc belts
     (dim, about 1–2 mm/yr) and the collision plateaus (bright, up to
     6.4 mm/yr at about 13° S, 180°).
   - `trench.png`: lines one cell wide along the convergent boundaries of
     `boundary_type.png`, on one side only, continuous over thousands of
     kilometres (from about 45° W 15° S to 55° E 35° N, and from 150° W
     35° S to 75° W 70° S), absent along divergent and transform
     boundaries.
   - `out/M2-04/10k/` (the golden case): the same features at 2° cells:
     arcs two to three cells wide, one 70 km plateau near 12° S 5° E with
     its uplift peak, trench lines with gaps where ocean meets ocean and
     the older side alternates between the plates.
   - `out/M2-04/650k/` (650,000 cells, seed 1, another partition): smooth
     arc belts along several active margins, one collision plateau near
     35° S 130° W and a weak, rate-scaled one near 15° N 25° W.
6. Golden fingerprints (`tests/cli.rs::generate_golden_fingerprints`,
   10,000 cells, seed 1): `crust_thickness` `0xAB89_127D_E598_AF4E`,
   `crust_age` `0x9B68_EB58_6AF8_700A`, `uplift_rate`
   `0x6725_720D_C775_FE38`, `trench` `0x0DC3_352B_B86C_E68A`, recorded
   after the tests above passed and the maps were looked at; the eight
   earlier values unchanged. `generate_writes_the_same_bytes_twice` covers
   the 25 files.
7. `mise run ci` passes: fmt, clippy `-D warnings`, the tests in debug and
   release (366 in debug, 357 in release), MSRV 1.85.0, cargo-deny
   (advisories, bans, licences, sources ok). `cargo doc -p kuniumi-sim`
   gives no warnings. Test-first: against stubs (no distances, zero
   fields) 11 of the first 12 new tests failed (the constants passed).

Timing (release, this machine): `kuniumi generate --seed 1 --cells 650000`
takes 4.8 s and 199 MB (4.2 s and 178 MB before; the peak is not in the
crust, so starting the heaps small, R1-4, left it unchanged); at 100,000
cells `generate` takes 0.98 s and 49 MB.

Review:
- Round 1 (APPROVE): the reviewer re-ran `mise run ci`, regenerated the
  maps at 10,000, 100,000 and 650,000 cells and confirmed the
  descriptions, recomputed the four fingerprints independently, confirmed
  the test-first claim at the stub commit, ran extreme valid parameters
  (radius 1 m to 2e7 m, fractions 0 to 1, 2 and 64 plates, 1,000 cells:
  finite, no −0), and mutated the belts (paths through any cell, sum
  instead of maximum), the cell across, the ages compared, the uplift and
  the continental age (each caught). R1-1 minor, no test covered ridge
  paths through continental cells (fixed:
  `ages_cross_continents_on_their_plate`, which fails under that
  mutation). R1-2 minor, the age seams inside plates were understated as
  "radial streaks" (fixed: their size, measured by the reviewer, is in
  evidence 5 and in M2-11's Scope as an input to its decision; the method
  is the owner's). R1-3 minor, Earth values without sources (fixed:
  Müller et al. 2008 for the oldest ocean floor, DeMets, Gordon & Argus
  2010 for the India–Eurasia rate, Dickinson 1973 for arc–trench gaps,
  Lavé & Avouac 2001 for Himalayan uplift; the arc ratio and the other
  profile values are marked as order-of-magnitude choices made with the
  owner). R1-4 nit, each search reserved a heap of `2n` entries (fixed:
  it starts at the number of sources).

Follow-ups: M2-11 (draft, ocean-floor age, passive margins and sea level).

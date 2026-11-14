# M3-01 Drainage graph: receivers, stack, accumulation

- Status: done
- Milestone: M3
- Depends on: M2-05
- Crates: kuniumi-sim (new `drainage` module)
- Design: 侵食と地形発達 / 1ステップの処理
- Decisions: none

## Goal

A reusable drainage graph on the global grid: where each land cell sends its
water, the processing order from outlets upstream, and accumulated drainage
area and discharge.

## Scope

In:
- `Drainage::build(grid, height_m, sea_level_m, radius_m)`:
  - receiver of each land cell = the neighbour with the steepest positive slope
    `(h_i − h_j) / (d_ij · R)`, ties to the lower index; a cell without a lower
    neighbour is a pit (its own receiver); ocean cells (h ≤ sea level) are
    outlets;
  - donors (the inverse of receivers) in CSR form;
  - stack order: breadth-first from outlets and pits through donors
    (Braun & Willett 2013), deterministic.
- `accumulate(&self, per_cell_input) -> Vec<f64>` in reverse stack order, for
  drainage area (m², input = cell area) and discharge (m³/s, input =
  precipitation rate × area).

Out:
- Resolving pits (M3-02).
- World fields and maps of drainage area and discharge (M3-05), and their
  golden fingerprints in the CLI tests.

## Implementation notes

Gaps found while refining (the ticket was `ready`), settled inside the
ticket's scope: the accessors besides `build` and `accumulate`, the checks
on the inputs, what flats do, the order of the roots and donors in the
stack, how the synthetic surfaces are measured, and how the maps and the
timing are produced.

### Module

`crates/kuniumi-sim/src/drainage.rs`, `pub mod drainage` in `lib.rs`, tests
in `drainage/tests.rs` (so `cargo test -p kuniumi-sim drainage` selects
them). A module of the crate, not a stage: the Erosion stage (M3-03 to
M3-05) and the Hydrology stage (M5-01) both route water on it. No world
field, no random numbers, no I/O.

### API

```rust
#[derive(Clone, Debug, PartialEq)]
pub struct Drainage { /* receivers, donor CSR, stack, base levels, pits */ }

impl Drainage {
    pub fn build(grid: &SphereGrid, height_m: &[f64], sea_level_m: f64,
                 radius_m: f64) -> Self;
    pub fn len(&self) -> usize;               // number of cells
    pub fn is_empty(&self) -> bool;
    pub fn receiver(&self, i: u32) -> u32;    // i itself for outlets and pits
    pub fn receivers(&self) -> &[u32];
    pub fn donors(&self, i: u32) -> &[u32];   // ascending cell index
    pub fn stack(&self) -> &[u32];            // every cell once, after its receiver
    pub fn base_levels(&self) -> &[u32];      // outlets and pits, index order
    pub fn pits(&self) -> &[u32];             // land base levels, index order
    pub fn accumulate(&self, per_cell_input: &[f64]) -> Vec<f64>;
}
```

`build` panics, documented, if `height_m` does not hold one value per cell,
if a height or the sea level is not finite, or if `radius_m` is not finite
and positive; `accumulate` if its input does not hold one value per cell.

### Algorithm

All in index order, integers and IEEE `−`, `×`, `÷` and comparisons only
(determinism rule 2); no floating-point equality (`float_cmp`):
`partial_cmp` decides ties.

1. Receivers. A cell with `h_i ≤ sea_level_m` is an outlet,
   `r(i) = i`. A land cell (`h_i > sea_level_m`) looks at the neighbours
   with `h_j < h_i` strictly and takes the one with the largest
   `s_ij = (h_i − h_j) / (d_ij · R)` (m/m, `d_ij` the centre distance in
   radians from `SphereGrid::centre_distances_rad`, `R = radius_m`); of
   equal slopes the lower cell index. Without such a neighbour the cell is
   a pit, `r(i) = i`. The receiver may be an ocean cell. Flats follow the
   same rule: neighbours of equal height are not lower, so a cell of a flat
   with no lower neighbour is a pit (the uniform +300 m interiors of the
   continents of M2-05 are pits until M3-02 fills them with a gradient).
   Every link goes strictly downhill, so following receivers ends at an
   outlet or a pit within `n − 1` steps and there are no cycles.
2. Donors: CSR by counting, `donor_offsets` (`n + 1`) and `donors`
   (`n −` the number of base levels); the donors of a cell in ascending
   index order.
3. Stack, breadth-first: the base levels (`r(i) = i`) in index order, then
   repeatedly the donors of the next cell of the stack appended, using the
   stack itself as the queue. Each cell is appended once, after its
   receiver; along the stack the number of links to the base level never
   decreases. Braun & Willett (2013) build their stack depth-first; the
   ticket asks for breadth-first, which has the property both their
   accumulation and their implicit solve need (each cell after its
   receiver) and needs no recursion.
4. `accumulate`: `a = input`, then for each cell of the stack after the
   base levels, in reverse, `a[r(i)] += a[i]`. Each cell's total is
   complete before it is passed on, because its donors come after it in
   the stack. Drainage area: `input_i = area_sr(i) · R²` (m²); discharge:
   `input_i = P_i · area_sr(i) · R²` with `P_i` the precipitation rate
   (m/s), giving m³/s.

Time and memory are O(n): each cell's neighbours are read once, and every
other step visits each cell once.

### Tests (`drainage/tests.rs`)

- By hand on the octahedron (all centre distances π/2, so the steepest
  neighbour is the lowest): an ocean case (an outlet at sea level exactly,
  an ocean cell with a lower ocean neighbour, land draining to the ocean),
  a case of land draining through another land cell to the ocean (with one
  ocean cell: every other cell neighbours an ocean cell when there are
  two), and land cases with ties and a flat (two pits); receivers,
  donors, stack, base levels and pits listed by hand; `accumulate` with
  inputs `2^i`, whose exact sums name the cells upstream.
- On a 1,000-cell grid, a cell whose nearer lower neighbour is steeper
  than its farther, lower one: the receiver is the steeper.
- Invariants on 10,000-cell surfaces: the Earth worlds of seeds 1–3
  through the Crust stage, uniform random heights, random heights from
  {0, 1, 2, 3} m (flats and ties everywhere) with sea level 0, and an
  all-flat land surface:
  - receivers: the definition above checked cell by cell;
  - no cycles: from every cell, following receivers reaches a base level
    within `n` steps;
  - donors: exactly the inverse of the receivers, ascending;
  - stack: every cell exactly once, each after its receiver, base levels
    first in index order, link depth non-decreasing (breadth-first);
    `base_levels` and `pits` as defined;
  - `accumulate`: with input 1 per cell, the exact number of cells
    upstream, checked against counts made by walking the receivers; drainage
    area summed over the base levels equals the sum of the inputs and 4πR²
    within 1e-12 relative; each cell's total equals its input plus its
    donors' totals within 1e-12 relative; discharge of random precipitation
    summed over the base levels equals `Σ P_i a_i R²` within 1e-12 relative.
- Synthetic surfaces on the 10,000-cell grid (`u` the unit tangent from a
  cell towards its receiver; the alignment bound is
  `acos(cos(G/2) − 2δ)`, `G` the widest angle between the directions of
  consecutive neighbours of the cell, because some neighbour lies within
  `G/2` of the downhill direction and the receiver is at least as steep,
  where the slope towards a neighbour at the angle `φ` from downhill is
  `|∇h| (cos φ + ε)`, `|ε| ≤ δ`, with `δ` derived from the geometry of
  each surface, `ℓ` the farthest neighbour; half the widest gap plus 2°,
  the first version of this plan, is reported too; review round 1):
  - tilted plane `h = 1,000 m · (p·d)`, sea level −500 m: no pit, every land
    cell drains to the ocean, every link goes down the plane, and for land
    cells 30° or more from `d`, `u` lies within the alignment bound of the
    downhill direction `−(d − (d·p)p)`, with `δ = ℓ²/6 + ℓ |cot ψ|/2`
    (`ψ` the angle from `d`; the slope is exactly
    `|∇h| (sin ℓ/ℓ · cos φ + (1 − cos ℓ)/ℓ · cot ψ)`);
  - cone `h = 2,000 m · (1 − θ/1 rad)`, `θ` the angle from the apex, sea
    level 0: no pit, every step goes outward (`θ` grows), `u` within the
    alignment bound of the outward radial direction for cells with
    `θ ≥ 0.2 rad`, with `δ = ℓ/2 · max(|cot(θ − ℓ)|, |cot(θ + ℓ)|)`
    (`d²θ/ds² = sin²φ cot θ` along a great circle), the water reaches the
    sea in each of 12 sectors of 30° around the apex, and no basin drains
    more than 1/12 of the cone's area;
  - bowl inside a rim (`h = 1,000 m · (1 − |θ − β|/β)`, `β = 0.5 rad`, sea
    level −100 m): exactly one pit, the cell nearest the centre; every cell
    more than two mean spacings inside the rim drains to it, every land cell
    more than two spacings outside drains to the ocean, and the pit's
    drainage area lies between the cap areas at `β ∓ 2` spacings and
    between the areas of those two sets of cells. A cell within `β − ℓ/2`
    drains inwards and one beyond `β + ℓ/2` outwards, and a Voronoi cell
    reaches at most `ρ` (the largest circumradius of the Delaunay
    triangles) from its centre, so the test also checks
    `L/2 + ρ < 2` spacings, `L` the farthest neighbour of all cells.
- Panics: wrong lengths, a non-finite height or sea level, a radius that is
  not finite and positive.
- Golden fingerprint: FNV-1a 64 over the receivers, the stack and the bits
  of the drainage area of the Earth world of seed 1 at 10,000 cells,
  recorded only after the tests above pass and the maps are viewed.

### Maps and timing

A scratch tool outside the repository (a Cargo package in the scratchpad
with path dependencies on `kuniumi-sim`, `kuniumi-io`, `kuniumi-sphere` and
`kuniumi-math`) draws equirectangular maps into `out/M3-01/`: the base-10
logarithm of the drainage area on the plane, the cone, the bowl, a smooth
random surface and the Earth world of seed 1 at 650,000 cells (ocean dark,
pits marked); and times `build` plus `accumulate` (area and discharge) in
release at 100,000, 650,000 and 2,600,000 cells, best of 5.

## Acceptance criteria

- [x] No cycles: following receivers from any cell reaches an outlet or a pit
      within n steps (on the surfaces of the invariant tests).
- [x] The stack contains every cell exactly once, each after its receiver
      (same surfaces), with link depth non-decreasing along it.
- [x] Accumulated area summed over outlets and pits equals the total area within
      1e-12 relative (Σ of the inputs and 4πR²).
- [x] Synthetic surfaces: a tilted plane drains one way, a cone drains radially,
      a bowl produces a single pit, as measured in the Tests notes.
- [x] O(n) time; 650,000 cells in under 1 s in release (time recorded, with
      100,000 and 2,600,000 cells for the scaling).
- [x] Maps of drainage area for the synthetic surfaces and a 650,000-cell
      world viewed and described.
- [x] `mise run ci` passes.

## Verification

```sh
mise exec -- cargo test -p kuniumi-sim drainage
mise exec -- cargo test -p kuniumi-sim --release drainage
```

Maps in `out/M3-01/` from the scratch tool (see Maps and timing).

## Result

What was done:
- `kuniumi-sim`: the new module `drainage` (`drainage.rs`), with
  `Drainage::build`, `len`, `is_empty`, `receiver`, `receivers`, `donors`,
  `stack`, `base_levels`, `pits` and `accumulate`, as planned; 16 tests in
  `drainage/tests.rs`. `lib.rs` lists the module and its documentation
  names it; the documentation of `stages::erosion` points to it.
- `CHANGELOG.md` (Added).

Deviations from the plan (the code follows it; the tests and maps differ
in these points, the Tests notes now say what the tests check; review
round 1, R1-1 and R1-2):
- The alignment of the links on the plane and the cone is checked against
  a bound derived from the exact geometry of each surface,
  `acos(cos(G/2) − 2δ)`, instead of half the widest gap plus 2°: a fixed
  2° is no bound where the sphere or the contours curve strongly (30° from
  `d` on the plane, near the apex of the cone), while the derived one
  holds for every cell. Half the widest gap plus 2° is reported as well:
  no checked cell lies beyond it (smallest margins 1.71° on the plane and
  1.91° on the cone).
- The bowl's cap bracket is checked together with the condition that
  makes it a consequence of the drain-to checks (`L/2 + ρ` < 2 spacings),
  and with the exact bracket by the areas of the two sets of cells.
- The case of land draining through land to the ocean is its own
  octahedron test: with two ocean cells, as the ocean case needs, every
  land cell of the octahedron neighbours an ocean cell and drains to it.
- Maps: the plane, the cone and the bowl as orthographic globes (centred
  in the band of the plane and on `d`, on the apex, on the centre), which
  show the one-way and radial patterns without the distortion of an
  equirectangular map near the poles; the smooth random surface and the
  Earth world equirectangular, with two close-ups of the Earth world. Land
  is coloured by the logarithm of the number of cells upstream (count, not
  area: the same picture on a grid of near-equal cells, and the scale of
  every map starts at one cell).

Evidence:

1. No cycles, stack, base levels, pits, donors and accumulation
   (`invariants_on_real_and_random_surfaces`): on the Earth worlds of seeds
   1–3 at 10,000 cells, uniform random heights, random heights from
   {0, 1, 2, 3} m (flats, ties, a quarter of the cells exactly at the sea
   level), flat land and ocean everywhere. For every cell: the receiver is
   the steepest lower neighbour by the formula, ties to the lower index, or
   the cell itself for an ocean cell or a land cell with no lower
   neighbour; following receivers reaches a base level within `n` links
   (longest paths 3, 6 and 5 links on the Earth worlds, 3 on the random
   surfaces); the donors are exactly the inverse, ascending; the stack
   holds every cell once, each after its receiver, the base levels first
   in index order, then ordered by the position of the receiver and the
   index (breadth-first), with the number of links never decreasing;
   `base_levels` and `pits` as defined. Counts (input 1) equal, exactly,
   the counts made by walking every path; summed over the base levels they
   are `n` exactly. Drainage area summed over the base levels equals the
   sum of the cell areas and 4πR² within 1e-12 relative, and each cell's
   total its own area plus its donors' totals within 1e-12 relative;
   discharge of random precipitation (0–3 m/yr) summed over the base
   levels equals `Σ P_i a_i R²` within 1e-12 relative. Earth worlds: 2,911,
   2,868 and 2,889 land cells, of which 1,832, 1,669 and 1,731 are pits
   (the flat +300 m continental interiors of M2-05); uniform random: 10
   pits among 5,055 land cells (about `5,000/(2⁶·7)` = 11 expected for six
   neighbours); random steps: 461.
2. By hand (`octahedron_ocean_by_hand`,
   `octahedron_land_through_land_to_the_ocean_by_hand`,
   `octahedron_land_by_hand`, `octahedron_flat_by_hand`): receivers,
   donors, stack, base levels, pits and the exact totals of the inputs
   `2^i` for an outlet exactly at sea level, an ocean cell next to a lower
   one, land draining through land to the ocean and to a pit, a flat of
   two pits, and ties that go to the lower index although the
   higher comes first in the ring. `the_receiver_is_the_steepest_not_the_lowest`:
   on a 1,000-cell grid the nearer, steeper neighbour wins over a farther,
   lower one, and the farther wins once it is steeper.
3. Synthetic surfaces (10,000 cells):
   - `a_tilted_plane_drains_one_way` (`h = 1,000 m · p·d`, sea level
     −500 m): no pit, every land cell drains to the ocean and every link
     goes down the plane; every land cell 30° or more from `d` (6,827
     cells) has its link within the alignment bound of the downhill
     direction (mean 15.97° from it; smallest margin to the bound 0.26°;
     none beyond half the widest gap plus 2°, smallest margin 1.71°).
   - `a_cone_drains_radially` (`h = 2,000 m · (1 − θ/1 rad)`, sea level
     0): no pit, every link goes outward, the coast receives water in all
     12 sectors of 30°, every cell at least 0.2 rad from the apex (2,197)
     has its link within the alignment bound of the outward direction
     (mean 16.65° from it; smallest margin 3.09°; none beyond half the
     widest gap plus 2°, smallest margin 1.91°), and the largest basin
     drains 2.92 % of the cone (bound 1/12, 8.33 %).
   - `a_bowl_produces_a_single_pit` (crater, rim at 0.5 rad, sea level
     −100 m): one pit, the cell nearest the centre (8,904);
     `L/2 + ρ` = 1.968 spacings; every cell more than two spacings inside
     the rim drains to the pit, every land cell more than two spacings
     outside to the ocean; the pit's area is 0.9999 of the cap inside the
     rim, between the caps at the rim ∓ 2 spacings (0.7406 and 1.2954 of
     it) and between the areas of the two sets of cells.
4. Panics: `build_needs_one_height_per_cell`, `build_needs_finite_heights`,
   `build_needs_finite_heights_not_infinite`,
   `build_needs_a_finite_sea_level`, `build_needs_a_finite_positive_radius`
   (0, −1, +∞, NaN), `accumulate_needs_one_value_per_cell`.
5. Golden value (`golden_fingerprint_of_an_earth_world`, the Earth world
   of seed 1 at 10,000 cells, FNV-1a 64 over receivers, stack and area
   bits): `0xB08B_3EE0_FB8F_353F`, recorded after 1–4 passed and the maps
   below were viewed; the same in debug and release.
6. Timing (release, this machine; the scratch tool below; `build` plus
   `accumulate` of area and of discharge, best of 5):

   | Cells | Earth world, seed 1 | Uniform random heights |
   | --- | --- | --- |
   | 100,000 | 1.6 ms (15.6 ns/cell) | 3.7 ms (36.6 ns/cell) |
   | 650,000 | 10.6 ms (16.2 ns/cell) | 24.6 ms (37.9 ns/cell) |
   | 2,600,000 | 42.0 ms (16.2 ns/cell) | 102.7 ms (39.5 ns/cell) |

   Linear in `n`; 650,000 cells far under 1 s. `build` is most of it
   (9.1 and 22.9 ms at 650,000 cells), one `accumulate` 0.8–0.9 ms.
7. Maps (`out/M3-01/`, a scratch tool outside the repository,
   `<scratchpad>/m3_01_maps`, a Cargo package with path dependencies on the
   workspace crates; ocean dark blue, pits red, land by the logarithm of
   the number of cells upstream in viridis), viewed:
   - `plane_100000.png` (globe centred in the band): parallel flow lines
     from the top of the disc down to the ocean cap at the bottom,
     converging into brighter streams near the coast; above, near the top
     of the plane, the lines fan out from `d`. Some sectors of parallel
     streaks follow the rows of the grid, the path bias of ticket M1-11.
     `plane_100000_top.png` (centred on `d`): flow radiates from the top.
   - `cone_100000.png` (the cone of the test, coast at 1 rad): streaks
     radiate from the apex in every direction to the coast ring, joining
     into brighter streams towards it; wedges of parallel streaks as in
     M1-11's cone images for jitter 0.35; no pit.
   - `bowl_100000.png`: a dark ring (the divide, count 1) on the rim; inside
     it the flow converges on the single red pit at the centre, the streams
     joining in branches near it; outside it the flow runs outward to the
     ocean ring.
   - `waves_100000.png` (a smooth random surface, a sum of 24 plane waves,
     sea level −150 m): divides radiate from the peaks, and the flow
     gathers into valleys and trunk streams (yellow) that reach the coast;
     a tree-like network, as expected of steepest descent on smooth relief.
     The same four at 10,000 cells look alike, coarser.
   - `earth_650000_1.png` (the Earth world of seed 1, 650,000 cells:
     195,865 land cells, 150,176 pits, 77 % of the land): the continental
     interiors are red, every cell a pit, being flat at +300 m; the margins
     drain to the sea in short streaks across the coast; the largest
     basin on land counts 162 cells. `earth_650000_1_arc.png` (15° S,
     130° E): the crest of an arc is a dark divide along the coast, its
     seaward flank draining to the sea and its landward flank into the
     pits of the flat interior. `earth_650000_1_plateau.png` (35° S,
     130° W): a flat plateau of pits, strips of sloping land around it
     draining to the sea and to the flats. Routing out of the pits is
     M3-02.
8. Test-first: against the stubs (commit "tests of the drainage graph
   against stubs") all 15 tests then present failed, the six panic tests
   because the stub's panic does not carry the expected message.
   Mutations, each caught by the final 16 tests (release): ties to the
   higher index (the two land octahedron tests with ties); the lowest
   instead of the steepest neighbour (five tests: the steepest-not-lowest
   test, the plane, the cone, the invariants and the golden value); ocean
   only strictly below the sea level (`octahedron_ocean_by_hand`);
   accumulation in forward order (five tests); donors in descending order
   (six tests). The round-2 reviewer ran 13 mutants of its own, all
   caught.
9. `mise run ci` passes: fmt, clippy `-D warnings`, rustdoc, the browser
   checks, the tests in debug and release, MSRV 1.85.0, cargo-deny.

Review:
- Round 1 (CHANGES_REQUESTED): R1-1 major, the cone test did not check
  what the plan asked (profile θ/0.5 rad instead of θ/1 rad, the 1/12
  basin bound only printed, the alignment check replaced by a check of
  the azimuth of the outlet); fixed: the plan's cone, the bound asserted,
  the alignment checked against a bound derived from the geometry, the
  plan's half gap + 2° reported. R1-2 major, "Deviations: none" while the
  tests differed from the plan (the plane only in a band, the bowl
  without the cap areas, no land-through-land case reaching the ocean, map
  projections); fixed: the plane on all land cells 30° or more from `d`,
  the bowl's cap bracket with the condition that proves it, a new
  octahedron test, the Tests notes describing the final checks and the
  remaining deviations listed with their reasons. R1-3 minor, the
  erosion docs listed stack and discharge before depression routing,
  against the design's order (fixed). R1-4 nit, a local FNV-1a duplicated
  `kuniumi_math::rng::fnv1a64` (fixed; golden value unchanged). The
  reviewer reproduced the timing, the maps (byte-identical) and the
  plan's measurements with its own checker.
- Round 2 (APPROVE): every round-1 finding resolved; the reviewer
  re-derived the plane, cone and bowl bounds, reproduced every quoted
  number, the timing and the maps, and caught 13 of 13 mutants. R2-1 nit,
  an unwrapped line and a grammar slip in the ticket; R2-2 nit, the
  mutation counts predated round 1 (both fixed).

Follow-ups: none.

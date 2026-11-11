# M2-03 Continents and boundary classification

- Status: done
- Milestone: M2
- Depends on: M2-02
- Crates: kuniumi-sim (stages/tectonics), kuniumi-cli (tests)
- Design: プレートテクトニクス / 境界ごとの処理; プレートテクトニクス / 初期条件の目安
- Decisions: none

## Goal

Continental crust covers the requested fraction of the surface, and every plate
boundary is classified as convergent, divergent or transform, with a rate.

## Scope

In:
- Continental cells: grow a few continental blobs inside plates (stream
  `"tectonics/continents"`) until the area-weighted continental fraction
  reaches `PlanetParams::continental_crust_fraction`. Some plates end up mostly
  oceanic.
- Boundary edges: for each neighbour pair on different plates, the relative
  velocity `v_j − v_i` projected on the edge normal (the great-circle direction
  from i to j) gives the convergence rate (m/s, positive when converging); the
  tangential part gives the shear. Convergent when the normal part dominates
  and is positive, divergent when it dominates and is negative, transform
  otherwise; the threshold angle is documented.
- Outputs: `crust_type` (class: 0 oceanic, 1 continental), `boundary_type`
  (class per cell: none, convergent, divergent, transform; the strongest of
  its edges) and `convergence_rate` (m/s).

Out:
- Thickness, age and elevation (M2-04, M2-05).

## Implementation notes

Gaps found while refining (the ticket was `ready`): the boundary method as
written cannot meet the "pure boundaries" criterion (below, put to the
owner), and the ticket leaves open how many blobs "a few" is, what
"strongest" means and the threshold angle. Measured with a throwaway
prototype in the scratchpad on the real grid (`SphereGrid::fibonacci`) and
on the real partitions and poles of M2-02.

### Question for the owner: the boundary normal

Two plates split by a great circle, with hand-chosen relative poles (pole
on the circle's axis: pure transform; pole on the circle: pure normal
motion, convergent on one side of the pole and divergent on the other;
cells within 10° of the relative pole, where the motion vanishes, left
out). Share of boundary cells classified as expected, 45° threshold, grid
seeds 1–5, four orientations of the circle:

| Normal of a boundary cell | 10,000 cells | 100,000 cells |
| --- | --- | --- |
| Literal: per edge, the direction from i to j; the cell takes its strongest edge | 67–74 % | 68–72 % |
| Per cell: sum of its Voronoi edges' normals towards the other plate, weighted by edge length | 97.7–99.7 % | 98.2–99.4 % |
| Recommended: that sum, added up along the boundary over 0.08 rad (about 510 km on Earth) | 100 % | 100 % |

A single Voronoi edge of a staircase boundary is often oblique to the
boundary as a whole (a third of the cells on these grids), so the literal
method mixes the three types at the scale of single cells along every
boundary. On the real 100,000- and 650,000-cell worlds of seed 1 the
literal map is speckled; with the recommended window convergent arcs,
ridges and transforms form segments hundreds to thousands of kilometres
long, at both resolutions. The window also gives the normal rate within
6.4 % (10,000 cells) and 0.6 % (100,000 cells) of the exact one in the
two-plate cases.

Recommended: the windowed normal below, with the 45° threshold and the
continents below. Alternatives: the literal method, with the criterion
relaxed to "about 70 % of the boundary cells"; or the per-cell normal
without a window, with the criterion relaxed to 97 %.

### Decision (owner, 2026-11-11)

The owner chose the recommended windowed normal: the per-cell,
edge-length-weighted normal summed along the boundary within 0.08 rad, and
approved the rest of this plan as written (blobs inside plates, plates taken
in random order until they cover twice the continental area, the 45°
threshold, the pair of largest relative speed for cells touching several
plates).

### Continents (`continents`)

```rust
pub const CONTINENTS_STREAM: &str = "tectonics/continents";
pub const CRUST_TYPE: &str = "crust_type";        // class: 0 oceanic, 1 continental
pub const CONTINENT_PLATE_COVER: f64 = 2.0;
pub fn continents(grid: &SphereGrid, plates: &[u32], plate_count: u32,
                  fraction: f64, seed: u64) -> Vec<bool>;   // true: continental
```

All random numbers from `CONTINENTS_STREAM`, in this order:

1. Slowness per cell, in cell order: `w_i = exp(3 (u − ½))`, as in the
   partition (`SLOWNESS_SPREAD`), so coastlines are irregular.
2. A random order of the plates: Fisher–Yates, for `i` from `K − 1` down to
   1, swap `i` with `rng.below(i + 1)`.
3. For each blob, when it starts: its seed, a uniformly random cell of its
   plate (`rng.below(plate size)` into the plate's cells in index order),
   and its growth rate `r_b = 1 + 2u` in `[1, 3)` (`MAX_GROWTH_RATE`).

Let `A = Σ a_i` (index order) and the target `A* = f A`. For `f > 0`,
blobs start on the plates of that order, one per plate, until the plates
taken cover at least `min(1, 2f) A` (`CONTINENT_PLATE_COVER = 2`: a
continent covers about half its plate on average; for Earth's `f = 0.4`
about 9–12 of 12 plates, measured over seeds 1–5). The blobs grow together
as in the partition (multi-source Dijkstra, priority `p + d_ij (w_i + w_j)
/ 2 / r_b`, heap of `(priority via total_cmp, cell, blob)`), each only into
cells of its own plate that are not continental yet, while the continental
area (summed in growth order) is below `A*`. If the heap empties first
(every blob has filled its plate), a blob starts on the next plate of the
order; when no plate is left every cell is continental. So the fraction is
within one cell's area above `f` (for `f = 0` no cell, for `f = 1` every
cell), each plate's continental cells are connected, small plates may be
filled, and the plates left out stay oceanic. Measured at `f = 0.4`, seeds
1–5, 10,000 and 100,000 cells: the most oceanic plate is 0–13 %
continental.

### Boundaries (`boundaries`)

```rust
#[repr(u8)]
pub enum BoundaryType { None = 0, Convergent = 1, Divergent = 2, Transform = 3 }
pub const BOUNDARY_TYPE: &str = "boundary_type";       // class
pub const CONVERGENCE_RATE: &str = "convergence_rate"; // m/s
pub const NORMAL_WINDOW_RAD: f64 = 0.08;
pub struct Boundaries { pub kind: Vec<BoundaryType>, pub convergence_rate_m_s: Vec<f64> }
pub fn boundaries(grid: &SphereGrid, plates: &[u32], poles: &[EulerPole],
                  radius_m: f64) -> Boundaries;
```

A cell `i` of plate `a` is a boundary cell of the pair `(a, b)` if a
neighbour is on plate `b`; the band of `(a, b)` is the boundary cells of
the pair on both sides.

1. Raw normal: `g_ib = Σ ℓ_ij u_ij` over the neighbours `j` on plate `b`,
   with `u_ij = (p_j − p_i)/|p_j − p_i|` the direction of the great circle
   from `i` to `j` at the edge's midpoint (the normal of the Voronoi edge)
   and `ℓ_ij` the edge length (rad). It is (up to the factor `−2 a_i`) the
   Green–Gauss gradient of plate `b`'s indicator on the cell, and points
   from `a` to `b`.
2. Window: `N_ib = Σ g_k` over the cells `k` of the band reached from `i`
   by a breadth-first search through neighbours in the band whose centres
   lie within `NORMAL_WINDOW_RAD` of `p_i`, each oriented from `a` to `b`
   (`g_kb` for `k` on `a`, `−g_ka` for `k` on `b`), summed in search order.
   `n̂` is `N_ib` projected on the tangent plane at `p_i` and normalised; if
   that is zero (a one-cell plate surrounded by `b`, as on the octahedron)
   the cell's unit direction to its first neighbour on `b` (CSR order) is
   used. `0.08 rad` is a few cells at 10,000 cells (spacing 0.035 rad) and
   a window in radians makes the types independent of the radius (the
   poles scale with `1/R`).
3. Relative velocity of `b` with respect to `a` at the cell:
   `Δv = (ω_b − ω_a) × (R p_i)`, `ω = angular_speed · axis` (m/s).
   Convergence rate `c = −Δv · n̂` (positive when the plates approach),
   shear `s = Δv · (p_i × n̂)`.
4. Type: convergent if `c > |s|`, divergent if `−c > |s|`, transform
   otherwise: the threshold is 45° between the relative velocity and the
   boundary normal, the angle at which the normal part stops dominating.
   No relative motion (`c = s = 0`) counts as transform.
5. A cell touching several plates takes the pair with the largest relative
   speed `|Δv|` ("the strongest"); a tie keeps the lower plate id. Cells
   that touch no other plate: `None`, rate `+0`.

`convergence_rate` holds `c` of the chosen pair (also for transform cells,
whose `|c| ≤ |s|`), so M2-04 can read the half spreading rate `−c/2` of
divergent cells.

### Stage

`run` adds, after the three fields of M2-02: `crust_type` (class),
`boundary_type` (class), `convergence_rate` (`"m/s"`), from
`world.params.continental_crust_fraction`.

### Tests

- `stages::tectonics::tests`: the new stream name; continental fraction
  within ±1 percentage point of 0.1, 0.4 and 0.7 for seeds 1–20 at 10,000
  cells (and exactly 0 and 1 for `f = 0` and `f = 1`); each plate's
  continental cells connected; at `f = 0.4` some plate at most 20 %
  continental for every seed; same seed, same bits.
- Two plates split by great circles of four orientations at 10,000 and
  100,000 cells: relative pole on the axis → every boundary cell
  transform; relative pole on the circle → every boundary cell farther
  than 10° from it convergent or divergent as the exact relative velocity
  says, both signs of the rotation, the rate within 10 % of `|Δv|` of the
  exact normal rate; swapping the two plates' ids gives the same bits;
  the octahedron fallback.
- Real partitions: `None` exactly for cells without a neighbour on another
  plate; the same types and rates (within 1e-12 relative) for Earth's
  radius and half of it; the stage adds its fields in order with their
  units.
- CLI: the six tectonics fields written; golden fingerprints of
  `crust_type`, `boundary_type` and `convergence_rate` for 10,000 cells,
  seed 1, recorded after the invariant tests and the visual check.

## Acceptance criteria

- [x] Continental area fraction within ±1 percentage point of targets 0.1, 0.4
      and 0.7 (seeds 1–20, 10,000 cells).
- [x] Synthetic two-plate tests with hand-chosen poles give pure convergent,
      pure divergent and pure transform boundaries: every boundary cell of
      the expected type (cells within 10° of the relative pole left out), at
      10,000 and 100,000 cells.
- [x] `boundary_type.png` viewed: boundaries are continuous lines along plate
      edges.
- [x] Golden fingerprints for 10,000 cells, seed 1.
- [x] `mise run ci` passes.

## Verification

```sh
cargo run -p kuniumi-cli -- generate --seed 1 --cells 100000 --out out/M2-03
```

## Result

What was done:
- `kuniumi-sim`, `stages::tectonics`: the submodules `continents`
  (`continents`, `CONTINENTS_STREAM`, `CRUST_TYPE`,
  `CONTINENT_PLATE_COVER`) and `boundaries` (`boundaries`, `Boundaries`,
  `BoundaryType`, `BOUNDARY_TYPE`, `CONVERGENCE_RATE`,
  `NORMAL_WINDOW_RAD`), re-exported from `tectonics`, each with its tests
  (`continents/tests.rs`, `boundaries/tests.rs`); the heap entry of the
  partition is shared (its third key renamed from `plate` to `source`, no
  change in behaviour). `run` adds `crust_type`, `boundary_type` and
  `convergence_rate` after the three plate fields; the module, `run` and
  `generate` documentation name them.
- `kuniumi-cli`: tests only (the eight fields written, the continental
  fraction and boundary types of the 1,000-cell world, the golden values).
- `CHANGELOG.md` (Added, Changed).

Deviations from the plan:
- `Boundaries` also holds `other_plate`, the plate across each boundary
  cell (`None` elsewhere): it lets the tests check that the strongest pair
  was taken, and M2-04 needs the plate (and so the crust) across.
- The fallback normal is used when the tangential part of the window sum
  is shorter than 1e-9 of the sum (the plan said "zero"), so that a radial
  sum with rounding noise does not give a random direction.
- The "swap the ids" test renumbers the twelve plates of real partitions
  (ids reversed, poles with them) instead of two hemispheres, where the
  renumbering is trivially the same computation.

Evidence:

1. Continental fraction (`continents::tests::the_continental_fraction_meets_the_target`,
   seeds 1–20 at 10,000 cells, 12 plates): every fraction lies in
   `[f, f + largest cell / total]` and within ±0.01 of `f`; the largest
   differences are 9.6e-5 (target 0.1), 1.4e-4 (0.4) and 1.1e-4 (0.7).
   `no_fraction_gives_no_continent_and_a_whole_fraction_every_cell` (0 and 1
   at 10,000 cells, and at 1,000 cells with 2 and 64 plates);
   `each_plate_has_one_connected_continent` (breadth-first search within
   each plate's continental cells, seeds 1–20 at 10,000 cells and 1,000
   cells with 2, 12 and 64 plates, fractions 0.1, 0.4, 0.7);
   `continents_start_on_plates_covering_twice_their_area` (plates holding
   continent cover at least `min(1, 2f)` of the area; for `f ≤ 0.25` some
   plate has none); `some_plates_stay_mostly_oceanic` (at 0.4 every seed
   has a plate without continent, and 1 to 6 plates at most 20 %
   continental); `the_same_seed_gives_the_same_continents`;
   `continents_rejects_*` (wrong length; −0.1, 1.5 and NaN). The 1,000-cell
   world of the CLI test is 40 % continental within one percentage point.
2. Pure boundaries (`boundaries::tests`, two plates split by four great
   circles, at 10,000 cells with grid seeds 1–3 and 100,000 cells with
   seed 1, each with and without a rotation common to both plates):
   `a_relative_pole_on_the_axis_gives_pure_transform` (every cell next to
   the other plate is transform, facing the other plate; no other cell is
   a boundary); `a_relative_pole_on_the_circle_gives_pure_convergence_and_divergence`
   (both signs of the rotation: every boundary cell farther than 10° from
   the relative pole is convergent or divergent as the exact relative
   velocity across the exact normal says, with the rate within 10 % of
   `|Δv|` of the exact normal rate; both types present).
   `a_plate_of_one_cell_takes_the_direction_of_a_neighbour` (octahedron,
   worked out by hand: types divergent, none, divergent, convergent,
   transform, transform and rates `−ωR, 0, −ωR, ωR, 0, 0`).
   On real partitions (seeds 1–5 at 10,000 cells):
   `boundary_cells_are_the_cells_next_to_another_plate` (boundary cells are
   exactly those with a neighbour on another plate; elsewhere rate `+0`,
   no plate across; the plate across is a neighbour's plate of largest
   relative speed; convergent rates positive, divergent negative; counts
   none 42,534, convergent 1,774, divergent 1,762, transform 3,930);
   `renumbering_the_plates_changes_nothing` (bits equal);
   `the_boundaries_do_not_depend_on_the_radius` (Earth's radius and half
   of it: same types and plates across, rates within 1e-12 of the largest).
   `tectonics::tests::the_stage_adds_its_fields` (six fields, names, units,
   kinds, order and values).
3. Maps (`kuniumi generate --seed 1 --cells 100000 --out out/M2-03/100k`,
   release, exit code 2 naming `crust`), viewed:
   - `boundary_type.png`: every plate outline of `plate_id.png` drawn as a
     continuous line two cells wide, closed around each plate and meeting
     at triple junctions, with no boundary cell inside a plate. Along the
     lines the types form segments hundreds to thousands of kilometres
     long: a long convergent arc (dark blue) from about 150° W to 70° W
     at 35° to 70° S, another along the north side of the large southern
     plate (about 15° W to 55° E),
     divergent ridges (light green) along the west side of the red plate
     and across the south, transform (dark green) along the north polar
     plate; at the scale of single cells some steps of transform remain
     inside ridges.
   - `convergence_rate.png` (viridis over its range, 0 is the teal
     interior): yellow where `boundary_type` is convergent, dark purple
     along the southern ridge, mid tones along transforms.
   - `crust_type.png`: about eight irregular continents (dark blue) of very
     different sizes, coastlines rough at the scale of a few cells; some
     are cut straight by plate boundaries (the dark blue oval plate is
     entirely continental), a south polar continent, and wide oceans
     (light blue) over two plates without continent.
   - `out/M2-03/10k/` (the golden case): the same features at 2° cells;
     boundary lines continuous, segments of each type several cells long;
     five or six landmasses, some with lakes of one cell.
   - `out/M2-03/650k/` (650,000 cells, seed 1, another partition):
     boundary segments of each type up to thousands of kilometres, as at
     100,000 cells; continents smoother relative to their size.
4. Golden fingerprints (`tests/cli.rs::generate_golden_fingerprints`,
   10,000 cells, seed 1): `crust_type` `0x791F_6B17_52CF_D2B8`,
   `boundary_type` `0x5163_2C85_7C3F_C400`, `convergence_rate`
   `0x04BB_B545_5B56_84CC`, recorded after the tests above passed and the
   maps were looked at; the five earlier values unchanged.
   `generate_writes_the_same_bytes_twice` covers the 17 files.
5. `mise run ci` passes: fmt, clippy `-D warnings`, the tests in debug and
   release (353 in debug, 344 in release), MSRV 1.85.0, cargo-deny
   (advisories, bans, licences, sources ok). `cargo doc -p kuniumi-sim`
   gives no warnings. Test-first: with `continents` and `boundaries`
   stubbed (no continent, no boundary) nine of the new tests failed.

Timing (release, this machine): `kuniumi generate --seed 1 --cells 650000`
takes 4.2 s and 177 MB (3.9 s and 175 MB after M2-02); `continents` takes
0.13 s and `boundaries` 0.07 s at 650,000 cells (0.012 s and 0.010 s at
100,000). At 100,000 cells `generate` takes 0.83 s and 44 MB.

Review:
- Round 1 (APPROVE): the reviewer re-ran `mise run ci`, regenerated the
  maps at 10,000, 100,000 and 650,000 cells and confirmed the descriptions,
  checked fractions from 1e-9 to 0.999999 with 2, 12 and 64 plates, ran
  the two-plate cases on untested grids (10,000 cells seeds 4–10, 100,000
  cells seeds 2–3, 650,000 cells seed 1) with random circles and common
  rotations (no cell misclassified), and mutated the plate restriction,
  the window and the orientation (each caught by the tests). R1-1 minor,
  continents enclose small oceanic pockets (single cells mostly; 12–15 at
  10,000 cells, 116–146 at 650,000) that M2-05 would turn into pits
  (recorded as the draft M2-10, and in the documentation of
  `continents`). R1-2 minor, below about 5,000 cells the window holds only
  the cell and the types are noisier (up to 2.5 % wrong at 1,000 cells;
  documented on `NORMAL_WINDOW_RAD`; a window floored at a few cells
  would be a design change for the owner). R1-3 minor, the `Crates:` line
  missed the CLI tests (fixed). R1-4 nit, the 10 % rate tolerance has
  little margin on untested grids (0.102 on one; the expected accuracy is
  now documented on `boundaries`). R1-5 nit, the "start another blob"
  branch is unreachable with `CONTINENT_PLATE_COVER = 2` and the tie-break
  of equal speeds is untested (the branch is documented as a guard; exact
  ties are practically impossible). R1-6 nit, a doc paragraph not rewrapped
  (fixed).

Follow-ups: M2-10 (draft, oceanic pockets enclosed by continents).

# M3-06 Convergent flow across flats

- Status: done
- Milestone: M3
- Depends on: M3-05, M3-08
- Crates: kuniumi-sim (stages::erosion), kuniumi-world (docs), kuniumi-cli (tests), docs
- Design: 侵食と地形発達 / 1ステップの処理
- Decisions: none

## Goal

Water crossing a flat (a drained flat or a filled depression) converges
into channels as on real low-gradient land, instead of running in the
straight parallel lines that the ε gradient of M3-02 gives.

## Scope

In:
- The decision, from the river maps of M3-05 and the measurements below,
  whether the parallel lines of the ε routing survive erosion and whether
  the flat resolution of Barnes, Lehman & Mulla (2014b) removes them;
  recorded with its evidence.
- Instead of that flat resolution, which does not remove them, an initial
  roughness of the land at the start of the Erosion stage (approved by the
  owner): a seeded perturbation of 1 m at most on every land cell, before
  the first erosion step.
- The design document (侵食と地形発達 / 1ステップの処理) states the initial
  roughness.
- Golden values of the Erosion stage and of `kuniumi generate` updated
  (world-breaking).

Out:
- The flat resolution of Barnes, Lehman & Mulla (2014b): measured below,
  it leaves the flow on the Earth worlds' flats parallel.
- The basin graph of Cordonnier, Bovy & Braun (2019).
- The roughness as a user setting (a constant here; the calibration of
  milestone 7 may make it one).
- Relief of the continental interiors from the crust itself (the uniform
  +300 m interiors of M2-05; milestones 6 and 7).

## Implementation notes

### Decision: what the parallel lines are, and what removes them

Measured with a scratch tool outside the repository
(`<scratchpad>/m3_06_proto`, a Cargo package with path dependencies on the
workspace crates) on the Earth world of seed 1 at 650,000 cells, default
settings; its erosion loop is that of `erode` (with the ε routing it
reproduces M3-05's world: the Hack exponent 0.5247028646816342 of M3-08
bit for bit, 180,665 land cells within 1 m of the sea level, 42,457 river
cells). Images in `out/M3-06/` (git-ignored).

1. **The parallel lines survive erosion.** After 1,000 steps the planed
   interiors carry trunks with pinnate tributaries parallel on each side
   (`out/M3-05/650k/closeups/east_rivers.png`, 4° N, 130° E) and long
   straight parallel rivers across the plains beside the collision belt
   (`out/M3-06/closeups_eps/plateau_rivers.png`, 35° S, 130° W); M3-08's
   Hack exponent is 0.51–0.52 for seeds 1–5, below Earth's 0.6.
2. **They are made in the first steps.** On the surface filled without ε,
   150,186 of the 195,865 land cells lie on 54 flats (the largest 28,259
   cells). The implicit incision lowers whole chains of receivers at once,
   so the flats vanish fast: 66,788 flat cells after 1 step, 22,611 after
   10, 649 after 30. The network the ε gradient draws in those steps is
   what the valleys then cut.
3. **The ε routing is Barnes et al.'s gradient towards lower terrain.** On
   that surface, the receivers of Priority-Flood+ε equal, all 650,000 of
   them, those of Barnes, Lehman & Mulla (2014b, "An efficient assignment
   of drainage direction over flat surfaces in raster digital elevation
   models", *Computers & Geosciences* 62, 128–135; arXiv:1511.04433) with
   the gradient towards lower terrain alone, which the paper (§2.1) says
   "results in parallel flows".
4. **Their full method does not converge the flow on these flats.**
   Prototyped from the paper (§2, Algorithms 1–7; the step on the surface
   filled without ε, receivers by the steepest drop of the mask): 49 of
   the 54 flats touch higher terrain (140,229 cells), the gradient away
   from it spans up to 232 cells and the one towards lower terrain up to
   90. Superposed, both are distance fields smooth over hundreds of cells,
   so the flow stays parallel except along the flats' medial lines: it
   changes the receivers of 35,360 flat cells (23.5 %) but not the
   pattern (`out/M3-06/flats/step0_{plateau,east}_{eps,barnes}_flow.png`).
   After 1,000 steps the river maps show the same pinnate tributaries and
   parallel rivers as those of ε (`out/M3-06/closeups_barnes/`,
   `closeups_s3_barnes/`), with Hack exponents 0.5187 (seed 1) and 0.5132
   (seed 3) against 0.5247 and 0.5093.
5. **A small roughness of the initial surface does.** Adding `A u_i`, `u_i`
   uniform in [0, 1) from a seeded stream, to every land cell before the
   first step (the ε routing unchanged) gives dendritic networks, trunks
   that wander and tributaries in every direction, without the pinnate
   pattern (`out/M3-06/closeups_eps_1/`, `closeups_s3_eps_1/`), and Hack
   exponents near Earth's:

   | Seed | ε (M3-08) | Barnes 2014b | roughness 1 cm | 1 m | 10 m |
   | --- | --- | --- | --- | --- | --- |
   | 1 | 0.5247 | 0.5187 | 0.6027 | 0.6110 | 0.6065 |
   | 2 | 0.5094 | | | 0.5819 | |
   | 3 | 0.5093 | 0.5132 | | 0.6126 | |
   | 4 | 0.5241 | | | 0.5982 | |
   | 5 | 0.5213 | | | 0.5892 | |

   The amplitude hardly matters: it only breaks the symmetry of an exactly
   level surface, which no real surface has (a cell of the default grid,
   28 km across, averages real relief of tens to hundreds of metres).
   It costs one pass over the cells, nothing per step.

Decided with the owner: no flat resolution; the initial roughness
instead, planned below. The alternatives put to the owner were to record
the decision only and move the roughness to a new ticket, to implement
Barnes et al. (2014b) as first planned, or to record the decision and
accept the parallel pattern until the interiors get relief from the
crust.

### Plan (initial roughness)

`kuniumi-sim`, `stages::erosion`:

```rust
/// The random stream of the initial roughness (docs/conventions.md,
/// "Determinism", rule 7).
pub const ROUGHNESS_STREAM: &str = "erosion/roughness";

/// The amplitude of the initial roughness (m): 1 m.
pub const INITIAL_ROUGHNESS_M: f64 = 1.0;

/// Adds `amplitude_m · u_i` to every land cell (`h_i > sea_level_m`), with
/// `u_i = rng.next_f64()` in [0, 1) drawn for every cell in index order,
/// land or ocean, so that a cell's perturbation does not depend on the
/// others. Ocean cells keep their heights bit for bit; land stays land.
pub fn roughen(height_m: &mut [f64], sea_level_m: f64, amplitude_m: f64, rng: &mut Rng);
```

- Panics if `amplitude_m` is negative or not finite, if a height or the
  sea level is not finite, or if a raised height is not finite.
- `run` (the stage): when `world.settings.erosion.steps > 0`, roughens
  the Crust stage's heights with `Rng::for_stream(world.seed,
  ROUGHNESS_STREAM)` and `INITIAL_ROUGHNESS_M` before `erode`; with 0
  steps the surface stays that of the Crust stage, as now. `erode` itself
  is unchanged (no randomness), so the golden values of M3-01 to M3-04
  stay.
- The module documentation and that of `run` describe the roughness and
  why (the decision above), citing a landscape-evolution source for
  starting from a slightly perturbed surface, checked in the original
  before it is cited (design document, 参考文献).
- `docs/design-ja.md`, 侵食と地形発達 / 1ステップの処理: one paragraph
  before the steps: the erosion starts from the surface of the Crust stage
  with a random roughness of at most 1 m on land, which decides the
  routing across flats that are exactly level; Priority-Flood+ε alone
  would draw parallel lines there (this ticket).

Tests (`stages/erosion/stage_tests.rs`, or a new `roughness_tests.rs`):

- `roughen_by_hand`: octahedron heights with ocean cells, land cells and a
  cell exactly at the sea level; the expected heights `h + A u` with `u`
  drawn in the test from a fresh `Rng::for_stream` of the same name; ocean
  and the cell at the sea level unchanged bit for bit; amplitude 0 leaves
  every height unchanged bit for bit.
- `roughen` panics: negative, NaN and infinite amplitudes; a NaN or
  infinite height or sea level; a land cell at `f64::MAX`.
- On the Earth worlds of seeds 1–3 at 10,000 cells: every land cell raised
  by at least 0 and less than `A`, the ocean unchanged, the land cells the
  same.
- `the_stage_erodes_and_routes_the_final_surface` and
  `the_stage_uses_every_erosion_setting` compare with `erode` of the
  roughened Crust surface; `zero_steps_keep_the_surface_of_the_crust_stage`
  unchanged.
- A synthetic flat: a cap of land at 300 m (radius 40°) in an ocean at
  −4,000 m, 40,000 cells, no uplift. Routed on the filled surface, the
  bare cap drains in parallel lines to the coast (largest basin 219 cells
  in the prototype); the roughened one converges (1,506). The test checks
  that the largest basin of the roughened cap is at least twice that of
  the bare one at step 0 and after 10 erosion steps (prototype: 1,104
  against 388), and the M3-02 invariants on the routed rough surface (no
  pit, every land cell reaches the ocean). (Review round 1, R1-2: over the
  roughness of five streams, every one at least twice at step 0; after 10
  steps the median at least twice and none below 1.5 times.)
- Golden values: the stage's four (`golden_fingerprints_of_the_erosion_stage`)
  and the CLI's four (`generate_golden_fingerprints`) recorded anew after
  the checks above and the visual check.

Evidence runs (release): `kuniumi generate --seed N --out out/M3-06/650k/seedN`
for seeds 1–5 (Hack exponents from `metrics.txt`, `rivers.png`), the
close-ups of item 5 (east, plateau, interior of seed 1; south pole of
seed 3), and the time of `kuniumi generate` for seed 1.

Size: about 60 lines of code, 250 of tests, 60 of documentation.

## Acceptance criteria

- [ ] The decision and its evidence recorded: the flow maps at step 0 (ε
      and Barnes et al. 2014b), the river maps after 1,000 steps (ε,
      Barnes, roughness), the Hack exponents above.
- [ ] `roughen` by hand and its panics: ocean cells unchanged bit for bit,
      land cells raised by at least 0 and less than the amplitude, nothing
      changed with 0 erosion steps.
- [ ] The flow over a synthetic flat converges with the roughness: on the
      40,000-cell cap, roughened from the streams of seeds 1–5, the
      largest basin at least twice that of the bare cap at step 0 for
      every stream, and after 10 steps at least twice for their median
      and 1.5 times for each (review round 1, R1-2); no pit, every land
      cell reaches the ocean (the M3-02 invariants).
- [ ] The Earth worlds of seeds 1–5 at 650,000 cells: Hack exponent from
      0.55 to 0.65 (Earth 0.6; 0.51–0.52 before); the close-ups show
      dendritic networks without the parallel pinnate tributaries, viewed.
- [ ] Golden values of the Erosion stage and of `kuniumi generate`
      updated, stated in the commit's `World-breaking:` line and in
      `CHANGELOG.md`; those of M3-01 to M3-04 unchanged.
- [ ] `kuniumi generate` for seed 1 at 650,000 cells within 5 % of the time
      of the same run without the roughness on this machine.
- [ ] `mise run ci` passes.

## Verification

```sh
mise exec -- cargo test -p kuniumi-sim erosion
mise exec -- cargo run --release -p kuniumi-cli -- generate --seed 1 --out out/M3-06/650k/seed1
```

## Result

The owner chose the initial roughness in this ticket (no flat
resolution), with the amplitude a constant 1 m, the stream
`erosion/roughness`, nothing with 0 erosion steps, the design document
updated and the golden values recorded anew (world-breaking).

What was done:
- `kuniumi-sim`, `stages::erosion`: `ROUGHNESS_STREAM`,
  `INITIAL_ROUGHNESS_M` and `roughen` as planned; `run` roughens the Crust
  stage's heights before `erode` when there are steps. The module
  documentation, that of `run` and of `roughen` (the decision and its
  sources: Barnes, Lehman & Mulla 2014b, §2.1; Cordonnier, Bovy & Braun
  2019, Fig. 1 and §3.3, checked in the paper), and the peak heights quoted
  there and in `ErosionSettings::DEFAULT_ERODIBILITY` (`kuniumi-world`).
- Tests: `stages/erosion/roughness_tests.rs` (4 tests); in
  `stage_tests.rs` the expected stage output erodes the roughened surface,
  and the stage's heights are checked to differ from the erosion of the
  bare surface; the golden values of the stage and of the CLI recorded
  anew.
- `docs/design-ja.md` (侵食と地形発達 / 1ステップの処理: the initial
  roughness and why; Barnes et al. 2014b in 参考文献), `docs/roadmap.md`
  (M3's tickets), `CHANGELOG.md` (Added, Changed, World-breaking).

Deviations from the plan:
- "Raised by less than the amplitude" holds only before rounding: `h +
  A u` with `u < 1` can round to `h + A`, or a hair above it, where a ULP
  of `h` exceeds `A (1 − u)`. The documentation says "less than the
  amplitude but for the rounding of the sum", and
  `roughness_of_the_earth_worlds` allows 1e-9 m above it.
- `roughen` raises a land cell only where `A u_i > 0`, so that amplitude 0
  (and a draw of 0) keeps the height's bits, the sign of a zero included;
  `roughen_by_hand` checks it with a land cell at −0.
- The cap test draws the roughness from five streams, not one (review
  round 1, R1-2), and criterion 3 says so.
- The roughness tests are in a new file, `roughness_tests.rs`;
  `kuniumi-world (docs)` added to the `Crates:` line for the peaks quoted
  in its documentation.

Evidence:

1. Criterion 1, the decision: "Decision" in the Implementation notes
   above, with the images in `out/M3-06/` (`flats/step0_*_flow.png`,
   `closeups_{eps,barnes,eps_1,eps_10}/`, `closeups_s3_*/`), all viewed.
2. Criterion 2 (`roughness_tests.rs`): `roughen_by_hand` (seven heights:
   ocean below and at the sea level of 0 m with both zeros, land from
   5e-324 m to 1,000 m; amplitudes 1, 2.5 and 1e-3 m; the expected `h + A
   u_i` with the i-th draw, one draw per cell, the generator advanced once
   per cell; a sea level of −698.5 m; amplitudes 0 and −0 unchanged bit
   for bit, a land cell at −0 raised to its draw);
   `roughen_needs_valid_inputs` (amplitudes −1, −5e-324, NaN, ±∞; sea
   levels and heights NaN, ±∞; a cell at `f64::MAX` raised by up to
   `f64::MAX`); `roughness_of_the_earth_worlds` (seeds 1–3 and the dry
   world at 10,000 cells: the ocean unchanged bit for bit, every land cell
   raised by 0 to 1 m and still land; mean rises 0.4991, 0.4981, 0.4969
   and 0.5017 m, shares below 0.1 m 0.100, 0.094, 0.107 and 0.101, from
   0.9 m 0.105, 0.099, 0.105 and 0.100). Nothing with 0 steps:
   `zero_steps_keep_the_surface_of_the_crust_stage` and the CLI's
   `generate_without_erosion_steps_keeps_the_relief_of_the_crust_stage`
   (the Crust stage's golden value), unchanged and passing.
   `the_stage_erodes_and_routes_the_final_surface` and
   `the_stage_uses_every_erosion_setting`: the stage equals `erode` of the
   surface roughened from `Rng::for_stream(seed, "erosion/roughness")`
   with 1 m, bit for bit, and differs from `erode` of the bare surface.
3. Criterion 3 (`a_rough_flat_drains_through_converging_rivers`): the cap
   of 4,677 land cells; largest basin 219 cells bare before erosion and
   388 after 10 steps. Roughened from the streams of seeds 1–5: 1,506,
   1,348, 1,569, 1,311 and 850 cells before erosion (3.88–7.16 times),
   1,104, 1,248, 1,332, 1,208 and 844 after (2.18–3.43 times, median
   3.11); seed 1 gives the prototype's numbers. No pit and every base
   level an ocean cell on all twelve surfaces.
4. Criterion 4 (`kuniumi generate --seed N --out out/M3-06/650k/seedN`,
   release, default settings; `metrics.txt`):

   | Seed | `hack_exponent` (M3-08 before) | `hack_basins` | highest (m) |
   | --- | --- | --- | --- |
   | 1 | 0.6110 (0.5247) | 10,614 | 14,103.7 |
   | 2 | 0.5819 (0.5094) | 11,020 | 8,034.0 |
   | 3 | 0.6126 (0.5093) | 10,115 | 12,642.4 |
   | 4 | 0.5982 (0.5241) | 9,811 | 5,606.0 |
   | 5 | 0.5892 (0.5213) | 10,556 | 8,885.2 |

   The same exponents as the prototype's roughness of 1 m (item 5 of the
   Decision). Viewed: `seed1/rivers.png` and `seed5/rivers.png`
   (dendritic networks on every continent, trunks to the coasts, none
   across the belts); `closeups_s1/east_rivers.png` (4° N, 130° E: trunks
   that wander, tributaries in every direction, none of the pinnate
   parallel tributaries of M3-05), `closeups_s1/interior_rivers.png`
   (31° N, 55° W: dendritic), `closeups_s1/plateau_both.png` (35° S,
   130° W: the belt's transverse valleys as before, the plains beside it
   dendritic instead of the long straight rivers);
   `closeups_s3/southpole_rivers.png` (seed 3: dendritic, no parallel
   coastal lines). Land within 1 m of the sea level: 180,666 of 195,865
   (seed 1; 180,665 before).
5. Criterion 5: golden values recorded after 2–4 and the images:
   the stage (`golden_fingerprints_of_the_erosion_stage`) `surface_height`
   `0xFDE1_B5E3_689C_E2AE`, `drainage_area` `0x381B_B3AE_6288_1477`,
   `discharge` `0x423E_19A2_C079_1331`, `depression_depth`
   `0xBF73_131B_22D5_2440`; the CLI (`generate_golden_fingerprints`)
   `0xDB1E_0D68_718E_F38C`, `0x74B7_97DE_863F_9C2A`,
   `0x714B_06C0_204A_48F6`, `0x2128_8731_7CB9_00EB`; the twelve earlier
   fields, the sea level, the Crust stage's height and the golden values
   of M3-01 to M3-04 unchanged; the same in debug and release.
6. Criterion 6, timing (release, this machine, seed 1 at 650,000 cells,
   one run at a time, `/usr/bin/time`): 112.81 and 113.25 s with the
   roughness, 113.19 and 113.67 s with `main`'s binary; 312 MB each. The
   outputs equal those of the runs above and of M3-05 respectively.
7. Criterion 7: `mise run ci` passes (fmt, clippy `-D warnings`, rustdoc,
   the browser checks, the tests in debug and release, MSRV 1.85.0,
   cargo-deny).
8. Test-first: against a stub that changed nothing, 5 of the new and
   changed tests failed (`roughen_by_hand`, `roughen_needs_valid_inputs`,
   `roughness_of_the_earth_worlds`,
   `a_rough_flat_drains_through_converging_rivers`,
   `the_stage_erodes_and_routes_the_final_surface`). Mutants (release),
   each caught: the ocean raised (4 tests), a cell at the sea level raised
   (3), no guard on a zero rise (1), draws for land cells only (2), the
   stage without roughness (4), roughness with 0 steps (1), another stream
   name (4), amplitude 2 m (5), no overflow check (1), amplitude 0 refused
   (1), half the draw (3).

Review:
- Round 1 (APPROVE): the reviewer ran `mise run ci`, re-ran the prototype
  (every number of the Decision reproduced), checked the prototype's
  flat resolution against Algorithms 4–7 of the paper, regenerated seeds
  1–5 at 650,000 cells (every `.kfield` byte-identical, the same Hack
  exponents), viewed the maps and close-ups of its own (also two sites of
  seed 4 of its choosing), re-timed seed 1 against `main` (112.96 and
  113.28 s against 113.84 and 113.32 s) and caught three mutants of its
  own. Two minor findings and three nits, all fixed: R1-1 (minor), the
  citation of Cordonnier et al. (2019) gave §4 instead of §3.3, and the
  documentation and the design document presented their setup as general
  practice: §3.3, worded as their runs; R1-2 (minor), the cap test's
  threshold after 10 steps rested on one realisation (1.64 times with the
  draws `1 − u`): five streams, median and floor; R1-3 (nit), an
  unneeded `#![allow(clippy::disallowed_methods)]` removed; R1-4 (nit),
  `run`'s "Panics" names `roughen`; R1-5 (nit), the Priority-Flood paper
  labelled 2014a in the design document.
- Round 2 (APPROVE, no findings): the reviewer ran `mise run ci`,
  confirmed R1-1 to R1-5 resolved (the cap test also passes with the
  draws `1 − u` that failed the old one: 1.64–3.81 times, median 3.11),
  re-ran the roughness tests (the printed numbers of items 2 and 3),
  regenerated seed 1 at 650,000 cells (every `.kfield` byte-identical,
  `hack_exponent` 0.6110), viewed its `rivers.png` and timed it (112.93 s,
  312 MB).

Follow-ups: none.

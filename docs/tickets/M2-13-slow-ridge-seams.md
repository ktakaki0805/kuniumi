# M2-13 Age steps around very slow ridge cells

- Status: done
- Milestone: M2
- Depends on: M2-11
- Crates: kuniumi-sim (stages/tectonics), kuniumi-cli (tests)
- Design: 地殻から標高へ / 海洋底：年代から水深; プレートテクトニクス / 境界ごとの処理
- Decisions: none

## Goal

Decide whether the remaining steps of ocean-floor age next to very slow
ridge cells are smoothed, so that the sea floor of M2-05 has no step of
kilometres beside a ridge.

## Scope

In:
- The finding of the review of M2-11: the time-nearest growth of M2-11
  removes most age seams (pairs of neighbouring oceanic cells of one plate
  more than 20 Myr apart: 82 at 100,000 cells, seed 1, against 1,955 with
  M2-04), but a ridge cell that spreads very slowly is still the first to
  reach its own neighbours when no faster front is near, as in clusters
  of slow cells near a relative pole. At 100,000 cells, seed 1, ridge cell
  39194 (0.096 cm/yr, 13° N 80° W) has a neighbour of ranked age 105.9 Myr,
  a step of about 3.2 km in M2-05's sea floor; 20 of the 82 seams touch a
  ridge cell, and the largest step is 106 Myr. At 650,000 cells, seeds
  1–10: no seam in six worlds; 3–110 seams in the others (seed 5: 110, 15
  next to a ridge cell), with largest steps of 25–75 Myr.
- Options to decide with the owner: a floor on the half spreading rate of
  a ridge cell used for the times (for example the plate's mean rate, or a
  fixed minimum); rates smoothed along each ridge; or keep the steps
  until M6.

Out:
- The time-evolving tectonics of M6.
- The field `convergence_rate` and the boundary types of M2-03: a floor
  applies only to the rate with which `ridge_times` grows its fronts.
- Steps that are not next to a slow ridge: at 10,000 cells (225 km cells)
  a front at 2 cm/yr takes 22 Myr per cell, so neighbours more than 20 Myr
  apart are the grid's resolution there, not an artefact.

## Implementation notes

Refined from `draft`. The options were measured with a throwaway
prototype in the scratchpad (not committed) on the real grids, plates,
continents and boundaries (Earth's parameters, 12 plates): it changes the
convergence rate of the divergent cells in a copy of `Boundaries` and runs
the committed `crust`, `surface_heights_m`, `sea_level_m` and
`land_fraction` on it, so each variant is exactly what `crust` would give
with those rates.

### Why the steps remain

A ridge cell's front needs `2 d R / v` to cross one cell (`v` the full
spreading rate `−c`, so `v/2` the half rate; `d R` the centre distance,
about 77 km at 100,000 cells and 30 km at 650,000). Near the relative
pole of two plates the rate falls to zero, and a cell next to a ridge cell
spreading at 0.1 cm/yr is first reached by that ridge's own front after
about 150 Myr when no faster front gets there earlier; ranked onto Earth's
distribution that is a step of about 100 Myr from the ridge's +0. Earth's
slowest ridges spread far faster than that: the ultraslow ridges, the
Gakkel and Southwest Indian ridges, at full rates below about 12 mm/yr,
their characteristics common up to about 20 mm/yr (Dick, Lin & Schouten
2003, "An ultraslow-spreading class of ocean ridge", Nature 426,
405–412); the Gakkel Ridge, the slowest, from about 13 mm/yr in the west
to about 6 mm/yr in the east (as summarised by O'Connor et al. 2021,
Nature Communications, doi:10.1038/s41467-021-27058-1).

### Prototype measurements

Seams (pairs of neighbouring oceanic cells of one plate whose ages differ
by more than 20 Myr, all of them and those touching a ridge cell) and the
largest step next to a ridge cell (Myr), on the current `main` (after
M2-10 and M2-12, so 81 seams at 100,000 cells, seed 1, where M2-11
counted 82). "Floor v" uses `max(−c, v)` as the full rate of each ridge
cell; "pair mean × f" uses `max(−c, f · mean)`, the mean over the
divergent cells of the same pair of plates; "smooth W" replaces each ridge
cell's rate with the mean over the divergent cells of its pair within
`W` along the boundary band.

| Option | 100k seed 1: seams (at a ridge) | step at a ridge | 100k seeds 2–5: seams | 650k seeds 1–10: seams | step at a ridge, 650k |
| --- | --- | --- | --- | --- | --- |
| as now | 81 (20) | 105.9 | 8, 13, 78, 38 | 0, 24, 8, 0, 110, 0, 0, 3, 0, 0 | 1.7–57.1 |
| floor 0.6 cm/yr | 13 (0) | 17.1 | 8, 0, 16, 8 | 0 in all | 1.7–7.7 |
| **floor 1 cm/yr** | **0 (0)** | **11.5** | **8, 0, 0, 0** | **0 in all** | **1.6–4.9** |
| floor 2 cm/yr | 0 (0) | 6.2 | 0, 0, 0, 0 | 0 in all | 1.1–3.0 |
| pair mean × 0.5 | 0 (0) | 13.1 | — | 0, except seed 2: 24 | 1.7–12.7 |
| smooth 1,000 km | 20 (2) | 22.0 | — | 0, except seed 2: 14 | 1.3–13.3 |
| smooth 500 km | 75 (9) | 36.4 | — | 0, except seeds 2: 15, 5: 46 | 1.5–16.6 |

Other effects of the fixed floors:

- Ridge cells raised: 1 cm/yr raises 80 of 1,561 at 100,000 cells, seed 1
  (0.1–7.7 % in seeds 1–5) and 0–8.2 % at 650,000 cells (seed 9: none);
  2 cm/yr raises 19–39 % (100,000 cells) and 0.6–23 % (650,000 cells).
- Ages (mean 60 Myr by the ranking, unchanged): at 100,000 cells, seed 1,
  the floor of 1 cm/yr changes 21 % of the ocean floor by more than 1 Myr
  and 2.3 % by more than 5 Myr (0.6 cm/yr: 3.8 % and 0.5 %; 2 cm/yr: 77 %
  and 32 %, since it changes the times along a quarter of the ridges); the
  largest change, 147 Myr, is the floor beside the slow ridge cells, which
  turns from about 100 Myr to a few Myr.
- Land fraction unchanged to 0.01 points in every world measured (100,000
  cells, seeds 1–5: 28.14–29.94 %; 650,000 cells, seeds 1–10:
  28.55–30.84 %), sea level within 1 m; the number of trench cells
  changes by at most 4 (1 cm/yr) or 5 (2 cm/yr) of 300–2,000.
- The eight seams that the floors of 0.6 and 1 cm/yr leave at 100,000
  cells, seed 2 (20–22 Myr, at 9–51 Myr of age, seven at 35–50° S
  130–145° W and one at 32° N 148° E) touch no ridge cell; the slowest ridge cell there spreads at 0.99 cm/yr,
  and fronts at 1.0–1.3 cm/yr take 12–15 Myr per cell at that resolution,
  so these are the steep but real age gradients of slow spreading, not
  the artefact of this ticket.

Smoothing along the ridges keeps the clusters of slow cells near a
relative pole slow (their neighbours along the ridge are slow too), so it
removes few seams. A floor from each pair's mean rate depends on the shape
of the whole ridge and leaves the 24 seams of seed 2 at 650,000 cells,
whose pair spreads slowly everywhere.

### Plan (recommended option, a fixed floor of 1 cm/yr)

World-breaking: `crust_age` changes wherever a ridge cell spreads more
slowly than the floor (also in the golden world, 9 of 413 ridge cells at
10,000 cells, seed 1), and with it `surface_height`, the sea level and
perhaps `trench`. Nothing else changes.

1. `stages::tectonics::crust`: a constant
   `pub const MIN_SPREADING_RATE_M_S: f64 = 0.01 / SECONDS_PER_YEAR;`
   (1 cm/yr full rate, 0.5 cm/yr half rate), documented with the sources
   above as an order-of-magnitude choice made with the owner: "the rate of
   Earth's ultraslow ridges; a ridge cell spreading more slowly, as near
   the relative pole of its plates, grows its front as if it spread at
   this rate". Re-exported from `stages::tectonics`.
2. `ridge_times`: the front of ridge cell `s` takes the time
   `2ℓ / max(−c_s, v_min)` (by the explicit comparison `larger`, not
   `f64::max`). Doc comments of `ridge_times`, `RidgeTimes::time_s` (no
   more "or if a rate underflows"), `earth_area_age` and `crust` (step 1)
   updated; the trench rule keeps comparing these times.
3. Tests (`crust/age_tests.rs`, `crust/tests.rs`), written first:
   - `best_offer` (the reference of `ridge_times`) uses
     `max(−c_s, 1 cm/yr)`, with 1 cm/yr written in the test, so
     `ridge_times_match_a_reference_search` checks the floor on all its
     worlds;
   - `ridge_times_on_the_octahedron`: a ridge at 0.2 cm/yr reaches its
     neighbours at `2q / 1 cm/yr` (by hand), and one at exactly 1 cm/yr
     and one faster keep their own rates;
   - `slow_ridges_reach_their_neighbours_within_a_step` (new): on the
     two-plate cases, seeds 1–3 at 10,000 cells and seed 1 at 100,000
     cells, every cell `j` of the plate next to a ridge cell `s` has
     `t_j ≤ 2 d_sj R / max(−c_s, 1 cm/yr)` (1e-12 relative), and at least
     one world has a ridge cell slower than 1 cm/yr;
   - `seams_almost_vanish` becomes `no_seams_at_100000_cells`: at 100,000
     cells, seed 1, no pair of neighbouring oceanic cells of one plate
     more than 20 Myr apart (81 before), and every step next to a ridge
     cell below 20 Myr (105.9 Myr before).
4. Golden values of `kuniumi-cli` (`tests/cli.rs`): `crust_age`,
   `surface_height`, the sea level and the reference metrics of
   `tools/hypsometry_reference.py` (and `trench` if it changes),
   re-recorded after the tests above pass and the maps are viewed.
5. `docs/design-ja.md` (地殻から標高へ / 海洋底): the sentence on the steps
   that remain replaced by the floor, with the owner's agreement;
   `CHANGELOG.md` (Changed, World-breaking).

Measured outside the tests with a scratch tool on the committed code: the
seams and the largest step next to a ridge cell at 100,000 cells (seeds
1–5) and 650,000 cells (seeds 1–10), the land fraction and sea level, and
the timing of `kuniumi generate` at 650,000 cells.

### Question for the owner

How should the age steps next to very slow ridge cells be handled?

- Recommended: **a fixed floor of 1 cm/yr** on the full spreading rate
  (0.5 cm/yr half rate) used for the times, as planned above: 0 seams at
  100,000 cells (seed 1) and at 650,000 cells (seeds 1–10), steps next to
  a ridge at most 11.5 Myr (100,000) and 4.9 Myr (650,000); land fraction
  unchanged.
- A floor of 0.6 cm/yr (Earth's slowest observed spreading, the eastern
  Gakkel Ridge): 13 seams remain at 100,000 cells, seed 1 (none at a
  ridge, the largest 24 Myr); none at 650,000 cells.
- A floor of 2 cm/yr (the upper end of ultraslow spreading): no seams
  anywhere measured, but it changes the rates of 10–40 % of the ridge
  cells and, at 100,000 cells, seed 1, a third of the ocean floor by more
  than 5 Myr.
- Keep the steps until M6: no change; this ticket records the
  measurements and closes.

The bound for the first criterion follows the choice; for the
recommended option it is the one written below.

### Decision (owner, 2026-11-12)

The owner chose the recommended option, **a fixed floor of 1 cm/yr** on
the full spreading rate used for the times, approving the plan above: the
bounds of the second criterion (100,000 cells, seed 1: no seam and every
step next to a ridge cell below 20 Myr; 650,000 cells, seeds 1–10: no
seam and every such step below 10 Myr), the world-breaking change, and
the update of the design document (地殻から標高へ / 海洋底) in the same
commit.

## Acceptance criteria

- [x] The floor: `MIN_SPREADING_RATE_M_S` is 1 cm/yr; on the octahedron a
      ridge slower than the floor reaches its neighbours at
      `2q / 1 cm/yr` and faster ridges keep their rates (by hand);
      `ridge_times` equals the reference search with the floored rate on
      the two-plate cases and on real worlds (seeds 1–3 at 10,000 cells,
      seed 1 at 100,000) bit for bit; every neighbour of a ridge cell on
      its plate is reached within `2 d R / max(−c, 1 cm/yr)` (1e-12
      relative).
- [x] Seams and the largest step next to a ridge cell: at 100,000 cells,
      seed 1, no seam (81 before) and every step next to a ridge cell
      below 20 Myr (105.9 Myr before), in a test; at 650,000 cells, seeds
      1–10, no seam and every step next to a ridge cell below 10 Myr
      (scratch tool); also reported for 100,000 cells, seeds 2–5.
- [x] Land fraction with Earth's parameters stays within 0.1 points of
      `main` at 100,000 cells (seeds 1–5) and 650,000 cells (seeds 1–10),
      and the existing land-fraction tests at 10,000 cells pass.
- [x] Maps of `crust_age` and `surface_height` of the affected seeds
      (100,000 cells, seed 1, near 13° N 80° W; 650,000 cells, seeds 3
      and 5) viewed and described, before and after.
- [x] Golden values re-recorded after the tests above pass and the maps
      are viewed; world-breaking line and `CHANGELOG.md` entry; design
      document updated.
- [x] `mise run ci` passes.

## Verification

```sh
cargo run --release -p kuniumi-cli -- generate --seed 1 --cells 100000 --out out/M2-13/100k
cargo run --release -p kuniumi-cli -- generate --seed 5 --out out/M2-13/seed5
cargo run --release -p kuniumi-cli -- generate --seed 3 --out out/M2-13/seed3
```

## Result

What was done:
- `kuniumi-sim`, `stages::tectonics::crust`: the constant
  `MIN_SPREADING_RATE_M_S` (1 cm/yr, with Dick, Lin & Schouten 2003 and
  O'Connor et al. 2021 as its order-of-magnitude sources), re-exported from
  `stages::tectonics`; `ridge_times` grows the front of ridge cell `s` at
  `larger(−c_s, MIN_SPREADING_RATE_M_S)`. Documentation of the module,
  `RidgeTimes::time_s`, `ridge_times`, `earth_area_age` and `crust`.
- Tests: `crust/tests.rs` (`SPREADING_FLOOR_M_S` and `front_rate_m_s`,
  written independently of the constant, used by `best_offer` and so by
  `reference_times` and every test built on them; the constant in
  `the_field_names_and_constants`); `crust/age_tests.rs` (the floor cases
  of `ridge_times_on_the_octahedron`, the floored rate in
  `ridge_times_match_a_reference_search`, the new
  `slow_ridges_reach_their_neighbours_within_a_step`, and
  `seams_almost_vanish` replaced by `no_seams_at_100000_cells`);
  `kuniumi-cli` tests: the golden values and the reference metrics.
- `docs/design-ja.md` (地殻から標高へ / 海洋底, the two sources in the
  references), `CHANGELOG.md` (Changed, World-breaking).

Deviations from the plan:
- Trench cells: the plan's measurements counted trench cells only; cell
  by cell, the side that sinks where ocean meets ocean changes in 12, 15
  and 2 cells at 100,000 cells (seeds 1, 3 and 5) and in 112, 72 and 66
  cells in seeds 3, 5 and 6 at 650,000 cells (none in the other seeds),
  near the slow ridges whose times drop; none in the golden world (its
  `trench` fingerprint is unchanged). Measured cell by cell against
  `main`'s binary, and by the reviewer in round 1.
- Worlds that change: the plan said `crust_age` changes wherever a ridge
  cell spreads more slowly than the floor; it changes only where such a
  ridge cell is the first to reach one of its neighbours. At 650,000 cells
  seeds 1–8 change (165,000–384,000 cells each, since the ranking moves
  nearly all floor a little); seed 9 has no ridge cell below the floor,
  and in seed 10 faster fronts reach every neighbour of its 102 slow ridge
  cells first, so both are bit for bit as on `main` (found in review
  round 1).

Evidence:

1. The floor: `the_field_names_and_constants` (1 cm/yr = 3.1688e-10 m/s);
   `ridge_times_on_the_octahedron` (by hand: a ridge at 0.2 cm/yr and one
   at exactly 1 cm/yr both give `2q / 1 cm/yr` and `2 · 2q / 1 cm/yr`; at
   0.2 cm/yr against 1.5 cm/yr the faster keeps its own rate and wins the
   four other cells at `2q / 1.5 cm/yr`; at 0.2 against 0.5 cm/yr both
   grow at the floor and the lower ridge wins the ties, where the faster
   one would win at its own rate); `ridge_times_match_a_reference_search`
   (the 32 two-plate cases, seeds 1–3 at 10,000 cells and seed 1 at
   100,000: times, paths and ridges bit for bit against `reference_times`
   with `max(−c_s, 1 cm/yr)`, every time `2ℓ/max(−c_s, 1 cm/yr)` bit for
   bit, and every cell the least offer of its neighbours' labels);
   `slow_ridges_reach_their_neighbours_within_a_step` (same worlds: every
   neighbour of a ridge cell on its plate reached by
   `2 d R / max(−c, 1 cm/yr)` within 1e-12; 314 of the 12,026 ridge cells
   spread more slowly than the floor). Test-first: against a stub (the
   constant defined, `ridge_times` unchanged) 7 of the 29 crust tests
   failed (commit "tests of the spreading floor against a stub"), among
   them the three above and `no_seams_at_100000_cells`.
2. Seams (pairs of neighbouring oceanic cells of one plate more than
   20 Myr apart) and steps next to a ridge cell: `no_seams_at_100000_cells`
   (100,000 cells, seed 1): 0 pairs (81 on `main`), largest step 15.1 Myr,
   11.5 Myr next to a ridge cell (105.9 Myr on `main`). Scratch tool on the
   committed code (the real `crust` stage functions), 650,000 cells, seeds
   1–10: no pair in any world (0, 24, 8, 0, 110, 0, 0, 3, 0, 0 on `main`);
   largest step next to a ridge cell 4.4, 3.7, 2.9, 1.6, 3.0, 2.6, 3.9, 4.9,
   2.8 and 1.7 Myr (8.0, 12.7, 50.6, 3.4, 57.1, 3.6, 7.1, 25.5, 2.8, 1.7 on
   `main`), largest step anywhere 2.2–9.1 Myr. 100,000 cells, seeds 2–5:
   8, 0, 0, 0 pairs (8, 13, 78, 38 on `main`), largest step next to a
   ridge cell 17.8, 8.4, 7.3, 7.7 Myr; the eight pairs of seed 2 touch no
   ridge cell and are the real gradients of ridges at 1.0–1.3 cm/yr
   (Implementation notes). The numbers equal those of the prototype's
   "floor 1 cm/yr" in every world, at the precision printed.
3. Land fraction (scratch tool, M2-05's stage): 100,000 cells, seeds 1–5,
   29.62, 29.94, 28.14, 29.52, 29.74 %; 650,000 cells, seeds 1–10, 30.13,
   30.76, 28.95, 28.71, 29.97, 29.09, 29.24, 30.74, 30.84, 28.55 %: equal
   to `main` at 0.01 points in every world; sea level within 1 m of
   `main` (−716.63 → −716.78 m at 100,000 cells, seed 1; −722.59 →
   −722.67 m and −701.00 → −701.03 m in seeds 3 and 5 at 650,000 cells). At
   10,000 cells (seeds 1–10) both land-fraction tests pass with the same
   means as after M2-12: `stages::crust::tests::land_fraction_of_earth_like_worlds`
   28.14–30.04 %, mean 28.90 %; `margin_tests::land_fraction_of_earth_like_worlds`
   27.8–29.6 %, mean 28.58 %.
4. Maps (`kuniumi generate`, this branch into `out/M2-13/100k`, `seed3`,
   `seed5` and `10k`, `main`'s binary into `out/M2-13/*-main`; crops
   magnified three times in `out/M2-13/crop-*.png`), viewed:
   - 100,000 cells, seed 1, `crust_age` near 13° N 80° W: on `main` the
     dark young band of the meridional ridge between the western plates
     ends abruptly at the plate's northern tip against a yellow, nearly
     180 Myr tongue of ocean floor (the step of 106 Myr); on the branch
     the dark band runs on to the tip and the floor west of it turns
     smoothly from dark through teal to green, with no tongue. The sharp
     edges left round the tip are plate boundaries, with other plates'
     floor across. `surface_height`: on `main` the pale ridge line ends
     against the deep floor of the tongue; on the branch the ridge fades
     into the deep floor over several cells. The rest of both maps is the
     same to the eye.
   - 650,000 cells, seed 5, near 1° N 127° E: on `main` a dark ridge
     ending in a hook beside a pale yellow-green wedge (57 Myr steps); on
     the branch the ridge continues as a smooth dark band, the wedge gone
     and replaced by a teal halo; `surface_height` shows the light ridge
     line widening smoothly into the floor instead of ending against deep
     floor. The trench line next to it (north-east) is the same.
   - 650,000 cells, seed 3, near 6° N 109° W: the young band along a
     coast-parallel ridge stops abruptly on `main` where the ridge slows,
     with older teal floor right beside the ridge cells; on the branch the
     young band continues along the coast and fades gradually;
     `surface_height` slightly shallower there, the arc and trench along
     the convergent coast to the north-east unchanged.
   - 10,000 cells, seed 1 (golden world), `crust_age`: dark ridges, ages
     rising away from them, the ridgeless plate north of the central
     continent light green as after M2-12; nothing unexpected.
5. Golden values (`tests/cli.rs`, 10,000 cells, seed 1), recorded after
   the tests above passed and the maps were looked at: `crust_age`
   `0x8B16_21E7_DFBF_D241`, `surface_height` `0xA90E_8D61_5D3A_E286`, sea
   level −747.218 732 917 737 4 m (was −747.142 m), the reference metrics
   from `tools/hypsometry_reference.py out/M2-13/10k` (land fraction
   0.291 092, unchanged); the other eleven fingerprints unchanged, `trench`
   among them. World-breaking line, `CHANGELOG.md` entry and design
   document.
6. `mise run ci` passes: fmt, clippy `-D warnings`, rustdoc `-D warnings`,
   the tests in debug and release, MSRV 1.85.0, cargo-deny (advisories,
   bans, licences, sources ok), before review and again after the fixes of
   round 1.

Timing (release, this machine, alternating runs): `kuniumi generate
--seed 1` (650,000 cells) 6.4–6.8 s and 208 MB on this branch, 6.1–6.6 s
and 208 MB with `main`'s binary measured the same way (the difference is
within the run-to-run spread; the reviewer measured 5.45 s against
5.34 s, 208 MB both).

Review:
- Round 1 (APPROVE): the reviewer re-ran `mise run ci`, the stub commit
  (7 of 29 crust tests failing), mutants of the floor (0.5 cm/yr fails 7
  tests, `no_seams_at_100000_cells` among them; 2 cm/yr fails 5), the
  seams, steps and land fractions at 100,000 cells (seeds 1–5) and
  650,000 cells (seeds 1–10) with its own tool over the whole pipeline on
  the branch and on `main` (the same numbers as evidence 2 and 3), the
  maps and crops (`out/review/M2-13/`), the golden world (only
  `crust_age` and `surface_height` differ from `main`) and the reference
  metrics. R1-1 minor: the World-breaking entry said every world with a
  ridge cell below the floor changes, but seed 10 at 650,000 cells does
  not, and the trench list left out seed 6 (650,000 cells) and seeds 3
  and 5 (100,000 cells) (fixed: `CHANGELOG.md` and the Deviations above,
  checked cell by cell against `main`'s binary). R1-2 nit: the `crust`
  doc did not name M2-13 (fixed). R1-3 nit: the bound on the step next to
  a ridge cell in `no_seams_at_100000_cells` added little to `seams == 0`
  (fixed: 15 Myr, against 11.5 Myr measured).

Follow-ups: none.

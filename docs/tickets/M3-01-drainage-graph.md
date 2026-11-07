# M3-01 Drainage graph: receivers, stack, accumulation

- Status: ready
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

## Acceptance criteria

- [ ] No cycles: following receivers from any cell reaches an outlet or a pit
      within n steps.
- [ ] The stack contains every cell exactly once, each after its receiver.
- [ ] Accumulated area summed over outlets and pits equals the total area within
      1e-12 relative.
- [ ] Synthetic surfaces: a tilted plane drains one way, a cone drains radially,
      a bowl produces a single pit.
- [ ] O(n) time; 650,000 cells in under 1 s in release (time recorded).
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-sim drainage
```

## Result

_(filled in when the ticket is implemented)_

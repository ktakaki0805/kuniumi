# M5-02 Lakes by water balance

- Status: draft
- Milestone: M5
- Depends on: M5-01
- Crates: kuniumi-sim (stages/hydrology)
- Design: 水文（川・湖）
- Decisions: none

## Goal

Depressions fill with water according to their water balance: lakes that spill
continue as rivers (fresh), lakes that balance by evaporation before spilling
become terminal salt lakes.

## Scope

In:
- Depression hierarchy from the Priority-Flood of M3-02 (which depression spills
  into which, at what height).
- Fill–Spill–Merge (Barnes et al. 2021): process the hierarchy from the leaves;
  each depression receives inflow plus precipitation on the lake surface and
  loses lake evaporation (≈ PET over open water); find the level where it
  balances or spills, and pass the surplus downstream.
- Outputs: `lake_depth` (m), `lake_kind` (class field: none, fresh, salt),
  `water_surface_height` (m); discharge updated to flow through lakes.

Out:
- Wetlands (M5-05 decides them from slope and runoff).

## Implementation notes

To be refined: data structure of the hierarchy, handling of depressions that
touch the ocean, volume–level relation from the cell areas.

## Acceptance criteria

- [ ] Water is conserved through the hierarchy (inflow = outflow + evaporation
      for every lake, 1e-9 relative).
- [ ] Synthetic cases: a humid bowl fills and spills; an arid bowl forms a salt
      lake below its rim; nested bowls merge when the inner one overflows.
- [ ] Lake maps for a generated world viewed and described (endorheic lakes
      appear in dry interiors).
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-sim lakes
```

## Result

_(filled in when the ticket is implemented)_

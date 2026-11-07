# M6-02 Boundary detection by overlap and gaps

- Status: draft
- Milestone: M6
- Depends on: M6-01
- Crates: kuniumi-sim (stages/tectonics)
- Design: プレートテクトニクス / 表現; プレートテクトニクス / 境界ごとの処理
- Decisions: none

## Goal

Each step, find where plates overlap (convergent boundaries) and where no plate
covers the surface (divergent gaps), and classify each convergent contact.

## Scope

In:
- Per grid cell after rotation: the set of plates whose samples cover it.
- Overlaps of two plates: convergent; the contact type (ocean–continent,
  ocean–ocean, continent–continent) from the crust types; transform where the
  relative motion is parallel to the boundary.
- Gaps: cells covered by no plate (divergent).
- A per-step boundary record (cells, plates, contact type, relative speed)
  consumed by M6-03 to M6-05.

## Implementation notes

To be refined: handling of triple overlaps, minimum gap size, the transform
angle threshold (shared with M2-03).

## Acceptance criteria

- [ ] Two plates moving towards each other produce an overlap band whose width
      grows with step × speed; moving apart produces a gap of matching width.
- [ ] Contact types are correct in synthetic two-plate tests for each pair.
- [ ] Boundary maps over several steps viewed and described.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-sim boundaries
```

## Result

_(filled in when the ticket is implemented)_

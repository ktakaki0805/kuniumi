# M6-04 Ridge spreading and ocean-floor age

- Status: draft
- Milestone: M6
- Depends on: M6-02
- Crates: kuniumi-sim (stages/tectonics)
- Design: プレートテクトニクス / 境界ごとの処理; 地殻から標高へ / 海洋底：年代から水深
- Decisions: none

## Goal

Gaps between diverging plates fill with new oceanic crust of age zero, so ocean
floor ages and depths emerge from the history instead of a distance formula.

## Scope

In:
- New oceanic samples in gap cells, assigned to the nearer of the two diverging
  plates, with age 0 and standard oceanic thickness.
- Ages advance by the step every step; depth from age with the Parsons &
  Sclater function of M2-05.
- Gaps inside a continent (continental rifting) start as a rift valley before
  becoming ocean (shared with M6-06).

## Implementation notes

To be refined: ridge placement inside wide gaps (midline), handling of
asymmetric spreading.

## Acceptance criteria

- [ ] Two diverging plates produce an age field symmetric about the ridge that
      increases linearly with distance (age ≈ distance / half rate).
- [ ] The ridge lies at about 2.5 km depth and ocean floor deepens away from it.
- [ ] Age maps viewed and described.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-sim spreading
```

## Result

_(filled in when the ticket is implemented)_

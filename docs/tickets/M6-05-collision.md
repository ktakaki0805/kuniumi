# M6-05 Continental collision, suturing and plate merging

- Status: draft
- Milestone: M6
- Depends on: M6-03
- Crates: kuniumi-sim (stages/tectonics)
- Design: プレートテクトニクス / 境界ごとの処理
- Decisions: none

## Goal

When continental crust reaches a subduction front it cannot sink: the continents
collide, the incoming block sutures to the other plate, a wide region thickens
into Himalayan ranges and plateaus, and the plates may merge.

## Scope

In:
- Detection: continental samples of the subducting plate arriving at the front.
- Suturing: the continental block (connected continental samples) transfers to
  the overriding plate; the contact region thickens towards about 70 km over a
  documented width; orogeny type "Himalayan" and age recorded; the suture line
  stored for later rifting (M6-06).
- Plate merging when the remaining relative motion or the contact length passes
  a documented threshold; the merged plate gets a combined Euler rotation.
- Uplift rate recorded for the erosion input.

## Implementation notes

To be refined: the thickening profile, merge criterion, conservation of
continental volume during suturing.

## Acceptance criteria

- [ ] Continental crust volume (thickness × area) is conserved through
      collision within a documented tolerance.
- [ ] A synthetic collision produces a thickened belt about 70 km thick and a
      recorded suture.
- [ ] Maps viewed and described.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-sim collision
```

## Result

_(filled in when the ticket is implemented)_

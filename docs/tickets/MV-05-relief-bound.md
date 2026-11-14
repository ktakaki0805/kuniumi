# MV-05 Relief exaggeration bounded by the deepest point

- Status: draft
- Milestone: MV
- Depends on: MV-02
- Crates: kuniumi-viewer
- Design: 実装ロードマップ / 技術スタック
- Decisions: D-19

## Goal

No relief exaggeration the viewer offers turns a valid world inside out:
the drawn radius `1 + e·h/R` of the deepest point stays clearly above 0 on
every planet within the bounds of D-12.

## Scope

In (review R3-1 of MV-02):
- The exaggeration of `kuniumi-viewer` is capped by its world, not only by
  the fixed ×100 of `RELIEF_STEPS`: the deepest height above the sea level
  is stored beside the highest (`view::Relief`), and the Up key, `--relief`
  and the start view keep `1 + e·h_min/R` at or above a floor (for example
  0.5), or the vertex shader puts a floor on the drawn radius; the choice is
  made when the ticket is refined.
- The doc comment of `MAX_RELIEF` (`args.rs`) gives Earth's trenches as
  7 km deep; Earth's deepest point is about 11 km below the sea level
  (7.2 km is the lowest cell of MV-02's generated world).

Out:
- Everything else of the viewer (MV-03, MV-04).

## Implementation notes

Found by the third review of MV-02: a world of the smallest radius
(1,500 km, Earth's density, Earth's ocean volume, seed 1, 100,000 cells)
has an ocean 46 km deep on average and its lowest cell 50.7 km below the
sea level, so the drawn radius falls below 0 for `e > 29.6`; at ×100 the
inverted sea floor fills the view around the sea sphere, at ×20 it shrinks
to a ball of 0.3–0.5 radii under it. The radius of the deepest point at
×100 measured on 20,000-cell worlds: 0.45 at 3,000 km, 0.78 at 4,500 km,
0.89 at Earth's size, 0.97 at 13,000 km. No crash; Down recovers the view.
The near and far planes come from the highest point only, so they need the
deepest one too if the floor is in the shader.

## Acceptance criteria

- [ ] On the 1,500 km world above and on Earth, the drawn radius of every
      vertex stays at or above the chosen floor for every exaggeration the
      keys and `--relief` reach (test), and a screenshot of the small world
      at the largest exaggeration is saved and described.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo run --release -p kuniumi-viewer -- <small world> --relief 100 --screenshot out/MV-05/small.png
```

## Result

_(filled in when the ticket is implemented)_

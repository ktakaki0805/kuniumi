# M6-03 Subduction, trenches and arcs

- Status: draft
- Milestone: M6
- Depends on: M6-02
- Crates: kuniumi-sim (stages/tectonics)
- Design: プレートテクトニクス / 境界ごとの処理
- Decisions: none

## Goal

At ocean–continent and ocean–ocean overlaps the denser plate sinks: its samples
are removed, a trench forms, and the overriding plate is uplifted into Andean
arcs or island arcs.

## Scope

In:
- Ocean–continent: the oceanic samples subduct (are removed); the overriding
  continental crust is thickened and uplifted as a function of distance from
  the front (profile documented, after Cortial et al. 2019).
- Ocean–ocean: the older (denser) ocean floor subducts; an island arc forms on
  the overriding plate.
- Trench marking on the subducting side; orogeny type and age recorded on the
  uplifted samples ("Andean", "island arc").
- Uplift rate recorded per sample for the final erosion input.

## Implementation notes

To be refined: uplift profile and its parameters, arc volcanism as thickness
added, subduction speed limits.

## Acceptance criteria

- [ ] Synthetic ocean–continent convergence builds a trench and a parallel arc
      at a documented distance; ocean–ocean picks the older plate to subduct.
- [ ] Oceanic area removed equals the overlap area each step.
- [ ] Maps of thickness and orogeny type viewed and described.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-sim subduction
```

## Result

_(filled in when the ticket is implemented)_

# M3-06 Convergent flow across flats

- Status: draft
- Milestone: M3
- Depends on: M3-05
- Crates: kuniumi-sim (drainage)
- Design: 侵食と地形発達 / 1ステップの処理
- Decisions: none

## Goal

Water crossing a flat (a drained flat or a filled depression) converges
into channels as on real low-gradient land, instead of running in the
straight parallel lines that the ε gradient of M3-02 gives.

## Scope

In:
- Decide from the river maps of M3-05 whether the parallel lines of the ε
  routing survive erosion. On the Earth world of seed 1 at 650,000 cells,
  150,162 of 195,865 land cells are drained flats (the uniform +300 m
  continental interiors of M2-05), and `earth_650000_1_plateau_flow.png`
  of M3-02 shows the lines.
- If they do: a flat-resolution step on the filled surface, after Barnes,
  Lehman & Mulla (2014, "An efficient assignment of drainage direction over
  flat surfaces in raster digital elevation models", *Computers &
  Geosciences* 62, 128–135; arXiv:1511.04433; cited by their Priority-Flood
  paper, §7.1, as the alternative to ε): a gradient away from higher
  terrain superimposed on one towards lower terrain, so that flow
  converges.

Out:
- The basin graph of Cordonnier, Bovy & Braun (2019).

## Implementation notes

To be refined once M3-05 has landed: whether the step is needed at all,
its place in the erosion step (on `fill_depressions`' output, before
`Drainage::build`), how it combines with the ε gradient on the global
grid, and its cost at 650,000 cells. Written from the paper alone: the
authors' library, RichDEM, is not a source (docs/conventions.md,
"Dependencies and tools", allows porting only from permissively licensed
code).

## Acceptance criteria

- [ ] The decision and its evidence (M3-05's maps) recorded.
- [ ] If implemented: routing on flats still leaves no pit and reaches the
      ocean (the M3-02 invariants), and the flow over a synthetic flat
      converges (measured on a test surface).
- [ ] `mise run ci` passes.

## Verification

```sh
mise exec -- cargo test -p kuniumi-sim drainage
```

## Result

_(filled in when the ticket is implemented)_

# M2-10 Oceanic pockets enclosed by continents

- Status: draft
- Milestone: M2
- Depends on: M2-03
- Crates: kuniumi-sim (stages/tectonics)
- Design: プレートテクトニクス / 初期条件の目安; 地殻から標高へ
- Decisions: none

## Goal

Decide what becomes of the small oceanic pockets that the continents of
M2-03 enclose; until then M2-04 and M2-05 turn each into ocean floor, a pit
of ocean depth inside a continent.

## Scope

In:
- The continental growth of M2-03 (`stages::tectonics::continents`) stops
  partway through its front, and the slowness of the cells spreads over a
  factor of 20, so a continent can enclose oceanic cells. Measured in the
  review of M2-03 (fraction 0.4, 12 plates, seeds 1–3): oceanic components
  of at most 5 cells number 17–19 at 10,000 cells, 70–81 at 100,000 and
  157–181 at 650,000; 12–15, 57–71 and 116–146 of them are single cells
  whose neighbours are all continental. Without a change, each becomes an
  isolated pit of ocean depth inside a continent in M2-05.
- Options to decide with the owner: fill enclosed oceanic components below
  a size (in cells or in area), keeping the continental fraction within
  M2-03's tolerance (for example by trimming the same area from the
  continent's front); or keep them and let M2-04/M2-05 treat them as
  inland basins (for example thinned continental crust instead of oceanic
  crust).

Out:
- Larger inland seas, which are a feature, not an artefact.

## Implementation notes

To be refined: the size limit (cells or steradians, so that it does not
depend on resolution), how the fraction is kept, and whether it changes the
golden values of `crust_type` (world-breaking before 1.0, allowed with a
`World-breaking:` line).

## Acceptance criteria

- [ ] No enclosed oceanic component below the chosen size for seeds 1–20
      at 10,000 and 100,000 cells (or the owner's alternative is met).
- [ ] Continental fraction still within ±1 percentage point of the target.
- [ ] `crust_type.png` viewed.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo run --release -p kuniumi-cli -- generate --seed 1 --cells 100000 --out out/M2-10
```

## Result

_(filled in when the ticket is implemented)_

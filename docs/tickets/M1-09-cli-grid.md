# M1-09 Grid stage, `kuniumi grid` and golden tests

- Status: ready
- Milestone: M1
- Depends on: M1-08
- Crates: kuniumi-world, kuniumi-sim, kuniumi-cli
- Design: 球面グリッドとデータ構造; 実装ロードマップ
- Decisions: D-07, D-10, D-14

## Goal

The pipeline builds the grid as its first stage, and a CLI command shows the
grid's quality as images and statistics. M1 is complete.

## Scope

In:
- `kuniumi_world::Settings { cells: u32 }`: generation settings, separate from
  the planet parameters; default 650,000; at least 1,000.
- `generate(params, &settings, seed, progress)` runs the Grid stage (random
  streams under `"grid"`), stores the grid in the `World`, then returns
  `NotImplemented { stage: Tectonics }`.
- CLI argument parsing moves to `clap` (derive, with subcommands); the
  hand-written parser is removed and the existing commands keep working.
- `kuniumi grid --cells N --seed S --out DIR` writes `cell_area.png`
  (equirectangular), `cell_area_globe.png`, `neighbors.png` (ring size,
  categorical) and `stats.txt` (cell count, minimum/mean/maximum area and their
  ratio, ring-size histogram, build time).
- A golden test of the Grid stage output for 10,000 cells, seed 1.
- `docs/roadmap.md`: M1 marked as met.

Out:
- Later stages.

## Acceptance criteria

- [ ] `kuniumi grid --cells 650000 --seed 1 --out out/M1-09` succeeds; the
      implementer views all three images and describes them (an even area map
      with no seams at the antimeridian or the poles).
- [ ] `stats.txt` agrees with the invariants of M1-06.
- [ ] `kuniumi generate --seed 1 --out out/gen` reports the Grid stage, then
      a clear "tectonics not implemented yet" error with exit code 2.
- [ ] `kuniumi --help` and the help of every subcommand list all options.
- [ ] `clap` is listed in the D-10 table.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo run --release -p kuniumi-cli -- grid --cells 650000 --seed 1 --out out/M1-09
```

## Result

_(filled in when the ticket is implemented)_

# MV-01 Read world directories back

- Status: ready
- Milestone: MV
- Depends on: M2-01
- Crates: kuniumi-io
- Design: 高解像度化と出力 / ゲームへの出力
- Decisions: D-19

## Goal

A world written by `kuniumi generate` can be loaded back into a `World`, so the
viewer and later tools can open it.

## Scope

In:
- `kuniumi_io::world_dir::read(dir) -> Result<World, ReadError>`: parse
  `world.json` with `serde_json` into the mirror structs of M2-01, rebuild or
  load the grid, and read every listed `.kfield`.
- Grid storage: either rebuild the grid from the stored settings and seed, or
  store it (positions and CSR) as `.kfield`-style files when writing. Measure
  the load time for 650,000 cells and choose; document the choice.
- Version checks: unknown `world.json` or `.kfield` versions give clear errors.

Out:
- Reading worlds produced by generator versions with different grids (error
  instead).

## Acceptance criteria

- [ ] `read(write(world)) == world` (fields bit-exact, identical grid) for a
      20,000-cell world.
- [ ] Malformed or truncated `world.json` and `.kfield` files give clear errors
      (tests).
- [ ] Load time for a 650,000-cell world recorded.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-io world_dir
```

## Result

_(filled in when the ticket is implemented)_

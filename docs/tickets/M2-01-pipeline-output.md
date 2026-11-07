# M2-01 Pipeline runner and world output

- Status: ready
- Milestone: M2
- Depends on: M1-09
- Crates: kuniumi-world, kuniumi-sim, kuniumi-io, kuniumi-cli
- Design: パイプライン全体像
- Decisions: D-07, D-10

## Goal

Stages plug into one runner, every stage output lands in the `World` as a named
SI field, and `kuniumi generate` writes whatever exists so far as `.kfield`
files, PNG maps and a `world.json` metadata file.

## Scope

In:
- A stage runner in `kuniumi-sim`: runs `StageId::PIPELINE` in order, reports
  progress, stops with `NotImplemented` at the first stage without an
  implementation, and (behind a `Settings` flag) returns the partial world
  instead of an error so the CLI can write it.
- `World` helpers: `add_field(name, unit, values)` (rejects duplicate names and
  wrong lengths) and `field_mut`; a way to mark a field as a class field so it
  is rendered with the categorical palette.
- `kuniumi_io::world_dir::write(dir, &World)`: one `<name>.kfield` and one
  `<name>.png` per field, and `world.json` with the generator version, seed,
  settings, planet parameters, sea level and the field list with units.
- `world.json` is written with `serde` (derive) and `serde_json`, through
  mirror structs that live in `kuniumi-io`, so `kuniumi-world` and
  `kuniumi-sim` do not depend on serde. Field order is fixed (struct order,
  never a hash map).

Out:
- Reading worlds back (MV-01).

## Acceptance criteria

- [ ] Runner test with fake stages covers order, progress callbacks and the stop
      behaviour.
- [ ] `world.json` parses back with `serde_json` into the same values, every
      `f64` bit-exactly; writing the same world twice gives identical bytes.
- [ ] `kuniumi generate --seed 1 --cells 20000 --out out/M2-01` writes the grid
      outputs and `world.json`, then exits with code 2 naming the next stage.
- [ ] `serde` and `serde_json` are listed in the D-10 table;
      `cargo tree -p kuniumi-sim` shows no serde.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo run -p kuniumi-cli -- generate --seed 1 --cells 20000 --out out/M2-01
cargo tree -p kuniumi-sim
```

## Result

_(filled in when the ticket is implemented)_

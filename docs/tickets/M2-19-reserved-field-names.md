# M2-19 Field names that collide with the files of a world directory

- Status: draft
- Milestone: M2
- Depends on: MV-01
- Crates: kuniumi-world (Field::is_valid_name) or kuniumi-io (world_dir)
- Design: パイプライン全体像
- Decisions: none

## Goal

No field name can make `world_dir::write` overwrite one of its own files.

## Scope

In:
- A field named `hypsometry` passes `Field::is_valid_name`; in a world with
  a `surface_height`, `write` writes its map as `hypsometry.png` and then
  overwrites it with the chart of the metrics (`HYPSOMETRY_FILE`), so the
  directory lists a map that is not there. Found while implementing MV-01;
  no stage uses the name. Reject such names (in `is_valid_name`, or in
  `write` and `read` of `world_dir`), or rename the chart's file.
- Check the other fixed names of a world directory (`world.json`,
  `grid.kgrid`, `metrics.txt`) the same way (none collides with
  `<name>.kfield` or `<name>.png` today).

Out:
- Other naming rules for fields.

## Implementation notes

To be refined: where the rule lives (`kuniumi-world` knows nothing of
files; `world_dir` does), and whether `read` rejects a listed `hypsometry`
too.

## Acceptance criteria

- [ ] A world with a field named `hypsometry` and a `surface_height` is
      rejected before anything is written, or its map and the chart are
      both kept (test).
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-io world_dir
```

## Result

_(filled in when the ticket is implemented)_

# M3-07 Erosion settings in `Settings`, `world.json` and settings files

- Status: draft
- Milestone: M3
- Depends on: M3-03
- Crates: kuniumi-world, kuniumi-io, kuniumi-config, kuniumi-cli
- Design: 侵食と地形発達 / 基本方程式
- Decisions: D-15

## Goal

The erosion settings that M3-03 defines (`kuniumi_world::ErosionSettings`:
erodibility K, discharge exponent m, time step and step count) become part
of a world's generation settings: checked with the others, recorded in
`world.json` and written to and read from settings files, so that the
Erosion stage of M3-05 can take them from `World::settings`.

## Scope

In:
- `Settings::erosion: ErosionSettings`, defaulting to
  `ErosionSettings::default()`; `Settings::validate` checks it after the
  existing fields (`ErosionSettings::validate`).
- `world.json`: the erosion settings under `settings` (SI, every `f64`
  written so that it reads back bit for bit), format version 3; how a
  version-2 file is treated (refused, as other versions are now, or read
  with the default erosion settings) is decided when the ticket is
  refined.
- Settings files (D-15): keys for the four settings, with the human-unit
  key for the time step (years) beside the SI one, the checks and messages
  of the existing keys, and `load(&write(p, s)) == (p, s)` bit for bit.
- The CLI's tests that list the lines of `world.json` and of settings files.

Out:
- Wiring the Erosion stage into the pipeline (M3-05).
- Physical bounds tighter than the domain of `ErosionSettings::validate`
  (M3-03), unless the refinement shows that user input needs them (ask the
  owner).

## Implementation notes

Split off M3-03 by the owner's decision (2026-11-15): with the solver, the
settings plumbing would have made one commit of about 1,800 lines in five
crates. To be refined: the TOML table and key names (the unit of K,
m^(1−3m) s^(m−1), depends on m, so its key cannot spell a fixed unit),
the format version history in `kuniumi_io::world_dir`, and the error names
of `ParamError::OutOfRange` for the nested settings.

## Acceptance criteria

- [ ] `Settings::validate` rejects invalid erosion settings by name; the
      defaults are valid.
- [ ] `world_dir::read(write(world)) == world` for non-default erosion
      settings, bit for bit; `world.json` documents format version 3.
- [ ] `config::load(&config::write(p, s)) == (p, s)` for non-default
      erosion settings, bit for bit; each new key rejected by key when out
      of range, and the time step's two keys together rejected.
- [ ] `mise run ci` passes.

## Verification

```sh
mise exec -- cargo test -p kuniumi-world -p kuniumi-io -p kuniumi-config -p kuniumi-cli
```

## Result

_(filled in when the ticket is implemented)_

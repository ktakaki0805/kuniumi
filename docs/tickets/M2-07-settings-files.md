# M2-07 Planet and generation settings from TOML files

- Status: ready
- Milestone: M2
- Depends on: M2-01
- Crates: new `kuniumi-config`, kuniumi-math, kuniumi-cli
- Design: プロジェクト方針 / 規約の要点
- Decisions: D-09, D-10, D-15

## Goal

Planets and generation settings are described in TOML files, so non-Earth
planets can be generated without changing code.

## Scope

In:
- New crate `kuniumi-config` using `toml` and `serde` (derive). Deserialise
  into raw structs whose fields are `Option`s, one per accepted key, with
  `#[serde(deny_unknown_fields)]`, then convert to `PlanetParams` and
  `Settings`. `toml` and `serde` stay out of `kuniumi-world` and `kuniumi-sim`.
- `load(&str) -> Result<(PlanetParams, Settings), ConfigError>` with tables
  `[planet]` and `[generation]`. Each quantity accepts its SI key or one
  human-unit key; giving both is an error:

  | Quantity | SI key | Human key |
  | --- | --- | --- |
  | radius | `radius_m` | `radius_km` |
  | surface gravity | `surface_gravity_m_s2` | — |
  | rotation period | `rotation_period_s` | `rotation_period_h` |
  | obliquity | `obliquity_rad` | `obliquity_deg` |
  | orbital period | `orbital_period_s` | `orbital_period_days` |
  | eccentricity | `orbital_eccentricity` | — |
  | stellar flux | `stellar_flux_w_m2` | — |
  | ocean volume | `ocean_volume_m3` | `ocean_volume_km3` |
  | surface pressure | `surface_pressure_pa` | `surface_pressure_bar` |
  | continental crust fraction | `continental_crust_fraction` | — |
  | rotation regime | `rotation` = `"free"` or `"tidally_locked"` | — |
  | tectonic regime | `tectonics` = `"plate_tectonics"` or `"stagnant_lid"` | — |

  Missing keys take Earth's value. Unknown keys, wrong types and invalid values
  (via `PlanetParams::validate`) are errors that name the key.
- Unit conversions through `kuniumi_math::units` (add the missing helpers).
- `write(&PlanetParams, &Settings) -> String` using SI keys.
- CLI: `kuniumi params` prints the Earth settings file; `kuniumi generate
  --planet FILE` reads one; `--seed` and `--cells` on the command line override
  the file.
- `CLAUDE.md`: `kuniumi-config` moves from planned to existing in the workspace
  list. `docs/roadmap.md`: if this is the last open M2 ticket, mark M2 as met.

Out:
- Other file formats; presets other than Earth.

## Acceptance criteria

- [ ] Round trip: `load(write(p, s))` returns `p` and `s` bit-exactly for Earth
      and for a random valid planet.
- [ ] Human-unit keys convert correctly (for example `obliquity_deg = 23.44`,
      `radius_km = 6371`), checked against hand-computed SI values.
- [ ] Errors for both keys given, an unknown key, a wrong type, an out-of-range
      value and malformed TOML, each tested, each naming the key or position.
- [ ] `toml` is listed in the D-10 table; `cargo tree -p kuniumi-sim` shows no
      `toml` or `serde`.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo run -p kuniumi-cli -- params > out/earth.toml
cargo run -p kuniumi-cli -- generate --planet out/earth.toml --seed 1 --cells 20000 --out out/M2-07
cargo tree -p kuniumi-sim
```

## Result

_(filled in when the ticket is implemented)_

# M4-16 Warnings of the stages in the world: an unconverged climate

- Status: draft
- Milestone: M4
- Depends on: M4-07
- Crates: kuniumi-world, kuniumi-sim (stages/climate), kuniumi-io, kuniumi-cli
- Design: 気候モデル / 計算の流れ
- Decisions: D-22

## Goal

A world records the warnings of its stages, its world directory keeps them,
and `kuniumi generate` prints them: first those of a climate that has not
reached its periodic state.

## Scope

In (split from M4-07 with the owner's approval when M4-07 was refined):
- `World::warnings` (human-readable messages of the stages, in the order
  raised, each naming its stage), written to `world.json` (format version
  5) and read back; `kuniumi generate` prints each on the standard error.
- The Climate stage's warnings (decided with the owner in M4-02): a planet
  whose energy balance has not converged within its 100 years
  (`energy_balance::Temperature::converged`, M4-02; and the sea surface's,
  M4-05) keeps the last year's temperatures and is reported with the last
  change; the warning says that the climate may be in transition (near an
  ice–albedo bifurcation the last year can lie mid-way into a snowball:
  Earth at 1,296 W/m² changed by 11.5 K in its 100th year). The moisture's
  `converged`, `unsettled_steps` and `failed_solves` (M4-06) and the
  currents' solves (`OceanCurrents::solves`, M4-04) likewise.
- Measured in M4-07 (1,000 cells, seed 3, the 128 corners of the v1 box):
  the energy balance of 4 corners and the sea surface of 2 not periodic
  after 100 years, the moisture of 11 after its 10 years; no unsettled
  sub-step, no failed solve.

Out:
- The viewer's display of the warnings (later, if wanted).

## Implementation notes

To be refined: the form of a warning (a string or a struct with the stage),
the wording, the format version's read path.

## Acceptance criteria

- [ ] A world with warnings writes and reads them back bit for bit;
      `world.json` format version 5.
- [ ] The corners of M4-07 that do not converge get the warnings, Earth none.
- [ ] `kuniumi generate` prints them.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-sim climate
cargo test -p kuniumi-io world_dir
cargo test -p kuniumi-cli
```

## Result

_(filled in when the ticket is implemented)_

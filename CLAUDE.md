# kuniumi — working guide for Claude Code

kuniumi is a physically based planet generator for games, written in Rust. It
builds a planet from causes rather than noise: plate tectonics → crust →
climate ⇄ erosion → hydrology → biomes. Realism comes first; generation runs
offline and may take tens of minutes.

## Read before working

- `docs/conventions.md` — units, coordinates, determinism, dependencies,
  commits. **Mandatory.**
- `docs/design-ja.md` — the design: every stage, its model and its parameters
  (Japanese). Tickets and stage modules point to its sections.
- `docs/decisions.md` — the decisions (D-01 …) and the open questions.
- `docs/tickets/README.md` — the work, as tickets of one commit each, in order.
- `docs/roadmap.md` — milestones and their exit criteria.

## Commands

The tools are pinned in `mise.toml` (D-23): the development Rust toolchain,
Rust 1.85.0 for the MSRV check, and cargo-deny.

```sh
mise install          # once per machine, and after mise.toml changes
mise run format       # cargo fmt --all
mise run check        # fmt check, clippy -D warnings, tests in debug and release
mise run ci           # check + MSRV build + cargo-deny: the gate for every change
cargo run -p kuniumi-cli -- help
```

Run checks through `mise run` (or `mise exec -- cargo …`) so they use the
pinned toolchains; with `mise activate` in the shell, plain `cargo` does too. A
change is done only when `mise run ci` passes. If `mise` is missing or
`mise install` fails, stop and tell the owner; do not fall back to another
toolchain.

Development is local (D-21): nothing is pushed until the owner decides, so CI
does not run yet. When it does, its jobs run the same mise tasks, plus the
release tests on Linux ARM64, macOS and Windows weekly (D-25).

## Workspace

```
kuniumi-math    deterministic math (libm), Rng, Vec3, SI units
kuniumi-sphere  spherical grid: points, Delaunay/Voronoi, CSR neighbours
kuniumi-world   PlanetParams, Settings, World, Field
kuniumi-config  TOML settings files: load and write PlanetParams, Settings
kuniumi-sim     the pipeline and its stages (no I/O)
kuniumi-io      .kfield files, world directories, PNG, fingerprints
kuniumi-cli     the `kuniumi` binary
```

Dependencies point downwards. `kuniumi-sim` never depends on `kuniumi-io` and
never touches the file system, clock or environment. Planned: `kuniumi-viewer`
(wgpu, native and web), `kuniumi-gpu`.

## Non-negotiable rules

1. **Determinism.** Same seed + same generator version ⇒ bit-identical output on
   every platform. Follow `docs/conventions.md#determinism`. Clippy blocks
   `f32`, std transcendental functions, `powi`, `HashMap`, `HashSet` and
   `available_parallelism`; never silence these lints in library code (only
   the viewer and GPU crates may allow `f32`).
2. **SI units** everywhere inside the code; unit suffixes on quantity names.
3. **Dependencies** follow D-10: well known, licence passes cargo-deny, pinned
   to an exact version when on the generation path, kept in the crate that
   needs it, listed in the D-10 table. Ask the owner before adding anything
   little known or borderline. Python only under `tools/` (D-11).
4. **Golden values** change only for an intentional world-breaking change,
   stated in the commit's `World-breaking:` line and in `CHANGELOG.md`.
5. **Scope.** No human geography in v1. v2 features (tidally locked, stagnant
   lid) keep their parameters and return `ParamError::Unsupported`.
6. **Git stays local.** One squashed commit per ticket on the local `main`.
   Never run `git push`, `git pull` or `git fetch`, never change remotes, never
   rewrite the history of `main`. The owner decides when to push.

## How to work

- Implement a ticket with the `implement-ticket` skill (`/implement-ticket
  M1-01`, or `/implement-ticket` for the next one). It refines a draft ticket
  and asks the owner to approve the plan, implements test-first on a local
  branch, runs `mise run ci`, loops with the independent `ticket-reviewer`
  subagent until it approves, then squashes the work into one commit on the
  local `main`.
- To run several tickets in a row, the owner starts `implement-all-tickets`
  (`/implement-all-tickets`, or with a milestone such as `M3`, a ticket ID such
  as `M2-05`, or a count such as `5`). It hands each ticket to a fresh
  `ticket-implementer` subagent that follows `implement-ticket`, relays
  questions to the owner, verifies every landed commit, and stops at the
  target, at a question nobody answers, or at the first ticket that cannot
  land.
- A ticket's status lives only in its `- Status:` line; the ticket's own commit
  sets it to `done`.
- Work discovered along the way becomes a new `draft` ticket, not a bigger
  commit.
- If the design document and the code disagree, or a design choice is unclear,
  ask the owner. If the design changes, update `docs/design-ja.md` in the same
  commit.

## How to check your own work

- Prefer invariant tests: CSR validity, conservation of mass and water, no
  cycles in drainage graphs, areas summing to 4π.
- Give every stage output a golden fingerprint test on a small grid with a
  fixed seed.
- Render anything spatial to a PNG under `./out/` (git-ignored) and look at it
  before calling the work done.
- Compare against Earth with the metrics of the design document's "検証方法"
  section.

Code, comments and technical documents are in English; reply to the owner in
Japanese.

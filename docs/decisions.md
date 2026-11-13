# Decisions

The current decisions of the project and the reasons for them, as of
2026-11-08. When a decision changes, edit it in place and append "(changed
YYYY-MM-DD: what changed)". A new decision takes the next free number. Open
questions are at the end.

## Project

**D-01 Name.** *kuniumi* (国生み, the Japanese creation myth in which the
islands are born). Crates are named `kuniumi-*`; the binary is `kuniumi`.

**D-02 Licence.** Dual-licensed `MIT OR Apache-2.0`, the Rust ecosystem
convention; Apache-2.0 adds an explicit patent grant. Copyright holder: Kousuke
Takaki.

**D-03 Repository and publication.** The repository is private until 1.0. The
owner creates the GitHub repository; its URL is an open question. Nothing is
published to crates.io before 1.0 (`publish = false`); crate names for
publication are decided then.

**D-04 Disclosure and community files.** The README does not describe how the
code is written. No CONTRIBUTING, CODE_OF_CONDUCT or SECURITY files while the
repository is private; revisit at 1.0.

**D-05 Documentation language.** Code, comments, README, `CLAUDE.md`,
conventions, decisions, tickets, skills and agents are in English. The design
document (`docs/design-ja.md`) is in Japanese. Claude replies to the owner in
Japanese.

## Architecture

**D-06 Workspace.** A cargo workspace; dependencies point downwards:

| Crate | Role |
| --- | --- |
| `kuniumi-math` | deterministic math, random numbers, vectors, SI units |
| `kuniumi-sphere` | the spherical cell grid |
| `kuniumi-world` | planet parameters, settings, the world data model |
| `kuniumi-config` | settings files (TOML, D-15) |
| `kuniumi-sim` | the generation pipeline (no I/O) |
| `kuniumi-io` | file formats, images, fingerprints |
| `kuniumi-cli` | the `kuniumi` command |

`kuniumi-sim` never depends on `kuniumi-io`. Planned crates: `kuniumi-viewer`
(D-19), `kuniumi-gpu` (D-08), `kuniumi-ffi` (engine bindings, D-07).
(changed 2026-11-13: `kuniumi-config`, built by M2-07, moved from the planned
crates into the table)

**D-07 Generation API and output.** `kuniumi_sim::generate` does no I/O and
never reads global state, the clock or the machine's core count, so it can be
embedded anywhere. v1 writes files only; engine integrations (C ABI for Unity
and Unreal, gdext for Godot, Bevy directly) can be added later without changing
the API.

**D-08 Reproducibility.** The CPU implementation in `f64` is the reference:
the same seed and generator version produce bit-identical output on Linux,
macOS and Windows, on x86_64 and ARM64, in debug and release builds. A GPU path
(`kuniumi-gpu`, `wgpu`, `f32`) may come later and is only statistically
equivalent to the reference. Reasons: wgpu's shader language has no portable
`f64` (the `SHADER_F64` feature is Vulkan-only, native-only and much slower),
GPU reductions and math functions differ between vendors, and some algorithms
(the implicit stream-power solve) need a different formulation on the GPU.

**D-09 Units and coordinates.** SI base units and radians everywhere inside the
code and in file formats; planet-centred right-handed frame. Details in
`docs/conventions.md`.

**D-10 Dependencies.** Few dependencies, not zero. A well-known, actively
maintained crate may be used when it clearly saves work or improves quality:

1. Its licence passes cargo-deny.
2. On the generation path (anything whose output reaches world data) it must
   be pure Rust arithmetic: no platform libm, no CPU-feature dispatch that
   changes floating-point results, no exposed hash-order iteration. Such crates
   are pinned to an exact version (`=x.y.z`); upgrading one is treated as
   potentially world-breaking and the golden tests decide.
3. It stays in the crate that needs it (for example, `serde` does not enter
   `kuniumi-sim` because I/O uses it).
4. It is listed in the table below.

Clear cases may be added by Claude Code without asking; little-known or
borderline crates need the owner's approval first. Kept in-house on purpose:
the random number generator (frozen algorithm), the spherical Delaunay
triangulation (the owner's choice), the FFT (D-17) and the `.kfield` and
`.kgrid` formats (MV-01).

| Crate | Used in | Purpose | Pinned | Added by |
| --- | --- | --- | --- | --- |
| `libm` | kuniumi-math | deterministic math functions | exact | M1-01 |
| `robust` | kuniumi-sphere | exact orient3d predicate | exact | M1-03 |
| `png` | kuniumi-io | compressed PNG output | no | M1-08 |
| `clap` | kuniumi-cli | command-line parsing | no | M1-09 |
| `serde`, `serde_json` | kuniumi-io | `world.json` | no | M2-01 |
| `toml`, `serde` | kuniumi-config | settings files | no | M2-07 |
| `wgpu`, `winit`, `wasm-bindgen` family | kuniumi-viewer | viewer | no | MV-02, MV-03 |

Also acceptable when a ticket needs them: `rayon` (parallel loops, following
the determinism rules for parallel code), `thiserror` (error types) and
`proptest` (property tests, development only).

**D-11 Developer tools in Python.** Python scripts are allowed under `tools/`
for development only: converting reference data (NetCDF, GeoTIFF) into simple
formats the Rust code reads, and computing independent reference values for
tests. They are never part of the build. The Python version is pinned in
`mise.toml` (D-23) and packages in `tools/requirements.txt`, both added with
the first script.

## Scope

**D-12 Planet range.** v1 generates Earth-like planets with variable radius,
gravity, rotation, obliquity, orbit and ocean volume. Tidally locked rotation
and stagnant-lid tectonics are planned for v2: their parameters exist from v1
and are rejected with `ParamError::Unsupported`. The climate model must not
hard-code three circulation cells; their number scales with the rotation rate.
`PlanetParams::validate` bounds the planet (chosen with the owner in M2-08;
the reasons are with the constants of `PlanetParams`): radius 1,500 to
13,000 km, surface gravity 1 to 40 m/s², rotation period 8 to 240 h, orbital
period 90 to 3,650 days, stellar flux 400 to 2,500 W/m², ocean volume 0 to
10¹¹ km³, surface pressure 0.1 to 10 bar, and a bulk density `3g / (4πGR)`
of 3,000 to 16,000 kg/m³ (`ParamError::Inconsistent` otherwise). Within these
bounds every planet rotates slower than break-up, which a compile-time
assertion over the bounds keeps true instead of a run-time check.

**D-13 Out of scope.** Human geography (civilisations, borders, cities) is not
part of v1.

## Technical choices

**D-14 Grid.** Jittered Fibonacci points (about 650,000 cells by default, about
30 km apart on an Earth-sized planet), triangulated as their 3D convex hull by
code written in-house, with exact predicates from `robust`; the Voronoi dual
gives the cells.

**D-15 Settings files.** Human-written planet and generation settings are TOML,
parsed with `toml` and `serde` in `kuniumi-config`. Keys spell their unit; each
quantity accepts its SI key or one human-unit key.

**D-16 Reference data.** Calibration starts from ETOPO 2022 topography and the
Köppen–Geiger maps of Beck et al. (2018); observed temperature and
precipitation datasets are chosen in M4-08. Each dataset's licence and citation
is recorded here before use. Data is downloaded by scripts and never committed.

**D-17 FFT.** Written in-house (radix-2), so results are identical on every CPU;
external FFT crates choose different code paths per instruction set.

**D-18 Biomes.** Vegetation biomes use the Whittaker scheme (annual mean
temperature and precipitation); special biomes are layered on top.

**D-19 Viewer.** `kuniumi-viewer`, built on `wgpu`, runs natively and in the
browser (WebAssembly with WebGPU, WebGL2 fallback). It may use `f32` for GPU
data (crate-level lint allowance); nothing it computes feeds back into
generation.

## Process

**D-20 Tickets and review loop.** Work is split into tickets of one commit each
in `docs/tickets/`, each with measurable acceptance criteria. Tickets of the
near milestones are written in full; later ones start as `draft` and are
refined just before implementation, with the owner approving each refined
plan. The `implement-ticket` skill implements one ticket and loops with the
independent `ticket-reviewer` subagent (a fresh context each round, at most
three rounds). The `implement-all-tickets` skill, started only by the owner,
runs tickets in order, each in a fresh `ticket-implementer` subagent that
follows `implement-ticket`; it relays questions to the owner and stops at the
first ticket that cannot land. A ticket's status lives only in its own file.

**D-21 Integration.** Development is local for now. A ticket is built and
reviewed on a local `ticket/<id>-*` branch and squashed into one commit on the
local `main`; there are no pull requests. Skills and agents never push, pull or
fetch and never change remotes; the owner decides when to push. Nobody
force-pushes or rewrites the history of `main`. Until the first push CI does
not run: `mise run ci` (D-23) covers its Linux checks, and agreement of golden
values across platforms is confirmed by CI after the first push.

**D-22 Changelog.** `CHANGELOG.md` in Keep a Changelog format, with a
*World-breaking* section listing every change after which an existing seed
produces a different world. Each ticket's commit updates it.

**D-23 Development environment.** [mise](https://mise.jdx.dev) pins the tools
in `mise.toml`: the development Rust toolchain (with rustfmt and clippy), Rust
1.85.0 for the MSRV check, and cargo-deny. `mise install` sets up a machine.
`mise run check` runs fmt, clippy, the documentation build (`mise run doc`:
`cargo doc --workspace --no-deps` with warnings as errors, added by M2-09) and
the tests in debug and release; `mise run ci` adds the MSRV build and
cargo-deny and is the gate for every change. The development toolchain is
raised on purpose, in its own commit, because new clippy and rustdoc lints and
rustfmt changes arrive with it; a toolchain change must not change any golden
value (D-08).

**D-24 MSRV.** `rust-version = "1.85"` (the edition 2024 minimum), checked by
`mise run msrv`. Raising it is recorded here.

**D-25 CI.** In a private repository Actions minutes are billed against the
account's included amount (GitHub Free: USD 12 per month; macOS minutes cost
about ten times Linux minutes). Every push runs the Linux checks (fmt, clippy,
rustdoc, tests in debug and release, MSRV, cargo-deny) on Linux x86_64; the
release tests on Linux ARM64, macOS and Windows run weekly and on demand. Every
job installs its tools with `jdx/mise-action` and runs the mise tasks of D-23,
so CI and local checks use the same toolchains and commands.

## Open questions

- **Game engine and map use.** Which engine the game uses and how it uses the
  map (globe, world map, hex board, regional maps). Needed before M8-03 and
  before any engine bindings.
- **Repository URL.** To be filled into `Cargo.toml` and the README once the
  owner has created the repository.

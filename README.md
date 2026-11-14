# kuniumi

*kuniumi* (国生み, "the birth of the land") is a physically based planet
generator for games. Instead of painting terrain with noise, it simulates the
causes: plate tectonics build crust, isostasy turns crust into relief, a climate
model brings wind and rain, erosion carves rivers and valleys, and hydrology and
biomes follow from the result.

> **Status: pre-alpha.** The workspace, data model and determinism tooling are
> in place; the generation stages are implemented ticket by ticket. See
> [docs/roadmap.md](docs/roadmap.md).

## Principles

- **Causes, not noise.** Plates → crust → climate ⇄ erosion → hydrology →
  biomes, each stage a simplified model from the Earth sciences.
- **Whole planets.** A spherical Voronoi grid of about 650,000 cells (about
  30 km apart on an Earth-sized planet), refined regionally on demand.
- **Reproducible.** The same seed and version produce bit-identical worlds on
  every supported platform.
- **Beyond Earth.** Radius, gravity, rotation, obliquity, orbit and ocean volume
  are parameters; tidally locked and stagnant-lid planets are planned.
- **Few dependencies.** Only well-known crates, pinned to exact versions
  wherever they affect the numbers.

## Building

Requires Rust 1.85 or later.

```sh
cargo build --release
cargo run --release -p kuniumi-cli -- help
cargo run --release -p kuniumi-cli -- params
cargo run --release -p kuniumi-cli -- grid --cells 650000 --seed 1 --out out/grid
cargo run --release -p kuniumi-cli -- generate --seed 1 --out out/world
cargo run --release -p kuniumi-cli -- params > planet.toml
cargo run --release -p kuniumi-cli -- generate --planet planet.toml --seed 1 --out out/planet
```

`params` prints the settings file of Earth with the default settings (TOML):
the planet's parameters under `[planet]` and the generation settings under
`[generation]`. Keys spell their SI unit (`radius_m`); some quantities may be
given in a human unit instead (`radius_km`, `obliquity_deg`,
`rotation_period_h`, `orbital_period_days`, `ocean_volume_km3`,
`surface_pressure_bar`), and missing keys take Earth's values. `generate
--planet FILE` reads such a file; `--cells` on the command line replaces the
file's.

`generate` runs the stages that exist so far and writes the world into the
output directory: the cell grid (`grid.kgrid`); a `.kfield` file and a PNG
map per field; once the Crust stage has run, `metrics.txt` (land fraction
and hypsometry against Earth's values) and `hypsometry.png` (the elevation
histogram); and `world.json` with the metadata. While later stages are
missing it then exits with code 2, naming the first of them. The directory
can be read back into the same world (`kuniumi_io::world_dir::read`).

A world directory can be explored as a 3D globe (Rust 1.87 or later, and a
GPU with Vulkan, Metal or Direct3D 12):

```sh
cargo run --release -p kuniumi-viewer -- out/world
cargo run --release -p kuniumi-viewer -- out/world --layer plate_id --screenshot out/plates.png
```

Drag to turn the globe and use the wheel to zoom; Left and Right choose the
field shown, in the colours of its PNG map; Up and Down change the relief
exaggeration (up to ×100, less on a world whose deepest point it would draw
closer to the centre than half the radius); S shows or hides the sea
surface; R resets the view; Esc or Q quits. The cell under the cursor is
highlighted and its values are shown in a panel at the top left with the
legend of the field shown (P hides and shows it); a click prints them to the
terminal. `--help` lists the options
(start view, screenshots without a window, with `--cursor` the cell under a
pixel and the panel, a frame-time benchmark).

The same viewer runs in a web browser (WebGPU, or WebGL2 where the browser
has no WebGPU) and opens a world directory chosen on the machine:
`mise run web-serve` builds the page and serves it on
`http://127.0.0.1:8000/`; see
[crates/kuniumi-viewer/README.md](crates/kuniumi-viewer/README.md).

For development, the toolchains and tools are pinned in `mise.toml` for
[mise](https://mise.jdx.dev):

```sh
mise install    # set up the pinned tools
mise run ci     # formatting, lints, docs, tests, MSRV build and licence checks
```

## Layout

| Crate | Purpose |
| --- | --- |
| `kuniumi-math` | Deterministic math, random numbers, vectors, SI units |
| `kuniumi-sphere` | The spherical cell grid |
| `kuniumi-world` | Planet parameters, settings and the world data model |
| `kuniumi-config` | Settings files (TOML) |
| `kuniumi-sim` | The generation pipeline (no I/O) |
| `kuniumi-io` | File formats and images |
| `kuniumi-cli` | The `kuniumi` command |
| `kuniumi-viewer` | The `kuniumi-viewer` command and its browser page: a world as a 3D globe |

Documentation: [design](docs/design-ja.md) (Japanese) ·
[conventions](docs/conventions.md) · [decisions](docs/decisions.md) ·
[roadmap](docs/roadmap.md) · [tickets](docs/tickets/README.md) ·
[changelog](CHANGELOG.md)

## License

Licensed under either of

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE))
- MIT license ([LICENSE-MIT](LICENSE-MIT))

at your option.

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in the work by you, as defined in the Apache-2.0 license, shall be
dual licensed as above, without any additional terms or conditions.

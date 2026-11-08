# Roadmap

Milestones follow the design document (`docs/design-ja.md`, "実装ロードマップ"):
first connect every stage in a simple form (M5 produces the first complete
planet), then make each stage realistic. The work itself is in
[docs/tickets](tickets/README.md); this file records each milestone's exit
criterion, and the ticket that meets it sets "Met" to yes.

| Milestone | Exit criterion | Tickets | Met |
| --- | --- | --- | --- |
| M0 Scaffold | Workspace, determinism lints, Rng, grid data structure, parameters, PNG and `.kfield`, development environment | — | yes |
| M1 Spherical grid and map rendering | Any scalar field on a 650k-cell grid renders to PNG; the grid validates; golden fingerprints are recorded | M1-01 … M1-09 | yes |
| M2 One-shot plates and initial relief | Land fraction and the bimodal hypsometry appear; planets load from TOML files | M2-01 … M2-07 | no |
| MV Viewer | A 650k-cell world can be explored as a 3D globe, natively and in the browser | MV-01 … MV-03 | no |
| M3 Erosion and drainage | River networks follow Hack's law | M3-01 … M3-05 | no |
| M4 Climate | With Earth's real topography, zonal-mean temperature and precipitation broadly match observations | M4-01 … M4-08 | no |
| M5 Hydrology and biomes | A complete planet is generated and exported with a validation report; with Earth's real topography, the Köppen map broadly matches observations | M5-01 … M5-06 | no |
| M6 Plate tectonics through time | Ocean-floor age gradients, young and old orogens, passive margins | M6-01 … M6-07 | no |
| M7 Coupling and extra processes | The validation metrics fall within Earth's range | M7-01 … M7-05 | no |
| M8 Refinement and output | Regional detail at 100 m–1 km consistent with the global world; game-facing exports | M8-01 … M8-04 | no |

Golden values are checked on the local machine by `mise run ci`. Their
agreement across Linux, macOS and Windows on x86_64 and ARM64 is confirmed by
CI once the owner pushes the repository (D-21, D-25).

## Later

- GPU path (`kuniumi-gpu`, `wgpu`, `f32`, statistical-equivalence tests; D-08)
- v2: tidally locked planets and stagnant-lid tectonics (D-12)
- Engine integrations (`kuniumi-ffi`, Godot, Bevy; D-07)

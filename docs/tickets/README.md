# Tickets

The work, as tickets of one commit each (D-20). A ticket is implemented with
the `implement-ticket` skill: it refines the ticket if needed, implements it on
a local branch, runs `mise run ci`, loops with the independent
`ticket-reviewer` subagent until it approves, and squashes the work into one
commit on the local `main` (D-21). `implement-all-tickets` runs tickets in the
order below, one after another, each in a fresh `ticket-implementer` subagent.

## Format

Every ticket starts with the same header lines:

| Line | Content |
| --- | --- |
| `- Status:` | `draft`, `ready`, `blocked` or `done` (see below) |
| `- Milestone:` | `M1` … `M8`, or `MV` for the viewer |
| `- Depends on:` | ticket IDs that must be `done` first, or `none` |
| `- Crates:` | crates the ticket touches or creates |
| `- Design:` | sections of `docs/design-ja.md` (`section / subsection`, `;`-separated) |
| `- Decisions:` | entries of `docs/decisions.md` (`D-NN`), or `none` |

Then: Goal, Scope, Implementation notes (required for `draft`), Acceptance
criteria (the last one is always "`mise run ci` passes."), Verification and
Result. Copy [TEMPLATE.md](TEMPLATE.md) for new tickets.

## Status

The status of a ticket lives in one place only: its `- Status:` line.

| Status | Meaning |
| --- | --- |
| `draft` | Goal and acceptance criteria are set; the plan is not. It is refined just before implementation, and the owner approves the plan before coding. |
| `ready` | Detailed enough to implement directly. |
| `blocked` | Waiting on a decision; the reason is written in the ticket. |
| `done` | On `main`; set in the ticket's own commit. |

Work in progress is visible as a local `ticket/<id>-*` branch.

```sh
grep -H "^- Status:" docs/tickets/M*.md
```

## Order

Tickets are listed in the recommended order. The next ticket is the first one
that is `ready` or `draft` and whose dependencies are all `done`. Tickets of
M1–M3 and MV-01 are `ready`; the rest are `draft`, because they are best
planned once the code they build on exists.

### M1 — Spherical grid and map rendering

| ID | Title | Depends on |
| --- | --- | --- |
| [M1-01](M1-01-det-math.md) | Deterministic math functions via libm | — |
| [M1-02](M1-02-sphere-coordinates.md) | Spherical coordinates, rotations and distances | M1-01 |
| [M1-03](M1-03-exact-orient3d.md) | Exact orient3d predicate | — |
| [M1-04](M1-04-fibonacci-points.md) | Jittered Fibonacci point set | M1-02 |
| [M1-05](M1-05-spherical-delaunay.md) | Spherical Delaunay triangulation | M1-03, M1-04 |
| [M1-06](M1-06-voronoi-grid.md) | Voronoi dual and `SphereGrid::fibonacci` | M1-05 |
| [M1-07](M1-07-nearest-cell.md) | Nearest-cell lookup | M1-06 |
| [M1-08](M1-08-render-maps.md) | Colour maps and map rendering | M1-07 |
| [M1-09](M1-09-cli-grid.md) | Grid stage, `kuniumi grid` and golden tests | M1-08 |

### M2 — One-shot plates and initial relief

| ID | Title | Depends on |
| --- | --- | --- |
| [M2-01](M2-01-pipeline-output.md) | Pipeline runner and world output | M1-09 |
| [M2-02](M2-02-plate-partition.md) | Plate partition and Euler poles | M2-01 |
| [M2-03](M2-03-continents-boundaries.md) | Continents and boundary classification | M2-02 |
| [M2-04](M2-04-one-shot-crust.md) | One-shot crust thickness, age and uplift | M2-03 |
| [M2-05](M2-05-isostasy-sea-level.md) | Elevation from crust and sea level | M2-04 |
| [M2-06](M2-06-hypsometry-metrics.md) | Land fraction and hypsometry metrics | M2-05 |
| [M2-07](M2-07-settings-files.md) | Planet and generation settings from TOML files | M2-01 |

### MV — Viewer

Placed after M2 so there is relief to look at; it can also be worked on in
parallel with later milestones.

| ID | Title | Depends on |
| --- | --- | --- |
| [MV-01](MV-01-world-reader.md) | Read world directories back | M2-01 |
| [MV-02](MV-02-native-viewer.md) | Native globe viewer | MV-01, M2-06 |
| [MV-03](MV-03-web-viewer.md) | Viewer in the browser | MV-02 |

### M3 — Erosion and drainage

| ID | Title | Depends on |
| --- | --- | --- |
| [M3-01](M3-01-drainage-graph.md) | Drainage graph: receivers, stack, accumulation | M2-05 |
| [M3-02](M3-02-depression-routing.md) | Depression routing (Priority-Flood) | M3-01 |
| [M3-03](M3-03-stream-power.md) | Implicit stream-power erosion with uplift | M3-02 |
| [M3-04](M3-04-hillslope-diffusion.md) | Hillslope diffusion | M3-03 |
| [M3-05](M3-05-erosion-stage.md) | Erosion stage, river maps and Hack's law | M3-04, M2-06 |

### M4 — Climate

| ID | Title | Depends on |
| --- | --- | --- |
| [M4-01](M4-01-insolation.md) | Insolation by latitude and month | M2-01 |
| [M4-02](M4-02-energy-balance.md) | Energy-balance temperature model | M4-01 |
| [M4-03](M4-03-pressure-winds.md) | Pressure belts and surface winds | M4-02 |
| [M4-04](M4-04-ocean-gyres.md) | Ocean basins and wind-driven gyres | M4-03 |
| [M4-05](M4-05-sst-upwelling.md) | Sea-surface temperature transport and upwelling | M4-04 |
| [M4-06](M4-06-moisture-precipitation.md) | Moisture, precipitation and rain shadows | M4-05 |
| [M4-07](M4-07-climate-stage.md) | Climate stage and monthly outputs | M4-06 |
| [M4-08](M4-08-earth-calibration.md) | Earth calibration harness | M4-07 |

### M5 — Hydrology and biomes (first complete planet)

| ID | Title | Depends on |
| --- | --- | --- |
| [M5-01](M5-01-runoff.md) | Evapotranspiration and runoff | M4-07, M3-05 |
| [M5-02](M5-02-lakes.md) | Lakes by water balance | M5-01 |
| [M5-03](M5-03-rivers.md) | Rivers, seasonal rivers and vector export | M5-02 |
| [M5-04](M5-04-koppen.md) | Köppen–Geiger classification | M4-07, M4-08 |
| [M5-05](M5-05-biomes.md) | Biomes, NPP and special biomes | M5-04, M5-02 |
| [M5-06](M5-06-first-planet.md) | Full world export and validation report | M5-03, M5-05 |

### M6 — Plate tectonics through time

| ID | Title | Depends on |
| --- | --- | --- |
| [M6-01](M6-01-crust-samples.md) | Crust samples, plate rotation and resampling | M5-06 |
| [M6-02](M6-02-boundary-detection.md) | Boundary detection by overlap and gaps | M6-01 |
| [M6-03](M6-03-subduction.md) | Subduction, trenches and arcs | M6-02 |
| [M6-04](M6-04-ridge-spreading.md) | Ridge spreading and ocean-floor age | M6-02 |
| [M6-05](M6-05-collision.md) | Continental collision, suturing and plate merging | M6-03 |
| [M6-06](M6-06-rifting-forces-hotspots.md) | Rifting, plate speeds and hotspots | M6-04, M6-05 |
| [M6-07](M6-07-tectonics-switch.md) | Orogeny decay, uplift recording and pipeline switch | M6-06 |

### M7 — Coupling and extra processes

| ID | Title | Depends on |
| --- | --- | --- |
| [M7-01](M7-01-climate-erosion-loop.md) | Climate–erosion alternation | M6-07 |
| [M7-02](M7-02-deposition.md) | Sediment deposition and marine sedimentation | M7-01 |
| [M7-03](M7-03-isostatic-rebound.md) | Isostatic rebound | M7-02 |
| [M7-04](M7-04-lithology-glacial.md) | Lithology-dependent erodibility and glacial erosion | M7-03 |
| [M7-05](M7-05-calibration.md) | Calibration against the validation metrics | M7-04 |

### M8 — Refinement and output

| ID | Title | Depends on |
| --- | --- | --- |
| [M8-01](M8-01-regional-patch.md) | Regional patch extraction and boundary conditions | M7-05 |
| [M8-02](M8-02-regional-simulation.md) | Regional re-simulation and climate downscaling | M8-01 |
| [M8-03](M8-03-projections-export.md) | Projections and game-facing exports | M8-02 |
| [M8-04](M8-04-real-dem-detail.md) | Detail from real elevation data (optional) | M8-02 |

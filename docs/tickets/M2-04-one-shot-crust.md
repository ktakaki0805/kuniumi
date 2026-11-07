# M2-04 One-shot crust thickness, age and uplift

- Status: ready
- Milestone: M2
- Depends on: M2-03
- Crates: kuniumi-sim (stages/tectonics)
- Design: プレートテクトニクス / 境界ごとの処理; 地殻から標高へ
- Decisions: none

## Goal

Without time evolution, give each cell a plausible crust thickness, ocean-floor
age and uplift rate from its distance to the boundaries, so later stages have
realistic input until M6 replaces this.

## Scope

In:
- Distances (m) from each boundary type by multi-source Dijkstra on the cell
  graph, with centre distances × radius as edge weights (deterministic heap,
  ties by index).
- Ocean-floor age = distance to the nearest divergent boundary ÷ that
  boundary's half spreading rate, capped at 180 Myr (stored in seconds).
- Continental thickness: 35 km base; continent–continent convergence thickens
  the crust up to about 70 km over a few hundred kilometres (plateau);
  ocean–continent convergence builds an arc (about 45–50 km) inland of the
  trench. Profiles are documented with their parameters in SI.
- A trench flag on the oceanic side of subduction boundaries.
- `uplift_rate` (m/s): non-zero near convergent boundaries, scaled by the
  convergence rate; it feeds the erosion stage until M6.
- Outputs: `crust_thickness` (m), `crust_age` (s), `uplift_rate` (m/s),
  `trench` (class).

Out:
- Conversion to elevation (M2-05).

## Acceptance criteria

- [ ] On a synthetic two-plate divergent setup, age increases monotonically with
      distance from the ridge and equals distance ÷ half rate within 1e-12
      relative.
- [ ] Thickness ranges: oceanic 6–8 km, continental 25–75 km for the default
      seeds.
- [ ] Uplift rates between 0 and 10 mm/yr (expressed in m/s).
- [ ] Images of the four fields viewed and described.
- [ ] Golden fingerprints for 10,000 cells, seed 1.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo run -p kuniumi-cli -- generate --seed 1 --cells 100000 --out out/M2-04
```

## Result

_(filled in when the ticket is implemented)_

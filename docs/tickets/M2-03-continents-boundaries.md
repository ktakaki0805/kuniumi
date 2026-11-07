# M2-03 Continents and boundary classification

- Status: ready
- Milestone: M2
- Depends on: M2-02
- Crates: kuniumi-sim (stages/tectonics)
- Design: プレートテクトニクス / 境界ごとの処理; プレートテクトニクス / 初期条件の目安
- Decisions: none

## Goal

Continental crust covers the requested fraction of the surface, and every plate
boundary is classified as convergent, divergent or transform, with a rate.

## Scope

In:
- Continental cells: grow a few continental blobs inside plates (stream
  `"tectonics/continents"`) until the area-weighted continental fraction
  reaches `PlanetParams::continental_crust_fraction`. Some plates end up mostly
  oceanic.
- Boundary edges: for each neighbour pair on different plates, the relative
  velocity `v_j − v_i` projected on the edge normal (the great-circle direction
  from i to j) gives the convergence rate (m/s, positive when converging); the
  tangential part gives the shear. Convergent when the normal part dominates
  and is positive, divergent when it dominates and is negative, transform
  otherwise; the threshold angle is documented.
- Outputs: `crust_type` (class: 0 oceanic, 1 continental), `boundary_type`
  (class per cell: none, convergent, divergent, transform; the strongest of
  its edges) and `convergence_rate` (m/s).

Out:
- Thickness, age and elevation (M2-04, M2-05).

## Acceptance criteria

- [ ] Continental area fraction within ±1 percentage point of targets 0.1, 0.4
      and 0.7.
- [ ] Synthetic two-plate tests with hand-chosen poles give pure convergent,
      pure divergent and pure transform boundaries.
- [ ] `boundary_type.png` viewed: boundaries are continuous lines along plate
      edges.
- [ ] Golden fingerprints for 10,000 cells, seed 1.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo run -p kuniumi-cli -- generate --seed 1 --cells 100000 --out out/M2-03
```

## Result

_(filled in when the ticket is implemented)_

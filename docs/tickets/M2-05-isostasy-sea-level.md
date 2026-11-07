# M2-05 Elevation from crust and sea level

- Status: ready
- Milestone: M2
- Depends on: M2-04
- Crates: kuniumi-sim (stages/crust)
- Design: 地殻から標高へ
- Decisions: none

## Goal

The Crust stage turns crust into relief and solves the sea level from the ocean
volume, giving the first real elevation map.

## Scope

In:
- Ocean depth from age after Parsons & Sclater (1977), with age converted to
  Myr only inside the formula: `d = 2500 + 350·√t` for young and
  `d = 6400 − 3200·exp(−t/62.8)` for old ocean floor (m). The two published
  curves overlap and switching at 70 Myr would leave a 78 m step, so switch
  where they intersect (t ≈ 26.3 Myr, d ≈ 4294 m), computed once by bisection.
- Continental elevation by Airy isostasy:
  `h = (T − T_ref)·(ρ_m − ρ_c)/ρ_m + h_ref` with ρ_c = 2800 kg/m³,
  ρ_m = 3300 kg/m³, T_ref = 35 km, h_ref = 300 m. Document that gravity
  cancels in Airy isostasy.
- Trenches deepened to about 7–10 km.
- Sea level `z` by bisection so that `Σ max(0, z − h_i) · a_i · R² =
  ocean_volume_m3`; ocean volume 0 (no ocean) and more water than any basin
  holds (water world) behave as documented.
- Outputs: `surface_height` (m) and `World::sea_level_m`.

Out:
- Metrics and histograms (M2-06).

## Acceptance criteria

- [ ] Formula tests: depth(0) = 2500 m; depth is continuous across the switch
      (difference below 1e-6 m) and increasing in t; thickness 35 km → 300 m,
      70 km → 5603 m.
- [ ] Bisection: `|volume(z) − V| / V < 1e-9`; monotonicity test; V = 0 and a
      huge V behave as documented.
- [ ] Earth defaults, seeds 1–10: land fraction between 20% and 40%.
- [ ] `surface_height.png` viewed and described (continents, ridges and
      trenches visible).
- [ ] Golden fingerprints for 10,000 cells, seed 1.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo run --release -p kuniumi-cli -- generate --seed 1 --out out/M2-05
```

## Result

_(filled in when the ticket is implemented)_

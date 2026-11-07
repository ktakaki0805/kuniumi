# M5-05 Biomes, NPP and special biomes

- Status: draft
- Milestone: M5
- Depends on: M5-04, M5-02
- Crates: kuniumi-sim (stages/biomes)
- Design: 気候区分とバイオーム; 気候区分とバイオーム / 特殊バイオームの規則
- Decisions: D-18

## Goal

The vegetation map a player sees: Whittaker biomes from annual temperature and
precipitation, vegetation density from net primary production, and special
biomes from terrain, water and sea temperature on top.

## Scope

In:
- Whittaker biomes with a continuous suitability per biome; classify by the
  maximum and blend near boundaries (the blend weight is an output).
- Miami model NPP (g/m²/yr):
  `min(3000 / (1 + exp(1.315 − 0.119 T)), 3000 (1 − exp(−0.000664 P)))` with T
  in °C and P in mm/yr at the boundary of the formula only.
- Special biomes overriding the base: alpine above the tree line, mangroves,
  salt flats on dry lake beds, dune fields, coral reefs (depth < 50 m, sea
  surface above about 20 °C all year), monthly sea ice, wetlands on flat, wet
  lowlands.
- Outputs: `biome` (class field), `biome_blend`, `npp`.

## Implementation notes

To be refined: the suitability functions, the tree-line rule, thresholds for
mangroves, dunes and wetlands.

## Acceptance criteria

- [ ] NPP unit tests at reference points of the Miami model.
- [ ] Every special-biome rule has a synthetic unit test.
- [ ] Biome map for Earth input and a generated world viewed and described
      (deserts, rainforests, taiga and tundra in plausible places).
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-sim biomes
```

## Result

_(filled in when the ticket is implemented)_

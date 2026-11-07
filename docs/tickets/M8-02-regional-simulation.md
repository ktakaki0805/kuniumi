# M8-02 Regional re-simulation and climate downscaling

- Status: draft
- Milestone: M8
- Depends on: M8-01
- Crates: kuniumi-sim (regional), kuniumi-math (new `fft` module), kuniumi-cli
- Design: 高解像度化と出力 / 3つの方法; 侵食と地形発達
- Decisions: D-17

## Goal

Re-run erosion on the fine regional grid driven by the global uplift and the
boundary conditions (after Schott et al. 2023), with climate downscaled by
lapse rate and a linear orographic precipitation model, so regional detail
stays consistent with the global world.

## Scope

In:
- The erosion solver (stream power, diffusion, deposition, rebound) generalised
  to the regional grid (shared code with the global solver where possible).
- Temperature by lapse rate on the fine elevation.
- Orographic precipitation with the linear theory of Smith & Barstad (2004),
  using an in-house radix-2 FFT in `kuniumi-math` (D-17): deterministic
  butterflies, `det` twiddle factors, power-of-two sizes with padding.
- A `kuniumi region` command writing the regional world directory.

## Implementation notes

To be refined: FFT API, padding and windowing, iteration count and time step of
the regional erosion.

## Acceptance criteria

- [ ] FFT: forward then inverse reproduces the input within 1e-12; results match
      a naive DFT on random inputs; a golden bit-pattern test pins the output.
- [ ] Smith & Barstad: a Gaussian hill gives windward enhancement and a lee
      rain shadow matching the published shape.
- [ ] Regional rivers connect to the global inflows and outlets; ridges and
      major valleys agree with the global world at its resolution (compared and
      described).
- [ ] `mise run ci` passes.

## Verification

```sh
cargo run --release -p kuniumi-cli -- generate --seed 1 --out out/M8-02-world
cargo run --release -p kuniumi-cli -- region --world out/M8-02-world --center 35,135 --size-km 1000 --spacing-m 500 --out out/M8-02
```

## Result

_(filled in when the ticket is implemented)_

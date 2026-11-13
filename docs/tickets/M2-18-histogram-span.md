# M2-18 Elevation spans beyond the histogram in world directories

- Status: draft
- Milestone: M2
- Depends on: M2-15
- Crates: kuniumi-sim (metrics), kuniumi-io (world_dir)
- Design: 検証方法
- Decisions: none

## Goal

`world_dir::write` rejects with `InvalidInput`, before writing anything,
every world whose elevation histogram it cannot build, and never needs
gigabytes for it.

## Scope

In:
- `world_dir::write` panics (`Hypsometry::new`, "more than 2^31") for a
  world whose `surface_height` spans more than 2³¹ bins of 100 m (about
  2.1e11 m), documented under Panics since M2-15. Return `InvalidInput`
  instead, without duplicating the binning arithmetic of
  `kuniumi_sim::metrics` (for example a public check there, or a fallible
  `ElevationMetrics` constructor).
- Below that limit the histogram still allocates a share per bin (up to
  2³¹ × 8 B = 16 GB) and `metrics.txt` writes a line per bin: decide a
  practical bound (for example 10⁶ bins, 100,000 km of relief) for world
  directories, or for `Hypsometry` itself.

Out:
- Bounds on the planet parameters (M2-08), which keep the pipeline's
  relief to tens of kilometres.

## Implementation notes

Found while implementing M2-15: `ElevationMetrics::new(&[1.0, 1.0], &[0.0,
1e12], 0.0)` panics (test `metrics_need_at_most_2_pow_31_bins`), and a
world with such heights makes `write` panic after `check_fields` and
`world_json` accepted it. The bin index of a cell is `floor((h − z) / w)`,
or −1 for an ocean cell whose quotient rounds to −0 (`Hypsometry`'s
documentation).

## Acceptance criteria

- [ ] `write` returns `InvalidInput` for a world whose elevations span
      more bins than the bound, and writes nothing; a test.
- [ ] The bound and its reason are documented where it is checked.
- [ ] `mise run ci` passes.

## Verification

```sh
mise exec -- cargo test -p kuniumi-io world_dir
mise exec -- cargo test -p kuniumi-sim metrics
```

## Result

_(filled in when the ticket is implemented)_

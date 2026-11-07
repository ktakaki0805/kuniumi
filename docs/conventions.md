# Conventions

Rules for every crate and every change. Some are enforced by `clippy.toml` and
`mise run ci`; the rest by the review loop and tests. Decisions referenced as
D-NN are in `docs/decisions.md`.

## Units

- Inside the code and in every file format, quantities are SI base units:
  metres, seconds, kilograms, kelvin. Angles are radians.
- Geological time is in seconds too (1 Myr ≈ 3.156 × 10¹³ s; `f64` is fine).
- Struct fields and parameters that carry a physical quantity end with their
  unit: `radius_m`, `rotation_period_s`, `temperature_k`,
  `surface_gravity_m_s2`, `obliquity_rad`. Dimensionless values have no suffix.
- Per-cell fields of a `World` are named without a unit (`surface_height`); the
  unit is stored with the values (`Field::unit`).
- Human units (Myr, °C, mm/year, degrees, km) appear only at the boundary:
  settings files, whose keys spell the unit (`obliquity_deg`), and printed
  output. Convert with `kuniumi_math::units`.

## Coordinates

- Planet-centred, right-handed Cartesian frame. +Z is the rotation axis towards
  the north pole; +X passes through longitude 0 on the equator; +Y through
  longitude +90° (east). The planet rotates counter-clockwise seen from +Z.
- The canonical position is a unit vector (`Vec3`). Latitude and longitude are
  derived: radians, latitude in [−π/2, π/2], longitude east-positive in (−π, π].
- For tidally locked planets (v2) the substellar point is longitude 0,
  latitude 0.

## Elevation, sea level and time

- Solid-surface height is measured from the reference radius, in metres.
- Sea level is a separate value on the same reference (`World::sea_level_m`).
  Height above sea level and water depth are derived quantities.
- Ice thickness is a separate field.
- Simulation time starts at 0 s; the end of the simulation is "the present".

## Cells and fields

- Cells are indexed by `u32`. Neighbour lists are CSR (`offsets`, `neighbors`),
  ordered counter-clockwise as seen from outside the planet.
- Cell areas are stored on the unit sphere in steradians; multiply by
  `radius_m²` for square metres.
- Per-cell data is struct-of-arrays: one `Vec<f64>` per quantity.
- Grid iteration is in index order unless an algorithm needs another order (a
  drainage stack, for example), and that order must itself be deterministic.

## Determinism

**Contract (D-08):** on the CPU path, the same seed and the same generator
version produce bit-identical output on every supported platform.

1. The CPU path uses `f64` only; `f32` is banned by clippy. Only
   `kuniumi-viewer` and `kuniumi-gpu` may opt out, with a crate-level
   `#![allow(clippy::disallowed_types)]`, and nothing they compute feeds back
   into generation.
2. Only IEEE 754 correctly rounded operations are used directly: `+ − × ÷`,
   `sqrt`, `mul_add`. Transcendental functions (`sin`, `exp`, `ln`, `atan2`,
   `pow`, …) come from `kuniumi_math::det`, which wraps the pure-Rust `libm`
   crate pinned to an exact version. The std versions are banned by clippy
   because they call the platform's libm.
3. `powi` is banned (unspecified precision); multiply explicitly.
4. `HashMap` and `HashSet` are banned: their iteration order is random per
   process. Use `Vec`, a sorted `Vec`, `BTreeMap` or `BTreeSet`.
5. Sorting by floating-point keys uses `total_cmp` and breaks ties by cell
   index.
6. Floating-point reductions (sums, minima, maxima) happen in a fixed order.
   Parallel code splits work into chunks of a constant size, independent of
   the number of threads, and combines partial results in chunk order. With
   `rayon`, use only order-preserving operations (`par_chunks(N)` with a
   constant `N`, `map`, `collect`, `for_each` on disjoint data), never its
   `sum` or `reduce` on floats. `std::thread::available_parallelism` must not
   influence results (banned).
7. Randomness comes only from `Rng::for_stream(seed, name)`. Stage streams are
   named after `StageId::name()`; sub-streams use a slash
   (`"tectonics/rifting"`). Stream names are part of the output contract.
8. Results never depend on the clock, environment variables, file-system order
   or memory addresses.
9. Stages never produce NaN or infinity; check with `debug_assert!` at stage
   boundaries.

### Golden tests

Every stage gets tests that run it on a small grid with a fixed seed and
compare `kuniumi_io::field_file::fingerprint` of its outputs with recorded
values. A failing golden test means the same seed would produce a different
world: fix the cause, not the expected value. `mise run ci` runs them on the
local machine; once the repository is pushed, CI also runs them on Linux
ARM64, macOS and Windows (D-25).

### World-breaking changes

A change that alters output for an existing seed is *world-breaking*. Before
1.0 such changes are allowed, but each one updates the golden values in the
same commit, explains itself in the commit's `World-breaking:` line, and is
listed under *World-breaking* in `CHANGELOG.md`. Renaming a random stream,
changing `kuniumi_math::rng` and upgrading a pinned dependency on the
generation path are always treated as potentially world-breaking. Raising the
Rust toolchain must not change golden values; if it does, find out why before
going on.

## GPU path (later)

The GPU path uses `f32` and is never used where seed reproducibility is
promised. Its tests check statistical equivalence with the CPU reference: the
validation metrics of the design document ("検証方法") must agree within stated
tolerances.

## Dependencies and tools

- Crates follow D-10: few and well known; licence passes cargo-deny; crates on
  the generation path are pure Rust and pinned to an exact version; each stays
  in the crate that needs it; every one is listed in the D-10 table. Ask the
  owner before adding anything little known or borderline.
- Development tools are pinned in `mise.toml` (D-23). Run checks through
  `mise run` so they use the pinned toolchains.
- Python is allowed only under `tools/`, for development (D-11).
- Ported code must come from permissively licensed sources (MIT, Apache-2.0,
  BSD, public domain), with its origin noted in a comment at the top of the
  ported code and in `THIRD_PARTY.md`. Never port from GPL or non-commercial
  sources.

## Commits and changelog

- A ticket lands on the local `main` as one squashed commit titled
  `<ticket id>: <ticket title>` (D-21). Other commits use a short prefix:
  `docs:`, `ci:`, `fix:`, `chore:`.
- Every commit that changes behaviour updates `CHANGELOG.md` under
  *Unreleased* (D-22).
- Every commit message has a `World-breaking: no` or
  `World-breaking: <reason>` line after its body.
- Nothing is pushed, pulled or fetched except by the owner; nobody force-pushes
  or rewrites the history of `main` (D-21).

## Code style

- Formatted by `mise run format`; `mise run ci` passes (rustfmt check, clippy
  with `-D warnings`, tests in debug and release, MSRV, cargo-deny).
- `unsafe` is forbidden.
- Every public item has a doc comment. Physics code cites its source (paper and
  equation) and documents each parameter with its SI unit and Earth value.
- Tests may use std transcendental functions as a reference; put
  `#[allow(clippy::disallowed_methods)]` on the test module, never on library
  code.

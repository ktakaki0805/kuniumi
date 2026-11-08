# Changelog

All notable changes to kuniumi are recorded here, following
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Nothing is released
before 1.0, so everything sits under *Unreleased*.

Besides Added, Changed, Fixed and Removed, each version has a **World-breaking**
section: every change after which an existing seed produces a different world
(`docs/conventions.md`, "World-breaking changes").

## [Unreleased]

### Added

- Cargo workspace: `kuniumi-math`, `kuniumi-sphere`, `kuniumi-world`,
  `kuniumi-sim`, `kuniumi-io`, `kuniumi-cli`.
- Deterministic random numbers (xoshiro256\*\* with named streams), `Vec3`, SI
  unit helpers.
- `SphereGrid` data structure with full validation and an octahedron test grid.
- `PlanetParams` with Earth defaults; tidally locked and stagnant-lid regimes
  reserved for v2.
- PNG writer and the `.kfield` field format with fingerprints.
- CLI commands `version`, `params`, `generate` (stub) and `test-image`.
- Determinism lints (`clippy.toml`) and cargo-deny configuration.
- Development environment in `mise.toml`: pinned Rust toolchains (development
  and MSRV 1.85.0), cargo-deny, and the tasks `format`, `check` and `ci`.
- CI workflow running the same mise tasks, with a weekly cross-platform
  determinism check.
- Design document, conventions, decisions, roadmap and tickets.
- Claude Code skills `implement-ticket` and `implement-all-tickets`, and the
  subagents `ticket-reviewer` and `ticket-implementer`.
- `kuniumi_math::det`: deterministic `sin`, `cos`, `sin_cos`, `tan`, `atan`,
  `atan2`, `asin`, `acos`, `exp`, `ln` and `pow` over the pure-Rust `libm`
  0.2.16 (pinned exactly), with accuracy, reference golden, special-value and
  output-fingerprint tests (M1-01).
- `tools/det_reference.py` and `tools/det_error.py`, the first development
  scripts (mpmath reference values and error measurement of `det`); Python
  3.14.10 pinned in `mise.toml`, packages in `tools/requirements.txt` (M1-01).
- `kuniumi_sphere::coords`: `from_lat_lon`, `lat_lon` (longitude in (−π, π],
  0 at the poles), the local frame `east` and `north` (defined at the poles as
  the limit along longitude 0), and `great_circle_rad`; and
  `Vec3::rotate_about` (Rodrigues' formula, right-hand rule, length kept
  within 1e-15 relative), all built on `det` and checked by output
  fingerprints; the example `rotate_length` measures the rotation's length
  change (M1-02).
- `kuniumi_sphere::predicates::orient3d` and `Orientation`: the exact sign of
  `det[b − a; c − a; d − a]` (`Positive` for a right-handed tetrahedron, the
  opposite of Shewchuk's sign), over Shewchuk's adaptive predicates from the
  `robust` crate 1.2.0 (pinned exactly); exact for coordinates that are zero or
  of magnitude in [2⁻²⁰⁰, 2²⁰⁰], checked against exact integer evaluations
  (M1-03).

### World-breaking

- None yet.

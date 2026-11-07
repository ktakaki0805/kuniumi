# M1-02 Spherical coordinates, rotations and distances

- Status: ready
- Milestone: M1
- Depends on: M1-01
- Crates: kuniumi-math, kuniumi-sphere
- Design: プロジェクト方針 / 規約の要点
- Decisions: D-09

## Goal

One well-tested place for converting between unit vectors and latitude and
longitude, for local east and north directions, great-circle distances and
rotations about an axis (used by plates, winds and rendering).

## Scope

In:
- `kuniumi_sphere::coords`: `from_lat_lon(lat_rad, lon_rad) -> Vec3`,
  `lat_lon(Vec3) -> (f64, f64)`, `east(Vec3)`, `north(Vec3)`,
  `great_circle_rad(a, b)`.
- `Vec3::rotate_about(axis, angle_rad)` in `kuniumi-math` (Rodrigues formula).

Out:
- Map projections (M1-08, M8-03).

## Implementation notes

- Latitude in [−π/2, π/2]; longitude east-positive in (−π, π]; −π maps to π.
- `lat_lon` uses `atan2(z, √(x² + y²))` and `atan2(y, x)`, which stay accurate
  near the poles (never `asin(z)`).
- At the poles longitude is defined as 0; document how `east` and `north` are
  defined there. The frame must be continuous away from the poles.
- `great_circle_rad(a, b) = atan2(|a × b|, a · b)`: accurate for tiny and
  near-antipodal angles, exactly 0 for `a == b`.
- `rotate_about` requires a unit axis (debug assertion) and calls
  `det::sin_cos` once.

## Acceptance criteria

- [ ] Round trip `lat_lon(from_lat_lon(φ, λ))` within 1e-15 rad over 10⁵ random
      points, including |φ| within 1e-12 of ±π/2 and λ near ±π.
- [ ] `east` and `north` are unit vectors, orthogonal to each other and to the
      position within 1e-15, and `east × north == position` (right-handed).
- [ ] Great-circle distance: symmetric; 0 for identical inputs; π/2 for
      orthogonal axes and π for antipodes within 1e-15; relative error below
      1e-6 for separations of 1e-9 rad.
- [ ] Rotating `Vec3::X` about `Vec3::Z` by π/2 gives `Vec3::Y` within 1e-16;
      rotation preserves length within 1e-15.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo test -p kuniumi-sphere coords
cargo test -p kuniumi-math vec3
```

## Result

_(filled in when the ticket is implemented)_

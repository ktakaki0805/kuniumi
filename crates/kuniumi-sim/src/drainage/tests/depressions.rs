//! Tests of the depression routing (ticket M3-02):
//!
//! - `next_up_by_hand`, the least double above a value;
//! - by hand on the octahedron: `octahedron_pit_by_hand`,
//!   `octahedron_negative_pit_by_hand`, `octahedron_flat_without_ocean_by_hand`
//!   and `octahedron_zero_signs_by_hand`;
//! - synthetic surfaces on the 10,000-cell grid:
//!   `a_bowl_fills_to_its_spill_height`, `the_outlet_is_the_lowest_saddle`
//!   and `nested_bowls_fill_correctly`;
//! - `invariants_on_real_and_random_surfaces`: the fill against the spill
//!   heights of an independent algorithm, unchanged heights where the water
//!   already descends, and routing without pits;
//! - the panics of `fill_depressions` and `depression_depth_m`;
//! - `golden_fingerprint_of_a_filled_earth_world`.
//!
//! Reference values: the spill height `S` of a cell, the lowest over paths
//! to the ocean of the highest height on the path, comes from the
//! iteration of Planchon & Darboux (2002, "A fast, simple and versatile
//! algorithm to fill the depressions of digital elevation models", *Catena*
//! 46, 159–176) with ε = 0 ([`spill_heights_m`]), not from Priority-Flood;
//! the synthetic bowls and the octahedron cases have their levels by
//! construction.
#![allow(clippy::disallowed_methods)]

use super::super::depressions::next_up;
use super::*;

/// The mean spacing of the 10,000-cell grid (rad).
fn spacing_rad() -> f64 {
    mean_spacing_rad(CELLS)
}

/// The unit in the last place of `|x|`: the distance from `|x|` to the next
/// double above it, by its bit pattern (not by [`next_up`]).
fn ulp(x: f64) -> f64 {
    let a = x.abs();
    f64::from_bits(a.to_bits() + 1) - a
}

/// The bound on the ε gradient of a fill of `h`: `2n ULP(max |h|)`. A cell
/// is raised by at most one ULP per cell of its flood path (at most `n`),
/// and the ULP of any filled height is at most twice that of the largest
/// height.
fn epsilon_bound_m(h: &[f64]) -> f64 {
    let largest = h.iter().fold(0.0, |m: f64, &x| m.max(x.abs()));
    2.0 * h.len() as f64 * ulp(largest)
}

/// The lowest cell of `h` by `(height, index)`, heights by `total_cmp`.
fn lowest(h: &[f64]) -> u32 {
    (0..h.len() as u32)
        .min_by(|&a, &b| h[a as usize].total_cmp(&h[b as usize]).then(a.cmp(&b)))
        .expect("cells")
}

/// The spill height of every cell (m): the lowest, over the paths from
/// the cell to an ocean cell (or, without ocean, to the lowest cell), of
/// the highest height on the path, by the iteration of Planchon & Darboux
/// (2002) with ε = 0: `W = h` on the ocean, `+∞` elsewhere; then sweeps,
/// forwards and backwards, of `W_i ← max(h_i, min_j W_j)` while any value
/// falls. `W` only falls and is bounded below by the spill height, which it
/// reaches along every optimal path.
fn spill_heights_m(grid: &SphereGrid, h: &[f64], sea_level_m: f64) -> Vec<f64> {
    let n = grid.len();
    let mut w: Vec<f64> = h
        .iter()
        .map(|&x| if x <= sea_level_m { x } else { f64::INFINITY })
        .collect();
    if w.iter().all(|x| x.is_infinite()) {
        let s = lowest(h) as usize;
        w[s] = h[s];
    }
    let mut sweeps = 0;
    loop {
        let mut changed = false;
        for forward in [true, false] {
            for k in 0..n {
                let i = if forward { k } else { n - 1 - k };
                let m = grid
                    .neighbors(i as u32)
                    .iter()
                    .map(|&j| w[j as usize])
                    .fold(f64::INFINITY, f64::min);
                let v = if h[i] > m { h[i] } else { m };
                if v < w[i] {
                    w[i] = v;
                    changed = true;
                }
            }
        }
        sweeps += 1;
        assert!(sweeps <= n, "no convergence");
        if !changed {
            break;
        }
    }
    w
}

/// Whether each cell has a strictly descending path to an ocean cell (an
/// ocean cell has, trivially) or, without ocean, to the lowest cell, by
/// visiting the cells in ascending height: a land cell has one if a lower
/// neighbour has.
fn strictly_drained(grid: &SphereGrid, h: &[f64], sea_level_m: f64) -> Vec<bool> {
    let n = grid.len();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| h[a].total_cmp(&h[b]).then(a.cmp(&b)));
    let seed = (!h.iter().any(|&x| x <= sea_level_m)).then(|| lowest(h) as usize);
    let mut strict = vec![false; n];
    for i in order {
        strict[i] = h[i] <= sea_level_m
            || seed == Some(i)
            || grid
                .neighbors(i as u32)
                .iter()
                .any(|&j| h[j as usize] < h[i] && strict[j as usize]);
    }
    strict
}

/// What [`check_fill`] found on a surface.
struct Fill {
    filled_m: Vec<f64>,
    spill_m: Vec<f64>,
    eps_m: f64,
    drainage: Drainage,
}

/// Fills `s`, checks every property of the fill and of the routing on it,
/// prints a summary, and returns the fill.
fn check_fill(s: &Surface<'_>) -> Fill {
    let (grid, h, sea, name) = (s.grid, &s.heights_m, s.sea_level_m, &s.name);
    let n = grid.len();
    let filled = fill_depressions(grid, h, sea);
    assert_eq!(filled.len(), n, "{name}");
    let again = fill_depressions(grid, h, sea);
    assert!(
        filled.iter().zip(&again).all(|(&a, &b)| same(a, b)),
        "{name}: not reproducible"
    );
    let spill = spill_heights_m(grid, h, sea);
    let strict = strictly_drained(grid, h, sea);
    let eps = epsilon_bound_m(h);
    let ocean = h.iter().any(|&x| x <= sea);
    let seed = (!ocean).then(|| lowest(h));

    let (mut raised, mut in_depressions, mut on_flats) = (0, 0, 0);
    let (mut flat_raise_m, mut excess_m): (f64, f64) = (0.0, 0.0);
    for i in 0..n {
        let (f, hi) = (filled[i], h[i]);
        assert!(f.is_finite() && f >= hi, "{name}: cell {i}: {f} below {hi}");
        if hi <= sea || seed == Some(i as u32) {
            assert!(same(f, hi), "{name}: base level {i} changed");
        }
        // Criterion 2's levels in general: the fill is the spill height
        // plus at most the ε gradient.
        assert!(
            spill[i] <= f && f - spill[i] <= eps,
            "{name}: cell {i}: filled {f}, spill height {}",
            spill[i]
        );
        // Criterion 3: unchanged, bit for bit, where the water descends.
        if strict[i] {
            assert!(same(f, hi), "{name}: drained cell {i} raised to {f}");
        }
        excess_m = excess_m.max(f - spill[i]);
        if !same(f, hi) {
            raised += 1;
            if spill[i] > hi {
                in_depressions += 1;
            } else {
                on_flats += 1;
                flat_raise_m = flat_raise_m.max(f - hi);
            }
        }
        // Every land cell but the lone seed has a strictly lower neighbour.
        if hi > sea && seed != Some(i as u32) {
            assert!(
                grid.neighbors(i as u32)
                    .iter()
                    .any(|&j| filled[j as usize] < f),
                "{name}: cell {i} has no lower neighbour on the fill"
            );
        }
    }

    let depth = depression_depth_m(h, &filled);
    assert_eq!(depth.len(), n, "{name}");
    for i in 0..n {
        assert!(same(depth[i], filled[i] - h[i]), "{name}: depth {i}");
        assert!(depth[i] >= 0.0, "{name}: depth {i}");
        if same(filled[i], h[i]) {
            assert_eq!(depth[i].to_bits(), 0, "{name}: depth {i} is not +0");
        }
    }

    // Criterion 1: routing on the filled surface leaves no pit (without
    // ocean, the seed), and every land cell reaches an ocean cell (or the
    // seed); the graph is valid on the filled surface (ticket M3-01).
    let d = Drainage::build(grid, &filled, sea, s.radius_m);
    let expected_pits: Vec<u32> = seed.into_iter().collect();
    assert_eq!(d.pits(), expected_pits, "{name}");
    for i in 0..n as u32 {
        let (b, _) = base_of(&d, i);
        assert!(
            h[b as usize] <= sea || Some(b) == seed,
            "{name}: cell {i} ends at {b} on land"
        );
    }
    let routed = Surface {
        name: format!("{name}, filled"),
        grid,
        heights_m: filled.clone(),
        sea_level_m: sea,
        radius_m: s.radius_m,
    };
    check_invariants(&routed, &d);

    let land = h.iter().filter(|&&x| x > sea).count();
    let pits_before = Drainage::build(grid, h, sea, s.radius_m).pits().len();
    let depressed = (0..n).filter(|&i| spill[i] > h[i]).count();
    println!(
        "{name}: {land} land cells, {pits_before} pits before, {} after; {depressed} cells in \
         depressions; raised {raised}: {in_depressions} in depressions, {on_flats} on drained \
         flats (largest raise {flat_raise_m:e} m); largest fill above the spill height \
         {excess_m:e} m, bound {eps:e} m",
        d.pits().len()
    );
    Fill {
        filled_m: filled,
        spill_m: spill,
        eps_m: eps,
        drainage: d,
    }
}

// ---------------------------------------------------------------------------
// `next_up`.

#[test]
fn next_up_by_hand() {
    // IEEE 754 binary64: sign, 11 exponent bits, 52 fraction bits; the
    // next double above a positive value adds one to its pattern, above a
    // negative one subtracts one, and above ±0 it is the least subnormal.
    for (x, expected) in [
        (0x3FF0_0000_0000_0000_u64, 0x3FF0_0000_0000_0001_u64), // 1 → 1 + 2⁻⁵²
        (0xBFF0_0000_0000_0000, 0xBFEF_FFFF_FFFF_FFFF),         // −1 → −(1 − 2⁻⁵³)
        (0x0000_0000_0000_0000, 0x0000_0000_0000_0001),         // +0 → 2⁻¹⁰⁷⁴
        (0x8000_0000_0000_0000, 0x0000_0000_0000_0001),         // −0 → 2⁻¹⁰⁷⁴
        (0x8000_0000_0000_0001, 0x8000_0000_0000_0000),         // −2⁻¹⁰⁷⁴ → −0
        (0x000F_FFFF_FFFF_FFFF, 0x0010_0000_0000_0000),         // subnormal → MIN_POSITIVE
        (0x8010_0000_0000_0000, 0x800F_FFFF_FFFF_FFFF),         // −MIN_POSITIVE → subnormal
        (0x3FFF_FFFF_FFFF_FFFF, 0x4000_0000_0000_0000),         // 2 − 2⁻⁵² → 2
        (0xC000_0000_0000_0000, 0xBFFF_FFFF_FFFF_FFFF),         // −2 → −(2 − 2⁻⁵²)
        (0x7FEF_FFFF_FFFF_FFFF, 0x7FF0_0000_0000_0000),         // f64::MAX → +∞
    ] {
        let got = next_up(f64::from_bits(x)).to_bits();
        assert_eq!(got, expected, "{x:#018X}: {got:#018X}");
    }
    // 200 m lies in [2⁷, 2⁸): one ULP is 2⁻⁴⁵.
    assert!(same(next_up(200.0), 200.0 + 1.0 / 35_184_372_088_832.0));
    assert!(same(next_up(f64::MAX), f64::INFINITY));
}

// ---------------------------------------------------------------------------
// By hand on the octahedron (rings in the M3-01 tests above: 0: [2, 4, 3,
// 5], 1: [2, 5, 3, 4], 2: [4, 0, 5, 1], 3: [4, 1, 5, 0], 4: [0, 2, 1, 3],
// 5: [0, 3, 1, 2]; all centre distances π/2, so the steepest neighbour is
// the lowest).

/// 2⁻⁴⁴, 2⁻⁴⁵ and 2⁻⁴⁶: one ULP in [256, 512), [128, 256) and [64, 128).
const TWO_POW_M44: f64 = 1.0 / 17_592_186_044_416.0;
const TWO_POW_M45: f64 = 1.0 / 35_184_372_088_832.0;
const TWO_POW_M46: f64 = 1.0 / 70_368_744_177_664.0;

/// Fills the octahedron and checks the filled heights, the depths and the
/// receivers on the fill against the values worked out by hand, bit for
/// bit.
fn check_octahedron(
    heights_m: [f64; 6],
    sea_level_m: f64,
    filled_m: [f64; 6],
    depth_m: [f64; 6],
    receivers: [u32; 6],
    pits: &[u32],
) {
    let grid = SphereGrid::octahedron();
    let filled = fill_depressions(&grid, &heights_m, sea_level_m);
    for i in 0..6 {
        assert!(same(filled[i], filled_m[i]), "cell {i}: {}", filled[i]);
    }
    let depth = depression_depth_m(&heights_m, &filled);
    for i in 0..6 {
        assert!(same(depth[i], depth_m[i]), "depth {i}: {}", depth[i]);
    }
    let d = Drainage::build(&grid, &filled, sea_level_m, EARTH_RADIUS_M);
    assert_eq!(d.receivers(), receivers);
    assert_eq!(d.pits(), pits);
}

#[test]
fn octahedron_pit_by_hand() {
    // Sea level 0; cell 5 (−100 m) is the ocean. Cell 4 (10 m) neighbours
    // 0 (300 m), 2 (400 m), 1 (200 m) and 3 (250 m), but not 5: a pit on
    // the original surface. The flood pops 5 and closes 0, 3, 1, 2 at
    // their heights; then pops 1 (200 m, the lowest), which closes 4 one
    // ULP above itself: 200 + 2⁻⁴⁵ m, a depth of 190 + 2⁻⁴⁵ m (exact, both
    // in [128, 256)). On the fill, 4 drains to 1 (its only lower
    // neighbour) and every other land cell to the ocean.
    let h = [300.0, 200.0, 400.0, 250.0, 10.0, -100.0];
    assert_eq!(
        Drainage::build(&SphereGrid::octahedron(), &h, 0.0, EARTH_RADIUS_M).pits(),
        [4]
    );
    check_octahedron(
        h,
        0.0,
        [300.0, 200.0, 400.0, 250.0, 200.0 + TWO_POW_M45, -100.0],
        [0.0, 0.0, 0.0, 0.0, 190.0 + TWO_POW_M45, 0.0],
        [5, 5, 5, 5, 1, 5],
        &[],
    );
}

#[test]
fn octahedron_negative_pit_by_hand() {
    // The same below the datum: sea level −1,000 m, cell 5 (−2,000 m) the
    // ocean, the pit 4 (−600 m) enclosed by 0 (−300 m), 2 (−200 m), 1
    // (−350 m) and 3 (−320 m). It fills to one ULP above −350 m, towards
    // zero: −350 + 2⁻⁴⁴ m, a depth of 250 + 2⁻⁴⁴ m (exact).
    check_octahedron(
        [-300.0, -350.0, -200.0, -320.0, -600.0, -2000.0],
        -1000.0,
        [
            -300.0,
            -350.0,
            -200.0,
            -320.0,
            -350.0 + TWO_POW_M44,
            -2000.0,
        ],
        [0.0, 0.0, 0.0, 0.0, 250.0 + TWO_POW_M44, 0.0],
        [5, 5, 5, 5, 1, 5],
        &[],
    );
}

#[test]
fn octahedron_flat_without_ocean_by_hand() {
    // Every cell at 100 m above a sea level of 0: no ocean, so the lowest
    // cell by (height, index), 0, is the one seed and stays. It closes its
    // ring 2, 4, 3, 5 at 100 + 2⁻⁴⁶ m (one ULP in [64, 128)); 2 pops first
    // of those (lowest index) and closes the antipode 1 one ULP higher, at
    // 100 + 2⁻⁴⁵ m. On the fill, 2, 3, 4, 5 drain to 0, the one pit, and 1
    // ties its four neighbours and takes the lowest index, 2.
    let up1 = 100.0 + TWO_POW_M46;
    let up2 = 100.0 + TWO_POW_M45;
    check_octahedron(
        [100.0; 6],
        0.0,
        [100.0, up2, up1, up1, up1, up1],
        [
            0.0,
            TWO_POW_M45,
            TWO_POW_M46,
            TWO_POW_M46,
            TWO_POW_M46,
            TWO_POW_M46,
        ],
        [0, 2, 0, 0, 0, 0],
        &[0],
    );
}

#[test]
fn octahedron_zero_signs_by_hand() {
    // Sea level −1 m, cell 5 (−2 m) the ocean, cells 0–3 at −2⁻¹⁰⁷⁴ m (the
    // negative subnormal nearest zero) and cell 4 at +0 m. The flood closes
    // 4 from cell 0, whose next double is −0: as +0 is not below −0, cell 4
    // keeps its height, +0 and not −0, and its depth is +0. On the fill,
    // its slopes towards 0–3 underflow to zero, a tie that goes to 0.
    let tiny = -f64::from_bits(1);
    let h = [tiny, tiny, tiny, tiny, 0.0, -2.0];
    check_octahedron(h, -1.0, h, [0.0; 6], [5, 5, 5, 5, 0, 5], &[]);
    let filled = fill_depressions(&SphereGrid::octahedron(), &h, -1.0);
    assert_eq!(filled[4].to_bits(), 0, "the sign of zero changed");
    let depth = depression_depth_m(&h, &filled);
    assert!(depth.iter().all(|x| x.to_bits() == 0), "{depth:?}");
}

// ---------------------------------------------------------------------------
// Synthetic surfaces: bowls around the centre `c`, with θ the angle from c
// and φ the azimuth around it. Rims are three mean spacings wide, more
// than the farthest neighbour of any cell (checked), so no link jumps a
// rim; notches are wedges of the rim, several spacings wide. Beyond the
// outer rim, the land falls from 400 m at 400 m per 0.3 rad of θ to the
// sea at 0 m and below; there, the height falls with θ, so every cell has
// a lower neighbour farther out and the water descends to the ocean.

/// The centre of the bowls.
fn bowl_centre() -> Vec3 {
    Vec3::new(0.6, 0.1, -0.75).normalized().unwrap()
}

/// The angle from the bowls' centre and the azimuth around it of every
/// cell of `grid` (rad).
fn polar(grid: &SphereGrid) -> (Vec<f64>, Vec<f64>) {
    let c = bowl_centre();
    grid.positions()
        .iter()
        .map(|&p| (angle(c, p), azimuth(c, p)))
        .unzip()
}

/// Whether the azimuth `phi` lies within `half_width` of `centre` (rad).
fn in_wedge(phi: f64, centre: f64, half_width: f64) -> bool {
    ((phi - centre + PI).rem_euclid(2.0 * PI) - PI).abs() < half_width
}

/// The land outside the outer rim, whose outer edge is at `edge` (rad).
fn outside_m(theta: f64, edge: f64) -> f64 {
    400.0 * (1.0 - (theta - edge) / 0.3)
}

/// Asserts that the rims hold: the farthest neighbour of every cell is
/// nearer than a rim is wide.
fn check_rims_hold(grid: &SphereGrid) {
    let farthest = (0..grid.len() as u32)
        .map(|i| farthest_neighbour_rad(grid, i))
        .fold(0.0, f64::max);
    assert!(
        farthest < 3.0 * spacing_rad(),
        "farthest neighbour {farthest} rad"
    );
}

/// A surface of the bowl tests on the 10,000-cell grid, sea level 0.
fn bowl_surface(name: &str, heights_m: Vec<f64>) -> Surface<'static> {
    Surface {
        name: name.into(),
        grid: grid(),
        heights_m,
        sea_level_m: 0.0,
        radius_m: EARTH_RADIUS_M,
    }
}

/// The path of receivers from cell `i` to its base level, `i` included.
fn path(d: &Drainage, i: u32) -> Vec<u32> {
    let mut out = vec![i];
    let mut c = i;
    while d.receiver(c) != c {
        c = d.receiver(c);
        out.push(c);
    }
    out
}

/// Asserts that the cells `floor` fill to `level_m` (m): their spill
/// height is `level_m` exactly, and their fill lies within the ε bound
/// above it; returns the largest fill above `level_m` (m).
fn check_level(fill: &Fill, floor: &[u32], level_m: f64, name: &str) -> f64 {
    assert!(!floor.is_empty(), "{name}");
    let mut above_m: f64 = 0.0;
    for &i in floor {
        let (s, f) = (fill.spill_m[i as usize], fill.filled_m[i as usize]);
        assert!(same(s, level_m), "{name}: cell {i}: spill height {s}");
        assert!(
            level_m <= f && f - level_m <= fill.eps_m,
            "{name}: cell {i}: filled to {f}"
        );
        above_m = above_m.max(f - level_m);
    }
    above_m
}

#[test]
fn a_bowl_fills_to_its_spill_height() {
    // Floor 200 m · θ/0.4 for θ < 0.4 rad (0–200 m); a rim of 1,000 m for
    // 0.4 ≤ θ < 0.4 + 3s, lowered to 500 m in a notch |φ| < 0.2 rad (4.4
    // spacings wide); outside it the land falls from 400 m. The only way
    // out of the bowl crosses the rim, at best through the notch, so the
    // spill height of the floor is 500 m.
    let grid = grid();
    check_rims_hold(grid);
    let (theta, phi) = polar(grid);
    let (rim, edge) = (0.4, 0.4 + 3.0 * spacing_rad());
    let notch = |i: usize| theta[i] >= rim && theta[i] < edge && in_wedge(phi[i], 0.0, 0.2);
    let h: Vec<f64> = (0..grid.len())
        .map(|i| match theta[i] {
            t if t < rim => 200.0 * t / rim,
            _ if notch(i) => 500.0,
            t if t < edge => 1000.0,
            t => outside_m(t, edge),
        })
        .collect();
    let s = bowl_surface("bowl", h);
    let fill = check_fill(&s);
    let floor: Vec<u32> = (0..grid.len() as u32)
        .filter(|&i| theta[i as usize] < rim)
        .collect();
    let above_m = check_level(&fill, &floor, 500.0, &s.name);
    for &i in &floor {
        let p = path(&fill.drainage, i);
        assert!(p.iter().any(|&c| notch(c as usize)), "cell {i}: no notch");
        for &c in &p {
            let t = theta[c as usize];
            assert!(
                t < rim || t >= edge || notch(c as usize),
                "cell {i} climbs the rim"
            );
        }
        assert!(s.heights_m[p[p.len() - 1] as usize] <= 0.0, "cell {i}");
    }
    let depth_m3: f64 = floor
        .iter()
        .map(|&i| (fill.filled_m[i as usize] - s.heights_m[i as usize]) * grid.area_sr(i))
        .sum::<f64>()
        * EARTH_RADIUS_M
        * EARTH_RADIUS_M;
    println!(
        "bowl: {} floor cells filled to 500 m + at most {above_m:e} m ({:.0} ULP of 500 m); \
         {depth_m3:.4e} m³ of depressions",
        floor.len(),
        above_m / ulp(500.0)
    );
}

#[test]
fn the_outlet_is_the_lowest_saddle() {
    // The bowl of `a_bowl_fills_to_its_spill_height` with two notches on
    // opposite sides: 500 m at φ = 0 and 700 m at φ = π. The floor fills to
    // the lower, and its water leaves through it, never through the 700 m
    // notch, which lies above the lake.
    let grid = grid();
    check_rims_hold(grid);
    let (theta, phi) = polar(grid);
    let (rim, edge) = (0.4, 0.4 + 3.0 * spacing_rad());
    let in_rim = |i: usize| theta[i] >= rim && theta[i] < edge;
    let low = |i: usize| in_rim(i) && in_wedge(phi[i], 0.0, 0.2);
    let high = |i: usize| in_rim(i) && in_wedge(phi[i], PI, 0.2);
    let h: Vec<f64> = (0..grid.len())
        .map(|i| match theta[i] {
            t if t < rim => 200.0 * t / rim,
            _ if low(i) => 500.0,
            _ if high(i) => 700.0,
            t if t < edge => 1000.0,
            t => outside_m(t, edge),
        })
        .collect();
    let s = bowl_surface("two notches", h);
    let fill = check_fill(&s);
    let floor: Vec<u32> = (0..grid.len() as u32)
        .filter(|&i| theta[i as usize] < rim)
        .collect();
    let above_m = check_level(&fill, &floor, 500.0, &s.name);
    let mut through_low = 0;
    for &i in &floor {
        let p = path(&fill.drainage, i);
        assert!(
            !p.iter().any(|&c| high(c as usize)),
            "cell {i}: 700 m notch"
        );
        assert!(
            p.iter().any(|&c| low(c as usize)),
            "cell {i}: no 500 m notch"
        );
        for &c in &p {
            let t = theta[c as usize];
            assert!(
                t < rim || t >= edge || low(c as usize),
                "cell {i} climbs the rim"
            );
        }
        assert!(s.heights_m[p[p.len() - 1] as usize] <= 0.0, "cell {i}");
        through_low += 1;
    }
    println!(
        "two notches: {through_low} floor cells filled to 500 m + at most {above_m:e} m, all \
         leaving through the 500 m notch"
    );
}

#[test]
fn nested_bowls_fill_correctly() {
    // An inner floor 50 m + 100 m · θ/0.15 for θ < 0.15 rad; an inner rim
    // of `H_in` for 0.15 ≤ θ < 0.15 + 3s, with an optional notch
    // |φ| < 0.5 rad (4.2 spacings); an annulus rising from 100 m to 200 m
    // out to 0.45 rad; the outer rim of 1,000 m for 0.45 ≤ θ < 0.45 + 3s
    // with a 500 m notch on the far side, |φ − π| < 0.2 rad; then the land
    // falling from 400 m. Levels by construction (inner floor, annulus):
    // - inner rim 300 m, no notch: the annulus fills to 500 m and drowns
    //   the inner rim, so both lie at 500 m;
    // - inner rim 800 m, notch 600 m: the inner floor spills at 600 m into
    //   the annulus, which spills at 500 m: 600 m and 500 m, perched;
    // - inner rim 800 m, notch 400 m: the annulus at 500 m drowns the
    //   inner notch, merging the two lakes at 500 m.
    let grid = grid();
    check_rims_hold(grid);
    let (theta, phi) = polar(grid);
    let s3 = 3.0 * spacing_rad();
    let (inner, inner_edge, outer, edge) = (0.15, 0.15 + s3, 0.45, 0.45 + s3);
    let inner_notch =
        |i: usize| theta[i] >= inner && theta[i] < inner_edge && in_wedge(phi[i], 0.0, 0.5);
    let outer_notch = |i: usize| theta[i] >= outer && theta[i] < edge && in_wedge(phi[i], PI, 0.2);
    for (rim_m, notch_m, inner_level_m, annulus_level_m) in [
        (300.0, None, 500.0, 500.0),
        (800.0, Some(600.0), 600.0, 500.0),
        (800.0, Some(400.0), 500.0, 500.0),
    ] {
        let h: Vec<f64> = (0..grid.len())
            .map(|i| match (theta[i], notch_m) {
                (t, _) if t < inner => 50.0 + 100.0 * t / inner,
                (_, Some(n)) if inner_notch(i) => n,
                (t, _) if t < inner_edge => rim_m,
                (t, _) if t < outer => 100.0 + 100.0 * (t - inner_edge) / (outer - inner_edge),
                _ if outer_notch(i) => 500.0,
                (t, _) if t < edge => 1000.0,
                (t, _) => outside_m(t, edge),
            })
            .collect();
        let name = match notch_m {
            None => format!("nested, inner rim {rim_m} m"),
            Some(n) => format!("nested, inner rim {rim_m} m, notch {n} m"),
        };
        let s = bowl_surface(&name, h);
        let fill = check_fill(&s);
        let floor: Vec<u32> = (0..grid.len() as u32)
            .filter(|&i| theta[i as usize] < inner)
            .collect();
        let annulus: Vec<u32> = (0..grid.len() as u32)
            .filter(|&i| (inner_edge..outer).contains(&theta[i as usize]))
            .collect();
        let inner_above_m = check_level(&fill, &floor, inner_level_m, &s.name);
        let annulus_above_m = check_level(&fill, &annulus, annulus_level_m, &s.name);
        for &i in &floor {
            let p = path(&fill.drainage, i);
            assert!(
                p.iter().any(|&c| outer_notch(c as usize)),
                "{}: cell {i}",
                s.name
            );
            if notch_m.is_some() {
                assert!(
                    p.iter().any(|&c| inner_notch(c as usize)),
                    "{}: cell {i}",
                    s.name
                );
            }
            assert!(
                s.heights_m[p[p.len() - 1] as usize] <= 0.0,
                "{}: cell {i}",
                s.name
            );
        }
        println!(
            "{}: {} inner cells at {inner_level_m} m + at most {inner_above_m:e} m, {} annulus \
             cells at {annulus_level_m} m + at most {annulus_above_m:e} m",
            s.name,
            floor.len(),
            annulus.len()
        );
    }
}

// ---------------------------------------------------------------------------
// Invariants.

/// A smooth random surface: a sum of 24 plane waves of random direction,
/// wavenumber 2–12 per radian and amplitude 1,000 m over the wavenumber.
fn waves_m(grid: &SphereGrid) -> Vec<f64> {
    let mut rng = Rng::for_stream(5, "test/drainage/waves");
    let waves: Vec<(Vec3, f64, f64, f64)> = (0..24)
        .map(|_| {
            let z = rng.range_f64(-1.0, 1.0);
            let lon = rng.range_f64(-PI, PI);
            let r = (1.0 - z * z).sqrt();
            let dir = Vec3::new(r * lon.cos(), r * lon.sin(), z);
            let k = rng.range_f64(2.0, 12.0);
            (dir, k, rng.range_f64(0.0, 2.0 * PI), 1000.0 / k)
        })
        .collect();
    grid.positions()
        .iter()
        .map(|&p| {
            waves
                .iter()
                .map(|&(w, k, phase, a)| a * (k * p.dot(w) + phase).sin())
                .sum()
        })
        .collect()
}

#[test]
fn invariants_on_real_and_random_surfaces() {
    // The surfaces of the M3-01 invariants (the Earth worlds of seeds 1–3,
    // uniform random heights, random steps, flat land without ocean and
    // ocean everywhere), a smooth random surface with sea level −150 m, and
    // the same without ocean (sea level 1 m below its lowest cell), whose
    // one seed is not cell 0: `check_fill` asserts that the lowest cell is
    // the one pit and every cell drains to it.
    let mut all = surfaces();
    let waves = waves_m(grid());
    let lowest_m = waves.iter().fold(f64::INFINITY, |m, &h| m.min(h));
    for (name, sea_level_m) in [
        ("smooth random", -150.0),
        ("smooth random without ocean", lowest_m - 1.0),
    ] {
        all.push(Surface {
            name: name.into(),
            grid: grid(),
            heights_m: waves.clone(),
            sea_level_m,
            radius_m: EARTH_RADIUS_M,
        });
    }
    for s in &all {
        let fill = check_fill(s);
        match s.name.as_str() {
            "flat land" => assert_eq!(fill.drainage.pits(), [0]),
            "smooth random without ocean" => {
                let seed = lowest(&s.heights_m);
                assert_ne!(seed, 0);
                assert_eq!(fill.drainage.pits(), [seed]);
            }
            "ocean everywhere" => assert!(
                fill.filled_m
                    .iter()
                    .zip(&s.heights_m)
                    .all(|(&a, &b)| same(a, b))
            ),
            _ => {}
        }
    }
}

// ---------------------------------------------------------------------------
// Panics.

#[test]
#[should_panic(expected = "one height per cell")]
fn fill_needs_one_height_per_cell() {
    let _ = fill_depressions(&SphereGrid::octahedron(), &[0.0; 7], 0.0);
}

#[test]
fn fill_needs_finite_heights() {
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let h = [0.0, 1.0, 2.0, bad, 4.0, 5.0];
        let result =
            std::panic::catch_unwind(|| fill_depressions(&SphereGrid::octahedron(), &h, 0.0));
        let message = result.expect_err("a panic");
        let text = message
            .downcast_ref::<String>()
            .map(String::as_str)
            .unwrap_or_default();
        assert!(text.contains("finite heights"), "{bad}: {text}");
    }
}

#[test]
#[should_panic(expected = "finite sea level")]
fn fill_needs_a_finite_sea_level() {
    let _ = fill_depressions(&SphereGrid::octahedron(), &[0.0; 6], f64::NAN);
}

#[test]
#[should_panic(expected = "next double is infinite")]
fn fill_needs_room_above_its_flats() {
    // Cell 5 is the ocean; 0–3 are closed at f64::MAX, and 4, behind them,
    // would rise one ULP above it.
    let h = [f64::MAX, f64::MAX, f64::MAX, f64::MAX, f64::MAX, -1.0];
    let _ = fill_depressions(&SphereGrid::octahedron(), &h, 0.0);
}

#[test]
#[should_panic(expected = "one filled height per height")]
fn depression_depth_needs_one_filled_height_per_height() {
    let _ = depression_depth_m(&[0.0; 6], &[0.0; 5]);
}

// ---------------------------------------------------------------------------
// Golden values.

#[test]
fn golden_fingerprint_of_a_filled_earth_world() {
    // The Earth world of seed 1 at 10,000 cells: FNV-1a over the bits of
    // the filled heights (f64, little-endian) and the receivers and the
    // stack of the graph on them (u32, little-endian). A change means that
    // the same world fills or drains differently. Recorded after the tests
    // above passed and the maps of out/M3-02 were viewed (ticket M3-02,
    // Result).
    let world = &worlds()[0];
    let heights = &world.field(SURFACE_HEIGHT).expect("heights").values;
    let filled = fill_depressions(&world.grid, heights, world.sea_level_m);
    let d = Drainage::build(
        &world.grid,
        &filled,
        world.sea_level_m,
        world.params.radius_m,
    );
    let bytes: Vec<u8> = filled
        .iter()
        .flat_map(|h| h.to_bits().to_le_bytes())
        .chain(
            d.receivers()
                .iter()
                .chain(d.stack())
                .flat_map(|v| v.to_le_bytes()),
        )
        .collect();
    let fingerprint = fnv1a64(&bytes);
    println!("fingerprint {fingerprint:#018X}");
    assert_eq!(fingerprint, GOLDEN_FILLED_FINGERPRINT);
}

/// The fingerprint of `golden_fingerprint_of_a_filled_earth_world`.
const GOLDEN_FILLED_FINGERPRINT: u64 = 0x5119_A7A4_1C16_BDAA;

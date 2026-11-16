//! Tests of the multigrid (ticket M4-09): the aggregates, the prolongator
//! and the Galerkin products against their definitions computed here with
//! dense matrices; solutions of random M-matrices and of the grid's systems
//! against Gaussian elimination done here; the iterations against the size
//! of the grid and over the planet range of `D`; the symmetry of the
//! V-cycle; the exact and degenerate cases; and the panics.
#![allow(clippy::disallowed_methods)]

use std::sync::OnceLock;

use kuniumi_math::Rng;
use kuniumi_sphere::SphereGrid;

use super::*;
use crate::operators::Laplacian;

/// The dense form of `m`.
fn dense_of(m: &SparseMatrix) -> Vec<Vec<f64>> {
    let n = m.len();
    let mut dense = vec![vec![0.0; n]; n];
    for (i, row) in dense.iter_mut().enumerate() {
        row[i] = m.diagonal()[i];
        let (columns, values) = m.row(i);
        for (&j, &v) in columns.iter().zip(values) {
            row[j as usize] = v;
        }
    }
    dense
}

/// The dense form of the rows `p` (`p.len()` rows of `p.width()` columns).
fn dense_rows(p: &Rows) -> Vec<Vec<f64>> {
    let mut dense = vec![vec![0.0; p.width()]; p.len()];
    for (i, row) in dense.iter_mut().enumerate() {
        let (columns, values) = p.row(i);
        for (&j, &v) in columns.iter().zip(values) {
            row[j as usize] = v;
        }
    }
    dense
}

/// Gaussian elimination with partial pivoting of the dense `a x = b`.
fn dense_solve(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Vec<f64> {
    let n = b.len();
    for col in 0..n {
        let pivot = (col..n)
            .max_by(|&p, &q| a[p][col].abs().total_cmp(&a[q][col].abs()))
            .expect("a row");
        a.swap(col, pivot);
        b.swap(col, pivot);
        for row in col + 1..n {
            let factor = a[row][col] / a[col][col];
            if factor == 0.0 {
                continue;
            }
            let (upper, lower) = a.split_at_mut(row);
            for (x, y) in lower[0][col..].iter_mut().zip(&upper[col][col..]) {
                *x -= factor * y;
            }
            b[row] -= factor * b[col];
        }
    }
    let mut x = vec![0.0; n];
    for row in (0..n).rev() {
        let mut s = b[row];
        for k in row + 1..n {
            s -= a[row][k] * x[k];
        }
        x[row] = s / a[row][row];
    }
    x
}

/// A random symmetric M-matrix of size `n`: a ring and about three random
/// neighbours per row with negative weights over two decades, and a
/// diagonal of the weights' sum plus a mass between `mass_lo` and
/// `mass_hi` times the row's typical weight.
fn random_spd(n: usize, mass_lo: f64, mass_hi: f64, rng: &mut Rng) -> SparseMatrix {
    let mut pairs: Vec<(usize, usize, f64)> = Vec::new();
    for i in 0..n {
        let mut columns = vec![(i + 1) % n];
        for _ in 0..2 {
            columns.push(rng.below(n as u64) as usize);
        }
        for j in columns {
            if j != i {
                pairs.push((
                    i.min(j),
                    i.max(j),
                    (10.0_f64).powf(rng.range_f64(-1.0, 1.0)),
                ));
            }
        }
    }
    let mut dense = vec![vec![0.0; n]; n];
    for &(i, j, w) in &pairs {
        dense[i][j] -= w;
        dense[j][i] -= w;
    }
    let mut diagonal = vec![0.0; n];
    let mut rows = vec![Vec::new(); n];
    for i in 0..n {
        let mut sum = 0.0;
        for (j, &v) in dense[i].iter().enumerate() {
            if v != 0.0 {
                rows[i].push((j as u32, v));
                sum -= v;
            }
        }
        diagonal[i] = sum + rng.range_f64(mass_lo, mass_hi);
    }
    SparseMatrix::from_rows(diagonal, &rows)
}

/// `symmetric` plus a random upwind advection: on each edge one direction
/// carries a random flow `q ≥ 0` out of the cell, `+q` on its diagonal and
/// `−q` on the receiving row, so the result is a non-symmetric M-matrix
/// whose columns are dominant (as the sea surface's, ticket M4-05).
fn with_advection(symmetric: &SparseMatrix, strength: f64, rng: &mut Rng) -> SparseMatrix {
    let n = symmetric.len();
    let mut dense = dense_of(symmetric);
    for i in 0..n {
        let (columns, _) = symmetric.row(i);
        for &j in columns {
            let j = j as usize;
            if i < j {
                let q = strength * rng.next_f64();
                let (from, to) = if rng.below(2) == 0 { (i, j) } else { (j, i) };
                dense[from][from] += q;
                dense[to][from] -= q;
            }
        }
    }
    let mut diagonal = vec![0.0; n];
    let mut rows = vec![Vec::new(); n];
    for i in 0..n {
        diagonal[i] = dense[i][i];
        let (columns, _) = symmetric.row(i);
        for &j in columns {
            rows[i].push((j, dense[i][j as usize]));
        }
    }
    SparseMatrix::from_rows(diagonal, &rows)
}

/// The year (s) of the Earth-like sub-step: a sidereal year.
const YEAR_S: f64 = 3.155_814_954e7;

/// The heat capacities of land and ocean (J m⁻² K⁻¹) and `B`
/// (W m⁻² K⁻¹), Earth's values of ticket M4-02.
const LAND_J_M2_K: f64 = 1.037e7;
const OCEAN_J_M2_K: f64 = 3.171e8;
const B_W_M2_K: f64 = 2.09;

/// Earth's `D` (W m⁻² K⁻¹) of ticket M4-02.
const D_EARTH: f64 = 0.583;

/// Whether the cell at `p` is land on the test planet: two continents.
fn land(p: kuniumi_math::Vec3) -> bool {
    p.x * p.y + 0.6 * p.z * p.x > 0.15
}

/// The masses `a_i (C_i/Δt + B)` of the energy balance's sub-step on
/// `grid`, with the land of [`land`].
fn substep_mass(grid: &SphereGrid) -> Vec<f64> {
    let dt = YEAR_S / 48.0;
    (0..grid.len() as u32)
        .map(|i| {
            let c = if land(grid.position(i)) {
                LAND_J_M2_K
            } else {
                OCEAN_J_M2_K
            };
            grid.area_sr(i) * (c / dt + B_W_M2_K)
        })
        .collect()
}

/// The energy balance's sub-step matrix on `grid` with the diffusion `d`.
fn substep_matrix(grid: &SphereGrid, d: f64) -> SparseMatrix {
    Laplacian::new(grid).matrix(&substep_mass(grid), d)
}

/// A right-hand side like the sub-step's: `a_i` times a smooth field of
/// a few thousand W/m² with a land–sea step.
fn substep_rhs(grid: &SphereGrid) -> Vec<f64> {
    (0..grid.len() as u32)
        .map(|i| {
            let p = grid.position(i);
            let t = 288.0 - 28.0 * (1.5 * p.z * p.z - 0.5) - if land(p) { 5.0 } else { 0.0 };
            grid.area_sr(i) * (500.0 * t + 100.0 * p.x)
        })
        .collect()
}

/// The grid of `cells` cells, built once per size.
fn grid(cells: u32) -> &'static SphereGrid {
    const SIZES: [u32; 5] = [1_000, 2_000, 2_500, 10_000, 40_000];
    static GRIDS: [OnceLock<SphereGrid>; 5] = [const { OnceLock::new() }; 5];
    let k = SIZES.iter().position(|&s| s == cells).expect("a test size");
    GRIDS[k].get_or_init(|| SphereGrid::fibonacci(cells, 1).expect("grid"))
}

/// The largest `|a_i − b_i|` and the largest `|b_i|`.
fn error(a: &[f64], b: &[f64]) -> (f64, f64) {
    let mut worst = 0.0_f64;
    let mut largest = 0.0_f64;
    for (x, y) in a.iter().zip(b) {
        worst = worst.max((x - y).abs());
        largest = largest.max(y.abs());
    }
    (worst, largest)
}

#[test]
fn laplacian_matrix_is_the_system_of_solve() {
    // `Laplacian::matrix(m, c)` is `M + c L`: its product with a random
    // vector equals `m_i u_i + c Σ_j w_ij (u_i − u_j)` computed here from
    // the weights, its diagonal is `m_i + c Σ_j w_ij` and it is symmetric
    // bit for bit.
    let g = grid(1_000);
    let lap = Laplacian::new(g);
    let mut rng = Rng::for_stream(1, "test/multigrid/matrix");
    let mass: Vec<f64> = (0..g.len()).map(|_| rng.range_f64(0.1, 1.0)).collect();
    let u: Vec<f64> = (0..g.len()).map(|_| rng.range_f64(-1.0, 1.0)).collect();
    let m = lap.matrix(&mass, 0.7);
    assert!(is_symmetric(&m));
    let mut y = vec![0.0; g.len()];
    m.apply(&u, &mut y);
    for i in 0..g.len() {
        let w = lap.weights(i as u32);
        let mut flux = 0.0;
        let mut sum = 0.0;
        for (k, &j) in g.neighbors(i as u32).iter().enumerate() {
            flux += w[k] * (u[i] - u[j as usize]);
            sum += w[k];
        }
        assert_eq!(m.diagonal()[i].to_bits(), (mass[i] + 0.7 * sum).to_bits());
        let expected = mass[i] * u[i] + 0.7 * flux;
        assert!((y[i] - expected).abs() <= 1e-14 * (1.0 + expected.abs()));
    }
}

/// Whether the entry of row `i` at column `j` of `m` is strong.
fn strong(m: &SparseMatrix, i: usize, j: usize) -> bool {
    let (columns, values) = m.row(i);
    let k = columns.binary_search(&(j as u32)).expect("a neighbour");
    let d = m.diagonal();
    values[k] * values[k] > STRENGTH_THRESHOLD * STRENGTH_THRESHOLD * d[i] * d[j]
}

#[test]
fn aggregates_cover_the_strongly_connected_cells() {
    // The Earth-like sub-step on 10,000 cells: diffusion dominates every
    // row, so every cell is in an aggregate, each aggregate is connected
    // by strong edges and has at least two cells, and there are between
    // n/15 and n/3 of them. With a coupling of 10⁻⁹ the masses dominate:
    // no cell has a strong neighbour and there is no aggregate.
    let g = grid(10_000);
    let m = substep_matrix(g, D_EARTH);
    let aggregates = aggregate(&m, STRENGTH_THRESHOLD);
    let n = m.len();
    assert!(
        aggregates.count > n / 15 && aggregates.count < n / 3,
        "{}",
        aggregates.count
    );
    let mut members = vec![Vec::new(); aggregates.count];
    for (i, &a) in aggregates.cell.iter().enumerate() {
        assert_ne!(a, NO_AGGREGATE, "cell {i} is in no aggregate");
        members[a as usize].push(i);
    }
    for (a, cells) in members.iter().enumerate() {
        assert!(cells.len() >= 2, "aggregate {a} has {} cells", cells.len());
        // Connected through strong edges inside the aggregate.
        let mut seen = vec![cells[0]];
        let mut k = 0;
        while k < seen.len() {
            let i = seen[k];
            for &j in m.row(i).0 {
                let j = j as usize;
                if aggregates.cell[j] == a as u32 && strong(&m, i, j) && !seen.contains(&j) {
                    seen.push(j);
                }
            }
            k += 1;
        }
        assert_eq!(seen.len(), cells.len(), "aggregate {a} is not connected");
    }
    let weak = substep_matrix(g, 1e-9);
    let none = aggregate(&weak, STRENGTH_THRESHOLD);
    assert_eq!(none.count, 0);
    assert!(none.cell.iter().all(|&a| a == NO_AGGREGATE));
}

#[test]
fn the_prolongator_is_the_smoothed_aggregation() {
    // For the graph Laplacian without mass (every entry strong, nothing
    // filtered) `L 1 = 0`, so the smoothed prolongator keeps the constant:
    // every row of P sums to 1 (to rounding). For the Earth-like sub-step
    // on 2,000 cells, where the ocean's mass makes some entries weak, P is
    // `(I − ω D_F⁻¹ A_F) P₀` computed here row by row: A_F keeps the strong
    // entries and adds the weak ones to the diagonal, D_F its diagonal,
    // `ω = 4/(3ρ)` with ρ the largest `Σ_j |a^F_ij| / a^F_ii`; its rows sum
    // to between 0 and 1.
    let g = grid(2_000);
    let lap = Laplacian::new(g);
    let zero = lap.matrix(&vec![0.0; g.len()], 1.0);
    let aggregates = aggregate(&zero, STRENGTH_THRESHOLD);
    let p = prolongator(&zero, STRENGTH_THRESHOLD, &aggregates);
    assert_eq!(p.len(), g.len());
    assert_eq!(p.width(), aggregates.count);
    for i in 0..p.len() {
        let sum: f64 = p.row(i).1.iter().sum();
        assert!((sum - 1.0).abs() <= 1e-14, "row {i} sums to {sum}");
    }
    let m = substep_matrix(g, D_EARTH);
    let n = m.len();
    let aggregates = aggregate(&m, STRENGTH_THRESHOLD);
    let p = prolongator(&m, STRENGTH_THRESHOLD, &aggregates);
    let d = m.diagonal();
    let mut weak_entries = 0;
    let mut filtered = vec![0.0; n];
    let mut rho = 0.0_f64;
    for i in 0..n {
        let (columns, values) = m.row(i);
        let (mut f, mut s) = (d[i], 0.0);
        for (&j, &v) in columns.iter().zip(values) {
            if strong(&m, i, j as usize) {
                s += v.abs();
            } else {
                f += v;
                weak_entries += 1;
            }
        }
        filtered[i] = f;
        rho = rho.max((f + s) / f);
    }
    assert!(weak_entries > 0, "the test needs weak entries");
    let omega = 4.0 / (3.0 * rho);
    for i in 0..n {
        let mut expected = vec![0.0; p.width()];
        if aggregates.cell[i] != NO_AGGREGATE {
            expected[aggregates.cell[i] as usize] += 1.0 - omega;
        }
        let (columns, values) = m.row(i);
        for (&j, &v) in columns.iter().zip(values) {
            let a = aggregates.cell[j as usize];
            if strong(&m, i, j as usize) && a != NO_AGGREGATE {
                expected[a as usize] -= omega * v / filtered[i];
            }
        }
        let mut got = vec![0.0; p.width()];
        let (columns, values) = p.row(i);
        for (&c, &v) in columns.iter().zip(values) {
            got[c as usize] = v;
        }
        for c in 0..p.width() {
            assert!(
                (got[c] - expected[c]).abs() <= 1e-14,
                "P[{i}][{c}] = {} against {}",
                got[c],
                expected[c]
            );
        }
        let sum: f64 = values.iter().sum();
        assert!((0.0..=1.0 + 1e-14).contains(&sum), "row {i} sums to {sum}");
    }
}

#[test]
fn coarse_operators_are_galerkin() {
    // Pᵀ A P against the dense product, for a random symmetric M-matrix
    // and for the same with advection (non-symmetric, the prolongator of
    // the symmetric one), within 10⁻¹⁴ of the largest entry; the
    // transpose of P is P's transpose exactly; the symmetric product is
    // symmetric to 10⁻¹⁴ of the largest entry (rounding).
    let mut rng = Rng::for_stream(2, "test/multigrid/galerkin");
    let s = random_spd(300, 0.01, 0.1, &mut rng);
    let a = with_advection(&s, 2.0, &mut rng);
    let aggregates = aggregate(&s, STRENGTH_THRESHOLD);
    let p = prolongator(&s, STRENGTH_THRESHOLD, &aggregates);
    let r = p.transpose();
    let dp = dense_rows(&p);
    let dr = dense_rows(&r);
    assert_eq!(r.len(), p.width());
    assert_eq!(r.width(), p.len());
    for (i, row) in dp.iter().enumerate() {
        for (c, &v) in row.iter().enumerate() {
            assert_eq!(v.to_bits(), dr[c][i].to_bits());
        }
    }
    for (matrix, symmetric) in [(&s, true), (&a, false)] {
        let coarse = galerkin(matrix, &p, &r);
        let da = dense_of(matrix);
        let nc = p.width();
        let n = matrix.len();
        let mut ap = vec![vec![0.0; nc]; n];
        for i in 0..n {
            for k in 0..n {
                if da[i][k] != 0.0 {
                    for c in 0..nc {
                        ap[i][c] += da[i][k] * dp[k][c];
                    }
                }
            }
        }
        let dc = dense_of(&coarse);
        let mut largest = 0.0_f64;
        let mut expected = vec![vec![0.0; nc]; nc];
        for c in 0..nc {
            for d in 0..nc {
                for i in 0..n {
                    expected[c][d] += dp[i][c] * ap[i][d];
                }
                largest = largest.max(expected[c][d].abs());
            }
        }
        for c in 0..nc {
            for d in 0..nc {
                assert!(
                    (dc[c][d] - expected[c][d]).abs() <= 1e-14 * largest,
                    "({c}, {d}): {} against {}",
                    dc[c][d],
                    expected[c][d]
                );
                if symmetric {
                    assert!((dc[c][d] - dc[d][c]).abs() <= 1e-14 * largest);
                }
            }
        }
    }
}

#[test]
fn solves_match_a_dense_solution() {
    // Random symmetric M-matrices of 5 to 900 rows (one level up to 400,
    // several above) and the Earth-like sub-step on 1,000 cells, by the
    // conjugate gradient method; the same with advection by BiCGSTAB, the
    // prolongators of the symmetric matrix. To a relative residual of
    // 10⁻¹³: within 10⁻⁹ of the largest value of the dense solution.
    let mut rng = Rng::for_stream(3, "test/multigrid/dense");
    let g = grid(1_000);
    let mut systems: Vec<SparseMatrix> = [5, 50, 400, 900]
        .iter()
        .map(|&n| random_spd(n, 0.001, 0.1, &mut rng))
        .collect();
    systems.push(substep_matrix(g, D_EARTH));
    for s in systems {
        let n = s.len();
        let b: Vec<f64> = (0..n).map(|_| rng.range_f64(-1.0, 1.0)).collect();
        let mg = Multigrid::new(s.clone());
        let mut x = vec![0.0; n];
        let solve = conjugate_gradient(&s, &mg, &b, 1e-13, 1_000, &mut x);
        assert!(solve.converged(), "n = {n}: {solve:?}");
        let (worst, largest) = error(&x, &dense_solve(dense_of(&s), b.clone()));
        assert!(worst <= 1e-9 * largest, "n = {n}: {worst:e} of {largest:e}");

        let a = with_advection(&s, 1.0, &mut rng);
        let ma = Multigrid::with_prolongators_of(a.clone(), &mg);
        assert_eq!(ma.levels(), mg.levels());
        if let Some(level) = mg.levels.first() {
            // The coarse operator of the non-symmetric hierarchy is
            // `Pᵀ A P` of the non-symmetric matrix, not the symmetric one's.
            let coarse = ma.levels.get(1).map_or(&ma.coarsest, |l| &l.matrix);
            let expected = galerkin(&a, &level.prolongator, &level.restriction);
            assert_eq!(coarse, &expected, "n = {n}");
            assert!(!is_symmetric(coarse), "n = {n}");
        }
        let mut x = vec![0.0; n];
        let solve = bicgstab(&a, &ma, &b, 1e-13, 1_000, &mut x);
        assert!(solve.converged(), "n = {n}: {solve:?}");
        let (worst, largest) = error(&x, &dense_solve(dense_of(&a), b.clone()));
        assert!(
            worst <= 1e-9 * largest,
            "n = {n} with advection: {worst:e} of {largest:e}"
        );
    }
}

/// The iterations of the conjugate gradient method with the multigrid for
/// the matrix `m` and the right-hand side `b` to 10⁻¹⁰.
fn iterations(m: &SparseMatrix, b: &[f64]) -> u32 {
    let mg = Multigrid::new(m.clone());
    let mut x = vec![0.0; m.len()];
    let solve = conjugate_gradient(m, &mg, b, 1e-10, 1_000, &mut x);
    assert!(solve.converged(), "{solve:?}");
    solve.iterations
}

#[test]
fn iterations_do_not_grow_with_the_grid() {
    // The energy balance's sub-step with Earth's D on 2,500, 10,000 and
    // 40,000 cells: at most 20 iterations, within 3 of each other (the
    // Jacobi preconditioner needs about twice as many each time the cells
    // quadruple; ticket M4-02).
    let counts: Vec<u32> = [2_500, 10_000, 40_000]
        .iter()
        .map(|&cells| {
            let g = grid(cells);
            iterations(&substep_matrix(g, D_EARTH), &substep_rhs(g))
        })
        .collect();
    println!("iterations on 2,500, 10,000 and 40,000 cells: {counts:?}");
    assert!(counts.iter().all(|&c| c <= 20), "{counts:?}");
    let (lo, hi) = (counts.iter().min().unwrap(), counts.iter().max().unwrap());
    assert!(hi - lo <= 3, "{counts:?}");
}

#[test]
fn robust_over_the_planet_range() {
    // 10,000 cells: the sub-step with D from D⊕/91 (8 h, 0.1 bar) to
    // 990 D⊕ (240 h, 10 bar; ticket M4-02) and the warm start's
    // `B a + D L` within 30 iterations; the moisture's potential
    // `10⁻⁹ a + L` (ticket M4-06) within 40.
    let g = grid(10_000);
    let b = substep_rhs(g);
    let mut counts = Vec::new();
    for d in [D_EARTH / 91.0, D_EARTH, 990.0 * D_EARTH] {
        counts.push(iterations(&substep_matrix(g, d), &b));
    }
    let lap = Laplacian::new(g);
    let areas: Vec<f64> = (0..g.len() as u32).map(|i| g.area_sr(i)).collect();
    let warm: Vec<f64> = areas.iter().map(|a| a * B_W_M2_K).collect();
    counts.push(iterations(&lap.matrix(&warm, D_EARTH), &b));
    println!("iterations, D⊕/91, D⊕, 990 D⊕ and the warm start: {counts:?}");
    assert!(counts.iter().all(|&c| c <= 30), "{counts:?}");
    let potential: Vec<f64> = areas.iter().map(|a| 1e-9 * a).collect();
    // A divergence: zero mean up to rounding.
    let raw: Vec<f64> = (0..g.len() as u32)
        .map(|i| {
            let p = g.position(i);
            g.area_sr(i) * (1.5 * p.z * p.z - 0.5 + 0.3 * p.x * p.y)
        })
        .collect();
    let total: f64 = raw.iter().sum();
    let area: f64 = areas.iter().sum();
    let rhs: Vec<f64> = raw
        .iter()
        .zip(&areas)
        .map(|(r, a)| r - a * total / area)
        .collect();
    let count = iterations(&lap.matrix(&potential, 1.0), &rhs);
    println!("iterations of the potential: {count}");
    assert!(count <= 40, "{count}");
}

#[test]
fn a_hierarchy_serves_a_matrix_close_to_its_own() {
    // The sea surface solves twelve months with the hierarchy of the
    // annual-mean matrix (ticket M4-05): the Earth-like sub-step on 10,000
    // cells with its diagonal raised by up to 10 % at random, solved with
    // the hierarchy of the unchanged matrix, by the conjugate gradient
    // method and by BiCGSTAB, within 30 iterations to 10⁻¹⁰ (12 to 15 with
    // its own hierarchy, `iterations_do_not_grow_with_the_grid`) and within
    // 10⁻⁹ of the solution with its own hierarchy. Without levels (a
    // diagonal matrix of 2,000 rows) the sweeps of the matrix being solved
    // give its solution in one iteration (to rounding).
    let g = grid(10_000);
    let m = substep_matrix(g, D_EARTH);
    let mut rng = Rng::for_stream(13, "test/multigrid/close");
    let diagonal: Vec<f64> = m
        .diagonal()
        .iter()
        .map(|d| d * (1.0 + 0.1 * rng.next_f64()))
        .collect();
    let close = m.with_diagonal(diagonal);
    let b = substep_rhs(g);
    let own = Multigrid::new(close.clone());
    let mut exact = vec![0.0; g.len()];
    conjugate_gradient(&close, &own, &b, 1e-12, 100, &mut exact);
    let mg = Multigrid::new(m.clone());
    let mut x = vec![0.0; g.len()];
    let solve = conjugate_gradient(&close, &mg, &b, 1e-10, 100, &mut x);
    println!("conjugate gradient: {solve:?}");
    assert!(solve.converged() && solve.iterations <= 30, "{solve:?}");
    let (worst, largest) = error(&x, &exact);
    assert!(worst <= 1e-9 * largest, "{worst:e}");
    let ma = Multigrid::with_prolongators_of(m.clone(), &mg);
    let mut x = vec![0.0; g.len()];
    let solve = bicgstab(&close, &ma, &b, 1e-10, 100, &mut x);
    println!("BiCGSTAB: {solve:?}");
    assert!(solve.converged() && solve.iterations <= 30, "{solve:?}");
    let (worst, largest) = error(&x, &exact);
    assert!(worst <= 1e-9 * largest, "{worst:e}");

    let n = 2_000;
    let d = SparseMatrix::from_rows(vec![2.0; n], &vec![Vec::new(); n]);
    let md = Multigrid::new(d);
    assert_eq!(md.levels(), 1);
    let other = SparseMatrix::from_rows(
        (0..n).map(|i| 1.0 + i as f64).collect(),
        &vec![Vec::new(); n],
    );
    let mut x = vec![0.0; n];
    let solve = conjugate_gradient(&other, &md, &vec![1.0; n], 1e-12, 10, &mut x);
    assert_eq!(solve.iterations, 1, "{solve:?}");
    assert!(
        x.iter()
            .enumerate()
            .all(|(i, &v)| (v * (1.0 + i as f64) - 1.0).abs() <= 1e-14)
    );
}

#[test]
fn the_cycle_smooths_the_finest_level_with_the_matrix_being_solved() {
    // A hierarchy of several levels used for a matrix whose diagonal
    // differs (as the sea surface's months do): the cycle equals, bit for
    // bit, a forward sweep of that matrix from zero, its residual
    // restricted, the cycle of the hierarchy's coarser levels, the
    // correction prolonged and a backward sweep of that matrix, composed
    // here; and it differs from the cycle of the hierarchy's own matrix.
    let g = grid(10_000);
    let m = substep_matrix(g, D_EARTH);
    let mg = Multigrid::new(m.clone());
    assert!(mg.levels() >= 3, "{}", mg.levels());
    let mut rng = Rng::for_stream(14, "test/multigrid/finest");
    let diagonal: Vec<f64> = m
        .diagonal()
        .iter()
        .map(|d| d * (1.0 + 0.1 * rng.next_f64()))
        .collect();
    let close = m.with_diagonal(diagonal);
    let n = g.len();
    let r: Vec<f64> = (0..n).map(|_| rng.range_f64(-1.0, 1.0)).collect();
    let mut work = mg.cycle_space();
    let mut z = vec![0.0; n];
    mg.cycle(&mut work, &close, &r, &mut z);
    // By hand.
    let level = &mg.levels[0];
    let mut x = vec![0.0; n];
    forward_sweep(&close, &r, &mut x);
    let mut residual = vec![0.0; n];
    residual_of(&close, &r, &x, &mut residual);
    let mut coarse_rhs = vec![0.0; level.restriction.len()];
    level.restriction.apply(&residual, &mut coarse_rhs);
    let below = Multigrid {
        levels: mg.levels[1..].to_vec(),
        coarsest: mg.coarsest.clone(),
        solver: mg.solver.clone(),
    };
    let mut coarse = vec![0.0; coarse_rhs.len()];
    below.apply(&coarse_rhs, &mut coarse);
    level.prolongator.add_product(&coarse, &mut x);
    backward_sweep(&close, &r, &mut x);
    assert!(z.iter().zip(&x).all(|(a, b)| a.to_bits() == b.to_bits()));
    let mut own = vec![0.0; n];
    mg.apply(&r, &mut own);
    assert!(own.iter().zip(&z).any(|(a, b)| a.to_bits() != b.to_bits()));
}

#[test]
fn the_v_cycle_is_symmetric() {
    // For a symmetric matrix the V-cycle (a forward sweep, the coarse
    // correction, a backward sweep) is a symmetric operator B, as the
    // conjugate gradient method needs: x·By = y·Bx to 10⁻¹².
    let g = grid(10_000);
    let mg = Multigrid::new(substep_matrix(g, D_EARTH));
    assert!(mg.levels() >= 3, "{}", mg.levels());
    let mut rng = Rng::for_stream(4, "test/multigrid/symmetric");
    let n = g.len();
    for _ in 0..3 {
        let x: Vec<f64> = (0..n).map(|_| rng.range_f64(-1.0, 1.0)).collect();
        let y: Vec<f64> = (0..n).map(|_| rng.range_f64(-1.0, 1.0)).collect();
        let (mut bx, mut by) = (vec![0.0; n], vec![0.0; n]);
        mg.apply(&x, &mut bx);
        mg.apply(&y, &mut by);
        let xby: f64 = x.iter().zip(&by).map(|(a, b)| a * b).sum();
        let ybx: f64 = y.iter().zip(&bx).map(|(a, b)| a * b).sum();
        assert!(
            (xby - ybx).abs() <= 1e-12 * xby.abs().max(ybx.abs()),
            "{xby} against {ybx}"
        );
    }
}

#[test]
fn the_hierarchy_is_reported() {
    // The Earth-like sub-step on 10,000 cells: several levels, each at
    // most 80 % of the one before, the coarsest at most 400 cells, and an
    // operator complexity between 1 and 2.
    let g = grid(10_000);
    let mg = Multigrid::new(substep_matrix(g, D_EARTH));
    let sizes = mg.level_sizes();
    println!("levels {sizes:?}, complexity {}", mg.operator_complexity());
    assert_eq!(sizes.len(), mg.levels());
    assert_eq!(sizes[0], 10_000);
    assert_eq!(mg.len(), 10_000);
    assert!(!mg.is_empty());
    assert!(sizes.len() >= 3);
    for w in sizes.windows(2) {
        assert!((w[1] as f64) <= COARSENING_LIMIT * w[0] as f64, "{sizes:?}");
    }
    assert!(*sizes.last().unwrap() <= COARSEST_CELLS);
    let c = mg.operator_complexity();
    assert!(c > 1.0 && c < 2.0, "{c}");
}

#[test]
fn small_systems_are_solved_by_the_coarsest_level() {
    // 300 rows: no coarsening, the dense factorisation solves the system,
    // so one iteration reaches the tolerance; the solution matches the
    // dense one to 10⁻¹².
    let mut rng = Rng::for_stream(5, "test/multigrid/small");
    let s = random_spd(300, 0.001, 0.1, &mut rng);
    let mg = Multigrid::new(s.clone());
    assert_eq!(mg.levels(), 1);
    let b: Vec<f64> = (0..300).map(|_| rng.range_f64(-1.0, 1.0)).collect();
    let mut x = vec![0.0; 300];
    let solve = conjugate_gradient(&s, &mg, &b, 1e-12, 100, &mut x);
    assert_eq!(solve.iterations, 1, "{solve:?}");
    let (worst, largest) = error(&x, &dense_solve(dense_of(&s), b));
    assert!(worst <= 1e-12 * largest, "{worst:e}");
}

#[test]
fn a_level_that_cannot_coarsen_is_smoothed() {
    // A coupling of 10⁻⁹ on 2,000 cells: no strong connection, so no
    // aggregate; the only level is smoothed by Gauss–Seidel sweeps, which
    // solve the nearly diagonal system within two iterations.
    let g = grid(2_000);
    let m = substep_matrix(g, 1e-9);
    let mg = Multigrid::new(m.clone());
    assert_eq!(mg.levels(), 1);
    let b = substep_rhs(g);
    let mut x = vec![0.0; g.len()];
    let solve = conjugate_gradient(&m, &mg, &b, 1e-12, 100, &mut x);
    assert!(solve.converged() && solve.iterations <= 2, "{solve:?}");
}

#[test]
fn repeatable_bit_for_bit() {
    // Two hierarchies of the same matrix, two solves: the same bits.
    let g = grid(2_500);
    let b = substep_rhs(g);
    let run = || {
        let m = substep_matrix(g, D_EARTH);
        let mg = Multigrid::new(m.clone());
        let mut x = vec![0.0; g.len()];
        let solve = conjugate_gradient(&m, &mg, &b, 1e-10, 100, &mut x);
        (x, solve.iterations, solve.relative_residual)
    };
    let (x, i, r) = run();
    let (y, j, s) = run();
    assert_eq!(i, j);
    assert_eq!(r.to_bits(), s.to_bits());
    assert!(x.iter().zip(&y).all(|(a, b)| a.to_bits() == b.to_bits()));
}

#[test]
fn the_stops() {
    // A zero right-hand side gives zero at once; a start that solves a
    // diagonal system exactly is returned unchanged; an iteration limit of
    // 2 stops a hard system unconverged after 2 iterations; for both
    // solvers. (BiCGSTAB's breakdowns and restarts are those of
    // `sparse::bicgstab`, the same code, tested there.)
    let g = grid(2_500);
    let m = substep_matrix(g, 990.0 * D_EARTH);
    let mg = Multigrid::new(m.clone());
    let n = g.len();
    let mut x = vec![1.0; n];
    let solve = conjugate_gradient(&m, &mg, &vec![0.0; n], 1e-10, 100, &mut x);
    assert!(solve.converged() && solve.iterations == 0);
    assert!(x.iter().all(|&v| v == 0.0));
    let mut x = vec![0.0; n];
    let solve = conjugate_gradient(&m, &mg, &substep_rhs(g), 1e-14, 2, &mut x);
    assert!(!solve.converged() && solve.iterations == 2, "{solve:?}");
    let ma = Multigrid::with_prolongators_of(m.clone(), &mg);
    let mut x = vec![0.0; n];
    let solve = bicgstab(&m, &ma, &substep_rhs(g), 1e-14, 2, &mut x);
    assert!(!solve.converged() && solve.iterations == 2, "{solve:?}");
    let mut x = vec![1.0; n];
    let solve = bicgstab(&m, &ma, &vec![0.0; n], 1e-10, 100, &mut x);
    assert!(solve.converged() && solve.iterations == 0);
    assert!(x.iter().all(|&v| v == 0.0));
    let d = SparseMatrix::from_rows(vec![2.0, 4.0, 8.0], &[vec![], vec![], vec![]]);
    let md = Multigrid::new(d.clone());
    let mut x = vec![0.5, 0.25, 0.125];
    let solve = conjugate_gradient(&d, &md, &[1.0, 1.0, 1.0], 1e-12, 100, &mut x);
    assert_eq!(solve.iterations, 0);
    assert_eq!(x, vec![0.5, 0.25, 0.125]);
    let solve = bicgstab(&d, &md, &[1.0, 1.0, 1.0], 1e-12, 100, &mut x);
    assert_eq!(solve.iterations, 0);
    assert_eq!(x, vec![0.5, 0.25, 0.125]);
}

#[test]
#[should_panic(expected = "symmetric")]
fn new_needs_a_symmetric_matrix() {
    let m = SparseMatrix::from_rows(vec![2.0, 2.0], &[vec![(1, -1.0)], vec![(0, -0.5)]]);
    Multigrid::new(m);
}

#[test]
#[should_panic(expected = "positive diagonal")]
fn new_needs_a_positive_diagonal() {
    let m = SparseMatrix::from_rows(vec![2.0, 0.0], &[vec![(1, -1.0)], vec![(0, -1.0)]]);
    Multigrid::new(m);
}

#[test]
#[should_panic(expected = "the size of the symmetric hierarchy")]
fn with_prolongators_of_needs_the_same_size() {
    let mut rng = Rng::for_stream(6, "test/multigrid/sizes");
    let mg = Multigrid::new(random_spd(20, 0.1, 1.0, &mut rng));
    Multigrid::with_prolongators_of(random_spd(21, 0.1, 1.0, &mut rng), &mg);
}

#[test]
#[should_panic(expected = "one right-hand side per row")]
fn conjugate_gradient_needs_one_rhs_per_row() {
    let mut rng = Rng::for_stream(7, "test/multigrid/rhs");
    let m = random_spd(20, 0.1, 1.0, &mut rng);
    let mg = Multigrid::new(m.clone());
    conjugate_gradient(&m, &mg, &[1.0; 19], 1e-10, 100, &mut [0.0; 20]);
}

#[test]
#[should_panic(expected = "the multigrid of a matrix of its size")]
fn conjugate_gradient_needs_the_hierarchy_of_its_size() {
    let mut rng = Rng::for_stream(8, "test/multigrid/hierarchy");
    let m = random_spd(20, 0.1, 1.0, &mut rng);
    let mg = Multigrid::new(random_spd(21, 0.1, 1.0, &mut rng));
    conjugate_gradient(&m, &mg, &[1.0; 20], 1e-10, 100, &mut [0.0; 20]);
}

#[test]
#[should_panic(expected = "finite right-hand side")]
fn conjugate_gradient_needs_a_finite_rhs() {
    let mut rng = Rng::for_stream(9, "test/multigrid/finite");
    let m = random_spd(20, 0.1, 1.0, &mut rng);
    let mg = Multigrid::new(m.clone());
    let mut b = vec![1.0; 20];
    b[3] = f64::NAN;
    conjugate_gradient(&m, &mg, &b, 1e-10, 100, &mut [0.0; 20]);
}

#[test]
#[should_panic(expected = "tolerance in (0, 1)")]
fn conjugate_gradient_needs_a_tolerance_below_one() {
    let mut rng = Rng::for_stream(10, "test/multigrid/tolerance");
    let m = random_spd(20, 0.1, 1.0, &mut rng);
    let mg = Multigrid::new(m.clone());
    conjugate_gradient(&m, &mg, &[1.0; 20], 1.0, 100, &mut [0.0; 20]);
}

#[test]
#[should_panic(expected = "finite positive diagonal")]
fn conjugate_gradient_needs_a_positive_diagonal() {
    let mut rng = Rng::for_stream(15, "test/multigrid/diagonal");
    let m = random_spd(20, 0.1, 1.0, &mut rng);
    let mg = Multigrid::new(m.clone());
    let mut diagonal = m.diagonal().to_vec();
    diagonal[5] = 0.0;
    conjugate_gradient(
        &m.with_diagonal(diagonal),
        &mg,
        &[1.0; 20],
        1e-10,
        100,
        &mut [0.0; 20],
    );
}

#[test]
#[should_panic(expected = "one value per row")]
fn bicgstab_needs_one_value_per_row() {
    let mut rng = Rng::for_stream(11, "test/multigrid/values");
    let m = random_spd(20, 0.1, 1.0, &mut rng);
    let mg = Multigrid::new(m.clone());
    bicgstab(&m, &mg, &[1.0; 20], 1e-10, 100, &mut [0.0; 19]);
}

#[test]
#[should_panic(expected = "one value per row")]
fn apply_needs_one_value_per_row() {
    let mut rng = Rng::for_stream(12, "test/multigrid/apply");
    let mg = Multigrid::new(random_spd(20, 0.1, 1.0, &mut rng));
    mg.apply(&[1.0; 20], &mut [0.0; 21]);
}

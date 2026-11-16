//! Tests of the sparse matrix and BiCGSTAB (ticket M4-04): the product by
//! hand, solutions of random non-symmetric M-matrices against Gaussian
//! elimination done here, the stops on the tolerance and on the iteration
//! limit, the breakdowns (at the first step, and after it with a restart),
//! the exact cases, and the panics.
#![allow(clippy::disallowed_methods)]

use kuniumi_math::Rng;

use super::*;

/// The 3 × 3 matrix
///
/// ```text
/// [ 4 −1  0 ]
/// [−2  5 −1 ]
/// [ 0 −3  6 ]
/// ```
fn small() -> SparseMatrix {
    SparseMatrix::from_rows(
        vec![4.0, 5.0, 6.0],
        &[vec![(1, -1.0)], vec![(0, -2.0), (2, -1.0)], vec![(1, -3.0)]],
    )
}

#[test]
fn apply_by_hand() {
    let m = small();
    assert_eq!(m.len(), 3);
    assert!(!m.is_empty());
    assert_eq!(m.diagonal(), &[4.0, 5.0, 6.0]);
    let mut y = vec![0.0; 3];
    m.apply(&[1.0, 2.0, 3.0], &mut y);
    // 4 − 2 = 2; −2 + 10 − 3 = 5; −6 + 18 = 12.
    assert_eq!(y, vec![2.0, 5.0, 12.0]);
    let empty = SparseMatrix::from_rows(Vec::new(), &[]);
    assert!(empty.is_empty());
}

#[test]
fn with_diagonal_keeps_the_rest() {
    // Ticket M4-05: the diagonal 1, 2, 3 instead of 4, 5, 6 changes only
    // the diagonal terms of the product: 1 − 2 = −1; −2 + 4 − 3 = −1;
    // −6 + 9 = 3.
    let m = small().with_diagonal(vec![1.0, 2.0, 3.0]);
    assert_eq!(m.diagonal(), &[1.0, 2.0, 3.0]);
    assert_eq!(m.row(1), small().row(1));
    let mut y = vec![0.0; 3];
    m.apply(&[1.0, 2.0, 3.0], &mut y);
    assert_eq!(y, vec![-1.0, -1.0, 3.0]);
}

#[test]
#[should_panic(expected = "one diagonal entry per row")]
fn with_diagonal_needs_one_entry_per_row() {
    small().with_diagonal(vec![1.0, 2.0]);
}

#[test]
#[should_panic(expected = "finite entries on its diagonal")]
fn with_diagonal_needs_finite_entries() {
    small().with_diagonal(vec![1.0, f64::INFINITY, 3.0]);
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

/// A random non-symmetric M-matrix of size `n` like those of the ocean
/// (ticket M4-04): about four negative off-diagonal entries per row, of
/// random sizes over three decades and unrelated in the two directions,
/// and a diagonal equal to the row's sum of their magnitudes times
/// `1 + slack`, so weakly diagonally dominant where `slack` is 0, plus a
/// strictly dominant first row (a Dirichlet boundary) to make it
/// non-singular. Returned sparse and dense.
fn m_matrix(n: usize, slack: f64, rng: &mut Rng) -> (SparseMatrix, Vec<Vec<f64>>) {
    let mut dense = vec![vec![0.0; n]; n];
    let mut rows = vec![Vec::new(); n];
    let mut diagonal = vec![0.0; n];
    for i in 0..n {
        // A ring, so that the graph is connected, and three random columns.
        let mut columns = vec![(i + 1) % n, (i + n - 1) % n];
        for _ in 0..3 {
            columns.push(rng.below(n as u64) as usize);
        }
        columns.sort_unstable();
        columns.dedup();
        let mut sum = 0.0;
        for &j in &columns {
            if j == i {
                continue;
            }
            let v = -(10.0_f64).powf(rng.range_f64(-1.5, 1.5));
            rows[i].push((j as u32, v));
            dense[i][j] = v;
            sum -= v;
        }
        diagonal[i] = sum * (1.0 + slack) + if i == 0 { 1.0 } else { 0.0 };
        dense[i][i] = diagonal[i];
    }
    (SparseMatrix::from_rows(diagonal, &rows), dense)
}

#[test]
fn bicgstab_matches_a_dense_solution() {
    // Sizes 5 to 400, diagonally dominant by 0 to 10 % (weakly dominant
    // rows are the hard case: the matrix is then nearly singular apart
    // from the boundary row). The solution within 1e-8 of the largest
    // |x| of the dense one; the reported residual is the true one.
    let mut rng = Rng::for_stream(11, "test/sparse/m-matrix");
    for &(n, slack) in &[(5, 0.1), (20, 0.0), (60, 0.01), (150, 0.0), (400, 0.05)] {
        let (m, dense) = m_matrix(n, slack, &mut rng);
        let b: Vec<f64> = (0..n).map(|_| rng.range_f64(-1.0, 1.0)).collect();
        let exact = dense_solve(dense, b.clone());
        let mut x = vec![0.0; n];
        let solve = bicgstab(&m, &b, 1e-12, 10_000, &mut x);
        assert!(solve.converged(), "n = {n}: {solve:?}");
        assert!(solve.relative_residual <= 1e-12, "n = {n}: {solve:?}");
        let largest = exact.iter().fold(0.0_f64, |a, v| a.max(v.abs()));
        for (k, (a, e)) in x.iter().zip(&exact).enumerate() {
            assert!(
                (a - e).abs() <= 1e-8 * largest,
                "n = {n}, x[{k}] = {a} against {e} (largest {largest})"
            );
        }
        // The true residual.
        let mut ax = vec![0.0; n];
        m.apply(&x, &mut ax);
        let (mut rr, mut bb) = (0.0, 0.0);
        for k in 0..n {
            rr += (b[k] - ax[k]) * (b[k] - ax[k]);
            bb += b[k] * b[k];
        }
        assert!(
            (rr / bb).sqrt() <= 1e-11,
            "n = {n}: true residual {}",
            (rr / bb).sqrt()
        );
    }
}

#[test]
fn bicgstab_starts_from_the_values() {
    // A start at the solution of a diagonal system is returned bit for bit
    // after no iteration; a start near the solution needs fewer
    // iterations than a start at zero.
    let d = SparseMatrix::from_rows(vec![2.0, 4.0, 8.0], &[vec![], vec![], vec![]]);
    let mut x = vec![0.5, 0.25, 0.125];
    let solve = bicgstab(&d, &[1.0, 1.0, 1.0], 1e-12, 100, &mut x);
    assert_eq!(solve.iterations, 0);
    assert!(solve.converged());
    assert_eq!(solve.relative_residual.to_bits(), 0.0_f64.to_bits());
    assert_eq!(x, vec![0.5, 0.25, 0.125]);
    // The Jacobi preconditioner solves a diagonal system in one step.
    let mut x = vec![0.0; 3];
    let solve = bicgstab(&d, &[1.0, 1.0, 1.0], 1e-12, 100, &mut x);
    assert_eq!(solve.iterations, 1);
    assert_eq!(x, vec![0.5, 0.25, 0.125]);

    let mut rng = Rng::for_stream(12, "test/sparse/start");
    let (m, dense) = m_matrix(200, 0.0, &mut rng);
    let b: Vec<f64> = (0..200).map(|_| rng.range_f64(-1.0, 1.0)).collect();
    let exact = dense_solve(dense, b.clone());
    let mut cold = vec![0.0; 200];
    let from_zero = bicgstab(&m, &b, 1e-12, 10_000, &mut cold);
    let mut warm: Vec<f64> = exact.iter().map(|v| v * (1.0 + 1e-6)).collect();
    let from_near = bicgstab(&m, &b, 1e-12, 10_000, &mut warm);
    assert!(from_zero.converged() && from_near.converged());
    assert!(
        from_near.iterations < from_zero.iterations,
        "{from_near:?} against {from_zero:?}"
    );
}

#[test]
fn bicgstab_stops_at_the_iteration_limit() {
    let mut rng = Rng::for_stream(13, "test/sparse/limit");
    let (m, _) = m_matrix(300, 0.0, &mut rng);
    let b: Vec<f64> = (0..300).map(|_| rng.range_f64(-1.0, 1.0)).collect();
    let mut x = vec![0.0; 300];
    let solve = bicgstab(&m, &b, 1e-12, 3, &mut x);
    assert_eq!(solve.iterations, 3);
    assert!(!solve.converged());
    assert!(solve.relative_residual > 1e-12 && solve.relative_residual.is_finite());
    // With no iteration allowed the start stays.
    let mut x = vec![0.0; 300];
    let solve = bicgstab(&m, &b, 1e-12, 0, &mut x);
    assert_eq!(solve.iterations, 0);
    assert!(!solve.converged());
    assert_eq!(solve.relative_residual.to_bits(), 1.0_f64.to_bits());
    assert!(x.iter().all(|&v| v == 0.0));
}

#[test]
fn bicgstab_stops_on_a_breakdown_at_the_first_step() {
    // A = [[1, 3], [−1, 1]], b = (1, −1), the solution (1, 0). From zero,
    // r = b and v = A D⁻¹ r = (−2, −2), so r̃·v = 0: the method breaks down
    // before its first step and returns the start, unconverged. From
    // (0, 0.1) it solves the system.
    let m = SparseMatrix::from_rows(vec![1.0, 1.0], &[vec![(1, 3.0)], vec![(0, -1.0)]]);
    let mut x = vec![0.0; 2];
    let solve = bicgstab(&m, &[1.0, -1.0], 1e-12, 100, &mut x);
    assert_eq!(solve.iterations, 0);
    assert!(!solve.converged());
    assert_eq!(solve.relative_residual.to_bits(), 1.0_f64.to_bits());
    assert!(x.iter().all(|v| v.to_bits() == 0));
    let mut x = vec![0.0, 0.1];
    let solve = bicgstab(&m, &[1.0, -1.0], 1e-12, 100, &mut x);
    assert!(solve.converged(), "{solve:?}");
    assert!((x[0] - 1.0).abs() < 1e-11 && x[1].abs() < 1e-11, "{x:?}");
}

#[test]
fn bicgstab_restarts_after_a_breakdown() {
    // A = [[1, 1, −1], [1, 1, 0.5], [1, −0.5, 1]] (unit diagonal, so the
    // preconditioner is the identity), b = e₁; the solution, by hand,
    // (5/9, −2/9, −2/3). From zero: r₀ = r̃ = e₁, v = A e₁ = (1, 1, 1),
    // α = 1, s = (0, −1, −1), t = A s = (0, −1.5, −0.5), ω = t·s / t·t =
    // 2 / 2.5 = 0.8, r₁ = s − ω t = (0, 0.2, −0.6): r̃·r₁ = 0 exactly, a
    // breakdown after the first step. Restarted with the true residual as
    // the new shadow residual the method solves the system within n = 3
    // further steps (BiCG's finite termination); keeping the old shadow
    // residual, or stopping, would leave it at r₁ after one step.
    let m = SparseMatrix::from_rows(
        vec![1.0, 1.0, 1.0],
        &[
            vec![(1, 1.0), (2, -1.0)],
            vec![(0, 1.0), (2, 0.5)],
            vec![(0, 1.0), (1, -0.5)],
        ],
    );
    let mut x = vec![0.0; 3];
    let solve = bicgstab(&m, &[1.0, 0.0, 0.0], 1e-12, 100, &mut x);
    assert!(solve.converged(), "{solve:?}");
    assert!(solve.iterations > 1 && solve.iterations <= 4, "{solve:?}");
    for (a, e) in x.iter().zip([5.0 / 9.0, -2.0 / 9.0, -2.0 / 3.0]) {
        assert!((a - e).abs() < 1e-12, "{x:?}");
    }
}

#[test]
fn a_zero_right_hand_side_gives_zero() {
    let m = small();
    let mut x = vec![3.0, -2.0, 1.0];
    let solve = bicgstab(&m, &[0.0, 0.0, 0.0], 1e-12, 100, &mut x);
    assert!(solve.converged());
    assert_eq!(solve.iterations, 0);
    assert_eq!(solve.relative_residual.to_bits(), 0.0_f64.to_bits());
    assert!(x.iter().all(|&v| v == 0.0 && v.is_sign_positive()));
}

#[test]
fn bicgstab_is_bit_for_bit_repeatable() {
    let mut rng = Rng::for_stream(14, "test/sparse/repeat");
    let (m, _) = m_matrix(250, 0.02, &mut rng);
    let b: Vec<f64> = (0..250).map(|_| rng.range_f64(-1.0, 1.0)).collect();
    let mut x = vec![0.0; 250];
    let mut y = vec![0.0; 250];
    let s = bicgstab(&m, &b, 1e-10, 1_000, &mut x);
    let t = bicgstab(&m, &b, 1e-10, 1_000, &mut y);
    assert_eq!(s, t);
    assert!(x.iter().zip(&y).all(|(a, b)| a.to_bits() == b.to_bits()));
}

#[test]
#[should_panic(expected = "one diagonal entry per row")]
fn from_rows_needs_one_row_per_diagonal_entry() {
    SparseMatrix::from_rows(vec![1.0, 1.0], &[vec![]]);
}

#[test]
#[should_panic(expected = "strictly increasing columns")]
fn from_rows_needs_sorted_columns() {
    SparseMatrix::from_rows(
        vec![1.0, 1.0, 1.0],
        &[vec![(2, -1.0), (1, -1.0)], vec![], vec![]],
    );
}

#[test]
#[should_panic(expected = "off the diagonal")]
fn from_rows_keeps_the_diagonal_apart() {
    SparseMatrix::from_rows(vec![1.0, 1.0], &[vec![(0, -1.0)], vec![]]);
}

#[test]
#[should_panic(expected = "column 2 of row 0")]
fn from_rows_needs_columns_in_range() {
    SparseMatrix::from_rows(vec![1.0, 1.0], &[vec![(2, -1.0)], vec![]]);
}

#[test]
#[should_panic(expected = "finite entries")]
fn from_rows_needs_finite_entries() {
    SparseMatrix::from_rows(vec![1.0, f64::NAN], &[vec![], vec![]]);
}

#[test]
#[should_panic(expected = "one right-hand side per row")]
fn bicgstab_needs_one_right_hand_side_per_row() {
    bicgstab(&small(), &[1.0, 1.0], 1e-10, 10, &mut [0.0; 3]);
}

#[test]
#[should_panic(expected = "one value per row")]
fn bicgstab_needs_one_value_per_row() {
    bicgstab(&small(), &[1.0; 3], 1e-10, 10, &mut [0.0; 2]);
}

#[test]
#[should_panic(expected = "tolerance in (0, 1)")]
fn bicgstab_needs_a_tolerance_below_one() {
    bicgstab(&small(), &[1.0; 3], 1.0, 10, &mut [0.0; 3]);
}

#[test]
#[should_panic(expected = "non-zero diagonal")]
fn bicgstab_needs_a_non_zero_diagonal() {
    let m = SparseMatrix::from_rows(vec![1.0, 0.0], &[vec![(1, -1.0)], vec![(0, -1.0)]]);
    bicgstab(&m, &[1.0, 1.0], 1e-10, 10, &mut [0.0; 2]);
}

#[test]
#[should_panic(expected = "finite starting values")]
fn bicgstab_needs_finite_starting_values() {
    bicgstab(&small(), &[1.0; 3], 1e-10, 10, &mut [0.0, f64::NAN, 0.0]);
}

#[test]
#[should_panic(expected = "finite right-hand side")]
fn bicgstab_needs_a_finite_right_hand_side() {
    bicgstab(
        &small(),
        &[1.0, f64::INFINITY, 1.0],
        1e-10,
        10,
        &mut [0.0; 3],
    );
}

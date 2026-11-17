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

// ---------------------------------------------------------------------------
// The downwind order and the Gauss–Seidel preconditioner (ticket M4-14).

/// Whether `order` holds every index below `n` once.
fn is_permutation(order: &[u32], n: usize) -> bool {
    let mut seen = vec![false; n];
    order.len() == n
        && order.iter().all(|&i| {
            let fresh = (i as usize) < n && !seen[i as usize];
            if fresh {
                seen[i as usize] = true;
            }
            fresh
        })
}

/// The position of each index in `order`.
fn positions(order: &[u32]) -> Vec<usize> {
    let mut position = vec![0; order.len()];
    for (q, &i) in order.iter().enumerate() {
        position[i as usize] = q;
    }
    position
}

/// A random upwind matrix of size `n` without cycles: each cell receives
/// from up to three cells of lower rank in a random ranking (random
/// inflows over three decades), and its diagonal is its outflows plus a
/// mass of 0.01 to 1, as the moisture's upwind advection with its mass
/// term (column sums positive). Returned sparse and dense, with the
/// ranking.
fn acyclic_flow(n: usize, rng: &mut Rng) -> (SparseMatrix, Vec<Vec<f64>>) {
    // A random ranking: cell `rank_order[r]` has rank r.
    let mut keys: Vec<(u64, usize)> = (0..n).map(|i| (rng.next_u64(), i)).collect();
    keys.sort_unstable();
    let rank_order: Vec<usize> = keys.iter().map(|&(_, i)| i).collect();
    let mut dense = vec![vec![0.0; n]; n];
    let mut outflow = vec![0.0; n];
    for r in 1..n {
        let i = rank_order[r];
        for _ in 0..3 {
            let j = rank_order[rng.below(r as u64) as usize];
            let w = (10.0_f64).powf(rng.range_f64(-1.5, 1.5));
            dense[i][j] -= w;
            outflow[j] += w;
        }
    }
    let mut rows = vec![Vec::new(); n];
    let mut diagonal = vec![0.0; n];
    for i in 0..n {
        for (j, &v) in dense[i].iter().enumerate() {
            if v != 0.0 {
                rows[i].push((j as u32, v));
            }
        }
        diagonal[i] = outflow[i] + rng.range_f64(0.01, 1.0);
        dense[i][i] = diagonal[i];
    }
    (SparseMatrix::from_rows(diagonal, &rows), dense)
}

#[test]
fn the_downwind_order_is_a_permutation() {
    // Random M-matrices with cycles everywhere (rings and random columns)
    // and an acyclic flow: every index once, and the same order bit for bit
    // when computed again.
    let mut rng = Rng::for_stream(15, "test/sparse/order");
    for &(n, slack) in &[(1, 0.1), (5, 0.1), (60, 0.0), (400, 0.05)] {
        let (m, _) = m_matrix(n, slack, &mut rng);
        let order = downwind_order(&m);
        assert!(is_permutation(&order, n), "n = {n}: {order:?}");
        assert_eq!(order, downwind_order(&m));
    }
    let (m, _) = acyclic_flow(300, &mut rng);
    assert!(is_permutation(&downwind_order(&m), 300));
    let empty = SparseMatrix::from_rows(Vec::new(), &[]);
    assert!(downwind_order(&empty).is_empty());
}

#[test]
fn an_acyclic_flow_is_ordered_downwind() {
    // Without cycles every cell comes after the cells that flow into it,
    // so the matrix renumbered in the order is lower triangular, one
    // Gauss–Seidel sweep solves it exactly and BiCGSTAB preconditioned by
    // it converges in one iteration, to the dense solution (within 10⁻¹²
    // of the largest value). In index order the same system needs more
    // iterations.
    let mut rng = Rng::for_stream(16, "test/sparse/acyclic");
    for n in [10, 100, 400] {
        let (m, dense) = acyclic_flow(n, &mut rng);
        let order = downwind_order(&m);
        let position = positions(&order);
        for i in 0..n {
            let (columns, values) = m.row(i);
            for (&j, &v) in columns.iter().zip(values) {
                assert!(
                    v >= 0.0 || position[j as usize] < position[i],
                    "n = {n}: {j} → {i}"
                );
            }
        }
        let p = m.permuted(&order);
        for q in 0..n {
            assert!(
                p.row(q).0.iter().all(|&c| (c as usize) < q),
                "row {q} has an upper entry"
            );
        }
        let b: Vec<f64> = (0..n).map(|_| rng.range_f64(0.0, 1.0)).collect();
        let exact = dense_solve(dense, b.clone());
        let pb: Vec<f64> = order.iter().map(|&i| b[i as usize]).collect();
        let mut px = vec![0.0; n];
        let solve = bicgstab_gauss_seidel(&p, &pb, 1e-12, 100, &mut px);
        assert!(solve.converged(), "n = {n}: {solve:?}");
        assert_eq!(solve.iterations, 1, "n = {n}: {solve:?}");
        let largest = exact.iter().fold(0.0_f64, |a, v| a.max(v.abs()));
        for (q, &i) in order.iter().enumerate() {
            let e = exact[i as usize];
            assert!(
                (px[q] - e).abs() <= 1e-12 * largest,
                "n = {n}: {} against {e}",
                px[q]
            );
        }
        let mut x = vec![0.0; n];
        let unordered = bicgstab_gauss_seidel(&m, &b, 1e-12, 100, &mut x);
        println!(
            "n = {n}: in index order {} iterations",
            unordered.iterations
        );
        assert!(unordered.converged() && unordered.iterations > 1);
    }
}

#[test]
fn a_cycle_is_broken_at_its_weakest_inflow() {
    // A ring 0 → 1 → … → 5 → 0 (cell i receives 1 from cell i − 1, except
    // cell 4, which receives only 0.01 from cell 3), a cell 6 that
    // receives from cell 2, and a cell 7 that receives from nothing; every
    // diagonal 2. By hand: cell 7 is ready first; then no cell is ready, and
    // cell 4 has the least inflow from unordered cells over its diagonal
    // (0.005 against 0.5), so the ring is cut before it and followed
    // downwind: 5, 0, 1, 2; then cells 3 and 6 are ready, by index.
    let mut rows = vec![Vec::new(); 8];
    for i in 0..6_u32 {
        let from = (i + 5) % 6;
        rows[i as usize].push((from, if i == 4 { -0.01 } else { -1.0 }));
    }
    rows[6].push((2, -1.0));
    let m = SparseMatrix::from_rows(vec![2.0; 8], &rows);
    assert_eq!(downwind_order(&m), vec![7, 4, 5, 0, 1, 2, 3, 6]);
    // Entries that are not inflows (positive, or zero) do not order.
    let mut rows = rows.clone();
    rows[7].push((6, 0.0));
    rows[3].push((7, 0.5));
    let m = SparseMatrix::from_rows(vec![2.0; 8], &rows);
    assert_eq!(downwind_order(&m), vec![7, 4, 5, 0, 1, 2, 3, 6]);
}

#[test]
fn a_cut_weighs_only_the_inflow_left_behind() {
    // Two rings, 0 ⇄ 1 (inflows 1) and 2 ⇄ 3, where cell 2 also receives
    // 1.5 from cell 0 and only 0.1 from cell 3, and cell 3 receives 1 from
    // cell 2; every diagonal 2. By hand: no cell is ready; the shares of
    // the inflow from unordered cells are 0.5, 0.5, 0.8 and 0.5, so ring 0
    // ⇄ 1 is cut first, before cell 0 (the lowest index), and cell 1 follows.
    // Then cell 2 has only its 0.1 from cell 3 left (share 0.05, against
    // cell 3's 0.5), so the second ring is cut before cell 2, although cell
    // 2's whole inflow (share 0.8) is the largest: 0, 1, 2, 3. Weighing the
    // whole inflow would cut before cell 3 instead: 0, 1, 3, 2.
    let m = SparseMatrix::from_rows(
        vec![2.0; 4],
        &[
            vec![(1, -1.0)],
            vec![(0, -1.0)],
            vec![(0, -1.5), (3, -0.1)],
            vec![(2, -1.0)],
        ],
    );
    assert_eq!(downwind_order(&m), vec![0, 1, 2, 3]);
}

#[test]
fn the_permuted_matrix_is_the_same_operator() {
    // Row and column q of the renumbered matrix are row and column
    // order[q] of the matrix: the same entries (bit for bit) at the new
    // places, each row's columns increasing, and the same product up to
    // the order of the sums (within 10⁻¹⁴ of the row's magnitude).
    let mut rng = Rng::for_stream(17, "test/sparse/permuted");
    let n = 200;
    let (m, _) = m_matrix(n, 0.02, &mut rng);
    let mut keys: Vec<(u64, u32)> = (0..n as u32).map(|i| (rng.next_u64(), i)).collect();
    keys.sort_unstable();
    let order: Vec<u32> = keys.iter().map(|&(_, i)| i).collect();
    let position = positions(&order);
    let p = m.permuted(&order);
    assert_eq!(p.len(), n);
    let x: Vec<f64> = (0..n).map(|_| rng.range_f64(-1.0, 1.0)).collect();
    let px: Vec<f64> = order.iter().map(|&i| x[i as usize]).collect();
    let (mut ax, mut pax) = (vec![0.0; n], vec![0.0; n]);
    m.apply(&x, &mut ax);
    p.apply(&px, &mut pax);
    for (q, &i) in order.iter().enumerate() {
        let i = i as usize;
        assert_same_bits(p.diagonal()[q], m.diagonal()[i]);
        let (columns, values) = m.row(i);
        let mut expected: Vec<(u32, u64)> = columns
            .iter()
            .zip(values)
            .map(|(&j, v)| (position[j as usize] as u32, v.to_bits()))
            .collect();
        expected.sort_unstable();
        let (pc, pv) = p.row(q);
        let actual: Vec<(u32, u64)> = pc.iter().zip(pv).map(|(&c, v)| (c, v.to_bits())).collect();
        assert_eq!(actual, expected, "row {q}");
        let scale = m.diagonal()[i] * 2.0;
        assert!((pax[q] - ax[i]).abs() <= 1e-14 * scale, "row {q}");
    }
}

/// Asserts that two doubles are the same bits.
fn assert_same_bits(a: f64, b: f64) {
    assert_eq!(a.to_bits(), b.to_bits(), "{a} against {b}");
}

#[test]
fn bicgstab_gauss_seidel_matches_a_dense_solution() {
    // The random M-matrices of `bicgstab_matches_a_dense_solution`, in index
    // order and in their downwind order: the solution within 10⁻¹⁰ of the
    // largest |x| of the dense one (the weakly dominant ones, nearly
    // singular, within 10⁻⁸ as for `bicgstab`), the reported residual the
    // true one; fewer iterations than the Jacobi preconditioner's.
    let mut rng = Rng::for_stream(11, "test/sparse/m-matrix");
    for &(n, slack) in &[(5, 0.1), (20, 0.0), (60, 0.01), (150, 0.0), (400, 0.05)] {
        let (m, dense) = m_matrix(n, slack, &mut rng);
        let b: Vec<f64> = (0..n).map(|_| rng.range_f64(-1.0, 1.0)).collect();
        let exact = dense_solve(dense, b.clone());
        let largest = exact.iter().fold(0.0_f64, |a, v| a.max(v.abs()));
        let bound = if slack == 0.0 { 1e-8 } else { 1e-10 };
        let mut jacobi = vec![0.0; n];
        let reference = bicgstab(&m, &b, 1e-12, 10_000, &mut jacobi);
        for order in [(0..n as u32).collect::<Vec<u32>>(), downwind_order(&m)] {
            let p = m.permuted(&order);
            let pb: Vec<f64> = order.iter().map(|&i| b[i as usize]).collect();
            let mut px = vec![0.0; n];
            let solve = bicgstab_gauss_seidel(&p, &pb, 1e-12, 10_000, &mut px);
            assert!(solve.converged(), "n = {n}: {solve:?}");
            assert!(solve.relative_residual <= 1e-12, "n = {n}: {solve:?}");
            assert!(
                solve.iterations <= reference.iterations,
                "n = {n}: {solve:?} against Jacobi {reference:?}"
            );
            for (q, &i) in order.iter().enumerate() {
                let e = exact[i as usize];
                assert!(
                    (px[q] - e).abs() <= bound * largest,
                    "n = {n}, x[{i}] = {} against {e} (largest {largest})",
                    px[q]
                );
            }
            let mut ax = vec![0.0; n];
            p.apply(&px, &mut ax);
            let (mut rr, mut bb) = (0.0, 0.0);
            for q in 0..n {
                rr += (pb[q] - ax[q]) * (pb[q] - ax[q]);
                bb += pb[q] * pb[q];
            }
            assert!(
                (rr / bb).sqrt() <= 1e-11,
                "n = {n}: true residual {}",
                (rr / bb).sqrt()
            );
        }
    }
}

#[test]
fn bicgstab_gauss_seidel_stops() {
    // A zero right-hand side gives zero; a start that solves the system is
    // returned bit for bit after no iteration; the iteration limit; the
    // same bits when solved again.
    let m = small();
    let mut x = vec![3.0, -2.0, 1.0];
    let solve = bicgstab_gauss_seidel(&m, &[0.0; 3], 1e-12, 100, &mut x);
    assert!(solve.converged() && solve.iterations == 0);
    assert!(x.iter().all(|&v| v == 0.0 && v.is_sign_positive()));
    let d = SparseMatrix::from_rows(vec![2.0, 4.0, 8.0], &[vec![], vec![], vec![]]);
    let mut x = vec![0.5, 0.25, 0.125];
    let solve = bicgstab_gauss_seidel(&d, &[1.0; 3], 1e-12, 100, &mut x);
    assert!(solve.converged() && solve.iterations == 0);
    assert_eq!(x, vec![0.5, 0.25, 0.125]);
    let mut rng = Rng::for_stream(13, "test/sparse/limit");
    let (m, _) = m_matrix(300, 0.0, &mut rng);
    let b: Vec<f64> = (0..300).map(|_| rng.range_f64(-1.0, 1.0)).collect();
    let mut x = vec![0.0; 300];
    let solve = bicgstab_gauss_seidel(&m, &b, 1e-12, 2, &mut x);
    assert_eq!(solve.iterations, 2);
    assert!(!solve.converged());
    assert!(solve.relative_residual > 1e-12 && solve.relative_residual.is_finite());
    let mut y = vec![0.0; 300];
    let again = bicgstab_gauss_seidel(&m, &b, 1e-12, 2, &mut y);
    assert_eq!(solve, again);
    assert!(x.iter().zip(&y).all(|(a, b)| a.to_bits() == b.to_bits()));
}

#[test]
#[should_panic(expected = "every row once")]
fn permuted_needs_a_permutation() {
    small().permuted(&[0, 0, 2]);
}

#[test]
#[should_panic(expected = "every row once")]
fn permuted_needs_one_entry_per_row() {
    small().permuted(&[0, 1]);
}

#[test]
#[should_panic(expected = "positive finite diagonal")]
fn the_downwind_order_needs_a_positive_diagonal() {
    let m = SparseMatrix::from_rows(vec![1.0, 0.0], &[vec![(1, -1.0)], vec![(0, -1.0)]]);
    downwind_order(&m);
}

#[test]
#[should_panic(expected = "non-zero diagonal")]
fn bicgstab_gauss_seidel_needs_a_non_zero_diagonal() {
    let m = SparseMatrix::from_rows(vec![1.0, 0.0], &[vec![(1, -1.0)], vec![(0, -1.0)]]);
    bicgstab_gauss_seidel(&m, &[1.0, 1.0], 1e-10, 10, &mut [0.0; 2]);
}

#[test]
#[should_panic(expected = "finite right-hand side")]
fn bicgstab_gauss_seidel_needs_a_finite_right_hand_side() {
    bicgstab_gauss_seidel(&small(), &[1.0, f64::NAN, 1.0], 1e-10, 10, &mut [0.0; 3]);
}

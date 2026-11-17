//! Sparse linear systems that are not symmetric (ticket M4-04): a square
//! matrix stored by rows ([`SparseMatrix`]) and the BiCGSTAB method with
//! the Jacobi preconditioner ([`bicgstab`]); for upwind advection, a
//! downwind order of the rows ([`downwind_order`]), the matrix renumbered
//! in it ([`SparseMatrix::permuted`]) and BiCGSTAB preconditioned by a
//! Gauss–Seidel sweep in that order ([`bicgstab_gauss_seidel`], ticket
//! M4-14).
//!
//! The ocean's Stommel model ([`crate::stages::climate::ocean`]) solves
//! with [`bicgstab`] the advection–diffusion system of its streamfunction,
//! whose matrix is a non-symmetric M-matrix; the moisture
//! ([`crate::stages::climate::moisture`]) solves its upwind advection with
//! [`bicgstab_gauss_seidel`] in the downwind order of each month; the
//! symmetric systems of the diffusion and the sea surface's system are
//! solved with the multigrid of [`crate::multigrid`] (ticket M4-09), whose
//! BiCGSTAB is this one with the V-cycle for the preconditioner.
//!
//! # The downwind order
//!
//! An upwind discretisation of advection has off its diagonal only the
//! inflows of each cell (`a_ij < 0`: cell `i` receives from cell `j`).
//! Numbered so that every cell comes after the cells that flow into it,
//! the matrix is lower triangular where the flow has no cycles, and one
//! Gauss–Seidel sweep in that order solves it: the "downwind numbering"
//! of W. Hackbusch & T. Probst (1997, "Downwind Gauß–Seidel smoothing for
//! convection dominated problems", Numer. Linear Algebra Appl. 4, 85–102)
//! and J. Bey & G. Wittum (1997, "Downwind numbering: robust multigrid for
//! convection-diffusion problems", Appl. Numer. Math. 23, 177–192). Where
//! streamlines close the flow has cycles, and the order must cut each of
//! them somewhere; [`downwind_order`] cuts where the inflow left behind is
//! the smallest share of a row's diagonal, so that the entries above the
//! diagonal in the order, which the sweep does not solve, are as small as
//! its greedy choice makes them (on the moisture's systems of the
//! Earth-parameter world, 6 % of a row's diagonal on average on 10,000 to
//! 650,000 cells, against 36 % in index order on 10,000; ticket M4-14,
//! prototype).
//!
//! Everything runs sequentially in row order with IEEE `+ − × ÷` and
//! `sqrt`: every product and dot product sums in a fixed order
//! (docs/conventions.md, "Determinism", rules 2 and 6), and the downwind
//! order compares integers and the bits of non-negative doubles, ties
//! broken by index (rule 5), so the results are the same on every
//! platform.

use core::cmp::Reverse;
use std::collections::BinaryHeap;

use crate::operators::Solve;

/// A square sparse matrix: the diagonal, and the entries off the diagonal
/// by row (compressed sparse rows), each row's columns strictly
/// increasing.
#[derive(Clone, Debug, PartialEq)]
pub struct SparseMatrix {
    /// The diagonal entries.
    diagonal: Vec<f64>,
    /// The entries of row `i` are `offsets[i]..offsets[i + 1]`.
    offsets: Vec<usize>,
    /// The column of each entry.
    columns: Vec<u32>,
    /// The value of each entry.
    values: Vec<f64>,
}

impl SparseMatrix {
    /// The matrix with the diagonal `diagonal` and, for each row `i`, the
    /// entries `rows[i]` off the diagonal as `(column, value)` pairs.
    ///
    /// # Panics
    ///
    /// If `rows` does not hold one row per diagonal entry, if the rows have
    /// more than `u32::MAX` entries in all, if a row's columns are not
    /// strictly increasing, are not below the size or include the row
    /// itself, or if an entry is not finite.
    pub fn from_rows(diagonal: Vec<f64>, rows: &[Vec<(u32, f64)>]) -> Self {
        let n = diagonal.len();
        assert_eq!(
            rows.len(),
            n,
            "a sparse matrix needs one diagonal entry per row: {n} for {} rows",
            rows.len()
        );
        assert!(
            n <= u32::MAX as usize,
            "a sparse matrix has at most u32::MAX rows, got {n}"
        );
        assert!(
            diagonal.iter().all(|d| d.is_finite()),
            "a sparse matrix needs finite entries on its diagonal"
        );
        let mut offsets = Vec::with_capacity(n + 1);
        let mut columns = Vec::new();
        let mut values = Vec::new();
        offsets.push(0);
        for (i, row) in rows.iter().enumerate() {
            let mut previous: Option<u32> = None;
            for &(column, value) in row {
                assert!(
                    (column as usize) < n,
                    "a sparse matrix needs columns below its size: column {column} of row {i}, \
                     size {n}"
                );
                assert!(
                    column as usize != i,
                    "a sparse matrix keeps the diagonal apart: row {i} lists column {column} off \
                     the diagonal"
                );
                assert!(
                    previous.is_none_or(|p| column > p),
                    "a sparse matrix needs strictly increasing columns in each row: row {i}"
                );
                assert!(
                    value.is_finite(),
                    "a sparse matrix needs finite entries: row {i}, column {column} has {value}"
                );
                previous = Some(column);
                columns.push(column);
                values.push(value);
            }
            offsets.push(columns.len());
        }
        Self {
            diagonal,
            offsets,
            columns,
            values,
        }
    }

    /// The number of rows (and of columns).
    pub fn len(&self) -> usize {
        self.diagonal.len()
    }

    /// Whether the matrix has no rows.
    pub fn is_empty(&self) -> bool {
        self.diagonal.is_empty()
    }

    /// The diagonal entries.
    pub fn diagonal(&self) -> &[f64] {
        &self.diagonal
    }

    /// The same entries off the diagonal with the diagonal `diagonal`: the
    /// systems of the sea-surface temperature, which differ month by month
    /// only on the diagonal (ticket M4-05).
    ///
    /// # Panics
    ///
    /// If `diagonal` does not hold one finite entry per row.
    pub fn with_diagonal(&self, diagonal: Vec<f64>) -> Self {
        assert_eq!(
            diagonal.len(),
            self.len(),
            "with_diagonal needs one diagonal entry per row: {} for {} rows",
            diagonal.len(),
            self.len()
        );
        assert!(
            diagonal.iter().all(|d| d.is_finite()),
            "a sparse matrix needs finite entries on its diagonal"
        );
        Self {
            diagonal,
            offsets: self.offsets.clone(),
            columns: self.columns.clone(),
            values: self.values.clone(),
        }
    }

    /// The entries of row `i` off the diagonal: their columns and values.
    pub fn row(&self, i: usize) -> (&[u32], &[f64]) {
        let range = self.offsets[i]..self.offsets[i + 1];
        (&self.columns[range.clone()], &self.values[range])
    }

    /// `y = A x`: per row, the diagonal term, then the others in column
    /// order.
    ///
    /// # Panics
    ///
    /// If `x` or `y` does not hold one value per row.
    pub fn apply(&self, x: &[f64], y: &mut [f64]) {
        let n = self.len();
        assert!(
            x.len() == n && y.len() == n,
            "apply needs one value per row: {} and {} values for {n} rows",
            x.len(),
            y.len()
        );
        for (i, out) in y.iter_mut().enumerate() {
            let mut sum = self.diagonal[i] * x[i];
            for k in self.offsets[i]..self.offsets[i + 1] {
                sum += self.values[k] * x[self.columns[k] as usize];
            }
            *out = sum;
        }
    }

    /// The matrix renumbered by `order`: row and column `q` of the result
    /// are row and column `order[q]` of this one, `P A Pᵀ` with
    /// `(P x)_q = x_{order[q]}`. Each row keeps its entries, bit for bit,
    /// with their columns renumbered and sorted (ticket M4-14).
    ///
    /// # Panics
    ///
    /// If `order` does not list every row once.
    pub fn permuted(&self, order: &[u32]) -> Self {
        let n = self.len();
        let mut position = vec![u32::MAX; n];
        let mut valid = order.len() == n;
        for (q, &i) in order.iter().enumerate() {
            let i = i as usize;
            if !valid || i >= n || position[i] != u32::MAX {
                valid = false;
                break;
            }
            // At most u32::MAX rows (`from_rows`).
            position[i] = q as u32;
        }
        assert!(
            valid,
            "permuted needs an order that lists every row once: {} entries for {n} rows",
            order.len()
        );
        let mut offsets = Vec::with_capacity(n + 1);
        let mut columns = Vec::with_capacity(self.columns.len());
        let mut values = Vec::with_capacity(self.values.len());
        let mut row: Vec<(u32, f64)> = Vec::new();
        offsets.push(0);
        for &i in order {
            let (cols, vals) = self.row(i as usize);
            row.clear();
            row.extend(
                cols.iter()
                    .zip(vals)
                    .map(|(&j, &v)| (position[j as usize], v)),
            );
            // The columns of a row are distinct, so the order is unique.
            row.sort_unstable_by_key(|e| e.0);
            for &(c, v) in &row {
                columns.push(c);
                values.push(v);
            }
            offsets.push(columns.len());
        }
        Self {
            diagonal: order.iter().map(|&i| self.diagonal[i as usize]).collect(),
            offsets,
            columns,
            values,
        }
    }
}

/// The downwind order of the rows of `matrix` (module documentation, "The
/// downwind order"): `order[q]` is the row in place `q`, for
/// [`SparseMatrix::permuted`] (ticket M4-14).
///
/// The entries off the diagonal that are negative are the inflows of
/// their row (`a_ij < 0`: row `i` receives from row `j`); the others do
/// not order. Kahn's topological sort of the inflow graph (A. B. Kahn 1962,
/// "Topological sorting of large networks", Commun. ACM 5, 558–562): a row
/// is ready when every row it receives from is ordered, and the ready row
/// of the lowest index comes next. Where no row is ready (the flow has
/// cycles), the next is the unordered row whose inflow from unordered rows
/// `Σ_j (−a_ij)` is the smallest share of its diagonal `a_ii`, ties by
/// index. A binary heap of `(ready or not, the share's bits, index)` with
/// lazy deletion: `O(nnz log n)`.
///
/// # Panics
///
/// If a diagonal entry is not positive and finite.
pub fn downwind_order(matrix: &SparseMatrix) -> Vec<u32> {
    let n = matrix.len();
    let diagonal = &matrix.diagonal;
    // Written so that NaN fails.
    assert!(
        diagonal.iter().all(|&d| d > 0.0 && d.is_finite()),
        "the downwind order needs a positive finite diagonal"
    );
    // The inflows of each row (their number and sum), and the rows that
    // each row flows into (by increasing index).
    let mut count = vec![0u32; n];
    let mut remaining = vec![0.0; n];
    let mut starts = vec![0usize; n + 1];
    for (i, (c, r)) in count.iter_mut().zip(&mut remaining).enumerate() {
        let (cols, vals) = matrix.row(i);
        for (&j, &v) in cols.iter().zip(vals) {
            if v < 0.0 {
                *c += 1;
                *r += -v;
                starts[j as usize + 1] += 1;
            }
        }
    }
    for j in 0..n {
        starts[j + 1] += starts[j];
    }
    let mut fill = starts.clone();
    let mut receivers = vec![(0u32, 0.0); starts[n]];
    for i in 0..n {
        let (cols, vals) = matrix.row(i);
        for (&j, &v) in cols.iter().zip(vals) {
            if v < 0.0 {
                receivers[fill[j as usize]] = (i as u32, -v);
                fill[j as usize] += 1;
            }
        }
    }
    // The key of a row: ready rows first, then the smallest share; the
    // share is non-negative (a remaining sum that rounding takes below
    // zero, or to −0, counts as +0), so its bits order as its value.
    let key = |count: u32, remaining: f64, d: f64| -> (u8, u64) {
        if count == 0 {
            (0, 0)
        } else {
            let r = if remaining > 0.0 { remaining } else { 0.0 };
            (1, (r / d).to_bits())
        }
    };
    let mut keys: Vec<(u8, u64)> = (0..n)
        .map(|i| key(count[i], remaining[i], diagonal[i]))
        .collect();
    let mut heap: BinaryHeap<Reverse<(u8, u64, u32)>> = keys
        .iter()
        .enumerate()
        .map(|(i, k)| Reverse((k.0, k.1, i as u32)))
        .collect();
    let mut ordered = vec![false; n];
    let mut order = Vec::with_capacity(n);
    while let Some(Reverse((ready, share, i))) = heap.pop() {
        let i = i as usize;
        if ordered[i] || keys[i] != (ready, share) {
            // A stale entry: the row is ordered or its key has changed.
            continue;
        }
        ordered[i] = true;
        order.push(i as u32);
        for &(j, w) in &receivers[starts[i]..starts[i + 1]] {
            let j = j as usize;
            if !ordered[j] {
                count[j] -= 1;
                remaining[j] -= w;
                keys[j] = key(count[j], remaining[j], diagonal[j]);
                heap.push(Reverse((keys[j].0, keys[j].1, j as u32)));
            }
        }
    }
    order
}

/// `Σ_i x_i y_i`, sequentially in index order.
fn dot(x: &[f64], y: &[f64]) -> f64 {
    let mut sum = 0.0;
    for (a, b) in x.iter().zip(y) {
        sum += a * b;
    }
    sum
}

/// Solves `A x = b` by BiCGSTAB (H. A. van der Vorst 1992, "Bi-CGSTAB: a
/// fast and smoothly converging variant of Bi-CG for the solution of
/// nonsymmetric linear systems", SIAM J. Sci. Stat. Comput. 13, 631–644),
/// preconditioned by the diagonal of `A`, in the form of R. Barrett et al.
/// (1994, *Templates for the Solution of Linear Systems*, SIAM, §2.3.8,
/// Fig. 2.10), with the shadow residual the starting residual; `values`
/// hold the start and become the solution, `A = matrix`, `b = rhs`.
///
/// It stops when the residual it updates falls to
/// `‖b − A x‖₂ ≤ tolerance · ‖b‖₂` (checked at the start, at the half step
/// and at the end of each iteration) or after `max_iterations` iterations.
/// Where the method breaks down (`r̃·r`, `r̃·v`, `t·t` or `ω` zero or not
/// finite), it restarts from the current iterate, with its true residual
/// `b − A x` as the residual and the new shadow residual; a breakdown at
/// the first step of a start stops it (on a nearly singular, strongly
/// non-normal system rounding can make `r̃·r` vanish near the attainable
/// accuracy: a random coast of 2,000 cells on a planet of 1,500 km broke
/// down at a relative residual of 1.05 × 10⁻⁸ after 374 iterations and,
/// restarted, reached 7 × 10⁻¹¹ after 487; ticket M4-04). The returned
/// [`Solve`] says whether it converged, with the iterations (each of two
/// products with `A`) and the relative residual as the method updated it.
/// A start that solves the system exactly is returned unchanged, bit for
/// bit; a zero right-hand side gives zero, the solution of a non-singular
/// system. Convergence needs a matrix the method suits, such as a
/// non-singular M-matrix.
///
/// # Panics
///
/// If `rhs` or `values` does not hold one entry per row, if a diagonal
/// entry is zero, if a right-hand side or a value is not finite, if
/// `tolerance` is not in `(0, 1)`.
pub fn bicgstab(
    matrix: &SparseMatrix,
    rhs: &[f64],
    tolerance: f64,
    max_iterations: u32,
    values: &mut [f64],
) -> Solve {
    let n = matrix.len();
    assert_eq!(
        rhs.len(),
        n,
        "bicgstab needs one right-hand side per row: {} for {n} rows",
        rhs.len()
    );
    assert_eq!(
        values.len(),
        n,
        "bicgstab needs one value per row: {} values for {n} rows",
        values.len()
    );
    assert!(
        matrix.diagonal.iter().all(|&d| d != 0.0),
        "bicgstab needs a non-zero diagonal (the Jacobi preconditioner)"
    );
    check_system("bicgstab", rhs, values, tolerance);
    let d = &matrix.diagonal;
    let mut jacobi = |r: &[f64], z: &mut [f64]| {
        for k in 0..n {
            z[k] = r[k] / d[k];
        }
    };
    bicgstab_with(matrix, rhs, tolerance, max_iterations, values, &mut jacobi)
}

/// Solves `A x = b` by the BiCGSTAB method of [`bicgstab`] (the same
/// iteration, stops, breakdowns and restarts), right-preconditioned by one
/// forward Gauss–Seidel sweep in row order, `M = D + L` with `D` the
/// diagonal of `A` and `L` its entries below the diagonal: for a matrix
/// renumbered in its [`downwind_order`] ([`SparseMatrix::permuted`]), so
/// that the sweep follows the flow and solves an acyclic upwind system
/// exactly, in one iteration (ticket M4-14). `values` hold the start and
/// become the solution, `A = matrix`, `b = rhs`; the returned [`Solve`] as
/// [`bicgstab`]'s.
///
/// # Panics
///
/// As [`bicgstab`]: if `rhs` or `values` does not hold one entry per row,
/// if a diagonal entry is zero, if a right-hand side or a value is not
/// finite, if `tolerance` is not in `(0, 1)`.
pub fn bicgstab_gauss_seidel(
    matrix: &SparseMatrix,
    rhs: &[f64],
    tolerance: f64,
    max_iterations: u32,
    values: &mut [f64],
) -> Solve {
    let n = matrix.len();
    assert_eq!(
        rhs.len(),
        n,
        "bicgstab_gauss_seidel needs one right-hand side per row: {} for {n} rows",
        rhs.len()
    );
    assert_eq!(
        values.len(),
        n,
        "bicgstab_gauss_seidel needs one value per row: {} values for {n} rows",
        values.len()
    );
    assert!(
        matrix.diagonal.iter().all(|&d| d != 0.0),
        "bicgstab_gauss_seidel needs a non-zero diagonal (the Gauss–Seidel sweep divides by it)"
    );
    check_system("bicgstab_gauss_seidel", rhs, values, tolerance);
    let mut sweep = |r: &[f64], z: &mut [f64]| {
        for q in 0..n {
            let mut sum = r[q];
            for k in matrix.offsets[q]..matrix.offsets[q + 1] {
                let c = matrix.columns[k] as usize;
                if c > q {
                    // The columns increase: the rest of the row is above.
                    break;
                }
                sum -= matrix.values[k] * z[c];
            }
            z[q] = sum / matrix.diagonal[q];
        }
    };
    bicgstab_with(matrix, rhs, tolerance, max_iterations, values, &mut sweep)
}

/// Checks the right-hand side, the start and the tolerance of a solve by
/// the solver `what`.
pub(crate) fn check_system(what: &str, rhs: &[f64], values: &[f64], tolerance: f64) {
    assert!(
        rhs.iter().all(|b| b.is_finite()),
        "{what} needs a finite right-hand side"
    );
    assert!(
        values.iter().all(|v| v.is_finite()),
        "{what} needs finite starting values"
    );
    // Written so that NaN fails.
    assert!(
        tolerance > 0.0 && tolerance < 1.0,
        "{what} needs a tolerance in (0, 1), got {tolerance}"
    );
}

/// The BiCGSTAB method of [`bicgstab`] with the right preconditioner
/// `precondition`, which sets its second argument to `M⁻¹` times its first
/// (`p̂ = M⁻¹ p`, `ŝ = M⁻¹ s` in Barrett et al.); with the diagonal of `A`
/// for `M` it is [`bicgstab`], bit for bit. The inputs are checked by the
/// caller (one entry per row, finite, `tolerance` in `(0, 1)`).
pub(crate) fn bicgstab_with(
    matrix: &SparseMatrix,
    rhs: &[f64],
    tolerance: f64,
    max_iterations: u32,
    values: &mut [f64],
    precondition: &mut dyn FnMut(&[f64], &mut [f64]),
) -> Solve {
    let n = matrix.len();
    let bb = dot(rhs, rhs);
    if bb == 0.0 {
        values.fill(0.0);
        return Solve::new(0, 0.0, true);
    }
    let limit = tolerance * tolerance * bb;
    let x = values;
    let mut r = vec![0.0; n];
    matrix.apply(x, &mut r);
    for (ri, bi) in r.iter_mut().zip(rhs) {
        *ri = bi - *ri;
    }
    let mut shadow = r.clone();
    let mut rr = dot(&r, &r);
    let mut converged = rr <= limit;
    let mut iterations = 0;
    let (mut rho, mut alpha, mut omega) = (1.0_f64, 1.0_f64, 1.0_f64);
    let mut p = vec![0.0; n];
    let mut v = vec![0.0; n];
    let mut p_hat = vec![0.0; n];
    let mut s = vec![0.0; n];
    let mut s_hat = vec![0.0; n];
    let mut t = vec![0.0; n];
    // Whether no step has been taken since the start or the last restart.
    let mut fresh = true;
    while !converged && iterations < max_iterations {
        // Written so that NaN breaks down too.
        let breakdown;
        let rho_next = dot(&shadow, &r);
        if !(rho_next != 0.0 && rho_next.is_finite()) {
            breakdown = true;
        } else {
            if fresh {
                p.copy_from_slice(&r);
            } else {
                let beta = (rho_next / rho) * (alpha / omega);
                for k in 0..n {
                    p[k] = r[k] + beta * (p[k] - omega * v[k]);
                }
            }
            rho = rho_next;
            precondition(&p, &mut p_hat);
            matrix.apply(&p_hat, &mut v);
            let sv = dot(&shadow, &v);
            if !(sv != 0.0 && sv.is_finite()) {
                breakdown = true;
            } else {
                alpha = rho / sv;
                for k in 0..n {
                    s[k] = r[k] - alpha * v[k];
                }
                iterations += 1;
                fresh = false;
                let ss = dot(&s, &s);
                if ss <= limit {
                    for k in 0..n {
                        x[k] += alpha * p_hat[k];
                    }
                    rr = ss;
                    converged = true;
                    break;
                }
                precondition(&s, &mut s_hat);
                matrix.apply(&s_hat, &mut t);
                let tt = dot(&t, &t);
                if !(tt > 0.0 && tt.is_finite()) {
                    for k in 0..n {
                        x[k] += alpha * p_hat[k];
                    }
                    breakdown = true;
                } else {
                    omega = dot(&t, &s) / tt;
                    for k in 0..n {
                        x[k] += alpha * p_hat[k] + omega * s_hat[k];
                        r[k] = s[k] - omega * t[k];
                    }
                    rr = dot(&r, &r);
                    converged = rr <= limit;
                    breakdown = !(omega != 0.0 && omega.is_finite());
                }
            }
        }
        if breakdown && !converged {
            if fresh {
                // No step since the last start: nothing more to gain.
                break;
            }
            // Restart from the current iterate, with its true residual as
            // the residual and the new shadow residual.
            matrix.apply(x, &mut r);
            for (ri, bi) in r.iter_mut().zip(rhs) {
                *ri = bi - *ri;
            }
            rr = dot(&r, &r);
            converged = rr <= limit;
            shadow.copy_from_slice(&r);
            fresh = true;
        }
    }
    let relative_residual = if rr == 0.0 { 0.0 } else { (rr / bb).sqrt() };
    Solve::new(iterations, relative_residual, converged)
}

#[cfg(test)]
mod tests;

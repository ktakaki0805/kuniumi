//! Sparse linear systems that are not symmetric (ticket M4-04): a square
//! matrix stored by rows ([`SparseMatrix`]) and the BiCGSTAB method with
//! the Jacobi preconditioner ([`bicgstab`]).
//!
//! The ocean's Stommel model ([`crate::stages::climate::ocean`]) solves
//! with it the advection–diffusion system of its streamfunction, whose
//! matrix is a non-symmetric M-matrix, and the moisture
//! ([`crate::stages::climate::moisture`]) its upwind advection; the
//! symmetric systems of the diffusion and the sea surface's system are
//! solved with the multigrid of [`crate::multigrid`] (ticket M4-09), whose
//! BiCGSTAB is this one with the V-cycle for the preconditioner.
//!
//! Everything runs sequentially in row order with IEEE `+ − × ÷` and
//! `sqrt`: every product and dot product sums in a fixed order
//! (docs/conventions.md, "Determinism", rules 2 and 6), so the results are
//! the same on every platform.

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

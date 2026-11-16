//! Smoothed-aggregation algebraic multigrid (ticket M4-09): a hierarchy of
//! ever coarser systems that preconditions the conjugate gradient method
//! for the symmetric systems of the climate model ([`conjugate_gradient`])
//! and BiCGSTAB for its non-symmetric one ([`bicgstab`]).
//!
//! # Why
//!
//! The implicit diffusion of the energy balance
//! ([`crate::stages::climate::energy_balance`], ticket M4-02) is stiff at
//! the grid's scale: on 650,000 cells `D L` reaches 10⁵ W m⁻² K⁻¹ per
//! cell against a heat capacity and emission `C/Δt + B` of about 18 for
//! land, and the iterations of the conjugate gradient method with the
//! Jacobi preconditioner ([`crate::operators::Laplacian::solve`]) grow with
//! the resolution (78, 153, 333 and 662 per solve on 10,000, 40,000,
//! 160,000 and 650,000 cells; ticket M4-02). Point smoothers such as Jacobi
//! and Gauss–Seidel remove the error that varies from cell to cell in a
//! few sweeps, but not the smooth error, which diffusion barely damps. A
//! multigrid removes the smooth error on coarser systems, where it varies
//! from cell to cell again, so that its iterations hardly depend on the
//! grid: 12, 13 and 15 per sub-step on 40,000, 160,000 and 650,000 cells
//! (ticket M4-09, prototype).
//!
//! # The method
//!
//! Smoothed aggregation (P. Vaněk, J. Mandel & M. Brezina 1996, "Algebraic
//! multigrid by smoothed aggregation for second and fourth order elliptic
//! problems", Computing 56, 179–196), built from the matrix alone
//! ([`Multigrid::new`]), for the symmetric M-matrices of the climate model
//! (positive diagonal, non-positive entries off it, rows weakly diagonally
//! dominant; on coarser levels small positive entries may appear):
//!
//! - **Strength.** The entry `a_ij` is strong if `a_ij² > θ² a_ii a_jj`,
//!   with `θ` = [`STRENGTH_THRESHOLD`] on the finest level, halved on each
//!   coarser one (Vaněk et al. 1996, §5).
//! - **Aggregation**, in index order: pass 1, a cell not yet aggregated
//!   whose strong neighbours are all free starts an aggregate of itself and
//!   them; pass 2, each remaining cell joins the aggregate of pass 1 of its
//!   strongest neighbour (the first in column order on a tie); pass 3, each
//!   remaining cell with a strong neighbour starts an aggregate of itself
//!   and its free strong neighbours. A cell without strong neighbours,
//!   whose mass dominates its row, belongs to no aggregate: the smoother
//!   alone solves it.
//! - **Prolongator.** The tentative prolongator `P₀` is 1 where a cell
//!   belongs to an aggregate (the constants, the near-null space of a
//!   diffusion), smoothed by one damped Jacobi step of the filtered matrix
//!   `A_F`, which keeps the strong entries and adds the weak ones to the
//!   diagonal: `P = (I − ω D_F⁻¹ A_F) P₀`, `ω = 4/(3ρ)` with `ρ` the
//!   Gershgorin bound `max_i Σ_j |a^F_ij| / a^F_ii` of the spectral radius
//!   of `D_F⁻¹ A_F` (Vaněk et al. 1996, eq. 3.3 with the bound instead of
//!   an eigenvalue estimate, so that nothing is iterated).
//! - **Coarse operators**, `A_c = Pᵀ A P` (Galerkin). Coarsening stops
//!   when a level has at most [`COARSEST_CELLS`] cells, or when its
//!   aggregates would keep more than [`COARSENING_LIMIT`] of them (as where
//!   the masses dominate and nothing is strong).
//! - **V-cycle** ([`Multigrid::apply`]): from zero, one forward
//!   Gauss–Seidel sweep, the residual restricted by `Pᵀ`, the coarse
//!   correction, prolonged by `P`, and one backward Gauss–Seidel sweep. The
//!   backward sweep is the adjoint of the forward one, so for a symmetric
//!   matrix the cycle is a symmetric operator, as the conjugate gradient
//!   method needs. The coarsest level is solved exactly by Gaussian
//!   elimination with partial pivoting where it has at most
//!   [`COARSEST_CELLS`] cells, otherwise (a level that cannot coarsen,
//!   nearly diagonal) by [`COARSEST_SWEEPS`] pairs of forward and backward
//!   sweeps.
//!
//! # Non-symmetric systems
//!
//! The sea surface with currents ([`crate::stages::climate::sea_surface`],
//! ticket M4-05) adds an upwind advection to the energy balance's matrix.
//! [`Multigrid::with_prolongators_of`] takes the aggregates and
//! prolongators of the hierarchy of a symmetric matrix (the energy
//! balance's) and forms the coarse operators `Pᵀ A P` of the non-symmetric
//! one, with the same smoothers: a preconditioner for [`bicgstab`]
//! (10 iterations per sub-step on 650,000 cells against 771 with the
//! Jacobi preconditioner; ticket M4-09, prototype). It does not suit
//! systems dominated by advection, such as the moisture's
//! ([`crate::stages::climate::moisture`]; ticket M4-14).
//!
//! # Determinism
//!
//! Everything runs sequentially in index order with IEEE `+ − × ÷` and
//! comparisons: the aggregation visits the cells in order, each row's
//! columns are sorted as integers, and every product and dot product sums
//! in a fixed order (docs/conventions.md, "Determinism", rules 2 and 6), so
//! the hierarchy and the solutions are the same on every platform.

use crate::operators::Solve;
use crate::sparse::{SparseMatrix, bicgstab_with, check_system};

/// `θ`, the strength of connection on the finest level (Vaněk et al. 1996,
/// §5): 0.08, halved on each coarser level.
pub const STRENGTH_THRESHOLD: f64 = 0.08;

/// The most cells of a level that is not coarsened further and is solved
/// exactly: 400.
pub const COARSEST_CELLS: usize = 400;

/// The largest fraction of a level's cells that its aggregates may number
/// for it to be coarsened: 0.8. Beyond it coarsening gains too little, and
/// the level is the coarsest.
pub const COARSENING_LIMIT: f64 = 0.8;

/// The pairs of forward and backward Gauss–Seidel sweeps that solve a
/// coarsest level of more than [`COARSEST_CELLS`] cells: 4.
pub const COARSEST_SWEEPS: usize = 4;

/// The aggregate of a cell that belongs to none.
pub(crate) const NO_AGGREGATE: u32 = u32::MAX;

/// The aggregates of a level: each cell's, or [`NO_AGGREGATE`], and their
/// number.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Aggregates {
    /// The aggregate of each cell.
    pub(crate) cell: Vec<u32>,
    /// The number of aggregates.
    pub(crate) count: usize,
}

/// A rectangular sparse matrix stored by rows, each row's columns strictly
/// increasing: a prolongator or a restriction.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Rows {
    /// The entries of row `i` are `offsets[i]..offsets[i + 1]`.
    offsets: Vec<usize>,
    /// The column of each entry.
    columns: Vec<u32>,
    /// The value of each entry.
    values: Vec<f64>,
    /// The number of columns.
    width: usize,
}

impl Rows {
    /// The number of rows.
    pub(crate) fn len(&self) -> usize {
        self.offsets.len() - 1
    }

    /// The number of columns.
    pub(crate) fn width(&self) -> usize {
        self.width
    }

    /// The columns and values of row `i`.
    pub(crate) fn row(&self, i: usize) -> (&[u32], &[f64]) {
        let range = self.offsets[i]..self.offsets[i + 1];
        (&self.columns[range.clone()], &self.values[range])
    }

    /// `y = R x`, each row summed in column order.
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        for (i, out) in y.iter_mut().enumerate() {
            let mut sum = 0.0;
            for k in self.offsets[i]..self.offsets[i + 1] {
                sum += self.values[k] * x[self.columns[k] as usize];
            }
            *out = sum;
        }
    }

    /// `y += R x`, each row summed in column order before it is added.
    fn add_product(&self, x: &[f64], y: &mut [f64]) {
        for (i, out) in y.iter_mut().enumerate() {
            let mut sum = 0.0;
            for k in self.offsets[i]..self.offsets[i + 1] {
                sum += self.values[k] * x[self.columns[k] as usize];
            }
            *out += sum;
        }
    }

    /// The transpose, each of its rows in increasing column order.
    pub(crate) fn transpose(&self) -> Rows {
        let mut offsets = vec![0usize; self.width + 1];
        for &c in &self.columns {
            offsets[c as usize + 1] += 1;
        }
        for c in 0..self.width {
            offsets[c + 1] += offsets[c];
        }
        let mut next = offsets.clone();
        let mut columns = vec![0u32; self.columns.len()];
        let mut values = vec![0.0; self.values.len()];
        for i in 0..self.len() {
            for k in self.offsets[i]..self.offsets[i + 1] {
                let c = self.columns[k] as usize;
                columns[next[c]] = i as u32;
                values[next[c]] = self.values[k];
                next[c] += 1;
            }
        }
        Rows {
            offsets,
            columns,
            values,
            width: self.len(),
        }
    }
}

/// Accumulates sparse rows: the sum of each column touched since the last
/// [`start`](Self::start), in the order of the additions. Rows are told
/// apart by their index, so each row index is used once.
struct Accumulator {
    sums: Vec<f64>,
    /// The row each column was last touched for.
    marks: Vec<usize>,
    touched: Vec<u32>,
}

impl Accumulator {
    fn new(width: usize) -> Self {
        Self {
            sums: vec![0.0; width],
            marks: vec![usize::MAX; width],
            touched: Vec::new(),
        }
    }

    /// Starts a row: forgets the columns touched so far (their marks tell
    /// the next row's from the last).
    fn start(&mut self) {
        self.touched.clear();
    }

    /// Adds `value` to column `column` of row `row`.
    fn add(&mut self, row: usize, column: u32, value: f64) {
        let c = column as usize;
        if self.marks[c] != row {
            self.marks[c] = row;
            self.sums[c] = value;
            self.touched.push(column);
        } else {
            self.sums[c] += value;
        }
    }

    /// The touched columns in increasing order with their sums.
    fn finish(&mut self) -> impl Iterator<Item = (u32, f64)> + '_ {
        self.touched.sort_unstable();
        self.touched.iter().map(|&c| (c, self.sums[c as usize]))
    }
}

/// Whether the entry `value` of a row with the diagonal `diagonal_i`, in
/// the column whose diagonal is `diagonal_j`, is strong for `theta`:
/// `a_ij² > θ² a_ii a_jj`.
fn is_strong(value: f64, diagonal_i: f64, diagonal_j: f64, theta: f64) -> bool {
    value * value > theta * theta * diagonal_i * diagonal_j
}

/// Whether `matrix` is symmetric bit for bit: every entry off the
/// diagonal has the bits of its mirror.
pub(crate) fn is_symmetric(matrix: &SparseMatrix) -> bool {
    (0..matrix.len()).all(|i| {
        let (columns, values) = matrix.row(i);
        columns.iter().zip(values).all(|(&j, &v)| {
            let (mirror, mirror_values) = matrix.row(j as usize);
            match mirror.binary_search(&(i as u32)) {
                Ok(k) => mirror_values[k].to_bits() == v.to_bits(),
                Err(_) => false,
            }
        })
    })
}

/// The aggregates of `matrix` for the strength `theta` (module
/// documentation, "The method"): three passes in index order.
pub(crate) fn aggregate(matrix: &SparseMatrix, theta: f64) -> Aggregates {
    let n = matrix.len();
    let d = matrix.diagonal();
    let strong = |i: usize| {
        let (columns, values) = matrix.row(i);
        columns
            .iter()
            .zip(values)
            .filter(move |&(&j, &v)| is_strong(v, d[i], d[j as usize], theta))
            .map(|(&j, &v)| (j as usize, v))
    };
    let mut cell = vec![NO_AGGREGATE; n];
    let mut count = 0u32;
    // Pass 1: cells whose strong neighbours are all free.
    for i in 0..n {
        if cell[i] != NO_AGGREGATE {
            continue;
        }
        let mut any = false;
        let mut free = true;
        for (j, _) in strong(i) {
            any = true;
            if cell[j] != NO_AGGREGATE {
                free = false;
                break;
            }
        }
        if any && free {
            cell[i] = count;
            for (j, _) in strong(i) {
                cell[j] = count;
            }
            count += 1;
        }
    }
    // Pass 2: the aggregate of pass 1 of the strongest neighbour.
    let first = cell.clone();
    for i in 0..n {
        if first[i] != NO_AGGREGATE {
            continue;
        }
        let mut best: Option<(f64, u32)> = None;
        for (j, v) in strong(i) {
            if first[j] != NO_AGGREGATE {
                let size = v.abs();
                if best.is_none_or(|(b, _)| size > b) {
                    best = Some((size, first[j]));
                }
            }
        }
        if let Some((_, a)) = best {
            cell[i] = a;
        }
    }
    // Pass 3: the rest, with their free strong neighbours.
    for i in 0..n {
        if cell[i] != NO_AGGREGATE || strong(i).next().is_none() {
            continue;
        }
        cell[i] = count;
        for (j, _) in strong(i) {
            if cell[j] == NO_AGGREGATE {
                cell[j] = count;
            }
        }
        count += 1;
    }
    Aggregates {
        cell,
        count: count as usize,
    }
}

/// The smoothed prolongator `(I − ω D_F⁻¹ A_F) P₀` of `matrix` for its
/// aggregates `aggregates` and the strength `theta` (module documentation,
/// "The method"): one row per cell, one column per aggregate.
pub(crate) fn prolongator(matrix: &SparseMatrix, theta: f64, aggregates: &Aggregates) -> Rows {
    let n = matrix.len();
    let d = matrix.diagonal();
    // The filtered diagonal and the Gershgorin bound ρ of D_F⁻¹ A_F.
    let mut filtered = Vec::with_capacity(n);
    let mut rho = 0.0_f64;
    for i in 0..n {
        let (columns, values) = matrix.row(i);
        let (mut f, mut s) = (d[i], 0.0);
        for (&j, &v) in columns.iter().zip(values) {
            if is_strong(v, d[i], d[j as usize], theta) {
                s += v.abs();
            } else {
                f += v;
            }
        }
        // Written so that NaN falls back too; not met by the matrices of
        // the climate model, whose rows are dominated by their diagonal.
        if !(f > 0.0 && f.is_finite()) {
            f = d[i];
        }
        let r = (f + s) / f;
        if r > rho {
            rho = r;
        }
        filtered.push(f);
    }
    let omega = 4.0 / (3.0 * rho);
    let mut sums = Accumulator::new(aggregates.count);
    let mut offsets = Vec::with_capacity(n + 1);
    let mut columns = Vec::new();
    let mut values = Vec::new();
    offsets.push(0);
    for i in 0..n {
        sums.start();
        let own = aggregates.cell[i];
        if own != NO_AGGREGATE {
            sums.add(i, own, 1.0 - omega);
        }
        let (row_columns, row_values) = matrix.row(i);
        for (&j, &v) in row_columns.iter().zip(row_values) {
            let a = aggregates.cell[j as usize];
            if a != NO_AGGREGATE && is_strong(v, d[i], d[j as usize], theta) {
                sums.add(i, a, -(omega * v / filtered[i]));
            }
        }
        for (c, v) in sums.finish() {
            columns.push(c);
            values.push(v);
        }
        offsets.push(columns.len());
    }
    Rows {
        offsets,
        columns,
        values,
        width: aggregates.count,
    }
}

/// The Galerkin product `R A P` of `matrix` with the prolongator `p` and
/// the restriction `r` (`Pᵀ`): first `A P` row by row (the diagonal's term,
/// then the others in column order), then `R (A P)`.
pub(crate) fn galerkin(matrix: &SparseMatrix, p: &Rows, r: &Rows) -> SparseMatrix {
    let n = matrix.len();
    let width = p.width();
    let d = matrix.diagonal();
    let mut sums = Accumulator::new(width);
    let mut ap_offsets = Vec::with_capacity(n + 1);
    let mut ap_columns = Vec::new();
    let mut ap_values = Vec::new();
    ap_offsets.push(0);
    for (i, &di) in d.iter().enumerate() {
        sums.start();
        let (columns, values) = p.row(i);
        for (&c, &v) in columns.iter().zip(values) {
            sums.add(i, c, di * v);
        }
        let (row_columns, row_values) = matrix.row(i);
        for (&j, &a) in row_columns.iter().zip(row_values) {
            let (columns, values) = p.row(j as usize);
            for (&c, &v) in columns.iter().zip(values) {
                sums.add(i, c, a * v);
            }
        }
        for (c, v) in sums.finish() {
            ap_columns.push(c);
            ap_values.push(v);
        }
        ap_offsets.push(ap_columns.len());
    }
    let ap = Rows {
        offsets: ap_offsets,
        columns: ap_columns,
        values: ap_values,
        width,
    };
    let mut sums = Accumulator::new(width);
    let mut diagonal = vec![0.0; width];
    let mut rows = Vec::with_capacity(width);
    for (c, slot) in diagonal.iter_mut().enumerate() {
        sums.start();
        let (cells, weights) = r.row(c);
        for (&i, &w) in cells.iter().zip(weights) {
            let (columns, values) = ap.row(i as usize);
            for (&e, &v) in columns.iter().zip(values) {
                sums.add(c, e, w * v);
            }
        }
        let mut row = Vec::new();
        for (e, v) in sums.finish() {
            if e as usize == c {
                *slot = v;
            } else {
                row.push((e, v));
            }
        }
        rows.push(row);
    }
    SparseMatrix::from_rows(diagonal, &rows)
}

/// Panics unless every diagonal entry of `matrix` is finite and positive.
fn check_diagonal(matrix: &SparseMatrix, what: &str) {
    // Written so that NaN fails.
    if let Some(i) = matrix
        .diagonal()
        .iter()
        .position(|d| !(*d > 0.0 && d.is_finite()))
    {
        panic!(
            "{what} needs a finite positive diagonal: row {i} has {}",
            matrix.diagonal()[i]
        );
    }
}

/// Gaussian elimination with partial pivoting of a dense coarsest level.
#[derive(Clone, Debug, PartialEq)]
struct Dense {
    /// The size.
    n: usize,
    /// The factors `L` (unit diagonal, below) and `U` by rows.
    factors: Vec<f64>,
    /// The row swapped with each row in turn.
    pivots: Vec<usize>,
}

impl Dense {
    /// The factorisation of `matrix`, the pivot of each column the largest
    /// magnitude below the diagonal (the first on a tie).
    ///
    /// # Panics
    ///
    /// If a pivot is zero (the level is singular).
    fn new(matrix: &SparseMatrix) -> Self {
        let n = matrix.len();
        let mut a = vec![0.0; n * n];
        for i in 0..n {
            a[i * n + i] = matrix.diagonal()[i];
            let (columns, values) = matrix.row(i);
            for (&j, &v) in columns.iter().zip(values) {
                a[i * n + j as usize] = v;
            }
        }
        let mut pivots = Vec::with_capacity(n);
        for k in 0..n {
            let mut p = k;
            let mut best = a[k * n + k].abs();
            for r in k + 1..n {
                let v = a[r * n + k].abs();
                if v > best {
                    best = v;
                    p = r;
                }
            }
            // Written so that NaN fails.
            assert!(
                best > 0.0 && best.is_finite(),
                "the coarsest level of a multigrid is singular (column {k})"
            );
            if p != k {
                for c in 0..n {
                    a.swap(k * n + c, p * n + c);
                }
            }
            pivots.push(p);
            let pivot = a[k * n + k];
            for r in k + 1..n {
                let factor = a[r * n + k] / pivot;
                a[r * n + k] = factor;
                if factor != 0.0 {
                    for c in k + 1..n {
                        a[r * n + c] -= factor * a[k * n + c];
                    }
                }
            }
        }
        Self {
            n,
            factors: a,
            pivots,
        }
    }

    /// `x = A⁻¹ b`.
    fn solve(&self, b: &[f64], x: &mut [f64]) {
        let n = self.n;
        let a = &self.factors;
        x.copy_from_slice(b);
        for (k, &p) in self.pivots.iter().enumerate() {
            x.swap(k, p);
        }
        for i in 0..n {
            let mut s = x[i];
            for k in 0..i {
                s -= a[i * n + k] * x[k];
            }
            x[i] = s;
        }
        for i in (0..n).rev() {
            let mut s = x[i];
            for k in i + 1..n {
                s -= a[i * n + k] * x[k];
            }
            x[i] = s / a[i * n + i];
        }
    }
}

/// How the coarsest level is solved.
#[derive(Clone, Debug, PartialEq)]
enum Coarsest {
    /// Exactly, at most [`COARSEST_CELLS`] cells.
    Dense(Dense),
    /// By [`COARSEST_SWEEPS`] pairs of Gauss–Seidel sweeps.
    Sweeps,
}

/// A level above the coarsest: its matrix, the prolongator from the next
/// level and the restriction to it.
#[derive(Clone, Debug, PartialEq)]
struct Level {
    matrix: SparseMatrix,
    prolongator: Rows,
    restriction: Rows,
}

/// The vectors of a V-cycle on every level, allocated once per solve.
struct Cycle {
    rhs: Vec<Vec<f64>>,
    solution: Vec<Vec<f64>>,
    residual: Vec<Vec<f64>>,
}

/// The multigrid hierarchy of a sparse matrix (module documentation): its
/// levels from the matrix itself to the coarsest, the prolongators and
/// restrictions between them, and the coarsest level's solver.
#[derive(Clone, Debug, PartialEq)]
pub struct Multigrid {
    levels: Vec<Level>,
    coarsest: SparseMatrix,
    solver: Coarsest,
}

impl Multigrid {
    /// The hierarchy of the symmetric M-matrix `matrix` (module
    /// documentation, "The method"). About 0.3 s for the energy balance's
    /// sub-step on 650,000 cells (ticket M4-09).
    ///
    /// # Panics
    ///
    /// If `matrix` is not symmetric, or a diagonal entry of it or of a
    /// coarse level is not finite and positive.
    pub fn new(matrix: SparseMatrix) -> Self {
        check_diagonal(&matrix, "Multigrid::new");
        assert!(
            is_symmetric(&matrix),
            "Multigrid::new needs a symmetric matrix (bit for bit)"
        );
        let mut levels = Vec::new();
        let mut current = matrix;
        let mut theta = STRENGTH_THRESHOLD;
        while current.len() > COARSEST_CELLS {
            let n = current.len();
            let aggregates = aggregate(&current, theta);
            if aggregates.count == 0 || aggregates.count as f64 > COARSENING_LIMIT * n as f64 {
                break;
            }
            let prolongator = prolongator(&current, theta, &aggregates);
            let restriction = prolongator.transpose();
            let coarse = galerkin(&current, &prolongator, &restriction);
            check_diagonal(&coarse, "a coarse level of the multigrid");
            levels.push(Level {
                matrix: current,
                prolongator,
                restriction,
            });
            current = coarse;
            theta *= 0.5;
        }
        Self::finish(levels, current)
    }

    /// The hierarchy of `matrix` with the aggregates, prolongators and
    /// restrictions of `symmetric` (module documentation, "Non-symmetric
    /// systems"): the coarse operators are `Pᵀ A P` of `matrix`, an
    /// M-matrix whose symmetric part is dominated by that of `symmetric`
    /// (the sea surface's matrix with the energy balance's).
    ///
    /// # Panics
    ///
    /// If `matrix` has not the size of the symmetric hierarchy, or a
    /// diagonal entry of it or of a coarse level is not finite and
    /// positive.
    pub fn with_prolongators_of(matrix: SparseMatrix, symmetric: &Multigrid) -> Self {
        assert!(
            matrix.len() == symmetric.len(),
            "with_prolongators_of needs a matrix of the size of the symmetric hierarchy: {} rows \
             for {}",
            matrix.len(),
            symmetric.len()
        );
        check_diagonal(&matrix, "Multigrid::with_prolongators_of");
        let mut levels = Vec::with_capacity(symmetric.levels.len());
        let mut current = matrix;
        for level in &symmetric.levels {
            let coarse = galerkin(&current, &level.prolongator, &level.restriction);
            check_diagonal(&coarse, "a coarse level of the multigrid");
            levels.push(Level {
                matrix: current,
                prolongator: level.prolongator.clone(),
                restriction: level.restriction.clone(),
            });
            current = coarse;
        }
        Self::finish(levels, current)
    }

    /// The hierarchy of the levels `levels` above the coarsest level
    /// `coarsest`, with its solver.
    fn finish(levels: Vec<Level>, coarsest: SparseMatrix) -> Self {
        let solver = if coarsest.len() <= COARSEST_CELLS {
            Coarsest::Dense(Dense::new(&coarsest))
        } else {
            Coarsest::Sweeps
        };
        Self {
            levels,
            coarsest,
            solver,
        }
    }

    /// The matrix of the finest level, the one the hierarchy was built
    /// from.
    pub fn matrix(&self) -> &SparseMatrix {
        self.levels.first().map_or(&self.coarsest, |l| &l.matrix)
    }

    /// The number of rows of the finest level.
    pub fn len(&self) -> usize {
        self.matrix().len()
    }

    /// Whether the finest level has no rows.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The number of levels, the finest and the coarsest included.
    pub fn levels(&self) -> usize {
        self.levels.len() + 1
    }

    /// The number of rows of each level, from the finest.
    pub fn level_sizes(&self) -> Vec<usize> {
        self.levels
            .iter()
            .map(|l| l.matrix.len())
            .chain(core::iter::once(self.coarsest.len()))
            .collect()
    }

    /// The entries (the diagonal and the others) of every level over those
    /// of the finest, the operator complexity: the memory of the hierarchy
    /// and the work of a V-cycle relative to those of its finest level;
    /// 1.35 to 1.41 for the energy balance's sub-step on 650,000 to 10,000
    /// cells (ticket M4-09).
    pub fn operator_complexity(&self) -> f64 {
        let entries = |m: &SparseMatrix| -> usize {
            m.len() + (0..m.len()).map(|i| m.row(i).0.len()).sum::<usize>()
        };
        let total: usize = self
            .levels
            .iter()
            .map(|l| entries(&l.matrix))
            .sum::<usize>()
            + entries(&self.coarsest);
        total as f64 / entries(self.matrix()) as f64
    }

    /// One V-cycle for `A z = r` from `z = 0` (module documentation, "The
    /// method"): `correction` becomes the cycle's approximation of
    /// `A⁻¹ residual`.
    ///
    /// # Panics
    ///
    /// If `residual` or `correction` does not hold one value per row.
    pub fn apply(&self, residual: &[f64], correction: &mut [f64]) {
        let n = self.len();
        assert!(
            residual.len() == n && correction.len() == n,
            "apply needs one value per row: {} and {} for {n} rows",
            residual.len(),
            correction.len()
        );
        let mut work = self.cycle_space();
        self.cycle(&mut work, self.matrix(), residual, correction);
    }

    /// The vectors of a V-cycle.
    fn cycle_space(&self) -> Cycle {
        let sizes = self.level_sizes();
        Cycle {
            rhs: sizes.iter().map(|&s| vec![0.0; s]).collect(),
            solution: sizes.iter().map(|&s| vec![0.0; s]).collect(),
            residual: sizes[..self.levels.len()]
                .iter()
                .map(|&s| vec![0.0; s])
                .collect(),
        }
    }

    /// The V-cycle of [`apply`](Self::apply) in the vectors `work`, with
    /// the sweeps of the finest level on `finest` (the matrix being solved,
    /// of the hierarchy's size): the hierarchy's own finest matrix, or one
    /// close to it, such as a month's matrix of the sea surface, whose
    /// diagonal differs from the hierarchy's by the month's upwelling.
    fn cycle(
        &self,
        work: &mut Cycle,
        finest: &SparseMatrix,
        residual: &[f64],
        correction: &mut [f64],
    ) {
        let depth = self.levels.len();
        work.rhs[0].copy_from_slice(residual);
        for (l, level) in self.levels.iter().enumerate() {
            let a = if l == 0 { finest } else { &level.matrix };
            let x = &mut work.solution[l];
            x.fill(0.0);
            forward_sweep(a, &work.rhs[l], x);
            residual_of(a, &work.rhs[l], x, &mut work.residual[l]);
            level
                .restriction
                .apply(&work.residual[l], &mut work.rhs[l + 1]);
        }
        let (b, x) = (&work.rhs[depth], &mut work.solution[depth]);
        match &self.solver {
            Coarsest::Dense(dense) => dense.solve(b, x),
            Coarsest::Sweeps => {
                let a = if depth == 0 { finest } else { &self.coarsest };
                x.fill(0.0);
                for _ in 0..COARSEST_SWEEPS {
                    forward_sweep(a, b, x);
                    backward_sweep(a, b, x);
                }
            }
        }
        for (l, level) in self.levels.iter().enumerate().rev() {
            let a = if l == 0 { finest } else { &level.matrix };
            let (fine, coarse) = work.solution.split_at_mut(l + 1);
            level.prolongator.add_product(&coarse[0], &mut fine[l]);
            backward_sweep(a, &work.rhs[l], &mut fine[l]);
        }
        correction.copy_from_slice(&work.solution[0]);
    }
}

/// One Gauss–Seidel sweep of `A x = b` in index order, in place.
fn forward_sweep(a: &SparseMatrix, b: &[f64], x: &mut [f64]) {
    let d = a.diagonal();
    for i in 0..a.len() {
        let (columns, values) = a.row(i);
        let mut s = b[i];
        for (&j, &v) in columns.iter().zip(values) {
            s -= v * x[j as usize];
        }
        x[i] = s / d[i];
    }
}

/// One Gauss–Seidel sweep of `A x = b` in reverse index order, in place.
fn backward_sweep(a: &SparseMatrix, b: &[f64], x: &mut [f64]) {
    let d = a.diagonal();
    for i in (0..a.len()).rev() {
        let (columns, values) = a.row(i);
        let mut s = b[i];
        for (&j, &v) in columns.iter().zip(values) {
            s -= v * x[j as usize];
        }
        x[i] = s / d[i];
    }
}

/// `r = b − A x`, each product summed as [`SparseMatrix::apply`] does.
fn residual_of(a: &SparseMatrix, b: &[f64], x: &[f64], r: &mut [f64]) {
    let d = a.diagonal();
    for i in 0..a.len() {
        let (columns, values) = a.row(i);
        let mut sum = d[i] * x[i];
        for (&j, &v) in columns.iter().zip(values) {
            sum += v * x[j as usize];
        }
        r[i] = b[i] - sum;
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

/// Panics unless `multigrid`, `rhs` and `values` fit `matrix` and the
/// diagonal of `matrix` is finite and positive (the Gauss–Seidel sweeps
/// divide by it), as `what`.
fn check_sizes(
    what: &str,
    matrix: &SparseMatrix,
    multigrid: &Multigrid,
    rhs: &[f64],
    values: &[f64],
) {
    let n = matrix.len();
    assert!(
        multigrid.len() == n,
        "{what} needs the multigrid of a matrix of its size: {} rows for {n}",
        multigrid.len()
    );
    assert!(
        rhs.len() == n,
        "{what} needs one right-hand side per row: {} for {n} rows",
        rhs.len()
    );
    assert!(
        values.len() == n,
        "{what} needs one value per row: {} values for {n} rows",
        values.len()
    );
    check_diagonal(matrix, what);
}

/// Solves the symmetric positive definite system `A x = b` by the
/// conjugate gradient method preconditioned by one V-cycle of `multigrid`
/// (J. R. Shewchuk 1994, "An Introduction to the Conjugate Gradient Method
/// Without the Agonizing Pain", Carnegie Mellon University, appendix B3,
/// with the V-cycle for `M⁻¹`); `A = matrix`, `b = rhs`, `values` hold the
/// start and become the solution. `multigrid` is the hierarchy of `matrix`
/// or of a symmetric matrix close to it: the cycle's sweeps on the finest
/// level are those of `matrix` itself, its coarser levels the hierarchy's
/// (a hierarchy of a single level of at most [`COARSEST_CELLS`] rows is
/// solved exactly, so its cycle is the inverse of its own matrix).
///
/// It stops as [`crate::operators::Laplacian::solve`] does: when the
/// residual it updates falls to `‖b − A x‖₂ ≤ tolerance · ‖b‖₂`, after
/// `max_iterations` iterations, or if rounding leaves no direction to go
/// (`p·Ap` or `r·z` not positive and finite); the returned [`Solve`] says
/// which. A start that already meets the tolerance is returned unchanged,
/// bit for bit; a zero right-hand side gives zero.
///
/// # Panics
///
/// If `multigrid`, `rhs` or `values` does not fit `matrix`, if a diagonal
/// entry of `matrix` is not finite and positive, if a right-hand side or a
/// value is not finite, if `tolerance` is not in `(0, 1)`.
pub fn conjugate_gradient(
    matrix: &SparseMatrix,
    multigrid: &Multigrid,
    rhs: &[f64],
    tolerance: f64,
    max_iterations: u32,
    values: &mut [f64],
) -> Solve {
    check_sizes("conjugate_gradient", matrix, multigrid, rhs, values);
    check_system("conjugate_gradient", rhs, values, tolerance);
    let n = matrix.len();
    let bb = dot(rhs, rhs);
    if bb == 0.0 {
        values.fill(0.0);
        return Solve::new(0, 0.0, true);
    }
    let limit = tolerance * tolerance * bb;
    let mut r = vec![0.0; n];
    residual_of(matrix, rhs, values, &mut r);
    let mut rr = dot(&r, &r);
    let mut converged = rr <= limit;
    let mut iterations = 0;
    if !converged {
        let mut work = multigrid.cycle_space();
        let mut z = vec![0.0; n];
        multigrid.cycle(&mut work, matrix, &r, &mut z);
        let mut rz = dot(&r, &z);
        let mut p = z.clone();
        let mut q = vec![0.0; n];
        // Written so that NaN stops too.
        while rz > 0.0 && rz.is_finite() && iterations < max_iterations {
            matrix.apply(&p, &mut q);
            let pq = dot(&p, &q);
            if !(pq > 0.0 && pq.is_finite()) {
                break;
            }
            let alpha = rz / pq;
            for k in 0..n {
                values[k] += alpha * p[k];
                r[k] -= alpha * q[k];
            }
            iterations += 1;
            rr = dot(&r, &r);
            converged = rr <= limit;
            if converged {
                break;
            }
            multigrid.cycle(&mut work, matrix, &r, &mut z);
            let rz_next = dot(&r, &z);
            let beta = rz_next / rz;
            rz = rz_next;
            for k in 0..n {
                p[k] = z[k] + beta * p[k];
            }
        }
    }
    let relative_residual = if rr == 0.0 { 0.0 } else { (rr / bb).sqrt() };
    Solve::new(iterations, relative_residual, converged)
}

/// Solves `A x = b` by BiCGSTAB ([`crate::sparse::bicgstab`], the same
/// method, stops, breakdowns and restarts) right-preconditioned by one
/// V-cycle of `multigrid` instead of the diagonal; `A = matrix`,
/// `b = rhs`, `values` hold the start and become the solution.
/// `multigrid` is the hierarchy of `matrix` or of a matrix close to it
/// ([`Multigrid::with_prolongators_of`]): the cycle's sweeps on the finest
/// level are those of `matrix` itself, its coarser levels the hierarchy's
/// (the sea surface uses one hierarchy for the twelve months, whose
/// matrices differ on the diagonal by the month's upwelling; a hierarchy of
/// a single level of at most [`COARSEST_CELLS`] rows is solved exactly, so
/// its cycle is the inverse of its own matrix).
///
/// # Panics
///
/// If `multigrid`, `rhs` or `values` does not fit `matrix`, if a diagonal
/// entry of `matrix` is not finite and positive, if a right-hand side or a
/// value is not finite, if `tolerance` is not in `(0, 1)`.
pub fn bicgstab(
    matrix: &SparseMatrix,
    multigrid: &Multigrid,
    rhs: &[f64],
    tolerance: f64,
    max_iterations: u32,
    values: &mut [f64],
) -> Solve {
    check_sizes("bicgstab", matrix, multigrid, rhs, values);
    check_system("bicgstab", rhs, values, tolerance);
    let mut work = multigrid.cycle_space();
    let mut precondition = |r: &[f64], z: &mut [f64]| multigrid.cycle(&mut work, matrix, r, z);
    bicgstab_with(
        matrix,
        rhs,
        tolerance,
        max_iterations,
        values,
        &mut precondition,
    )
}

#[cfg(test)]
mod tests;

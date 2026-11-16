//! Finite-volume operators on the spherical grid (design: "侵食と地形発達 /
//! 基本方程式", the term `K_d ∇²h`; ticket M3-04).
//!
//! [`Laplacian`] is the finite-volume Laplacian of a [`SphereGrid`]: the
//! flux across the Voronoi edge between cells `i` and `j` is
//! `κ (u_j − u_i) / (d_ij R) · ℓ_ij R`, with `ℓ_ij` the length of the edge
//! and `d_ij` the distance between the two centres (radians on the unit
//! sphere, [`SphereGrid::edge_lengths_rad`] and
//! [`SphereGrid::centre_distances_rad`]) and `R` the radius, and cell `i`
//! changes by the sum of its fluxes over its area `a_i R²`
//! ([`SphereGrid::area_sr`]):
//!
//! ```text
//! (∇²u)_i = Σ_j w_ij (u_j − u_i) / (a_i R²),   w_ij = ℓ_ij / d_ij.
//! ```
//!
//! This is the two-point flux of an orthogonal (Voronoi) mesh: the Voronoi
//! edge is perpendicular to the arc between the two centres, so the
//! difference of the centre values is the normal derivative across it.
//! The weights are the same, bit for bit, in both directions, so every
//! flux leaves one cell as it enters the other.
//!
//! Its accuracy (ticket M3-04, measured on the jittered grid of
//! [`SphereGrid::fibonacci`] at 10⁴, 10⁵ and 10⁶ cells):
//!
//! - Pointwise, [`Laplacian::apply`] is not consistent in general. Where
//!   the centre is not the centroid of its cell, the two-point flux misses
//!   part of the curvature of a field whose second derivatives differ by
//!   direction, and that error does not shrink as the grid is refined: for
//!   the degree-2 harmonics `3z² − 1` and `xy` it stays at about 3 % (RMS)
//!   and 17–24 % (largest) of the largest `|∇²u|` at every size. Only
//!   where the Hessian is isotropic, as for the degree-1 harmonic `z`
//!   (`−z` times the identity on the sphere), does the error fall with the
//!   spacing: 1.2 %, 0.40 % and 0.13 % of `2/R²` for `∇²z = −2z/R²`.
//! - The solutions of [`Laplacian::implicit_step`] converge nonetheless:
//!   for both degree-2 harmonics, one step of `c = Δt κ / R² = 1` lies
//!   within 1e-4, 1e-5 and 1e-6 of `u / (1 + 6c)`, relative to the largest
//!   `|u|`, so at second order in the spacing, as measured (the error
//!   bounds that R. Eymard, T. Gallouët and R. Herbin 2000, "Finite volume
//!   methods", *Handbook of Numerical Analysis* VII, prove for the
//!   two-point flux scheme on admissible, orthogonal meshes are of first
//!   order).
//!
//! So the operator serves the implicit diffusion; a stage that needs `∇²u`
//! cell by cell, as a source term, needs a consistent operator or a
//! centroidal grid.
//!
//! [`Laplacian::implicit_step`] advances `∂u/∂t = κ ∇²u` by one
//! backward-Euler step, solved by the conjugate gradient method, on some
//! cells while the others hold fixed values (a Dirichlet boundary): the
//! hillslope diffusion of the Erosion stage
//! ([`crate::stages::erosion::diffuse`]), and later the diffusion of
//! marine sediment. [`Laplacian::solve`] solves the more general
//! `(M + c L) u = b` with a mass per cell by the same method, over every
//! cell: the diffusion of heat in the energy balance of the Climate stage
//! ([`crate::stages::climate::energy_balance`], ticket M4-02).
//!
//! [`Gradient`] is the least-squares gradient of a cell field in each
//! cell's east/north frame (ticket M4-03): the pressure gradient that
//! drives the surface winds of the Climate stage
//! ([`crate::stages::climate::circulation`]), and the slope that lifts
//! moist air ([`crate::stages::climate::moisture`], ticket M4-06). (The
//! curl of the wind stress of the ocean, [`crate::stages::climate::ocean`],
//! is the circulation around each cell instead, ticket M4-04.) The
//! non-symmetric systems of the ocean are solved by
//! [`crate::sparse::bicgstab`].
//!
//! Everything here depends only on its inputs, through IEEE `+ − × ÷` and
//! `sqrt` in a fixed order: every sum and dot product runs sequentially over
//! the cells in index order (docs/conventions.md, "Determinism", rules 2
//! and 6), so the results are the same on every platform.

use kuniumi_sphere::{SphereGrid, coords};

/// The relative residual at which [`Laplacian::implicit_step`] stops:
/// `‖b − A u'‖₂ ≤ 10⁻¹² ‖b‖₂`. [`Laplacian::solve`] takes its tolerance
/// from the caller.
pub const SOLVER_TOLERANCE: f64 = 1.0e-12;

/// The number of iterations after which [`Laplacian::implicit_step`] and
/// [`Laplacian::solve`] stop whether or not they have reached their
/// tolerance: 10,000. The step needs one or two iterations for the
/// hillslope diffusion of the Erosion stage, about 15 at the explicit limit
/// and a few hundred to a few thousand at 10⁴ to 10¹² times it (10⁴ and
/// 10⁵ cells), 3,600 at 7 × 10⁵ times it on 10⁶ cells (ticket M3-04).
pub const SOLVER_MAX_ITERATIONS: u32 = 10_000;

/// The finite-volume Laplacian of a grid (module documentation): per
/// neighbour pair the weight `w_ij = ℓ_ij / d_ij` (dimensionless), per cell
/// its area `a_i` (sr) and `Σ_j w_ij`, in the neighbour order of the grid.
#[derive(Clone, Debug, PartialEq)]
pub struct Laplacian {
    /// CSR offsets: the neighbours of cell `i` are
    /// `offsets[i]..offsets[i + 1]`.
    offsets: Vec<u32>,
    /// The neighbours, as in the grid.
    neighbors: Vec<u32>,
    /// `w_ij` per CSR entry.
    weights: Vec<f64>,
    /// `a_i` (sr) per cell.
    areas_sr: Vec<f64>,
    /// `Σ_j w_ij` per cell, summed in neighbour order.
    weight_sums: Vec<f64>,
}

/// What one [`Laplacian::implicit_step`], [`Laplacian::solve`] or
/// [`crate::sparse::bicgstab`] did.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Solve {
    /// The iterations of the conjugate gradient method (of BiCGSTAB for
    /// [`crate::sparse::bicgstab`]).
    pub iterations: u32,
    /// The relative residual `‖b − A u'‖₂ / ‖b‖₂` at the end, as the
    /// method updates it (0 where nothing was solved, `+∞` where `b` is 0
    /// but the residual is not).
    pub relative_residual: f64,
    /// Whether the solve stopped on the tolerance (or had nothing to
    /// solve), as the loop decided it.
    converged: bool,
}

impl Solve {
    /// The report of a solve, for the other solvers of the crate
    /// ([`crate::sparse::bicgstab`], ticket M4-04).
    pub(crate) fn new(iterations: u32, relative_residual: f64, converged: bool) -> Self {
        Self {
            iterations,
            relative_residual,
            converged,
        }
    }

    /// Whether the solve stopped because the residual reached its
    /// tolerance ([`SOLVER_TOLERANCE`] for
    /// [`implicit_step`](Laplacian::implicit_step), the caller's for
    /// [`solve`](Laplacian::solve) and [`crate::sparse::bicgstab`]), or had
    /// nothing to solve; `false` where it stopped on its iteration limit
    /// ([`SOLVER_MAX_ITERATIONS`], or the caller's for
    /// [`crate::sparse::bicgstab`]) or for lack of a direction.
    pub fn converged(&self) -> bool {
        self.converged
    }
}

impl Laplacian {
    /// The Laplacian of `grid`: the weights `ℓ_ij / d_ij` of every
    /// neighbour pair (`d_ij > 0` on a valid grid), the areas and the sums
    /// of the weights. O(n).
    pub fn new(grid: &SphereGrid) -> Self {
        let n = grid.len();
        let mut offsets = Vec::with_capacity(n + 1);
        let mut neighbors = Vec::new();
        let mut weights = Vec::new();
        let mut areas_sr = Vec::with_capacity(n);
        let mut weight_sums = Vec::with_capacity(n);
        offsets.push(0);
        for i in 0..n as u32 {
            let mut sum = 0.0;
            let lengths = grid.edge_lengths_rad(i);
            let distances = grid.centre_distances_rad(i);
            for (k, &j) in grid.neighbors(i).iter().enumerate() {
                let w = lengths[k] / distances[k];
                neighbors.push(j);
                weights.push(w);
                sum += w;
            }
            // The grid's own CSR offsets are u32 (at most 6n − 12 entries,
            // `kuniumi_sphere::MAX_CELLS`), so these fit.
            offsets.push(neighbors.len() as u32);
            areas_sr.push(grid.area_sr(i));
            weight_sums.push(sum);
        }
        Self {
            offsets,
            neighbors,
            weights,
            areas_sr,
            weight_sums,
        }
    }

    /// The number of cells.
    pub fn len(&self) -> usize {
        self.areas_sr.len()
    }

    /// Whether there are no cells.
    pub fn is_empty(&self) -> bool {
        self.areas_sr.is_empty()
    }

    /// The CSR range of cell `i`.
    fn range(&self, i: usize) -> core::ops::Range<usize> {
        self.offsets[i] as usize..self.offsets[i + 1] as usize
    }

    /// The weights `w_ij = ℓ_ij / d_ij` of the neighbours of cell `i`, in
    /// the order of [`SphereGrid::neighbors`]; the same, bit for bit, as the
    /// weight listed for `i` among the neighbours of `j`.
    pub fn weights(&self, i: u32) -> &[f64] {
        &self.weights[self.range(i as usize)]
    }

    /// `Σ_j w_ij (u_j − u_i)` of cell `i` with `u = values`, in neighbour
    /// order.
    fn flux_sum(&self, i: usize, values: &[f64]) -> f64 {
        let u = values[i];
        let mut sum = 0.0;
        for k in self.range(i) {
            sum += self.weights[k] * (values[self.neighbors[k] as usize] - u);
        }
        sum
    }

    /// The Laplacian of `values` (one per cell, in a unit `[u]`) on a
    /// sphere of radius `radius_m`, in `[u]`/m²:
    /// `(∇²u)_i = Σ_j w_ij (u_j − u_i) / (a_i R²)` (module documentation),
    /// the sum in neighbour order. A constant gives exactly 0.
    ///
    /// Cell by cell this is not a consistent approximation of `∇²u` on an
    /// irregular grid: for a field whose second derivatives differ by
    /// direction its error does not shrink with refinement (about 3 % RMS
    /// on the jittered global grid; module documentation). Sums of it over
    /// many cells, and the solutions of [`implicit_step`](Self::implicit_step),
    /// do converge.
    ///
    /// # Panics
    ///
    /// If `values` does not hold one value per cell, or if `radius_m` is
    /// not finite and positive.
    pub fn apply(&self, radius_m: f64, values: &[f64]) -> Vec<f64> {
        let n = self.len();
        assert_eq!(
            values.len(),
            n,
            "apply needs one value per cell: {} values for {n} cells",
            values.len()
        );
        // Written so that NaN fails.
        assert!(
            radius_m > 0.0 && radius_m.is_finite(),
            "apply needs a finite positive radius, got {radius_m}"
        );
        let radius_sq_m2 = radius_m * radius_m;
        (0..n)
            .map(|i| self.flux_sum(i, values) / (self.areas_sr[i] * radius_sq_m2))
            .collect()
    }

    /// One backward-Euler step of the diffusion `∂u/∂t = κ ∇²u` on the
    /// cells where `free[i]`, in place: `values` (one per cell, in any unit)
    /// are the values before the step and become those after it; the other
    /// cells are fixed (a Dirichlet boundary) and keep their values, which
    /// their free neighbours see. `κ = diffusivity_m2_s` (m²/s),
    /// `Δt = time_step_s` (s), `R = radius_m` (m).
    ///
    /// With `c = Δt κ / R²` (dimensionless), each free cell satisfies
    ///
    /// ```text
    /// a_i u'_i + c Σ_j w_ij (u'_i − u'_j) = a_i u_i,
    /// ```
    ///
    /// with `u'_j` the fixed value of a fixed neighbour: the system
    /// `A u' = b`, `A = M + c L`, with `M` the diagonal of the areas and `L`
    /// the weighted graph Laplacian of the free cells (its diagonal
    /// `Σ_j w_ij` over every neighbour), and
    /// `b_i = a_i u_i + c Σ_{j fixed} w_ij u_j`. `A` is symmetric and
    /// positive definite (the weights are symmetric bit for bit), and an
    /// M-matrix: in exact arithmetic `u'` lies within the range of `u` and
    /// the fixed values, for every `Δt` (a discrete maximum principle), so
    /// the step is stable however large `c` is.
    ///
    /// The solver is the conjugate gradient method preconditioned by the
    /// diagonal of `A` (J. R. Shewchuk 1994, "An Introduction to the
    /// Conjugate Gradient Method Without the Agonizing Pain", Carnegie
    /// Mellon University, appendix B3), started from `u`, over the free
    /// cells in index order. The system is divided by `1 + c` so that its
    /// entries stay of the order of the areas and weights for any finite
    /// `c`. It stops when the residual it updates falls to
    /// `‖b − A u'‖₂ ≤` [`SOLVER_TOLERANCE`]` · ‖b‖₂`, after
    /// [`SOLVER_MAX_ITERATIONS`] iterations, or if rounding leaves no
    /// direction to go; the returned [`Solve`] says which.
    ///
    /// Where no free cell has a fixed neighbour (on a connected grid: all
    /// cells free) the step conserves `Σ a_i u_i` in exact arithmetic, but
    /// the solver's tolerance alone does not keep it within 1e-12 for every
    /// field (a lone spike loses 1e-11 of it; ticket M3-04). There the
    /// solved values are shifted by the constant
    /// `(Σ a_i u_i − Σ a_i u'_i) / Σ a_i` (compensated sums in index order),
    /// which restores the sum up to rounding. The constant vector is the
    /// eigenvector of `M⁻¹A` (of the pencil `(A, M)`) with the smallest
    /// eigenvalue, 1, along which the solver's error collects; the shift
    /// removes the error's component along it in the `M` inner product, so
    /// it can only reduce the error's `M`-norm, and it changes no
    /// difference between cells.
    ///
    /// With `κ = 0`, `Δt = 0`, no free cell, or a `c` that rounds to 0, the
    /// values are unchanged bit for bit. Values far beyond 10¹⁵⁰ in
    /// magnitude can overflow the dot products.
    ///
    /// # Panics
    ///
    /// If `values` or `free` does not hold one entry per cell; if
    /// `radius_m` is not finite and positive; if the diffusivity or the
    /// time step is negative or not finite; if a value is not finite; if
    /// `c = Δt κ / R²` is not finite.
    pub fn implicit_step(
        &self,
        radius_m: f64,
        diffusivity_m2_s: f64,
        time_step_s: f64,
        free: &[bool],
        values: &mut [f64],
    ) -> Solve {
        let n = self.len();
        assert_eq!(
            values.len(),
            n,
            "implicit_step needs one value per cell: {} values for {n} cells",
            values.len()
        );
        assert_eq!(
            free.len(),
            n,
            "implicit_step needs one flag per cell: {} flags for {n} cells",
            free.len()
        );
        // Written so that NaN fails.
        assert!(
            radius_m > 0.0 && radius_m.is_finite(),
            "implicit_step needs a finite positive radius, got {radius_m}"
        );
        assert!(
            diffusivity_m2_s >= 0.0 && diffusivity_m2_s.is_finite(),
            "implicit_step needs a finite, non-negative diffusivity, got {diffusivity_m2_s}"
        );
        assert!(
            time_step_s >= 0.0 && time_step_s.is_finite(),
            "implicit_step needs a finite, non-negative time step, got {time_step_s}"
        );
        if let Some(i) = values.iter().position(|v| !v.is_finite()) {
            panic!(
                "implicit_step needs finite values: cell {i} has {}",
                values[i]
            );
        }
        let nothing = Solve {
            iterations: 0,
            relative_residual: 0.0,
            converged: true,
        };
        let cells: Vec<usize> = (0..n).filter(|&i| free[i]).collect();
        if diffusivity_m2_s == 0.0 || time_step_s == 0.0 || cells.is_empty() {
            return nothing;
        }
        let c = time_step_s * diffusivity_m2_s / (radius_m * radius_m);
        assert!(
            c.is_finite(),
            "implicit_step needs a finite step Δt κ / R², got {c}"
        );
        if c == 0.0 {
            return nothing;
        }

        // The system divided by 1 + c: A' = M/(1 + c) + c/(1 + c) L.
        let scale = 1.0 + c;
        let coupling = c / scale;
        let mass: Vec<f64> = self.areas_sr.iter().map(|&a| a / scale).collect();
        let diagonal: Vec<f64> = (0..n)
            .map(|i| mass[i] + coupling * self.weight_sums[i])
            .collect();

        // The initial residual b − A'u = c/(1 + c) Σ_j w_ij (u_j − u_i),
        // computed as fluxes so that nothing cancels, and ‖b‖₂².
        let u: &[f64] = values;
        let mut residual = vec![0.0; n];
        let (mut bb, mut closed) = (0.0, true);
        for &i in &cells {
            residual[i] = coupling * self.flux_sum(i, u);
            let mut boundary = 0.0;
            for k in self.range(i) {
                let j = self.neighbors[k] as usize;
                if !free[j] {
                    boundary += self.weights[k] * u[j];
                    closed = false;
                }
            }
            let b = mass[i] * u[i] + coupling * boundary;
            bb += b * b;
        }
        let limit = SOLVER_TOLERANCE * SOLVER_TOLERANCE * bb;

        let Descent {
            change,
            iterations,
            rr,
            converged,
        } = self.conjugate_gradient(&cells, &mass, coupling, &diagonal, residual, limit);
        let relative_residual = if rr == 0.0 { 0.0 } else { (rr / bb).sqrt() };

        if closed {
            let before = compensated_sum(cells.iter().map(|&i| self.areas_sr[i] * values[i]));
            for &i in &cells {
                values[i] += change[i];
            }
            let after = compensated_sum(cells.iter().map(|&i| self.areas_sr[i] * values[i]));
            let area = compensated_sum(cells.iter().map(|&i| self.areas_sr[i]));
            let shift = (before - after) / area;
            for &i in &cells {
                values[i] += shift;
            }
        } else {
            for &i in &cells {
                values[i] += change[i];
            }
        }
        Solve {
            iterations,
            relative_residual,
            converged,
        }
    }

    /// Solves `(M + c L) u = b` over every cell, in place: `M` the
    /// diagonal of `mass` (one positive value per cell, in a unit `[m]`),
    /// `c = coupling` (in `[m]`, as each weight `w_ij` is dimensionless),
    /// `L` the weighted graph Laplacian, `(L u)_i = Σ_j w_ij (u_i − u_j)`
    /// in neighbour order, and `b = rhs` (in `[m][u]`); `values` hold the
    /// starting guess and become the solution. With `mass` the areas `a_i`
    /// and `c = Δt κ / R²` this is the step of
    /// [`implicit_step`](Self::implicit_step) with every cell free; with
    /// other masses it is the backward-Euler step of a diffusion with a
    /// heat capacity and a linear loss per cell, as in the energy balance
    /// of the Climate stage
    /// ([`crate::stages::climate::energy_balance`], ticket M4-02).
    ///
    /// `A = M + c L` is symmetric and positive definite (the weights are
    /// symmetric bit for bit and every mass is positive). The solver is the
    /// conjugate gradient method preconditioned by the diagonal of `A`,
    /// the same code as in [`implicit_step`](Self::implicit_step), over the
    /// cells in index order, started from `values`. It stops when the
    /// residual it updates falls to `‖b − A u‖₂ ≤ tolerance · ‖b‖₂`, after
    /// [`SOLVER_MAX_ITERATIONS`] iterations, or if rounding leaves no
    /// direction to go; the returned [`Solve`] says which. The initial
    /// residual is `b_i − m_i u_i + c Σ_j w_ij (u_j − u_i)`, so that a
    /// starting guess that solves the system exactly (a constant `u` with
    /// `b = M u`, for one) is returned unchanged, bit for bit.
    ///
    /// # Panics
    ///
    /// If `mass`, `rhs` or `values` does not hold one entry per cell; if a
    /// mass is not finite and positive; if `coupling` is negative or not
    /// finite; if a right-hand side or a value is not finite; if
    /// `tolerance` is not in `(0, 1)`.
    pub fn solve(
        &self,
        mass: &[f64],
        coupling: f64,
        rhs: &[f64],
        tolerance: f64,
        values: &mut [f64],
    ) -> Solve {
        let n = self.len();
        assert_eq!(
            mass.len(),
            n,
            "solve needs one mass per cell: {} masses for {n} cells",
            mass.len()
        );
        assert_eq!(
            rhs.len(),
            n,
            "solve needs one right-hand side per cell: {} for {n} cells",
            rhs.len()
        );
        assert_eq!(
            values.len(),
            n,
            "solve needs one value per cell: {} values for {n} cells",
            values.len()
        );
        // Written so that NaN fails.
        if let Some(i) = mass.iter().position(|m| !(*m > 0.0 && m.is_finite())) {
            panic!(
                "solve needs finite positive masses: cell {i} has {}",
                mass[i]
            );
        }
        assert!(
            coupling >= 0.0 && coupling.is_finite(),
            "solve needs a finite, non-negative coupling, got {coupling}"
        );
        if let Some(i) = rhs.iter().position(|b| !b.is_finite()) {
            panic!(
                "solve needs a finite right-hand side: cell {i} has {}",
                rhs[i]
            );
        }
        if let Some(i) = values.iter().position(|v| !v.is_finite()) {
            panic!("solve needs finite values: cell {i} has {}", values[i]);
        }
        assert!(
            tolerance > 0.0 && tolerance < 1.0,
            "solve needs a tolerance in (0, 1), got {tolerance}"
        );
        let cells: Vec<usize> = (0..n).collect();
        let diagonal: Vec<f64> = (0..n)
            .map(|i| mass[i] + coupling * self.weight_sums[i])
            .collect();
        let u: &[f64] = values;
        let mut residual = vec![0.0; n];
        let mut bb = 0.0;
        for i in 0..n {
            residual[i] = rhs[i] - mass[i] * u[i] + coupling * self.flux_sum(i, u);
            bb += rhs[i] * rhs[i];
        }
        let limit = tolerance * tolerance * bb;
        let Descent {
            change,
            iterations,
            rr,
            converged,
        } = self.conjugate_gradient(&cells, mass, coupling, &diagonal, residual, limit);
        let relative_residual = if rr == 0.0 { 0.0 } else { (rr / bb).sqrt() };
        for i in 0..n {
            values[i] += change[i];
        }
        Solve {
            iterations,
            relative_residual,
            converged,
        }
    }

    /// The preconditioned conjugate gradient method of
    /// [`implicit_step`](Self::implicit_step) and [`solve`](Self::solve)
    /// on `cells` (in their order) for `A = diag(mass) + coupling L`,
    /// whose entries outside `cells` are held fixed: from the residual
    /// `residual` of the starting guess, with `diagonal` the diagonal of
    /// `A` (the Jacobi preconditioner), until `‖r‖₂² ≤ limit`, after
    /// [`SOLVER_MAX_ITERATIONS`] iterations, or when rounding leaves no
    /// direction to go (`p·Ap` or `r·z` not positive and finite).
    fn conjugate_gradient(
        &self,
        cells: &[usize],
        mass: &[f64],
        coupling: f64,
        diagonal: &[f64],
        mut residual: Vec<f64>,
        limit: f64,
    ) -> Descent {
        let n = self.len();
        let mut change = vec![0.0; n];
        let mut z = vec![0.0; n];
        for &i in cells {
            z[i] = residual[i] / diagonal[i];
        }
        let mut direction = z.clone();
        let mut image = vec![0.0; n];
        let mut rz = dot(cells, &residual, &z);
        let mut rr = dot(cells, &residual, &residual);
        let mut iterations = 0;
        let mut converged = rr <= limit;
        while !converged && iterations < SOLVER_MAX_ITERATIONS {
            // image = A direction; the entries of `direction` outside
            // `cells` (fixed cells) are 0.
            for &i in cells {
                let p = direction[i];
                let mut sum = 0.0;
                for k in self.range(i) {
                    sum += self.weights[k] * (p - direction[self.neighbors[k] as usize]);
                }
                image[i] = mass[i] * p + coupling * sum;
            }
            let pq = dot(cells, &direction, &image);
            // Written so that NaN stops too.
            if !(pq > 0.0 && pq.is_finite()) {
                break;
            }
            let alpha = rz / pq;
            for &i in cells {
                change[i] += alpha * direction[i];
                residual[i] -= alpha * image[i];
            }
            iterations += 1;
            rr = dot(cells, &residual, &residual);
            converged = rr <= limit;
            if converged {
                break;
            }
            for &i in cells {
                z[i] = residual[i] / diagonal[i];
            }
            let rz_next = dot(cells, &residual, &z);
            if !(rz_next > 0.0 && rz_next.is_finite()) {
                break;
            }
            let beta = rz_next / rz;
            rz = rz_next;
            for &i in cells {
                direction[i] = z[i] + beta * direction[i];
            }
        }
        Descent {
            change,
            iterations,
            rr,
            converged,
        }
    }
}

/// What [`Laplacian::conjugate_gradient`] returns: the change to add to the
/// starting guess, the iterations, the last `‖r‖₂²` and whether it stopped
/// on the limit.
struct Descent {
    change: Vec<f64>,
    iterations: u32,
    rr: f64,
    converged: bool,
}

/// `Σ_k x_k y_k` over `cells`, sequentially in their order.
fn dot(cells: &[usize], x: &[f64], y: &[f64]) -> f64 {
    let mut sum = 0.0;
    for &i in cells {
        sum += x[i] * y[i];
    }
    sum
}

/// The sum of `terms` in their order, compensated (A. Neumaier 1974,
/// "Rundungsfehleranalyse einiger Verfahren zur Summation endlicher
/// Summen", *ZAMM* 54, 39–51): the rounding error of each addition is
/// collected exactly (`+ −` and comparisons only) and added at the end, so
/// the result is accurate to about one rounding of the sum.
fn compensated_sum(terms: impl Iterator<Item = f64>) -> f64 {
    let (mut sum, mut compensation) = (0.0_f64, 0.0_f64);
    for x in terms {
        let t = sum + x;
        compensation += if sum.abs() >= x.abs() {
            (sum - t) + x
        } else {
            (x - t) + sum
        };
        sum = t;
    }
    sum + compensation
}

/// The gradient of a cell field on the unit sphere, by least squares over
/// each cell's neighbours (module documentation; ticket M4-03).
///
/// For cell `i` at `p_i`, with the east and north unit vectors `e_i` and
/// `n_i` of [`coords::east`] and [`coords::north`], each neighbour `j`
/// gives the displacement `d_ij` in the tangent plane: the components
/// `((p_j − p_i)·e_i, (p_j − p_i)·n_i)`, the direction of `p_j` from `p_i`,
/// scaled to the great-circle distance `θ_ij`
/// ([`SphereGrid::centre_distances_rad`]), so that `|d_ij| = θ_ij`. The
/// gradient `g_i` minimises `Σ_j (u_j − u_i − g_i·d_ij)² / θ_ij²`, the
/// weighted least-squares fit of a plane through the cell's value: the
/// normal equations `(Σ_j w_ij d_ij d_ijᵀ) g_i = Σ_j w_ij d_ij (u_j − u_i)`
/// with `w_ij = 1/θ_ij²`, solved once per cell into one coefficient pair
/// per neighbour, so that
///
/// ```text
/// g_i = Σ_j c_ij (u_j − u_i)        (east and north, per radian of arc)
/// ```
///
/// in neighbour order. Divide by the radius for the gradient per metre.
///
/// A constant gives exactly 0. The fit is exact for a field linear in
/// `d`, so the gradient is first-order accurate on the jittered grid of
/// [`SphereGrid::fibonacci`]: for the harmonics `z`, `3z² − 1` and `xy`
/// its error is about 0.5–0.9 % (RMS) and 2–5 % (largest) of the largest
/// `|∇u|` at 10⁴ cells and halves with each fourfold refinement (ticket
/// M4-03, `operators::tests`). On the octahedron, whose neighbours lie a
/// quarter turn away, the gradient of `z` at the equator is `2/π` instead
/// of 1: the chords of a coarse grid.
#[derive(Clone, Debug, PartialEq)]
pub struct Gradient {
    /// CSR offsets: the neighbours of cell `i` are
    /// `offsets[i]..offsets[i + 1]`.
    offsets: Vec<u32>,
    /// The neighbours, as in the grid.
    neighbors: Vec<u32>,
    /// `c_ij` east, per CSR entry (per radian).
    east: Vec<f64>,
    /// `c_ij` north, per CSR entry (per radian).
    north: Vec<f64>,
}

impl Gradient {
    /// The gradient of `grid`: the coefficients of every neighbour pair.
    /// O(n).
    ///
    /// # Panics
    ///
    /// If the neighbours of a cell do not span its tangent plane (not on a
    /// valid grid, where every cell is surrounded by at least three
    /// neighbours).
    pub fn new(grid: &SphereGrid) -> Self {
        let n = grid.len();
        let mut offsets = Vec::with_capacity(n + 1);
        let mut neighbors = Vec::new();
        let mut east = Vec::new();
        let mut north = Vec::new();
        offsets.push(0);
        let mut displacements = Vec::new();
        for i in 0..n as u32 {
            let p = grid.position(i);
            let (e, nv) = (coords::east(p), coords::north(p));
            let distances = grid.centre_distances_rad(i);
            displacements.clear();
            let (mut xx, mut xy, mut yy) = (0.0, 0.0, 0.0);
            for (k, &j) in grid.neighbors(i).iter().enumerate() {
                let q = grid.position(j);
                let d = kuniumi_math::Vec3::new(q.x - p.x, q.y - p.y, q.z - p.z);
                let (x, y) = (d.dot(e), d.dot(nv));
                let theta = distances[k];
                // Scaled to the great-circle distance; `x² + y² > 0` for
                // distinct neighbours off the antipode, as on a valid grid.
                let scale = theta / (x * x + y * y).sqrt();
                let (x, y, w) = (x * scale, y * scale, 1.0 / (theta * theta));
                xx += w * x * x;
                xy += w * x * y;
                yy += w * y * y;
                displacements.push((j, x, y, w));
            }
            let determinant = xx * yy - xy * xy;
            // Written so that NaN fails.
            assert!(
                determinant > 0.0 && determinant.is_finite(),
                "Gradient needs neighbours that surround each cell: cell {i} has a \
                 singular normal matrix"
            );
            for &(j, x, y, w) in &displacements {
                neighbors.push(j);
                east.push(w * (yy * x - xy * y) / determinant);
                north.push(w * (xx * y - xy * x) / determinant);
            }
            // At most 6n − 12 entries (`kuniumi_sphere::MAX_CELLS`), as in
            // the grid's own u32 offsets.
            offsets.push(neighbors.len() as u32);
        }
        Self {
            offsets,
            neighbors,
            east,
            north,
        }
    }

    /// The number of cells.
    pub fn len(&self) -> usize {
        self.offsets.len() - 1
    }

    /// Whether there are no cells.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The gradient of `values` (one per cell, in a unit `[u]`): its east
    /// and north components at every cell, in `[u]` per radian of arc.
    ///
    /// # Panics
    ///
    /// If `values` does not hold one value per cell.
    pub fn apply(&self, values: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let n = self.len();
        assert_eq!(
            values.len(),
            n,
            "apply needs one value per cell: {} values for {n} cells",
            values.len()
        );
        let mut east = Vec::with_capacity(n);
        let mut north = Vec::with_capacity(n);
        for i in 0..n {
            let u = values[i];
            let (mut ge, mut gn) = (0.0, 0.0);
            for k in self.offsets[i] as usize..self.offsets[i + 1] as usize {
                let difference = values[self.neighbors[k] as usize] - u;
                ge += self.east[k] * difference;
                gn += self.north[k] * difference;
            }
            east.push(ge);
            north.push(gn);
        }
        (east, north)
    }
}

#[cfg(test)]
mod tests;

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
//! ([`crate::stages::erosion::diffuse`]), and later the diffusion of heat
//! and of marine sediment.
//!
//! Everything here depends only on its inputs, through IEEE `+ − × ÷` and
//! `sqrt` in a fixed order: every sum and dot product runs sequentially over
//! the cells in index order (docs/conventions.md, "Determinism", rules 2
//! and 6), so the results are the same on every platform.

use kuniumi_sphere::SphereGrid;

/// The relative residual at which [`Laplacian::implicit_step`] stops:
/// `‖b − A u'‖₂ ≤ 10⁻¹² ‖b‖₂`.
pub const SOLVER_TOLERANCE: f64 = 1.0e-12;

/// The number of iterations after which [`Laplacian::implicit_step`] stops
/// whether or not it has reached [`SOLVER_TOLERANCE`]: 10,000. The step
/// needs one or two iterations for the hillslope diffusion of the Erosion
/// stage, about 15 at the explicit limit and a few hundred to a few
/// thousand at 10⁴ to 10¹² times it (10⁴ and 10⁵ cells), 3,600 at
/// 7 × 10⁵ times it on 10⁶ cells (ticket M3-04).
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

/// What one [`Laplacian::implicit_step`] did.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Solve {
    /// The iterations of the conjugate gradient method.
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
    /// Whether the solve stopped because the residual reached
    /// [`SOLVER_TOLERANCE`], or had nothing to solve; `false` where it
    /// stopped on [`SOLVER_MAX_ITERATIONS`] or for lack of a direction.
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

        let mut change = vec![0.0; n];
        let mut z = vec![0.0; n];
        for &i in &cells {
            z[i] = residual[i] / diagonal[i];
        }
        let mut direction = z.clone();
        let mut image = vec![0.0; n];
        let mut rz = dot(&cells, &residual, &z);
        let mut rr = dot(&cells, &residual, &residual);
        let mut iterations = 0;
        let mut converged = rr <= limit;
        while !converged && iterations < SOLVER_MAX_ITERATIONS {
            // image = A' direction; the fixed cells' entries of `direction`
            // are 0.
            for &i in &cells {
                let p = direction[i];
                let mut sum = 0.0;
                for k in self.range(i) {
                    sum += self.weights[k] * (p - direction[self.neighbors[k] as usize]);
                }
                image[i] = mass[i] * p + coupling * sum;
            }
            let pq = dot(&cells, &direction, &image);
            // Written so that NaN stops too.
            if !(pq > 0.0 && pq.is_finite()) {
                break;
            }
            let alpha = rz / pq;
            for &i in &cells {
                change[i] += alpha * direction[i];
                residual[i] -= alpha * image[i];
            }
            iterations += 1;
            rr = dot(&cells, &residual, &residual);
            converged = rr <= limit;
            if converged {
                break;
            }
            for &i in &cells {
                z[i] = residual[i] / diagonal[i];
            }
            let rz_next = dot(&cells, &residual, &z);
            if !(rz_next > 0.0 && rz_next.is_finite()) {
                break;
            }
            let beta = rz_next / rz;
            rz = rz_next;
            for &i in &cells {
                direction[i] = z[i] + beta * direction[i];
            }
        }
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

#[cfg(test)]
mod tests;

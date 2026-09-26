//! Numerical field solver for OpenRDSon.
//!
//! v0 solves the steady current-continuity equation `∇·(σ∇φ) = 0` with the
//! finite-element method on the structured hexahedral mesh, applies Dirichlet
//! terminal voltages and computes the effective resistance from the dissipated
//! power (`R = ΔV² / P`, with `P = uᵀKu`).
//!
//! Element matrices are the standard trilinear 8-node hex stiffness evaluated
//! with 2×2×2 Gauss quadrature. Triangular prisms (used for the few
//! non-rectangular geometry regions) are not yet assembled.

use openrdson_core::mesh::{ElementType, Mesh, Point3};
use std::collections::BTreeMap;

/// Row-wise sparse matrix (`rows[i]` = sorted `(column, value)` pairs).
#[derive(Debug, Clone, Default)]
pub struct Sparse {
    pub n: usize,
    pub rows: Vec<Vec<(u32, f64)>>,
}

impl Sparse {
    pub fn new(n: usize) -> Self {
        Self {
            n,
            rows: vec![Vec::new(); n],
        }
    }

    pub fn matvec(&self, x: &[f64], y: &mut [f64]) {
        for (i, row) in self.rows.iter().enumerate() {
            let mut acc = 0.0;
            for &(j, v) in row {
                acc += v * x[j as usize];
            }
            y[i] = acc;
        }
    }

    pub fn diag(&self, i: usize) -> f64 {
        self.rows[i]
            .iter()
            .find(|(j, _)| *j as usize == i)
            .map(|(_, v)| *v)
            .unwrap_or(0.0)
    }

    fn add(&mut self, i: usize, j: u32, v: f64) {
        self.rows[i].push((j, v));
    }

    /// Add a conductance `g` between nodes `i` and `j` (a channel/device edge).
    /// Call [`Sparse::finalize`] afterwards.
    pub fn add_edge(&mut self, i: u32, j: u32, g: f64) {
        let (i, j) = (i as usize, j as usize);
        self.rows[i].push((j as u32, -g));
        self.rows[j].push((i as u32, -g));
        self.rows[i].push((i as u32, g));
        self.rows[j].push((j as u32, g));
    }

    /// Sum duplicate entries and sort each row (deterministic).
    pub fn finalize(&mut self) {
        for row in &mut self.rows {
            row.sort_by_key(|(j, _)| *j);
            let mut merged: Vec<(u32, f64)> = Vec::with_capacity(row.len());
            for &(j, v) in row.iter() {
                if let Some(last) = merged.last_mut() {
                    if last.0 == j {
                        last.1 += v;
                        continue;
                    }
                }
                merged.push((j, v));
            }
            *row = merged;
        }
    }
}

const SIGNS: [(f64, f64, f64); 8] = [
    (-1.0, -1.0, -1.0),
    (1.0, -1.0, -1.0),
    (1.0, 1.0, -1.0),
    (-1.0, 1.0, -1.0),
    (-1.0, -1.0, 1.0),
    (1.0, -1.0, 1.0),
    (1.0, 1.0, 1.0),
    (-1.0, 1.0, 1.0),
];

fn det3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

fn inv3(m: &[[f64; 3]; 3]) -> Option<[[f64; 3]; 3]> {
    let d = det3(m);
    if d.abs() < 1e-300 {
        return None;
    }
    let id = 1.0 / d;
    Some([
        [
            (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * id,
            (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * id,
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * id,
        ],
        [
            (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * id,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * id,
            (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * id,
        ],
        [
            (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * id,
            (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * id,
            (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * id,
        ],
    ])
}

fn comp(p: Point3, a: usize) -> f64 {
    match a {
        0 => p.x,
        1 => p.y,
        _ => p.z,
    }
}

/// 8×8 trilinear hex stiffness (conductance) matrix for a uniform `sigma`.
fn hex_stiffness(coords: &[Point3; 8], sigma: f64) -> [[f64; 8]; 8] {
    let g = 1.0 / 3.0f64.sqrt();
    let gp = [-g, g];
    let mut ke = [[0.0f64; 8]; 8];

    for &xi in &gp {
        for &eta in &gp {
            for &zeta in &gp {
                let mut dn = [[0.0f64; 3]; 8];
                for i in 0..8 {
                    let (si, ei, zi) = SIGNS[i];
                    dn[i][0] = 0.125 * si * (1.0 + ei * eta) * (1.0 + zi * zeta);
                    dn[i][1] = 0.125 * ei * (1.0 + si * xi) * (1.0 + zi * zeta);
                    dn[i][2] = 0.125 * zi * (1.0 + si * xi) * (1.0 + ei * eta);
                }
                let mut j = [[0.0f64; 3]; 3];
                for (i, c) in coords.iter().enumerate() {
                    for a in 0..3 {
                        for b in 0..3 {
                            j[a][b] += comp(*c, a) * dn[i][b];
                        }
                    }
                }
                let detj = det3(&j);
                if detj <= 0.0 {
                    continue;
                }
                let Some(inv) = inv3(&j) else { continue };

                // Physical gradients: grad[i][a] = sum_b inv[b][a] * dn[i][b].
                let mut grad = [[0.0f64; 3]; 8];
                for i in 0..8 {
                    for a in 0..3 {
                        let mut s = 0.0;
                        for b in 0..3 {
                            s += inv[b][a] * dn[i][b];
                        }
                        grad[i][a] = s;
                    }
                }
                for i in 0..8 {
                    for k in 0..8 {
                        let dot = grad[i][0] * grad[k][0]
                            + grad[i][1] * grad[k][1]
                            + grad[i][2] * grad[k][2];
                        ke[i][k] += sigma * dot * detj;
                    }
                }
            }
        }
    }
    ke
}

/// Assemble the global conductance matrix for hex elements with uniform `sigma`.
pub fn assemble_conductance(mesh: &Mesh, sigma: f64) -> (Sparse, Vec<String>) {
    assemble_conductance_by(mesh, &|_| sigma)
}

/// Assemble the global conductance matrix, looking up conductivity per element.
///
/// `sigma_of` receives the element's `material_id` (0 = unassigned) and returns
/// its conductivity in S/m. Returns the matrix plus diagnostics for skipped
/// elements.
pub fn assemble_conductance_by(
    mesh: &Mesh,
    sigma_of: &dyn Fn(u32) -> f64,
) -> (Sparse, Vec<String>) {
    let n = mesh.nodes.len();
    let mut a = Sparse::new(n);
    let mut diagnostics = Vec::new();
    let mut skipped = 0usize;

    for e in mesh.elements.iter() {
        if e.kind != ElementType::Hex || e.nodes.len() != 8 {
            skipped += 1;
            continue;
        }
        let mut coords = [Point3::new(0.0, 0.0, 0.0); 8];
        for (k, &ni) in e.nodes.iter().enumerate() {
            coords[k] = mesh.nodes[ni as usize];
        }
        let sigma = sigma_of(e.material_id);
        if sigma <= 0.0 {
            continue; // insulating region
        }
        let ke = hex_stiffness(&coords, sigma);
        for i in 0..8 {
            for k in 0..8 {
                a.add(e.nodes[i] as usize, e.nodes[k], ke[i][k]);
            }
        }
    }
    if skipped > 0 {
        diagnostics.push(format!(
            "{skipped} non-hex elements skipped during assembly (prism support pending)"
        ));
    }
    openrdson_core::log_debug!(
        "assembled conductance: {} nodes, {} elements ({} skipped)",
        n,
        mesh.elements.len().saturating_sub(skipped),
        skipped
    );
    a.finalize();
    (a, diagnostics)
}

/// Incomplete Cholesky IC(0) preconditioner: `L Lᵀ ≈ A` with `L` sharing A's
/// lower-triangular sparsity. Much stronger than Jacobi on 2D/3D Laplacians.
pub struct Ico {
    /// Rows of L: `(column, value)` with `column <= row`.
    l: Vec<Vec<(u32, f64)>>,
    /// Columns of L: `col[j] = [(row > j, value)]`.
    col: Vec<Vec<(u32, f64)>>,
    diag: Vec<f64>,
}

impl Ico {
    pub fn new(a: &Sparse) -> Option<Self> {
        let n = a.n;
        let mut l: Vec<Vec<(u32, f64)>> = vec![Vec::new(); n];
        let mut col: Vec<Vec<(u32, f64)>> = vec![Vec::new(); n];
        let mut diag = vec![0.0f64; n];
        for i in 0..n {
            for &(jc, aij) in &a.rows[i] {
                let j = jc as usize;
                if j > i {
                    continue;
                }
                let mut s = aij;
                let (mut p, mut q) = (0usize, 0usize);
                while p < l[i].len() && q < l[j].len() {
                    let (ci, vi) = l[i][p];
                    let (cj, vj) = l[j][q];
                    if ci == cj {
                        s -= vi * vj;
                        p += 1;
                        q += 1;
                    } else if ci < cj {
                        p += 1;
                    } else {
                        q += 1;
                    }
                }
                if i == j {
                    if s <= 0.0 {
                        return None;
                    }
                    let lii = s.sqrt();
                    diag[i] = lii;
                    l[i].push((i as u32, lii));
                } else {
                    let ljj = diag[j];
                    if ljj == 0.0 {
                        return None;
                    }
                    let lij = s / ljj;
                    l[i].push((j as u32, lij));
                    col[j].push((i as u32, lij));
                }
            }
        }
        Some(Ico { l, col, diag })
    }

    /// Apply `z = M⁻¹ r` with `M = L Lᵀ`.
    pub fn apply(&self, r: &[f64], z: &mut [f64]) {
        let n = r.len();
        let mut y = vec![0.0f64; n];
        for i in 0..n {
            let mut s = r[i];
            for &(j, v) in &self.l[i] {
                if (j as usize) < i {
                    s -= v * y[j as usize];
                }
            }
            y[i] = s / self.diag[i];
        }
        for i in (0..n).rev() {
            let mut s = y[i];
            for &(k, v) in &self.col[i] {
                s -= v * z[k as usize];
            }
            z[i] = s / self.diag[i];
        }
    }
}

/// IC(0)-preconditioned conjugate gradient.
pub fn pcg_ic0(
    a: &Sparse,
    pre: &Ico,
    b: &[f64],
    x: &mut [f64],
    tol: f64,
    max_iter: usize,
) -> Result<usize, String> {
    let n = a.n;
    if n == 0 {
        return Ok(0);
    }
    let mut r = vec![0.0f64; n];
    let mut ax = vec![0.0f64; n];
    a.matvec(x, &mut ax);
    for i in 0..n {
        r[i] = b[i] - ax[i];
    }
    let mut z = vec![0.0f64; n];
    let mut p = vec![0.0f64; n];
    let mut ap = vec![0.0f64; n];
    pre.apply(&r, &mut z);
    p.copy_from_slice(&z);
    let mut rz: f64 = r.iter().zip(&z).map(|(a, b)| a * b).sum();
    let bnorm = b.iter().map(|v| v * v).sum::<f64>().sqrt().max(1e-300);

    for it in 0..max_iter {
        let rnorm = r.iter().map(|v| v * v).sum::<f64>().sqrt();
        if rnorm / bnorm < tol {
            return Ok(it);
        }
        a.matvec(&p, &mut ap);
        let pap: f64 = p.iter().zip(&ap).map(|(a, b)| a * b).sum();
        if pap.abs() < 1e-300 {
            return Err("IC(0)-CG breakdown".into());
        }
        let alpha = rz / pap;
        for i in 0..n {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        pre.apply(&r, &mut z);
        let rz_new: f64 = r.iter().zip(&z).map(|(a, b)| a * b).sum();
        let beta = rz_new / rz.max(1e-300);
        for i in 0..n {
            p[i] = z[i] + beta * p[i];
        }
        rz = rz_new;
    }
    Err(format!("IC(0)-CG did not converge in {max_iter} iterations"))
}

/// Solve `K u = 0` with Dirichlet fixed values using Jacobi-preconditioned CG.
pub fn solve_dirichlet(
    a: &Sparse,
    fixed: &BTreeMap<u32, f64>,
    tol: f64,
    max_iter: usize,
) -> Result<Vec<f64>, String> {
    let x0 = vec![0.0f64; a.n];
    solve_dirichlet_warm(a, fixed, &x0, tol, max_iter)
}

/// Like [`solve_dirichlet`], but seeds the CG with the initial guess `x0`
/// (length `a.n`). Fixed nodes are still forced to their Dirichlet values; only
/// the free-node block of `x0` is used as the CG starting point.
pub fn solve_dirichlet_warm(
    a: &Sparse,
    fixed: &BTreeMap<u32, f64>,
    x0: &[f64],
    tol: f64,
    max_iter: usize,
) -> Result<Vec<f64>, String> {
    let n = a.n;
    let mut u = vec![0.0f64; n];
    for (&i, &v) in fixed {
        u[i as usize] = v;
    }

    // Map free nodes to reduced indices.
    let mut free_index = vec![usize::MAX; n];
    let mut free_nodes = Vec::new();
    for i in 0..n {
        if !fixed.contains_key(&(i as u32)) {
            free_index[i] = free_nodes.len();
            free_nodes.push(i);
        }
    }
    let m = free_nodes.len();
    if m == 0 {
        return Ok(u);
    }

    // Reduced matrix and rhs: K_ff u_f = -K_fc u_c.
    let mut red = Sparse::new(m);
    let mut b = vec![0.0f64; m];
    for (ri, &i) in free_nodes.iter().enumerate() {
        for &(j, v) in &a.rows[i] {
            let j = j as usize;
            if fixed.contains_key(&(j as u32)) {
                b[ri] -= v * u[j];
            } else {
                red.add(ri, free_index[j] as u32, v);
            }
        }
    }
    red.finalize();

    let mut x = vec![0.0f64; m];
    for (ri, &i) in free_nodes.iter().enumerate() {
        x[ri] = x0[i];
    }
    let iters = if let Some(pre) = Ico::new(&red) {
        pcg_ic0(&red, &pre, &b, &mut x, tol, max_iter)?
    } else {
        pcg(&red, &b, &mut x, tol, max_iter)?
    };
    openrdson_core::log_debug!(
        "linear solve: {m} free DOF converged in {iters} iterations"
    );
    for (ri, &i) in free_nodes.iter().enumerate() {
        u[i] = x[ri];
    }
    Ok(u)
}

/// Direct Dirichlet solve via `faer`'s sparse Cholesky. Intended for the **2D
/// sheet network**, where a direct factorization is far faster than PCG (fill
/// is O(n log n) in 2D); it falls back to [`solve_dirichlet`] (IC(0)-PCG) if the
/// factorization fails or the feature is disabled.
#[cfg(feature = "faer-sheet")]
pub fn solve_dirichlet_direct(
    a: &Sparse,
    fixed: &BTreeMap<u32, f64>,
    tol: f64,
    max_iter: usize,
) -> Result<Vec<f64>, String> {
    use faer::prelude::*;
    use faer::sparse::{SparseColMat, Triplet};
    use faer::Side;

    let n = a.n;
    let mut u = vec![0.0f64; n];
    for (&i, &v) in fixed {
        u[i as usize] = v;
    }
    let mut free_index = vec![usize::MAX; n];
    let mut free_nodes = Vec::new();
    for i in 0..n {
        if !fixed.contains_key(&(i as u32)) {
            free_index[i] = free_nodes.len();
            free_nodes.push(i);
        }
    }
    let m = free_nodes.len();
    if m == 0 {
        return Ok(u);
    }
    let mut trips: Vec<Triplet<usize, usize, f64>> = Vec::new();
    let mut b = vec![0.0f64; m];
    for (ri, &i) in free_nodes.iter().enumerate() {
        for &(j, v) in &a.rows[i] {
            let j = j as usize;
            if fixed.contains_key(&(j as u32)) {
                b[ri] -= v * u[j];
            } else {
                trips.push(Triplet::new(ri, free_index[j], v));
            }
        }
    }
    let Ok(mat) = SparseColMat::<usize, f64>::try_new_from_triplets(m, m, &trips) else {
        return solve_dirichlet(a, fixed, tol, max_iter);
    };
    let Ok(llt) = mat.sp_cholesky(Side::Lower) else {
        return solve_dirichlet(a, fixed, tol, max_iter);
    };
    let bm = faer::Mat::from_fn(m, 1, |i, _| b[i]);
    let x = llt.solve(&bm);
    for (ri, &i) in free_nodes.iter().enumerate() {
        u[i] = x[(ri, 0)];
    }
    Ok(u)
}

/// Jacobi-preconditioned conjugate gradient.
pub fn pcg(
    a: &Sparse,
    b: &[f64],
    x: &mut [f64],
    tol: f64,
    max_iter: usize,
) -> Result<usize, String> {
    let n = a.n;
    if n == 0 {
        return Ok(0);
    }
    let mut r = vec![0.0f64; n];
    let mut ax = vec![0.0f64; n];
    a.matvec(x, &mut ax);
    for i in 0..n {
        r[i] = b[i] - ax[i];
    }
    let mut z = vec![0.0f64; n];
    let mut p = vec![0.0f64; n];
    let mut ap = vec![0.0f64; n];
    let mut rz = 0.0f64;
    for i in 0..n {
        let d = a.diag(i);
        let minv = if d.abs() > 1e-300 { 1.0 / d } else { 1.0 };
        z[i] = minv * r[i];
        p[i] = z[i];
        rz += r[i] * z[i];
    }
    let bnorm = b.iter().map(|v| v * v).sum::<f64>().sqrt().max(1e-300);

    for it in 0..max_iter {
        let rnorm = r.iter().map(|v| v * v).sum::<f64>().sqrt();
        if rnorm / bnorm < tol {
            return Ok(it);
        }
        a.matvec(&p, &mut ap);
        let pap: f64 = p.iter().zip(&ap).map(|(a, b)| a * b).sum();
        if pap.abs() < 1e-300 {
            return Err("CG breakdown (pᵀAp ≈ 0)".into());
        }
        let alpha = rz / pap;
        for i in 0..n {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        let mut rz_new = 0.0;
        for i in 0..n {
            let d = a.diag(i);
            let minv = if d.abs() > 1e-300 { 1.0 / d } else { 1.0 };
            z[i] = minv * r[i];
            rz_new += r[i] * z[i];
        }
        let beta = rz_new / rz.max(1e-300);
        for i in 0..n {
            p[i] = z[i] + beta * p[i];
        }
        rz = rz_new;
    }
    Err(format!("CG did not converge in {max_iter} iterations"))
}

/// Total dissipated power `P = uᵀKu` for a solution `u` (watts, for σ in S/m).
pub fn total_power(a: &Sparse, u: &[f64]) -> f64 {
    let mut ku = vec![0.0f64; a.n];
    a.matvec(u, &mut ku);
    u.iter().zip(&ku).map(|(a, b)| a * b).sum()
}

/// Effective resistance from a two-terminal solve with `ΔV = 1`: `R = 1 / P`.
pub fn effective_resistance_from_power(a: &Sparse, u: &[f64]) -> f64 {
    let p = total_power(a, u);
    if p.abs() < 1e-300 {
        f64::INFINITY
    } else {
        1.0 / p
    }
}

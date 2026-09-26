//! Head-to-head benchmark: the in-house IC(0)-PCG solver vs `faer`'s sparse
//! Cholesky, on structured 3D hex-grid Laplacians (the same SPD form as the
//! FEM conductance matrix).
//!
//! Run: `cargo run --release --example solve_bench`

use std::collections::BTreeMap;
use std::time::Instant;

use openrdson_core::mesh::{Element, ElementType, Mesh, Point3};
use openrdson_solver::{assemble_conductance_by, solve_dirichlet, Sparse};

fn grid(n: usize) -> Mesh {
    let idx = |i: usize, j: usize, k: usize| (i + n * (j + n * k)) as u32;
    let mut nodes = Vec::with_capacity(n * n * n);
    for k in 0..n {
        for j in 0..n {
            for i in 0..n {
                nodes.push(Point3::new(i as f64, j as f64, k as f64));
            }
        }
    }
    let mut elements = Vec::new();
    for k in 0..n - 1 {
        for j in 0..n - 1 {
            for i in 0..n - 1 {
                elements.push(Element {
                    kind: ElementType::Hex,
                    nodes: vec![
                        idx(i, j, k),
                        idx(i + 1, j, k),
                        idx(i + 1, j + 1, k),
                        idx(i, j + 1, k),
                        idx(i, j, k + 1),
                        idx(i + 1, j, k + 1),
                        idx(i + 1, j + 1, k + 1),
                        idx(i, j + 1, k + 1),
                    ],
                    material_id: 1,
                    net: None,
                    layer: None,
                    device_ref: None,
                    source_polygon: None,
                });
            }
        }
    }
    Mesh {
        nodes,
        elements,
        sets: Default::default(),
    }
}

/// faer sparse-Cholesky solve of the same Dirichlet system.
fn solve_faer(a: &Sparse, fixed: &BTreeMap<u32, f64>) -> Vec<f64> {
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
    let mut trips = Vec::new();
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
    let mat = SparseColMat::<usize, f64>::try_new_from_triplets(m, m, &trips).unwrap();
    let llt = mat.sp_cholesky(Side::Lower).unwrap();
    let bm = faer::Mat::from_fn(m, 1, |i, _| b[i]);
    let x = llt.solve(&bm);
    for (ri, &i) in free_nodes.iter().enumerate() {
        u[i] = x[(ri, 0)];
    }
    u
}

fn run(n: usize) {
    let mesh = grid(n);
    let (a, _) = assemble_conductance_by(&mesh, &|_| 1.0);

    // Dirichlet: z=0 face at 1 V, z=(n-1) face at 0 V.
    let mut fixed = BTreeMap::new();
    for k in [0usize, n - 1] {
        for j in 0..n {
            for i in 0..n {
                let id = (i + n * (j + n * k)) as u32;
                fixed.insert(id, if k == 0 { 1.0 } else { 0.0 });
            }
        }
    }

    let t = Instant::now();
    let u_home = solve_dirichlet(&a, &fixed, 1e-10, 20000).unwrap();
    let dt_home = t.elapsed();

    let t = Instant::now();
    let u_faer = solve_faer(&a, &fixed);
    let dt_faer = t.elapsed();

    let diff = u_home
        .iter()
        .zip(&u_faer)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0f64, f64::max);

    println!(
        "n={n:3}  dof={:8}  home(IC0-PCG)={:>10.3?}  faer(Cholesky)={:>10.3?}  speedup={:>6.1}x  maxdiff={:.2e}",
        mesh.nodes.len(),
        dt_home,
        dt_faer,
        dt_home.as_secs_f64() / dt_faer.as_secs_f64().max(1e-9),
        diff
    );
}

/// 2D 5-point Laplacian (the sheet network's form), built directly.
fn grid2d(n: usize) -> Sparse {
    let idx = |i: usize, j: usize| (i + n * j) as u32;
    let mut a = Sparse::new(n * n);
    for j in 0..n {
        for i in 0..n {
            let id = idx(i, j);
            if i + 1 < n {
                a.add_edge(id, idx(i + 1, j), 1.0);
            }
            if j + 1 < n {
                a.add_edge(id, idx(i, j + 1), 1.0);
            }
        }
    }
    a.finalize();
    a
}

fn run2d(n: usize) {
    let a = grid2d(n);
    let mut fixed = BTreeMap::new();
    for j in 0..n {
        for i in 0..n {
            let id = (i + n * j) as u32;
            if i == 0 || i == n - 1 || j == 0 || j == n - 1 {
                fixed.insert(id, if i == 0 { 1.0 } else { 0.0 });
            }
        }
    }

    let t = Instant::now();
    let u_home = solve_dirichlet(&a, &fixed, 1e-10, 20000).unwrap();
    let dt_home = t.elapsed();

    let t = Instant::now();
    let u_faer = solve_faer(&a, &fixed);
    let dt_faer = t.elapsed();

    let diff = u_home
        .iter()
        .zip(&u_faer)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0f64, f64::max);
    println!(
        "2D n={n:3} dof={:8}  home(IC0-PCG)={:>10.3?}  faer(Cholesky)={:>10.3?}  speedup={:>6.1}x  maxdiff={:.2e}",
        n * n,
        dt_home,
        dt_faer,
        dt_home.as_secs_f64() / dt_faer.as_secs_f64().max(1e-9),
        diff
    );
}

fn main() {
    println!("in-house IC(0)-PCG  vs  faer sparse Cholesky");
    println!("-- 3D hex-grid Laplacian (full-stack FEM form) --");
    for n in [20usize, 30, 40, 50] {
        run(n);
    }
    println!("-- 2D 5-point Laplacian (sheet-network form) --");
    for n in [100usize, 200, 300, 500, 800] {
        run2d(n);
    }
}

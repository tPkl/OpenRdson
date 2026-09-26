//! Analytic acceptance test: a uniform bar's FEM resistance must equal `ρL/A`
//! (equivalently `L/(σA)`), per the solver goal acceptance criterion #1.

use openrdson_core::geometry::Point;
use openrdson_geometry::extrude_polygon;
use openrdson_meshing::{mesh_region, MeshConfig};
use openrdson_solver::{assemble_conductance, effective_resistance_from_power, solve_dirichlet};
use std::collections::BTreeMap;

fn bar(l: f64, w: f64, t: f64) -> openrdson_core::geometry::SolidRegion {
    let pts = vec![
        Point::new(0.0, 0.0),
        Point::new(l, 0.0),
        Point::new(l, w),
        Point::new(0.0, w),
    ];
    extrude_polygon(&pts, 0.0, t, 1, None, Some("Metal1".into()), None)
}

#[test]
fn fem_resistance_matches_rho_l_over_a() {
    let (l, w, t) = (10e-6, 2e-6, 0.325e-6);
    let sigma = 5.8e7; // S/m (copper)
    let region = bar(l, w, t);
    let mesh = mesh_region(
        &region,
        &MeshConfig {
            lateral_cell: 1e-6,
            z_cell: 0.1e-6,
        },
    )
    .unwrap();

    let (a, diags) = assemble_conductance(&mesh, sigma);
    assert!(diags.is_empty(), "unexpected assembly diagnostics: {diags:?}");

    // Dirichlet terminals: x = 0 face at 1 V, x = L face at 0 V.
    let mut fixed = BTreeMap::new();
    let eps = 1e-12;
    for (i, p) in mesh.nodes.iter().enumerate() {
        if p.x.abs() < eps {
            fixed.insert(i as u32, 1.0);
        } else if (p.x - l).abs() < eps {
            fixed.insert(i as u32, 0.0);
        }
    }
    assert!(!fixed.is_empty());

    let u = solve_dirichlet(&a, &fixed, 1e-12, 5000).unwrap();
    let r = effective_resistance_from_power(&a, &u);
    let analytic = l / (sigma * w * t);

    assert!(
        (r - analytic).abs() / analytic < 1e-6,
        "FEM R = {r:e} Ω vs analytic {analytic:e} Ω"
    );
}

#[test]
fn resistance_scales_linearly_with_length() {
    let sigma = 1e6;
    let w = 1e-6;
    let t = 0.2e-6;
    let mut r = [0.0; 2];
    for (k, l) in [2e-6f64, 4e-6].iter().enumerate() {
        let region = bar(*l, w, t);
        let mesh = mesh_region(
            &region,
            &MeshConfig {
                lateral_cell: 0.5e-6,
                z_cell: 0.1e-6,
            },
        )
        .unwrap();
        let (a, _) = assemble_conductance(&mesh, sigma);
        let mut fixed = BTreeMap::new();
        for (i, p) in mesh.nodes.iter().enumerate() {
            if p.x.abs() < 1e-12 {
                fixed.insert(i as u32, 1.0);
            } else if (p.x - *l).abs() < 1e-12 {
                fixed.insert(i as u32, 0.0);
            }
        }
        let u = solve_dirichlet(&a, &fixed, 1e-12, 5000).unwrap();
        r[k] = effective_resistance_from_power(&a, &u);
    }
    // Doubling length doubles resistance.
    assert!((r[1] / r[0] - 2.0).abs() < 1e-6, "R ratio {}", r[1] / r[0]);
}

#[test]
fn resistance_scales_inversely_with_cross_section() {
    let sigma = 1e6;
    let l = 4e-6;
    let t = 0.2e-6;
    let mut r = [0.0; 2];
    for (k, w) in [1e-6f64, 2e-6].iter().enumerate() {
        let region = bar(l, *w, t);
        let mesh = mesh_region(
            &region,
            &MeshConfig {
                lateral_cell: 0.5e-6,
                z_cell: 0.1e-6,
            },
        )
        .unwrap();
        let (a, _) = assemble_conductance(&mesh, sigma);
        let mut fixed = BTreeMap::new();
        for (i, p) in mesh.nodes.iter().enumerate() {
            if p.x.abs() < 1e-12 {
                fixed.insert(i as u32, 1.0);
            } else if (p.x - l).abs() < 1e-12 {
                fixed.insert(i as u32, 0.0);
            }
        }
        let u = solve_dirichlet(&a, &fixed, 1e-12, 5000).unwrap();
        r[k] = effective_resistance_from_power(&a, &u);
    }
    // Doubling cross-section halves resistance.
    assert!((r[1] / r[0] - 0.5).abs() < 1e-6, "R ratio {}", r[1] / r[0]);
}

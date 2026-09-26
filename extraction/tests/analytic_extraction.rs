//! Analytic acceptance tests for resistive network extraction.

use openrdson_core::geometry::{Point, SolidModel, SolidRegion};
use openrdson_extraction::{terminal_network, two_terminal_resistance};
use openrdson_geometry::extrude_polygon;
use openrdson_meshing::{mesh_model, weld_nodes, MeshConfig};

fn bar_region(x0: f64, x1: f64, w: f64, t: f64, mat: u32) -> SolidRegion {
    let pts = vec![
        Point::new(x0, 0.0),
        Point::new(x1, 0.0),
        Point::new(x1, w),
        Point::new(x0, w),
    ];
    extrude_polygon(&pts, 0.0, t, mat, None, None, None)
}

fn terminals_at_x(mesh: &openrdson_core::mesh::Mesh, x: f64) -> Vec<u32> {
    mesh.nodes
        .iter()
        .enumerate()
        .filter(|(_, p)| (p.x - x).abs() < 1e-12)
        .map(|(i, _)| i as u32)
        .collect()
}

fn cfg() -> MeshConfig {
    MeshConfig {
        lateral_cell: 0.5e-6,
        z_cell: 0.1e-6,
    }
}

#[test]
fn single_bar_resistance_matches_analytic() {
    let (l, w, t) = (4e-6, 1e-6, 0.2e-6);
    let sigma = 1e6;
    let model = SolidModel {
        regions: vec![bar_region(0.0, l, w, t, 1)],
        material_names: vec!["A".into()],
    };
    let (mesh, _) = mesh_model(&model, &cfg());
    let mesh = weld_nodes(&mesh, 1e-12);

    let a = terminals_at_x(&mesh, 0.0);
    let b = terminals_at_x(&mesh, l);
    let r = two_terminal_resistance(&mesh, &|_| sigma, &a, &b, 1e-12, 5000).unwrap();
    let analytic = l / (sigma * w * t);
    assert!((r - analytic).abs() / analytic < 1e-6, "R={r:e} vs {analytic:e}");
}

#[test]
fn two_materials_in_series_add_resistance() {
    let (l1, l2) = (4e-6, 4e-6);
    let (w, t) = (1e-6, 0.2e-6);
    let (sigma1, sigma2) = (1e6, 2e6);
    let model = SolidModel {
        regions: vec![
            bar_region(0.0, l1, w, t, 1),
            bar_region(l1, l1 + l2, w, t, 2),
        ],
        material_names: vec!["A".into(), "B".into()],
    };
    let (mesh, diags) = mesh_model(&model, &cfg());
    assert!(diags.is_empty());
    // Weld the shared interface at x = l1 so the two materials are in series.
    let mesh = weld_nodes(&mesh, 1e-12);

    let a = terminals_at_x(&mesh, 0.0);
    let b = terminals_at_x(&mesh, l1 + l2);
    let r = two_terminal_resistance(
        &mesh,
        &|id| match id {
            1 => sigma1,
            2 => sigma2,
            _ => 0.0,
        },
        &a,
        &b,
        1e-12,
        5000,
    )
    .unwrap();

    let area = w * t;
    let analytic = l1 / (sigma1 * area) + l2 / (sigma2 * area);
    assert!(
        (r - analytic).abs() / analytic < 1e-6,
        "series R={r:e} vs analytic {analytic:e}"
    );
}

#[test]
fn terminal_network_is_symmetric_with_zero_diagonal() {
    let (l, w, t) = (4e-6, 1e-6, 0.2e-6);
    let model = SolidModel {
        regions: vec![bar_region(0.0, l, w, t, 1)],
        material_names: vec!["A".into()],
    };
    let (mesh, _) = mesh_model(&model, &cfg());
    let mesh = weld_nodes(&mesh, 1e-12);

    let terminals = vec![
        ("D".to_string(), terminals_at_x(&mesh, 0.0)),
        ("S".to_string(), terminals_at_x(&mesh, l)),
    ];
    let net = terminal_network(&mesh, &|_| 1e6, &terminals, 1e-12, 5000).unwrap();
    assert_eq!(net.terminals, vec!["D", "S"]);
    assert_eq!(net.rmatrix[0][0], 0.0);
    assert_eq!(net.rmatrix[1][1], 0.0);
    assert!(net.rmatrix[0][1] > 0.0);
    assert!((net.rmatrix[0][1] - net.rmatrix[1][0]).abs() < 1e-18);
    assert_eq!(net.resistance("D", "S"), Some(net.rmatrix[0][1]));
}

//! Analytic acceptance test: a rectangular bar's mesh volume must match the
//! exact geometric volume (meshing `goal.md` M2 / acceptance criterion #2).

use openrdson_core::geometry::{Point, Solid};
use openrdson_geometry::extrude_polygon;
use openrdson_meshing::{mesh_region, mesh_volume, MeshConfig};

fn bar_footprint(l: f64, w: f64) -> Vec<Point> {
    vec![
        Point::new(0.0, 0.0),
        Point::new(l, 0.0),
        Point::new(l, w),
        Point::new(0.0, w),
    ]
}

#[test]
fn rectangle_extrudes_to_a_box() {
    let pts = bar_footprint(10e-6, 2e-6);
    let region = extrude_polygon(
        &pts,
        0.0,
        0.325e-6,
        1,
        Some("net_M1".into()),
        Some("Metal1".into()),
        None,
    );
    assert!(matches!(region.solid, Solid::Box { .. }));
    // 10 um * 2 um * 0.325 um = 6.5e-18 m^3
    assert!((region.volume() - 6.5e-18).abs() / 6.5e-18 < 1e-12);
}

#[test]
fn mesh_volume_matches_analytic_bar_within_0_1_percent() {
    let pts = bar_footprint(10e-6, 2e-6);
    let region = extrude_polygon(&pts, 0.0, 0.325e-6, 1, None, None, None);
    let analytic = region.volume();

    let cfg = MeshConfig {
        lateral_cell: 1e-6,
        z_cell: 0.1e-6,
    };
    let mesh = mesh_region(&region, &cfg).expect("box meshes");

    // 10 x 2 lateral cells, 4 z layers (ceil(0.325/0.1)).
    assert_eq!(mesh.elements.len(), 80);
    assert_eq!(mesh.nodes.len(), 11 * 3 * 5);

    let vol = mesh_volume(&mesh);
    assert!(
        (vol - analytic).abs() / analytic < 1e-3,
        "mesh volume {vol:e} vs analytic {analytic:e}"
    );
}

#[test]
fn mesh_is_tagged_with_material_and_net() {
    let pts = bar_footprint(2e-6, 2e-6);
    let region = extrude_polygon(
        &pts,
        0.0,
        0.1e-6,
        7,
        Some("net_M2".into()),
        Some("Metal2".into()),
        Some("X0".into()),
    );
    let mesh = mesh_region(&region, &MeshConfig::default()).unwrap();
    let e = &mesh.elements[0];
    assert_eq!(e.material_id, 7);
    assert_eq!(e.net.as_deref(), Some("net_M2"));
    assert_eq!(e.layer.as_deref(), Some("Metal2"));
    assert_eq!(e.device_ref.as_deref(), Some("X0"));
}

#[test]
fn non_rectangular_footprint_meshes_as_triangular_prism() {
    // A triangle is not an axis-aligned rectangle.
    let pts = vec![
        Point::new(0.0, 0.0),
        Point::new(1e-6, 0.0),
        Point::new(0.0, 1e-6),
    ];
    let region = extrude_polygon(&pts, 0.0, 1e-6, 1, None, None, None);
    assert!(matches!(region.solid, Solid::Prism { .. }));

    let mesh = mesh_region(&region, &MeshConfig::default()).expect("prism meshes");
    assert_eq!(mesh.elements.len(), 1);
    assert_eq!(mesh.elements[0].kind, openrdson_core::mesh::ElementType::Prism);
    assert_eq!(mesh.nodes.len(), 6);
    // Area 0.5e-12 m^2 * height 1e-6 m = 0.5e-18 m^3.
    assert!((mesh_volume(&mesh) - 0.5e-18).abs() / 0.5e-18 < 1e-9);
}

#[test]
fn grid_sizes_round_up_to_cover_the_region() {
    let pts = bar_footprint(2.5e-6, 0.5e-6);
    let region = extrude_polygon(&pts, 0.0, 0.25e-6, 1, None, None, None);
    let cfg = MeshConfig {
        lateral_cell: 1e-6,
        z_cell: 0.1e-6,
    };
    let mesh = mesh_region(&region, &cfg).unwrap();
    // nx=3, ny=1, nz=3 -> 9 elements; nodes 4*2*4 = 32.
    assert_eq!(mesh.elements.len(), 9);
    assert_eq!(mesh.nodes.len(), 32);
}

#[test]
fn layer_refinement_overrides_cell_size() {
    use openrdson_core::geometry::SolidModel;
    use openrdson_meshing::{mesh_model_with_options, MeshOptions};
    let a = extrude_polygon(&bar_footprint(4.0, 2.0), 0.0, 1.0, 1, None, Some("Metal1".into()), None);
    let b = extrude_polygon(&bar_footprint(4.0, 2.0), 0.0, 1.0, 2, None, Some("Poly1".into()), None);
    let model = SolidModel {
        regions: vec![a, b],
        material_names: vec!["Metal1".into(), "Poly1".into()],
    };
    let opts = MeshOptions {
        default: MeshConfig {
            lateral_cell: 2.0,
            z_cell: 1.0,
        },
        per_layer: Default::default(),
    }
    .with_layer("Poly1", 0.5, 1.0);
    let (mesh, diags) = mesh_model_with_options(&model, &opts);
    assert!(diags.is_empty());
    // Metal1 (2 um cells): 2x1 = 2 hexes; Poly1 (0.5 um cells): 8x4 = 32 hexes.
    let metal1 = mesh.elements.iter().filter(|e| e.material_id == 1).count();
    let poly1 = mesh.elements.iter().filter(|e| e.material_id == 2).count();
    assert_eq!(metal1, 2);
    assert_eq!(poly1, 32);
}

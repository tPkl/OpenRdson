//! The coarse-mesh short fix: nodes shared across connectivity components must
//! be split so distinct nets never share a node.

use openrdson_core::geometry::{Point, SolidModel};
use openrdson_extraction::split_nodes_by_component;
use openrdson_geometry::extrude_polygon;
use openrdson_meshing::{mesh_model_global, MeshConfig};
use std::collections::HashSet;

fn box_region(x0: f64, y0: f64, x1: f64, y1: f64, source: usize) -> openrdson_core::geometry::SolidRegion {
    let pts = vec![
        Point::new(x0, y0),
        Point::new(x1, y0),
        Point::new(x1, y1),
        Point::new(x0, y1),
    ];
    let mut r = extrude_polygon(&pts, 0.0, 1.0, 1, None, None, None);
    r.source_polygon = Some(source);
    r
}

#[test]
fn split_isolates_two_nets_that_share_an_interface() {
    // Two abutting boxes belonging to different nets (components 0 and 1).
    let model = SolidModel {
        regions: vec![
            box_region(0.0, 0.0, 1.0, 1.0, 0),
            box_region(1.0, 0.0, 2.0, 1.0, 1),
        ],
        material_names: vec!["A".into()],
    };
    let cfg = MeshConfig {
        lateral_cell: 0.5,
        z_cell: 0.5,
    };
    let (shared, _) = mesh_model_global(&model, &cfg);
    // The global grid shares the interface nodes (conforming)...
    assert_eq!(shared.nodes.len(), 45);

    let component_of = vec![0usize, 1];
    let split = split_nodes_by_component(&shared, &component_of);

    // ...and after splitting, no node is used by both components.
    let mut nodes0: HashSet<u32> = HashSet::new();
    let mut nodes1: HashSet<u32> = HashSet::new();
    for e in &split.elements {
        let c = e.source_polygon.unwrap();
        for &n in &e.nodes {
            if c == 0 {
                nodes0.insert(n);
            } else {
                nodes1.insert(n);
            }
        }
    }
    assert!(
        nodes0.intersection(&nodes1).next().is_none(),
        "nets still share nodes after split"
    );
    assert!(split.nodes.len() > shared.nodes.len());
}

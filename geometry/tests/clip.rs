//! Tests for clipping a solid model to a window.

use openrdson_core::geometry::{Bbox, Point, Solid, SolidModel};
use openrdson_geometry::{clip_model, extrude_polygon};

fn box_region(x0: f64, y0: f64, x1: f64, y1: f64, z0: f64, z1: f64) -> SolidModel {
    let pts = vec![
        Point::new(x0, y0),
        Point::new(x1, y0),
        Point::new(x1, y1),
        Point::new(x0, y1),
    ];
    SolidModel {
        regions: vec![extrude_polygon(&pts, z0, z1, 1, None, None, None)],
        material_names: vec!["A".into()],
    }
}

#[test]
fn clipping_a_box_to_a_smaller_window() {
    let model = box_region(0.0, 0.0, 10.0, 10.0, 0.0, 1.0);
    let window = Bbox {
        min_x: 2.0,
        min_y: 2.0,
        max_x: 4.0,
        max_y: 4.0,
    };
    let clipped = clip_model(&model, window);
    assert_eq!(clipped.regions.len(), 1);
    assert!((clipped.total_volume() - 4.0).abs() < 1e-12);
    match &clipped.regions[0].solid {
        Solid::Box {
            min_x,
            min_y,
            max_x,
            max_y,
            ..
        } => {
            assert_eq!((*min_x, *min_y, *max_x, *max_y), (2.0, 2.0, 4.0, 4.0));
        }
        _ => panic!("expected a box"),
    }
}

#[test]
fn clipping_splits_a_polygon_crossing_the_window() {
    // An L-shaped polygon crossing the window edge gets split into rectangles.
    let pts = vec![
        Point::new(0.0, 0.0),
        Point::new(6.0, 0.0),
        Point::new(6.0, 2.0),
        Point::new(2.0, 2.0),
        Point::new(2.0, 6.0),
        Point::new(0.0, 6.0),
    ];
    let model = SolidModel {
        regions: vec![extrude_polygon(&pts, 0.0, 1.0, 1, None, None, None)],
        material_names: vec!["A".into()],
    };
    let window = Bbox {
        min_x: 1.0,
        min_y: 1.0,
        max_x: 3.0,
        max_y: 3.0,
    };
    let clipped = clip_model(&model, window);
    // Inside [1,3]x[1,3] the L covers [1,3]x[1,2] (area 2) plus [1,2]x[2,3]
    // (area 1) = 3 (the top-right [2,3]x[2,3] is outside the L).
    assert!((clipped.total_volume() - 3.0).abs() < 1e-12, "vol={}", clipped.total_volume());
    assert!(!clipped.regions.is_empty());
}

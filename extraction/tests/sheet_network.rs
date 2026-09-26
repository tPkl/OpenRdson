//! Correctness of the 2.5D sheet-resistance network on a known square.

use openrdson_core::geometry::Point;
use openrdson_extraction::{SheetNetwork, SheetPolygon};

fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Vec<Point> {
    vec![
        Point::new(x0, y0),
        Point::new(x1, y0),
        Point::new(x1, y1),
        Point::new(x0, y1),
        Point::new(x0, y0),
    ]
}

#[test]
fn square_sheet_has_one_ohm_per_square() {
    // 10x10 unit square, sheet resistance 1 ohm/sq -> R = 1 ohm.
    let poly = SheetPolygon {
        physical: "M".into(),
        points: rect(0.0, 0.0, 10.0, 10.0),
        component: 0,
        r_sheet: 1.0,
    };
    let net = SheetNetwork::build(&[poly], &[], &[], &[], 0.5, 200);
    assert!(net.node_count() > 0);

    let left: Vec<u32> = (0..net.node_count() as u32)
        .filter(|&n| net.x[n as usize].abs() < 1e-12)
        .collect();
    let right: Vec<u32> = (0..net.node_count() as u32)
        .filter(|&n| (net.x[n as usize] - 10.0).abs() < 1e-12)
        .collect();
    assert!(!left.is_empty() && !right.is_empty());

    let r = net.resistance(&left, &right, 1e-10, 5000).unwrap();
    assert!((r - 1.0).abs() / 1.0 < 0.05, "R = {r}, expected ~1 ohm/sq");
}

#[test]
fn doubling_length_doubles_resistance() {
    let poly = SheetPolygon {
        physical: "M".into(),
        points: rect(0.0, 0.0, 20.0, 10.0),
        component: 0,
        r_sheet: 1.0,
    };
    let net = SheetNetwork::build(&[poly], &[], &[], &[], 0.5, 200);
    let left: Vec<u32> = (0..net.node_count() as u32)
        .filter(|&n| net.x[n as usize].abs() < 1e-12)
        .collect();
    let right: Vec<u32> = (0..net.node_count() as u32)
        .filter(|&n| (net.x[n as usize] - 20.0).abs() < 1e-12)
        .collect();
    let r = net.resistance(&left, &right, 1e-10, 5000).unwrap();
    // 2 squares -> ~2 ohm.
    assert!((r - 2.0).abs() / 2.0 < 0.05, "R = {r}, expected ~2 ohm/sq");
}

#[test]
fn two_isolated_nets_are_open() {
    let a = SheetPolygon { physical: "M".into(), points: rect(0.0, 0.0, 10.0, 10.0), component: 0, r_sheet: 1.0 };
    let b = SheetPolygon { physical: "M".into(), points: rect(20.0, 0.0, 30.0, 10.0), component: 1, r_sheet: 1.0 };
    let net = SheetNetwork::build(&[a, b], &[], &[], &[], 1.0, 100);
    let left: Vec<u32> = (0..net.node_count() as u32)
        .filter(|&n| net.component[n as usize] == 0 && net.x[n as usize].abs() < 1e-12)
        .collect();
    let right: Vec<u32> = (0..net.node_count() as u32)
        .filter(|&n| net.component[n as usize] == 1 && (net.x[n as usize] - 20.0).abs() < 1e-12)
        .collect();
    assert!(!left.is_empty() && !right.is_empty());
    let r = net.resistance_with_edges(&left, &right, &[], 1e-8, 2000).unwrap();
    println!("isolated nets R = {r:e}");
    assert!(r.is_infinite() || r > 1e9, "expected open, got {r:e}");
}

#[test]
fn via_grouping_combines_parallel_vias() {
    use openrdson_extraction::{group_vias, SheetVia};
    // Four 4-ohm vias close together -> one equivalent 1-ohm via.
    let vias: Vec<SheetVia> = (0..4)
        .map(|i| SheetVia {
            bottom: "M1".into(),
            top: "M2".into(),
            component: 0,
            center: Point::new(i as f64 * 0.1, 0.0),
            resistance: 4.0,
            half_width: 0.0,
            half_height: 0.0,
        })
        .collect();
    let g = group_vias(&vias, 1.0);
    assert_eq!(g.len(), 1, "expected a single grouped via");
    assert!((g[0].resistance - 1.0).abs() < 1e-12, "R={}", g[0].resistance);

    // Far-apart vias stay separate.
    let far: Vec<SheetVia> = (0..3)
        .map(|i| SheetVia {
            bottom: "M1".into(),
            top: "M2".into(),
            component: 0,
            center: Point::new(i as f64 * 100.0, 0.0),
            resistance: 1.0,
            half_width: 0.0,
            half_height: 0.0,
        })
        .collect();
    assert_eq!(group_vias(&far, 1.0).len(), 3);
}

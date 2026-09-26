use openrdson_core::mesh::{Element, ElementType, Mesh, Point3};
use openrdson_io::{
    cell_current_density, cell_power_current, write_vtu, VtkCellType, VtuMesh,
};

#[test]
fn writes_unstructured_grid_with_point_and_cell_data() {
    let path = std::env::temp_dir().join("openrdson_vtu_test.vtu");
    let vm = VtuMesh {
        points: vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.5, 0.5, 0.0],
        ],
        cells: vec![
            (VtkCellType::Quad, vec![0, 1, 2, 3]),
            (VtkCellType::Triangle, vec![0, 1, 4]),
        ],
        point_data: vec![("potential_V".into(), vec![0.0, 1.0, 2.0, 3.0, 0.5])],
        cell_data: vec![("current_density_A_m".into(), vec![10.0, 20.0])],
    };
    write_vtu(&path, &vm).unwrap();
    let text = std::fs::read_to_string(&path).unwrap();
    assert!(text.contains("NumberOfPoints=\"5\""));
    assert!(text.contains("NumberOfCells=\"2\""));
    assert!(text.contains("Name=\"connectivity\""));
    assert!(text.contains("Name=\"offsets\""));
    assert!(text.contains("Name=\"potential_V\""));
    assert!(text.contains("Name=\"current_density_A_m\""));

    // Cell types: Quad = 9, Triangle = 5.
    let seg = text.split("Name=\"types\"").nth(1).unwrap();
    let seg = seg.split("</DataArray>").next().unwrap();
    let nums: Vec<i32> = seg
        .split(|c: char| !c.is_ascii_digit())
        .filter(|s| !s.is_empty())
        .filter_map(|s| s.parse().ok())
        .collect();
    assert_eq!(nums, vec![9, 5]);
    std::fs::remove_file(&path).ok();
}

#[test]
fn subset_cells_remaps_points_and_data() {
    let vm = VtuMesh {
        points: vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [9.0, 9.0, 0.0],
        ],
        cells: vec![
            (VtkCellType::Quad, vec![0, 1, 2, 3]),
            (VtkCellType::Quad, vec![1, 4, 2, 0]),
        ],
        point_data: vec![("potential_V".into(), vec![0.0, 1.0, 2.0, 3.0, 4.0])],
        cell_data: vec![("layer_id".into(), vec![0.0, 1.0])],
    };
    let sub = vm.subset_cells(&[false, true]);
    assert_eq!(sub.cells.len(), 1);
    // Cell 1 references points 1, 4, 2, 0 -> 4 unique, in encounter order.
    assert_eq!(sub.points.len(), 4);
    assert_eq!(sub.points[0], [1.0, 0.0, 0.0]);
    assert_eq!(sub.point_data[0].1, vec![1.0, 4.0, 2.0, 0.0]);
    assert_eq!(sub.cells[0].1, vec![0, 1, 2, 3]);
    assert_eq!(sub.cell_data[0].1, vec![1.0]);
}

#[test]
fn current_density_of_linear_field_on_a_cube() {
    // Unit cube, u = x -> grad u = (1, 0, 0), so |J| = sigma.
    let nodes = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(1.0, 1.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(0.0, 0.0, 1.0),
        Point3::new(1.0, 0.0, 1.0),
        Point3::new(1.0, 1.0, 1.0),
        Point3::new(0.0, 1.0, 1.0),
    ];
    let mesh = Mesh {
        nodes: nodes.clone(),
        elements: vec![Element {
            kind: ElementType::Hex,
            nodes: vec![0, 1, 2, 3, 4, 5, 6, 7],
            material_id: 1,
            net: None,
            layer: None,
            device_ref: None,
            source_polygon: None,
        }],
        sets: Default::default(),
    };
    let potential: Vec<f64> = nodes.iter().map(|p| p.x).collect();
    let sigma = vec![0.0, 2.0]; // material id 1 -> 2 S/m
    let j = cell_current_density(&mesh, &potential, &sigma);
    assert_eq!(j.len(), 1);
    assert!((j[0] - 2.0).abs() < 1e-12, "got {}", j[0]);
}

#[test]
fn power_and_current_of_linear_field_on_a_cube() {
    // Unit cube, u = x, sigma = 2: power = sigma|grad|^2 V = 2 W, dV = 1 -> I = 2 A.
    let nodes = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(1.0, 1.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(0.0, 0.0, 1.0),
        Point3::new(1.0, 0.0, 1.0),
        Point3::new(1.0, 1.0, 1.0),
        Point3::new(0.0, 1.0, 1.0),
    ];
    let mesh = Mesh {
        nodes: nodes.clone(),
        elements: vec![Element {
            kind: ElementType::Hex,
            nodes: vec![0, 1, 2, 3, 4, 5, 6, 7],
            material_id: 1,
            net: None,
            layer: None,
            device_ref: None,
            source_polygon: None,
        }],
        sets: Default::default(),
    };
    let potential: Vec<f64> = nodes.iter().map(|p| p.x).collect();
    let sigma = vec![0.0, 2.0];
    let (p, i) = cell_power_current(&mesh, &potential, &sigma);
    assert_eq!(p.len(), 1);
    assert!((p[0] - 2.0).abs() < 1e-12, "power = {}", p[0]);
    assert!((i[0] - 2.0).abs() < 1e-12, "current = {}", i[0]);
}

#[test]
fn rejects_mismatched_point_data() {
    let path = std::env::temp_dir().join("openrdson_vtu_bad.vtu");
    let vm = VtuMesh {
        points: vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        cells: vec![],
        point_data: vec![("potential_V".into(), vec![0.0])],
        cell_data: vec![],
    };
    assert!(write_vtu(&path, &vm).is_err());
}

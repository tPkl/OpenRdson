//! GDSII writer + KLayout colormap round-trip tests.

use openrdson_io::{
    read_gds_file, write_colormap, write_gds_file, write_grouped_colormap, GdsBoundary, VizItem,
};
use std::path::PathBuf;

fn tmp(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("openrdson_test_{}_{name}", std::process::id()))
}

#[test]
fn gds_write_read_round_trip() {
    let path = tmp("roundtrip.gds");
    let cells = vec![(
        "TOP".to_string(),
        vec![GdsBoundary {
            layer: 5,
            datatype: 0,
            points: vec![(0, 0), (100, 0), (100, 50), (0, 50)],
        }],
    )];
    write_gds_file(&path, "LIB", 1e-9, &cells).unwrap();

    let lib = read_gds_file(&path).unwrap();
    assert_eq!(lib.name, "LIB");
    assert_eq!(lib.structs.len(), 1);
    assert_eq!(lib.structs[0].name, "TOP");
    assert_eq!(lib.structs[0].elements.len(), 1);
    assert_eq!(lib.structs[0].elements[0].layer, 5);
    assert_eq!(lib.structs[0].elements[0].xy.len(), 5); // closed ring
    assert!((lib.db_unit_meters - 1e-9).abs() < 1e-24);
    std::fs::remove_file(&path).ok();
}

#[test]
fn colormap_writes_gds_and_lyp() {
    let path = tmp("cmap.gds");
    let items = vec![
        VizItem {
            points: vec![(0, 0), (10, 0), (10, 10), (0, 10)],
            value: 0.0,
        },
        VizItem {
            points: vec![(10, 0), (20, 0), (20, 10), (10, 10)],
            value: 1.0,
        },
        VizItem {
            points: vec![(20, 0), (30, 0), (30, 10), (20, 10)],
            value: 0.5,
        },
    ];
    let bins = write_colormap(&path, "LIB", "TOP", 1e-9, &items, 8, 100, 0.0, 1.0).unwrap();
    assert_eq!(bins.len(), 8);

    let lib = read_gds_file(&path).unwrap();
    let layers: std::collections::BTreeSet<i16> =
        lib.structs[0].elements.iter().map(|e| e.layer).collect();
    // min value -> layer 100, max -> 107, mid -> 103/104.
    assert!(layers.contains(&100));
    assert!(layers.contains(&107));
    assert_eq!(lib.structs[0].elements.len(), 3);

    // KLayout auto-loads <base>.lyp (same name, no .gds) and matches GDS layers
    // via the `<layer>/<datatype>@<view>` source format.
    let lyp = path.with_extension("lyp");
    assert!(lyp.exists(), "expected {}", lyp.display());
    let text = std::fs::read_to_string(&lyp).unwrap();
    assert!(text.contains("<source>100/0@*</source>"), "bad source: {text}");
    assert!(text.contains("<source>107/0@*</source>"));
    std::fs::remove_file(&path).ok();
    std::fs::remove_file(&lyp).ok();
}

#[test]
fn grouped_colormap_has_one_group_per_layer_with_bins() {
    let path = tmp("grouped.gds");
    let layers = vec![
        (
            "Metal1".to_string(),
            vec![
                VizItem { points: vec![(0, 0), (10, 0), (10, 10), (0, 10)], value: 0.0 },
                VizItem { points: vec![(10, 0), (20, 0), (20, 10), (10, 10)], value: 1.0 },
            ],
        ),
        (
            "Poly1".to_string(),
            vec![VizItem { points: vec![(0, 0), (5, 0), (5, 5), (0, 5)], value: 0.5 }],
        ),
    ];
    write_grouped_colormap(&path, "LIB", "TOP", 1e-9, &layers, 32, 1000, 32).unwrap();

    let lib = read_gds_file(&path).unwrap();
    let ls: std::collections::BTreeSet<i16> =
        lib.structs[0].elements.iter().map(|e| e.layer).collect();
    // Metal1 bin 0 -> 1000, bin 31 -> 1031; Poly1 -> 1032..1063.
    assert!(ls.contains(&1000));
    assert!(ls.contains(&1031));
    assert!(ls.contains(&1032));

    let lyp = path.with_extension("lyp");
    let text = std::fs::read_to_string(&lyp).unwrap();
    // Two groups (Metal1, Poly1); each child is its own <group-members> block.
    // Metal1: values 0.0 and 1.0 -> bins 0 and 31 -> layers 1000, 1031.
    // Poly1:  value 0.5         -> bin 0         -> layer 1032.
    assert_eq!(text.matches("<properties>").count(), 2, "one group per layer");
    // 2 layers x 32 bins, each bin its own <group-members> block.
    assert_eq!(text.matches("<group-members>").count(), 64, "one block per bin");
    assert!(text.contains("<name>Metal1</name>"));
    assert!(text.contains("<name>Poly1</name>"));
    assert!(text.contains("<source>1000/0</source>"));
    assert!(text.contains("<source>1031/0</source>"));
    assert!(text.contains("<source>1032/0</source>"));
    std::fs::remove_file(&path).ok();
    std::fs::remove_file(&lyp).ok();
}

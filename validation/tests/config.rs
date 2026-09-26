//! Meshing-settings config parsing.

use openrdson_validation::MeshSettings;

#[test]
fn mesh_config_parses_layer_overrides() {
    let mut s = MeshSettings::default();
    s.load_config_str("# layer-based refinement (microns)\nPoly1 0.5\nMetal3 1.0\n");
    assert_eq!(s.per_layer_cell_m.get("Poly1").copied(), Some(0.5e-6));
    assert_eq!(s.per_layer_cell_m.get("Metal3").copied(), Some(1.0e-6));
    // Comments and blank lines ignored; unknown/garbage lines skipped.
    s.load_config_str("\n#comment\nbadline\n");
    assert_eq!(s.per_layer_cell_m.len(), 2);
}

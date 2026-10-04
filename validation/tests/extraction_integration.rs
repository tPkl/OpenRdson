//! End-to-end extraction: minimal ICT + GDS fixtures run through the full
//! `SheetFullArray` pipeline, checked against analytic resistance values.
//!
//! These are the "integration" tests: they exercise the ICT parser, the GDS
//! layer map, connectivity extraction, the 2.5D sheet mesh and the solve as one
//! path, so a regression anywhere in that chain fails here.

use openrdson_io::{write_gds_file, GdsBoundary};
use openrdson_validation::{BiasTerminal, MeshSettings, SheetFullArray};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

const UM: f64 = 1e-6;
static DIR_N: AtomicU64 = AtomicU64::new(0);

/// A minimal two-metal process: both metals 1 Ω/sq, and a reference via of
/// 1 Ω at 0.25 µm² (so a via of area `A` µm² has resistance `0.25 / A` Ω).
const ICT: &str = r#"process "test" {
    temp_reference 25
}
conductor "M1" {
    resistivity 1.0
}
conductor "M2" {
    resistivity 1.0
}
via "V1" {
    top_layer "M2"
    bottom_layer "M1"
    area_resistance 1.0 0.25
}
"#;

const GDS_MAP: &str = "net_M1 10 0\nnet_M2 20 0\nvia1 30 0\n";
const CCI_MAP: &str = "conducting_layers\nnet_M1 M1\nnet_M2 M2\nvia_layers\nvia1 V1\n";

fn temp_dir() -> PathBuf {
    let d = std::env::temp_dir().join(format!(
        "openrdson_it_{}_{}",
        std::process::id(),
        DIR_N.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn write_fixture(dir: &Path) {
    std::fs::write(dir.join("stack.ict"), ICT).unwrap();
    std::fs::write(dir.join("layout.gds.map"), GDS_MAP).unwrap();
    std::fs::write(dir.join("sft_cci.map"), CCI_MAP).unwrap();
    std::fs::write(dir.join("layout.ports"), "").unwrap();
    std::fs::write(dir.join("layout.devtab"), "").unwrap();
    std::fs::write(dir.join("layout.spi"), "").unwrap();
    std::fs::write(dir.join("model.csv"), "").unwrap();
}

/// A rectangular GDS boundary on `layer`; coordinates in µm (db unit = 1 µm).
fn rect(layer: i16, x0: f64, y0: f64, x1: f64, y1: f64) -> GdsBoundary {
    GdsBoundary {
        layer,
        datatype: 0,
        points: vec![
            (x0 as i32, y0 as i32),
            (x1 as i32, y0 as i32),
            (x1 as i32, y1 as i32),
            (x0 as i32, y1 as i32),
        ],
    }
}

/// A full-width contact on `layer` at `x_um`. The footprint spans the whole 4 µm
/// strip height, so the terminal is a bar contact (no point-contact spreading).
fn terminal(name: &str, layer: &str, x_um: f64) -> BiasTerminal {
    BiasTerminal {
        name: name.into(),
        voltage: if name.eq_ignore_ascii_case("d") { 1.0 } else { 0.0 },
        x_um,
        y_um: 2.0,
        // Narrow in x (single edge column) but full height, so the terminal is
        // a bar contact with no point-contact spreading and no length loss.
        dx_um: 0.2,
        dy_um: 4.0,
        layer: Some(layer.into()),
    }
}

fn settings(dir: &Path) -> MeshSettings {
    let mut s = MeshSettings::default();
    s.paths.cci_dir = dir.to_path_buf();
    s.paths.layout = "layout.gds".into();
    s.paths.gds_map = "layout.gds.map".into();
    s.paths.ports = "layout.ports".into();
    s.paths.devtab = "layout.devtab".into();
    s.paths.spi = "layout.spi".into();
    s.paths.tech_ict = dir.join("stack.ict");
    s.paths.layer_map = dir.join("sft_cci.map");
    s.paths.model_csv = dir.join("model.csv");
    s.temperature_c = 27.0;
    s.sheet_cell_m = 0.5 * UM;
    s.sheet_cap = 200;
    s.via_group_radius_m = 0.0; // keep vias distinct so parallel combination is exact
    s
}

fn extract_rds(s: &MeshSettings) -> f64 {
    let fa = SheetFullArray::build_with(s).unwrap_or_else(|e| panic!("build failed: {e}"));
    fa.rds_open()
        .unwrap_or_else(|e| panic!("rds_open failed: {e}"))
}

fn assert_close(actual: f64, expected: f64, tol: f64, ctx: &str) {
    let err = (actual - expected).abs() / expected;
    assert!(
        err < tol,
        "{ctx}: R = {actual:.6} ohm, expected {expected:.6} ohm ({:.2}% off, > {:.1}% tol)",
        err * 100.0,
        tol * 100.0
    );
}

#[test]
fn metal_line_resistance_matches_sheet_resistance() {
    let dir = temp_dir();
    write_fixture(&dir);
    // 40 x 4 µm M1 strip, 1 Ω/sq -> R = Rs * L/W = 10 Ω.
    write_gds_file(
        dir.join("layout.gds"),
        "TOP",
        UM,
        &[("TOP".into(), vec![rect(10, 0.0, 0.0, 40.0, 4.0)])],
    )
    .unwrap();

    let mut s = settings(&dir);
    s.terminals = vec![terminal("D", "M1", 0.0), terminal("S", "M1", 40.0)];
    assert_close(extract_rds(&s), 10.0, 0.02, "metal line");
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn m1_via_m2_series_resistance_matches_analytic() {
    let dir = temp_dir();
    write_fixture(&dir);
    // M1 and M2 overlap; a full-width 1 x 4 µm via at x=20 connects them.
    write_gds_file(
        dir.join("layout.gds"),
        "TOP",
        UM,
        &[(
            "TOP".into(),
            vec![
                rect(10, 0.0, 0.0, 40.0, 4.0),
                rect(20, 0.0, 0.0, 40.0, 4.0),
                rect(30, 20.0, 0.0, 21.0, 4.0),
            ],
        )],
    )
    .unwrap();

    let mut s = settings(&dir);
    s.terminals = vec![terminal("D", "M1", 0.0), terminal("S", "M2", 40.0)];
    // Metal: 1 Ω/sq * 40/4 = 10 Ω. Via: 0.25 Ω·µm² ref, area 4 µm² -> 0.0625 Ω.
    // The via's finite 1 µm x-footprint overlaps ~1 µm of metal, so the series
    // answer carries ~1% modeling uncertainty; a 3% tolerance is still far below
    // any single-parser/connectivity regression (which shows up as open or ~4x).
    assert_close(extract_rds(&s), 10.0625, 0.03, "M1-via-M2");
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn parallel_vias_combine_like_one_wide_via() {
    let dir = temp_dir();
    write_fixture(&dir);
    let metals = vec![rect(10, 0.0, 0.0, 40.0, 4.0), rect(20, 0.0, 0.0, 40.0, 4.0)];
    let mut s = settings(&dir);
    s.terminals = vec![terminal("D", "M1", 0.0), terminal("S", "M2", 40.0)];

    // One full-width via (area 4 µm² -> 0.0625 Ω).
    let mut one = metals.clone();
    one.push(rect(30, 20.0, 0.0, 21.0, 4.0));
    write_gds_file(dir.join("layout.gds"), "TOP", UM, &[("TOP".into(), one)]).unwrap();
    let r_wide = extract_rds(&s);

    // Two half-width vias in parallel (each 1 x 2 µm -> 0.125 Ω -> 0.0625 Ω).
    let mut two = metals.clone();
    two.push(rect(30, 20.0, 0.0, 21.0, 2.0));
    two.push(rect(30, 20.0, 2.0, 21.0, 4.0));
    write_gds_file(dir.join("layout.gds"), "TOP", UM, &[("TOP".into(), two)]).unwrap();
    let r_parallel = extract_rds(&s);

    assert_close(r_parallel, r_wide, 0.02, "parallel vias");
    std::fs::remove_dir_all(&dir).ok();
}

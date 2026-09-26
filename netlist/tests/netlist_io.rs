//! Round-trip and format tests for netlist IO.

use openrdson_netlist::{
    parse_spice_resistors, write_spef, write_spice_resistors, write_spice_two_terminal,
    NamedResistor,
};

#[test]
fn spice_resistor_round_trip() {
    let rs = vec![
        NamedResistor {
            a: "D".into(),
            b: "S".into(),
            resistance: 0.265e-3,
        },
        NamedResistor {
            a: "G".into(),
            b: "S".into(),
            resistance: 1.0e9,
        },
    ];
    let text = write_spice_resistors("test", &rs);
    assert!(text.contains(".SUBCKT PARASITIC_DUT"));
    let parsed = parse_spice_resistors(&text);
    assert_eq!(parsed.len(), 2);
    assert_eq!(parsed[0].a, "D");
    assert!((parsed[0].resistance - 0.265e-3).abs() < 1e-15);
}

#[test]
fn spef_contains_res_and_values() {
    let rs = vec![NamedResistor {
        a: "D".into(),
        b: "S".into(),
        resistance: 1.23,
    }];
    let spef = write_spef("design", &rs);
    assert!(spef.contains("*SPEF"));
    assert!(spef.contains("*R_UNIT 1 OHM"));
    assert!(spef.contains("*RES"));
    assert!(spef.contains("1.230000e0"));
}

#[test]
fn two_terminal_spice_emits_single_resistor() {
    use openrdson_extraction::ParasiticNetwork;
    let net = ParasiticNetwork {
        terminals: vec!["D".into(), "S".into()],
        rmatrix: vec![vec![0.0, 0.5], vec![0.5, 0.0]],
        resistors: Vec::new(),
    };
    let text = write_spice_two_terminal(&net, "dut");
    assert!(text.contains("R1 D S 5.000000e-1"));
    assert!(text.contains("D <-> S = 5.000000e-1"));
}

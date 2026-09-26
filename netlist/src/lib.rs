//! Netlist & parasitic IO for OpenRDSon.
//!
//! v0 writes/reads two-terminal resistor networks in SPICE and SPEF. The
//! reduced `ParasiticNetwork` from extraction is emitted as a single resistor
//! for a two-terminal pair (with the full matrix preserved as comments).

use openrdson_extraction::ParasiticNetwork;
use std::fmt::Write as _;

/// A named two-terminal resistor.
#[derive(Debug, Clone, PartialEq)]
pub struct NamedResistor {
    pub a: String,
    pub b: String,
    pub resistance: f64,
}

/// Write a SPICE resistor subcircuit.
pub fn write_spice_resistors(title: &str, resistors: &[NamedResistor]) -> String {
    let mut s = String::new();
    let _ = writeln!(s, "* {title}");
    let _ = writeln!(s, ".SUBCKT PARASITIC_DUT");
    for (i, r) in resistors.iter().enumerate() {
        let _ = writeln!(
            s,
            "R{i} {} {} {:.6e}",
            sanitize(&r.a),
            sanitize(&r.b),
            r.resistance
        );
    }
    let _ = writeln!(s, ".ENDS PARASITIC_DUT");
    s
}

/// Write a SPICE deck for a reduced two-terminal network.
pub fn write_spice_two_terminal(net: &ParasiticNetwork, title: &str) -> String {
    let mut s = String::new();
    let _ = writeln!(s, "* {title}");
    let _ = writeln!(s, "* reduced terminal resistance matrix (ohms)");
    for (i, ti) in net.terminals.iter().enumerate() {
        for (j, tj) in net.terminals.iter().enumerate() {
            if j > i {
                let _ = writeln!(s, "*   {} <-> {} = {:.6e}", ti, tj, net.rmatrix[i][j]);
            }
        }
    }
    if net.terminals.len() == 2 {
        let _ = writeln!(
            s,
            "R1 {} {} {:.6e}",
            sanitize(&net.terminals[0]),
            sanitize(&net.terminals[1]),
            net.rmatrix[0][1]
        );
    }
    s
}

/// Parse `R<name> <a> <b> <value>` resistor lines from SPICE text.
pub fn parse_spice_resistors(text: &str) -> Vec<NamedResistor> {
    let mut out = Vec::new();
    for line in text.lines() {
        let line = line.trim();
        if line.starts_with('*') || line.starts_with('.') {
            continue;
        }
        let t: Vec<&str> = line.split_whitespace().collect();
        if t.len() < 4 {
            continue;
        }
        if !t[0].starts_with('R') && !t[0].starts_with('r') {
            continue;
        }
        if let Ok(v) = t[3].parse::<f64>() {
            out.push(NamedResistor {
                a: t[1].to_string(),
                b: t[2].to_string(),
                resistance: v,
            });
        }
    }
    out
}

/// Write a minimal SPEF file for a set of named resistors.
pub fn write_spef(design: &str, resistors: &[NamedResistor]) -> String {
    let mut s = String::new();
    let _ = writeln!(s, "*SPEF \"IEEE 1481-1998\"");
    let _ = writeln!(s, "*DESIGN \"{design}\"");
    let _ = writeln!(s, "*DATE \"today\"");
    let _ = writeln!(s, "*VERSION \"1.0\"");
    let _ = writeln!(s, "*DESIGN_FLOW \"NAME_SCOPE LOCAL\" \"PIN_CAP NONE\"");
    let _ = writeln!(s, "*DIVIDER /");
    let _ = writeln!(s, "*DELIMITER :");
    let _ = writeln!(s, "*BUS_DELIMITER [ ]");
    let _ = writeln!(s, "*T_UNIT 1 NS");
    let _ = writeln!(s, "*C_UNIT 1 PF");
    let _ = writeln!(s, "*R_UNIT 1 OHM");
    let _ = writeln!(s, "*L_UNIT 1 HENRY");
    let _ = writeln!(s);
    let _ = writeln!(s, "*NAME_MAP");
    for r in resistors {
        let _ = writeln!(s, "*{} {}", sanitize(&r.a), sanitize(&r.a));
    }
    let _ = writeln!(s);
    let _ = writeln!(s, "*D_NET net");
    let _ = writeln!(s, "*CONN");
    let _ = writeln!(s, "*RES");
    for (i, r) in resistors.iter().enumerate() {
        let _ = writeln!(
            s,
            "{i} {} {} {:.6e}",
            sanitize(&r.a),
            sanitize(&r.b),
            r.resistance
        );
    }
    let _ = writeln!(s, "*END");
    s
}

fn sanitize(name: &str) -> String {
    name.chars()
        .map(|c| if c.is_alphanumeric() || c == '_' { c } else { '_' })
        .collect()
}

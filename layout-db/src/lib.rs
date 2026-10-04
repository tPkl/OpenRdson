//! Layout database ingestion for OpenRDSon.
//!
//! Turns an annotated GDSII/AGF file plus its `*.gds.map` layer map, ports and
//! `devtab` device templates into the shared [`LayoutIR`] / [`DeviceTemplate`]
//! model. Generic device *recognition* lives in `device-recognition`.

pub mod connectivity;

pub use connectivity::{
    apply_channel_cut, channel_regions, extract_connectivity, Component, Connectivity,
};

use openrdson_core::device::{DeviceKind, DeviceTemplate, SpiInstance};
use openrdson_core::diag::Diag;
use openrdson_core::geometry::Point;
use openrdson_core::layout::{LayoutIR, LayoutPolygon, Port};
use openrdson_io::{read_gds_file, read_gds_map_file, ElementKind, GdsLibrary};
use std::collections::BTreeMap;
use std::io;
use std::path::Path;

/// Build a [`LayoutIR`] from an AGF/GDSII file and its `*.gds.map`.
pub fn load_layout_ir<P: AsRef<Path>, Q: AsRef<Path>>(
    agf_path: P,
    gds_map_path: Q,
) -> io::Result<LayoutIR> {
    let lib = read_gds_file(agf_path)?;
    let layer_names = read_gds_map_file(gds_map_path)?;
    Ok(library_to_layout_ir(&lib, &layer_names))
}

/// Convert a parsed [`GdsLibrary`] into a [`LayoutIR`] using a layer map.
pub fn library_to_layout_ir(
    lib: &GdsLibrary,
    layer_names: &BTreeMap<i16, String>,
) -> LayoutIR {
    // Use the first structure as the top cell.
    let top = lib.structs.first();
    let cell = top.map(|s| s.name.clone()).unwrap_or_default();
    let mut polygons = Vec::new();

    if let Some(top) = top {
        for e in &top.elements {
            // v0 handles polygonal conducting geometry.
            if !matches!(e.kind, ElementKind::Boundary | ElementKind::Path) {
                continue;
            }
            if e.xy.len() < 3 {
                continue;
            }
            let points: Vec<Point> = e
                .xy
                .iter()
                .map(|&(x, y)| {
                    Point::new(x as f64 * lib.db_unit_meters, y as f64 * lib.db_unit_meters)
                })
                .collect();
            let logical = layer_names.get(&e.layer);
            let net = logical
                .filter(|n| n.starts_with("net_"))
                .cloned();
            polygons.push(LayoutPolygon {
                layer: openrdson_core::geometry::LayerKey::new(e.layer, e.datatype),
                net,
                points,
                props: e.props.clone(),
            });
        }
    }

    LayoutIR {
        cell,
        db_unit_meters: lib.db_unit_meters,
        polygons,
        ports: Vec::new(),
        layer_names: layer_names.clone(),
    }
}

/// Scale every coordinate in a layout by `factor`, updating `db_unit_meters`.
///
/// Used to correct a malformed GDS/AGF `UNITS` record (e.g. both doubles equal,
/// implying a 1 m user unit) once the true database unit is known.
pub fn rescale_layout(layout: &mut openrdson_core::layout::LayoutIR, factor: f64) {
    if !(factor.is_finite() && factor > 0.0) || (factor - 1.0).abs() < 1e-15 {
        return;
    }
    for p in &mut layout.polygons {
        for pt in &mut p.points {
            pt.x *= factor;
            pt.y *= factor;
        }
    }
    for port in &mut layout.ports {
        port.x *= factor;
        port.y *= factor;
    }
    layout.db_unit_meters *= factor;
}

/// Parse a CCI `*.ports` file. Coordinates are database units and are converted
/// to meters with `db_unit_meters`. Returns the ports plus diagnostics for
/// malformed lines.
pub fn parse_ports(text: &str, db_unit_meters: f64) -> (Vec<Port>, Vec<Diag>) {
    let mut ports = Vec::new();
    let mut diags = Vec::new();
    for (i, raw) in text.lines().enumerate() {
        let line_no = i + 1;
        let line = raw.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let t: Vec<&str> = line.split_whitespace().collect();
        // <name> <id> <direction> <x> <y> <net>
        if t.len() < 6 {
            diags.push(Diag::at(
                line_no,
                format!(
                    "expected 6 fields \"<name> <id> <direction> <x> <y> <net>\", got {}",
                    t.len()
                ),
            ));
            continue;
        }
        let (Ok(x), Ok(y)) = (t[3].parse::<f64>(), t[4].parse::<f64>()) else {
            diags.push(Diag::at(
                line_no,
                format!("non-numeric coordinate (x=\"{}\", y=\"{}\")", t[3], t[4]),
            ));
            continue;
        };
        // The trailing token is the port's logical layer/net name in this DB.
        ports.push(Port {
            name: t[0].to_string(),
            id: t[1].parse::<i64>().ok(),
            direction: Some(t[2].to_string()),
            x: x * db_unit_meters,
            y: y * db_unit_meters,
            net: t[5].to_string(),
            layer: t[5].to_string(),
        });
    }
    (ports, diags)
}

pub fn read_ports_file<P: AsRef<Path>>(path: P, db_unit_meters: f64) -> io::Result<Vec<Port>> {
    let path = path.as_ref();
    let text = std::fs::read_to_string(path)?;
    let (ports, diags) = parse_ports(&text, db_unit_meters);
    openrdson_core::diag::log_diags(&path.display().to_string(), &diags);
    Ok(ports)
}

/// Parse a CCI `devtab` device table into device *definitions* (templates),
/// plus diagnostics for entries that could not be parsed.
pub fn parse_devtab(text: &str) -> (Vec<DeviceTemplate>, Vec<Diag>) {
    let lines: Vec<&str> = text.lines().collect();
    let starts: Vec<usize> = lines
        .iter()
        .enumerate()
        .filter(|(_, l)| l.trim_start().starts_with("Device Entry"))
        .map(|(i, _)| i)
        .collect();

    let mut templates = Vec::new();
    let mut diags = Vec::new();
    for (n, &start) in starts.iter().enumerate() {
        let end = starts.get(n + 1).copied().unwrap_or(lines.len());
        // Skip the "Device Entry N" line and the following timestamp line.
        let body = &lines[(start + 2).min(end)..end];
        match parse_device_entry(body) {
            Some(t) => templates.push(t),
            None => diags.push(Diag::at(
                start + 1,
                "malformed 'Device Entry' block (missing name/kind/model or terminal/property/param counts)",
            )),
        }
    }
    (templates, diags)
}

fn parse_device_entry(lines: &[&str]) -> Option<DeviceTemplate> {
    let mut it = lines.iter().map(|s| s.trim()).filter(|s| !s.is_empty());
    let name = it.next()?.to_string();
    let kind = DeviceKind::parse(it.next()?);
    let model_raw = it.next()?;
    let model = if model_raw == "(null)" {
        None
    } else {
        Some(model_raw.to_string())
    };
    // Three reserved "(null)" fields.
    it.next();
    it.next();
    it.next();

    let term_count: usize = it.next()?.parse().ok()?;
    let mut terminals = Vec::new();
    for _ in 0..term_count {
        let l = it.next()?;
        let t: Vec<&str> = l.split_whitespace().collect();
        if t.len() >= 2 {
            terminals.push((t[0].to_string(), t[1].to_string()));
        }
    }

    let prop_count: usize = it.next()?.parse().ok()?;
    let mut property_layers = Vec::new();
    for _ in 0..prop_count {
        property_layers.push(it.next()?.to_string());
    }

    let param_count: usize = it.next()?.parse().ok()?;
    let mut params = Vec::new();
    for _ in 0..param_count {
        params.push(it.next()?.to_string());
    }

    Some(DeviceTemplate {
        name,
        kind,
        model,
        terminals,
        property_layers,
        params,
    })
}

pub fn read_devtab_file<P: AsRef<Path>>(path: P) -> io::Result<Vec<DeviceTemplate>> {
    let path = path.as_ref();
    let text = std::fs::read_to_string(path)?;
    let (templates, diags) = parse_devtab(&text);
    openrdson_core::diag::log_diags(&path.display().to_string(), &diags);
    Ok(templates)
}

/// Parse instance lines from a CCI/LVS SPICE netlist (`*.spi`).
///
/// An instance line looks like:
/// ```text
/// XX0 1 3 2 1 device_model w=1e-05 weff=0.001 nf=100 nx=10 ny=10 m=1 mlay=1 avnx=10 $X=23680 $Y=10280 $D=202
/// ```
///
/// Returns the instances plus diagnostics for instance lines that could not be
/// resolved (e.g. no model name found).
pub fn parse_spi_instances(text: &str) -> (Vec<SpiInstance>, Vec<Diag>) {
    let mut out = Vec::new();
    let mut diags = Vec::new();
    for (i, raw) in text.lines().enumerate() {
        let line_no = i + 1;
        let line = raw.trim();
        if !line.starts_with('X') || line.starts_with("XX_") {
            continue;
        }
        let t: Vec<&str> = line.split_whitespace().collect();
        if t.len() < 3 {
            diags.push(Diag::at(
                line_no,
                format!("instance line has too few fields: \"{line}\""),
            ));
            continue;
        }
        // Find the model: first token after the name that is not a node/number.
        let mut model_idx = None;
        for (i, tok) in t.iter().enumerate().skip(1) {
            if tok.contains('=') {
                break;
            }
            if tok.parse::<f64>().is_err() {
                model_idx = Some(i);
                break;
            }
        }
        let Some(mi) = model_idx else {
            diags.push(Diag::at(
                line_no,
                format!("could not find a device model on instance line: \"{line}\""),
            ));
            continue;
        };
        let nodes = t[1..mi].iter().map(|s| s.to_string()).collect();
        let model = t[mi].to_string();
        let mut params = std::collections::BTreeMap::new();
        let mut x = None;
        let mut y = None;
        for tok in &t[mi + 1..] {
            let Some((k, v)) = tok.split_once('=') else {
                continue;
            };
            match k {
                "$X" => x = v.parse::<f64>().ok(),
                "$Y" => y = v.parse::<f64>().ok(),
                _ if k.starts_with('$') => {}
                _ => {
                    if let Ok(val) = v.parse::<f64>() {
                        params.insert(k.to_string(), val);
                    }
                }
            }
        }
        out.push(SpiInstance {
            name: t[0].to_string(),
            nodes,
            model,
            params,
            x,
            y,
        });
    }
    (out, diags)
}

pub fn read_spi_instances<P: AsRef<Path>>(path: P) -> io::Result<Vec<SpiInstance>> {
    let path = path.as_ref();
    let text = std::fs::read_to_string(path)?;
    let (instances, diags) = parse_spi_instances(&text);
    openrdson_core::diag::log_diags(&path.display().to_string(), &diags);
    Ok(instances)
}

/// Parse a CCI `*.lnn` layout net-name table (`<node id> <net name>`), plus
/// diagnostics for lines whose leading id is not a number.
pub fn parse_net_names(text: &str) -> (std::collections::BTreeMap<i64, String>, Vec<Diag>) {
    let mut out = std::collections::BTreeMap::new();
    let mut diags = Vec::new();
    for (i, raw) in text.lines().enumerate() {
        let line_no = i + 1;
        let line = raw.trim();
        if line.is_empty()
            || line.starts_with('#')
            || line.starts_with('%')
            || line.starts_with('*')
        {
            continue;
        }
        let t: Vec<&str> = line.split_whitespace().collect();
        if t.len() < 2 {
            continue;
        }
        match t[0].parse::<i64>() {
            Ok(id) => {
                out.insert(id, t[1].to_string());
            }
            Err(_) => diags.push(Diag::at(
                line_no,
                format!("expected \"<node id> <net name>\", got \"{line}\""),
            )),
        }
    }
    (out, diags)
}

pub fn read_net_names<P: AsRef<Path>>(
    path: P,
) -> io::Result<std::collections::BTreeMap<i64, String>> {
    let path = path.as_ref();
    let text = std::fs::read_to_string(path)?;
    let (names, diags) = parse_net_names(&text);
    openrdson_core::diag::log_diags(&path.display().to_string(), &diags);
    Ok(names)
}

/// A device template from the CCI netlist `.DEVTMPLT` dictionary: the model,
/// seed layer, and the ordered `(terminal name, layer)` list.
#[derive(Debug, Clone, PartialEq)]
pub struct CciDeviceTemplate {
    pub model: String,
    pub seed: String,
    pub terminals: Vec<(String, String)>,
}

/// Parse the `*.DEVTMPLT` lines of a CCI netlist (`*.pin_xy_spi`), plus
/// diagnostics for malformed template lines.
pub fn parse_device_templates(
    text: &str,
) -> (std::collections::BTreeMap<i64, CciDeviceTemplate>, Vec<Diag>) {
    let mut out = std::collections::BTreeMap::new();
    let mut diags = Vec::new();
    for (i, raw) in text.lines().enumerate() {
        let line_no = i + 1;
        let Some(rest) = raw.trim().strip_prefix("*.DEVTMPLT ") else {
            continue;
        };
        let t: Vec<&str> = rest.split_whitespace().collect();
        if t.len() < 3 {
            diags.push(Diag::at(
                line_no,
                format!("expected \"*.DEVTMPLT <id> <model> <seed> [term(layer) ...]\", got \"{rest}\""),
            ));
            continue;
        }
        let Ok(num) = t[0].parse::<i64>() else {
            diags.push(Diag::at(
                line_no,
                format!("non-numeric .DEVTMPLT id \"{}\"", t[0]),
            ));
            continue;
        };
        let model = t[1].trim_end_matches("()").to_string();
        let seed = t[2].to_string();
        let terminals = t[3..]
            .iter()
            .filter_map(|tok| {
                let (layer, term) = tok.split_once('(')?;
                Some((term.trim_end_matches(')').to_string(), layer.to_string()))
            })
            .collect();
        out.insert(
            num,
            CciDeviceTemplate {
                model,
                seed,
                terminals,
            },
        );
    }
    (out, diags)
}

pub fn read_cci_device_templates<P: AsRef<Path>>(
    path: P,
) -> io::Result<std::collections::BTreeMap<i64, CciDeviceTemplate>> {
    let path = path.as_ref();
    let text = std::fs::read_to_string(path)?;
    let (templates, diags) = parse_device_templates(&text);
    openrdson_core::diag::log_diags(&path.display().to_string(), &diags);
    Ok(templates)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_lnn_net_names() {
        let (m, diags) = parse_net_names("# header\n% cell 0\n1 D\n2 S\n3 G\n");
        assert!(diags.is_empty(), "{diags:?}");
        assert_eq!(m.get(&1).map(String::as_str), Some("D"));
        assert_eq!(m.get(&3).map(String::as_str), Some("G"));
    }

    #[test]
    fn flags_a_lnn_line_without_a_numeric_id() {
        let (_m, diags) = parse_net_names("1 D\nnotanid S\n");
        assert!(
            diags
                .iter()
                .any(|d| d.message.contains("expected \"<node id> <net name>\"")),
            "{diags:?}"
        );
        assert_eq!(diags[0].line, 2);
    }

    #[test]
    fn parses_devtmplt_terminal_order() {
        let text = "*.DEVTMPLT 202 model_a() seed_a net_d(d) net_g(g) net_s(s) net_sub(sub)\n";
        let (m, diags) = parse_device_templates(text);
        assert!(diags.is_empty(), "{diags:?}");
        let t = m.get(&202).unwrap();
        assert_eq!(t.model, "model_a");
        assert_eq!(t.seed, "seed_a");
        assert_eq!(t.terminals[0], ("d".to_string(), "net_d".to_string()));
        assert_eq!(t.terminals[1].0, "g");
        assert_eq!(t.terminals[2], ("s".to_string(), "net_s".to_string()));
    }

    #[test]
    fn parses_ports_and_scales_to_meters() {
        // <name> <id> <direction> <x> <y> <net>
        let (ports, diags) = parse_ports("D 1 inout 1000 2000 net_M1\nS 2 inout 0 2000 net_M2\n", 1e-9);
        assert!(diags.is_empty(), "{diags:?}");
        assert_eq!(ports.len(), 2);
        assert_eq!(ports[0].name, "D");
        assert_eq!(ports[0].x, 1000.0 * 1e-9);
        assert_eq!(ports[0].net, "net_M1");
        assert_eq!(ports[1].y, 2000.0 * 1e-9);
    }

    #[test]
    fn flags_ports_with_too_few_fields_and_bad_coordinates() {
        let (_p, diags) = parse_ports("D 1 inout 1000\n", 1e-9);
        assert!(
            diags.iter().any(|d| d.message.contains("expected 6 fields")),
            "{diags:?}"
        );
        let (_p, diags) = parse_ports("D 1 inout abc 2000 net_M1\n", 1e-9);
        assert!(
            diags
                .iter()
                .any(|d| d.message.contains("non-numeric coordinate")),
            "{diags:?}"
        );
    }

    #[test]
    fn flags_a_malformed_devtab_entry() {
        // A "Device Entry" block whose body is missing the required fields.
        let (_t, diags) = parse_devtab("Device Entry 1\nAug 12 2026\n\n");
        assert!(
            diags.iter().any(|d| d.message.contains("malformed 'Device Entry'")),
            "{diags:?}"
        );
    }

    #[test]
    fn parses_a_spi_instance_and_flags_a_modelless_line() {
        let (inst, diags) = parse_spi_instances(
            "X0 1 3 2 dev w=1e-6 nf=10 $X=100 $Y=200\nX1 1 3 2\n",
        );
        assert_eq!(inst.len(), 1);
        assert_eq!(inst[0].name, "X0");
        assert_eq!(inst[0].model, "dev");
        assert_eq!(inst[0].params.get("nf").copied(), Some(10.0));
        assert_eq!(inst[0].x, Some(100.0));
        assert!(
            diags
                .iter()
                .any(|d| d.message.contains("could not find a device model")),
            "{diags:?}"
        );
        assert_eq!(diags[0].line, 2);
    }
}

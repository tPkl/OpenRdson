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
/// to meters with `db_unit_meters`.
pub fn parse_ports(text: &str, db_unit_meters: f64) -> Vec<Port> {
    let mut ports = Vec::new();
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let t: Vec<&str> = line.split_whitespace().collect();
        // <name> <id> <direction> <x> <y> <net>
        if t.len() < 6 {
            continue;
        }
        let (Ok(x), Ok(y)) = (t[3].parse::<f64>(), t[4].parse::<f64>()) else {
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
    ports
}

pub fn read_ports_file<P: AsRef<Path>>(path: P, db_unit_meters: f64) -> io::Result<Vec<Port>> {
    let text = std::fs::read_to_string(path)?;
    Ok(parse_ports(&text, db_unit_meters))
}

/// Parse a CCI `devtab` device table into device *definitions* (templates).
pub fn parse_devtab(text: &str) -> Vec<DeviceTemplate> {
    let lines: Vec<&str> = text.lines().collect();
    let starts: Vec<usize> = lines
        .iter()
        .enumerate()
        .filter(|(_, l)| l.trim_start().starts_with("Device Entry"))
        .map(|(i, _)| i)
        .collect();

    let mut templates = Vec::new();
    for (n, &start) in starts.iter().enumerate() {
        let end = starts.get(n + 1).copied().unwrap_or(lines.len());
        // Skip the "Device Entry N" line and the following timestamp line.
        let body = &lines[(start + 2).min(end)..end];
        if let Some(t) = parse_device_entry(body) {
            templates.push(t);
        }
    }
    templates
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
    let text = std::fs::read_to_string(path)?;
    Ok(parse_devtab(&text))
}

/// Parse instance lines from a CCI/LVS SPICE netlist (`*.spi`).
///
/// An instance line looks like:
/// ```text
/// XX0 1 3 2 1 device_model w=1e-05 weff=0.001 nf=100 nx=10 ny=10 m=1 mlay=1 avnx=10 $X=23680 $Y=10280 $D=202
/// ```
pub fn parse_spi_instances(text: &str) -> Vec<SpiInstance> {
    let mut out = Vec::new();
    for line in text.lines() {
        let line = line.trim();
        if !line.starts_with('X') || line.starts_with("XX_") {
            continue;
        }
        let t: Vec<&str> = line.split_whitespace().collect();
        if t.len() < 3 {
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
        let Some(mi) = model_idx else { continue };
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
    out
}

pub fn read_spi_instances<P: AsRef<Path>>(path: P) -> io::Result<Vec<SpiInstance>> {
    let text = std::fs::read_to_string(path)?;
    Ok(parse_spi_instances(&text))
}

/// Parse a CCI `*.lnn` layout net-name table (`<node id> <net name>`).
pub fn parse_net_names(text: &str) -> std::collections::BTreeMap<i64, String> {
    let mut out = std::collections::BTreeMap::new();
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty()
            || line.starts_with('#')
            || line.starts_with('%')
            || line.starts_with('*')
        {
            continue;
        }
        let t: Vec<&str> = line.split_whitespace().collect();
        if t.len() >= 2 {
            if let Ok(id) = t[0].parse::<i64>() {
                out.insert(id, t[1].to_string());
            }
        }
    }
    out
}

pub fn read_net_names<P: AsRef<Path>>(
    path: P,
) -> io::Result<std::collections::BTreeMap<i64, String>> {
    Ok(parse_net_names(&std::fs::read_to_string(path)?))
}

/// A device template from the CCI netlist `.DEVTMPLT` dictionary: the model,
/// seed layer, and the ordered `(terminal name, layer)` list.
#[derive(Debug, Clone, PartialEq)]
pub struct CciDeviceTemplate {
    pub model: String,
    pub seed: String,
    pub terminals: Vec<(String, String)>,
}

/// Parse the `*.DEVTMPLT` lines of a CCI netlist (`*.pin_xy_spi`).
pub fn parse_device_templates(text: &str) -> std::collections::BTreeMap<i64, CciDeviceTemplate> {
    let mut out = std::collections::BTreeMap::new();
    for line in text.lines() {
        let Some(rest) = line.trim().strip_prefix("*.DEVTMPLT ") else {
            continue;
        };
        let t: Vec<&str> = rest.split_whitespace().collect();
        if t.len() < 3 {
            continue;
        }
        let Ok(num) = t[0].parse::<i64>() else { continue };
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
    out
}

pub fn read_cci_device_templates<P: AsRef<Path>>(
    path: P,
) -> io::Result<std::collections::BTreeMap<i64, CciDeviceTemplate>> {
    Ok(parse_device_templates(&std::fs::read_to_string(path)?))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_lnn_net_names() {
        let m = parse_net_names("# header\n% cell 0\n1 D\n2 S\n3 G\n");
        assert_eq!(m.get(&1).map(String::as_str), Some("D"));
        assert_eq!(m.get(&3).map(String::as_str), Some("G"));
    }

    #[test]
    fn parses_devtmplt_terminal_order() {
        let text = "*.DEVTMPLT 202 model_a() seed_a net_d(d) net_g(g) net_s(s) net_sub(sub)\n";
        let m = parse_device_templates(text);
        let t = m.get(&202).unwrap();
        assert_eq!(t.model, "model_a");
        assert_eq!(t.seed, "seed_a");
        assert_eq!(t.terminals[0], ("d".to_string(), "net_d".to_string()));
        assert_eq!(t.terminals[1].0, "g");
        assert_eq!(t.terminals[2], ("s".to_string(), "net_s".to_string()));
    }
}

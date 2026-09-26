//! GDSII writer and KLayout visualization helpers.
//!
//! Writes computed quantities (potential, current density, resistance, mesh
//! regions) as GDSII boundaries binned onto layers, plus a KLayout
//! `<gds>.lyp` layer-properties file that colours the bins into a colormap.
//! Opening the `.gds` in KLayout auto-loads the `.lyp`.

use std::io;
use std::path::Path;

/// A boundary element to write.
#[derive(Debug, Clone)]
pub struct GdsBoundary {
    pub layer: i16,
    pub datatype: i16,
    /// Closed ring in database units (closing point added automatically).
    pub points: Vec<(i32, i32)>,
}

/// A visualization item: a footprint (database units) and its scalar value.
#[derive(Debug, Clone)]
pub struct VizItem {
    pub points: Vec<(i32, i32)>,
    pub value: f64,
}

fn push_record(buf: &mut Vec<u8>, rtype: u8, dtype: u8, payload: &[u8]) {
    let len = (4 + payload.len()) as u16;
    buf.extend_from_slice(&len.to_be_bytes());
    buf.push(rtype);
    buf.push(dtype);
    buf.extend_from_slice(payload);
}

fn push_string(buf: &mut Vec<u8>, rtype: u8, s: &str) {
    let mut b = s.as_bytes().to_vec();
    if b.len() % 2 == 1 {
        b.push(0);
    }
    push_record(buf, rtype, 0x06, &b);
}

/// Write a GDSII library with one or more cells of boundary elements.
pub fn write_gds_file<P: AsRef<Path>>(
    path: P,
    lib_name: &str,
    db_unit_meters: f64,
    cells: &[(String, Vec<GdsBoundary>)],
) -> io::Result<()> {
    let mut buf: Vec<u8> = Vec::new();
    push_record(&mut buf, 0x00, 0x02, &600i16.to_be_bytes()); // HEADER v600
    push_record(&mut buf, 0x01, 0x02, &[0u8; 24]); // BGNLIB
    push_string(&mut buf, 0x02, lib_name);
    // GDSII UNITS = (database unit in user units, database unit in meters).
    // Use a 1 µm user unit (the conventional choice) so viewers such as KLayout
    // display the layout in microns rather than interpreting user unit = 1 m.
    let mut units = Vec::with_capacity(16);
    units.extend_from_slice(&(db_unit_meters / 1e-6).to_be_bytes());
    units.extend_from_slice(&db_unit_meters.to_be_bytes());
    push_record(&mut buf, 0x03, 0x05, &units); // UNITS

    for (name, elements) in cells {
        push_record(&mut buf, 0x05, 0x02, &[0u8; 24]); // BGNSTR
        push_string(&mut buf, 0x06, name);
        for e in elements {
            if e.points.len() < 3 {
                continue;
            }
            push_record(&mut buf, 0x08, 0x00, &[]); // BOUNDARY
            push_record(&mut buf, 0x0D, 0x02, &e.layer.to_be_bytes()); // LAYER
            push_record(&mut buf, 0x0E, 0x02, &e.datatype.to_be_bytes()); // DATATYPE
            let mut xy = Vec::with_capacity((e.points.len() + 1) * 8);
            for &(x, y) in &e.points {
                xy.extend_from_slice(&x.to_be_bytes());
                xy.extend_from_slice(&y.to_be_bytes());
            }
            let (x0, y0) = e.points[0];
            xy.extend_from_slice(&x0.to_be_bytes());
            xy.extend_from_slice(&y0.to_be_bytes());
            push_record(&mut buf, 0x10, 0x03, &xy); // XY
            push_record(&mut buf, 0x11, 0x00, &[]); // ENDEL
        }
        push_record(&mut buf, 0x07, 0x00, &[]); // ENDSTR
    }
    push_record(&mut buf, 0x04, 0x00, &[]); // ENDLIB
    std::fs::write(path, buf)
}

/// A blue→cyan→green→yellow→red colour for `t` in `[0, 1]`.
fn jet(t: f64) -> (u8, u8, u8) {
    let t = t.clamp(0.0, 1.0);
    let (r, g, b) = if t < 0.25 {
        (0.0, 4.0 * t, 1.0)
    } else if t < 0.5 {
        (0.0, 1.0, 1.0 - 4.0 * (t - 0.25))
    } else if t < 0.75 {
        (4.0 * (t - 0.5), 1.0, 0.0)
    } else {
        (1.0, 1.0 - 4.0 * (t - 0.75), 0.0)
    };
    ((r * 255.0) as u8, (g * 255.0) as u8, (b * 255.0) as u8)
}

/// Write a value map as a GDSII colormap: each item's footprint goes on the
/// layer `layer_base + bin`, where `bin` is the value's quantile in
/// `[vmin, vmax]`. Also writes `<path>.lyp` for KLayout.
///
/// Returns the `(layer, value)` bin centres for reference.
pub fn write_colormap<P: AsRef<Path>>(
    path: P,
    lib_name: &str,
    cell_name: &str,
    db_unit_meters: f64,
    items: &[VizItem],
    n_bins: usize,
    layer_base: i16,
    vmin: f64,
    vmax: f64,
) -> io::Result<Vec<(i16, f64)>> {
    let n_bins = n_bins.max(1);
    let span = (vmax - vmin).max(1e-300);
    let mut boundaries = Vec::with_capacity(items.len());
    for it in items {
        if it.points.len() < 3 || !it.value.is_finite() {
            continue;
        }
        let t = ((it.value - vmin) / span).clamp(0.0, 1.0);
        let bin = ((t * (n_bins as f64 - 1.0)).round() as usize).min(n_bins - 1);
        boundaries.push(GdsBoundary {
            layer: layer_base + bin as i16,
            datatype: 0,
            points: it.points.clone(),
        });
    }
    write_gds_file(&path, lib_name, db_unit_meters, &[(cell_name.to_string(), boundaries)])?;

    // KLayout auto-loads a layer-properties file with the same base name
    // (e.g. potential.gds -> potential.lyp).
    let mut lyp_path = path.as_ref().to_path_buf();
    lyp_path.set_extension("lyp");
    let mut bins = Vec::with_capacity(n_bins);
    for b in 0..n_bins {
        let t = if n_bins <= 1 { 0.0 } else { b as f64 / (n_bins - 1) as f64 };
        bins.push((layer_base + b as i16, vmin + t * span));
    }
    write_lyp(&lyp_path, &bins, vmin, vmax)?;
    Ok(bins)
}

fn lyp_entry(s: &mut String, layer: i16, name: &str, color: (u8, u8, u8), indent: &str) {
    let (r, g, b) = color;
    let color = format!("#{r:02X}{g:02X}{b:02X}");
    s.push_str(&format!("{indent}<properties>\n"));
    s.push_str(&format!("{indent} <frame-color>{color}</frame-color>\n"));
    s.push_str(&format!("{indent} <fill-color>{color}</fill-color>\n"));
    s.push_str(&format!("{indent} <frame-brightness>0</frame-brightness>\n"));
    s.push_str(&format!("{indent} <fill-brightness>0</fill-brightness>\n"));
    s.push_str(&format!("{indent} <dither-pattern>I0</dither-pattern>\n"));
    s.push_str(&format!("{indent} <valid>true</valid>\n"));
    s.push_str(&format!("{indent} <visible>true</visible>\n"));
    s.push_str(&format!("{indent} <transparent>false</transparent>\n"));
    s.push_str(&format!("{indent} <width>1</width>\n"));
    s.push_str(&format!("{indent} <marked>false</marked>\n"));
    s.push_str(&format!("{indent} <animation>0</animation>\n"));
    s.push_str(&format!("{indent} <name>{name}</name>\n"));
    // KLayout GDS layer source format is `<layer>/<datatype>@<view>`; a
    // `layer/...` prefix is parsed as a *name* source and matches nothing.
    s.push_str(&format!("{indent} <source>{layer}/0@*</source>\n"));
    s.push_str(&format!("{indent}</properties>\n"));
}

/// A child of a KLayout group: a `<group-members>` element holding the entry's
/// properties directly. KLayout's format nests *each* child as its own
/// `<group-members>` block under the group `<properties>` — a single wrapping
/// `<group-members>` containing `<properties>` children is parsed wrongly and
/// collapses every group into the last one.
fn lyp_group_child(s: &mut String, layer: i16, name: &str, color: (u8, u8, u8), indent: &str) {
    let (r, g, b) = color;
    let color = format!("#{r:02X}{g:02X}{b:02X}");
    s.push_str(&format!("{indent}<group-members>\n"));
    s.push_str(&format!("{indent} <source>{layer}/0</source>\n"));
    s.push_str(&format!("{indent} <name>{name}</name>\n"));
    s.push_str(&format!("{indent} <frame-color>{color}</frame-color>\n"));
    s.push_str(&format!("{indent} <fill-color>{color}</fill-color>\n"));
    s.push_str(&format!("{indent} <width>1</width>\n"));
    s.push_str(&format!("{indent} <frame-brightness>0</frame-brightness>\n"));
    s.push_str(&format!("{indent} <fill-brightness>0</fill-brightness>\n"));
    s.push_str(&format!("{indent} <dither-pattern>I0</dither-pattern>\n"));
    s.push_str(&format!("{indent} <visible>true</visible>\n"));
    s.push_str(&format!("{indent} <transparent>false</transparent>\n"));
    s.push_str(&format!("{indent} <marked>false</marked>\n"));
    s.push_str(&format!("{indent} <animation>0</animation>\n"));
    s.push_str(&format!("{indent}</group-members>\n"));
}

/// Write a single GDS containing a **grouped** colormap: each physical layer
/// becomes a KLayout group whose children are the value bins.
///
/// `layers` is `(physical_layer_name, cells)`; `n_bins` sets the colormap depth.
/// Bin `b` of physical layer `p` is written on layer
/// `layer_base + p_index*stride + b`, datatype 0. One `<name>.lyp` is written
/// with a `<properties>` group per physical layer whose children are each a
/// `<group-members>` block (KLayout's canonical grouped format).
pub fn write_grouped_colormap<P: AsRef<Path>>(
    path: P,
    lib_name: &str,
    cell_name: &str,
    db_unit_meters: f64,
    layers: &[(String, Vec<VizItem>)],
    n_bins: usize,
    layer_base: i16,
    stride: i16,
) -> io::Result<()> {
    let n_bins = n_bins.max(1);
    let mut boundaries = Vec::new();
    let mut groups: Vec<(String, Vec<(i16, f64)>)> = Vec::new();

    for (p_idx, (phys, items)) in layers.iter().enumerate() {
        let vals: Vec<f64> = items.iter().map(|i| i.value).filter(|v| v.is_finite()).collect();
        if vals.is_empty() {
            continue;
        }
        let vmin = vals.iter().cloned().fold(f64::INFINITY, f64::min);
        let vmax = vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let span = (vmax - vmin).max(1e-300);
        let base = layer_base + p_idx as i16 * stride;
        let mut bins = Vec::with_capacity(n_bins);
        for it in items {
            if it.points.len() < 3 || !it.value.is_finite() {
                continue;
            }
            let t = ((it.value - vmin) / span).clamp(0.0, 1.0);
            let bin = ((t * (n_bins as f64 - 1.0)).round() as usize).min(n_bins - 1);
            boundaries.push(GdsBoundary {
                layer: base + bin as i16,
                datatype: 0,
                points: it.points.clone(),
            });
        }
        for b in 0..n_bins {
            let t = if n_bins <= 1 { 0.0 } else { b as f64 / (n_bins - 1) as f64 };
            bins.push((base + b as i16, vmin + t * span));
        }
        groups.push((phys.clone(), bins));
    }

    write_gds_file(
        &path,
        lib_name,
        db_unit_meters,
        &[(cell_name.to_string(), boundaries)],
    )?;

    let mut lyp_path = path.as_ref().to_path_buf();
    lyp_path.set_extension("lyp");
    let mut s = String::from("<?xml version=\"1.0\" encoding=\"utf-8\"?>\n<layer-properties>\n");
    for (phys, bins) in &groups {
        s.push_str(" <properties>\n");
        s.push_str(&format!("  <name>{phys}</name>\n"));
        s.push_str("  <expanded>true</expanded>\n");
        for &(layer, value) in bins {
            let t = if n_bins <= 1 {
                0.0
            } else {
                (value - bins[0].1) / (bins[bins.len() - 1].1 - bins[0].1).max(1e-300)
            };
            lyp_group_child(
                &mut s,
                layer,
                &format!("{value:.4e}"),
                jet(t.clamp(0.0, 1.0)),
                "  ",
            );
        }
        s.push_str(" </properties>\n");
    }
    s.push_str("</layer-properties>\n");
    std::fs::write(lyp_path, s)
}

/// Write a KLayout `.lyp` with one named entry per `(layer, name)`.
pub fn write_named_lyp<P: AsRef<Path>>(path: P, entries: &[(i16, String)]) -> io::Result<()> {
    let mut s = String::from("<?xml version=\"1.0\" encoding=\"utf-8\"?>\n<layer-properties>\n");
    for (i, (layer, name)) in entries.iter().enumerate() {
        let t = (i % 16) as f64 / 15.0;
        lyp_entry(&mut s, *layer, name, jet(t), " ");
    }
    s.push_str("</layer-properties>\n");
    std::fs::write(path, s)
}

/// Write a KLayout layer-properties (`.lyp`) file colouring each bin.
pub fn write_lyp<P: AsRef<Path>>(
    path: P,
    bins: &[(i16, f64)],
    vmin: f64,
    vmax: f64,
) -> io::Result<()> {
    let span = (vmax - vmin).max(1e-300);
    let mut s = String::from("<?xml version=\"1.0\" encoding=\"utf-8\"?>\n<layer-properties>\n");
    for &(layer, value) in bins {
        let t = ((value - vmin) / span).clamp(0.0, 1.0);
        lyp_entry(&mut s, layer, &format!("{value:.4e}"), jet(t), " ");
    }
    s.push_str("</layer-properties>\n");
    std::fs::write(path, s)
}

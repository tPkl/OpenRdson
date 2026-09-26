//! A minimal, dependency-free GDSII / annotated-GDSII (AGF) reader.
//!
//! It parses the records needed to reconstruct layout polygons and their LVS
//! element properties: `HEADER, BGNLIB, LIBNAME, UNITS, BGNSTR, STRNAME,
//! BOUNDARY, PATH, SREF, AREF, TEXT, LAYER, DATATYPE, XY, ENDEL, ENDSTR,
//! PROPATTR, PROPVALUE, ENDLIB`.
//!
//! Coordinates are kept as raw database units (`i32`) here; conversion to
//! meters happens when building [`openrdson_core::layout::LayoutIR`].

use std::io::{self, Read};
use std::path::Path;

// Record types.
const R_HEADER: u8 = 0x00;
const R_BGNLIB: u8 = 0x01;
const R_LIBNAME: u8 = 0x02;
const R_UNITS: u8 = 0x03;
const R_ENDLIB: u8 = 0x04;
const R_BGNSTR: u8 = 0x05;
const R_STRNAME: u8 = 0x06;
const R_ENDSTR: u8 = 0x07;
const R_BOUNDARY: u8 = 0x08;
const R_PATH: u8 = 0x09;
const R_SREF: u8 = 0x0A;
const R_AREF: u8 = 0x0B;
const R_TEXT: u8 = 0x0C;
const R_LAYER: u8 = 0x0D;
const R_DATATYPE: u8 = 0x0E;
const R_XY: u8 = 0x10;
const R_ENDEL: u8 = 0x11;
const R_TEXTTYPE: u8 = 0x16;
const R_PROPATTR: u8 = 0x2B;
const R_PROPVALUE: u8 = 0x2C;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ElementKind {
    Boundary,
    Path,
    Text,
    SRef,
    ARef,
    Other,
}

#[derive(Debug, Clone, PartialEq)]
pub struct GdsElement {
    pub kind: ElementKind,
    pub layer: i16,
    pub datatype: i16,
    /// Raw database-unit coordinate pairs.
    pub xy: Vec<(i32, i32)>,
    /// `(PROPATTR, PROPVALUE)` pairs attached to this element.
    pub props: Vec<(i16, String)>,
}

impl Default for GdsElement {
    fn default() -> Self {
        Self {
            kind: ElementKind::Other,
            layer: 0,
            datatype: 0,
            xy: Vec::new(),
            props: Vec::new(),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct GdsStruct {
    pub name: String,
    pub elements: Vec<GdsElement>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct GdsLibrary {
    pub name: String,
    /// Size of one database unit in meters (`UNITS[1]`).
    pub db_unit_meters: f64,
    /// Size of the user unit in meters (`UNITS[1] / UNITS[0]`).
    pub user_unit_meters: f64,
    pub structs: Vec<GdsStruct>,
}

fn err(msg: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg.into())
}

fn read_f64s(payload: &[u8], off: usize) -> Vec<f64> {
    payload[off..]
        .chunks_exact(8)
        .map(|c| f64::from_be_bytes([c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7]]))
        .collect()
}

fn read_i16(payload: &[u8], off: usize) -> Option<i16> {
    if payload.len() >= off + 2 {
        Some(i16::from_be_bytes([payload[off], payload[off + 1]]))
    } else {
        None
    }
}

fn read_string(payload: &[u8], off: usize) -> String {
    let bytes = &payload[off..];
    let end = bytes.iter().position(|&b| b == 0).unwrap_or(bytes.len());
    String::from_utf8_lossy(&bytes[..end]).trim_end().to_string()
}

/// Read a GDSII/AGF file from disk.
pub fn read_gds_file<P: AsRef<Path>>(path: P) -> io::Result<GdsLibrary> {
    let mut f = std::fs::File::open(path)?;
    let mut buf = Vec::new();
    f.read_to_end(&mut buf)?;
    read_gds(&buf)
}

/// Parse an in-memory GDSII/AGF byte stream.
pub fn read_gds(buf: &[u8]) -> io::Result<GdsLibrary> {
    let mut lib = GdsLibrary {
        name: String::new(),
        db_unit_meters: 0.0,
        user_unit_meters: 0.0,
        structs: Vec::new(),
    };

    let mut cur_struct: Option<GdsStruct> = None;
    let mut cur_elem: Option<GdsElement> = None;
    let mut i = 0usize;
    let mut seen_endlib = false;

    while i + 4 <= buf.len() {
        let len = u16::from_be_bytes([buf[i], buf[i + 1]]) as usize;
        if len < 4 {
            return Err(err(format!("record at offset {i} has invalid length {len}")));
        }
        if i + len > buf.len() {
            return Err(err(format!(
                "record at offset {i} claims length {len} but only {} bytes remain",
                buf.len() - i
            )));
        }
        let rtype = buf[i + 2];
        let payload = &buf[i + 4..i + len];

        match rtype {
            R_HEADER => {}
            R_BGNLIB => {}
            R_LIBNAME => lib.name = read_string(payload, 0),
            R_UNITS => {
                let vals = read_f64s(payload, 0);
                if vals.len() >= 2 {
                    let db = vals[1];
                    lib.db_unit_meters = db;
                    lib.user_unit_meters = if vals[0] != 0.0 { db / vals[0] } else { db };
                } else {
                    return Err(err("UNITS record with fewer than 2 doubles"));
                }
            }
            R_BGNSTR => cur_struct = Some(GdsStruct { name: String::new(), elements: Vec::new() }),
            R_STRNAME => {
                if let Some(s) = cur_struct.as_mut() {
                    s.name = read_string(payload, 0);
                }
            }
            R_ENDSTR => {
                if let Some(s) = cur_struct.take() {
                    lib.structs.push(s);
                }
            }
            R_BOUNDARY | R_PATH | R_TEXT | R_SREF | R_AREF => {
                let kind = match rtype {
                    R_BOUNDARY => ElementKind::Boundary,
                    R_PATH => ElementKind::Path,
                    R_TEXT => ElementKind::Text,
                    R_SREF => ElementKind::SRef,
                    _ => ElementKind::ARef,
                };
                let mut e = GdsElement::default();
                e.kind = kind;
                cur_elem = Some(e);
            }
            R_LAYER => {
                if let (Some(e), Some(v)) = (cur_elem.as_mut(), read_i16(payload, 0)) {
                    e.layer = v;
                }
            }
            R_DATATYPE | R_TEXTTYPE => {
                if let (Some(e), Some(v)) = (cur_elem.as_mut(), read_i16(payload, 0)) {
                    e.datatype = v;
                }
            }
            R_XY => {
                if let Some(e) = cur_elem.as_mut() {
                    let ints: Vec<i32> = payload
                        .chunks_exact(4)
                        .map(|c| i32::from_be_bytes([c[0], c[1], c[2], c[3]]))
                        .collect();
                    e.xy = ints.chunks_exact(2).map(|p| (p[0], p[1])).collect();
                }
            }
            R_PROPATTR => {
                if let (Some(e), Some(v)) = (cur_elem.as_mut(), read_i16(payload, 0)) {
                    e.props.push((v, String::new()));
                }
            }
            R_PROPVALUE => {
                if let Some(e) = cur_elem.as_mut() {
                    let val = read_string(payload, 0);
                    match e.props.last_mut() {
                        Some((_, s)) if s.is_empty() => *s = val,
                        _ => e.props.push((0, val)),
                    }
                }
            }
            R_ENDEL => {
                if let Some(e) = cur_elem.take() {
                    if let Some(s) = cur_struct.as_mut() {
                        s.elements.push(e);
                    }
                }
            }
            R_ENDLIB => {
                seen_endlib = true;
                break;
            }
            _ => {
                // SNAME, COLROW, WIDTH, STRING, etc. are not needed to build
                // conducting polygons; ignored deliberately.
            }
        }
        i += len;
    }

    if cur_elem.is_some() || cur_struct.is_some() {
        return Err(err("file ended inside an element/structure (truncated?)"));
    }
    if !seen_endlib {
        return Err(err("missing ENDLIB record"));
    }
    Ok(lib)
}

impl GdsLibrary {
    pub fn element_count(&self) -> usize {
        self.structs.iter().map(|s| s.elements.len()).sum()
    }
}

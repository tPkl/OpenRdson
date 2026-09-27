//! Device IR: templates (definitions from `devtab`) and recognized instances.

use crate::geometry::Bbox;
use std::collections::BTreeMap;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DeviceKind {
    Mos,
    Resistor,
    Capacitor,
    Diode,
    Pnp,
    Npn,
    User,
    Unknown,
}

impl DeviceKind {
    pub fn parse(s: &str) -> Self {
        match s.trim().to_ascii_uppercase().as_str() {
            "MOS" => DeviceKind::Mos,
            "RESISTOR" | "RES" => DeviceKind::Resistor,
            "CAPACITOR" | "CAP" => DeviceKind::Capacitor,
            "DIODE" => DeviceKind::Diode,
            "PNP" => DeviceKind::Pnp,
            "NPN" => DeviceKind::Npn,
            "USER" => DeviceKind::User,
            _ => DeviceKind::Unknown,
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            DeviceKind::Mos => "MOS",
            DeviceKind::Resistor => "RESISTOR",
            DeviceKind::Capacitor => "CAPACITOR",
            DeviceKind::Diode => "DIODE",
            DeviceKind::Pnp => "PNP",
            DeviceKind::Npn => "NPN",
            DeviceKind::User => "USER",
            DeviceKind::Unknown => "UNKNOWN",
        }
    }
}

/// A device *definition* from the CCI device table (`devtab`).
#[derive(Debug, Clone, PartialEq)]
pub struct DeviceTemplate {
    /// Seed layer / rule name.
    pub name: String,
    pub kind: DeviceKind,
    /// Subcircuit / model name.
    pub model: Option<String>,
    /// Terminal name -> layer name, e.g. `d -> net_nwell_O_600`.
    pub terminals: Vec<(String, String)>,
    /// Property layers carrying instance parameters.
    pub property_layers: Vec<String>,
    /// Instance parameter names, e.g. `w nx ns ny m mlay avnx`.
    pub params: Vec<String>,
}

/// A recognized device *instance* (produced by `device-recognition`).
#[derive(Debug, Clone, PartialEq)]
pub struct DeviceInstance {
    pub name: String,
    pub kind: DeviceKind,
    pub model: Option<String>,
    /// Name of the template (seed rule) that produced this instance.
    pub template: String,
    /// Terminal name -> logical layer (from the template).
    pub terminal_layers: Vec<(String, String)>,
    /// Terminal name -> resolved electrical net. Empty until connectivity is
    /// traced; the layer mapping in `terminal_layers` is always available.
    pub terminals: Vec<(String, String)>,
    /// Unit vector of the channel/current-flow axis, if determined.
    pub channel_axis: Option<(f64, f64)>,
    /// Channel orientation in degrees, if determined.
    pub channel_angle_deg: Option<f64>,
    /// Axis-aligned bounds of the recognized device region (meters).
    pub bbox: Option<Bbox>,
    /// Reference location (lower-left of the seed region, meters).
    pub location: Option<(f64, f64)>,
    pub params: BTreeMap<String, f64>,
    /// Recognition confidence in `[0, 1]`.
    pub confidence: f64,
}

/// One instance line parsed from a golden LVS SPICE netlist.
#[derive(Debug, Clone, PartialEq)]
pub struct SpiInstance {
    pub name: String,
    /// Connected node names, in order.
    pub nodes: Vec<String>,
    pub model: String,
    pub params: BTreeMap<String, f64>,
    /// `$X` location in database units, if present.
    pub x: Option<f64>,
    /// `$Y` location in database units, if present.
    pub y: Option<f64>,
}

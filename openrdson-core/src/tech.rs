//! Process technology / stack IR, produced by `techfile`.

use std::collections::BTreeMap;

/// Mapping between logical (rulefile) layer names and physical (techfile) layers,
/// parsed from `sft_cci.map`.
#[derive(Debug, Clone, Default)]
pub struct LayerMap {
    /// logical layer -> physical layer for conductors.
    pub conducting: BTreeMap<String, String>,
    /// logical layer -> physical layer for vias/contacts.
    pub via: BTreeMap<String, String>,
}

/// A conducting layer (metal, poly, diffusion-contact metal, pad, ...).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Conductor {
    pub name: String,
    /// Absolute bottom height of the layer in meters (ICT is µm; converted on load).
    pub height: Option<f64>,
    pub thickness: Option<f64>,
    /// Resistivity table. ICT encodes width-dependent tables as
    /// `rho width rho width ...`; a single value means bulk resistivity.
    pub resistivity: Vec<f64>,
    pub layer_type: Option<String>,
    pub gate_forming_layer: bool,
    pub temp_tc1: Option<f64>,
    pub temp_tc2: Option<f64>,
}

impl Conductor {
    pub fn z_bottom(&self) -> Option<f64> {
        self.height
    }
    pub fn z_top(&self) -> Option<f64> {
        match (self.height, self.thickness) {
            (Some(h), Some(t)) => Some(h + t),
            _ => None,
        }
    }
    /// Width-dependent sheet resistance (Ω/sq). The ICT encodes `rho width`
    /// pairs (µm); `width_um` selects/interpolates. A single value is used as-is.
    pub fn sheet_resistance(&self, width_um: f64) -> Option<f64> {
        interp_rho_width(&self.resistivity, width_um)
    }
    /// Temperature multiplier `1 + tc1·ΔT + tc2·ΔT²`.
    pub fn temp_factor(&self, t: f64, t_ref: f64) -> f64 {
        temp_factor(self.temp_tc1, self.temp_tc2, t, t_ref)
    }
}

impl Diffusion {
    /// Width-dependent sheet resistance (Ω/sq), as for [`Conductor`].
    pub fn sheet_resistance(&self, width_um: f64) -> Option<f64> {
        interp_rho_width(&self.resistivity, width_um)
    }
}

impl Via {
    /// Temperature multiplier for the via resistance.
    pub fn temp_factor(&self, t: f64, t_ref: f64) -> f64 {
        temp_factor(self.temp_tc1, self.temp_tc2, t, t_ref)
    }
}

fn temp_factor(tc1: Option<f64>, tc2: Option<f64>, t: f64, t_ref: f64) -> f64 {
    let dt = t - t_ref;
    1.0 + tc1.unwrap_or(0.0) * dt + tc2.unwrap_or(0.0) * dt * dt
}

/// Interpolate an ICT `rho width` table at `width_um`.
fn interp_rho_width(vals: &[f64], width_um: f64) -> Option<f64> {
    if vals.is_empty() {
        return None;
    }
    if vals.len() < 2 {
        return Some(vals[0]);
    }
    let mut pts: Vec<(f64, f64)> = vals
        .chunks(2)
        .filter(|c| c.len() == 2)
        .map(|c| (c[1], c[0])) // (width, rho)
        .collect();
    if pts.is_empty() {
        return Some(vals[0]);
    }
    pts.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    if width_um <= pts[0].0 {
        return Some(pts[0].1);
    }
    let last = pts[pts.len() - 1];
    if width_um >= last.0 {
        return Some(last.1);
    }
    for w in pts.windows(2) {
        let ((x0, v0), (x1, v1)) = (w[0], w[1]);
        if width_um >= x0 && width_um <= x1 {
            let t = if (x1 - x0).abs() < 1e-15 {
                0.0
            } else {
                (width_um - x0) / (x1 - x0)
            };
            return Some(v0 + t * (v1 - v0));
        }
    }
    Some(last.1)
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct Dielectric {
    pub name: String,
    pub height: Option<f64>,
    pub thickness: Option<f64>,
    pub dielectric_constant: Option<f64>,
    pub conformal: bool,
    pub expanded_from: Option<String>,
    pub top_thickness: Option<f64>,
    pub side_expand: Option<f64>,
    pub theta: Option<f64>,
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct Diffusion {
    pub name: String,
    pub height: Option<f64>,
    pub thickness: Option<f64>,
    pub resistivity: Vec<f64>,
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct Via {
    pub name: String,
    pub top_layer: Option<String>,
    pub bottom_layer: Option<String>,
    /// `area_resistance` table: `resistance area` pairs (ICT convention).
    pub area_resistance: Vec<f64>,
    pub min_width: Option<f64>,
    pub min_spacing: Option<f64>,
    pub min_top_encl: Option<f64>,
    pub min_bot_encl: Option<f64>,
    pub temp_tc1: Option<f64>,
    pub temp_tc2: Option<f64>,
}

/// A parsed ICT process stack.
#[derive(Debug, Clone, Default)]
pub struct TechStack {
    pub process: String,
    /// Reference temperature (°C) for the `temp_tc1`/`temp_tc2` coefficients.
    pub temp_reference: Option<f64>,
    pub conductors: Vec<Conductor>,
    pub dielectrics: Vec<Dielectric>,
    pub diffusions: Vec<Diffusion>,
    pub vias: Vec<Via>,
}

impl TechStack {
    /// Reference temperature for the temperature coefficients (default 25 °C).
    pub fn temp_ref(&self) -> f64 {
        self.temp_reference.unwrap_or(25.0)
    }
}

impl TechStack {
    pub fn conductor(&self, name: &str) -> Option<&Conductor> {
        self.conductors.iter().find(|c| c.name == name)
    }
    pub fn dielectric(&self, name: &str) -> Option<&Dielectric> {
        self.dielectrics.iter().find(|d| d.name == name)
    }
    pub fn diffusion(&self, name: &str) -> Option<&Diffusion> {
        self.diffusions.iter().find(|d| d.name == name)
    }
    pub fn via(&self, name: &str) -> Option<&Via> {
        self.vias.iter().find(|v| v.name == name)
    }
    pub fn conductor_names(&self) -> Vec<&str> {
        self.conductors.iter().map(|c| c.name.as_str()).collect()
    }
    pub fn via_names(&self) -> Vec<&str> {
        self.vias.iter().map(|v| v.name.as_str()).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn width_dependent_sheet_resistance() {
        let c = Conductor {
            resistivity: vec![0.08, 0.16, 0.06, 0.4, 0.1, 1.0],
            ..Default::default()
        };
        assert!((c.sheet_resistance(0.16).unwrap() - 0.08).abs() < 1e-12);
        assert!((c.sheet_resistance(0.4).unwrap() - 0.06).abs() < 1e-12);
        // interpolated between (0.16, 0.08) and (0.4, 0.06)
        assert!((c.sheet_resistance(0.28).unwrap() - 0.07).abs() < 1e-12);
        // clamped outside the table
        assert!((c.sheet_resistance(0.01).unwrap() - 0.08).abs() < 1e-12);
        assert!((c.sheet_resistance(10.0).unwrap() - 0.1).abs() < 1e-12);
    }

    #[test]
    fn single_value_resistivity_is_constant() {
        let c = Conductor {
            resistivity: vec![0.09],
            ..Default::default()
        };
        assert_eq!(c.sheet_resistance(0.5), Some(0.09));
    }

    #[test]
    fn temp_factor_applies_coefficients() {
        let c = Conductor {
            temp_tc1: Some(0.003),
            temp_tc2: Some(0.0),
            ..Default::default()
        };
        assert!((c.temp_factor(125.0, 25.0) - 1.3).abs() < 1e-9);
        assert!((c.temp_factor(25.0, 25.0) - 1.0).abs() < 1e-9);
    }
}

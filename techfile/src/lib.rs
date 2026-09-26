//! Process technology handling for OpenRDSon.
//!
//! Parses the ICT stack (via `openrdson-io`), parses the logical->physical layer
//! map (`sft_cci.map`) and provides z-stack validation.

use openrdson_core::tech::{Conductor, LayerMap, TechStack};
use openrdson_io::read_ict_file;
use std::io;
use std::path::Path;

pub use openrdson_io::{parse_ict, parse_sft_cci_map, read_sft_cci_map};

pub fn load_tech_stack<P: AsRef<Path>>(ict_path: P) -> io::Result<TechStack> {
    read_ict_file(ict_path)
}

/// Return conductors sorted by bottom z-height (for cross-section assembly).
pub fn conductors_by_height(stack: &TechStack) -> Vec<&Conductor> {
    let mut c: Vec<&Conductor> = stack.conductors.iter().collect();
    c.sort_by(|a, b| {
        a.height
            .partial_cmp(&b.height)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    c
}

/// Basic consistency diagnostics for the z-stack. Returns human-readable issues.
pub fn validate_z_stack(stack: &TechStack) -> Vec<String> {
    let mut issues = Vec::new();
    if stack.conductors.is_empty() {
        issues.push("no conductors parsed".to_string());
    }
    for c in &stack.conductors {
        match (c.height, c.thickness) {
            (None, _) => issues.push(format!("conductor {} has no height", c.name)),
            (_, None) => issues.push(format!("conductor {} has no thickness", c.name)),
            (Some(_), Some(t)) if t <= 0.0 => {
                issues.push(format!("conductor {} has non-positive thickness {t}", c.name))
            }
            _ => {}
        }
    }
    for v in &stack.vias {
        if v.top_layer.is_none() || v.bottom_layer.is_none() {
            issues.push(format!("via {} missing top/bottom layer", v.name));
        }
        if v.area_resistance.is_empty() {
            issues.push(format!("via {} has no area_resistance", v.name));
        }
    }
    issues
}

/// Resolve a logical layer to a physical techfile layer, if mapped.
pub fn physical_layer<'a>(map: &'a LayerMap, logical: &str) -> Option<&'a str> {
    map.conducting
        .get(logical)
        .or_else(|| map.via.get(logical))
        .map(|s| s.as_str())
}

//! Generic, template-driven device recognition from layout + seed layers.
//!
//! v0 scope (M0–M2 of `goal.md`):
//! - find each template's **seed layer** by name,
//! - one device instance per seed polygon,
//! - recover the instance name from the seed polygon's LVS property
//!   (`PROPATTR = 7` carries the instance name, e.g. `X0`),
//! - attach the template's terminal→layer map,
//! - derive the **channel/current-flow axis** from the seed region geometry.
//!
//! Net resolution (tracing terminal layers to electrical nets) and parameter
//! extraction from `*_prop` layers are deferred to later milestones.

use openrdson_core::device::{DeviceInstance, DeviceTemplate};
use openrdson_core::geometry::Bbox;
use openrdson_core::layout::{LayoutIR, LayoutPolygon};
use std::collections::BTreeMap;

/// `PROPATTR` value that carries the instance name.
pub const PROP_INSTANCE_NAME: i16 = 7;

#[derive(Debug, Clone)]
pub struct RecognitionConfig {
    /// Skip seed polygons smaller than this area (m²) to reject slivers.
    pub min_seed_area_m2: f64,
}

impl Default for RecognitionConfig {
    fn default() -> Self {
        // ~0.01 um² in m²
        Self {
            min_seed_area_m2: 1e-14,
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct RecognitionResult {
    pub instances: Vec<DeviceInstance>,
    pub diagnostics: Vec<String>,
    /// Number of seed polygons found per template.
    pub seeds_per_template: BTreeMap<String, usize>,
}

fn layer_number(layout: &LayoutIR, logical: &str) -> Option<i16> {
    layout
        .layer_names
        .iter()
        .find(|(_, n)| n.as_str() == logical)
        .map(|(k, _)| *k)
}

fn instance_name(poly: &LayoutPolygon) -> Option<String> {
    poly.props
        .iter()
        .find(|(a, _)| *a == PROP_INSTANCE_NAME)
        .map(|(_, v)| v.clone())
}

/// Channel/current-flow axis from a seed region's bounds.
///
/// For an LDMOS the seed is a tall thin rectangle whose long axis is
/// the gate/finger direction; current flows along the short axis. We therefore
/// return the short-axis direction. This is deliberately simple for v0 and is
/// cross-checked by `device-recognition-verify`.
fn channel_axis_from_bbox(bbox: &Bbox) -> Option<((f64, f64), f64)> {
    let w = bbox.width();
    let h = bbox.height();
    if w <= 0.0 || h <= 0.0 || (w - h).abs() < f64::EPSILON {
        return None;
    }
    if h >= w {
        // Finger direction is Y; current flows along X.
        Some(((1.0, 0.0), 0.0))
    } else {
        // Finger direction is X; current flows along Y.
        Some(((0.0, 1.0), 90.0))
    }
}

/// Recognize device instances from the layout using the parsed `devtab`
/// templates.
pub fn recognize(
    layout: &LayoutIR,
    templates: &[DeviceTemplate],
    cfg: &RecognitionConfig,
) -> RecognitionResult {
    let mut result = RecognitionResult::default();

    for t in templates {
        let Some(seed_layer) = layer_number(layout, &t.name) else {
            continue; // This template's seed layer is not in this layout.
        };
        let seeds: Vec<&LayoutPolygon> = layout
            .polygons
            .iter()
            .filter(|p| p.layer.layer == seed_layer)
            .collect();
        if seeds.is_empty() {
            continue;
        }
        result.seeds_per_template.insert(t.name.clone(), seeds.len());

        for (i, seed) in seeds.iter().enumerate() {
            let bbox = Bbox::from_points(&seed.points);
            if bbox.width() * bbox.height() < cfg.min_seed_area_m2 {
                result.diagnostics.push(format!(
                    "{}: seed {} too small ({:.3e} m²), skipped",
                    t.name,
                    i,
                    bbox.width() * bbox.height()
                ));
                continue;
            }
            let name = instance_name(seed).unwrap_or_else(|| format!("{}_{}", t.name, i));
            let (axis, angle) = match channel_axis_from_bbox(&bbox) {
                Some((a, ang)) => (Some(a), Some(ang)),
                None => (None, None),
            };
            // Confidence: seed present + all terminal layers resolvable.
            let all_layers = t
                .terminals
                .iter()
                .all(|(_, l)| layer_number(layout, l).is_some());
            let confidence = if all_layers { 1.0 } else { 0.5 };
            if !all_layers {
                result.diagnostics.push(format!(
                    "{}: instance {} has unresolved terminal layers",
                    t.name, name
                ));
            }

            result.instances.push(DeviceInstance {
                name,
                kind: t.kind,
                model: t.model.clone(),
                template: t.name.clone(),
                terminal_layers: t.terminals.clone(),
                terminals: Vec::new(),
                channel_axis: axis,
                channel_angle_deg: angle,
                bbox: Some(bbox),
                location: Some((bbox.min_x, bbox.min_y)),
                params: BTreeMap::new(),
                confidence,
            });
        }
    }

    result
}

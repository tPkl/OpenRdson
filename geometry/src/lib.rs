//! 3D geometry construction for OpenRDSon.
//!
//! v0 implements the 2.5D path: extrude axis-aligned rectangular footprints
//! into boxes, and general footprints into prisms. Polygon boolean cleanup and
//! conformal dielectrics are later milestones.

use openrdson_core::geometry::{
    is_axis_aligned_rectangle, Bbox, Point, Solid, SolidModel, SolidRegion,
};
use openrdson_core::layout::LayoutIR;
use openrdson_core::tech::{LayerMap, TechStack, Via};
use std::collections::BTreeMap;

/// Extrude a polygon footprint between `z_bottom` and `z_top`.
///
/// Axis-aligned rectangles become [`Solid::Box`] (fast structured meshing);
/// everything else becomes a [`Solid::Prism`].
pub fn extrude_polygon(
    points: &[Point],
    z_bottom: f64,
    z_top: f64,
    material_id: u32,
    net: Option<String>,
    layer: Option<String>,
    device_ref: Option<String>,
) -> SolidRegion {
    let solid = if is_axis_aligned_rectangle(points) {
        let b = Bbox::from_points(points);
        Solid::Box {
            min_x: b.min_x,
            min_y: b.min_y,
            max_x: b.max_x,
            max_y: b.max_y,
            z_bottom,
            z_top,
        }
    } else {
        Solid::Prism {
            footprint: points.to_vec(),
            z_bottom,
            z_top,
        }
    };
    SolidRegion {
        material_id,
        net,
        layer,
        device_ref,
        source_polygon: None,
        solid,
    }
}

/// Total volume of a set of regions (m³).
pub fn total_volume(regions: &[SolidRegion]) -> f64 {
    regions.iter().map(|r| r.volume()).sum()
}

/// Clip a solid model to a rectangular window, splitting polygons that cross the
/// window boundary. Used to extract a fine local mesh around a device without
/// pulling in full-die routing polygons.
pub fn clip_model(model: &SolidModel, window: Bbox) -> SolidModel {
    use openrdson_core::geomops::polygon_intersection_rects;
    let ring = vec![
        Point::new(window.min_x, window.min_y),
        Point::new(window.max_x, window.min_y),
        Point::new(window.max_x, window.max_y),
        Point::new(window.min_x, window.max_y),
        Point::new(window.min_x, window.min_y),
    ];
    let mut regions = Vec::new();
    for r in &model.regions {
        let push_box = |x0: f64, y0: f64, x1: f64, y1: f64, z0: f64, z1: f64, out: &mut Vec<SolidRegion>| {
            if x1 > x0 && y1 > y0 && z1 > z0 {
                out.push(SolidRegion {
                    material_id: r.material_id,
                    net: r.net.clone(),
                    layer: r.layer.clone(),
                    device_ref: r.device_ref.clone(),
                    source_polygon: r.source_polygon,
                    solid: Solid::Box {
                        min_x: x0,
                        min_y: y0,
                        max_x: x1,
                        max_y: y1,
                        z_bottom: z0,
                        z_top: z1,
                    },
                });
            }
        };
        match &r.solid {
            Solid::Box {
                min_x,
                min_y,
                max_x,
                max_y,
                z_bottom,
                z_top,
            } => push_box(
                min_x.max(window.min_x),
                min_y.max(window.min_y),
                max_x.min(window.max_x),
                max_y.min(window.max_y),
                *z_bottom,
                *z_top,
                &mut regions,
            ),
            Solid::Prism {
                footprint,
                z_bottom,
                z_top,
            } => {
                for rect in polygon_intersection_rects(footprint, &ring) {
                    push_box(rect.x0, rect.y0, rect.x1, rect.y1, *z_bottom, *z_top, &mut regions);
                }
            }
        }
    }
    SolidModel {
        regions,
        material_names: model.material_names.clone(),
    }
}

/// Diagnostics from [`assemble_stack`].
#[derive(Debug, Clone, Default)]
pub struct StackAssemblyReport {
    pub regions: usize,
    pub boxes: usize,
    pub prisms: usize,
    /// Logical layers present in the layout but not in the layer map.
    pub unmatched_layers: BTreeMap<String, usize>,
    /// Physical layers that had no resolvable z-range.
    pub no_z_range: BTreeMap<String, usize>,
    /// Logical layers intentionally skipped as non-geometry markers.
    pub skipped_markers: BTreeMap<String, usize>,
}

fn conductor_z(stack: &TechStack, name: &str) -> Option<(f64, f64)> {
    if let Some(c) = stack.conductor(name) {
        return Some((c.height?, c.height? + c.thickness?));
    }
    if let Some(d) = stack.diffusion(name) {
        return Some((d.height?, d.height? + d.thickness?));
    }
    None
}

fn via_z(stack: &TechStack, via: &Via) -> Option<(f64, f64)> {
    let bottom = via.bottom_layer.as_deref()?;
    let top = via.top_layer.as_deref()?;
    let (_, z_bot_top) = conductor_z(stack, bottom)?;
    let (z_top_bot, _) = conductor_z(stack, top)?;
    Some((z_bot_top, z_top_bot))
}

fn is_marker(logical: &str) -> bool {
    logical.starts_with("seed_")
        || logical.contains("_prop")
        || logical.contains("marker")
        || logical.starts_with("BLOCK_")
        || logical.contains("drawing")
}

/// Assemble tagged 3D solids for every conducting/via polygon in the layout,
/// using the ICT z-stack and the logical→physical layer map.
///
/// v0 policy: polygons whose logical layer starts with `net_` and is present in
/// the layer map are extruded; `seed_*`, `*_prop`, `*marker*`, `BLOCK_*` and
/// `*drawing*` layers are treated as non-geometry markers and skipped (they are
/// consumed by device recognition instead).
pub fn assemble_stack(
    layout: &LayoutIR,
    stack: &TechStack,
    map: &LayerMap,
) -> (SolidModel, StackAssemblyReport) {
    let mut model = SolidModel::default();
    let mut report = StackAssemblyReport::default();

    // Build the material table: conductors, then diffusions, then vias.
    for c in &stack.conductors {
        model.material_names.push(c.name.clone());
    }
    for d in &stack.diffusions {
        model.material_names.push(d.name.clone());
    }
    for v in &stack.vias {
        model.material_names.push(v.name.clone());
    }
    let material_id = |name: &str| -> Option<u32> {
        model
            .material_names
            .iter()
            .position(|n| n == name)
            .map(|i| (i + 1) as u32)
    };

    for (poly_index, poly) in layout.polygons.iter().enumerate() {
        let Some(logical) = layout.layer_name(poly.layer.layer) else {
            report
                .unmatched_layers
                .entry(format!("layer:{}", poly.layer.layer))
                .or_insert(0);
            continue;
        };
        if is_marker(logical) {
            *report.skipped_markers.entry(logical.to_string()).or_insert(0) += 1;
            continue;
        }
        let physical = map
            .conducting
            .get(logical)
            .or_else(|| map.via.get(logical));
        let Some(physical) = physical else {
            *report.unmatched_layers.entry(logical.to_string()).or_insert(0) += 1;
            continue;
        };

        let z = if let Some(v) = stack.via(physical) {
            via_z(stack, v)
        } else {
            conductor_z(stack, physical)
        };
        let Some((z0, z1)) = z else {
            *report.no_z_range.entry(physical.clone()).or_insert(0) += 1;
            continue;
        };
        if z1 <= z0 {
            *report.no_z_range.entry(physical.clone()).or_insert(0) += 1;
            continue;
        }

        let mut region = extrude_polygon(
            &poly.points,
            z0,
            z1,
            material_id(physical).unwrap_or(0),
            poly.net.clone(),
            Some(physical.clone()),
            None,
        );
        region.source_polygon = Some(poly_index);
        match region.solid {
            Solid::Box { .. } => report.boxes += 1,
            Solid::Prism { .. } => report.prisms += 1,
        }
        report.regions += 1;
        model.regions.push(region);
    }

    (model, report)
}

/// Convenience: extrude a rectangular region given explicit bounds.
#[allow(clippy::too_many_arguments)]
pub fn extrude_box(
    bbox: Bbox,
    z_bottom: f64,
    z_top: f64,
    material_id: u32,
    net: Option<String>,
    layer: Option<String>,
    device_ref: Option<String>,
) -> SolidRegion {
    SolidRegion {
        material_id,
        net,
        layer,
        device_ref,
        source_polygon: None,
        solid: Solid::Box {
            min_x: bbox.min_x,
            min_y: bbox.min_y,
            max_x: bbox.max_x,
            max_y: bbox.max_y,
            z_bottom,
            z_top,
        },
    }
}

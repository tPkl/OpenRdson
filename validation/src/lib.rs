//! Aggregated acceptance validation for OpenRDSon.
//!
//! Runs the generic (database-independent) analytic checks and returns a
//! machine-readable pass/fail list.

use openrdson_channel::ModelTable;
use openrdson_core::geometry::{Point, SolidModel};
use openrdson_device_recognition::{recognize, RecognitionConfig};
use openrdson_extraction::two_terminal_resistance;
use openrdson_geometry::extrude_polygon;
use openrdson_layout_db::{load_layout_ir, read_devtab_file, read_spi_instances};
use openrdson_meshing::{mesh_model, weld_nodes, MeshConfig};
use openrdson_techfile::{load_tech_stack, read_sft_cci_map};
use std::path::PathBuf;

#[derive(Debug, Clone)]
pub struct Check {
    pub name: String,
    pub passed: bool,
    pub detail: String,
}

/// Input file locations for a project (configurable via the YAML config).
#[derive(Debug, Clone, Default)]
pub struct ProjectPaths {
    pub cci_dir: PathBuf,
    pub layout: String,
    pub gds_map: String,
    pub ports: String,
    pub devtab: String,
    pub spi: String,
    /// CCI layout net-name table (`*.lnn`): node id -> net name.
    pub net_names: String,
    /// CCI netlist with `.DEVTMPLT` terminal-order definitions (`*.pin_xy_spi`).
    pub pin_xy: String,
    pub tech_ict: PathBuf,
    pub layer_map: PathBuf,
    pub model_csv: PathBuf,
}

impl ProjectPaths {
    pub fn layout_path(&self) -> PathBuf {
        self.cci_dir.join(&self.layout)
    }
    pub fn gds_map_path(&self) -> PathBuf {
        self.cci_dir.join(&self.gds_map)
    }
    pub fn ports_path(&self) -> PathBuf {
        self.cci_dir.join(&self.ports)
    }
    pub fn devtab_path(&self) -> PathBuf {
        self.cci_dir.join(&self.devtab)
    }
    pub fn spi_path(&self) -> PathBuf {
        self.cci_dir.join(&self.spi)
    }
    pub fn net_names_path(&self) -> PathBuf {
        self.cci_dir.join(&self.net_names)
    }
    pub fn pin_xy_path(&self) -> PathBuf {
        self.cci_dir.join(&self.pin_xy)
    }

    /// Every input file the extraction reads, in a stable order. Used to key
    /// the field cache so a changed database (or a regenerated `model.csv`)
    /// invalidates a cached solve.
    pub fn input_files(&self) -> Vec<PathBuf> {
        vec![
            self.layout_path(),
            self.gds_map_path(),
            self.ports_path(),
            self.devtab_path(),
            self.spi_path(),
            self.net_names_path(),
            self.pin_xy_path(),
            self.tech_ict.clone(),
            self.layer_map.clone(),
            self.model_csv.clone(),
        ]
    }
}

/// A user-defined voltage terminal (config input).
#[derive(Debug, Clone)]
pub struct BiasTerminal {
    pub name: String,
    /// Absolute voltage (V).
    pub voltage: f64,
    /// Contact centre (µm, relative to the layout origin).
    pub x_um: f64,
    pub y_um: f64,
    /// Contact footprint (µm); 0 binds only the nearest node.
    pub dx_um: f64,
    pub dy_um: f64,
    /// Optional contact level (physical ICT name, else logical gds.map name).
    pub layer: Option<String>,
}

/// A terminal resolved onto the sheet network: the nodes tied to its voltage.
#[derive(Debug, Clone)]
pub struct ResolvedTerminal {
    pub name: String,
    pub voltage: f64,
    pub layer: String,
    pub nodes: Vec<u32>,
}

/// What to do when an active device has no channel model available.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MissingModel {
    Error,
    Open,
}

/// User-tunable meshing/extraction accuracy options (R3D-style).
#[derive(Debug, Clone)]
pub struct MeshSettings {
    /// Default 2D sheet cell size for the full-array extraction (m).
    pub sheet_cell_m: f64,
    /// Max cells per rectangle dimension (perf cap).
    pub sheet_cap: usize,
    /// Radius within which same-net vias are grouped (m); 0 disables grouping.
    pub via_group_radius_m: f64,
    /// Per-physical-layer sheet cell overrides (m).
    pub per_layer_cell_m: std::collections::BTreeMap<String, f64>,
    /// Explicit database-unit override (nm). The AGF `UNITS` record can be
    /// malformed (both doubles equal -> a 1 m user unit), so when this is set it
    /// wins; otherwise the unit is inferred from the golden SPI finger width vs
    /// the recognized seed footprint.
    pub db_unit_nm: Option<f64>,
    /// Input file locations.
    pub paths: ProjectPaths,
    /// User voltage terminals (absolute bias).
    pub terminals: Vec<BiasTerminal>,
    /// Behavior when an active device has no model.
    pub missing_model: MissingModel,
    /// Operating temperature (°C) for the interconnect temperature coefficients.
    pub temperature_c: f64,
    /// Adaptive (quadtree) sheet-meshing controls.
    pub adaptive: openrdson_extraction::adaptive::AdaptiveConfig,
}

impl Default for MeshSettings {
    fn default() -> Self {
        Self {
            sheet_cell_m: 2e-6,
            sheet_cap: 30,
            via_group_radius_m: 0.5e-6,
            per_layer_cell_m: std::collections::BTreeMap::new(),
            db_unit_nm: None,
            paths: ProjectPaths::default(),
            terminals: Vec::new(),
            missing_model: MissingModel::Error,
            temperature_c: 27.0,
            adaptive: openrdson_extraction::adaptive::AdaptiveConfig::default(),
        }
    }
}

impl MeshSettings {
    /// Parse a mesh config: lines `LayerName lateral_um` (`#` starts a comment).
    pub fn load_config_str(&mut self, text: &str) {
        for line in text.lines() {
            let line = line.split('#').next().unwrap_or("").trim();
            if line.is_empty() {
                continue;
            }
            let t: Vec<&str> = line.split_whitespace().collect();
            if t.len() >= 2 {
                if let Ok(um) = t[1].parse::<f64>() {
                    if um > 0.0 {
                        self.per_layer_cell_m.insert(t[0].to_string(), um * 1e-6);
                    }
                }
            }
        }
    }

    pub fn load_config_file(&mut self, path: &std::path::Path) -> Result<(), String> {
        let text = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
        self.load_config_str(&text);
        Ok(())
    }
}

/// Infer the true database unit (m) from the golden SPI finger width vs the
/// recognized seed footprint. The AGF `UNITS` record can be malformed (both
/// doubles equal, implying a 1 m user unit); the seed's finger-direction extent
/// should equal the SPI `w`, so `dbu = w · dbu_loaded / extent_loaded`.
fn infer_db_unit(
    layout: &openrdson_core::layout::LayoutIR,
    devices: &[openrdson_core::device::DeviceInstance],
    spi: &[openrdson_core::device::SpiInstance],
    override_nm: Option<f64>,
) -> f64 {
    if let Some(nm) = override_nm {
        if nm > 0.0 {
            return nm * 1e-9;
        }
    }
    let dbu = layout.db_unit_meters;
    let mut used = vec![false; spi.len()];
    let mut implied: Vec<f64> = Vec::new();
    for inst in devices {
        let (Some(bb), Some(loc)) = (inst.bbox, inst.location) else {
            continue;
        };
        let mut best: Option<(usize, f64)> = None;
        for (i, g) in spi.iter().enumerate() {
            if used[i] {
                continue;
            }
            let (Some(gx), Some(gy)) = (g.x, g.y) else { continue };
            let d = ((gx * dbu - loc.0).powi(2) + (gy * dbu - loc.1).powi(2)).sqrt();
            if best.map(|(_, bd)| d < bd).unwrap_or(true) {
                best = Some((i, d));
            }
        }
        let Some((i, _)) = best else { continue };
        used[i] = true;
        let Some(w) = spi[i].params.get("w").copied() else {
            continue;
        };
        let extent = bb.width().max(bb.height());
        if w > 0.0 && extent > 0.0 {
            implied.push(w * dbu / extent);
        }
    }
    if implied.is_empty() {
        return dbu;
    }
    implied.sort_by(|a, b| a.partial_cmp(b).unwrap());
    implied[implied.len() / 2]
}

/// Resolve a terminal's contact layer to a physical layer name: an explicit
/// physical (ICT) name, an explicit logical (gds.map) name, or the port's
/// logical layer mapped through the layer map.
fn resolve_contact_layer(
    stack: &openrdson_core::tech::TechStack,
    map: &openrdson_core::tech::LayerMap,
    explicit: Option<&str>,
    port_layer: Option<&str>,
) -> Option<String> {
    if let Some(l) = explicit {
        if stack.conductor(l).is_some() || stack.diffusion(l).is_some() || stack.via(l).is_some() {
            return Some(l.to_string());
        }
        if let Some(ph) = map.conducting.get(l).or_else(|| map.via.get(l)) {
            return Some(ph.clone());
        }
        return Some(l.to_string());
    }
    port_layer.and_then(|l| map.conducting.get(l).or_else(|| map.via.get(l)).cloned())
}

/// Correct a layout's database unit in place using the golden SPI vs the
/// recognized seed footprint (see [`infer_db_unit`]).
pub fn correct_db_unit(layout: &mut openrdson_core::layout::LayoutIR, settings: &MeshSettings) {
    let Ok(templates) = read_devtab_file(settings.paths.devtab_path()) else {
        return;
    };
    let spi = read_spi_instances(settings.paths.spi_path()).unwrap_or_default();
    let prelim = recognize(layout, &templates, &RecognitionConfig::default());
    let dbu = infer_db_unit(layout, &prelim.instances, &spi, settings.db_unit_nm);
    if (dbu - layout.db_unit_meters).abs() > 1e-15 {
        openrdson_core::log_info!(
            "database unit corrected {:.4} nm -> {:.4} nm (golden SPI vs seed)",
            layout.db_unit_meters * 1e9,
            dbu * 1e9
        );
        let factor = dbu / layout.db_unit_meters;
        openrdson_layout_db::rescale_layout(layout, factor);
    }
}

/// Load the layout with its database unit corrected, so every consumer
/// (sheet network, stack solids, recognition, visualization) shares one scale.
pub fn load_layout(
    settings: &MeshSettings,
) -> Result<openrdson_core::layout::LayoutIR, String> {
    let paths = &settings.paths;
    let mut layout = load_layout_ir(paths.layout_path(), paths.gds_map_path())
        .map_err(|e| e.to_string())?;
    correct_db_unit(&mut layout, settings);
    Ok(layout)
}

/// A device channel with its drain/source/gate nodes, connectivity components
/// (stable across mesh rebuilds), finger width, and recognized seed footprint.
#[derive(Debug, Clone)]
pub struct ChannelDevice {
    pub d: u32,
    pub s: u32,
    /// Gate node (for Vgs), if the device has a gate terminal.
    pub gate: Option<u32>,
    /// Connectivity components of the drain/source/gate terminals. These are
    /// layout-derived and stable across mesh rebuilds, so the channel can be
    /// re-projected onto an adaptive mesh (node indices are not stable).
    pub d_comp: usize,
    pub s_comp: usize,
    pub gate_comp: Option<usize>,
    /// Per-instance channel width (m): the drawn finger width `w` from the SPI.
    pub weff: f64,
    /// Recognized seed footprint (m).
    pub bbox: Option<openrdson_core::geometry::Bbox>,
    pub nx: usize,
    pub ny: usize,
    pub spi_name: Option<String>,
}

/// A device channel as solved with its local (IR-drop aware) bias.
#[derive(Debug, Clone)]
pub struct ChannelState {
    pub d: u32,
    pub s: u32,
    /// Local gate-source voltage (0 if the device has no gate).
    pub vgs: f64,
    /// Local drain-source voltage (model sign convention).
    pub vds: f64,
    /// Channel resistance at the local bias (Ω).
    pub resistance: f64,
    pub conductance: f64,
    /// Current through the channel (A).
    pub current: f64,
    /// Joule power in the channel (W).
    pub power: f64,
    pub weff: f64,
    /// Seed footprint (m), for placing the channel in the visualization.
    pub bbox: Option<openrdson_core::geometry::Bbox>,
    pub nx: usize,
    pub ny: usize,
}

/// Per-cell scalar quantities, for visualization/export (KLayout colormaps,
/// VTK cell data). All values are per mesh cell (quad).
#[derive(Debug, Clone, Default)]
pub struct CellScalars {
    /// Average cell potential (V).
    pub potential: f64,
    /// Current-density magnitude `|∇u|/r_sheet` (A/m).
    pub current_density: f64,
    /// Joule power in the cell (W).
    pub power: f64,
    /// Current through the cell (A).
    pub current: f64,
    /// Resistance contribution of the cell (Ω), `power / I_total²`.
    pub resistance: f64,
    /// Quadtree refinement level of the cell.
    pub level: u8,
    /// Net (connectivity component) id.
    pub net: usize,
}

/// Per-cell fields from the sheet solve, for VTK export.
#[derive(Debug, Clone, Default)]
pub struct SheetField {
    /// `(x, y, potential)` per node.
    pub nodes: Vec<(f64, f64, f64)>,
    pub quads: Vec<[u32; 4]>,
    pub dx: Vec<f64>,
    pub dy: Vec<f64>,
    pub gx: Vec<f64>,
    pub gy: Vec<f64>,
    /// Sheet resistance (Ω/sq) per cell.
    pub r_sheet: Vec<f64>,
    /// Physical layer name per cell.
    pub layer: Vec<String>,
    /// Net (connectivity component) id per cell.
    pub component: Vec<usize>,
    /// Net id per node.
    pub node_component: Vec<usize>,
    /// Resolved bias terminals (name, voltage, contact nodes).
    pub terminals: Vec<ResolvedTerminal>,
    /// Total current (A) summed over the terminals at the solve bias.
    pub total_current: f64,
    /// Per-device channels, each with its own local-bias resistance.
    pub channels: Vec<ChannelState>,
}

/// Full-array 2.5D **sheet-resistance** extraction.
///
/// Each conducting layer is a 2D sheet resistor mesh; vias/contacts are lumped
/// resistors; each device is a bias-dependent channel resistor. This scales to
/// the whole die where 3D FEM cannot.
///
/// The conductor net graph is merged by connectivity component (see
/// `SheetNetwork::merge_nets_by_component`) and the solve uses a tight tolerance
/// so the open-state resistance is recovered. Validated by the `full-array`
/// acceptance check: `Rds(0) ≈ 1.3e10 Ω`, `Rds(-2.5) ≈ 1.2e1 Ω`.
pub struct SheetFullArray {
    net: openrdson_extraction::SheetNetwork,
    terminals: Vec<ResolvedTerminal>,
    channels: Vec<ChannelDevice>,
    /// Channel model; `None` for purely resistive layouts.
    model: Option<ModelTable>,
    missing_model: MissingModel,
    db_unit_meters: f64,
    r_sheet: std::collections::BTreeMap<String, f64>,
    /// Mesh inputs, retained so the network can be rebuilt on an adaptive grid.
    polygons: Vec<openrdson_extraction::SheetPolygon>,
    connections: Vec<(usize, usize)>,
    vias: Vec<openrdson_extraction::SheetVia>,
    term_specs: Vec<TermSpec>,
    sheet_cell_m: f64,
    sheet_cap: usize,
    per_layer_cell_m: std::collections::BTreeMap<String, f64>,
}

/// A terminal resolved to a fixed geometric contact (layer + rectangle), so it
/// can be re-projected onto a rebuilt (adaptive) mesh.
#[derive(Debug, Clone)]
pub struct TermSpec {
    pub name: String,
    pub voltage: f64,
    pub layer: String,
    pub x: f64,
    pub y: f64,
    pub hx: f64,
    pub hy: f64,
}

/// Outcome of an adaptive refine→solve run.
#[derive(Debug, Clone)]
pub struct AdaptiveReport {
    pub rds: f64,
    pub nodes: usize,
    pub edges: usize,
    pub cells: usize,
    pub iters: usize,
    /// Leaf count per refinement level (index 0 = base grid).
    pub levels: Vec<usize>,
    /// `(nodes, Rds, indicator)` per iteration.
    pub history: Vec<(usize, f64, f64)>,
}

/// Re-project terminal specs onto a (possibly rebuilt) network.
fn resolve_terminals_on(
    net: &openrdson_extraction::SheetNetwork,
    specs: &[TermSpec],
) -> Result<Vec<ResolvedTerminal>, String> {
    let mut out = Vec::new();
    for t in specs {
        let seed = net.nearest_on(&t.layer, t.x, t.y).ok_or_else(|| {
            format!("terminal '{}': no mesh node on layer {} near ({},{})", t.name, t.layer, t.x, t.y)
        })?;
        let comp = net.component[seed as usize];
        let mut nodes: Vec<u32> = if t.hx > 0.0 || t.hy > 0.0 {
            net.physical_nodes
                .get(&t.layer)
                .map(|ns| {
                    ns.iter()
                        .copied()
                        .filter(|&n| {
                            net.component[n as usize] == comp
                                && (net.x[n as usize] - t.x).abs() <= t.hx + 1e-15
                                && (net.y[n as usize] - t.y).abs() <= t.hy + 1e-15
                        })
                        .collect()
                })
                .unwrap_or_default()
        } else {
            vec![seed]
        };
        if nodes.is_empty() {
            // Footprint smaller than the mesh cell: bind the nearest node.
            openrdson_core::log_warn!(
                "terminal '{}': footprint covers no mesh node (smaller than the cell); binding the nearest node",
                t.name
            );
            nodes = vec![seed];
        }
        out.push(ResolvedTerminal {
            name: t.name.clone(),
            voltage: t.voltage,
            layer: t.layer.clone(),
            nodes,
        });
    }
    Ok(out)
}

/// A-priori geometry refinement, two independent parts:
///
/// * **Terminal region** (radius-gated, off when `radius <= 0`): terminals are
///   the current entry points, so each is refined as a distance-tapered region
///   — cells within `radius` of the contact go to full `max_level` at the
///   contact, tapering to level 0 at the radius.
///
/// * **Via point-weld fix** (always on): each via splits its containing leaf
///   once, so a via is not welded to a single coarse node. Cheap, bounded by
///   `max_cells`, and consistently improves accuracy (resolves the point
///   contact), so it is independent of the terminal-refinement radius.
fn a_priori_refine(
    grid: &mut openrdson_extraction::adaptive::AdaptiveGrid,
    vias: &[openrdson_extraction::SheetVia],
    terms: &[TermSpec],
    radius: f64,
    max_cells: usize,
) {
    // Terminals first: they are the voltage contacts and what the user expects
    // to see refined around, and there are only a handful.
    if radius > 0.0 {
        for t in terms {
            if grid.leaf_count() >= max_cells {
                return;
            }
            for c in grid.leaves() {
                if grid.leaf_count() >= max_cells {
                    return;
                }
                let cell = &grid.cells[c as usize];
                // Distance from the contact to the cell rectangle: 0 for the
                // cell containing the contact (so it gets full refinement),
                // growing toward `radius` at the edge of the region.
                let dx = (cell.x0 - t.x).max(0.0).max(t.x - cell.x1);
                let dy = (cell.y0 - t.y).max(0.0).max(t.y - cell.y1);
                let d = (dx * dx + dy * dy).sqrt();
                if d > radius {
                    continue;
                }
                // Full depth at the contact, tapering to zero at `radius`.
                let target = (grid.max_level as f64 * (1.0 - d / radius)).round() as u8;
                refine_to_level(grid, c, target, t.x, t.y, max_cells);
            }
        }
    }
    // Via point-weld fix: split the leaf containing each via once (nearest-leaf
    // fallback), so the point-welded via gets a local node without refining a
    // whole disc.
    for v in vias {
        if grid.leaf_count() >= max_cells {
            return;
        }
        let mut best: Option<(u32, f64)> = None;
        for c in grid.leaves() {
            let cell = &grid.cells[c as usize];
            if v.center.x >= cell.x0
                && v.center.x <= cell.x1
                && v.center.y >= cell.y0
                && v.center.y <= cell.y1
            {
                best = Some((c, 0.0));
                break;
            }
            let d = (cell.cx() - v.center.x).powi(2) + (cell.cy() - v.center.y).powi(2);
            if best.map(|(_, bd)| d < bd).unwrap_or(true) {
                best = Some((c, d));
            }
        }
        if let Some((c, _)) = best {
            // Split the via's cell once (to level 1), matching the original
            // point-weld fix: a cell with several vias is split once, not once
            // per via (which would over-refine to level 2, 3, ...).
            if grid.cells[c as usize].level < 1 {
                grid.refine(c);
            }
        }
    }
}

/// Refine `c` toward `target` level, always following the child that contains
/// `(sx, sy)` so the refinement stays centred on the seed rather than a fixed
/// quadrant.
fn refine_to_level(
    grid: &mut openrdson_extraction::adaptive::AdaptiveGrid,
    c: u32,
    target: u8,
    sx: f64,
    sy: f64,
    max_cells: usize,
) {
    let mut c = c;
    while grid.cells[c as usize].level < target && grid.leaf_count() < max_cells {
        if !grid.refine(c) {
            break;
        }
        let children = grid.cells[c as usize].children;
        c = children
            .iter()
            .copied()
            .find(|&ch| {
                let cc = &grid.cells[ch as usize];
                sx >= cc.x0 && sx <= cc.x1 && sy >= cc.y0 && sy <= cc.y1
            })
            .unwrap_or(children[0]);
    }
}

impl SheetFullArray {
    pub fn build() -> Result<Self, String> {
        Self::build_with(&MeshSettings::default())
    }

    pub fn build_with(settings: &MeshSettings) -> Result<Self, String> {
        let paths = &settings.paths;
        let mut layout = load_layout_ir(paths.layout_path(), paths.gds_map_path())
            .map_err(|e| e.to_string())?;
        let stack = load_tech_stack(&paths.tech_ict).map_err(|e| e.to_string())?;
        let map = read_sft_cci_map(&paths.layer_map).map_err(|e| e.to_string())?;
        // Correct a malformed `UNITS` record using the golden SPI finger width
        // cross-referenced against the recognized seed footprint.
        correct_db_unit(&mut layout, settings);
        layout.ports =
            openrdson_layout_db::read_ports_file(paths.ports_path(), layout.db_unit_meters)
                .map_err(|e| e.to_string())?;
        let conn = openrdson_layout_db::extract_connectivity(&layout, &stack, &map);

        // Width-dependent sheet resistance (Ω/sq) at the operating temperature:
        // the ICT `rho width` table is interpolated, then the `temp_tc`
        // coefficients applied relative to the process reference temperature.
        let t_op = settings.temperature_c;
        let t_ref = stack.temp_ref();
        let r_sheet = |physical: &str, width_um: f64| -> Option<f64> {
            if let Some(c) = stack.conductor(physical) {
                c.sheet_resistance(width_um)
                    .map(|r| r * c.temp_factor(t_op, t_ref))
            } else if let Some(d) = stack.diffusion(physical) {
                d.sheet_resistance(width_um)
            } else {
                None
            }
        };
        let is_marker = |logical: &str| {
            logical.starts_with("seed_")
                || logical.contains("_prop")
                || logical.contains("marker")
                || logical.starts_with("BLOCK_")
                || logical.contains("drawing")
        };

        let mut polygons = Vec::new();
        let mut r_sheet_by_layer: std::collections::BTreeMap<String, f64> =
            std::collections::BTreeMap::new();
        // Per metal/diffusion layer: (min Rs, max Rs, polygon count) for the
        // consistency log.
        let mut metal_layers: std::collections::BTreeMap<String, (f64, f64, usize)> =
            std::collections::BTreeMap::new();
        let mut remap = vec![usize::MAX; layout.polygons.len()];
        for (i, poly) in layout.polygons.iter().enumerate() {
            let Some(logical) = layout.layer_name(poly.layer.layer) else { continue };
            if is_marker(logical) {
                continue;
            }
            if let Some(physical) = map.conducting.get(logical) {
                // Line width proxy: the smaller bbox dimension (µm).
                let pbb = openrdson_core::geometry::Bbox::from_points(&poly.points);
                let width_um = (pbb.width().min(pbb.height()) * 1e6).max(1e-6);
                if let Some(rs) = r_sheet(physical, width_um) {
                    remap[i] = polygons.len();
                    r_sheet_by_layer
                        .entry(physical.clone())
                        .or_insert(rs);
                    metal_layers
                        .entry(physical.clone())
                        .and_modify(|e| {
                            e.0 = e.0.min(rs);
                            e.1 = e.1.max(rs);
                            e.2 += 1;
                        })
                        .or_insert((rs, rs, 1));
                    polygons.push(openrdson_extraction::SheetPolygon {
                        physical: physical.clone(),
                        points: poly.points.clone(),
                        component: conn.component_of[i],
                        r_sheet: rs,
                    });
                }
            }
        }
        for (layer, (rmin, rmax, n)) in &metal_layers {
            openrdson_core::log_info!(
                "metal {layer}: {n} polygons, solver Rs {rmin:.4}..{rmax:.4} ohm/sq (width-dependent)"
            );
        }
        let connections: Vec<(usize, usize)> = conn
            .same_layer_connections
            .iter()
            .filter_map(|&(a, b)| {
                let ra = remap.get(a).copied().unwrap_or(usize::MAX);
                let rb = remap.get(b).copied().unwrap_or(usize::MAX);
                if ra != usize::MAX && rb != usize::MAX {
                    Some((ra, rb))
                } else {
                    None
                }
            })
            .collect();

        let mut vias = Vec::new();
        // Per via layer: (R_ref, min solver R, max solver R, count) for the
        // consistency log.
        let mut via_layers: std::collections::BTreeMap<String, (f64, f64, f64, usize)> =
            std::collections::BTreeMap::new();
        for (i, poly) in layout.polygons.iter().enumerate() {
            let Some(logical) = layout.layer_name(poly.layer.layer) else { continue };
            let Some(physical) = map.via.get(logical) else { continue };
            let Some(via) = stack.via(physical) else { continue };
            let (Some(bottom), Some(top)) =
                (via.bottom_layer.as_deref(), via.top_layer.as_deref())
            else {
                continue;
            };
            let bb = openrdson_core::geometry::Bbox::from_points(&poly.points);
            // ICT `area_resistance` is a `(R_ref, A_ref)` pair: the resistance of
            // a via whose area is the reference `A_ref` (µm²). Scale by
            // `A_ref / A_actual` so a drawn via larger than the reference gets a
            // proportionally lower resistance.
            let (r_ref, a_ref) = match via.area_resistance.as_slice() {
                [r, a, ..] => (*r, *a),
                [r] => (*r, 0.0),
                _ => (1.0, 0.0),
            };
            let a_actual = (bb.width() * 1e6) * (bb.height() * 1e6);
            let r_via = if a_ref > 0.0 && a_actual > 0.0 {
                r_ref * a_ref / a_actual
            } else {
                r_ref
            };
            via_layers
                .entry(physical.clone())
                .and_modify(|e| {
                    e.1 = e.1.min(r_via);
                    e.2 = e.2.max(r_via);
                    e.3 += 1;
                })
                .or_insert((r_ref, r_via, r_via, 1));
            vias.push(openrdson_extraction::SheetVia {
                bottom: bottom.to_string(),
                top: top.to_string(),
                center: bb.center(),
                resistance: r_via * via.temp_factor(t_op, t_ref),
                component: conn.component_of[i],
                half_width: bb.width() * 0.5,
                half_height: bb.height() * 0.5,
            });
        }
        for (layer, (r_ref, rmin, rmax, n)) in &via_layers {
            openrdson_core::log_info!(
                "via {layer}: {n} polygons, ICT R_ref={r_ref:.4} ohm -> solver R {rmin:.4}..{rmax:.4} ohm (area-scaled)"
            );
        }

        // Sparsify dense via arrays by grouping nearby same-net vias.
        let raw_vias = vias.len();
        let vias = if settings.via_group_radius_m > 0.0 {
            openrdson_extraction::group_vias(&vias, settings.via_group_radius_m)
        } else {
            vias
        };
        openrdson_core::log_info!(
            "via grouping: {raw_vias} -> {} equivalent vias",
            vias.len()
        );

        let net = openrdson_extraction::SheetNetwork::build_with_cells(
            &polygons,
            &connections,
            &vias,
            &[],
            settings.sheet_cell_m,
            settings.sheet_cap,
            &settings.per_layer_cell_m,
        );
        // Per-layer mesh size (the solver's actual degrees of freedom), which
        // the polygon counts alone do not convey.
        let layer_stats = net.layer_stats();
        for (layer, (n, e, c)) in &layer_stats {
            openrdson_core::log_info!("mesh {layer}: {n} nodes / {e} edges / {c} quads");
        }
        let same_layer: usize = layer_stats.values().map(|(_, e, _)| *e).sum();
        openrdson_core::log_info!(
            "mesh cross-layer (via/merge) edges: {} of {} total",
            net.edge_count().saturating_sub(same_layer),
            net.edge_count()
        );

        // ---- Resolve the user's voltage terminals onto the mesh ----
        let mut terminals: Vec<ResolvedTerminal> = Vec::new();
        let mut term_specs: Vec<TermSpec> = Vec::new();
        for t in &settings.terminals {
            let port = layout.ports.iter().find(|p| p.name == t.name);
            let layer = resolve_contact_layer(
                &stack,
                &map,
                t.layer.as_deref(),
                port.map(|p| p.layer.as_str()),
            )
            .ok_or_else(|| {
                format!(
                    "terminal '{}': no contact layer (specify `layer:` or use a known port)",
                    t.name
                )
            })?;
            let (x, y) = if t.x_um != 0.0 || t.y_um != 0.0 {
                (t.x_um * 1e-6, t.y_um * 1e-6)
            } else if let Some(p) = port {
                (p.x, p.y)
            } else {
                return Err(format!("terminal '{}': no position and no matching port", t.name));
            };
            let seed = net.nearest_on(&layer, x, y).ok_or_else(|| {
                format!("terminal '{}': no mesh node on layer {layer} near ({x},{y})", t.name)
            })?;
            let comp = net.component[seed as usize];
            let hx = (t.dx_um.abs() * 1e-6) / 2.0;
            let hy = (t.dy_um.abs() * 1e-6) / 2.0;
            let mut nodes: Vec<u32> = if hx > 0.0 || hy > 0.0 {
                net.physical_nodes
                    .get(&layer)
                    .map(|ns| {
                        ns.iter()
                            .copied()
                            .filter(|&n| {
                                net.component[n as usize] == comp
                                    && (net.x[n as usize] - x).abs() <= hx + 1e-15
                                    && (net.y[n as usize] - y).abs() <= hy + 1e-15
                            })
                            .collect()
                    })
                    .unwrap_or_default()
            } else {
                vec![seed]
            };
            if nodes.is_empty() {
                // Footprint smaller than the mesh cell: bind the nearest node.
                openrdson_core::log_warn!(
                    "terminal '{}': {}x{} um footprint covers no mesh node (smaller than the cell); binding the nearest node",
                    t.name, t.dx_um, t.dy_um
                );
                nodes = vec![seed];
            }
            terminals.push(ResolvedTerminal {
                name: t.name.clone(),
                voltage: t.voltage,
                layer: layer.clone(),
                nodes,
            });
            term_specs.push(TermSpec {
                name: t.name.clone(),
                voltage: t.voltage,
                layer,
                x,
                y,
                hx,
                hy,
            });
        }
        // Every top port must have a voltage.
        for p in &layout.ports {
            if !terminals.iter().any(|r| r.name == p.name) {
                return Err(format!(
                    "port '{}' has no terminal definition (every top port must have a voltage)",
                    p.name
                ));
            }
        }
        openrdson_core::log_info!(
            "terminals: {} resolved ({} contact nodes)",
            terminals.len(),
            terminals.iter().map(|t| t.nodes.len()).sum::<usize>()
        );

        let templates = read_devtab_file(paths.devtab_path()).map_err(|e| e.to_string())?;
        let devices = recognize(&layout, &templates, &RecognitionConfig::default());
        // Golden LVS/CCI instances: authoritative finger width (`w`), array
        // shape (`nx`/`ny`) and placement.
        let spi = read_spi_instances(paths.spi_path()).unwrap_or_default();
        // CCI net-name table (node id -> name) and the `.DEVTMPLT` dictionary
        // (template -> ordered terminal names/layers).
        let net_names =
            openrdson_layout_db::read_net_names(paths.net_names_path()).unwrap_or_default();
        let dev_templates =
            openrdson_layout_db::read_cci_device_templates(paths.pin_xy_path()).unwrap_or_default();
        let dbu = layout.db_unit_meters;
        openrdson_core::log_debug!(
            "cci: {} net names (lnn), {} device templates (pin_xy_spi)",
            net_names.len(),
            dev_templates.len()
        );
        let mut used_spi = vec![false; spi.len()];
        let mut channels = Vec::new();
        // Resolve a device terminal role (drain/source/gate). A terminal whose
        // name matches the conventional role wins (D/G/S, PVIN/POUT, ...); else
        // fall back to the connectivity component of the devtab terminal layer
        // nearest the device (the generalized orientation mapping).
        let term_comp = |roles: &[&str]| -> Option<usize> {
            terminals
                .iter()
                .find(|t| roles.iter().any(|r| t.name.eq_ignore_ascii_case(r)))
                .map(|t| net.component[t.nodes[0] as usize])
        };
        let devtab_comp = |inst: &openrdson_core::device::DeviceInstance, layer: &str| -> Option<usize> {
            let logical = inst
                .terminal_layers
                .iter()
                .find(|(n, _)| n.eq_ignore_ascii_case(layer))
                .map(|(_, l)| l.as_str())?;
            let bb = inst.bbox?;
            let mut best: Option<(usize, f64)> = None;
            for (i, poly) in layout.polygons.iter().enumerate() {
                if layout.layer_name(poly.layer.layer) != Some(logical) {
                    continue;
                }
                let pb = openrdson_core::geometry::Bbox::from_points(&poly.points);
                if pb.max_x < bb.min_x
                    || pb.min_x > bb.max_x
                    || pb.max_y < bb.min_y
                    || pb.min_y > bb.max_y
                {
                    continue;
                }
                let c = conn.component_of[i];
                if c == usize::MAX {
                    continue;
                }
                let d =
                    (pb.center().x - bb.center().x).powi(2) + (pb.center().y - bb.center().y).powi(2);
                if best.map(|(_, bd)| d < bd).unwrap_or(true) {
                    best = Some((c, d));
                }
            }
            best.map(|(c, _)| c)
        };
        let role_node = |inst: &openrdson_core::device::DeviceInstance,
                         roles: &[&str],
                         layer: &str|
         -> Option<u32> {
            let comp = term_comp(roles).or_else(|| devtab_comp(inst, layer))?;
            let c = inst.bbox?.center();
            let mut best: Option<(u32, f64)> = None;
            for n in 0..net.node_count() {
                if net.component[n] != comp {
                    continue;
                }
                let d = (net.x[n] - c.x).powi(2) + (net.y[n] - c.y).powi(2);
                if best.map(|(_, bd)| d < bd).unwrap_or(true) {
                    best = Some((n as u32, d));
                }
            }
            best.map(|(n, _)| n)
        };
        // SPI LVS node id -> port name, and port name -> terminal component.
        let node_port: std::collections::BTreeMap<i64, String> = layout
            .ports
            .iter()
            .filter_map(|p| p.id.map(|id| (id, p.name.clone())))
            .collect();
        let term_comp_by_name: std::collections::BTreeMap<&str, usize> = terminals
            .iter()
            .map(|t| (t.name.as_str(), net.component[t.nodes[0] as usize]))
            .collect();
        let nearest_on_comp = |comp: usize, c: Point| -> Option<u32> {
            let mut best: Option<(u32, f64)> = None;
            for n in 0..net.node_count() {
                if net.component[n] != comp {
                    continue;
                }
                let d = (net.x[n] - c.x).powi(2) + (net.y[n] - c.y).powi(2);
                if best.map(|(_, bd)| d < bd).unwrap_or(true) {
                    best = Some((n as u32, d));
                }
            }
            best.map(|(n, _)| n)
        };
        // Device import summary: layout/schematic names, cell type, and W/L.
        let mut device_summary: Vec<(String, String, String, f64, f64)> = Vec::new();
        for inst in &devices.instances {
            let Some(bb) = inst.bbox else { continue };
            let c = bb.center();
            // Match the recognized seed to its nearest unused SPI instance.
            let mut best: Option<(usize, f64)> = None;
            for (i, g) in spi.iter().enumerate() {
                if used_spi[i] {
                    continue;
                }
                let (Some(gx), Some(gy)) = (g.x, g.y) else { continue };
                let dist = ((gx * dbu - c.x).powi(2) + (gy * dbu - c.y).powi(2)).sqrt();
                if best.map(|(_, bd)| dist < bd).unwrap_or(true) {
                    best = Some((i, dist));
                }
            }
            let spi_inst = best.map(|(i, _)| {
                used_spi[i] = true;
                &spi[i]
            });
            // Import-summary fields: finger width W (SPI), gate length L (seed
            // bbox short axis), and the schematic model name.
            let w = spi_inst.and_then(|g| g.params.get("w").copied()).unwrap_or(0.0);
            let l = bb.width().min(bb.height());
            let model = spi_inst
                .map(|g| g.model.clone())
                .or_else(|| inst.model.clone())
                .unwrap_or_default();
            device_summary.push((
                inst.template.clone(),
                model,
                inst.kind.as_str().to_string(),
                w,
                l,
            ));
            // Orientation from the golden SPI: the instance's terminal nets are
            // listed in the model's terminal order (d, g, s, sub), so the SPI
            // tells us which net is drain/gate/source.
            let spi_comp = |idx: usize| -> Option<usize> {
                let g = spi_inst?;
                let num = g.nodes.get(idx)?.parse::<i64>().ok()?;
                // Prefer the CCI net name (covers internal nets too).
                if let Some(name) = net_names.get(&num) {
                    if let Some(comp) = term_comp_by_name.get(name.as_str()) {
                        return Some(*comp);
                    }
                }
                let pname = node_port.get(&num)?;
                term_comp_by_name.get(pname.as_str()).copied()
            };
            // Terminal order from the `.DEVTMPLT` dictionary (via the SPI `$D`
            // template id), so roles are found by terminal NAME, not position.
            let tmpl = spi_inst
                .and_then(|g| g.params.get("$D").copied())
                .and_then(|d| dev_templates.get(&(d as i64)));
            let role_idx = |names: &[&str]| -> Option<usize> {
                tmpl?.terminals
                    .iter()
                    .position(|(n, _)| names.iter().any(|r| n.eq_ignore_ascii_case(r)))
            };
            let mut d = role_idx(&["d", "drain"])
                .and_then(spi_comp)
                .or_else(|| spi_comp(0))
                .and_then(|comp| nearest_on_comp(comp, c));
            let mut gate = role_idx(&["g", "gate"])
                .and_then(spi_comp)
                .or_else(|| spi_comp(1))
                .and_then(|comp| nearest_on_comp(comp, c));
            let mut s = role_idx(&["s", "source"])
                .and_then(spi_comp)
                .or_else(|| spi_comp(2))
                .and_then(|comp| nearest_on_comp(comp, c));
            // Fall back to name/devtab role resolution when the SPI is absent.
            if d.is_none() {
                d = role_node(inst, &["d", "drain"], "d");
            }
            if s.is_none() {
                s = role_node(inst, &["s", "source"], "s");
            }
            if gate.is_none() {
                gate = role_node(inst, &["g", "gate"], "g");
            }
            let (Some(d), Some(s)) = (d, s) else {
                continue;
            };
            // Per-instance channel width = the drawn finger width `w`. The CCI
            // `weff` is the *device's* total combined finger width (`w · nf`),
            // so using it once per SPI entry over-counts the total width by the
            // finger count when the array is listed one entry per finger.
            let weff = spi_inst
                .and_then(|g| g.params.get("w").copied())
                .filter(|w| *w > 0.0)
                .unwrap_or(1e-3);
            let dim = |k: &str| -> usize {
                spi_inst
                    .and_then(|g| g.params.get(k).copied())
                    .map(|v| v.max(1.0) as usize)
                    .unwrap_or(1)
            };
            channels.push(ChannelDevice {
                d,
                s,
                gate,
                d_comp: net.component[d as usize],
                s_comp: net.component[s as usize],
                gate_comp: gate.map(|g| net.component[g as usize]),
                weff,
                bbox: Some(bb),
                nx: dim("nx"),
                ny: dim("ny"),
                spi_name: spi_inst.map(|g| g.name.clone()),
            });
        }
        openrdson_core::log_info!(
            "channels: {} devices, {} matched to golden SPI",
            channels.len(),
            channels.iter().filter(|c| c.spi_name.is_some()).count()
        );
        // CCI device import summary: group by seed template and report the
        // layout/schematic names, cell type, and the finger width W / gate
        // length L. Logged once per process (the import is identical across the
        // Vgs sweep and the validation checks).
        static CCI_LOGGED: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);
        if !CCI_LOGGED.swap(true, std::sync::atomic::Ordering::Relaxed) {
            let mut by_tmpl: std::collections::BTreeMap<String, (usize, String, String, f64, f64)> =
                std::collections::BTreeMap::new();
            for (tmpl, model, kind, w, l) in &device_summary {
                let e = by_tmpl
                    .entry(tmpl.clone())
                    .or_insert((0, model.clone(), kind.clone(), *w, *l));
                e.0 += 1;
            }
            for (tmpl, (count, model, kind, w, l)) in &by_tmpl {
                openrdson_core::log_info!(
                    "CCI device: seed={tmpl} (layout)  model={model} (schematic)  type={kind}  instances={count}  W={:.3} um  L={:.3} um",
                    *w * 1e6,
                    *l * 1e6
                );
            }
            // Nets found (CCI .lnn) and which are used as simulation terminals.
            openrdson_core::log_info!("CCI nets: {} nets", net_names.len());
            for (id, name) in &net_names {
                let used = terminals
                    .iter()
                    .any(|t| t.name.eq_ignore_ascii_case(name));
                openrdson_core::log_info!(
                    "  net {id} -> {name}{}",
                    if used { "  [terminal]" } else { "" }
                );
            }
        }
        for t in &terminals {
            openrdson_core::log_debug!(
                "terminal {} = {:.3} V on {} ({} nodes, comps {:?})",
                t.name,
                t.voltage,
                t.layer,
                t.nodes.len(),
                t.nodes
                    .iter()
                    .map(|&n| net.component[n as usize])
                    .collect::<std::collections::BTreeSet<_>>()
            );
        }
        if let Some(ch) = channels.first() {
            openrdson_core::log_debug!(
                "device0: d={} (comp {}), s={} (comp {}), gate={:?} (comp {:?})",
                ch.d,
                net.component[ch.d as usize],
                ch.s,
                net.component[ch.s as usize],
                ch.gate,
                ch.gate.map(|g| net.component[g as usize])
            );
        }

        // Channel model is optional: required only when active devices exist.
        let model = std::fs::read_to_string(&paths.model_csv)
            .ok()
            .and_then(|csv| ModelTable::from_csv(&csv).ok());
        if model.is_none() && !channels.is_empty() {
            match settings.missing_model {
                MissingModel::Error => {
                    return Err(format!(
                        "{} active device(s) recognized but no channel model at {} \
                         (set solver.missing_model: open to continue)",
                        channels.len(),
                        paths.model_csv.display()
                    ));
                }
                MissingModel::Open => openrdson_core::log_warn!(
                    "{} active device(s) have no model; binding them to an open resistance",
                    channels.len()
                ),
            }
        }

        Ok(Self {
            net,
            terminals,
            channels,
            model,
            missing_model: settings.missing_model,
            db_unit_meters: layout.db_unit_meters,
            r_sheet: r_sheet_by_layer,
            polygons,
            connections,
            vias,
            term_specs,
            sheet_cell_m: settings.sheet_cell_m,
            sheet_cap: settings.sheet_cap,
            per_layer_cell_m: settings.per_layer_cell_m.clone(),
        })
    }

    /// The extracted network (nodes, edges, per-cell data).
    pub fn net(&self) -> &openrdson_extraction::SheetNetwork {
        &self.net
    }

    /// Clone this extraction with a different network; terminals are
    /// re-projected onto the new mesh from their stored geometric contacts.
    fn rebuild(&self, net: openrdson_extraction::SheetNetwork) -> Result<Self, String> {
        let terminals = resolve_terminals_on(&net, &self.term_specs)?;
        // Re-project each device's d/s/gate node indices onto the rebuilt mesh.
        // Node numbering changes with refinement; the connectivity components
        // (stored on the device) do not, so resolve from the device footprint.
        let channels: Vec<ChannelDevice> = self
            .channels
            .iter()
            .filter_map(|ch| {
                let c = ch.bbox?.center();
                let d = net.nearest_on_comp(ch.d_comp, c.x, c.y)?;
                let s = net.nearest_on_comp(ch.s_comp, c.x, c.y)?;
                let gate = ch
                    .gate_comp
                    .and_then(|gc| net.nearest_on_comp(gc, c.x, c.y));
                Some(ChannelDevice {
                    d,
                    s,
                    gate,
                    ..ch.clone()
                })
            })
            .collect();
        Ok(Self {
            net,
            terminals,
            channels,
            model: self.model.clone(),
            missing_model: self.missing_model.clone(),
            db_unit_meters: self.db_unit_meters,
            r_sheet: self.r_sheet.clone(),
            polygons: self.polygons.clone(),
            connections: self.connections.clone(),
            vias: self.vias.clone(),
            term_specs: self.term_specs.clone(),
            sheet_cell_m: self.sheet_cell_m,
            sheet_cap: self.sheet_cap,
            per_layer_cell_m: self.per_layer_cell_m.clone(),
        })
    }

    /// Adaptive (quadtree) refine→solve loop.
    ///
    /// Builds the base grid, optionally refines near vias/terminals a priori,
    /// then repeatedly: build the network, solve, estimate the per-cell power
    /// indicator, Dörfler-mark and refine, and re-balance. Stops on indicator
    /// convergence, the iteration cap, or the node budget.
    pub fn adaptive_solve(
        &self,
        cfg: &openrdson_extraction::adaptive::AdaptiveConfig,
        temperature: f64,
    ) -> Result<AdaptiveReport, String> {
        Ok(self.adaptive_solve_inner(cfg, temperature)?.0)
    }

    /// As [`Self::adaptive_solve`], but also returns the final adaptive
    /// [`SheetFullArray`] (network, terminals and channel model all on the
    /// refined mesh), so the caller can render the adaptive field.
    pub fn adaptive_solve_full(
        &self,
        cfg: &openrdson_extraction::adaptive::AdaptiveConfig,
        temperature: f64,
    ) -> Result<(AdaptiveReport, SheetFullArray), String> {
        let (report, fa) = self.adaptive_solve_inner(cfg, temperature)?;
        let fa = fa.ok_or_else(|| "adaptive solve produced no network".to_string())?;
        Ok((report, fa))
    }

    fn adaptive_solve_inner(
        &self,
        cfg: &openrdson_extraction::adaptive::AdaptiveConfig,
        temperature: f64,
    ) -> Result<(AdaptiveReport, Option<SheetFullArray>), String> {
        use openrdson_extraction::adaptive::{build_network, cell_power, AdaptiveGrid};
        let mut grid = AdaptiveGrid::from_polygons(
            &self.polygons,
            self.sheet_cell_m,
            self.sheet_cap,
            &self.per_layer_cell_m,
        );
        grid.max_level = cfg.max_level;
        // Bound the a-priori refinement so a dense via array cannot explode the
        // mesh: at most 4x the base leaf count (or the node budget / 4).
        let a_priori_cap = grid.leaf_count().saturating_mul(4).min(cfg.max_nodes / 4);
        a_priori_refine(
            &mut grid,
            &self.vias,
            &self.term_specs,
            cfg.refine_near_vias_m,
            a_priori_cap,
        );
        grid.balance();

        let mut last_fa: Option<SheetFullArray> = None;
        let mut history: Vec<(usize, f64, f64)> = Vec::new();
        let mut prev_eta = f64::INFINITY;
        let mut last = (f64::NAN, 0usize, 0usize, 0usize);
        let mut used_iters = 0usize;
        // Warm-start state carried across refinement iterations.
        let mut prev_u: Option<Vec<f64>> = None;
        let mut prev_r: Option<Vec<f64>> = None;
        let mut prev_net: Option<openrdson_extraction::SheetNetwork> = None;
        for it in 0..cfg.iters.max(1) {
            used_iters = it + 1;
            let net = build_network(&grid, &self.polygons, &self.connections, &self.vias, &[]);
            let (nodes, edges) = (net.node_count(), net.edge_count());
            if nodes > cfg.max_nodes {
                break;
            }
            {
                let comps: std::collections::BTreeSet<usize> = net
                    .component
                    .iter()
                    .copied()
                    .filter(|&c| c != usize::MAX)
                    .collect();
                openrdson_core::log_debug!(
                    "adaptive iter {}: {nodes} nodes / {edges} edges / {} leaves, {} net components",
                    it + 1,
                    grid.leaf_count(),
                    comps.len()
                );
            }
            let fa = self.rebuild(net)?;
            // Warm-start the solve: map the previous potential onto the rebuilt
            // mesh (node numbering changes with refinement) and carry the
            // previous channel resistances forward.
            let u0 = match (&prev_net, &prev_u) {
                (Some(old_net), Some(old_u)) => Some(openrdson_extraction::map_solution(
                    old_net,
                    old_u,
                    fa.net(),
                    grid.base_cell,
                )),
                _ => None,
            };
            let (u, channels) = fa.solve_local_bias_warm(
                temperature,
                1e-3,
                20,
                prev_r.as_deref(),
                u0.as_deref(),
            )?;
            let r = fa.rds(temperature).unwrap_or(f64::NAN);
            let power = cell_power(fa.net(), &u);
            let eta: f64 = power.iter().sum();
            prev_u = Some(u.clone());
            prev_r = Some(channels.iter().map(|c| c.resistance).collect());
            prev_net = Some(fa.net().clone());
            last_fa = Some(fa);
            history.push((nodes, r, eta));
            last = (r, nodes, edges, grid.leaf_count());
            if cfg.tol > 0.0 && (prev_eta - eta).abs() / eta.max(1e-30) < cfg.tol {
                break;
            }
            prev_eta = eta;
            // Dörfler marking: refine the highest-power cells until the marked
            // set carries `mark_frac` of the total indicator.
            let order = grid.leaf_order();
            let mut idx: Vec<usize> = (0..power.len()).collect();
            idx.sort_by(|&a, &b| {
                power[b]
                    .partial_cmp(&power[a])
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            let target = cfg.mark_frac * eta;
            let mut acc = 0.0;
            let mut marked: Vec<u32> = Vec::new();
            for &i in &idx {
                if acc >= target && !marked.is_empty() {
                    break;
                }
                acc += power[i];
                if i < order.len() {
                    marked.push(order[i]);
                }
            }
            let mut any = false;
            for c in marked {
                if grid.refine(c) {
                    any = true;
                }
            }
            if !any {
                break;
            }
            grid.balance();
        }
        let (rds, nodes, edges, cells) = last;
        let mut levels = vec![0usize; grid.max_level as usize + 1];
        for c in grid.cells.iter().filter(|c| c.is_leaf()) {
            let l = c.level as usize;
            if l < levels.len() {
                levels[l] += 1;
            }
        }
        Ok((
            AdaptiveReport {
                rds,
                nodes,
                edges,
                cells,
                iters: used_iters,
                levels,
                history,
            },
            last_fa,
        ))
    }

    /// Per-cell scalar quantities over the whole layout for visualization:
    /// `(physical_layer, corners, CellScalars)` plus the total terminal current.
    /// Computes the field once (same solve as [`Self::sheet_fields`]).
    pub fn cell_scalars(
        &self,
        temperature: f64,
    ) -> Result<(Vec<(String, [(f64, f64); 4], CellScalars)>, f64), String> {
        let sf = self.sheet_fields(temperature)?;
        let it2 = (sf.total_current * sf.total_current).max(1e-300);
        let mut out = Vec::with_capacity(sf.quads.len());
        for (qi, q) in sf.quads.iter().enumerate() {
            let phys = sf.layer[qi].clone();
            let corners = [
                (sf.nodes[q[0] as usize].0, sf.nodes[q[0] as usize].1),
                (sf.nodes[q[1] as usize].0, sf.nodes[q[1] as usize].1),
                (sf.nodes[q[2] as usize].0, sf.nodes[q[2] as usize].1),
                (sf.nodes[q[3] as usize].0, sf.nodes[q[3] as usize].1),
            ];
            let rs = sf.r_sheet[qi].max(1e-30);
            let g2 = sf.gx[qi] * sf.gx[qi] + sf.gy[qi] * sf.gy[qi];
            let j = g2.sqrt() / rs;
            let power = g2 / rs * sf.dx[qi] * sf.dy[qi];
            let pot = (sf.nodes[q[0] as usize].2
                + sf.nodes[q[1] as usize].2
                + sf.nodes[q[2] as usize].2
                + sf.nodes[q[3] as usize].2)
                / 4.0;
            let us: Vec<f64> = q.iter().map(|&n| sf.nodes[n as usize].2).collect();
            let hi = us.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            let lo = us.iter().cloned().fold(f64::INFINITY, f64::min);
            let dv = (hi - lo).abs();
            let current = if dv > 1e-15 { power / dv } else { 0.0 };
            out.push((
                phys,
                corners,
                CellScalars {
                    potential: pot,
                    current_density: j,
                    power,
                    current,
                    resistance: power / it2,
                    level: self.net.cell_level.get(qi).copied().unwrap_or(0),
                    net: sf.component[qi],
                },
            ));
        }
        Ok((out, sf.total_current))
    }

    /// Per-cell visualization data over the whole layout: physical layer, corner
    /// coordinates (m), average potential, and current-density magnitude
    /// (`|∇φ| / R_sheet`). Use these quads to redraw the layout as a mesh mosaic.
    #[allow(clippy::type_complexity)]
    pub fn viz_cells(
        &self,
        temperature: f64,
    ) -> Result<Vec<(String, [(f64, f64); 4], f64, f64)>, String> {
        let (u, _) = self.solve_local_bias(temperature, 1e-3, 20)?;
        let mut out = Vec::with_capacity(self.net.cell_count());
        for &cell in &self.net.cells {
            let (phys, c, dx, dy) = self.net.cell_info(cell);
            let (u0, u1, u2, u3) = (
                u[cell[0] as usize],
                u[cell[1] as usize],
                u[cell[2] as usize],
                u[cell[3] as usize],
            );
            let pot = (u0 + u1 + u2 + u3) / 4.0;
            let gx = ((u1 + u2) - (u0 + u3)) / (2.0 * dx);
            let gy = ((u2 + u3) - (u0 + u1)) / (2.0 * dy);
            let rs = self.r_sheet.get(&phys).copied().unwrap_or(1.0);
            let j = (gx * gx + gy * gy).sqrt() / rs.max(1e-30);
            out.push((phys, c, pot, j));
        }
        Ok(out)
    }

    pub fn db_unit_meters(&self) -> f64 {
        self.db_unit_meters
    }

    /// Node potentials and edge currents for visualization, in meters:
    /// `(nodes: (x, y, potential), edges: (x0, y0, x1, y1, current))`.
    #[allow(clippy::type_complexity)]
    pub fn viz_at(
        &self,
        temperature: f64,
    ) -> Result<(Vec<(f64, f64, f64)>, Vec<(f64, f64, f64, f64, f64)>), String> {
        let (u, channels) = self.solve_local_bias(temperature, 1e-3, 20)?;
        let extra: Vec<(u32, u32, f64)> = channels
            .iter()
            .filter(|c| c.conductance > 0.0)
            .map(|c| (c.d, c.s, c.conductance))
            .collect();
        let nodes = (0..self.net.node_count())
            .map(|i| (self.net.x[i], self.net.y[i], u[i]))
            .collect();
        let edges = self
            .net
            .edges
            .iter()
            .chain(extra.iter())
            .map(|&(i, j, g)| {
                let (i, j) = (i as usize, j as usize);
                (
                    self.net.x[i],
                    self.net.y[i],
                    self.net.x[j],
                    self.net.y[j],
                    g * (u[i] - u[j]),
                )
            })
            .collect();
        Ok((nodes, edges))
    }

    /// Solve with **per-device, bias-dependent** channel resistances, iterating
    /// until each device's local Vgs/Vds (set by the terminal bias and the metal
    /// IR drop) is self-consistent. The bias comes from the resolved terminals;
    /// each device's width comes from the golden SPI (`weff`).
    pub fn solve_local_bias(
        &self,
        temperature: f64,
        tol: f64,
        max_iter: usize,
    ) -> Result<(Vec<f64>, Vec<ChannelState>), String> {
        self.solve_local_bias_warm(temperature, tol, max_iter, None, None)
    }

    /// [`Self::solve_local_bias`] with warm-started initial guesses: `r0` seeds
    /// the channel resistances (length = `channels.len()`), `u0` seeds the
    /// potential (length = node count, already mapped to the current mesh).
    pub fn solve_local_bias_warm(
        &self,
        temperature: f64,
        tol: f64,
        max_iter: usize,
        r0: Option<&[f64]>,
        u0: Option<&[f64]>,
    ) -> Result<(Vec<f64>, Vec<ChannelState>), String> {
        self.solve_local_bias_terms(&self.terminals, temperature, tol, max_iter, r0, u0)
    }

    /// [`Self::solve_local_bias_warm`] with an explicit terminal set: the mesh
    /// and channels are unchanged, only the Dirichlet bias comes from
    /// `terminals`. Used by the Vgs sweep to re-solve at a shifted gate bias
    /// without rebuilding the network.
    fn solve_local_bias_terms(
        &self,
        terminals: &[ResolvedTerminal],
        temperature: f64,
        tol: f64,
        max_iter: usize,
        r0: Option<&[f64]>,
        u0: Option<&[f64]>,
    ) -> Result<(Vec<f64>, Vec<ChannelState>), String> {
        // Dirichlet nodes from every terminal.
        let fixed: Vec<(u32, f64)> = terminals
            .iter()
            .flat_map(|t| t.nodes.iter().map(move |&n| (n, t.voltage)))
            .collect();
        // Applied component potentials, for the initial guess.
        let mut comp_v: std::collections::BTreeMap<usize, f64> = std::collections::BTreeMap::new();
        for t in terminals {
            for &n in &t.nodes {
                comp_v.insert(self.net.component[n as usize], t.voltage);
            }
        }
        let is_open = self.missing_model == MissingModel::Open;
        let r_of = |ch: &ChannelDevice, vgs: f64, vds: f64| -> f64 {
            self.model
                .as_ref()
                .and_then(|m| m.r_device(temperature, vgs, vds, ch.weff))
                .filter(|r| r.is_finite() && *r > 0.0)
                .unwrap_or(if is_open { 1e18 } else { f64::INFINITY })
        };
        // Initial channel resistances: warm-started from `r0` when present,
        // else from the applied terminal voltages.
        let mut r: Vec<f64> = match r0 {
            Some(r0) if r0.len() == self.channels.len() => r0.to_vec(),
            _ => self
                .channels
                .iter()
                .map(|ch| {
                    let vs = comp_v.get(&self.net.component[ch.s as usize]).copied().unwrap_or(0.0);
                    let vg = ch
                        .gate
                        .map(|g| comp_v.get(&self.net.component[g as usize]).copied().unwrap_or(0.0))
                        .unwrap_or(0.0);
                    let vd = comp_v.get(&self.net.component[ch.d as usize]).copied().unwrap_or(0.0);
                    r_of(ch, vg - vs, vd - vs)
                })
                .collect(),
        };
        let mut u_prev: Option<Vec<f64>> = u0.map(|x| x.to_vec());
        let mut u: Vec<f64> = Vec::new();
        for _ in 0..max_iter.max(1) {
            let extra: Vec<(u32, u32, f64)> = self
                .channels
                .iter()
                .zip(r.iter())
                .filter(|(_, ri)| ri.is_finite() && **ri > 0.0)
                .map(|(ch, &ri)| (ch.d, ch.s, 1.0 / ri))
                .collect();
            u = self
                .net
                .solve_fixed_warm(&fixed, &extra, u_prev.as_deref(), 1e-8, 20000)?;
            u_prev = Some(u.clone());
            let mut max_change = 0.0f64;
            for (i, ch) in self.channels.iter().enumerate() {
                let vgs = ch.gate.map(|g| u[g as usize]).unwrap_or(0.0) - u[ch.s as usize];
                let vds = u[ch.d as usize] - u[ch.s as usize];
                let ri = r_of(ch, vgs, vds);
                if ri.is_finite() && ri > 0.0 {
                    let change = (ri - r[i]).abs() / r[i].max(1e-30);
                    max_change = max_change.max(change);
                    r[i] = ri;
                }
            }
            if max_change < tol {
                break;
            }
        }
        let channels: Vec<ChannelState> = self
            .channels
            .iter()
            .zip(r.iter())
            .map(|(ch, &ri)| {
                let vgs = ch.gate.map(|g| u[g as usize]).unwrap_or(0.0) - u[ch.s as usize];
                let vds = u[ch.d as usize] - u[ch.s as usize];
                let g = if ri.is_finite() && ri > 0.0 { 1.0 / ri } else { 0.0 };
                let i = g * vds.abs();
                ChannelState {
                    d: ch.d,
                    s: ch.s,
                    vgs,
                    vds,
                    resistance: ri,
                    conductance: g,
                    current: i,
                    power: vds.abs() * i,
                    weff: ch.weff,
                    bbox: ch.bbox,
                    nx: ch.nx,
                    ny: ch.ny,
                }
            })
            .collect();
        for t in terminals {
            let comp = self.net.component[t.nodes[0] as usize];
            let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
            let mut cn = 0usize;
            for (n, &c) in self.net.component.iter().enumerate() {
                if c == comp {
                    cn += 1;
                    lo = lo.min(u[n]);
                    hi = hi.max(u[n]);
                }
            }
            openrdson_core::log_debug!(
                "term {} = {}V on {} -> comp {comp}: {cn} nodes, u in [{lo:.3}, {hi:.3}]",
                t.name,
                t.voltage,
                t.layer
            );
        }
        {
            let mut v: Vec<(usize, f64)> =
                channels.iter().enumerate().map(|(i, c)| (i, c.resistance)).collect();
            v.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
            for (i, r) in v.iter().rev().take(4).chain(v.iter().take(3)) {
                let ch = &self.channels[*i];
                let (x, y) = ch
                    .bbox
                    .map(|b| (b.center().x, b.center().y))
                    .unwrap_or((0.0, 0.0));
                openrdson_core::log_debug!(
                    "chan[{i:3}] R={r:.3e} d={}(c{},u={:.3}) s={}(c{},u={:.3}) g={:?}(c{:?},u={:?}) @({:.1},{:.1})um",
                    ch.d,
                    self.net.component[ch.d as usize],
                    u[ch.d as usize],
                    ch.s,
                    self.net.component[ch.s as usize],
                    u[ch.s as usize],
                    ch.gate,
                    ch.gate.map(|g| self.net.component[g as usize]),
                    ch.gate.map(|g| u[g as usize]),
                    x * 1e6,
                    y * 1e6
                );
            }
        }
        Ok((u, channels))
    }

    pub fn sheet_fields(&self, temperature: f64) -> Result<SheetField, String> {
        let (u, channels) = self.solve_local_bias(temperature, 1e-3, 20)?;
        let extra: Vec<(u32, u32, f64)> = channels
            .iter()
            .filter(|c| c.conductance > 0.0)
            .map(|c| (c.d, c.s, c.conductance))
            .collect();
        let nodes: Vec<(f64, f64, f64)> = (0..self.net.node_count())
            .map(|i| (self.net.x[i], self.net.y[i], u[i]))
            .collect();
        let quads = self.net.cells.clone();
        let n = quads.len();
        let mut dx = Vec::with_capacity(n);
        let mut dy = Vec::with_capacity(n);
        let mut gx = Vec::with_capacity(n);
        let mut gy = Vec::with_capacity(n);
        let mut r_sheet = Vec::with_capacity(n);
        let mut layer = Vec::with_capacity(n);
        let mut component = Vec::with_capacity(n);
        for (qi, &cell) in quads.iter().enumerate() {
            let (phys, _c, cdx, cdy) = self.net.cell_info(cell);
            let (u0, u1, u2, u3) = (
                u[cell[0] as usize],
                u[cell[1] as usize],
                u[cell[2] as usize],
                u[cell[3] as usize],
            );
            dx.push(cdx);
            dy.push(cdy);
            gx.push(((u1 + u2) - (u0 + u3)) / (2.0 * cdx));
            gy.push(((u2 + u3) - (u0 + u1)) / (2.0 * cdy));
            // Use the polygon's own width-dependent sheet resistance (the value
            // the solver used), falling back to the per-layer summary.
            r_sheet.push(
                self.net
                    .cell_r_sheet
                    .get(qi)
                    .copied()
                    .unwrap_or_else(|| self.r_sheet.get(&phys).copied().unwrap_or(1.0)),
            );
            layer.push(phys);
            component.push(self.net.component[cell[0] as usize]);
        }
        // Total current: the current entering every terminal.
        let term_nodes: std::collections::BTreeSet<u32> = self
            .terminals
            .iter()
            .flat_map(|t| t.nodes.iter().copied())
            .collect();
        let mut total_current = 0.0;
        for &(a, b, g) in self.net.edges.iter().chain(extra.iter()) {
            if term_nodes.contains(&a) {
                total_current += g * (u[a as usize] - u[b as usize]);
            } else if term_nodes.contains(&b) {
                total_current += g * (u[b as usize] - u[a as usize]);
            }
        }
        Ok(SheetField {
            nodes,
            quads,
            dx,
            dy,
            gx,
            gy,
            r_sheet,
            layer,
            component,
            node_component: self.net.component.clone(),
            terminals: self.terminals.clone(),
            total_current: total_current.abs(),
            channels,
        })
    }

    /// Nearest sheet node on a physical layer to `(x, y)`. Used to map a via
    /// onto the potentials of the two layers it connects.
    pub fn nearest_node_on(&self, layer: &str, x: f64, y: f64) -> Option<u32> {
        self.net.nearest_on(layer, x, y)
    }

    /// Via connections for visualization: `(x, y, bottom_node, top_node,
    /// conductance)` for every network edge that joins two different physical
    /// layers. These are the actual solved via resistors, so their current is
    /// `conductance · ΔV` rather than a nearest-node approximation.
    #[allow(clippy::type_complexity)]
    pub fn via_links(&self) -> Vec<(f64, f64, u32, u32, f64)> {
        let mut out = Vec::new();
        for &(a, b, g) in &self.net.edges {
            if self.net.physical[a as usize] != self.net.physical[b as usize] {
                let x = (self.net.x[a as usize] + self.net.x[b as usize]) / 2.0;
                let y = (self.net.y[a as usize] + self.net.y[b as usize]) / 2.0;
                out.push((x, y, a, b, g));
            }
        }
        out
    }

    /// Diagnostic: how well the sheet conductor graph matches the connectivity
    /// net components (detects fragmentation of the D/S nets).
    pub fn diagnostic_report(&self) -> String {
        let n = self.net.node_count();
        let mut parent: Vec<u32> = (0..n as u32).collect();
        fn find(parent: &mut [u32], mut x: u32) -> u32 {
            while parent[x as usize] != x {
                parent[x as usize] = parent[parent[x as usize] as usize];
                x = parent[x as usize];
            }
            x
        }
        for &(i, j, _) in &self.net.edges {
            let (ri, rj) = (find(&mut parent, i), find(&mut parent, j));
            if ri != rj {
                parent[ri as usize] = rj;
            }
        }
        let mut sizes: std::collections::HashMap<u32, usize> = std::collections::HashMap::new();
        for i in 0..n {
            *sizes.entry(find(&mut parent, i as u32)).or_insert(0) += 1;
        }
        format!(
            "nodes={n} graph_components={} terminals={}",
            sizes.len(),
            self.terminals.len()
        )
    }

    /// A human-readable resistance breakdown of the full extraction: per-layer
    /// contribution, total Rdson, current draw and power, plus terminal,
    /// device and mesh summary. Each layer's `R` is its share of the total
    /// (`P_layer / I²`), so the column sums to the total Rdson.
    pub fn resistance_report(&self, temperature: f64) -> Result<String, String> {
        let (u, channels) = self.solve_local_bias(temperature, 1e-3, 20)?;

        // Current entering each terminal.
        let mut node_term: std::collections::BTreeMap<u32, usize> =
            std::collections::BTreeMap::new();
        for (ti, t) in self.terminals.iter().enumerate() {
            for &n in &t.nodes {
                node_term.entry(n).or_insert(ti);
            }
        }
        let mut term_i = vec![0.0f64; self.terminals.len()];
        for &(a, b, g) in &self.net.edges {
            let cur = g * (u[a as usize] - u[b as usize]);
            if let Some(&ti) = node_term.get(&a) {
                term_i[ti] += cur;
            } else if let Some(&ti) = node_term.get(&b) {
                term_i[ti] -= cur;
            }
        }
        // Identify the drain (power) and source terminals by name, so the gate
        // (a high-impedance control terminal) is never mistaken for a power
        // terminal. Fall back to the largest voltage difference when the names
        // aren't the conventional "D"/"S".
        let (pa, pb, ds_dv) = {
            let d = self.terminals.iter().position(|t| t.name.eq_ignore_ascii_case("d"));
            let s = self.terminals.iter().position(|t| t.name.eq_ignore_ascii_case("s"));
            match (d, s) {
                (Some(d), Some(s)) => (
                    d,
                    s,
                    (self.terminals[d].voltage - self.terminals[s].voltage).abs(),
                ),
                _ => {
                    let (mut pa, mut pb) = (0usize, 0usize);
                    let mut maxdv = -1.0f64;
                    for i in 0..self.terminals.len() {
                        for j in (i + 1)..self.terminals.len() {
                            let dv = (self.terminals[i].voltage - self.terminals[j].voltage).abs();
                            if dv > maxdv {
                                maxdv = dv;
                                pa = i;
                                pb = j;
                            }
                        }
                    }
                    (pa, pb, maxdv)
                }
            }
        };
        // Total current drawn by the "power" (drain) terminal: sum over every
        // contact on that same net (a net may have several contacts).
        let v_pa = self.terminals[pa].voltage;
        let i_draw: f64 = self
            .terminals
            .iter()
            .zip(&term_i)
            .filter(|(t, _)| (t.voltage - v_pa).abs() < 1e-12)
            .map(|(_, &i)| i.abs())
            .sum();
        let it2 = (i_draw * i_draw).max(1e-300);

        // Per-layer power from the actual network edges (same-layer edges are
        // the metal sheet resistors; cross-layer edges are the vias).
        let mut layer_power: std::collections::BTreeMap<String, f64> =
            std::collections::BTreeMap::new();
        let mut via_power = 0.0f64;
        for &(a, b, g) in &self.net.edges {
            let dv = u[a as usize] - u[b as usize];
            let power = g * dv * dv;
            let (pa, pb) = (
                &self.net.physical[a as usize],
                &self.net.physical[b as usize],
            );
            if pa == pb {
                *layer_power.entry(pa.clone()).or_insert(0.0) += power;
            } else {
                via_power += power;
            }
        }
        let channel_power: f64 = channels.iter().map(|c| c.power).sum();
        let p_total = layer_power.values().sum::<f64>() + via_power + channel_power;

        let mut rows: Vec<(String, f64, f64)> = layer_power
            .iter()
            .map(|(k, &p)| (k.clone(), p / it2, p))
            .collect();
        rows.push(("Vias/contacts".into(), via_power / it2, via_power));
        rows.push(("Channel (devices)".into(), channel_power / it2, channel_power));
        rows.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(std::cmp::Ordering::Equal));

        let mut s = String::new();
        s.push_str("==================== OpenRDSon resistance summary ====================\n");
        let bias: Vec<String> = self
            .terminals
            .iter()
            .map(|t| format!("{}={:.3}V", t.name, t.voltage))
            .collect();
        s.push_str(&format!("bias: {}   T={:.1} C\n", bias.join("  "), temperature));
        s.push_str("----------------------------------------------------------------------\n");
        s.push_str(&format!(
            "{:<24}{:>15}{:>15}{:>9}\n",
            "Layer", "R (ohm)", "Power (W)", "Share"
        ));
        s.push_str("----------------------------------------------------------------------\n");
        for (name, r, p) in &rows {
            let share = if p_total.abs() > 1e-300 {
                100.0 * p / p_total
            } else {
                0.0
            };
            s.push_str(&format!(
                "{:<24}{:>15.4e}{:>15.4e}{:>8.1}%\n",
                name, r, p, share
            ));
        }
        s.push_str("----------------------------------------------------------------------\n");
        let r_total = p_total / it2;
        s.push_str(&format!(
            "{:<24}{:>15.4e}{:>15.4e}{:>8.1}%\n",
            "TOTAL", r_total, p_total, 100.0
        ));
        s.push_str("----------------------------------------------------------------------\n");

        let (n0, n1) = (
            self.terminals.get(pa).map(|t| t.name.as_str()).unwrap_or("?"),
            self.terminals.get(pb).map(|t| t.name.as_str()).unwrap_or("?"),
        );
        let rdson = if i_draw.abs() > 1e-300 {
            ds_dv / i_draw
        } else {
            f64::INFINITY
        };
        s.push_str(&format!("Total Rdson ({n0}-{n1})    : {rdson:.4e} ohm\n"));
        s.push_str(&format!("Total current draw      : {i_draw:.4e} A\n"));
        s.push_str(&format!("Total power             : {p_total:.4e} W\n"));
        s.push_str("Terminal currents:\n");
        for (t, &i) in self.terminals.iter().zip(&term_i) {
            s.push_str(&format!(
                "  {:<12}{:>10.3} V {:>15.4e} A\n",
                t.name, t.voltage, i
            ));
        }
        if channels.is_empty() {
            s.push_str("Devices                 : none (resistive-only network)\n");
        } else {
            let mut rs: Vec<f64> = channels
                .iter()
                .map(|c| c.resistance)
                .filter(|r| r.is_finite())
                .collect();
            rs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            let med = rs.get(rs.len() / 2).copied().unwrap_or(f64::NAN);
            s.push_str(&format!(
                "Devices                 : {} channels, Rds {:.3e}..{:.3e} ohm (median {:.3e})\n",
                channels.len(),
                rs.first().copied().unwrap_or(f64::NAN),
                rs.last().copied().unwrap_or(f64::NAN),
                med
            ));
        }
        s.push_str(&format!(
            "Sheet network           : {} nodes, {} edges, {} cells\n",
            self.net.node_count(),
            self.net.edge_count(),
            self.net.cell_count()
        ));
        s.push_str("======================================================================\n");
        Ok(s)
    }

    pub fn network_stats(&self) -> (usize, usize, usize) {
        (
            self.net.node_count(),
            self.net.edge_count(),
            self.channels.len(),
        )
    }

    /// All terminal nodes and voltage on one named net (case-insensitive).
    fn named_terminal(&self, name: &str) -> Option<(f64, std::collections::BTreeSet<u32>)> {
        let mut nodes = std::collections::BTreeSet::new();
        let mut v = f64::NAN;
        for t in &self.terminals {
            if t.name.eq_ignore_ascii_case(name) {
                v = t.voltage;
                nodes.extend(t.nodes.iter().copied());
            }
        }
        (!nodes.is_empty()).then(|| (v, nodes))
    }

    /// Drain and source terminal groups, identified by name (`d`/`s`),
    /// falling back to the two groups with the largest voltage difference.
    /// Returns `(drain_nodes, source_nodes, |Vd - Vs|)`.
    fn ds_groups(
        &self,
    ) -> Option<(std::collections::BTreeSet<u32>, std::collections::BTreeSet<u32>, f64)> {
        if let (Some((vd, d)), Some((vs, s))) =
            (self.named_terminal("d"), self.named_terminal("s"))
        {
            return Some((d, s, (vd - vs).abs()));
        }
        let mut map: std::collections::BTreeMap<String, (f64, std::collections::BTreeSet<u32>)> =
            std::collections::BTreeMap::new();
        for t in &self.terminals {
            let e = map
                .entry(t.name.to_lowercase())
                .or_insert((t.voltage, Default::default()));
            e.0 = t.voltage;
            e.1.extend(t.nodes.iter().copied());
        }
        let groups: Vec<(f64, std::collections::BTreeSet<u32>)> = map.into_values().collect();
        let mut best: Option<(usize, usize, f64)> = None;
        for i in 0..groups.len() {
            for j in (i + 1)..groups.len() {
                let d = (groups[i].0 - groups[j].0).abs();
                if best.map(|(_, _, bd)| d > bd).unwrap_or(true) {
                    best = Some((i, j, d));
                }
            }
        }
        let (i, j, dv) = best?;
        Some((groups[i].1.clone(), groups[j].1.clone(), dv))
    }

    /// Resistance with no channel, between the drain and source terminals.
    pub fn rds_open(&self) -> Result<f64, String> {
        let (d, s, _) = self
            .ds_groups()
            .ok_or_else(|| "need drain and source terminals for rds_open".to_string())?;
        let d: Vec<u32> = d.into_iter().collect();
        let s: Vec<u32> = s.into_iter().collect();
        self.net.resistance_with_edges(&d, &s, &[], 1e-9, 20000)
    }

    /// Effective drain–source resistance at the configured terminal bias:
    /// `ΔV / I`, with `I` the total current entering the drain.
    pub fn rds(&self, temperature: f64) -> Result<f64, String> {
        let (u, channels) = self.solve_local_bias(temperature, 1e-3, 20)?;
        self.rds_from_solution(&u, &channels)
    }

    /// Drain–source resistance at gate-source voltage `vgs`, reusing the current
    /// mesh and channel set: only the gate terminal's voltage is moved, so this
    /// is far cheaper than rebuilding the network (used by the Vgs sweep).
    pub fn rds_at_vgs(&self, vgs: f64, temperature: f64) -> Result<f64, String> {
        let vs = self
            .named_terminal("s")
            .map(|(v, _)| v)
            .ok_or_else(|| "need a source terminal for rds_at_vgs".to_string())?;
        let mut terms = self.terminals.clone();
        for t in &mut terms {
            if t.name.eq_ignore_ascii_case("g") {
                t.voltage = vs + vgs;
            }
        }
        let (u, channels) =
            self.solve_local_bias_terms(&terms, temperature, 1e-3, 20, None, None)?;
        self.rds_from_solution(&u, &channels)
    }

    /// `ΔV / I` from an already-solved potential field and channel set.
    fn rds_from_solution(&self, u: &[f64], channels: &[ChannelState]) -> Result<f64, String> {
        let extra: Vec<(u32, u32, f64)> = channels
            .iter()
            .filter(|c| c.conductance > 0.0)
            .map(|c| (c.d, c.s, c.conductance))
            .collect();
        let (d_nodes, _, dv) = self
            .ds_groups()
            .ok_or_else(|| "need drain and source terminals for rds".to_string())?;
        let mut i_total = 0.0;
        for &(a, b, g) in self.net.edges.iter().chain(extra.iter()) {
            if d_nodes.contains(&a) {
                i_total += g * (u[a as usize] - u[b as usize]);
            } else if d_nodes.contains(&b) {
                i_total += g * (u[b as usize] - u[a as usize]);
            }
        }
        if i_total.abs() < 1e-300 {
            return Ok(f64::INFINITY);
        }
        Ok(dv / i_total.abs())
    }
}

fn bar_region(x0: f64, x1: f64, w: f64, t: f64, mat: u32) -> openrdson_core::geometry::SolidRegion {
    let pts = vec![
        Point::new(x0, 0.0),
        Point::new(x1, 0.0),
        Point::new(x1, w),
        Point::new(x0, w),
    ];
    extrude_polygon(&pts, 0.0, t, mat, None, None, None)
}

fn terminals_at_x(mesh: &openrdson_core::mesh::Mesh, x: f64) -> Vec<u32> {
    mesh.nodes
        .iter()
        .enumerate()
        .filter(|(_, p)| (p.x - x).abs() < 1e-12)
        .map(|(i, _)| i as u32)
        .collect()
}

fn bar_resistance(l: f64, w: f64, t: f64, sigma: f64) -> f64 {
    let model = SolidModel {
        regions: vec![bar_region(0.0, l, w, t, 1)],
        material_names: vec!["A".into()],
    };
    let (mesh, _) = mesh_model(
        &model,
        &MeshConfig {
            lateral_cell: 0.5e-6,
            z_cell: 0.1e-6,
        },
    );
    let mesh = weld_nodes(&mesh, 1e-12);
    two_terminal_resistance(
        &mesh,
        &|_| sigma,
        &terminals_at_x(&mesh, 0.0),
        &terminals_at_x(&mesh, l),
        1e-12,
        5000,
    )
    .unwrap_or(f64::NAN)
}

/// Run the generic (database-independent) acceptance checks.
pub fn run_acceptance() -> Vec<Check> {
    let mut checks = Vec::new();

    // 1. Analytic bar R = L/(σA).
    let (l, w, t, sigma) = (4e-6, 1e-6, 0.2e-6, 1e6);
    let r = bar_resistance(l, w, t, sigma);
    let analytic = l / (sigma * w * t);
    checks.push(Check {
        name: "solver/extraction: bar R = rho*L/A".into(),
        passed: (r - analytic).abs() / analytic < 1e-6,
        detail: format!("R={r:.6e} analytic={analytic:.6e}"),
    });

    // 2. Series two-material R = R1 + R2.
    let model = SolidModel {
        regions: vec![
            bar_region(0.0, 4e-6, 1e-6, 0.2e-6, 1),
            bar_region(4e-6, 8e-6, 1e-6, 0.2e-6, 2),
        ],
        material_names: vec!["A".into(), "B".into()],
    };
    let (mesh, _) = mesh_model(
        &model,
        &MeshConfig {
            lateral_cell: 0.5e-6,
            z_cell: 0.1e-6,
        },
    );
    let mesh = weld_nodes(&mesh, 1e-12);
    let rseries = two_terminal_resistance(
        &mesh,
        &|id| match id {
            1 => 1e6,
            2 => 2e6,
            _ => 0.0,
        },
        &terminals_at_x(&mesh, 0.0),
        &terminals_at_x(&mesh, 8e-6),
        1e-12,
        5000,
    )
    .unwrap_or(f64::NAN);
    let area = 1e-6 * 0.2e-6;
    let analytic_series = 4e-6 / (1e6 * area) + 4e-6 / (2e6 * area);
    checks.push(Check {
        name: "extraction: series materials add".into(),
        passed: (rseries - analytic_series).abs() / analytic_series < 1e-6,
        detail: format!("R={rseries:.6e} analytic={analytic_series:.6e}"),
    });

    checks
}

#[cfg(test)]
mod apriori_tests {
    use super::*;
    use openrdson_extraction::adaptive::AdaptiveGrid;
    use openrdson_extraction::SheetPolygon;
    use std::collections::BTreeMap;

    fn grid_10x10(base: f64) -> AdaptiveGrid {
        let poly = SheetPolygon {
            physical: "M1".into(),
            points: vec![
                Point::new(0.0, 0.0),
                Point::new(10.0, 0.0),
                Point::new(10.0, 10.0),
                Point::new(0.0, 10.0),
            ],
            component: 0,
            r_sheet: 1.0,
        };
        AdaptiveGrid::from_polygons(&[poly], base, 4000, &BTreeMap::new())
    }

    fn term(x: f64, y: f64) -> TermSpec {
        TermSpec {
            name: "T".into(),
            voltage: 0.0,
            layer: "M1".into(),
            x,
            y,
            hx: 0.0,
            hy: 0.0,
        }
    }

    #[test]
    fn terminal_refines_region_and_stays_centered() {
        let mut g = grid_10x10(1.0);
        g.max_level = 4;
        // (7.3, 7.7) lies in the NW quadrant of base cell (7,7)-(8,8); a
        // fixed-SW recursion would refine the wrong corner.
        a_priori_refine(&mut g, &[], &[term(7.3, 7.7)], 3.0, 100_000);

        // Deep refinement happened (region, not a single split).
        let deepest = g
            .leaves()
            .iter()
            .map(|&c| g.cells[c as usize].level)
            .max()
            .unwrap();
        assert_eq!(deepest, 4, "expected max_level refinement, got level {deepest}");

        // The deepest leaf contains the terminal (not the SW corner).
        let l4_contains_terminal = g.leaves().iter().any(|&c| {
            let cell = &g.cells[c as usize];
            cell.level == 4
                && 7.3 >= cell.x0
                && 7.3 <= cell.x1
                && 7.7 >= cell.y0
                && 7.7 <= cell.y1
        });
        assert!(l4_contains_terminal, "deepest leaf is not centred on the terminal");

        // More than one leaf was refined (a region, not a single cell).
        assert!(g.leaves().iter().filter(|&c| g.cells[*c as usize].level > 0).count() > 1);
    }

    #[test]
    fn terminal_in_ne_corner_refines_ne_not_sw() {
        let mut g = grid_10x10(1.0);
        g.max_level = 4;
        a_priori_refine(&mut g, &[], &[term(9.8, 9.8)], 1.0, 100_000);
        // (9.8, 9.8) is in the NE quadrant of cell (9,9)-(10,10). The level-4
        // leaf must contain it.
        assert!(g.leaves().iter().any(|&c| {
            let cell = &g.cells[c as usize];
            cell.level == 4 && 9.8 >= cell.x0 && 9.8 <= cell.x1 && 9.8 >= cell.y0 && 9.8 <= cell.y1
        }));
    }

    #[test]
    fn field_cache_roundtrip() {
        let field = SheetField {
            nodes: vec![
                (0.0, 0.0, 1.0),
                (1.0, 0.0, 2.0),
                (1.0, 1.0, 3.0),
                (0.0, 1.0, 4.0),
            ],
            quads: vec![[0, 1, 2, 3]],
            dx: vec![1.0],
            dy: vec![1.0],
            gx: vec![1.0],
            gy: vec![2.0],
            r_sheet: vec![7.0],
            layer: vec!["Metal1".into()],
            component: vec![0],
            node_component: vec![0; 4],
            terminals: vec![ResolvedTerminal {
                name: "D".into(),
                voltage: 0.0,
                layer: "Metal1".into(),
                nodes: vec![0, 1],
            }],
            total_current: 0.5,
            channels: vec![ChannelState {
                d: 0,
                s: 3,
                vgs: -2.0,
                vds: -1.5,
                resistance: 100.0,
                conductance: 0.01,
                current: 0.015,
                power: 0.0225,
                weff: 10e-6,
                bbox: Some(openrdson_core::geometry::Bbox {
                    min_x: 0.0,
                    min_y: 0.0,
                    max_x: 1.0,
                    max_y: 1.0,
                }),
                nx: 2,
                ny: 2,
            }],
        };
        let c = FieldCache {
            db_unit_meters: 1e-9,
            field,
            edges: vec![(0, 3, 0.5)],
            physical: vec!["Metal1".into(); 4],
            cell_level: vec![1],
            report: "Total Rdson (D-S) = 1.0 ohm".into(),
        };
        let c2 = FieldCache::from_bytes(&c.to_bytes()).expect("roundtrip");
        assert_eq!(c.field.nodes.len(), c2.field.nodes.len());
        assert_eq!(c.field.quads, c2.field.quads);
        assert_eq!(c.field.layer, c2.field.layer);
        assert_eq!(c.edges, c2.edges);
        assert_eq!(c.cell_level, c2.cell_level);
        assert_eq!(c.report, c2.report);
        assert!((c.field.total_current - c2.field.total_current).abs() < 1e-15);
        assert!((c.field.channels[0].resistance - c2.field.channels[0].resistance).abs() < 1e-15);
        assert!((c.field.channels[0].bbox.unwrap().max_x - c2.field.channels[0].bbox.unwrap().max_x).abs() < 1e-15);
    }
}

// ============================================================================
// Field cache: persist the solved sheet field so `viz` can re-export the
// results (KLayout/ParaView) without re-solving, and so users can post-process
// the field themselves or keep several accuracy levels on disk.
// ============================================================================

/// A solved 2.5D sheet field plus the network data the visualization export
/// needs. Persisted to `field_cache_<fingerprint>.bin` (see
/// [`field_cache_fingerprint`]) keyed by the meshing/accuracy settings.
#[derive(Debug, Clone)]
pub struct FieldCache {
    pub db_unit_meters: f64,
    pub field: SheetField,
    /// Network edges `(a, b, conductance)` (for via links).
    pub edges: Vec<(u32, u32, f64)>,
    /// Physical layer name per node (for via-link layer mapping).
    pub physical: Vec<String>,
    /// Quadtree refinement level per cell.
    pub cell_level: Vec<u8>,
    /// Pre-computed resistance report (human-readable summary).
    pub report: String,
}

impl FieldCache {
    /// Per-cell scalars (the same quantities `cell_scalars` returns), derived
    /// from the stored field without re-solving.
    pub fn cell_scalars(&self) -> Vec<(String, [(f64, f64); 4], CellScalars)> {
        let sf = &self.field;
        let it2 = (sf.total_current * sf.total_current).max(1e-300);
        let mut out = Vec::with_capacity(sf.quads.len());
        for (qi, q) in sf.quads.iter().enumerate() {
            let phys = sf.layer[qi].clone();
            let corners = [
                (sf.nodes[q[0] as usize].0, sf.nodes[q[0] as usize].1),
                (sf.nodes[q[1] as usize].0, sf.nodes[q[1] as usize].1),
                (sf.nodes[q[2] as usize].0, sf.nodes[q[2] as usize].1),
                (sf.nodes[q[3] as usize].0, sf.nodes[q[3] as usize].1),
            ];
            let rs = sf.r_sheet[qi].max(1e-30);
            let g2 = sf.gx[qi] * sf.gx[qi] + sf.gy[qi] * sf.gy[qi];
            let j = g2.sqrt() / rs;
            let power = g2 / rs * sf.dx[qi] * sf.dy[qi];
            let pot = (sf.nodes[q[0] as usize].2
                + sf.nodes[q[1] as usize].2
                + sf.nodes[q[2] as usize].2
                + sf.nodes[q[3] as usize].2)
                / 4.0;
            let hi = q.iter().map(|&n| sf.nodes[n as usize].2).fold(f64::NEG_INFINITY, f64::max);
            let lo = q.iter().map(|&n| sf.nodes[n as usize].2).fold(f64::INFINITY, f64::min);
            let dv = (hi - lo).abs();
            let current = if dv > 1e-15 { power / dv } else { 0.0 };
            out.push((
                phys,
                corners,
                CellScalars {
                    potential: pot,
                    current_density: j,
                    power,
                    current,
                    resistance: power / it2,
                    level: self.cell_level.get(qi).copied().unwrap_or(0),
                    net: sf.component[qi],
                },
            ));
        }
        out
    }

    /// Via connections `(x, y, bottom_node, top_node, conductance)`, from the
    /// stored network edges + node positions.
    pub fn via_links(&self) -> Vec<(f64, f64, u32, u32, f64)> {
        let mut out = Vec::new();
        for &(a, b, g) in &self.edges {
            if self.physical[a as usize] != self.physical[b as usize] {
                let x = (self.field.nodes[a as usize].0 + self.field.nodes[b as usize].0) / 2.0;
                let y = (self.field.nodes[a as usize].1 + self.field.nodes[b as usize].1) / 2.0;
                out.push((x, y, a, b, g));
            }
        }
        out
    }

    pub fn save(&self, path: &std::path::Path) -> Result<(), String> {
        std::fs::write(path, self.to_bytes()).map_err(|e| e.to_string())
    }

    pub fn load(path: &std::path::Path) -> Result<Self, String> {
        let bytes = std::fs::read(path).map_err(|e| e.to_string())?;
        Self::from_bytes(&bytes)
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        let mut b: Vec<u8> = Vec::new();
        let mut w = |bytes: &[u8]| b.extend_from_slice(bytes);
        w(b"RDSONFC1");
        w(&1u32.to_le_bytes());
        w(&self.db_unit_meters.to_le_bytes());
        let f = &self.field;
        w(&(f.nodes.len() as u64).to_le_bytes());
        for &(x, y, u) in &f.nodes {
            w(&x.to_le_bytes());
            w(&y.to_le_bytes());
            w(&u.to_le_bytes());
        }
        w(&(f.quads.len() as u64).to_le_bytes());
        for q in &f.quads {
            for &n in q {
                w(&n.to_le_bytes());
            }
        }
        for v in [&f.dx, &f.dy, &f.gx, &f.gy, &f.r_sheet] {
            w(&(v.len() as u64).to_le_bytes());
            for &x in v.iter() {
                w(&x.to_le_bytes());
            }
        }
        w(&(f.layer.len() as u64).to_le_bytes());
        for s in &f.layer {
            w(&(s.len() as u64).to_le_bytes());
            w(s.as_bytes());
        }
        for v in [&f.component, &f.node_component] {
            w(&(v.len() as u64).to_le_bytes());
            for &x in v.iter() {
                w(&(x as u64).to_le_bytes());
            }
        }
        w(&(f.terminals.len() as u64).to_le_bytes());
        for t in &f.terminals {
            w(&(t.name.len() as u64).to_le_bytes());
            w(t.name.as_bytes());
            w(&t.voltage.to_le_bytes());
            w(&(t.layer.len() as u64).to_le_bytes());
            w(t.layer.as_bytes());
            w(&(t.nodes.len() as u64).to_le_bytes());
            for &n in &t.nodes {
                w(&n.to_le_bytes());
            }
        }
        w(&f.total_current.to_le_bytes());
        w(&(f.channels.len() as u64).to_le_bytes());
        for c in &f.channels {
            w(&c.d.to_le_bytes());
            w(&c.s.to_le_bytes());
            w(&c.vgs.to_le_bytes());
            w(&c.vds.to_le_bytes());
            w(&c.resistance.to_le_bytes());
            w(&c.conductance.to_le_bytes());
            w(&c.current.to_le_bytes());
            w(&c.power.to_le_bytes());
            w(&c.weff.to_le_bytes());
            match c.bbox {
                Some(bb) => {
                    w(&[1u8]);
                    w(&bb.min_x.to_le_bytes());
                    w(&bb.min_y.to_le_bytes());
                    w(&bb.max_x.to_le_bytes());
                    w(&bb.max_y.to_le_bytes());
                }
                None => w(&[0u8]),
            }
            w(&(c.nx as u64).to_le_bytes());
            w(&(c.ny as u64).to_le_bytes());
        }
        w(&(self.edges.len() as u64).to_le_bytes());
        for &(a, bb, g) in &self.edges {
            w(&a.to_le_bytes());
            w(&bb.to_le_bytes());
            w(&g.to_le_bytes());
        }
        w(&(self.physical.len() as u64).to_le_bytes());
        for s in &self.physical {
            w(&(s.len() as u64).to_le_bytes());
            w(s.as_bytes());
        }
        w(&(self.cell_level.len() as u64).to_le_bytes());
        w(&self.cell_level);
        w(&(self.report.len() as u64).to_le_bytes());
        w(self.report.as_bytes());
        b
    }

    pub fn from_bytes(bytes: &[u8]) -> Result<Self, String> {
        let mut r = Rdr::new(bytes);
        if r.take(8)? != b"RDSONFC1" {
            return Err("field cache: bad magic".into());
        }
        let version = r.u32()?;
        if version != 1 {
            return Err(format!("field cache: unsupported version {version}"));
        }
        let db_unit_meters = r.f64()?;
        let n_nodes = r.u64()? as usize;
        let mut nodes = Vec::with_capacity(n_nodes);
        for _ in 0..n_nodes {
            nodes.push((r.f64()?, r.f64()?, r.f64()?));
        }
        let n_quads = r.u64()? as usize;
        let mut quads = Vec::with_capacity(n_quads);
        for _ in 0..n_quads {
            quads.push([r.u32()?, r.u32()?, r.u32()?, r.u32()?]);
        }
        let dx = r.vec_f64()?;
        let dy = r.vec_f64()?;
        let gx = r.vec_f64()?;
        let gy = r.vec_f64()?;
        let r_sheet = r.vec_f64()?;
        let layer = r.vec_str()?;
        let component = r.vec_usize()?;
        let node_component = r.vec_usize()?;
        let n_terms = r.u64()? as usize;
        let mut terminals = Vec::with_capacity(n_terms);
        for _ in 0..n_terms {
            let name = r.str()?;
            let voltage = r.f64()?;
            let layer_t = r.str()?;
            let n_nodes = r.u64()? as usize;
            let mut tnodes = Vec::with_capacity(n_nodes);
            for _ in 0..n_nodes {
                tnodes.push(r.u32()?);
            }
            terminals.push(ResolvedTerminal {
                name,
                voltage,
                layer: layer_t,
                nodes: tnodes,
            });
        }
        let total_current = r.f64()?;
        let n_channels = r.u64()? as usize;
        let mut channels = Vec::with_capacity(n_channels);
        for _ in 0..n_channels {
            let d = r.u32()?;
            let s = r.u32()?;
            let vgs = r.f64()?;
            let vds = r.f64()?;
            let resistance = r.f64()?;
            let conductance = r.f64()?;
            let current = r.f64()?;
            let power = r.f64()?;
            let weff = r.f64()?;
            let has_bbox = r.u8()? == 1;
            let bbox = if has_bbox {
                Some(openrdson_core::geometry::Bbox {
                    min_x: r.f64()?,
                    min_y: r.f64()?,
                    max_x: r.f64()?,
                    max_y: r.f64()?,
                })
            } else {
                None
            };
            let nx = r.u64()? as usize;
            let ny = r.u64()? as usize;
            channels.push(ChannelState {
                d,
                s,
                vgs,
                vds,
                resistance,
                conductance,
                current,
                power,
                weff,
                bbox,
                nx,
                ny,
            });
        }
        let n_edges = r.u64()? as usize;
        let mut edges = Vec::with_capacity(n_edges);
        for _ in 0..n_edges {
            edges.push((r.u32()?, r.u32()?, r.f64()?));
        }
        let physical = r.vec_str()?;
        let n_levels = r.u64()? as usize;
        let cell_level = r.take(n_levels)?.to_vec();
        let report = r.str()?;
        Ok(FieldCache {
            db_unit_meters,
            field: SheetField {
                nodes,
                quads,
                dx,
                dy,
                gx,
                gy,
                r_sheet,
                layer,
                component,
                node_component,
                terminals,
                total_current,
                channels,
            },
            edges,
            physical,
            cell_level,
            report,
        })
    }
}

/// A minimal byte cursor for the field-cache format.
struct Rdr<'a> {
    b: &'a [u8],
    i: usize,
}

impl<'a> Rdr<'a> {
    fn new(b: &'a [u8]) -> Self {
        Self { b, i: 0 }
    }
    fn take(&mut self, n: usize) -> Result<&'a [u8], String> {
        if self.i + n > self.b.len() {
            return Err("field cache: truncated".into());
        }
        let s = &self.b[self.i..self.i + n];
        self.i += n;
        Ok(s)
    }
    fn u8(&mut self) -> Result<u8, String> {
        Ok(self.take(1)?[0])
    }
    fn u32(&mut self) -> Result<u32, String> {
        Ok(u32::from_le_bytes(self.take(4)?.try_into().unwrap()))
    }
    fn u64(&mut self) -> Result<u64, String> {
        Ok(u64::from_le_bytes(self.take(8)?.try_into().unwrap()))
    }
    fn f64(&mut self) -> Result<f64, String> {
        Ok(f64::from_le_bytes(self.take(8)?.try_into().unwrap()))
    }
    fn str(&mut self) -> Result<String, String> {
        let n = self.u64()? as usize;
        Ok(String::from_utf8_lossy(self.take(n)?).into_owned())
    }
    fn vec_f64(&mut self) -> Result<Vec<f64>, String> {
        let n = self.u64()? as usize;
        let mut v = Vec::with_capacity(n);
        for _ in 0..n {
            v.push(self.f64()?);
        }
        Ok(v)
    }
    fn vec_usize(&mut self) -> Result<Vec<usize>, String> {
        let n = self.u64()? as usize;
        let mut v = Vec::with_capacity(n);
        for _ in 0..n {
            v.push(self.u64()? as usize);
        }
        Ok(v)
    }
    fn vec_str(&mut self) -> Result<Vec<String>, String> {
        let n = self.u64()? as usize;
        let mut v = Vec::with_capacity(n);
        for _ in 0..n {
            v.push(self.str()?);
        }
        Ok(v)
    }
}

impl SheetFullArray {
    /// Solve and capture the field for caching / viz post-processing. This
    /// performs the (uniform or adaptive) sheet solve and packages the result
    /// into a [`FieldCache`] that can be re-exported without re-solving.
    pub fn field_cache(&self, temperature: f64) -> Result<FieldCache, String> {
        let field = self.sheet_fields(temperature)?;
        let report = self.resistance_report(temperature).unwrap_or_default();
        Ok(FieldCache {
            db_unit_meters: self.db_unit_meters,
            field,
            edges: self.net.edges.clone(),
            physical: self.net.physical.clone(),
            cell_level: self.net.cell_level.clone(),
            report,
        })
    }
}

/// A deterministic, filename-safe fingerprint of the accuracy/bias settings, so
/// the field cache can keep one entry per accuracy level.
/// Bump whenever the *meaning* of the cached data changes — e.g. the
/// resistance-report computation — so a stale cache is re-solved instead of
/// replayed. The binary format version inside [`FieldCache`] tracks layout
/// changes; this constant tracks content-semantics changes.
pub const FIELD_CACHE_CONTENT_VERSION: u32 = 3;

pub fn field_cache_fingerprint(settings: &MeshSettings) -> String {
    let mut s = format!(
        "v{FIELD_CACHE_CONTENT_VERSION}|cell={:.6e}|cap={}|via={:.6e}|temp={:.3}|model={:?}|dbu={:?}",
        settings.sheet_cell_m,
        settings.sheet_cap,
        settings.via_group_radius_m,
        settings.temperature_c,
        settings.missing_model,
        settings.db_unit_nm
    );
    for (l, c) in &settings.per_layer_cell_m {
        s.push_str(&format!("|{l}={:.6e}", c));
    }
    let a = &settings.adaptive;
    s.push_str(&format!(
        "|adapt={}:{}:{}:{:.4}:{:.6e}",
        a.enabled, a.max_level, a.iters, a.mark_frac, a.refine_near_vias_m
    ));
    for t in &settings.terminals {
        s.push_str(&format!(
            "|t{}:{:.4}:{:.4}:{:.4}:{:.4}:{:.4}:{}",
            t.name,
            t.voltage,
            t.x_um,
            t.y_um,
            t.dx_um,
            t.dy_um,
            t.layer.as_deref().unwrap_or("")
        ));
    }
    // Input database files: path + size + mtime, so switching databases or
    // regenerating a file (model.csv, layout, ICT, ...) invalidates the cache.
    // A missing file still contributes its path, so pointing at a different
    // (empty) database also changes the fingerprint.
    for p in settings.paths.input_files() {
        let meta = std::fs::metadata(&p)
            .ok()
            .and_then(|m| m.modified().ok().map(|mt| (m.len(), mt)));
        s.push_str(&format!("|f{}={:?}", p.display(), meta));
    }
    let mut h: u64 = 0xcbf29ce484222325;
    for b in s.bytes() {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    format!("{:016x}", h)
}

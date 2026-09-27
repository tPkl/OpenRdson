//! OpenRDSon pipeline driver.
//!
//! Usage: `openrdson [FLAGS] <command>`
//!
//! Commands: `sheet-rds`, `viz`, `all`.
//! Every extraction command reads its input database from the YAML config
//! (`--config <file>`); see `--print-default-config` for the template.
//!
//! Logging flags: `--verbose`/`-v` (DEBUG), `--trace` (TRACE), `--quiet`/`-q`
//! (WARN). Default is INFO.
//!
//! Adaptive (quadtree) 2.5D meshing: `--adaptive` enables the refine->solve
//! loop; `--adaptive-level <n>` (max quadtree level), `--adaptive-iters <n>`,
//! `--adaptive-mark <0..1>` (Doerfler fraction), `--refine-near-vias <um>`
//! (a-priori geometry refinement radius). Config keys mirror these under
//! `meshing:` (`adaptive`, `adaptive_max_level`, ...).
//!
//! `viz` writes KLayout GDS colormaps and VTU/ParaView files into the
//! configured output directory; the solved field is cached so a re-export does
//! not re-solve.

use openrdson_geometry::assemble_stack;
use openrdson_techfile::{load_tech_stack, read_sft_cci_map};
use std::collections::BTreeMap;

fn print_help() {
    println!(
        "OpenRDSon — on-resistance (Rdson) extraction pipeline

USAGE:
    openrdson [FLAGS] <command>

COMMANDS:
    all           Rdson solve + table + viz export [default]
    sheet-rds     Full-array 2.5D sheet extraction; prints the Rds(Vgs) table
    viz           Export the solved field (GDS/KLayout + VTU/ParaView); cached
                  to `field_cache_<hash>.bin`, re-solves only when stale

CONFIGURATION:
    --config <file>           YAML config (default: ./openrdson.yaml)
    --print-default-config    Emit the default config template and exit

LOGGING:
    -v, --verbose             DEBUG logging
    -q, --quiet               WARN logging (errors/warnings only)
        --trace               TRACE logging

MESHING ACCURACY (R3D-style):
    --sheet-cell <um>         Full-array sheet cell size (um)
    --sheet-cap <n>           Max cells per rectangle dimension
    --via-radius <um>         Group same-net vias within a radius (0 = off)

ADAPTIVE (QUADTREE) MESHING:
    --adaptive                Enable the refine→solve loop
    --adaptive-level <n>      Max quadtree refinement level
    --adaptive-iters <n>      Max refine→solve iterations
    --adaptive-mark <f>       Doerfler marking fraction (0..1)
    --refine-near-vias <um>   A-priori refinement radius around terminals/vias

VISUALIZATION:
    --viz-zscale <f>          Vertical exaggeration for the 2.5D stack
    --out-dbu-nm <n>          Output GDS database unit (nm)

EXAMPLES:
    openrdson --config openrdson.yaml sheet-rds
    openrdson --config openrdson.yaml sheet-rds --adaptive --sheet-cell 4
    openrdson --config openrdson.yaml viz --adaptive
    openrdson --print-default-config > openrdson.yaml"
    );
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();

    // Print usage and exit.
    if args.iter().any(|a| a == "-h" || a == "--help" || a == "help") {
        print_help();
        return;
    }

    // Emit the documented configuration template and exit.
    if args.iter().any(|a| a == "--print-default-config") {
        print!("{}", openrdson_config::SolverConfig::default_yaml());
        return;
    }

    let config = load_config(&args);

    // Logging: CLI flags override the config; default INFO.
    let cfg_level = config
        .as_ref()
        .map(|c| c.logging.level)
        .unwrap_or(openrdson_config::LogLevel::Info);
    let level = if args.iter().any(|a| a == "--quiet" || a == "-q") {
        openrdson_core::log::WARN
    } else if args.iter().any(|a| a == "--trace") {
        openrdson_core::log::TRACE
    } else if args.iter().any(|a| a == "--verbose" || a == "-v") {
        openrdson_core::log::DEBUG
    } else {
        match cfg_level {
            openrdson_config::LogLevel::Error => openrdson_core::log::ERROR,
            openrdson_config::LogLevel::Warn => openrdson_core::log::WARN,
            openrdson_config::LogLevel::Info => openrdson_core::log::INFO,
            openrdson_config::LogLevel::Debug => openrdson_core::log::DEBUG,
            openrdson_config::LogLevel::Trace => openrdson_core::log::TRACE,
        }
    };
    openrdson_core::log::set_level(level);

    const VALUE_FLAGS: [&str; 12] = [
        "--mesh-config",
        "--via-radius",
        "--sheet-cell",
        "--sheet-cap",
        "--db-unit-nm",
        "--out-dbu-nm",
        "--viz-zscale",
        "--config",
        "--adaptive-level",
        "--adaptive-iters",
        "--adaptive-mark",
        "--refine-near-vias",
    ];
    let mut cmd = None;
    let mut i = 0;
    while i < args.len() {
        let a = &args[i];
        if VALUE_FLAGS.contains(&a.as_str()) {
            i += 2;
            continue;
        }
        if a.starts_with('-') {
            i += 1;
            continue;
        }
        cmd = Some(a.clone());
        break;
    }
    let cmd = cmd.unwrap_or_else(|| "all".to_string());
    let settings = mesh_settings(&args, config.as_ref());
    if let Some(c) = &config {
        openrdson_core::log_info!(
            "config: cci={} out={} sheet_cell={:.3}um z_scale={} threads={}",
            c.project.cci_dir.display(),
            c.project.output_dir.display(),
            c.meshing.sheet_cell_um,
            c.visualization.z_scale,
            c.performance.threads
        );
    }
    match cmd.as_str() {
        "sheet-rds" => {
            let _ = run_sheet_rds(&settings);
        }
        "viz" => run_viz(&settings, config.as_ref(), None),
        "all" | _ => run_all(&settings, config.as_ref()),
    }
}

/// Load the YAML config named by `--config <file>`, if present.
fn load_config(args: &[String]) -> Option<openrdson_config::SolverConfig> {
    let i = args.iter().position(|a| a == "--config")?;
    let Some(p) = args.get(i + 1) else {
        eprintln!("--config requires a file path");
        std::process::exit(2);
    };
    match std::fs::read_to_string(p) {
        Ok(text) => match openrdson_config::SolverConfig::from_yaml(&text) {
            Ok(c) => Some(c),
            Err(e) => {
                eprintln!("config error in {p}: {e}");
                std::process::exit(2);
            }
        },
        Err(e) => {
            eprintln!("cannot read config {p}: {e}");
            std::process::exit(2);
        }
    }
}

/// Build `MeshSettings` from the config, then apply CLI flag overrides.
fn mesh_settings(
    args: &[String],
    config: Option<&openrdson_config::SolverConfig>,
) -> openrdson_validation::MeshSettings {
    let mut s = openrdson_validation::MeshSettings::default();
    if let Some(c) = config {
        s.sheet_cell_m = c.meshing.sheet_cell_um * 1e-6;
        s.sheet_cap = c.meshing.sheet_cap;
        s.via_group_radius_m = c.meshing.via_group_radius_um * 1e-6;
        s.per_layer_cell_m = c
            .meshing
            .per_layer_cell_um
            .iter()
            .map(|(k, v)| (k.clone(), v * 1e-6))
            .collect();
        s.db_unit_nm = c.units.db_unit_nm;
        s.paths.cci_dir = c.project.cci_dir.clone();
        s.paths.layout = c.project.layout.clone();
        s.paths.gds_map = c.project.gds_map.clone();
        s.paths.ports = c.project.ports.clone();
        s.paths.devtab = c.project.devtab.clone();
        s.paths.spi = c.project.spi.clone();
        s.paths.tech_ict = c.project.tech_ict.clone();
        s.paths.layer_map = c.project.layer_map.clone();
        s.paths.model_csv = c.project.model_csv.clone();
        s.terminals = c
            .terminals
            .iter()
            .map(|t| openrdson_validation::BiasTerminal {
                name: t.name.clone(),
                voltage: t.voltage,
                x_um: t.x_um,
                y_um: t.y_um,
                dx_um: t.dx_um,
                dy_um: t.dy_um,
                layer: t.layer.clone(),
            })
            .collect();
        s.missing_model = match c.solver.missing_model {
            openrdson_config::MissingModel::Error => openrdson_validation::MissingModel::Error,
            openrdson_config::MissingModel::Open => openrdson_validation::MissingModel::Open,
        };
        s.temperature_c = c.solver.temperature_c;
        s.adaptive = openrdson_extraction::adaptive::AdaptiveConfig {
            enabled: c.meshing.adaptive,
            tol: c.meshing.adaptive_tol,
            max_level: c.meshing.adaptive_max_level,
            max_nodes: c.meshing.adaptive_max_nodes,
            mark_frac: c.meshing.adaptive_mark_frac,
            iters: c.meshing.adaptive_iters,
            refine_near_vias_m: c.meshing.refine_near_vias_um * 1e-6,
        };
    }
    let get = |flag: &str| -> Option<String> {
        args.iter()
            .position(|a| a == flag)
            .and_then(|i| args.get(i + 1))
            .cloned()
    };
    if let Some(v) = get("--via-radius").and_then(|v| v.parse::<f64>().ok()) {
        s.via_group_radius_m = v * 1e-6;
    }
    if let Some(v) = get("--sheet-cell").and_then(|v| v.parse::<f64>().ok()) {
        s.sheet_cell_m = v * 1e-6;
    }
    if let Some(v) = get("--sheet-cap").and_then(|v| v.parse::<usize>().ok()) {
        s.sheet_cap = v;
    }
    if let Some(v) = get("--db-unit-nm").and_then(|v| v.parse::<f64>().ok()) {
        s.db_unit_nm = Some(v);
    }
    if let Some(p) = get("--mesh-config") {
        if let Err(e) = s.load_config_file(std::path::Path::new(&p)) {
            openrdson_core::log_warn!("could not read mesh config {p}: {e}");
        }
    }
    if args.iter().any(|a| a == "--adaptive") {
        s.adaptive.enabled = true;
    }
    if let Some(v) = get("--adaptive-level").and_then(|v| v.parse::<u8>().ok()) {
        s.adaptive.max_level = v;
        s.adaptive.enabled = true;
    }
    if let Some(v) = get("--adaptive-iters").and_then(|v| v.parse::<usize>().ok()) {
        s.adaptive.iters = v;
        s.adaptive.enabled = true;
    }
    if let Some(v) = get("--adaptive-mark").and_then(|v| v.parse::<f64>().ok()) {
        s.adaptive.mark_frac = v;
        s.adaptive.enabled = true;
    }
    if let Some(v) = get("--refine-near-vias").and_then(|v| v.parse::<f64>().ok()) {
        s.adaptive.refine_near_vias_m = v * 1e-6;
        s.adaptive.enabled = true;
    }
    s
}


fn run_sheet_rds(
    settings: &openrdson_validation::MeshSettings,
) -> Option<openrdson_validation::SheetFullArray> {
    let _section = openrdson_core::log::SectionGuard::new("sheet-rds");
    log_databases(settings);
    let temperature = settings.temperature_c;

    // Build the base (uniform) array once.
    let base = match openrdson_validation::SheetFullArray::build_with(settings) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("sheet-rds failed: {e}");
            return None;
        }
    };

    // Adaptive refinement: solve on the refined quadtree mesh and report the
    // *refined* result. The refined array replaces the base array for the
    // resistance summary, the Vgs sweep and the visualization export, so the
    // reported numbers all come from the adaptive mesh (not the base grid).
    let fa = if settings.adaptive.enabled {
        match base.adaptive_solve_full(&settings.adaptive, temperature) {
            Ok((rep, refined)) => {
                // The first history entry is the base (uniform) grid solve, so
                // it doubles as the "before refinement" reference without an
                // extra solve.
                let uniform = rep.history.first().map(|h| h.1).unwrap_or(f64::NAN);
                println!(
                    "adaptive: Rds={:.6e} ohm  nodes={} edges={} cells={} iters={} levels={:?}",
                    rep.rds, rep.nodes, rep.edges, rep.cells, rep.iters, rep.levels
                );
                println!("  uniform (base grid) Rds={uniform:.6e} ohm");
                for (n, r, eta) in &rep.history {
                    println!("   iter: nodes={n:>7} Rds={r:>12.6e} indicator={eta:.6e}");
                }
                refined
            }
            Err(e) => {
                eprintln!("adaptive failed: {e}");
                base
            }
        }
    } else {
        base
    };

    // Summary at the configured bias (on the refined mesh when adaptive).
    let (nodes, edges, channels) = fa.network_stats();
    println!("sheet-rds: {nodes} nodes / {edges} edges / {channels} channel connections");
    println!("  open (no channel) = {:?}", fa.rds_open());
    if let Ok(report) = fa.resistance_report(temperature) {
        print!("{report}");
    }

    // Vgs sweep: re-use the same mesh, only moving the gate bias, so each step
    // is a cheap re-solve instead of a full network rebuild.
    let has_g = settings
        .terminals
        .iter()
        .any(|t| t.name.eq_ignore_ascii_case("g"));
    let has_s = settings
        .terminals
        .iter()
        .any(|t| t.name.eq_ignore_ascii_case("s"));
    if !(has_g && has_s) {
        return Some(fa);
    }
    println!("\nVgs sweep (Rdson vs gate-source voltage):");
    println!("  {:>6}  {:>14}", "Vgs", "Rds(ohm)");
    for vgs in [0.0, -0.5, -1.0, -1.5, -2.0, -2.5] {
        match fa.rds_at_vgs(vgs, temperature) {
            Ok(r) => println!("  {vgs:>6.2}  {r:>14.4e}"),
            Err(e) => eprintln!("sheet-rds failed: {e}"),
        }
    }
    Some(fa)
}

fn to_du(v: f64, db_unit_meters: f64) -> i32 {
    (v / db_unit_meters).round() as i32
}

/// Make a physical-layer name safe for a file name.
fn sanitize(name: &str) -> String {
    name.chars()
        .map(|c| if c.is_ascii_alphanumeric() { c } else { '_' })
        .collect()
}

/// Absolute `(z_bottom, z_top)` of a physical layer, in meters, from the ICT
/// stack. Vias span from the top of their bottom layer to the bottom of their
/// top layer.
fn layer_z(stack: &openrdson_core::tech::TechStack, name: &str) -> Option<(f64, f64)> {
    if let Some(c) = stack.conductor(name) {
        if let (Some(h), Some(t)) = (c.height, c.thickness) {
            return Some((h, h + t));
        }
    }
    if let Some(d) = stack.diffusion(name) {
        if let (Some(h), Some(t)) = (d.height, d.thickness) {
            return Some((h, h + t));
        }
    }
    if let Some(v) = stack.via(name) {
        let zb = v.bottom_layer.as_deref().and_then(|l| layer_z(stack, l)).map(|(_, z1)| z1);
        let zt = v.top_layer.as_deref().and_then(|l| layer_z(stack, l)).map(|(z0, _)| z0);
        if let (Some(zb), Some(zt)) = (zb, zt) {
            return Some((zb, zt));
        }
    }
    None
}

/// Conductor/diffusion thickness (meters) for a physical layer.
fn layer_thickness(stack: &openrdson_core::tech::TechStack, name: &str) -> Option<f64> {
    if let Some(c) = stack.conductor(name) {
        return c.thickness;
    }
    if let Some(d) = stack.diffusion(name) {
        return d.thickness;
    }
    None
}

/// Load the tech stack and every extruded conducting/via region, for the 2.5D
/// visualization (correct per-polygon via geometry, unlike the 5 µm voxel mesh).
#[allow(clippy::type_complexity)]
/// Describe an ICT `rho width` resistivity table compactly (Ω/sq over the
/// width range it covers).
fn rho_table_summary(vals: &[f64]) -> String {
    if vals.is_empty() {
        return "-".to_string();
    }
    if vals.len() < 2 {
        return format!("{:.4} ohm/sq", vals[0]);
    }
    let mut pts: Vec<(f64, f64)> = vals
        .chunks(2)
        .filter(|c| c.len() == 2)
        .map(|c| (c[1], c[0])) // (width_um, rho)
        .collect();
    if pts.is_empty() {
        return format!("{:.4} ohm/sq", vals[0]);
    }
    pts.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    format!(
        "{:.4}..{:.4} ohm/sq @ {:.2}..{:.2}um",
        pts[0].1,
        pts[pts.len() - 1].1,
        pts[0].0,
        pts[pts.len() - 1].0
    )
}

/// Log the extracted ICT process stack (conductors, diffusions, vias) and the
/// CCI logical->physical layer map: the resistances the solver will use.
fn log_ict_summary(tech: &openrdson_core::tech::TechStack, map: &openrdson_core::tech::LayerMap) {
    openrdson_core::log_info!(
        "ICT: process={} temp_ref={:.1}C  {} conductors / {} diffusions / {} vias / {} dielectrics",
        tech.process,
        tech.temp_ref(),
        tech.conductors.len(),
        tech.diffusions.len(),
        tech.vias.len(),
        tech.dielectrics.len()
    );
    for c in openrdson_techfile::conductors_by_height(tech) {
        let z = match (c.height, c.thickness) {
            (Some(h), Some(t)) => format!("{:.3}-{:.3}", h * 1e6, (h + t) * 1e6),
            (Some(h), None) => format!("{:.3}-?", h * 1e6),
            _ => "?-?".to_string(),
        };
        let t = c
            .thickness
            .map(|t| format!("{:.3}", t * 1e6))
            .unwrap_or_else(|| "-".to_string());
        openrdson_core::log_info!(
            "  conductor {:<8} z={:<13} um t={:<7} um Rs={:<30} tc1={:.4} tc2={:.4}",
            c.name,
            z,
            t,
            rho_table_summary(&c.resistivity),
            c.temp_tc1.unwrap_or(0.0),
            c.temp_tc2.unwrap_or(0.0)
        );
    }
    for d in &tech.diffusions {
        let t = d
            .thickness
            .map(|t| format!("{:.3}", t * 1e6))
            .unwrap_or_else(|| "-".to_string());
        openrdson_core::log_info!(
            "  diffusion {:<8} t={:<7} um Rs={}",
            d.name,
            t,
            rho_table_summary(&d.resistivity)
        );
    }
    for v in &tech.vias {
        let (r_ref, a_ref) = match v.area_resistance.as_slice() {
            [r, a, ..] => (*r, *a),
            [r] => (*r, 0.0),
            _ => (f64::NAN, 0.0),
        };
        let ra = if a_ref > 0.0 { r_ref * a_ref } else { 0.0 };
        openrdson_core::log_info!(
            "  via {:<10} {} -> {:<8} R_ref={:.4} ohm @ A_ref={:.4} um^2  R*A={:.4} ohm.um^2  tc1={:.4}",
            v.name,
            v.bottom_layer.as_deref().unwrap_or("?"),
            v.top_layer.as_deref().unwrap_or("?"),
            r_ref,
            a_ref,
            ra,
            v.temp_tc1.unwrap_or(0.0)
        );
    }
    for (logical, physical) in &map.conducting {
        openrdson_core::log_debug!("  CCI map conductor {logical} -> {physical}");
    }
    for (logical, physical) in &map.via {
        openrdson_core::log_debug!("  CCI map via {logical} -> {physical}");
    }
}

/// Log the CCI connectivity data: net names, device templates, and the golden
/// SPI instances (per-instance finger width and summed channel width).
fn log_cci_summary(paths: &openrdson_validation::ProjectPaths) {
    let nets = openrdson_layout_db::read_net_names(paths.net_names_path()).unwrap_or_default();
    let templates =
        openrdson_layout_db::read_cci_device_templates(paths.pin_xy_path()).unwrap_or_default();
    let spi = openrdson_layout_db::read_spi_instances(paths.spi_path()).unwrap_or_default();
    openrdson_core::log_info!(
        "CCI: {} net names / {} device templates / {} SPI instances",
        nets.len(),
        templates.len(),
        spi.len()
    );
    for (id, name) in &nets {
        openrdson_core::log_info!("  CCI net {id} -> {name}");
    }
    let mut per_model: BTreeMap<String, usize> = BTreeMap::new();
    let mut total_w = 0.0;
    let mut total_weff = 0.0;
    for inst in &spi {
        *per_model.entry(inst.model.clone()).or_insert(0) += 1;
        total_w += inst.params.get("w").copied().unwrap_or(0.0);
        total_weff += inst.params.get("weff").copied().unwrap_or(0.0);
    }
    for (model, n) in &per_model {
        openrdson_core::log_info!("  CCI SPI model {model}: {n} instances");
    }
    openrdson_core::log_info!(
        "  CCI SPI channel width: sum(w)={:.3} um (used)  sum(weff)={:.3} um (combined, not used)",
        total_w * 1e6,
        total_weff * 1e6
    );
    for (id, t) in templates.iter().take(8) {
        let terms: Vec<String> = t
            .terminals
            .iter()
            .map(|(term, layer)| format!("{term}({layer})"))
            .collect();
        openrdson_core::log_info!(
            "  CCI template {id}: model={} seed={} terms=[{}]",
            t.model,
            t.seed,
            terms.join(" ")
        );
    }
}

fn load_tech_and_regions(
    settings: &openrdson_validation::MeshSettings,
) -> Option<(
    openrdson_core::tech::TechStack,
    Vec<openrdson_core::geometry::SolidRegion>,
)> {
    // Use the same database-unit correction as the sheet network, or the vias
    // and metals end up on different scales.
    let layout = openrdson_validation::load_layout(settings).ok()?;
    let tech = load_tech_stack(&settings.paths.tech_ict).ok()?;
    let map = read_sft_cci_map(&settings.paths.layer_map).ok()?;
    let (solid, _) = assemble_stack(&layout, &tech, &map);
    Some((tech, solid.regions))
}

/// Load the ICT stack + CCI layer map and log the extracted database summary
/// (the resistances and connectivity the solver will use).
fn log_databases(settings: &openrdson_validation::MeshSettings) {
    use openrdson_core::log_warn;
    match load_tech_stack(&settings.paths.tech_ict) {
        Ok(tech) => match read_sft_cci_map(&settings.paths.layer_map) {
            Ok(map) => {
                log_ict_summary(&tech, &map);
                log_cci_summary(&settings.paths);
            }
            Err(e) => log_warn!("layer map {}: {e}", settings.paths.layer_map.display()),
        },
        Err(e) => log_warn!("ICT {}: {e}", settings.paths.tech_ict.display()),
    }
}

/// Shoelace area of a polygon footprint.
fn polygon_area(pts: &[(f64, f64)]) -> f64 {
    let n = pts.len();
    if n < 3 {
        return 0.0;
    }
    let mut a = 0.0;
    for i in 0..n {
        let (x0, y0) = pts[i];
        let (x1, y1) = pts[(i + 1) % n];
        a += x0 * y1 - x1 * y0;
    }
    (a * 0.5).abs()
}

/// Build the 2.5D field stack: each layer's sheet-extraction quads extruded to
/// its ICT z-range, plus every via polygon extruded between the two layers it
/// connects.
///
/// Cell data: `current_density_A_m2`, `power_W`, `current_A`,
/// `resistance_ohm` (per-element contribution, `power / I_total²`),
/// `via_current_A`, `net_id`, `layer_id`. Returns the mesh and the layer table.
#[allow(clippy::type_complexity)]
fn build_2p5d_stack(
    tech: &openrdson_core::tech::TechStack,
    via_regions: &[openrdson_core::geometry::SolidRegion],
    links: &[(f64, f64, u32, u32, f64)],
    cell_level: &[u8],
    sf: &openrdson_validation::SheetField,
    z_scale: f64,
) -> (openrdson_io::VtuMesh, Vec<String>) {
    use openrdson_core::geometry::Solid;
    use openrdson_io::{VtkCellType, VtuMesh};

    let i_total = sf.total_current.max(1e-300);
    let mut all_layers: Vec<String> = Vec::new();
    let layer_id = |name: &str, all: &mut Vec<String>| -> f64 {
        match all.iter().position(|l| l == name) {
            Some(i) => i as f64,
            None => {
                all.push(name.to_string());
                (all.len() - 1) as f64
            }
        }
    };

    let mut points: Vec<[f64; 3]> = Vec::new();
    let mut potentials: Vec<f64> = Vec::new();
    let mut cells: Vec<(VtkCellType, Vec<u32>)> = Vec::new();
    let mut cj: Vec<f64> = Vec::new();
    let mut cp: Vec<f64> = Vec::new();
    let mut ci: Vec<f64> = Vec::new();
    let mut cr: Vec<f64> = Vec::new();
    let mut cv: Vec<f64> = Vec::new();
    let mut cn: Vec<f64> = Vec::new();
    let mut cl: Vec<f64> = Vec::new();
    let mut cs: Vec<f64> = Vec::new();
    let mut cch: Vec<f64> = Vec::new();
    let mut cvolt: Vec<f64> = Vec::new();

    // Conducting/diffusion layers: one slab per sheet quad.
    for (qi, q) in sf.quads.iter().enumerate() {
        let layer = &sf.layer[qi];
        let Some((z0, z1)) = layer_z(tech, layer) else {
            continue;
        };
        let base = points.len() as u32;
        for &n in q {
            let (x, y, u) = sf.nodes[n as usize];
            points.push([x, y, z0 * z_scale]);
            potentials.push(u);
        }
        for &n in q {
            let (x, y, u) = sf.nodes[n as usize];
            points.push([x, y, z1 * z_scale]);
            potentials.push(u);
        }
        cells.push((
            VtkCellType::Hex,
            vec![
                base,
                base + 1,
                base + 2,
                base + 3,
                base + 4,
                base + 5,
                base + 6,
                base + 7,
            ],
        ));
        let rs = sf.r_sheet[qi].max(1e-30);
        let g2 = sf.gx[qi] * sf.gx[qi] + sf.gy[qi] * sf.gy[qi];
        let power = g2 / rs * sf.dx[qi] * sf.dy[qi]; // W
        let us: Vec<f64> = q.iter().map(|&n| sf.nodes[n as usize].2).collect();
        let hi = us.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let lo = us.iter().cloned().fold(f64::INFINITY, f64::min);
        let dv = (hi - lo).abs();
        let t = layer_thickness(tech, layer).unwrap_or(1e-6).max(1e-12);
        cj.push(g2.sqrt() / rs / t); // A/m²
        cp.push(power);
        ci.push(if dv > 1e-15 { power / dv } else { 0.0 });
        cr.push(power / (i_total * i_total));
        cv.push(0.0);
        cn.push(sf.component[qi] as f64);
        cl.push(layer_id(layer, &mut all_layers));
        cs.push(1.0 / (rs * t)); // σ = 1/(R_sheet·thickness)
        cvolt.push(0.0);
        cch.push(0.0);
    }

    // Vias: extrude every via polygon between its bottom and top layer z.
    for r in via_regions {
        let Some(layer) = r.layer.as_deref() else {
            continue;
        };
        if tech.via(layer).is_none() {
            continue;
        }
        let (z0, z1) = (r.solid.z_bottom(), r.solid.z_top());
        if z1 <= z0 {
            continue;
        }
        let footprint: Vec<(f64, f64)> = match &r.solid {
            Solid::Box {
                min_x,
                min_y,
                max_x,
                max_y,
                ..
            } => vec![
                (*min_x, *min_y),
                (*max_x, *min_y),
                (*max_x, *max_y),
                (*min_x, *max_y),
            ],
            Solid::Prism { footprint, .. } => footprint.iter().map(|p| (p.x, p.y)).collect(),
        };
        let n = footprint.len();
        if n < 3 {
            continue;
        }
        let area = polygon_area(&footprint).max(1e-30);
        let (cx, cy) = (
            footprint.iter().map(|p| p.0).sum::<f64>() / n as f64,
            footprint.iter().map(|p| p.1).sum::<f64>() / n as f64,
        );
        // Match the via polygon to its solved network link (nearest by center).
        // The link carries the real solved current, `g·ΔV`.
        let mut best: Option<(u32, u32, f64)> = None;
        let mut bd = f64::INFINITY;
        for &(lx, ly, a, b, g) in links {
            let d = (lx - cx) * (lx - cx) + (ly - cy) * (ly - cy);
            if d < bd {
                bd = d;
                best = Some((a, b, g));
            }
        }
        let (ub, ut, g_via, net) = match best {
            Some((a, b, g)) => (
                sf.nodes[a as usize].2,
                sf.nodes[b as usize].2,
                g,
                sf.node_component[a as usize] as f64,
            ),
            None => (0.0, 0.0, 0.0, -1.0),
        };
        let dv = (ut - ub).abs();
        let i_via = g_via * dv;
        let power = dv * i_via;
        let jvia = i_via / area;
        // σ = height / (R·A) with the ICT area-resistance product.
        let ra_prod = tech
            .via(layer)
            .map(|v| {
                if v.area_resistance.len() >= 2 {
                    v.area_resistance[0] * v.area_resistance[1] * 1e-12
                } else {
                    0.0
                }
            })
            .unwrap_or(0.0);
        let sigma_via = if ra_prod > 0.0 { (z1 - z0) / ra_prod } else { 0.0 };
        let base = points.len() as u32;
        for &(x, y) in &footprint {
            points.push([x, y, z0 * z_scale]);
            potentials.push(ub);
        }
        for &(x, y) in &footprint {
            points.push([x, y, z1 * z_scale]);
            potentials.push(ut);
        }
        let lid = layer_id(layer, &mut all_layers);
        let nn = n as u32;
        let nsub = if n == 4 { 1.0 } else { (n - 2) as f64 };
        if n == 4 {
            cells.push((
                VtkCellType::Hex,
                vec![
                    base,
                    base + 1,
                    base + 2,
                    base + 3,
                    base + nn,
                    base + nn + 1,
                    base + nn + 2,
                    base + nn + 3,
                ],
            ));
            cj.push(jvia);
            cp.push(power);
            ci.push(i_via);
            cr.push(power / (i_total * i_total));
            cv.push(i_via);
            cn.push(net);
            cl.push(lid);
            cs.push(sigma_via);
            cvolt.push(0.0);
            cch.push(0.0);
        } else {
            for i in 1..nn - 1 {
                cells.push((
                    VtkCellType::Prism,
                    vec![
                        base,
                        base + i,
                        base + i + 1,
                        base + nn,
                        base + nn + i,
                        base + nn + i + 1,
                    ],
                ));
                cj.push(jvia);
                cp.push(power / nsub);
                ci.push(i_via / nsub);
                cr.push(power / nsub / (i_total * i_total));
                cv.push(i_via / nsub);
                cn.push(net);
                cl.push(lid);
                cs.push(sigma_via);
                cch.push(0.0);
                cvolt.push(0.0);
            }
        }
    }

    // Device channels, placed over the recognized seed footprint and subdivided
    // into the SPI's `nx` x `ny` finger grid. Each device carries its LOCAL
    // (IR-drop aware) bias and Rds; `channel_resistance_ohm` is its Rds.
    let (cz0, cz1) = layer_z(tech, "RX").unwrap_or((0.0, 0.1e-6));
    for ch in &sf.channels {
        let (xd, yd, ud) = sf.nodes[ch.d as usize];
        let (xs, ys, us) = sf.nodes[ch.s as usize];
        let (ddx, ddy) = (xs - xd, ys - yd);
        let len_ds = (ddx * ddx + ddy * ddy).sqrt();
        let l2 = (ddx * ddx + ddy * ddy).max(1e-30);
        let quads: Vec<[(f64, f64); 4]> = if let Some(bb) = ch.bbox {
            let nx = ch.nx.max(1);
            let ny = ch.ny.max(1);
            let cw = (bb.max_x - bb.min_x) / nx as f64;
            let chh = (bb.max_y - bb.min_y) / ny as f64;
            let mut v = Vec::with_capacity(nx * ny);
            for i in 0..nx {
                for j in 0..ny {
                    let (x0, y0) = (bb.min_x + i as f64 * cw, bb.min_y + j as f64 * chh);
                    let (x1, y1) = (x0 + cw, y0 + chh);
                    v.push([(x0, y0), (x1, y0), (x1, y1), (x0, y1)]);
                }
            }
            v
        } else {
            if len_ds < 1e-12 {
                continue;
            }
            let w = (0.5 * len_ds).max(0.2e-6);
            let (px, py) = (-ddy / len_ds * w / 2.0, ddx / len_ds * w / 2.0);
            vec![[
                (xd + px, yd + py),
                (xs + px, ys + py),
                (xs - px, ys - py),
                (xd - px, yd - py),
            ]]
        };
        let nc = quads.len() as f64;
        let lid = layer_id("Channel", &mut all_layers);
        let net = sf.node_component[ch.d as usize] as f64;
        // Device-level channel conductivity σ = L / (R·A).
        let sigma = if ch.resistance.is_finite() && ch.resistance > 0.0 && ch.weff > 0.0 {
            len_ds / (ch.resistance * ch.weff * (cz1 - cz0).max(1e-30))
        } else {
            0.0
        };
        // Channel-length interpolation axis: the bbox's shorter side (perpendicular
        // to the gate). Interpolating along the raw drain->source node vector
        // instead gives a diagonal/vertical gradient, because the nearest drain and
        // source mesh nodes are not exactly aligned -- a discretization artifact,
        // not the physical channel direction.
        let bbox_axis = ch.bbox.map(|b| {
            let length_is_x = b.width() < b.height();
            let (lo, hi) = if length_is_x { (b.min_x, b.max_x) } else { (b.min_y, b.max_y) };
            let ds_along = if length_is_x { ddx } else { ddy };
            (length_is_x, lo, hi, ds_along >= 0.0)
        });
        for quad in &quads {
            let (ccx, ccy) = (
                quad.iter().map(|p| p.0).sum::<f64>() / 4.0,
                quad.iter().map(|p| p.1).sum::<f64>() / 4.0,
            );
            let t = match bbox_axis {
                Some((length_is_x, lo, hi, drain_at_lo)) => {
                    let along = if length_is_x {
                        ((ccx - lo) / (hi - lo).max(1e-30)).clamp(0.0, 1.0)
                    } else {
                        ((ccy - lo) / (hi - lo).max(1e-30)).clamp(0.0, 1.0)
                    };
                    if drain_at_lo { along } else { 1.0 - along }
                }
                None => (((ccx - xd) * ddx + (ccy - yd) * ddy) / l2).clamp(0.0, 1.0),
            };
            let ucell = ud + t * (us - ud);
            let base = points.len() as u32;
            for &(x, y) in quad {
                points.push([x, y, cz0 * z_scale]);
                potentials.push(ucell);
            }
            for &(x, y) in quad {
                points.push([x, y, cz1 * z_scale]);
                potentials.push(ucell);
            }
            cells.push((
                VtkCellType::Hex,
                vec![
                    base,
                    base + 1,
                    base + 2,
                    base + 3,
                    base + 4,
                    base + 5,
                    base + 6,
                    base + 7,
                ],
            ));
            let area = polygon_area(quad).max(1e-30) * (cz1 - cz0).max(1e-30);
            cj.push(ch.current / nc / area);
            cp.push(ch.power / nc);
            ci.push(ch.current / nc);
            cr.push(ch.power / nc / (i_total * i_total));
            cv.push(0.0);
            cn.push(net);
            cl.push(lid);
            cs.push(sigma);
            cch.push(ch.resistance);
            cvolt.push(0.0);
        }
    }

    // Terminal contacts: a box over each terminal's mesh nodes at its contact
    // layer, coloured by its voltage.
    for t in &sf.terminals {
        let (mut x0, mut y0, mut x1, mut y1) = (
            f64::INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NEG_INFINITY,
        );
        for &n in &t.nodes {
            let (x, y, _) = sf.nodes[n as usize];
            x0 = x0.min(x);
            y0 = y0.min(y);
            x1 = x1.max(x);
            y1 = y1.max(y);
        }
        if !x0.is_finite() {
            continue;
        }
        let pad = 0.25e-6;
        if x1 - x0 < pad {
            x0 -= pad;
            x1 += pad;
        }
        if y1 - y0 < pad {
            y0 -= pad;
            y1 += pad;
        }
        let (_, z1) = layer_z(tech, &t.layer).unwrap_or((0.0, 0.2e-6));
        let corners = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)];
        let base = points.len() as u32;
        for &(x, y) in &corners {
            points.push([x, y, z1 * z_scale]);
            potentials.push(t.voltage);
        }
        for &(x, y) in &corners {
            points.push([x, y, (z1 + 0.3e-6) * z_scale]);
            potentials.push(t.voltage);
        }
        cells.push((
            VtkCellType::Hex,
            vec![
                base,
                base + 1,
                base + 2,
                base + 3,
                base + 4,
                base + 5,
                base + 6,
                base + 7,
            ],
        ));
        cj.push(0.0);
        cp.push(0.0);
        ci.push(0.0);
        cr.push(0.0);
        cv.push(0.0);
        cn.push(sf.node_component[t.nodes[0] as usize] as f64);
        cl.push(layer_id("Terminal", &mut all_layers));
        cs.push(0.0);
        cch.push(0.0);
        cvolt.push(t.voltage);
    }

    // Per-cell quadtree refinement level: the sheet quads are the first
    // `sf.quads.len()` cells (vias/channels/terminals follow), so the level
    // array is the network's `cell_level` padded with zeros for the rest.
    let mut crl: Vec<f64> = cell_level.iter().map(|&l| l as f64).collect();
    crl.resize(cells.len(), 0.0);

    let vm = VtuMesh {
        points,
        cells,
        point_data: vec![("potential_V".to_string(), potentials)],
        cell_data: vec![
            ("current_density_A_m2".to_string(), cj),
            ("power_W".to_string(), cp),
            ("current_A".to_string(), ci),
            ("resistance_ohm".to_string(), cr),
            ("via_current_A".to_string(), cv),
            ("net_id".to_string(), cn),
            ("layer_id".to_string(), cl),
            ("conductivity_S_m".to_string(), cs),
            ("channel_resistance_ohm".to_string(), cch),
            ("contact_voltage_V".to_string(), cvolt),
            ("refinement_level".to_string(), crl),
        ],
    };
    (vm, all_layers)
}

fn run_viz(
    settings: &openrdson_validation::MeshSettings,
    config: Option<&openrdson_config::SolverConfig>,
    prebuilt: Option<openrdson_validation::SheetFullArray>,
) {
    let _section = openrdson_core::log::SectionGuard::new("viz");
    if prebuilt.is_none() {
        log_databases(settings);
    }
    use openrdson_core::{log_info, log_warn};
    use openrdson_io::{
        write_gds_file, write_grouped_colormap,
        write_named_lyp, write_vtu, GdsBoundary, VizItem, VtkCellType, VtuMesh,
    };
    let cwd = std::env::current_dir().unwrap_or_else(|_| std::path::PathBuf::from("."));
    let out = config
        .map(|c| cwd.join(&c.project.output_dir))
        .unwrap_or_else(|| cwd.join("openrdson_viz"));
    if let Err(e) = std::fs::create_dir_all(&out) {
        eprintln!("viz: cannot create {}: {e}", out.display());
        return;
    }

    let temperature = config.map(|c| c.solver.temperature_c).unwrap_or(27.0);
    let fp = openrdson_validation::field_cache_fingerprint(settings);
    let cache_path = out.join(format!("field_cache_{fp}.bin"));
    // When we already hold a freshly-solved array (the `all` path), use it
    // directly rather than loading a cached field: a stale cache would export a
    // field that disagrees with the Rds table `run_sheet_rds` just printed.
    let standalone = prebuilt.is_none();
    let cache = match prebuilt {
        Some(fa) => {
            // `prebuilt` is already the final array `run_sheet_rds` solved and
            // reported (adaptive or uniform); just render its field. Re-running
            // `adaptive_solve_full` here would re-do the whole refinement loop.
            let c = match fa.field_cache(temperature) {
                Ok(c) => c,
                Err(e) => {
                    eprintln!("viz solve failed: {e}");
                    return;
                }
            };
            let _ = c.save(&cache_path);
            c
        }
        None => match openrdson_validation::FieldCache::load(&cache_path) {
            Ok(c) => {
                log_info!("loaded field cache {}", cache_path.display());
                c
            }
            Err(_) => {
                log_info!("building full-array sheet extraction for visualization...");
                let mut fa = match openrdson_validation::SheetFullArray::build_with(settings) {
                    Ok(f) => f,
                    Err(e) => {
                        eprintln!("viz failed: {e}");
                        return;
                    }
                };
                // When adaptive refinement is enabled, solve on the quadtree mesh.
                if settings.adaptive.enabled {
                    match fa.adaptive_solve_full(&settings.adaptive, temperature) {
                        Ok((report, fa_adaptive)) => {
                            log_info!(
                                "adaptive mesh: Rds={:.6e} ohm nodes={} cells={} iters={} levels={:?}",
                                report.rds,
                                report.nodes,
                                report.cells,
                                report.iters,
                                report.levels
                            );
                            fa = fa_adaptive;
                        }
                        Err(e) => log_warn!("adaptive solve failed (using uniform mesh): {e}"),
                    }
                }
                let c = match fa.field_cache(temperature) {
                    Ok(c) => c,
                    Err(e) => {
                        eprintln!("viz solve failed: {e}");
                        return;
                    }
                };
                match c.save(&cache_path) {
                    Ok(()) => log_info!("saved field cache {}", cache_path.display()),
                    Err(e) => log_warn!("failed to save field cache: {e}"),
                }
                c
            }
        },
    };
    let layout_db = cache.db_unit_meters;
    let args: Vec<String> = std::env::args().collect();
    let get_arg = |flag: &str| -> Option<String> {
        args.windows(2).find(|w| w[0] == flag).map(|w| w[1].clone())
    };
    // Output database unit (nm): CLI > config > layout unit.
    let out_dbu_nm = get_arg("--out-dbu-nm")
        .and_then(|v| v.parse::<f64>().ok())
        .or_else(|| config.and_then(|c| c.visualization.out_db_unit_nm));
    let db = out_dbu_nm.map(|nm| nm * 1e-9).unwrap_or(layout_db);
    log_info!("output database unit = {:.4} nm", db * 1e9);
    // Vertical exaggeration for the 3D stack: the die is much wider than it is
    // tall, so layers are invisible at true scale. CLI > config > 25.
    let z_scale = get_arg("--viz-zscale")
        .and_then(|v| v.parse::<f64>().ok())
        .or_else(|| config.map(|c| c.visualization.z_scale))
        .unwrap_or(25.0);
    let n_bins = config.map(|c| c.visualization.bins).unwrap_or(32).max(2);
    // A standalone `viz` prints the (cached or freshly-solved) resistance report;
    // the `all` path already printed it via `run_sheet_rds`, so it is skipped here.
    if standalone && !cache.report.is_empty() {
        print!("{}", cache.report);
    }
    let cells = cache.cell_scalars();
    log_info!(
        "redrawing the layout from {} mesh cells across {} layers",
        cells.len(),
        {
            let mut s = std::collections::BTreeSet::new();
            for (p, _, _) in &cells {
                s.insert(p.clone());
            }
            s.len()
        }
    );

    // One GDS layer per physical layer for the mesh mosaic.
    let mut phys_layer: std::collections::BTreeMap<String, i16> = std::collections::BTreeMap::new();
    for (phys, _, _) in &cells {
        let n = phys_layer.len() as i16;
        phys_layer.entry(phys.clone()).or_insert(3000 + n);
    }
    let mut mesh_bounds = Vec::with_capacity(cells.len());
    for (phys, c, _s) in &cells {
        let pts: Vec<(i32, i32)> = c
            .iter()
            .map(|&(x, y)| (to_du(x, db), to_du(y, db)))
            .collect();
        mesh_bounds.push(GdsBoundary {
            layer: *phys_layer.get(phys).unwrap(),
            datatype: 0,
            points: pts,
        });
    }

    let mesh_path = out.join("mesh.gds");
    if write_gds_file(&mesh_path, "OPENRDSON", db, &[("MESH".into(), mesh_bounds)]).is_ok() {
        let entries: Vec<(i16, String)> =
            phys_layer.iter().map(|(k, v)| (*v, k.clone())).collect();
        let mut lyp = mesh_path.clone();
        lyp.set_extension("lyp");
        let _ = write_named_lyp(&lyp, &entries);
        log_info!(
            "wrote {} ({} cells, {} layers)",
            mesh_path.display(),
            cells.len(),
            phys_layer.len()
        );
    }

    // One grouped GDS colormap per scalar quantity: one KLayout group per
    // physical layer, each group containing the colormap bins.
    let to_db = |(x, y): (f64, f64)| (to_du(x, db), to_du(y, db));
    let write_scalar = |name: &str,
                        base: i16,
                        pick: &dyn Fn(&openrdson_validation::CellScalars) -> f64| {
        let mut by_layer: std::collections::BTreeMap<String, Vec<VizItem>> =
            std::collections::BTreeMap::new();
        for (phys, c, s) in cells.iter() {
            let pts: Vec<(i32, i32)> = c.iter().map(|&(x, y)| to_db((x, y))).collect();
            by_layer
                .entry(phys.clone())
                .or_default()
                .push(VizItem { points: pts, value: pick(s) });
        }
        let layers: Vec<(String, Vec<VizItem>)> = by_layer.into_iter().collect();
        let path = out.join(format!("{name}.gds"));
        match write_grouped_colormap(
            &path,
            "OPENRDSON",
            &name.to_uppercase(),
            db,
            &layers,
            n_bins,
            base,
            32,
        ) {
            Ok(()) => log_info!(
                "wrote {} ({} groups x {n_bins} bins)",
                path.display(),
                layers.len()
            ),
            Err(_) => log_warn!("failed to write {name} map"),
        }
    };
    write_scalar("potential", 1000, &|s| s.potential);
    write_scalar("current_density", 2000, &|s| s.current_density);
    write_scalar("power", 4000, &|s| s.power);
    write_scalar("current", 5000, &|s| s.current);
    write_scalar("resistance", 6000, &|s| s.resistance);

    // VTK/ParaView export of the sheet network: per-vertex potential plus
    // per-cell current density, power, current, resistance contribution and
    // net id — all selectable quantities in the viewer.
    {
        let sf: &openrdson_validation::SheetField = &cache.field;
        let links = cache.via_links();
            // Device channels: one KLayout group coloured by channel resistance.
            {
                let mut ch_items: Vec<VizItem> = Vec::new();
                for ch in &sf.channels {
                    let Some(bb) = ch.bbox else { continue };
                    let (nx, ny) = (ch.nx.max(1), ch.ny.max(1));
                    let cw = (bb.max_x - bb.min_x) / nx as f64;
                    let chh = (bb.max_y - bb.min_y) / ny as f64;
                    for i in 0..nx {
                        for j in 0..ny {
                            let (x0, y0) =
                                (bb.min_x + i as f64 * cw, bb.min_y + j as f64 * chh);
                            let pts = vec![
                                (to_du(x0, db), to_du(y0, db)),
                                (to_du(x0 + cw, db), to_du(y0, db)),
                                (to_du(x0 + cw, db), to_du(y0 + chh, db)),
                                (to_du(x0, db), to_du(y0 + chh, db)),
                            ];
                            ch_items.push(VizItem { points: pts, value: ch.resistance });
                        }
                    }
                }
                if !ch_items.is_empty() {
                    let layers: Vec<(String, Vec<VizItem>)> =
                        vec![("Channel".to_string(), ch_items)];
                    let path = out.join("channel_resistance.gds");
                    match write_grouped_colormap(
                        &path, "OPENRDSON", "CHANNELR", db, &layers, n_bins, 7000, 32,
                    ) {
                        Ok(()) => log_info!(
                            "wrote {} ({} groups x {n_bins} bins)",
                            path.display(),
                            layers.len()
                        ),
                        Err(_) => log_warn!("failed to write channel resistance map"),
                    }
                }
            }

            // Vias: one KLayout group coloured by via current (conductance * dV).
            {
                let mut via_items: Vec<VizItem> = Vec::new();
                let s = 0.25e-6;
                for &(x, y, b, t, g) in &links {
                    let i_via = g * (sf.nodes[t as usize].2 - sf.nodes[b as usize].2);
                    let pts = vec![
                        (to_du(x - s, db), to_du(y - s, db)),
                        (to_du(x + s, db), to_du(y - s, db)),
                        (to_du(x + s, db), to_du(y + s, db)),
                        (to_du(x - s, db), to_du(y + s, db)),
                    ];
                    via_items.push(VizItem { points: pts, value: i_via });
                }
                if !via_items.is_empty() {
                    let layers: Vec<(String, Vec<VizItem>)> =
                        vec![("Vias".to_string(), via_items)];
                    let path = out.join("via_current.gds");
                    match write_grouped_colormap(
                        &path, "OPENRDSON", "VIACURRENT", db, &layers, n_bins, 8000, 32,
                    ) {
                        Ok(()) => log_info!(
                            "wrote {} ({} groups x {n_bins} bins)",
                            path.display(),
                            layers.len()
                        ),
                        Err(_) => log_warn!("failed to write via current map"),
                    }
                }
            }
            let mut layer_ids: std::collections::BTreeMap<String, usize> =
                std::collections::BTreeMap::new();
            let mut layer_names: Vec<String> = Vec::new();
            let cell_layer: Vec<f64> = sf
                .layer
                .iter()
                .map(|l| {
                    let next = layer_names.len();
                    *layer_ids.entry(l.clone()).or_insert_with(|| {
                        layer_names.push(l.clone());
                        next
                    }) as f64
                })
                .collect();
            let n = sf.quads.len();
            let it2 = (sf.total_current * sf.total_current).max(1e-300);
            let mut cell_j = Vec::with_capacity(n);
            let mut cell_power = Vec::with_capacity(n);
            let mut cell_current = Vec::with_capacity(n);
            let mut cell_res = Vec::with_capacity(n);
            let mut cell_net = Vec::with_capacity(n);
            for i in 0..n {
                let rs = sf.r_sheet[i].max(1e-30);
                let g2 = sf.gx[i] * sf.gx[i] + sf.gy[i] * sf.gy[i];
                let power = g2 / rs * sf.dx[i] * sf.dy[i];
                let us: Vec<f64> =
                    sf.quads[i].iter().map(|&nd| sf.nodes[nd as usize].2).collect();
                let hi = us.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                let lo = us.iter().cloned().fold(f64::INFINITY, f64::min);
                let dv = (hi - lo).abs();
                cell_j.push(g2.sqrt() / rs);
                cell_power.push(power);
                cell_current.push(if dv > 1e-15 { power / dv } else { 0.0 });
                cell_res.push(power / it2);
                cell_net.push(sf.component[i] as f64);
            }
            let vm = VtuMesh {
                points: sf.nodes.iter().map(|&(x, y, _)| [x, y, 0.0]).collect(),
                cells: sf
                    .quads
                    .iter()
                    .map(|q| (VtkCellType::Quad, q.to_vec()))
                    .collect(),
                point_data: vec![(
                    "potential_V".to_string(),
                    sf.nodes.iter().map(|&(_, _, u)| u).collect(),
                )],
                cell_data: vec![
                    ("current_density_A_m".to_string(), cell_j),
                    ("power_W".to_string(), cell_power),
                    ("current_A".to_string(), cell_current),
                    ("resistance_ohm".to_string(), cell_res),
                    ("net_id".to_string(), cell_net),
                    ("layer_id".to_string(), cell_layer.clone()),
                    (
                        "refinement_level".to_string(),
                        cache.cell_level.iter().map(|&l| l as f64).collect(),
                    ),
                ],
            };
            let path = out.join("sheet.vtu");
            match write_vtu(&path, &vm) {
                Ok(()) => log_info!(
                    "wrote {} ({} nodes, {} quads, {} layers)",
                    path.display(),
                    vm.points.len(),
                    vm.cells.len(),
                    layer_names.len()
                ),
                Err(e) => log_warn!("failed to write sheet.vtu: {e}"),
            }
            // One file per physical layer, e.g. `sheet_Metal1.vtu`.
            for (id, name) in layer_names.iter().enumerate() {
                let keep: Vec<bool> = cell_layer.iter().map(|&l| l as usize == id).collect();
                if !keep.iter().any(|&k| k) {
                    continue;
                }
                let sub = vm.subset_cells(&keep);
                let f = out.join(format!("sheet_{}.vtu", sanitize(name)));
                match write_vtu(&f, &sub) {
                    Ok(()) => log_info!(
                        "wrote {} ({} nodes, {} quads)",
                        f.display(),
                        sub.points.len(),
                        sub.cells.len()
                    ),
                    Err(e) => log_warn!("failed to write {}: {e}", f.display()),
                }
            }

            // 2.5D field stack: each layer's sheet mesh extruded to its ICT
            // z-range (correct per-polygon vias), carrying every quantity.
            if let Some((tech, via_regions)) = load_tech_and_regions(settings) {
                let (vm, all_layers) =
                    build_2p5d_stack(&tech, &via_regions, &links, &cache.cell_level, sf, z_scale);
                let path = out.join("stack.vtu");
                match write_vtu(&path, &vm) {
                    Ok(()) => log_info!(
                        "wrote {} ({} nodes, {} cells, {} layers; z x{})",
                        path.display(),
                        vm.points.len(),
                        vm.cells.len(),
                        all_layers.len(),
                        z_scale
                    ),
                    Err(e) => log_warn!("failed to write stack.vtu: {e}"),
                }
                for (id, name) in all_layers.iter().enumerate() {
                    let keep: Vec<bool> =
                        vm.cell_data[6].1.iter().map(|&l| l as usize == id).collect();
                    if !keep.iter().any(|&k| k) {
                        continue;
                    }
                    let sub = vm.subset_cells(&keep);
                    let f = out.join(format!("stack_{}.vtu", sanitize(name)));
                    match write_vtu(&f, &sub) {
                        Ok(()) => log_info!(
                            "wrote {} ({} nodes, {} cells)",
                            f.display(),
                            sub.points.len(),
                            sub.cells.len()
                        ),
                        Err(e) => log_warn!("failed to write {}: {e}", f.display()),
                    }
                }
            }
    }

    println!("visualization written to {}", out.display());
    println!("mesh.gds (per-layer mesh mosaic), potential.gds, current_density.gds");
    println!("open a .gds in KLayout; its .lyp is auto-loaded (same base name)");
    println!("sheet.vtu (VTK/ParaView: interpolated colormap, quantity dropdown)");
    println!("sheet_<Layer>.vtu (one file per physical layer)");
    println!("stack.vtu (2.5D field stack: all layers + every via, potential + current density)");
    println!("stack_<Layer>.vtu (per-layer 3D stacks, incl. Terminal contacts)");
}

fn run_all(
    settings: &openrdson_validation::MeshSettings,
    config: Option<&openrdson_config::SolverConfig>,
) {
    let _section = openrdson_core::log::SectionGuard::new("all");

    // Full-array Rdson extraction: solve and print the Rds(Vgs) table, then
    // reuse the solved field for the visualization export (no re-solve).
    println!("\n=== Rdson (full-array sheet extraction) ===");
    let summary_fa = run_sheet_rds(settings);

    // Export the solved field for viewing (KLayout + ParaView).
    run_viz(settings, config, summary_fa);
}


//! OpenRDSon solver configuration: a documented YAML schema plus a
//! dependency-free loader.
//!
//! Point the executable at a config file and the whole pipeline is driven from
//! it:
//!
//! ```text
//! openrdson --config openrdson.yaml viz
//! ```
//!
//! See [`SolverConfig::default_yaml`] for the fully commented template (also
//! emitted by `openrdson --print-default-config`).

pub mod yaml;

use std::collections::BTreeMap;
use std::path::PathBuf;
use yaml::Value;

/// Input files and output location.
#[derive(Debug, Clone)]
pub struct Project {
    /// Directory holding the CCI database (layout, map, ports, devtab, spi).
    pub cci_dir: PathBuf,
    pub layout: String,
    pub gds_map: String,
    pub ports: String,
    pub devtab: String,
    pub spi: String,
    /// ICT process stack.
    pub tech_ict: PathBuf,
    /// Logical -> physical layer map (`*.map`).
    pub layer_map: PathBuf,
    /// Bias-dependent channel model table (`model.csv`).
    pub model_csv: PathBuf,
    /// Directory the visualization database is written to.
    pub output_dir: PathBuf,
}

impl Default for Project {
    fn default() -> Self {
        Self {
            cci_dir: PathBuf::new(),
            layout: String::new(),
            gds_map: String::new(),
            ports: String::new(),
            devtab: String::new(),
            spi: String::new(),
            tech_ict: PathBuf::new(),
            layer_map: PathBuf::new(),
            model_csv: PathBuf::new(),
            output_dir: PathBuf::from("openrdson_viz"),
        }
    }
}

/// Unit handling. `db_unit_nm = None` means "infer from the golden SPI vs the
/// recognized seed footprint".
#[derive(Debug, Clone, Default)]
pub struct Units {
    pub db_unit_nm: Option<f64>,
}

/// Meshing / discretization settings.
#[derive(Debug, Clone)]
pub struct Meshing {
    pub sheet_cell_um: f64,
    pub sheet_cap: usize,
    pub via_group_radius_um: f64,
    /// Per-physical-layer sheet cell overrides (µm).
    pub per_layer_cell_um: BTreeMap<String, f64>,
    /// Adaptive (quadtree) sheet refinement.
    pub adaptive: bool,
    /// Relative indicator convergence tolerance.
    pub adaptive_tol: f64,
    /// Maximum quadtree refinement level.
    pub adaptive_max_level: u8,
    /// Node budget; refinement stops above it.
    pub adaptive_max_nodes: usize,
    /// Dörfler marking fraction (0..1).
    pub adaptive_mark_frac: f64,
    /// Maximum refine->solve iterations.
    pub adaptive_iters: usize,
    /// A-priori refinement radius around vias/terminals (µm); 0 disables.
    pub refine_near_vias_um: f64,
}

impl Default for Meshing {
    fn default() -> Self {
        Self {
            sheet_cell_um: 2.0,
            sheet_cap: 30,
            via_group_radius_um: 5.0,
            per_layer_cell_um: BTreeMap::new(),
            adaptive: false,
            adaptive_tol: 1e-3,
            adaptive_max_level: 3,
            adaptive_max_nodes: 400_000,
            adaptive_mark_frac: 0.5,
            adaptive_iters: 8,
            refine_near_vias_um: 0.0,
        }
    }
}

/// A user-defined voltage terminal (a contact on a net).
///
/// `name` should match a ports-file entry where one exists; extra "special
/// node" terminals are allowed. The contact ties every mesh node inside the
/// `dx_um × dy_um` rectangle, on the terminal's net, to `voltage`.
#[derive(Debug, Clone)]
pub struct Terminal {
    /// Port name (may repeat for multiple contacts on one net) or a free-form
    /// special-node name.
    pub name: String,
    /// Absolute terminal voltage (V).
    pub voltage: f64,
    /// Contact centre, in µm relative to the layout origin.
    pub x_um: f64,
    pub y_um: f64,
    /// Contact footprint (µm). Defaults to a single nearest node when 0.
    pub dx_um: f64,
    pub dy_um: f64,
    /// Contact level (physical ICT name, else logical gds.map name). Defaults
    /// to the port's layer.
    pub layer: Option<String>,
}

impl Default for Terminal {
    fn default() -> Self {
        Self {
            name: String::new(),
            voltage: 0.0,
            x_um: 0.0,
            y_um: 0.0,
            dx_um: 0.0,
            dy_um: 0.0,
            layer: None,
        }
    }
}

/// What to do when an active device has no channel model available.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MissingModel {
    /// Abort with an error (default).
    Error,
    /// Bind the device to a very large ("open") channel resistance.
    Open,
}

/// Bias point and solve controls.
#[derive(Debug, Clone)]
pub struct Solver {
    pub temperature_c: f64,
    /// Re-solve each device's channel at its own IR-dropped Vds.
    pub local_bias: bool,
    /// Relative convergence tolerance for the local-bias iteration.
    pub tolerance: f64,
    pub max_iterations: usize,
    /// Behavior when a device model is required but absent.
    pub missing_model: MissingModel,
}

impl Default for Solver {
    fn default() -> Self {
        Self {
            temperature_c: 27.0,
            local_bias: true,
            tolerance: 1.0e-3,
            max_iterations: 20,
            missing_model: MissingModel::Error,
        }
    }
}

/// Resource controls.
#[derive(Debug, Clone, Default)]
pub struct Performance {
    /// Worker threads; 0 = all available cores.
    pub threads: usize,
    /// Soft memory cap in MB; 0 = unlimited.
    pub max_memory_mb: u64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LogLevel {
    Error,
    Warn,
    Info,
    Debug,
    Trace,
}

impl LogLevel {
    fn parse(s: &str) -> Option<Self> {
        match s.trim().to_ascii_lowercase().as_str() {
            "error" => Some(Self::Error),
            "warn" | "warning" => Some(Self::Warn),
            "info" => Some(Self::Info),
            "debug" => Some(Self::Debug),
            "trace" => Some(Self::Trace),
            _ => None,
        }
    }
}

#[derive(Debug, Clone)]
pub struct Logging {
    pub level: LogLevel,
}

impl Default for Logging {
    fn default() -> Self {
        Self {
            level: LogLevel::Info,
        }
    }
}

/// Visualization output options.
#[derive(Debug, Clone)]
pub struct Visualization {
    pub enabled: bool,
    /// Vertical exaggeration for the 2.5D stack.
    pub z_scale: f64,
    /// Colormap bins for the KLayout GDS maps.
    pub bins: usize,
    /// Output database unit (nm); `None` = the layout's own unit.
    pub out_db_unit_nm: Option<f64>,
}

impl Default for Visualization {
    fn default() -> Self {
        Self {
            enabled: true,
            z_scale: 25.0,
            bins: 32,
            out_db_unit_nm: None,
        }
    }
}

/// The complete solver configuration.
#[derive(Debug, Clone, Default)]
pub struct SolverConfig {
    pub project: Project,
    pub units: Units,
    pub meshing: Meshing,
    pub solver: Solver,
    pub terminals: Vec<Terminal>,
    pub performance: Performance,
    pub logging: Logging,
    pub visualization: Visualization,
}

impl SolverConfig {
    /// The documented YAML template.
    pub fn default_yaml() -> &'static str {
        include_str!("default.yaml")
    }

    /// Parse a configuration from YAML text. Missing keys keep their defaults.
    pub fn from_yaml(text: &str) -> Result<Self, String> {
        let root = yaml::parse(text)?;
        let mut c = SolverConfig::default();

        if let Some(p) = root.get("project") {
            c.project.cci_dir = path_of(p.get("cci_dir"), &c.project.cci_dir);
            c.project.layout = str_of(p.get("layout"), &c.project.layout);
            c.project.gds_map = str_of(p.get("gds_map"), &c.project.gds_map);
            c.project.ports = str_of(p.get("ports"), &c.project.ports);
            c.project.devtab = str_of(p.get("devtab"), &c.project.devtab);
            c.project.spi = str_of(p.get("spi"), &c.project.spi);
            c.project.tech_ict = path_of(p.get("tech_ict"), &c.project.tech_ict);
            c.project.layer_map = path_of(p.get("layer_map"), &c.project.layer_map);
            c.project.model_csv = path_of(p.get("model_csv"), &c.project.model_csv);
            c.project.output_dir = path_of(p.get("output_dir"), &c.project.output_dir);
        }

        if let Some(u) = root.get("units") {
            if let Some(v) = u.get("db_unit_nm") {
                c.units.db_unit_nm = match v {
                    Value::Str(s) if s.eq_ignore_ascii_case("auto") => None,
                    _ => v.as_f64().filter(|x| *x > 0.0),
                };
            }
        }

        if let Some(m) = root.get("meshing") {
            c.meshing.sheet_cell_um = f64_of(m.get("sheet_cell_um"), c.meshing.sheet_cell_um);
            c.meshing.sheet_cap = usize_of(m.get("sheet_cap"), c.meshing.sheet_cap);
            c.meshing.via_group_radius_um =
                f64_of(m.get("via_group_radius_um"), c.meshing.via_group_radius_um);
            c.meshing.adaptive = bool_of(m.get("adaptive"), c.meshing.adaptive);
            c.meshing.adaptive_tol = f64_of(m.get("adaptive_tol"), c.meshing.adaptive_tol);
            c.meshing.adaptive_max_level = usize_of(
                m.get("adaptive_max_level"),
                c.meshing.adaptive_max_level as usize,
            ) as u8;
            c.meshing.adaptive_max_nodes =
                usize_of(m.get("adaptive_max_nodes"), c.meshing.adaptive_max_nodes);
            c.meshing.adaptive_mark_frac =
                f64_of(m.get("adaptive_mark_frac"), c.meshing.adaptive_mark_frac);
            c.meshing.adaptive_iters =
                usize_of(m.get("adaptive_iters"), c.meshing.adaptive_iters);
            c.meshing.refine_near_vias_um =
                f64_of(m.get("refine_near_vias_um"), c.meshing.refine_near_vias_um);
            if let Some(map) = m.get("per_layer_cell_um").and_then(|v| match v {
                Value::Map(mm) => Some(mm),
                _ => None,
            }) {
                for (k, v) in map {
                    if let Some(um) = v.as_f64() {
                        c.meshing.per_layer_cell_um.insert(k.clone(), um);
                    }
                }
            }
        }

        if let Some(s) = root.get("solver") {
            c.solver.temperature_c = f64_of(s.get("temperature_c"), c.solver.temperature_c);
            c.solver.local_bias = bool_of(s.get("local_bias"), c.solver.local_bias);
            c.solver.tolerance = f64_of(s.get("tolerance"), c.solver.tolerance);
            c.solver.max_iterations = usize_of(s.get("max_iterations"), c.solver.max_iterations);
            if let Some(m) = s.get("missing_model").and_then(|v| v.as_str()) {
                c.solver.missing_model = match m.to_ascii_lowercase().as_str() {
                    "open" => MissingModel::Open,
                    _ => MissingModel::Error,
                };
            }
        }

        if let Some(list) = root.get("terminals").and_then(|v| v.as_list()) {
            for item in list {
                let t = Terminal {
                    name: str_of(item.get("name"), ""),
                    voltage: f64_of(item.get("voltage"), 0.0),
                    x_um: f64_of(item.get("x_um"), 0.0),
                    y_um: f64_of(item.get("y_um"), 0.0),
                    dx_um: f64_of(item.get("dx_um"), 0.0),
                    dy_um: f64_of(item.get("dy_um"), 0.0),
                    layer: item
                        .get("layer")
                        .and_then(|v| v.as_str())
                        .map(|s| s.to_string()),
                };
                if !t.name.is_empty() {
                    c.terminals.push(t);
                }
            }
        }

        if let Some(p) = root.get("performance") {
            c.performance.threads = usize_of(p.get("threads"), c.performance.threads);
            c.performance.max_memory_mb = usize_of(p.get("max_memory_mb"), 0) as u64;
        }

        if let Some(l) = root.get("logging") {
            if let Some(level) = l.get("level").and_then(|v| v.as_str()).and_then(LogLevel::parse) {
                c.logging.level = level;
            }
        }

        if let Some(v) = root.get("visualization") {
            c.visualization.enabled = bool_of(v.get("enabled"), c.visualization.enabled);
            c.visualization.z_scale = f64_of(v.get("z_scale"), c.visualization.z_scale);
            c.visualization.bins = usize_of(v.get("bins"), c.visualization.bins);
            if let Some(o) = v.get("out_db_unit_nm") {
                c.visualization.out_db_unit_nm = match o {
                    Value::Str(s) if s.eq_ignore_ascii_case("auto") => None,
                    _ => o.as_f64().filter(|x| *x > 0.0),
                };
            }
        }

        Ok(c)
    }
}

fn str_of(v: Option<&Value>, default: &str) -> String {
    v.and_then(|v| v.as_str())
        .map(|s| s.to_string())
        .unwrap_or_else(|| default.to_string())
}

fn path_of(v: Option<&Value>, default: &PathBuf) -> PathBuf {
    v.and_then(|v| v.as_str())
        .map(PathBuf::from)
        .unwrap_or_else(|| default.clone())
}

fn f64_of(v: Option<&Value>, default: f64) -> f64 {
    v.and_then(|v| v.as_f64()).unwrap_or(default)
}

fn usize_of(v: Option<&Value>, default: usize) -> usize {
    v.and_then(|v| v.as_i64())
        .filter(|i| *i >= 0)
        .map(|i| i as usize)
        .unwrap_or(default)
}

fn bool_of(v: Option<&Value>, default: bool) -> bool {
    v.and_then(|v| v.as_bool()).unwrap_or(default)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_nested_config() {
        let cfg = SolverConfig::from_yaml(
            "project:\n  cci_dir: /tmp/cci\n  layout: foo.agf\n\
             meshing:\n  sheet_cell_um: 0.5\n  per_layer_cell_um:\n    Metal1: 0.25\n\
             solver:\n  temperature_c: 85.0\n  local_bias: false\n\
             logging:\n  level: debug\n",
        )
        .unwrap();
        assert_eq!(cfg.project.cci_dir, PathBuf::from("/tmp/cci"));
        assert_eq!(cfg.project.layout, "foo.agf");
        assert_eq!(cfg.meshing.sheet_cell_um, 0.5);
        assert_eq!(cfg.meshing.per_layer_cell_um.get("Metal1"), Some(&0.25));
        assert_eq!(cfg.solver.temperature_c, 85.0);
        assert!(!cfg.solver.local_bias);
        assert_eq!(cfg.logging.level, LogLevel::Debug);
        // untouched keys keep defaults
        assert_eq!(cfg.meshing.sheet_cap, 30);
    }

    #[test]
    fn parses_terminals() {
        let cfg = SolverConfig::from_yaml(
            "terminals:\n  - name: PVIN\n    voltage: 12.0\n    x_um: 10\n    y_um: 20\n    dx_um: 4\n    dy_um: 4\n    layer: Metal2\n  - name: POUT\n    voltage: 0.0\n    x_um: 30\n    y_um: 40\nsolver:\n  missing_model: open\n",
        )
        .unwrap();
        assert_eq!(cfg.terminals.len(), 2);
        assert_eq!(cfg.terminals[0].name, "PVIN");
        assert_eq!(cfg.terminals[0].voltage, 12.0);
        assert_eq!(cfg.terminals[0].layer.as_deref(), Some("Metal2"));
        assert_eq!(cfg.terminals[1].name, "POUT");
        assert_eq!(cfg.terminals[1].dx_um, 0.0);
        assert_eq!(cfg.solver.missing_model, MissingModel::Open);
    }

    #[test]
    fn auto_unit_is_none() {
        let cfg = SolverConfig::from_yaml("units:\n  db_unit_nm: auto\n").unwrap();
        assert_eq!(cfg.units.db_unit_nm, None);
        let cfg = SolverConfig::from_yaml("units:\n  db_unit_nm: 1\n").unwrap();
        assert_eq!(cfg.units.db_unit_nm, Some(1.0));
    }

    #[test]
    fn default_template_round_trips() {
        let cfg = SolverConfig::from_yaml(SolverConfig::default_yaml()).unwrap();
        assert_eq!(cfg.project.layout, "");
        assert_eq!(cfg.meshing.sheet_cell_um, 2.0);
        assert_eq!(cfg.visualization.z_scale, 25.0);
    }
}

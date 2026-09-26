//! Power-device channel modelling for OpenRDSon.
//!
//! - **Mode B (implemented here):** bias-dependent channel/drift resistance from
//!   the channel-model `Id(T, Vgs, Vds)` lookup table.
//! - **Mode A (harness here):** co-simulation with ngspice via a headless
//!   co-process. The deck generator emits the extracted parasitic network plus
//!   the device model, runs `ngspice -b`, and parses the operating point.
//!
//! See `goal.md` for the Mode A/B contract.

use openrdson_core::device::{ChannelElement, SpiInstance};
use std::collections::BTreeMap;
use std::io::Write;
use std::process::Command;

fn err(msg: impl Into<String>) -> String {
    msg.into()
}

/// A regular `Id(T, Vgs, Vds)` lookup table parsed from the model table.
#[derive(Debug, Clone)]
pub struct ModelTable {
    pub temperature_axis: Vec<f64>,
    pub vgs_axis: Vec<f64>,
    pub vds_axis: Vec<f64>,
    /// Reference channel length (m) from the CSV header.
    pub l_ref: f64,
    /// Reference finger width (m) from the CSV header.
    pub w_ref: f64,
    /// Flat `Id` table indexed `[i_vds][i_vgs][i_t]`.
    ids: Vec<Option<f64>>,
    n_t: usize,
    n_vgs: usize,
}

impl ModelTable {
    pub fn from_csv(text: &str) -> Result<Self, String> {
        let mut l_ref = None;
        let mut w_ref = None;
        let mut rows: Vec<(f64, f64, f64, f64)> = Vec::new();
        let mut t_axis: Vec<f64> = Vec::new();
        let mut vgs_axis: Vec<f64> = Vec::new();
        let mut vds_axis: Vec<f64> = Vec::new();

        for line in text.lines() {
            let line = line.trim();
            if line.is_empty() {
                continue;
            }
            if let Some(rest) = line.strip_prefix("L ") {
                l_ref = rest.trim().parse::<f64>().ok();
                continue;
            }
            if let Some(rest) = line.strip_prefix("Wfinger ") {
                w_ref = rest.trim().parse::<f64>().ok();
                continue;
            }
            if line.starts_with("Temperature")
                || line.starts_with("Vgsmax")
                || line.starts_with("Model")
                || line.starts_with("axis names")
            {
                continue;
            }
            let parts: Vec<&str> = line.split(',').collect();
            if parts.len() != 4 {
                continue;
            }
            let vals: Vec<f64> = parts.iter().filter_map(|p| p.trim().parse::<f64>().ok()).collect();
            if vals.len() != 4 {
                continue;
            }
            let (t, vgs, vds, id) = (vals[0], vals[1], vals[2], vals[3]);
            if !t_axis.iter().any(|v| (*v - t).abs() < 1e-12) {
                t_axis.push(t);
            }
            if !vgs_axis.iter().any(|v| (*v - vgs).abs() < 1e-12) {
                vgs_axis.push(vgs);
            }
            if !vds_axis.iter().any(|v| (*v - vds).abs() < 1e-12) {
                vds_axis.push(vds);
            }
            rows.push((t, vgs, vds, id));
        }

        let (n_t, n_vgs, n_vds) = (t_axis.len(), vgs_axis.len(), vds_axis.len());
        if n_t == 0 || n_vgs == 0 || n_vds == 0 {
            return Err("model table had no usable rows".into());
        }
        let mut ids = vec![None; n_t * n_vgs * n_vds];
        let idx_of = |axis: &[f64], v: f64| axis.iter().position(|a| (*a - v).abs() < 1e-12);
        for (t, vgs, vds, id) in rows {
            let it = idx_of(&t_axis, t).unwrap();
            let ig = idx_of(&vgs_axis, vgs).unwrap();
            let iv = idx_of(&vds_axis, vds).unwrap();
            ids[iv * n_t * n_vgs + ig * n_t + it] = Some(id);
        }

        Ok(Self {
            temperature_axis: t_axis,
            vgs_axis,
            vds_axis,
            l_ref: l_ref.unwrap_or(0.4e-6),
            w_ref: w_ref.unwrap_or(1e-6),
            ids,
            n_t,
            n_vgs,
        })
    }

    fn id_at_index(&self, iv: usize, ig: usize, it: usize) -> Option<f64> {
        self.ids[iv * self.n_t * self.n_vgs + ig * self.n_t + it]
    }

    /// Trilinear interpolation of `Id`. Returns `None` if the query is outside
    /// the table or a required corner is missing.
    pub fn id_at(&self, t: f64, vgs: f64, vds: f64) -> Option<f64> {
        let (t0, t1, ft) = bracket(&self.temperature_axis, t)?;
        let (g0, g1, fg) = bracket(&self.vgs_axis, vgs)?;
        let (v0, v1, fv) = bracket(&self.vds_axis, vds)?;
        let mut acc = 0.0;
        for (iv, fv_) in [(v0, 1.0 - fv), (v1, fv)] {
            for (ig, fg_) in [(g0, 1.0 - fg), (g1, fg)] {
                for (it, ft_) in [(t0, 1.0 - ft), (t1, ft)] {
                    let c = self.id_at_index(iv, ig, it)?;
                    acc += c * fv_ * fg_ * ft_;
                }
            }
        }
        Some(acc)
    }

    /// On-resistance `R = Vds / Id` (magnitude). Returns `None` if unavailable,
    /// `Some(inf)` if the device is off.
    pub fn r_on(&self, t: f64, vgs: f64, vds: f64) -> Option<f64> {
        let id = self.id_at(t, vgs, vds)?;
        if id.abs() < 1e-15 {
            return Some(f64::INFINITY);
        }
        Some((vds / id).abs())
    }

    /// Device resistance scaled from the reference finger width to `weff`
    /// (`R ∝ 1/width`).
    pub fn r_device(&self, t: f64, vgs: f64, vds: f64, weff: f64) -> Option<f64> {
        let r = self.r_on(t, vgs, vds)?;
        if weff <= 0.0 {
            return None;
        }
        Some(r * (self.w_ref / weff))
    }
}

/// Find the bracketing interval for `x` in `axis` (ascending or descending).
/// Returns `(i0, i1, frac)` where `frac` is the fraction from `i0` to `i1`.
fn bracket(axis: &[f64], x: f64) -> Option<(usize, usize, f64)> {
    if axis.len() == 1 {
        if (axis[0] - x).abs() < 1e-9 {
            return Some((0, 0, 0.0));
        }
        return None;
    }
    let ascending = axis[1] > axis[0];
    for i in 0..axis.len() - 1 {
        let (a, b) = (axis[i], axis[i + 1]);
        let within = if ascending {
            x >= a - 1e-12 && x <= b + 1e-12
        } else {
            x <= a + 1e-12 && x >= b - 1e-12
        };
        if within {
            let frac = if (b - a).abs() < 1e-300 {
                0.0
            } else {
                (x - a) / (b - a)
            };
            return Some((i, i + 1, frac));
        }
    }
    None
}

/// Parse `key = value` model parameters from a Spectre `.scs` file (e.g.
/// `rdw`, `rdw_inner`, `delta_rdlcw`). Values may use SPICE suffixes.
pub fn parse_scs_params(text: &str) -> BTreeMap<String, f64> {
    let mut out = BTreeMap::new();
    for raw in text.lines() {
        let line = raw.trim().trim_start_matches('+').trim();
        if let Some((k, v)) = line.split_once('=') {
            let key = k.trim().to_string();
            if let Some(val) = parse_spice_number(v.trim()) {
                out.insert(key, val);
            }
        }
    }
    out
}

/// Parse a SPICE-style number with an optional suffix (`k`, `m`, `u`, `n`, `p`).
pub fn parse_spice_number(s: &str) -> Option<f64> {
    let s = s.trim();
    let (num, mult) = if let Some(p) = s.strip_suffix('k').or_else(|| s.strip_suffix('K')) {
        (p, 1e3)
    } else if let Some(p) = s.strip_suffix('u').or_else(|| s.strip_suffix('U')) {
        (p, 1e-6)
    } else if let Some(p) = s.strip_suffix('n').or_else(|| s.strip_suffix('N')) {
        (p, 1e-9)
    } else if let Some(p) = s.strip_suffix('p').or_else(|| s.strip_suffix('P')) {
        (p, 1e-12)
    } else if let Some(p) = s.strip_suffix('m') {
        (p, 1e-3)
    } else {
        (s, 1.0)
    };
    num.trim().parse::<f64>().ok().map(|v| v * mult)
}

/// Equivalent conductivity of a channel of `length` and cross-section `area`
/// that presents resistance `r` (S/m). Used to feed the channel into the FEM.
pub fn conductivity_from_resistance(r: f64, length: f64, area: f64) -> f64 {
    if r <= 0.0 || area <= 0.0 {
        return 0.0;
    }
    length / (r * area)
}

/// Build per-instance channel elements from golden SPI instance parameters.
pub fn channel_elements_from_spi(
    spi: &[SpiInstance],
    table: &ModelTable,
    temperature: f64,
    vgs: f64,
    vds: f64,
) -> Vec<ChannelElement> {
    let mut out = Vec::new();
    for inst in spi {
        // Per-instance channel width = the drawn finger width `w`. The CCI
        // `weff` is the device's total combined finger width (`w · nf`).
        let weff = inst.params.get("w").copied().unwrap_or(table.w_ref);
        if let Some(r) = table.r_device(temperature, vgs, vds, weff) {
            out.push(ChannelElement {
                device_ref: inst.name.clone(),
                temperature,
                vgs,
                vds,
                r_channel: r,
            });
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Mode A: ngspice co-process harness
// ---------------------------------------------------------------------------

/// Headless ngspice runner (Mode A).
#[derive(Debug, Clone)]
pub struct Ngspice {
    pub binary: String,
}

impl Default for Ngspice {
    fn default() -> Self {
        Self {
            binary: "ngspice".to_string(),
        }
    }
}

impl Ngspice {
    /// Whether the configured ngspice binary is runnable.
    pub fn is_available(&self) -> bool {
        Command::new(&self.binary)
            .arg("--version")
            .output()
            .map(|o| o.status.success())
            .unwrap_or(false)
    }

    /// Run a SPICE deck in batch mode and return stdout.
    pub fn run(&self, deck: &str) -> Result<String, String> {
        use std::sync::atomic::{AtomicU64, Ordering};
        static COUNTER: AtomicU64 = AtomicU64::new(0);
        let id = COUNTER.fetch_add(1, Ordering::Relaxed);
        let dir = std::env::temp_dir()
            .join(format!("openrdson_ngspice_{}_{id}", std::process::id()));
        std::fs::create_dir_all(&dir).map_err(|e| err(format!("mkdir: {e}")))?;
        let path = dir.join("deck.cir");
        let mut f = std::fs::File::create(&path).map_err(|e| err(format!("create deck: {e}")))?;
        f.write_all(deck.as_bytes())
            .map_err(|e| err(format!("write deck: {e}")))?;
        let out = Command::new(&self.binary)
            .arg("-b")
            .arg(&path)
            .output()
            .map_err(|e| err(format!("spawn ngspice: {e}")))?;
        if !out.status.success() {
            return Err(format!(
                "ngspice failed ({}): {}",
                out.status,
                String::from_utf8_lossy(&out.stderr)
            ));
        }
        Ok(String::from_utf8_lossy(&out.stdout).to_string())
    }

    /// Run a deck and parse the `.op` node voltages / source currents.
    ///
    /// Parses the tabular section of the ngspice output: lines of the form
    /// `<name> <value>` after the `Node  Voltage` / `Source  Current` headers.
    pub fn operating_point(&self, deck: &str) -> Result<BTreeMap<String, f64>, String> {
        let text = self.run(deck)?;
        Ok(parse_operating_point(&text))
    }
}

/// Parse `name value` pairs from an ngspice batch operating-point output.
pub fn parse_operating_point(text: &str) -> BTreeMap<String, f64> {
    let mut out = BTreeMap::new();
    for line in text.lines() {
        let t: Vec<&str> = line.split_whitespace().collect();
        if t.len() == 2 {
            if let (true, Ok(v)) = (
                t[0].chars().next().map(|c| c.is_alphabetic()).unwrap_or(false),
                t[1].parse::<f64>(),
            ) {
                // Skip the divider rows and header labels.
                if t[0] != "Node" && t[0] != "Source" {
                    out.insert(t[0].to_string(), v);
                }
            }
        }
    }
    out
}

/// Fit a level-1-like channel conductance and threshold from a reference
/// on-resistance, for the ngspice behavioral channel resistor:
/// `R(Vgs) = 1 / (gch · |Vgs − vto|)`.
pub fn fit_behavioral_channel(r_ref: f64, vgs_ref: f64, vto: f64) -> Option<f64> {
    let overdrive = (vgs_ref - vto).abs();
    if r_ref <= 0.0 || overdrive <= 0.0 {
        return None;
    }
    Some(1.0 / (r_ref * overdrive))
}

/// Least-squares fit of `(gch, vto)` for the behavioral channel resistor
/// `R(Vgs) = 1 / (gch · |Vgs − vto|)` from `(vgs, R)` samples.
///
/// For each candidate `vto` the best `gch` is the linear least-squares solution
/// of `1/R = gch·|Vgs−vto|`; the `vto` with the smallest residual is returned.
pub fn fit_behavioral_channel_ls(samples: &[(f64, f64)]) -> Option<(f64, f64)> {
    let good: Vec<(f64, f64)> = samples
        .iter()
        .copied()
        .filter(|(_, r)| r.is_finite() && *r > 0.0)
        .collect();
    if good.len() < 2 {
        return None;
    }
    let vmin = good.iter().map(|(v, _)| *v).fold(f64::INFINITY, f64::min);
    let vmax = good.iter().map(|(v, _)| *v).fold(f64::NEG_INFINITY, f64::max);
    let mut best: Option<(f64, f64, f64)> = None; // (sse, gch, vto)
    let steps = 400;
    for k in 0..=steps {
        let vto = vmin - 0.5 + (vmax - vmin + 1.0) * (k as f64 / steps as f64);
        let mut sxy = 0.0;
        let mut sxx = 0.0;
        for &(vgs, r) in &good {
            let d = (vgs - vto).abs();
            let y = 1.0 / r;
            sxy += d * y;
            sxx += d * d;
        }
        if sxx <= 0.0 {
            continue;
        }
        let gch = sxy / sxx;
        if gch <= 0.0 {
            continue;
        }
        let sse: f64 = good
            .iter()
            .map(|&(vgs, r)| {
                let d = (vgs - vto).abs();
                let pred = 1.0 / (gch * d);
                (pred - r).powi(2)
            })
            .sum();
        if best.map(|(b, _, _)| sse < b).unwrap_or(true) {
            best = Some((sse, gch, vto));
        }
    }
    best.map(|(_, gch, vto)| (gch, vto))
}

/// Solve a 4x4 linear system with partial pivoting.
fn solve4(mut m: [[f64; 4]; 4], mut b: [f64; 4]) -> Option<[f64; 4]> {
    for col in 0..4 {
        let mut piv = col;
        for r in col + 1..4 {
            if m[r][col].abs() > m[piv][col].abs() {
                piv = r;
            }
        }
        if m[piv][col].abs() < 1e-300 {
            return None;
        }
        m.swap(col, piv);
        b.swap(col, piv);
        for r in col + 1..4 {
            let f = m[r][col] / m[col][col];
            for c in col..4 {
                m[r][c] -= f * m[col][c];
            }
            b[r] -= f * b[col];
        }
    }
    let mut x = [0.0f64; 4];
    for i in (0..4).rev() {
        let mut s = b[i];
        for j in i + 1..4 {
            s -= m[i][j] * x[j];
        }
        x[i] = s / m[i][i];
    }
    Some(x)
}

/// Least-squares cubic fit of `1/R = c0 + c1·Vgs + c2·Vgs² + c3·Vgs³` from
/// `(vgs, R)` samples. This matches the device `R(Vgs)` curve far better than a
/// single-overdrive level-1 fit, so the ngspice Mode A co-simulation tracks the
/// model over the whole strong-inversion range.
pub fn fit_channel_polynomial(samples: &[(f64, f64)]) -> Option<[f64; 4]> {
    let good: Vec<(f64, f64)> = samples
        .iter()
        .copied()
        .filter(|(_, r)| r.is_finite() && *r > 0.0)
        .collect();
    if good.len() < 4 {
        return None;
    }
    let mut m = [[0.0f64; 4]; 4];
    let mut b = [0.0f64; 4];
    for &(v, r) in &good {
        let y = 1.0 / r;
        let p = [1.0, v, v * v, v * v * v];
        for i in 0..4 {
            for j in 0..4 {
                m[i][j] += p[i] * p[j];
            }
            b[i] += p[i] * y;
        }
    }
    solve4(m, b)
}

/// Generate a Mode A deck using a cubic polynomial channel resistor.
pub fn mode_a_polynomial_deck(
    r_access: f64,
    c: [f64; 4],
    vds: f64,
    vgs: f64,
) -> String {
    format!(
        "* OpenRDSon Mode A co-simulation (parasitic R + polynomial channel)\n\
         Vd d 0 dc {vds}\n\
         Vg g 0 dc {vgs}\n\
         Racc d dint {racc}\n\
         Rch dint s R='1/(({c0})+({c1})*V(g)+({c2})*V(g)**2+({c3})*V(g)**3)'\n\
         Vs s 0 dc 0\n\
         .op\n\
         .end\n",
        vds = vds,
        vgs = vgs,
        racc = r_access,
        c0 = c[0],
        c1 = c[1],
        c2 = c[2],
        c3 = c[3],
    )
}

/// Generate a Mode A SPICE deck that co-simulates the extracted **parasitic
/// resistance** with a **bias-dependent behavioral channel** in ngspice.
///
/// `r_access` is the metal/contact access resistance (ohms), `gch` and `vto`
/// parameterise the channel resistor. The deck places the access resistance in
/// series with the behavioral channel, so `Rds = Vds / Id`.
pub fn mode_a_behavioral_deck(
    r_access: f64,
    gch: f64,
    vto: f64,
    vds: f64,
    vgs: f64,
) -> String {
    format!(
        "* OpenRDSon Mode A co-simulation (parasitic R + behavioral channel)\n\
         Vd d 0 dc {vds}\n\
         Vg g 0 dc {vgs}\n\
         Racc d dint {racc}\n\
         Rch dint s R='1/({gch}*abs(V(g)-({vto})))'\n\
         Vs s 0 dc 0\n\
         .op\n\
         .end\n",
        vds = vds,
        vgs = vgs,
        racc = r_access,
        gch = gch,
        vto = vto,
    )
}

/// Extract `Rds` from an ngspice operating point by driving `Vd`.
pub fn rds_from_op(op: &BTreeMap<String, f64>, vds: f64) -> Option<f64> {
    let i = op.get("vd#branch")?.abs();
    if i < 1e-30 {
        return None;
    }
    Some(vds.abs() / i)
}

/// Generate a Mode A SPICE deck for one device instance driven at a bias,
/// including a series parasitic resistance.
pub fn generate_mode_a_deck(
    device: &SpiInstance,
    model_include: &str,
    vgs: f64,
    vds: f64,
    r_parasitic: f64,
) -> String {
    let w = device.params.get("w").copied().unwrap_or(10e-6);
    let nf = device.params.get("nf").copied().unwrap_or(1.0);
    let nx = device.params.get("nx").copied().unwrap_or(1.0);
    let ny = device.params.get("ny").copied().unwrap_or(1.0);
    format!(
        "* OpenRDSon Mode A co-simulation deck\n\
         * device {name} model {model}\n\
         .include \"{include}\"\n\
         Vg g 0 dc {vgs}\n\
         Vd dint 0 dc {vds}\n\
         Rpar dint d {rpar}\n\
         Xd d g s s {model} w={w} nf={nf} nx={nx} ny={ny} m=1 mlay=1 avnx=1\n\
         Vs s 0 dc 0\n\
         .op\n\
         .end\n",
        name = device.name,
        model = device.model,
        include = model_include,
        vgs = vgs,
        vds = vds,
        rpar = r_parasitic,
        w = w,
        nf = nf,
        nx = nx,
        ny = ny,
    )
}

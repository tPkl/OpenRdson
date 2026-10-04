//! Power-device channel modelling for OpenRDSon.
//!
//! Bias-dependent channel/drift resistance from the channel-model
//! `Id(T, Vgs, Vds)` lookup table.

use openrdson_core::diag::Diag;

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
    /// Parse the model table. Returns the table plus diagnostics for malformed
    /// rows/headers. Errors only when the table has no usable rows at all.
    pub fn from_csv(text: &str) -> Result<(Self, Vec<Diag>), String> {
        let mut l_ref = None;
        let mut w_ref = None;
        let mut rows: Vec<(f64, f64, f64, f64)> = Vec::new();
        let mut t_axis: Vec<f64> = Vec::new();
        let mut vgs_axis: Vec<f64> = Vec::new();
        let mut vds_axis: Vec<f64> = Vec::new();
        let mut diags = Vec::new();

        for (i, raw) in text.lines().enumerate() {
            let line_no = i + 1;
            let line = raw.trim();
            if line.is_empty() {
                continue;
            }
            if let Some(rest) = line.strip_prefix("L ") {
                l_ref = match rest.trim().parse::<f64>() {
                    Ok(v) => Some(v),
                    Err(_) => {
                        diags.push(Diag::at(
                            line_no,
                            format!("\"L\" expects a number, got \"{}\"", rest.trim()),
                        ));
                        None
                    }
                };
                continue;
            }
            if let Some(rest) = line.strip_prefix("Wfinger ") {
                w_ref = match rest.trim().parse::<f64>() {
                    Ok(v) => Some(v),
                    Err(_) => {
                        diags.push(Diag::at(
                            line_no,
                            format!("\"Wfinger\" expects a number, got \"{}\"", rest.trim()),
                        ));
                        None
                    }
                };
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
                diags.push(Diag::at(
                    line_no,
                    format!(
                        "expected 4 comma-separated values (T,Vgs,Vds,Id), got {}: \"{line}\"",
                        parts.len()
                    ),
                ));
                continue;
            }
            let vals: Vec<f64> = parts.iter().filter_map(|p| p.trim().parse::<f64>().ok()).collect();
            if vals.len() != 4 {
                diags.push(Diag::at(line_no, format!("non-numeric data row: \"{line}\"")));
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

        Ok((
            Self {
                temperature_axis: t_axis,
                vgs_axis,
                vds_axis,
                l_ref: l_ref.unwrap_or(0.4e-6),
                w_ref: w_ref.unwrap_or(1e-6),
                ids,
                n_t,
                n_vgs,
            },
            diags,
        ))
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

#[cfg(test)]
mod tests {
    use super::*;

    const CSV: &str = "Temperature 25\n\
        Vgsmax 0\n\
        L 4.0e-7\n\
        Wfinger 1.0e-6\n\
        Model MY_NMOS\n\
        axis names: temperature,Vgs,Vds,Id\n\
        27.0,0.0,0.1,0.0\n\
        27.0,-1.0,0.1,1.0e-3\n\
        27.0,-2.0,0.1,1.0e-2\n";

    #[test]
    fn parses_a_clean_table_without_diagnostics() {
        let (t, diags) = ModelTable::from_csv(CSV).unwrap();
        assert!(diags.is_empty(), "{diags:?}");
        assert_eq!(t.l_ref, 4.0e-7);
        assert_eq!(t.w_ref, 1.0e-6);
        assert_eq!(t.vgs_axis.len(), 3);
        // Id is read back through the trilinear lookup.
        assert_eq!(t.id_at(27.0, -2.0, 0.1), Some(1.0e-2));
    }

    #[test]
    fn reports_a_row_with_the_wrong_column_count() {
        let (_t, diags) =
            ModelTable::from_csv("27.0,0.0,0.0,0.0\n27.0,0.0,0.0\n").unwrap();
        assert!(
            diags
                .iter()
                .any(|d| d.message.contains("expected 4 comma-separated values")),
            "{diags:?}"
        );
        assert_eq!(diags[0].line, 2);
    }

    #[test]
    fn reports_a_non_numeric_data_row() {
        let (_t, diags) =
            ModelTable::from_csv("27.0,0.0,0.0,0.0\n27.0,x,0.0,0.0\n").unwrap();
        assert!(
            diags.iter().any(|d| d.message.contains("non-numeric data row")),
            "{diags:?}"
        );
    }

    #[test]
    fn reports_a_bad_header_value() {
        let (_t, diags) =
            ModelTable::from_csv("L notanumber\n27.0,0.0,0.0,0.0\n").unwrap();
        assert!(
            diags
                .iter()
                .any(|d| d.message.contains("\"L\" expects a number")),
            "{diags:?}"
        );
    }

    #[test]
    fn errors_when_there_are_no_usable_rows() {
        assert!(ModelTable::from_csv("Temperature 25\nModel X\n").is_err());
    }
}

//! Independent verification of device recognition.
//!
//! Treats the recognizer as a black box and compares its [`DeviceInstance`]s to
//! the golden LVS SPICE netlist. Ground truth for this database
//! is the per-instance `$X/$Y` location (database units), which the recognition
//! stage never consumes.

use openrdson_core::device::{DeviceInstance, SpiInstance};

#[derive(Debug, Clone, PartialEq)]
pub struct InstanceMatch {
    pub golden_name: String,
    pub recognized_name: String,
    pub distance_m: f64,
}

#[derive(Debug, Clone, Default)]
pub struct VerificationReport {
    pub golden_count: usize,
    pub recognized_count: usize,
    pub matched: Vec<InstanceMatch>,
    /// Golden instances with no recognized counterpart.
    pub missing: Vec<String>,
    /// Recognized instances with no golden counterpart.
    pub extra: Vec<String>,
    /// Maximum match distance (m).
    pub max_distance_m: f64,
    /// Whether every matched instance has a determined channel orientation.
    pub all_oriented: bool,
}

impl VerificationReport {
    pub fn passed(&self, tolerance_m: f64) -> bool {
        self.missing.is_empty()
            && self.extra.is_empty()
            && self.matched.len() == self.golden_count
            && self.max_distance_m <= tolerance_m
            && self.all_oriented
    }
}

/// Match recognized instances to golden SPI instances by nearest location.
///
/// `db_unit_meters` converts the golden `$X/$Y` (database units) to meters.
pub fn verify_against_spi(
    recognized: &[DeviceInstance],
    golden: &[SpiInstance],
    db_unit_meters: f64,
) -> VerificationReport {
    let mut report = VerificationReport {
        golden_count: golden.len(),
        recognized_count: recognized.len(),
        ..Default::default()
    };

    let mut used = vec![false; recognized.len()];
    let mut max_d = 0.0f64;

    for g in golden {
        let (Some(gx), Some(gy)) = (g.x, g.y) else {
            report.missing.push(g.name.clone());
            continue;
        };
        let gx = gx * db_unit_meters;
        let gy = gy * db_unit_meters;

        let mut best: Option<(usize, f64)> = None;
        for (i, r) in recognized.iter().enumerate() {
            if used[i] {
                continue;
            }
            let Some((rx, ry)) = r.location else { continue };
            let d = ((rx - gx).powi(2) + (ry - gy).powi(2)).sqrt();
            if best.map(|(_, bd)| d < bd).unwrap_or(true) {
                best = Some((i, d));
            }
        }
        match best {
            Some((i, d)) => {
                used[i] = true;
                max_d = max_d.max(d);
                report.matched.push(InstanceMatch {
                    golden_name: g.name.clone(),
                    recognized_name: recognized[i].name.clone(),
                    distance_m: d,
                });
            }
            None => report.missing.push(g.name.clone()),
        }
    }

    for (i, r) in recognized.iter().enumerate() {
        if !used[i] {
            report.extra.push(r.name.clone());
        }
    }

    report.max_distance_m = max_d;
    report.all_oriented = report
        .matched
        .iter()
        .all(|_| recognized.iter().all(|r| r.channel_axis.is_some()));
    report
}

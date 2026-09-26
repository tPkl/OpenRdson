//! Parser for ICT process-stack files (`*.ict`).
//!
//! The ICT format is a brace-delimited key/value description of conductors,
//! dielectrics, diffusions and vias. Lengths in the file are in micrometres and
//! are converted to metres here (the boundary conversion). Resistivities are
//! kept as the raw ICT table (metals encode width-dependent sheet resistance as
//! `rho width rho width ...`).

use openrdson_core::tech::{Conductor, Dielectric, Diffusion, TechStack, Via};
use std::io;
use std::path::Path;

const UM: f64 = 1e-6;

fn um(v: f64) -> f64 {
    v * UM
}

fn strip_comment(line: &str) -> &str {
    match line.find('#') {
        Some(i) => &line[..i],
        None => line,
    }
}

#[derive(Default)]
struct Block {
    fields: Vec<(String, Vec<String>)>,
}

impl Block {
    fn add(&mut self, key: &str, vals: &[&str]) {
        self.fields
            .push((key.to_string(), vals.iter().map(|s| s.to_string()).collect()));
    }

    fn values(&self, key: &str) -> Option<&Vec<String>> {
        self.fields.iter().find(|(k, _)| k == key).map(|(_, v)| v)
    }

    fn first_float(&self, key: &str) -> Option<f64> {
        self.values(key)?.first()?.parse::<f64>().ok()
    }

    fn first_um(&self, key: &str) -> Option<f64> {
        self.first_float(key).map(um)
    }

    fn floats(&self, key: &str) -> Vec<f64> {
        self.values(key)
            .map(|v| v.iter().filter_map(|s| s.parse::<f64>().ok()).collect())
            .unwrap_or_default()
    }

    fn first_str(&self, key: &str) -> Option<String> {
        self.values(key)?.first().map(|s| s.trim_matches('"').to_string())
    }

    /// Try several spellings (the ICT mixes `expandedFrom` and `expanded_from`).
    fn first_str_any(&self, keys: &[&str]) -> Option<String> {
        keys.iter().find_map(|k| self.first_str(k))
    }

    fn first_um_any(&self, keys: &[&str]) -> Option<f64> {
        keys.iter().find_map(|k| self.first_um(k))
    }

    fn first_bool(&self, key: &str) -> bool {
        matches!(
            self.first_str(key).map(|s| s.to_ascii_lowercase()).as_deref(),
            Some("true") | Some("1") | Some("yes")
        )
    }
}

fn parse_block_start(line: &str) -> Option<(String, String)> {
    if !line.ends_with('{') {
        return None;
    }
    let head = line.trim_end_matches('{').trim();
    let toks: Vec<&str> = head.split_whitespace().collect();
    if toks.len() < 2 {
        return None;
    }
    let kind = toks[0];
    if !matches!(
        kind,
        "conductor" | "dielectric" | "diffusion" | "via" | "process"
    ) {
        return None;
    }
    Some((kind.to_string(), toks[1].trim_matches('"').to_string()))
}

fn finish(stack: &mut TechStack, kind: &str, name: String, b: &Block) {
    match kind {
        "process" => {
            stack.process = name;
            stack.temp_reference = b.first_float("temp_reference");
        }
        "conductor" => stack.conductors.push(Conductor {
            name,
            height: b.first_um("height"),
            thickness: b.first_um("thickness"),
            resistivity: b.floats("resistivity"),
            layer_type: b.first_str("layer_type"),
            gate_forming_layer: b.first_bool("gate_forming_layer"),
            temp_tc1: b.first_float("temp_tc1"),
            temp_tc2: b.first_float("temp_tc2"),
        }),
        "dielectric" => stack.dielectrics.push(Dielectric {
            name,
            height: b.first_um("height"),
            thickness: b.first_um("thickness"),
            dielectric_constant: b.first_float("dielectric_constant"),
            conformal: b.first_bool("conformal"),
            expanded_from: b.first_str_any(&["expanded_from", "expandedFrom"]),
            top_thickness: b.first_um_any(&["top_thickness", "topThickness"]),
            side_expand: b.first_um_any(&["side_expand", "sideExpand"]),
            theta: b.first_float("theta"),
        }),
        "diffusion" => stack.diffusions.push(Diffusion {
            name,
            height: b.first_um("height"),
            thickness: b.first_um("thickness"),
            resistivity: b.floats("resistivity"),
        }),
        "via" => stack.vias.push(Via {
            name,
            top_layer: b.first_str("top_layer"),
            bottom_layer: b.first_str("bottom_layer"),
            area_resistance: b.floats("area_resistance"),
            min_width: b.first_um("min_width"),
            min_spacing: b.first_um("min_spacing"),
            min_top_encl: b.first_um("min_top_encl"),
            min_bot_encl: b.first_um("min_bot_encl"),
            temp_tc1: b.first_float("temp_tc1"),
            temp_tc2: b.first_float("temp_tc2"),
        }),
        _ => {}
    }
}

/// Parse ICT text into a [`TechStack`].
pub fn parse_ict(text: &str) -> TechStack {
    let mut stack = TechStack::default();
    let mut block: Option<(String, String)> = None;
    let mut cur = Block::default();

    for raw in text.lines() {
        let line = strip_comment(raw).trim();
        if line.is_empty() {
            continue;
        }
        if line == "}" || line.starts_with('}') {
            if let Some((kind, name)) = block.take() {
                finish(&mut stack, &kind, name, &cur);
                cur = Block::default();
            }
            continue;
        }
        if let Some((kind, name)) = parse_block_start(line) {
            block = Some((kind, name));
            cur = Block::default();
            continue;
        }
        if block.is_some() {
            let toks: Vec<&str> = line.split_whitespace().collect();
            if !toks.is_empty() {
                cur.add(toks[0], &toks[1..]);
            }
        }
    }

    stack
}

pub fn read_ict_file<P: AsRef<Path>>(path: P) -> io::Result<TechStack> {
    let text = std::fs::read_to_string(path)?;
    Ok(parse_ict(&text))
}

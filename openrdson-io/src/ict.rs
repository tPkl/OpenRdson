//! Parser for ICT process-stack files (`*.ict`).
//!
//! The ICT format is a brace-delimited key/value description of conductors,
//! dielectrics, diffusions and vias. Lengths in the file are in micrometres and
//! are converted to metres here (the boundary conversion). Resistivities are
//! kept as the raw ICT table (metals encode width-dependent sheet resistance as
//! `rho width rho width ...`).
//!
//! [`parse_ict`] returns the parsed [`TechStack`] together with a list of
//! [`Diag`]s describing syntax/semantics problems found while reading, so the
//! `read_*` wrapper can surface them to the user instead of silently dropping
//! malformed input.

use openrdson_core::diag::Diag;
use openrdson_core::tech::{Conductor, Dielectric, Diffusion, TechStack, Via};
use std::collections::BTreeMap;
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
    /// 1-based line of the block header, for diagnostics.
    start_line: usize,
    kind: String,
    name: String,
    fields: Vec<(String, Vec<String>)>,
}

impl Block {
    fn new(start_line: usize, kind: String, name: String) -> Self {
        Self {
            start_line,
            kind,
            name,
            fields: Vec::new(),
        }
    }

    fn add(&mut self, key: &str, vals: &[&str]) {
        self.fields
            .push((key.to_string(), vals.iter().map(|s| s.to_string()).collect()));
    }

    fn values(&self, key: &str) -> Option<&Vec<String>> {
        self.fields.iter().find(|(k, _)| k == key).map(|(_, v)| v)
    }

    /// A numeric field; if the key is present but the value is not a number,
    /// record a diagnostic and return `None`.
    fn num(&self, key: &str, diags: &mut Vec<Diag>) -> Option<f64> {
        let first = self.values(key)?.first()?;
        match first.parse::<f64>() {
            Ok(v) => Some(v),
            Err(_) => {
                diags.push(Diag::at(
                    self.start_line,
                    format!(
                        "{} \"{}\": \"{}\" expects a number, got \"{}\"",
                        self.kind, self.name, key, first
                    ),
                ));
                None
            }
        }
    }

    fn num_um(&self, key: &str, diags: &mut Vec<Diag>) -> Option<f64> {
        self.num(key, diags).map(um)
    }

    fn num_um_any(&self, keys: &[&str], diags: &mut Vec<Diag>) -> Option<f64> {
        for k in keys {
            if let Some(v) = self.num_um(k, diags) {
                return Some(v);
            }
        }
        None
    }

    /// A numeric list; non-numeric entries are reported and dropped.
    fn num_list(&self, key: &str, diags: &mut Vec<Diag>) -> Vec<f64> {
        let Some(vals) = self.values(key) else {
            return Vec::new();
        };
        let mut out = Vec::with_capacity(vals.len());
        for v in vals {
            match v.parse::<f64>() {
                Ok(x) => out.push(x),
                Err(_) => diags.push(Diag::at(
                    self.start_line,
                    format!(
                        "{} \"{}\": \"{}\" expects a number, got \"{}\"",
                        self.kind, self.name, key, v
                    ),
                )),
            }
        }
        out
    }

    fn first_str(&self, key: &str) -> Option<String> {
        self.values(key)?.first().map(|s| s.trim_matches('"').to_string())
    }

    /// Try several spellings (the ICT mixes `expandedFrom` and `expanded_from`).
    fn first_str_any(&self, keys: &[&str]) -> Option<String> {
        keys.iter().find_map(|k| self.first_str(k))
    }

    fn first_bool(&self, key: &str) -> bool {
        matches!(
            self.first_str(key).map(|s| s.to_ascii_lowercase()).as_deref(),
            Some("true") | Some("1") | Some("yes")
        )
    }
}

/// Parse `<kind> "<name>" {`. Returns `None` when the head is not a recognised
/// block kind (the caller reports it).
fn parse_block_head(line: &str) -> Option<(String, String)> {
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

fn finish(stack: &mut TechStack, b: &Block, diags: &mut Vec<Diag>) {
    let ln = b.start_line;
    let name = b.name.clone();
    match b.kind.as_str() {
        "process" => {
            stack.process = name;
            stack.temp_reference = b.num("temp_reference", diags);
        }
        "conductor" => {
            let resistivity = b.num_list("resistivity", diags);
            if resistivity.is_empty() {
                diags.push(Diag::at(
                    ln,
                    format!("conductor \"{name}\" has no resistivity (sheet resistance is unknown)"),
                ));
            } else if resistivity.len() >= 3 && resistivity.len() % 2 == 1 {
                diags.push(Diag::at(
                    ln,
                    format!(
                        "conductor \"{name}\" resistivity has {} values; a width table needs rho/width pairs",
                        resistivity.len()
                    ),
                ));
            }
            stack.conductors.push(Conductor {
                name,
                height: b.num_um("height", diags),
                thickness: b.num_um("thickness", diags),
                resistivity,
                layer_type: b.first_str("layer_type"),
                gate_forming_layer: b.first_bool("gate_forming_layer"),
                temp_tc1: b.num("temp_tc1", diags),
                temp_tc2: b.num("temp_tc2", diags),
            });
        }
        "dielectric" => stack.dielectrics.push(Dielectric {
            name,
            height: b.num_um("height", diags),
            thickness: b.num_um("thickness", diags),
            dielectric_constant: b.num("dielectric_constant", diags),
            conformal: b.first_bool("conformal"),
            expanded_from: b.first_str_any(&["expanded_from", "expandedFrom"]),
            top_thickness: b.num_um_any(&["top_thickness", "topThickness"], diags),
            side_expand: b.num_um_any(&["side_expand", "sideExpand"], diags),
            theta: b.num("theta", diags),
        }),
        "diffusion" => {
            let resistivity = b.num_list("resistivity", diags);
            if resistivity.is_empty() {
                diags.push(Diag::at(
                    ln,
                    format!("diffusion \"{name}\" has no resistivity (sheet resistance is unknown)"),
                ));
            }
            stack.diffusions.push(Diffusion {
                name,
                height: b.num_um("height", diags),
                thickness: b.num_um("thickness", diags),
                resistivity,
            });
        }
        "via" => {
            let top_layer = b.first_str("top_layer");
            let bottom_layer = b.first_str("bottom_layer");
            let area_resistance = b.num_list("area_resistance", diags);
            if top_layer.is_none() {
                diags.push(Diag::at(ln, format!("via \"{name}\" has no top_layer")));
            }
            if bottom_layer.is_none() {
                diags.push(Diag::at(ln, format!("via \"{name}\" has no bottom_layer")));
            }
            if area_resistance.is_empty() {
                diags.push(Diag::at(ln, format!("via \"{name}\" has no area_resistance")));
            }
            stack.vias.push(Via {
                name,
                top_layer,
                bottom_layer,
                area_resistance,
                min_width: b.num_um("min_width", diags),
                min_spacing: b.num_um("min_spacing", diags),
                min_top_encl: b.num_um("min_top_encl", diags),
                min_bot_encl: b.num_um("min_bot_encl", diags),
                temp_tc1: b.num("temp_tc1", diags),
                temp_tc2: b.num("temp_tc2", diags),
            });
        }
        _ => {}
    }
}

/// Stack-level connectivity checks.
///
/// A real process stack must conduct: every conducting layer (metal, poly,
/// diffusion, pad) should be reachable through the via graph, and every via
/// should reference real layers. This flags isolated layers and broken
/// conduction paths, which otherwise show up only as a mysterious open circuit
/// during extraction.
fn check_stack_connectivity(
    stack: &TechStack,
    layer_lines: &BTreeMap<String, usize>,
    via_lines: &BTreeMap<String, usize>,
    diags: &mut Vec<Diag>,
) {
    use std::collections::BTreeSet;

    // Names that count as conducting nodes (metals + diffusions).
    let conducting: BTreeSet<&str> = stack
        .conductors
        .iter()
        .map(|c| c.name.as_str())
        .chain(stack.diffusions.iter().map(|d| d.name.as_str()))
        .collect();
    if conducting.is_empty() {
        return;
    }

    // (1) Every via endpoint must exist; a via cannot jump between equal layers.
    let mut touched: BTreeSet<String> = BTreeSet::new();
    for v in &stack.vias {
        let ln = via_lines.get(&v.name).copied().unwrap_or(0);
        for (label, layer) in [("top_layer", &v.top_layer), ("bottom_layer", &v.bottom_layer)] {
            if let Some(l) = layer {
                if conducting.contains(l.as_str()) {
                    touched.insert(l.clone());
                } else {
                    diags.push(Diag::at(
                        ln,
                        format!(
                            "via \"{}\" {} \"{}\" is not a conducting layer (conductor/diffusion) in this stack",
                            v.name, label, l
                        ),
                    ));
                }
            }
        }
        if let (Some(t), Some(b)) = (v.top_layer.as_deref(), v.bottom_layer.as_deref()) {
            if t == b {
                diags.push(Diag::at(
                    ln,
                    format!(
                        "via \"{}\" has the same top_layer and bottom_layer (\"{t}\")",
                        v.name
                    ),
                ));
            }
        }
    }

    // (2) Every conducting layer should be touched by at least one via.
    for c in &stack.conductors {
        if !touched.contains(&c.name) {
            let ln = layer_lines.get(&c.name).copied().unwrap_or(0);
            diags.push(Diag::at(
                ln,
                format!(
                    "conductor \"{}\" is not connected to any via/contact (no conduction path in or out of this layer)",
                    c.name
                ),
            ));
        }
    }
    for d in &stack.diffusions {
        if !touched.contains(&d.name) {
            let ln = layer_lines.get(&d.name).copied().unwrap_or(0);
            diags.push(Diag::at(
                ln,
                format!(
                    "diffusion \"{}\" is not connected to any via/contact",
                    d.name
                ),
            ));
        }
    }

    // (3) The via graph over the conducting layers must be a single component.
    let mut adj: BTreeMap<&str, Vec<&str>> = BTreeMap::new();
    for v in &stack.vias {
        if let (Some(t), Some(b)) = (v.top_layer.as_deref(), v.bottom_layer.as_deref()) {
            if t != b && conducting.contains(t) && conducting.contains(b) {
                adj.entry(t).or_default().push(b);
                adj.entry(b).or_default().push(t);
            }
        }
    }
    let mut comp: BTreeMap<&str, usize> = BTreeMap::new();
    let mut ncomp = 0usize;
    for name in &conducting {
        if comp.contains_key(name) {
            continue;
        }
        let cid = ncomp;
        ncomp += 1;
        comp.insert(name, cid);
        let mut queue = vec![*name];
        while let Some(x) = queue.pop() {
            if let Some(nbrs) = adj.get(x) {
                for &y in nbrs {
                    if !comp.contains_key(y) {
                        comp.insert(y, cid);
                        queue.push(y);
                    }
                }
            }
        }
    }
    if ncomp > 1 {
        let mut groups: Vec<Vec<&str>> = vec![Vec::new(); ncomp];
        for (name, &cid) in &comp {
            groups[cid].push(name);
        }
        for g in &mut groups {
            g.sort_unstable();
        }
        let desc = groups
            .iter()
            .map(|g| format!("[{}]", g.join(", ")))
            .collect::<Vec<_>>()
            .join("  ");
        diags.push(Diag::file(format!(
            "the conduction path is broken: the via graph splits the stack into {ncomp} disconnected groups: {desc}"
        )));
    }
}

/// Parse ICT text into a [`TechStack`] plus diagnostics.
pub fn parse_ict(text: &str) -> (TechStack, Vec<Diag>) {
    let mut stack = TechStack::default();
    let mut diags = Vec::new();
    let mut block: Option<Block> = None;
    let mut layer_lines: BTreeMap<String, usize> = BTreeMap::new();
    let mut via_lines: BTreeMap<String, usize> = BTreeMap::new();

    for (i, raw) in text.lines().enumerate() {
        let line_no = i + 1;
        let line = strip_comment(raw).trim();
        if line.is_empty() {
            continue;
        }
        if line == "}" || line.starts_with('}') {
            match block.take() {
                Some(b) => {
                    finish(&mut stack, &b, &mut diags);
                    match b.kind.as_str() {
                        "conductor" | "diffusion" => {
                            layer_lines.insert(b.name.clone(), b.start_line);
                        }
                        "via" => {
                            via_lines.insert(b.name.clone(), b.start_line);
                        }
                        _ => {}
                    }
                }
                None => diags.push(Diag::at(line_no, "unexpected '}' with no open block")),
            }
            continue;
        }
        if line.ends_with('{') {
            if block.is_some() {
                diags.push(Diag::at(
                    line_no,
                    "block opened before the previous block was closed",
                ));
            }
            match parse_block_head(line) {
                Some((kind, name)) => block = Some(Block::new(line_no, kind, name)),
                None => {
                    let kind = line
                        .trim_end_matches('{')
                        .split_whitespace()
                        .next()
                        .unwrap_or("?");
                    diags.push(Diag::at(
                        line_no,
                        format!(
                            "unrecognized block type \"{kind}\" (expected conductor/dielectric/diffusion/via/process)"
                        ),
                    ));
                    block = None;
                }
            }
            continue;
        }
        match block.as_mut() {
            Some(b) => {
                let toks: Vec<&str> = line.split_whitespace().collect();
                if !toks.is_empty() {
                    b.add(toks[0], &toks[1..]);
                }
            }
            None => diags.push(Diag::at(
                line_no,
                format!("field outside any block: \"{line}\""),
            )),
        }
    }
    if let Some(b) = block {
        diags.push(Diag::at(
            b.start_line,
            format!(
                "unterminated block {} \"{}\" (missing closing '}}')",
                b.kind, b.name
            ),
        ));
    }
    check_stack_connectivity(&stack, &layer_lines, &via_lines, &mut diags);
    (stack, diags)
}

pub fn read_ict_file<P: AsRef<Path>>(path: P) -> io::Result<TechStack> {
    let path = path.as_ref();
    let text = std::fs::read_to_string(path)?;
    let (stack, diags) = parse_ict(&text);
    openrdson_core::diag::log_diags(&path.display().to_string(), &diags);
    Ok(stack)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_a_minimal_stack_without_diagnostics() {
        let (stack, diags) = parse_ict(
            r#"process "p1" {
                temp_reference 25
            }
            conductor "M1" {
                height 0.5
                thickness 0.3
                resistivity 0.0822 0.16 0.06 40.0
                temp_tc1 0.003
            }
            conductor "M2" {
                resistivity 0.06
            }
            via "V1" {
                top_layer "M2"
                bottom_layer "M1"
                area_resistance 0.9 0.04
            }
            "#,
        );
        assert!(diags.is_empty(), "unexpected diagnostics: {diags:?}");
        assert_eq!(stack.process, "p1");
        assert_eq!(stack.temp_reference, Some(25.0));
        assert_eq!(stack.conductors.len(), 2);
        assert_eq!(stack.conductors[0].name, "M1");
        assert_eq!(stack.conductors[0].height, Some(0.5e-6));
        // width-dependent table preserved in raw order
        assert_eq!(stack.conductors[0].resistivity, vec![0.0822, 0.16, 0.06, 40.0]);
        assert_eq!(stack.vias.len(), 1);
        assert_eq!(stack.vias[0].top_layer.as_deref(), Some("M2"));
        assert_eq!(stack.vias[0].bottom_layer.as_deref(), Some("M1"));
        assert_eq!(stack.vias[0].area_resistance, vec![0.9, 0.04]);
    }

    #[test]
    fn flags_a_conductor_without_resistivity() {
        let (_s, diags) = parse_ict("conductor \"M1\" {\n  height 0.5\n}\n");
        assert!(
            diags.iter().any(|d| d.message.contains("no resistivity")),
            "expected a resistivity diagnostic, got {diags:?}"
        );
        assert_eq!(diags[0].line, 1);
    }

    #[test]
    fn flags_a_via_missing_its_layers_and_area() {
        let (_s, diags) = parse_ict("via \"V1\" {\n  min_width 0.1\n}\n");
        let msgs: Vec<_> = diags.iter().map(|d| d.message.as_str()).collect();
        assert!(msgs.iter().any(|m| m.contains("no top_layer")), "{diags:?}");
        assert!(msgs.iter().any(|m| m.contains("no bottom_layer")), "{diags:?}");
        assert!(msgs.iter().any(|m| m.contains("no area_resistance")), "{diags:?}");
    }

    #[test]
    fn flags_an_odd_resistivity_table() {
        let (_s, diags) = parse_ict("conductor \"M1\" {\n  resistivity 0.1 0.2 0.3\n}\n");
        assert!(
            diags.iter().any(|d| d.message.contains("rho/width pairs")),
            "{diags:?}"
        );
    }

    #[test]
    fn flags_non_numeric_values() {
        let (_s, diags) = parse_ict("conductor \"M1\" {\n  height abc\n  resistivity 0.1\n}\n");
        assert!(
            diags
                .iter()
                .any(|d| d.message.contains("\"height\" expects a number")),
            "{diags:?}"
        );
        assert_eq!(diags[0].line, 1);
    }

    #[test]
    fn flags_unterminated_and_unknown_and_stray_braces() {
        let (_s, diags) = parse_ict("conductor \"M1\" {\n  resistivity 0.1\n");
        assert!(
            diags.iter().any(|d| d.message.contains("unterminated block")),
            "{diags:?}"
        );

        let (_s, diags) = parse_ict("gizmo \"X\" {\n  x 1\n}\n");
        assert!(
            diags
                .iter()
                .any(|d| d.message.contains("unrecognized block type")),
            "{diags:?}"
        );

        let (_s, diags) = parse_ict("}\n");
        assert!(
            diags.iter().any(|d| d.message.contains("unexpected '}'")),
            "{diags:?}"
        );

        let (_s, diags) = parse_ict("stray 1\n");
        assert!(
            diags
                .iter()
                .any(|d| d.message.contains("field outside any block")),
            "{diags:?}"
        );
    }

    #[test]
    fn comments_and_blank_lines_are_ignored() {
        let (stack, diags) = parse_ict(
            "# a comment\n\n\
             conductor \"M1\" {  # trailing\n  resistivity 1.0\n}\n\
             conductor \"M2\" {\n  resistivity 1.0\n}\n\
             via \"V1\" {\n  top_layer \"M2\"\n  bottom_layer \"M1\"\n  area_resistance 1.0 0.25\n}\n",
        );
        assert!(diags.is_empty(), "{diags:?}");
        assert_eq!(stack.conductors[0].resistivity, vec![1.0]);
    }

    #[test]
    fn flags_a_conductor_without_a_via() {
        let (_s, diags) = parse_ict("conductor \"M1\" {\n  resistivity 1.0\n}\n");
        assert!(
            diags
                .iter()
                .any(|d| d.message.contains("not connected to any via/contact")),
            "{diags:?}"
        );
        assert_eq!(diags[0].line, 1);
    }

    #[test]
    fn flags_a_via_to_an_unknown_layer() {
        let (_s, diags) = parse_ict(
            "conductor \"M1\" {\n  resistivity 1.0\n}\n\
             via \"V1\" {\n  top_layer \"M9\"\n  bottom_layer \"M1\"\n  area_resistance 1.0 0.25\n}\n",
        );
        assert!(
            diags
                .iter()
                .any(|d| d.message.contains("\"M9\" is not a conducting layer")),
            "{diags:?}"
        );
    }

    #[test]
    fn flags_a_broken_conduction_path() {
        // Two internally-connected sub-stacks (M1-M2 and M3-M4) that are not
        // joined, so no current path crosses the whole stack.
        let (_s, diags) = parse_ict(
            "conductor \"M1\" {\n  resistivity 1.0\n}\n\
             conductor \"M2\" {\n  resistivity 1.0\n}\n\
             conductor \"M3\" {\n  resistivity 1.0\n}\n\
             conductor \"M4\" {\n  resistivity 1.0\n}\n\
             via \"V12\" {\n  top_layer \"M2\"\n  bottom_layer \"M1\"\n  area_resistance 1.0 0.25\n}\n\
             via \"V34\" {\n  top_layer \"M4\"\n  bottom_layer \"M3\"\n  area_resistance 1.0 0.25\n}\n",
        );
        assert!(
            diags
                .iter()
                .any(|d| d.message.contains("conduction path is broken")),
            "{diags:?}"
        );
        // The two groups are named in the message.
        let msg = diags
            .iter()
            .find(|d| d.message.contains("conduction path is broken"))
            .unwrap();
        assert!(msg.message.contains("M1") && msg.message.contains("M3"), "{msg}");
    }
}

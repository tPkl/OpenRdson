//! Parser for the CCI `sft_cci.map` logical→physical layer map.
//!
//! The file is sectioned: entries under `conducting_layers` map a logical
//! conductor to a physical ICT layer, entries under `via_layers` map a logical
//! via. [`parse_sft_cci_map`] returns the map plus diagnostics for malformed or
//! misplaced entries.

use openrdson_core::diag::Diag;
use openrdson_core::tech::LayerMap;
use std::io;
use std::path::Path;

/// Parse `sft_cci.map` into a [`LayerMap`] plus diagnostics.
pub fn parse_sft_cci_map(text: &str) -> (LayerMap, Vec<Diag>) {
    let mut map = LayerMap::default();
    let mut diags = Vec::new();
    #[derive(PartialEq)]
    enum Section {
        None,
        Conducting,
        Via,
    }
    let mut section = Section::None;
    for (i, raw) in text.lines().enumerate() {
        let line_no = i + 1;
        let line = raw.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        match line {
            "conducting_layers" => {
                section = Section::Conducting;
                continue;
            }
            "via_layers" => {
                section = Section::Via;
                continue;
            }
            _ => {}
        }
        let t: Vec<&str> = line.split_whitespace().collect();
        if t.len() < 2 {
            diags.push(Diag::at(
                line_no,
                format!("expected \"<logical> <physical>\", got \"{line}\""),
            ));
            continue;
        }
        let (logical, physical) = (t[0].to_string(), t[1].to_string());
        match section {
            Section::Conducting => {
                map.conducting.insert(logical, physical);
            }
            Section::Via => {
                map.via.insert(logical, physical);
            }
            Section::None => diags.push(Diag::at(
                line_no,
                format!(
                    "\"{logical} {physical}\" appears before a 'conducting_layers' or 'via_layers' header"
                ),
            )),
        }
    }
    (map, diags)
}

pub fn read_sft_cci_map<P: AsRef<Path>>(path: P) -> io::Result<LayerMap> {
    let path = path.as_ref();
    let text = std::fs::read_to_string(path)?;
    let (map, diags) = parse_sft_cci_map(&text);
    openrdson_core::diag::log_diags(&path.display().to_string(), &diags);
    Ok(map)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_sections_without_diagnostics() {
        let (map, diags) = parse_sft_cci_map(
            "conducting_layers\nnet_M1 M1\nnet_M2 M2\nvia_layers\nvia1 V1\n",
        );
        assert!(diags.is_empty(), "{diags:?}");
        assert_eq!(map.conducting.get("net_M1").map(String::as_str), Some("M1"));
        assert_eq!(map.conducting.get("net_M2").map(String::as_str), Some("M2"));
        assert_eq!(map.via.get("via1").map(String::as_str), Some("V1"));
        assert!(!map.conducting.contains_key("via1"));
    }

    #[test]
    fn flags_entries_before_a_section_header() {
        let (map, diags) = parse_sft_cci_map("net_M1 M1\nconducting_layers\nnet_M2 M2\n");
        assert!(
            diags
                .iter()
                .any(|d| d.message.contains("before a 'conducting_layers'")),
            "{diags:?}"
        );
        assert_eq!(diags[0].line, 1);
        // still parsed the valid entry after the header
        assert_eq!(map.conducting.get("net_M2").map(String::as_str), Some("M2"));
    }

    #[test]
    fn flags_short_lines() {
        let (_m, diags) = parse_sft_cci_map("conducting_layers\nonlyone\n");
        assert!(
            diags.iter().any(|d| d.message.contains("expected")),
            "{diags:?}"
        );
        assert_eq!(diags[0].line, 2);
    }

    #[test]
    fn ignores_comments_and_blanks() {
        let (map, diags) = parse_sft_cci_map("# c\n\nconducting_layers\nnet_M1 M1\n");
        assert!(diags.is_empty(), "{diags:?}");
        assert_eq!(map.conducting.len(), 1);
    }
}

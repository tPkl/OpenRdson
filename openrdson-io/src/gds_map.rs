//! Parser for the CCI `*.gds.map` layer-number -> logical-layer-name file.
//!
//! Lines:
//! ```text
//! Gds_Map 1000
//! Layers:
//! 0 0 65 Aug 12 13:11:17 2026
//! net_M1 41 0
//! seed_device 45 0
//! END OF RESPONSE
//! ```
//!
//! Data lines are exactly `<name> <layer> <datatype>`. Header/timestamp lines
//! have a different token count and are skipped. [`parse_gds_map`] returns the
//! map plus diagnostics for data-looking lines with a bad layer/datatype.

use openrdson_core::diag::Diag;
use std::collections::BTreeMap;
use std::io;
use std::path::Path;

/// Parse the map text into `layer number -> logical name` plus diagnostics.
pub fn parse_gds_map(text: &str) -> (BTreeMap<i16, String>, Vec<Diag>) {
    let mut map = BTreeMap::new();
    let mut diags = Vec::new();
    for (i, raw) in text.lines().enumerate() {
        let line_no = i + 1;
        let line = raw.trim();
        if line.is_empty() {
            continue;
        }
        let toks: Vec<&str> = line.split_whitespace().collect();
        // Data lines are exactly three tokens; everything else (the
        // `Gds_Map`/`Layers:`/timestamp/`END OF RESPONSE` lines) has a
        // different shape and is not validated.
        if toks.len() != 3 {
            continue;
        }
        if line.starts_with("END OF RESPONSE") {
            continue;
        }
        let (name, layer_s, datatype_s) = (toks[0], toks[1], toks[2]);
        match (layer_s.parse::<i16>(), datatype_s.parse::<i16>()) {
            (Ok(layer), Ok(_datatype)) => {
                map.insert(layer, name.to_string());
            }
            _ => diags.push(Diag::at(
                line_no,
                format!("expected \"<name> <layer> <datatype>\", got \"{line}\""),
            )),
        }
    }
    (map, diags)
}

pub fn read_gds_map_file<P: AsRef<Path>>(path: P) -> io::Result<BTreeMap<i16, String>> {
    let path = path.as_ref();
    let text = std::fs::read_to_string(path)?;
    let (map, diags) = parse_gds_map(&text);
    openrdson_core::diag::log_diags(&path.display().to_string(), &diags);
    Ok(map)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_a_typical_map_without_diagnostics() {
        let text = "Gds_Map 1000\nLayers:\n0 0 65 Aug 12 13:11:17 2026\n\
                    net_M1 41 0\nseed_device 45 0\nEND OF RESPONSE\n";
        let (map, diags) = parse_gds_map(text);
        assert!(diags.is_empty(), "{diags:?}");
        assert_eq!(map.get(&41).map(String::as_str), Some("net_M1"));
        assert_eq!(map.get(&45).map(String::as_str), Some("seed_device"));
    }

    #[test]
    fn flags_a_data_line_with_a_bad_layer() {
        let (_m, diags) = parse_gds_map("net_M1 notanumber 0\n");
        assert!(
            diags
                .iter()
                .any(|d| d.message.contains("expected \"<name> <layer> <datatype>\"")),
            "{diags:?}"
        );
        assert_eq!(diags[0].line, 1);
    }

    #[test]
    fn header_and_timestamp_lines_are_not_flagged() {
        let (_m, diags) =
            parse_gds_map("Gds_Map 1000\nLayers:\n0 0 65 Aug 12 13:11:17 2026\nEND OF RESPONSE\n");
        assert!(diags.is_empty(), "{diags:?}");
    }
}

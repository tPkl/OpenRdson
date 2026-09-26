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

use std::collections::BTreeMap;
use std::io;
use std::path::Path;

/// Parse the map text into `layer number -> logical name`.
pub fn parse_gds_map(text: &str) -> BTreeMap<i16, String> {
    let mut map = BTreeMap::new();
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let toks: Vec<&str> = line.split_whitespace().collect();
        // Data lines are exactly `<name> <layer> <datatype>`; header/timestamp
        // lines have a different token count and are skipped.
        if toks.len() != 3 {
            continue;
        }
        let (name, layer_s, datatype_s) = (toks[0], toks[1], toks[2]);
        let (Ok(layer), Ok(_datatype)) = (layer_s.parse::<i16>(), datatype_s.parse::<i16>()) else {
            continue;
        };
        map.insert(layer, name.to_string());
    }
    map
}

pub fn read_gds_map_file<P: AsRef<Path>>(path: P) -> io::Result<BTreeMap<i16, String>> {
    let text = std::fs::read_to_string(path)?;
    Ok(parse_gds_map(&text))
}

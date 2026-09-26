//! Parser for the CCI `sft_cci.map` logical→physical layer map.

use openrdson_core::tech::LayerMap;
use std::io;
use std::path::Path;

/// Parse `sft_cci.map` into a [`LayerMap`].
pub fn parse_sft_cci_map(text: &str) -> LayerMap {
    let mut map = LayerMap::default();
    #[derive(PartialEq)]
    enum Section {
        None,
        Conducting,
        Via,
    }
    let mut section = Section::None;
    for line in text.lines() {
        let line = line.trim();
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
            Section::None => {}
        }
    }
    map
}

pub fn read_sft_cci_map<P: AsRef<Path>>(path: P) -> io::Result<LayerMap> {
    let text = std::fs::read_to_string(path)?;
    Ok(parse_sft_cci_map(&text))
}

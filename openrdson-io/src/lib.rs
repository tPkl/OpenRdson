//! File-format readers/writers for OpenRDSon.
//!
//! All parsers are dependency-free and deterministic. Errors carry enough
//! context (record offset / line number) to diagnose malformed inputs.

pub mod gds_map;
pub mod gdsii;
pub mod ict;
pub mod sft_cci;
pub mod viz;
pub mod vtu;

pub use gdsii::{read_gds_file, GdsElement, GdsLibrary, GdsStruct, ElementKind};
pub use gds_map::{read_gds_map_file, parse_gds_map};
pub use ict::{read_ict_file, parse_ict};
pub use sft_cci::{parse_sft_cci_map, read_sft_cci_map};
pub use viz::{
    write_colormap, write_gds_file, write_grouped_colormap, write_lyp, write_named_lyp,
    GdsBoundary, VizItem,
};
pub use vtu::{write_vtu, VtkCellType, VtuMesh};

//! Shared intermediate-representation (IR) types for the OpenRDSon pipeline.
//!
//! Every crate exchanges data through these types. Keep this crate
//! dependency-free and stable.

pub mod device;
pub mod diag;
pub mod geometry;
pub mod geomops;
pub mod log;
pub mod layout;
pub mod mesh;
pub mod tech;

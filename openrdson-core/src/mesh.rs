//! Mesh IR: the interface between `meshing` and `solver`.
//!
//! v0 is intentionally minimal; it is extended as meshing lands.

use crate::geometry::Point;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ElementType {
    Tet,
    Hex,
    Prism,
    Pyramid,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Element {
    pub kind: ElementType,
    /// Indices into [`Mesh::nodes`].
    pub nodes: Vec<u32>,
    pub material_id: u32,
    pub net: Option<String>,
    pub layer: Option<String>,
    pub device_ref: Option<String>,
    /// Index of the source layout polygon (connectivity provenance).
    pub source_polygon: Option<usize>,
}

/// A conforming finite-element mesh.
#[derive(Debug, Clone, Default)]
pub struct Mesh {
    /// Node coordinates (meters).
    pub nodes: Vec<Point3>,
    pub elements: Vec<Element>,
    /// Named node/element sets (regions, boundaries, terminals).
    pub sets: std::collections::BTreeMap<String, Vec<u32>>,
}

/// A 3D point in meters.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Point3 {
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }
}

impl From<Point> for Point3 {
    fn from(p: Point) -> Self {
        Point3::new(p.x, p.y, 0.0)
    }
}

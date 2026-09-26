//! Layout IR: the connectivity-aware view of a layout database.

use crate::geometry::{Bbox, LayerKey, Point};
use std::collections::BTreeMap;

/// A single conducting (or annotated) polygon on a layer.
#[derive(Debug, Clone, PartialEq)]
pub struct LayoutPolygon {
    pub layer: LayerKey,
    /// Logical net this polygon belongs to, if resolved (from `*.gds.map`).
    pub net: Option<String>,
    /// Polygon outline in meters (SI).
    pub points: Vec<Point>,
    /// GDSII element properties (`PROPATTR`, `PROPVALUE`) attached to the polygon.
    pub props: Vec<(i16, String)>,
}

/// A named terminal / port of the top cell.
#[derive(Debug, Clone, PartialEq)]
pub struct Port {
    pub name: String,
    /// LVS node id from the ports table (matches the SPI instance node numbers).
    pub id: Option<i64>,
    pub net: String,
    pub layer: String,
    pub x: f64,
    pub y: f64,
    pub direction: Option<String>,
}

/// The unified layout model produced by `layout-db`.
#[derive(Debug, Clone, Default)]
pub struct LayoutIR {
    /// Top cell / structure name.
    pub cell: String,
    /// Size of one database unit in meters.
    pub db_unit_meters: f64,
    pub polygons: Vec<LayoutPolygon>,
    pub ports: Vec<Port>,
    /// GDS layer number -> logical layer name (from `*.gds.map`).
    pub layer_names: BTreeMap<i16, String>,
}

impl LayoutIR {
    pub fn bbox(&self) -> Bbox {
        let mut b = Bbox::empty();
        for p in &self.polygons {
            b.extend(&Bbox::from_points(&p.points));
        }
        b
    }

    pub fn polygon_count(&self) -> usize {
        self.polygons.len()
    }

    pub fn vertex_count(&self) -> usize {
        self.polygons.iter().map(|p| p.points.len()).sum()
    }

    /// Number of polygons per GDS layer number.
    pub fn layer_histogram(&self) -> BTreeMap<i16, usize> {
        let mut h = BTreeMap::new();
        for p in &self.polygons {
            *h.entry(p.layer.layer).or_insert(0) += 1;
        }
        h
    }

    /// Number of vertices per GDS layer number.
    pub fn vertex_histogram(&self) -> BTreeMap<i16, usize> {
        let mut h = BTreeMap::new();
        for p in &self.polygons {
            *h.entry(p.layer.layer).or_insert(0) += p.points.len();
        }
        h
    }

    pub fn layer_name(&self, layer: i16) -> Option<&str> {
        self.layer_names.get(&layer).map(|s| s.as_str())
    }

    /// All distinct logical net names referenced by polygons.
    pub fn nets(&self) -> BTreeMap<String, usize> {
        let mut m = BTreeMap::new();
        for p in &self.polygons {
            if let Some(net) = &p.net {
                *m.entry(net.clone()).or_insert(0) += 1;
            }
        }
        m
    }
}

//! Layout connectivity extraction: group conducting polygons into electrical
//! nets.
//!
//! Connectivity rules (v0):
//! - two polygons on the **same physical layer** are connected if their
//!   footprints overlap or touch;
//! - a **via/contact** polygon connects polygons on its `bottom_layer` and
//!   `top_layer` that overlap the via footprint.
//!
//! Overlap uses bounding boxes (exact for rectilinear shapes in a typical
//! database). Net assignment of electrical names (D/S/G) is done by associating
//! ports with the component that contains them.

use openrdson_core::geometry::{polygon_area, Bbox, Point};
use openrdson_core::geomops::{polygon_intersection_rects, polygon_minus_polygon, Rect};
use openrdson_core::layout::{LayoutIR, LayoutPolygon};
use openrdson_core::tech::{LayerMap, TechStack};
use std::collections::{BTreeMap, BTreeSet, HashMap};

/// A connected electrical net: the set of polygon indices and physical layers.
#[derive(Debug, Clone, Default)]
pub struct Component {
    pub polygons: Vec<usize>,
    pub physical_layers: BTreeSet<String>,
}

#[derive(Debug, Clone, Default)]
pub struct Connectivity {
    /// Root component index per layout polygon (`usize::MAX` if non-conducting).
    pub component_of: Vec<usize>,
    pub components: Vec<Component>,
    /// Physical layer of each conducting polygon.
    pub poly_physical: Vec<Option<String>>,
    /// Port name -> component index it lands on.
    pub port_component: BTreeMap<String, usize>,
    /// Port name -> polygon index it lands on.
    pub port_polygon: BTreeMap<String, usize>,
    /// Direct same-layer connections (overlap/touch) between polygon indices.
    pub same_layer_connections: Vec<(usize, usize)>,
}

struct UnionFind {
    parent: Vec<usize>,
    rank: Vec<u8>,
}

impl UnionFind {
    fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
            rank: vec![0; n],
        }
    }
    fn find(&mut self, x: usize) -> usize {
        let mut r = x;
        while self.parent[r] != r {
            r = self.parent[r];
        }
        let mut c = x;
        while self.parent[c] != c {
            let next = self.parent[c];
            self.parent[c] = r;
            c = next;
        }
        r
    }
    fn union(&mut self, a: usize, b: usize) {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra == rb {
            return;
        }
        if self.rank[ra] < self.rank[rb] {
            self.parent[ra] = rb;
        } else if self.rank[ra] > self.rank[rb] {
            self.parent[rb] = ra;
        } else {
            self.parent[rb] = ra;
            self.rank[ra] += 1;
        }
    }
}

fn bbox_overlaps(a: &Bbox, b: &Bbox, eps: f64) -> bool {
    a.min_x - eps <= b.max_x
        && b.min_x - eps <= a.max_x
        && a.min_y - eps <= b.max_y
        && b.min_y - eps <= a.max_y
}

/// Even-odd x-intervals of `ring` on the horizontal line `y`.
fn x_intervals_at(ring: &[Point], y: f64) -> Vec<(f64, f64)> {
    let n = ring.len();
    let mut xs: Vec<f64> = Vec::new();
    for i in 0..n {
        let a = ring[i];
        let b = ring[(i + 1) % n];
        if (a.y <= y && b.y > y) || (b.y <= y && a.y > y) {
            let x = a.x + (y - a.y) / (b.y - a.y) * (b.x - a.x);
            xs.push(x);
        }
    }
    xs.sort_by(|p, q| p.partial_cmp(q).unwrap());
    xs.chunks_exact(2).map(|c| (c[0], c[1])).collect()
}

/// True if two polygons overlap or touch (shared edge/corner counts as
/// connected). Uses exact scanline interval intersection, which is correct for
/// the interdigitated comb shapes where bounding-box overlap would be a false
/// positive.
fn polygons_connected(pa: &[Point], pb: &[Point], eps: f64) -> bool {
    if pa.len() < 3 || pb.len() < 3 {
        return false;
    }
    let ba = Bbox::from_points(pa);
    let bb = Bbox::from_points(pb);
    if !bbox_overlaps(&ba, &bb, eps) {
        return false;
    }
    let mut ys: Vec<f64> = pa.iter().chain(pb.iter()).map(|p| p.y).collect();
    ys.sort_by(|a, b| a.partial_cmp(b).unwrap());
    ys.dedup_by(|a, b| (*a - *b).abs() < 1e-15);
    for k in 0..ys.len().saturating_sub(1) {
        let y = 0.5 * (ys[k] + ys[k + 1]);
        let ia = x_intervals_at(pa, y);
        let ib = x_intervals_at(pb, y);
        for &(a0, a1) in &ia {
            for &(b0, b1) in &ib {
                if a1.min(b1) - a0.max(b0) > -eps {
                    return true;
                }
            }
        }
    }
    false
}

fn point_in_polygon(p: Point, ring: &[Point]) -> bool {
    let n = ring.len();
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let (xi, yi) = (ring[i].x, ring[i].y);
        let (xj, yj) = (ring[j].x, ring[j].y);
        if (yi > p.y) != (yj > p.y) {
            let x_int = xi + (p.y - yi) / (yj - yi) * (xj - xi);
            if p.x < x_int {
                inside = !inside;
            }
        }
        j = i;
    }
    inside
}

fn is_marker(logical: &str) -> bool {
    logical.starts_with("seed_")
        || logical.contains("_prop")
        || logical.contains("marker")
        || logical.starts_with("BLOCK_")
        || logical.contains("drawing")
}

/// Union all same-layer polygons whose footprints overlap or touch.
fn union_same_layer(
    uf: &mut UnionFind,
    items: &[(usize, Bbox)],
    layout: &LayoutIR,
    eps: f64,
    cell: f64,
    connections: &mut Vec<(usize, usize)>,
) {
    // Spatial hash over a uniform grid.
    let key = |x: f64, y: f64| -> (i64, i64) {
        ((x / cell).floor() as i64, (y / cell).floor() as i64)
    };
    let mut grid: HashMap<(i64, i64), Vec<usize>> = HashMap::new();
    for (k, (_, b)) in items.iter().enumerate() {
        let (x0, y0) = key(b.min_x, b.min_y);
        let (x1, y1) = key(b.max_x, b.max_y);
        for gx in x0..=x1 {
            for gy in y0..=y1 {
                grid.entry((gx, gy)).or_default().push(k);
            }
        }
    }
    let mut checked: BTreeSet<(usize, usize)> = BTreeSet::new();
    for (k, (_, b)) in items.iter().enumerate() {
        let (x0, y0) = key(b.min_x, b.min_y);
        let (x1, y1) = key(b.max_x, b.max_y);
        for gx in x0..=x1 {
            for gy in y0..=y1 {
                if let Some(bucket) = grid.get(&(gx, gy)) {
                    for &m in bucket {
                        if m <= k {
                            continue;
                        }
                        let pair = (k, m);
                        if !checked.insert(pair) {
                            continue;
                        }
                        if polygons_connected(
                            &layout.polygons[items[k].0].points,
                            &layout.polygons[items[m].0].points,
                            eps,
                        ) {
                            uf.union(items[k].0, items[m].0);
                            connections.push((items[k].0, items[m].0));
                        }
                    }
                }
            }
        }
    }
}

/// Extract electrical connectivity for the layout.
pub fn extract_connectivity(
    layout: &LayoutIR,
    stack: &TechStack,
    map: &LayerMap,
) -> Connectivity {
    let eps = 1e-12;
    let cell = 4e-6;
    let n = layout.polygons.len();
    let mut uf = UnionFind::new(n);
    let mut poly_physical: Vec<Option<String>> = vec![None; n];
    // Same-layer union is keyed by *logical* layer, because distinct logical
    // layers can share a physical layer (e.g. a source/drain net and a body net
    // that both map to `RX`) yet must not be shorted by overlap.
    let mut by_logical: BTreeMap<String, Vec<(usize, Bbox)>> = BTreeMap::new();
    let mut physical_of_logical: BTreeMap<String, String> = BTreeMap::new();
    let mut logicals_by_physical: BTreeMap<String, Vec<String>> = BTreeMap::new();
    let mut vias: Vec<(usize, String, Bbox)> = Vec::new();

    for (i, poly) in layout.polygons.iter().enumerate() {
        let Some(logical) = layout.layer_name(poly.layer.layer) else {
            continue;
        };
        if is_marker(logical) {
            continue;
        }
        let physical = match map
            .conducting
            .get(logical)
            .or_else(|| map.via.get(logical))
        {
            Some(p) => p.clone(),
            None => continue,
        };
        let bbox = Bbox::from_points(&poly.points);
        if bbox.is_empty() {
            continue;
        }
        poly_physical[i] = Some(physical.clone());
        physical_of_logical.insert(logical.to_string(), physical.clone());
        let e = logicals_by_physical.entry(physical.clone()).or_default();
        if !e.iter().any(|l| l == logical) {
            e.push(logical.to_string());
        }
        if map.via.contains_key(logical) {
            vias.push((i, physical, bbox));
        } else {
            by_logical
                .entry(logical.to_string())
                .or_default()
                .push((i, bbox));
        }
    }

    let mut same_layer_connections: Vec<(usize, usize)> = Vec::new();
    for items in by_logical.values() {
        union_same_layer(&mut uf, items, layout, eps, cell, &mut same_layer_connections);
    }

    // Via connectivity: weld all conducting polygons (of any logical layer that
    // shares the via's physical bottom/top layer) overlapping the via footprint.
    for (vi, physical, _vb) in &vias {
        let Some(via) = stack.via(physical) else {
            continue;
        };
        let (Some(bottom), Some(top)) = (via.bottom_layer.as_deref(), via.top_layer.as_deref())
        else {
            continue;
        };
        let via_pts = &layout.polygons[*vi].points;
        let mut group = vec![*vi];
        for phys in [bottom, top] {
            // A contact lands on the feature drawn on top of its layer. When the
            // layer has nested regions (e.g. source/drain diffusion inside a
            // full-die well, both on RX), connect only to the *smallest-area*
            // overlapping polygon, not the enclosing substrate.
            if let Some(logicals) = logicals_by_physical.get(phys) {
                let mut best: Option<(usize, f64)> = None;
                for logical in logicals {
                    if let Some(items) = by_logical.get(logical) {
                        for (j, _jb) in items {
                            if polygons_connected(via_pts, &layout.polygons[*j].points, eps) {
                                let a = polygon_area(&layout.polygons[*j].points);
                                if best.map(|(_, ba)| a < ba).unwrap_or(true) {
                                    best = Some((*j, a));
                                }
                            }
                        }
                    }
                }
                if let Some((j, _)) = best {
                    group.push(j);
                }
            }
        }
        for k in 1..group.len() {
            uf.union(group[0], group[k]);
        }
    }

    // Collect components.
    let mut root_to_comp: HashMap<usize, usize> = HashMap::new();
    let mut component_of = vec![usize::MAX; n];
    let mut components: Vec<Component> = Vec::new();
    for i in 0..n {
        if poly_physical[i].is_none() {
            continue;
        }
        let r = uf.find(i);
        let ci = *root_to_comp.entry(r).or_insert_with(|| {
            components.push(Component::default());
            components.len() - 1
        });
        component_of[i] = ci;
        components[ci].polygons.push(i);
        if let Some(p) = &poly_physical[i] {
            components[ci].physical_layers.insert(p.clone());
        }
    }

    // Assign ports to the component containing their landing polygon.
    let mut port_component = BTreeMap::new();
    let mut port_polygon = BTreeMap::new();
    for port in &layout.ports {
        let Some(physical) = map.conducting.get(&port.layer) else {
            continue;
        };
        let p = Point::new(port.x, port.y);
        let mut hit = None;
        for (i, poly) in layout.polygons.iter().enumerate() {
            if poly_physical[i].as_deref() != Some(physical.as_str()) {
                continue;
            }
            if point_in_polygon(p, &poly.points) {
                hit = Some(i);
                break;
            }
        }
        if let Some(i) = hit {
            port_polygon.insert(port.name.clone(), i);
            port_component.insert(port.name.clone(), component_of[i]);
        }
    }

    Connectivity {
        component_of,
        components,
        poly_physical,
        port_component,
        port_polygon,
        same_layer_connections,
    }
}

impl Connectivity {
    pub fn component_count(&self) -> usize {
        self.components.len()
    }
}

/// The channel regions of a device: `active ∩ gate`, as rectangles.
pub fn channel_regions(layout: &LayoutIR, gate_logical: &str, active_logical: &str) -> Vec<Rect> {
    let gates: Vec<&Vec<Point>> = layout
        .polygons
        .iter()
        .filter(|p| layout.layer_name(p.layer.layer) == Some(gate_logical))
        .map(|p| &p.points)
        .collect();
    let actives: Vec<&Vec<Point>> = layout
        .polygons
        .iter()
        .filter(|p| layout.layer_name(p.layer.layer) == Some(active_logical))
        .map(|p| &p.points)
        .collect();
    let mut out = Vec::new();
    for a in &actives {
        let ba = Bbox::from_points(a);
        for g in &gates {
            let bg = Bbox::from_points(g);
            if !bbox_overlaps(&ba, &bg, 0.0) {
                continue;
            }
            out.extend(polygon_intersection_rects(a, g));
        }
    }
    out
}

/// Return a copy of the layout with each `active_logical` polygon cut by the
/// union of `gate_logical` polygons (subtracting the channel region).
///
/// This splits a continuous active diffusion into source-side and drain-side
/// pieces, so connectivity no longer shorts them through the channel.
pub fn apply_channel_cut(
    layout: &LayoutIR,
    gate_logical: &str,
    active_logical: &str,
) -> LayoutIR {
    let gates: Vec<Vec<Point>> = layout
        .polygons
        .iter()
        .filter(|p| layout.layer_name(p.layer.layer) == Some(gate_logical))
        .map(|p| p.points.clone())
        .collect();
    if gates.is_empty() {
        return layout.clone();
    }

    let mut polygons = Vec::with_capacity(layout.polygons.len());
    for poly in &layout.polygons {
        if layout.layer_name(poly.layer.layer) == Some(active_logical) {
            // Decompose the active polygon into rectangles, then subtract each
            // gate polygon in turn.
            let mut rects = polygon_minus_polygon(&poly.points, &[]);
            for g in &gates {
                if rects.is_empty() {
                    break;
                }
                let mut next = Vec::new();
                for r in &rects {
                    next.extend(polygon_minus_polygon(&r.to_ring(), g));
                }
                rects = next;
            }
            for r in rects {
                polygons.push(LayoutPolygon {
                    layer: poly.layer,
                    net: poly.net.clone(),
                    points: r.to_ring(),
                    props: poly.props.clone(),
                });
            }
        } else {
            polygons.push(poly.clone());
        }
    }

    LayoutIR {
        cell: layout.cell.clone(),
        db_unit_meters: layout.db_unit_meters,
        polygons,
        ports: layout.ports.clone(),
        layer_names: layout.layer_names.clone(),
    }
}

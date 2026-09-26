//! Adaptive (quadtree) 2.5D sheet meshing.
//!
//! The uniform builder in [`crate::sheet`] grids every polygon rectangle with a
//! fixed `nx × ny` cell. This module keeps the same base grid but stores it as a
//! **quadtree**: each base cell is a level-0 leaf and can be split into four
//! children (level+1). The mesh is rebuilt from the leaves, so refinement is a
//! pure function of the grid state.
//!
//! Conformity uses a 2:1 balance plus a coarse-edge split: a coarse leaf edge
//! whose finer neighbour contributes a midpoint node is split into two half
//! resistors there. The network therefore stays a plain resistor graph and the
//! solver is unchanged.

use crate::sheet::{SheetChannel, SheetNetwork, SheetPolygon, SheetVia};
use openrdson_core::geomops::polygon_minus_polygon;
use std::collections::BTreeMap;

/// Sentinel for "no cell" in `QuadCell::parent` / `children`.
pub const NO_CELL: u32 = u32::MAX;

/// Adaptive-refinement controls.
#[derive(Debug, Clone)]
pub struct AdaptiveConfig {
    pub enabled: bool,
    /// Relative indicator convergence: stop when the total indicator changes by
    /// less than this between iterations.
    pub tol: f64,
    pub max_level: u8,
    /// Stop when the node count exceeds this.
    pub max_nodes: usize,
    /// Dörfler marking fraction: refine the smallest set of cells whose
    /// indicator sums to at least this fraction of the total.
    pub mark_frac: f64,
    /// Maximum refine→solve iterations.
    pub iters: usize,
    /// A-priori refinement radius around vias/terminals (m); 0 disables.
    pub refine_near_vias_m: f64,
}

impl Default for AdaptiveConfig {
    fn default() -> Self {
        Self {
            enabled: false,
            tol: 1e-3,
            max_level: 3,
            max_nodes: 400_000,
            mark_frac: 0.5,
            iters: 8,
            refine_near_vias_m: 0.0,
        }
    }
}

/// Per-cell power (W) implied by a node potential, used as the refinement
/// indicator: `(|∇u|²)·dx·dy / r_sheet`.
pub fn cell_power(net: &SheetNetwork, u: &[f64]) -> Vec<f64> {
    let mut out = Vec::with_capacity(net.cells.len());
    for (i, cell) in net.cells.iter().enumerate() {
        let (gx, gy, dx, dy) = cell_grad(net, cell, u);
        let rs = net.cell_r_sheet.get(i).copied().unwrap_or(1.0).max(1e-30);
        out.push((gx * gx + gy * gy) * dx * dy / rs);
    }
    out
}

/// Per-cell gradient `(∂u/∂x, ∂u/∂y)` and extents `(dx, dy)`.
fn cell_grad(net: &SheetNetwork, cell: &[u32; 4], u: &[f64]) -> (f64, f64, f64, f64) {
    let (u0, u1, u2, u3) = (
        u[cell[0] as usize],
        u[cell[1] as usize],
        u[cell[2] as usize],
        u[cell[3] as usize],
    );
    let dx = (net.x[cell[1] as usize] - net.x[cell[0] as usize]).abs().max(1e-30);
    let dy = (net.y[cell[3] as usize] - net.y[cell[0] as usize]).abs().max(1e-30);
    (
        ((u1 + u2) - (u0 + u3)) / (2.0 * dx),
        ((u2 + u3) - (u0 + u1)) / (2.0 * dy),
        dx,
        dy,
    )
}

/// ZZ-style gradient-jump error indicator: for each cell, the edge-length
/// weighted jump of the current-density gradient against its edge neighbours.
/// Unlike the raw power, this is a genuine error estimate and localises the
/// refinement at high-curvature features (necks, corners, contacts).
pub fn gradient_jump_indicator(net: &SheetNetwork, u: &[f64]) -> Vec<f64> {
    let n = net.cells.len();
    let mut gx = vec![0.0f64; n];
    let mut gy = vec![0.0f64; n];
    for (i, cell) in net.cells.iter().enumerate() {
        let (a, b, _, _) = cell_grad(net, cell, u);
        gx[i] = a;
        gy[i] = b;
    }
    // Map each (unordered) edge node pair to the cells that own it.
    let mut edge_cells: BTreeMap<(u32, u32), Vec<usize>> = BTreeMap::new();
    for (i, cell) in net.cells.iter().enumerate() {
        for k in 0..4 {
            let (a, b) = (cell[k], cell[(k + 1) % 4]);
            edge_cells
                .entry((a.min(b), a.max(b)))
                .or_default()
                .push(i);
        }
    }
    let mut out = vec![0.0f64; n];
    for ((a, b), cells) in &edge_cells {
        if cells.len() != 2 {
            continue;
        }
        let (i, j) = (cells[0], cells[1]);
        let dgx = gx[i] - gx[j];
        let dgy = gy[i] - gy[j];
        let jump = (dgx * dgx + dgy * dgy).sqrt();
        let len = ((net.x[*b as usize] - net.x[*a as usize]).powi(2)
            + (net.y[*b as usize] - net.y[*a as usize]).powi(2))
        .sqrt();
        out[i] += jump * len;
        out[j] += jump * len;
    }
    out
}



/// A quadtree cell: an axis-aligned rectangle at a refinement level.
#[derive(Debug, Clone)]
pub struct QuadCell {
    pub x0: f64,
    pub y0: f64,
    pub x1: f64,
    pub y1: f64,
    pub level: u8,
    pub parent: u32,
    pub children: [u32; 4],
    /// Index of the source polygon (for `r_sheet` and component lookup).
    pub poly: u32,
}

impl QuadCell {
    pub fn is_leaf(&self) -> bool {
        self.children[0] == NO_CELL
    }
    pub fn w(&self) -> f64 {
        self.x1 - self.x0
    }
    pub fn h(&self) -> f64 {
        self.y1 - self.y0
    }
    pub fn cx(&self) -> f64 {
        0.5 * (self.x0 + self.x1)
    }
    pub fn cy(&self) -> f64 {
        0.5 * (self.y0 + self.y1)
    }
}

/// Quadtree over the polygon rectangles.
#[derive(Debug, Clone, Default)]
pub struct AdaptiveGrid {
    pub cells: Vec<QuadCell>,
    pub roots: Vec<u32>,
    /// Nominal base cell size (m); leaf size is `base / 2^level`.
    pub base_cell: f64,
    pub max_level: u8,
}

impl AdaptiveGrid {
    /// Build the level-0 grid for `polygons` using the same cell sizing as the
    /// uniform builder (`cell_target`, `cap`, `per_layer_cell`).
    pub fn from_polygons(
        polygons: &[SheetPolygon],
        cell_target: f64,
        cap: usize,
        per_layer_cell: &BTreeMap<String, f64>,
    ) -> Self {
        let mut grid = AdaptiveGrid {
            base_cell: cell_target,
            max_level: 0,
            ..Default::default()
        };
        for (pi, poly) in polygons.iter().enumerate() {
            let cell = per_layer_cell
                .get(&poly.physical)
                .copied()
                .filter(|c| *c > 0.0)
                .unwrap_or(cell_target);
            for rect in polygon_minus_polygon(&poly.points, &[]) {
                let w = rect.x1 - rect.x0;
                let h = rect.y1 - rect.y0;
                if w <= 0.0 || h <= 0.0 {
                    continue;
                }
                let nx = ((w / cell).ceil() as usize).clamp(1, cap);
                let ny = ((h / cell).ceil() as usize).clamp(1, cap);
                let dx = w / nx as f64;
                let dy = h / ny as f64;
                for j in 0..ny {
                    for i in 0..nx {
                        let id = grid.cells.len() as u32;
                        grid.cells.push(QuadCell {
                            x0: rect.x0 + dx * i as f64,
                            y0: rect.y0 + dy * j as f64,
                            x1: rect.x0 + dx * (i + 1) as f64,
                            y1: rect.y0 + dy * (j + 1) as f64,
                            level: 0,
                            parent: NO_CELL,
                            children: [NO_CELL; 4],
                            poly: pi as u32,
                        });
                        grid.roots.push(id);
                    }
                }
            }
        }
        grid
    }

    /// Leaf cell indices, in deterministic arena order.
    pub fn leaves(&self) -> Vec<u32> {
        (0..self.cells.len() as u32)
            .filter(|&c| self.cells[c as usize].is_leaf())
            .collect()
    }

    /// Leaves in the same order [`build_network`] emits mesh cells (grouped by
    /// polygon, then arena order), so a per-cell indicator maps back to a leaf.
    pub fn leaf_order(&self) -> Vec<u32> {
        let mut by_poly: BTreeMap<u32, Vec<u32>> = BTreeMap::new();
        for leaf in self.leaves() {
            by_poly
                .entry(self.cells[leaf as usize].poly)
                .or_default()
                .push(leaf);
        }
        let mut out = Vec::new();
        for v in by_poly.values() {
            out.extend_from_slice(v);
        }
        out
    }

    pub fn leaf_count(&self) -> usize {
        self.cells.iter().filter(|c| c.is_leaf()).count()
    }

    /// Split a leaf into four children. Returns false if it is not a leaf or is
    /// already at `max_level`.
    pub fn refine(&mut self, c: u32) -> bool {
        let cell = self.cells[c as usize].clone();
        if !cell.is_leaf() || cell.level >= self.max_level {
            return false;
        }
        let (mx, my) = (cell.cx(), cell.cy());
        let (x0, y0, x1, y1) = (cell.x0, cell.y0, cell.x1, cell.y1);
        // Children in fixed order: SW, SE, NW, NE.
        let quads = [(x0, y0, mx, my), (mx, y0, x1, my), (x0, my, mx, y1), (mx, my, x1, y1)];
        let mut ids = [NO_CELL; 4];
        for (k, (qx0, qy0, qx1, qy1)) in quads.iter().enumerate() {
            let id = self.cells.len() as u32;
            self.cells.push(QuadCell {
                x0: *qx0,
                y0: *qy0,
                x1: *qx1,
                y1: *qy1,
                level: cell.level + 1,
                parent: c,
                children: [NO_CELL; 4],
                poly: cell.poly,
            });
            ids[k] = id;
        }
        self.cells[c as usize].children = ids;
        true
    }

    /// Enforce a 2:1 balance: no leaf may be more than one level finer than an
    /// edge-adjacent leaf. Repeats until stable.
    pub fn balance(&mut self) {
        loop {
            let mut changed = false;
            let leaves = self.leaves();
            for &a in &leaves {
                if !self.cells[a as usize].is_leaf() {
                    continue;
                }
                let ca = self.cells[a as usize].clone();
                for &b in &leaves {
                    if a == b || !self.cells[b as usize].is_leaf() {
                        continue;
                    }
                    let cb = self.cells[b as usize].clone();
                    if ca.level <= cb.level + 1 {
                        continue;
                    }
                    if shares_edge(&ca, &cb) {
                        if self.refine(b) {
                            changed = true;
                        }
                    }
                }
            }
            if !changed {
                break;
            }
        }
    }

    /// Refine coarse leaves until every edge-adjacent pair sits at the same
    /// level, making the mesh **conforming with no hanging nodes**. Then the
    /// node-centred finite-volume weighting (`apply_fv_weights`) applies
    /// verbatim (no coarse-edge split needed).
    pub fn conform(&mut self) {
        loop {
            let mut changed = false;
            let leaves = self.leaves();
            for &a in &leaves {
                if !self.cells[a as usize].is_leaf() {
                    continue;
                }
                let ca = self.cells[a as usize].clone();
                for &b in &leaves {
                    if a == b || !self.cells[b as usize].is_leaf() {
                        continue;
                    }
                    let cb = self.cells[b as usize].clone();
                    if ca.level < cb.level && shares_edge(&ca, &cb) {
                        if self.refine(a) {
                            changed = true;
                        }
                        break;
                    }
                }
            }
            if !changed {
                break;
            }
        }
    }

    /// Verify the mesh is conforming in the graph sense: for every coarse leaf
    /// edge whose midpoint node exists (a hanging node contributed by a finer
    /// neighbour), the network must contain the two half-edges through it.
    pub fn check_network_conforming(&self, net: &SheetNetwork) -> Result<(), String> {
        let key = |x: f64, y: f64| ((x / 1e-12).round() as i64, (y / 1e-12).round() as i64);
        let mut node_at: BTreeMap<(i64, i64), u32> = BTreeMap::new();
        for (n, (&x, &y)) in net.x.iter().zip(net.y.iter()).enumerate() {
            node_at.insert(key(x, y), n as u32);
        }
        let has_edge = |a: u32, b: u32| -> bool {
            net.edges
                .iter()
                .any(|&(x, y, _)| (x == a && y == b) || (x == b && y == a))
        };
        for leaf in self.leaves() {
            let c = &self.cells[leaf as usize];
            let corners = [
                ((c.x0, c.y0), (c.x1, c.y0), (c.cx(), c.y0)),
                ((c.x1, c.y0), (c.x1, c.y1), (c.x1, c.cy())),
                ((c.x1, c.y1), (c.x0, c.y1), (c.cx(), c.y1)),
                ((c.x0, c.y1), (c.x0, c.y0), (c.x0, c.cy())),
            ];
            for ((ax, ay), (bx, by), (mx, my)) in corners {
                let (Some(&a), Some(&b)) = (node_at.get(&key(ax, ay)), node_at.get(&key(bx, by)))
                else {
                    continue;
                };
                if let Some(&m) = node_at.get(&key(mx, my)) {
                    if m != a && m != b && (!has_edge(a, m) || !has_edge(m, b)) {
                        return Err(format!(
                            "hanging node {m} on leaf edge ({a},{b}) is not split"
                        ));
                    }
                }
            }
        }
        Ok(())
    }
}

/// Whether two axis-aligned rectangles share a positive-length edge.
fn shares_edge(a: &QuadCell, b: &QuadCell) -> bool {
    let x_overlap = a.x0.min(a.x1) < b.x1 && b.x0 < a.x1;
    let y_overlap = a.y0.min(a.y1) < b.y1 && b.y0 < a.y1;
    let touch_v = (a.x1 - b.x0).abs() < 1e-15 || (b.x1 - a.x0).abs() < 1e-15;
    let touch_h = (a.y1 - b.y0).abs() < 1e-15 || (b.y1 - a.y0).abs() < 1e-15;
    (touch_v && y_overlap) || (touch_h && x_overlap)
}

/// Build a [`SheetNetwork`] from the leaves of `grid`.
///
/// `connections`, `vias` and `channels` are handled exactly as in the uniform
/// builder. Coarse edges adjacent to a finer leaf are split at the midpoint
/// (hanging node) so the mesh is conforming.
pub fn build_network(
    grid: &AdaptiveGrid,
    polygons: &[SheetPolygon],
    connections: &[(usize, usize)],
    vias: &[SheetVia],
    channels: &[SheetChannel],
) -> SheetNetwork {
    let mut net = SheetNetwork::default();
    let mut key_map: BTreeMap<(String, usize, i64, i64), u32> = BTreeMap::new();
    let q = |v: f64| (v / 1e-12).round() as i64;

    let add_node = |net: &mut SheetNetwork,
                        key_map: &mut BTreeMap<(String, usize, i64, i64), u32>,
                        physical: &str,
                        component: usize,
                        x: f64,
                        y: f64|
     -> u32 {
        let key = (physical.to_string(), component, q(x), q(y));
        if let Some(&n) = key_map.get(&key) {
            if net.component[n as usize] == usize::MAX {
                net.component[n as usize] = component;
            }
            return n;
        }
        let n = net.x.len() as u32;
        net.physical.push(physical.to_string());
        net.x.push(x);
        net.y.push(y);
        net.component.push(component);
        net.physical_nodes
            .entry(physical.to_string())
            .or_default()
            .push(n);
        key_map.insert(key, n);
        n
    };

    // Group leaves by polygon so `poly_nodes` keeps the uniform semantics.
    let mut by_poly: BTreeMap<u32, Vec<u32>> = BTreeMap::new();
    for leaf in grid.leaves() {
        by_poly
            .entry(grid.cells[leaf as usize].poly)
            .or_default()
            .push(leaf);
    }

    for (pi, poly) in polygons.iter().enumerate() {
        let leaves = by_poly.get(&(pi as u32)).cloned().unwrap_or_default();
        let mut pnodes: Vec<u32> = Vec::new();
        // Nodes first (all leaf corners), so midpoint lookups work.
        for &lf in &leaves {
            let c = &grid.cells[lf as usize];
            for &(x, y) in &[
                (c.x0, c.y0),
                (c.x1, c.y0),
                (c.x1, c.y1),
                (c.x0, c.y1),
            ] {
                let n = add_node(&mut net, &mut key_map, &poly.physical, poly.component, x, y);
                pnodes.push(n);
            }
        }
        net.poly_nodes.push(pnodes);

        // Edges: 4 per leaf, split at a hanging midpoint when present.
        let node_at = |key_map: &BTreeMap<(String, usize, i64, i64), u32>,
                       x: f64,
                       y: f64|
         -> Option<u32> {
            key_map
                .get(&(poly.physical.clone(), poly.component, q(x), q(y)))
                .copied()
        };
        let edge = |net: &mut SheetNetwork,
                        seen: &mut std::collections::BTreeSet<(u32, u32)>,
                        a: u32,
                        b: u32,
                        r: f64,
                        mid: (f64, f64)| {
            if r <= 0.0 {
                return;
            }
            match node_at(&key_map, mid.0, mid.1) {
                Some(m) if m != a && m != b => {
                    emit_edge(net, seen, a, m, 1.0 / (r * 0.5));
                    emit_edge(net, seen, m, b, 1.0 / (r * 0.5));
                }
                _ => emit_edge(net, seen, a, b, 1.0 / r),
            }
        };
        let mut seen: std::collections::BTreeSet<(u32, u32)> = std::collections::BTreeSet::new();
        for &lf in &leaves {
            let c = &grid.cells[lf as usize];
            let (n00, n10, n11, n01) = (
                node_at(&key_map, c.x0, c.y0).unwrap(),
                node_at(&key_map, c.x1, c.y0).unwrap(),
                node_at(&key_map, c.x1, c.y1).unwrap(),
                node_at(&key_map, c.x0, c.y1).unwrap(),
            );
            let rh = poly.r_sheet * c.w() / c.h();
            let rv = poly.r_sheet * c.h() / c.w();
            edge(&mut net, &mut seen, n00, n10, rh, (c.cx(), c.y0));
            edge(&mut net, &mut seen, n01, n11, rh, (c.cx(), c.y1));
            edge(&mut net, &mut seen, n00, n01, rv, (c.x0, c.cy()));
            edge(&mut net, &mut seen, n10, n11, rv, (c.x1, c.cy()));
            net.cells.push([n00, n10, n11, n01]);
            net.cell_r_sheet.push(poly.r_sheet);
            net.cell_level.push(c.level);
            net.cell_parent.push(c.parent);
        }
    }

    // Same-layer overlap/touch connections (nearest node pair), as uniform.
    for &(a, b) in connections {
        if a >= net.poly_nodes.len() || b >= net.poly_nodes.len() {
            continue;
        }
        let (pa, pb) = (net.poly_nodes[a].clone(), net.poly_nodes[b].clone());
        if pa.is_empty() || pb.is_empty() {
            continue;
        }
        let (ax, ay) = (net.x[pa[0] as usize], net.y[pa[0] as usize]);
        if let Some(na) = nearest(&pa, &net.x, &net.y, ax, ay) {
            let (nx, ny) = (net.x[na as usize], net.y[na as usize]);
            if let Some(nb) = nearest(&pb, &net.x, &net.y, nx, ny) {
                net.edges.push((na, nb, 1.0));
            }
        }
    }

    // Vias: connect same-net nodes on the bottom and top layers. A point via
    // (zero footprint) connects the single nearest node on each layer; a via
    // with a footprint spreads its current over every node in the rect
    // [x±half_width] x [y±half_height] (the via is a uniform vertical resistor,
    // so the total conductance 1/R is split evenly across the node pairs).
    for via in vias {
        if !(via.resistance > 0.0) {
            continue;
        }
        let footprint = |physical: &str| -> Vec<u32> {
            let Some(nodes) = net.physical_nodes.get(physical) else {
                return Vec::new();
            };
            let cand: Vec<u32> = nodes
                .iter()
                .copied()
                .filter(|&n| net.component[n as usize] == via.component)
                .collect();
            if via.half_width > 0.0 && via.half_height > 0.0 {
                let inside: Vec<u32> = cand
                    .iter()
                    .copied()
                    .filter(|&n| {
                        (net.x[n as usize] - via.center.x).abs() <= via.half_width
                            && (net.y[n as usize] - via.center.y).abs() <= via.half_height
                    })
                    .collect();
                if !inside.is_empty() {
                    return inside;
                }
                // Footprint smaller than the mesh cell: fall back to the
                // nearest node so a sub-cell via stays connected.
            }
            nearest(&cand, &net.x, &net.y, via.center.x, via.center.y)
                .into_iter()
                .collect()
        };
        let b = footprint(&via.bottom);
        let t = footprint(&via.top);
        if b.is_empty() || t.is_empty() {
            continue;
        }
        let g_each = 1.0 / via.resistance / (b.len() * t.len()) as f64;
        for &bn in &b {
            for &tn in &t {
                net.edges.push((bn, tn, g_each));
            }
        }
    }

    // Channels.
    for ch in channels {
        if let Some((d, s)) =
            net.channel_nodes(ch.drain, ch.source, ch.drain_component, ch.source_component)
        {
            if ch.resistance > 0.0 && ch.resistance.is_finite() {
                net.edges.push((d, s, 1.0 / ch.resistance));
            }
        }
    }

    crate::sheet::apply_fv_weights(&mut net);
    net.merge_nets_by_component();
    net
}

/// Push an edge unless the (unordered) node pair was already emitted, so each
/// shared mesh edge contributes its conductance exactly once (matching the
/// uniform builder, which adds each interior edge once).
fn emit_edge(
    net: &mut SheetNetwork,
    seen: &mut std::collections::BTreeSet<(u32, u32)>,
    a: u32,
    b: u32,
    g: f64,
) {
    let key = (a.min(b), a.max(b));
    if seen.insert(key) {
        net.edges.push((a, b, g));
    }
}

fn nearest(nodes: &[u32], x: &[f64], y: &[f64], px: f64, py: f64) -> Option<u32> {
    let mut best = None;
    let mut bd = f64::INFINITY;
    for &n in nodes {
        let d = (x[n as usize] - px).powi(2) + (y[n as usize] - py).powi(2);
        if d < bd {
            bd = d;
            best = Some(n);
        }
    }
    best
}

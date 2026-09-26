//! 2.5D sheet-resistance network extraction.
//!
//! Full-die 3D FEM is infeasible (0.16 µm features across a 362×1075 µm die).
//! This module extracts a **2D sheet-resistance network per conducting layer**,
//! with lumped via/contact resistors between layers and bias-dependent channel
//! resistors between source/drain. It scales to the whole array.
//!
//! Each polygon is decomposed into rectangles (exact for rectilinear shapes) and
//! gridded adaptively (cell size = min(target, rectangle extent), capped per
//! dimension). Coincident grid nodes are welded, so adjacent rectangles/polygons
//! connect automatically.

use openrdson_core::geomops::polygon_minus_polygon;
use openrdson_core::geometry::Point;
use openrdson_solver::{solve_dirichlet, solve_dirichlet_warm, Sparse};
use std::collections::{BTreeMap, HashMap};

/// A conducting polygon to grid (physical layer, outline, net component).
#[derive(Debug, Clone)]
pub struct SheetPolygon {
    pub physical: String,
    pub points: Vec<Point>,
    pub component: usize,
    pub r_sheet: f64,
}

/// A via/contact connecting two physical layers. A point via (zero footprint)
/// connects the single nearest node on each layer; a via with a nonzero
/// footprint spreads its current over every node in the rect
/// `[x±half_width] × [y±half_height]`.
#[derive(Debug, Clone)]
pub struct SheetVia {
    pub bottom: String,
    pub top: String,
    pub center: Point,
    pub resistance: f64,
    /// Net component of the via, so it only connects same-net nodes.
    pub component: usize,
    /// Half-extent of the via footprint. Zero means a point via.
    pub half_width: f64,
    pub half_height: f64,
}

/// A bias-dependent channel connection between a drain and a source point.
#[derive(Debug, Clone)]
pub struct SheetChannel {
    pub drain: Point,
    pub source: Point,
    pub drain_component: usize,
    pub source_component: usize,
    pub resistance: f64,
}

/// Group nearby vias/contacts (same bottom/top layer and net, within `radius`)
/// into a single equivalent via: conductances add (parallel) and the position is
/// the centroid. This sparsifies dense via arrays without changing the
/// equivalent resistance.
pub fn group_vias(vias: &[SheetVia], radius: f64) -> Vec<SheetVia> {
    // BTreeMap (not HashMap): deterministic grouping order, so re-runs are
    // reproducible.
    use std::collections::BTreeMap;
    let cell = radius.max(1e-12);
    let q = |v: f64| (v / cell).floor() as i64;
    // (conductance, x0, x1, y0, y1, count): the footprint is the bbox of the
    // grouped vias' footprints and the resistance is the parallel combination.
    let mut groups: BTreeMap<
        (String, String, usize, i64, i64),
        (f64, f64, f64, f64, f64, usize),
    > = BTreeMap::new();
    for v in vias {
        if !(v.resistance > 0.0) {
            continue;
        }
        let key = (
            v.bottom.clone(),
            v.top.clone(),
            v.component,
            q(v.center.x),
            q(v.center.y),
        );
        let e = groups.entry(key).or_insert((
            0.0,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
            0,
        ));
        e.0 += 1.0 / v.resistance;
        e.1 = e.1.min(v.center.x - v.half_width);
        e.2 = e.2.max(v.center.x + v.half_width);
        e.3 = e.3.min(v.center.y - v.half_height);
        e.4 = e.4.max(v.center.y + v.half_height);
        e.5 += 1;
    }
    groups
        .into_iter()
        .map(|((bottom, top, component, _, _), (g, x0, x1, y0, y1, _))| SheetVia {
            bottom,
            top,
            component,
            center: Point::new((x0 + x1) * 0.5, (y0 + y1) * 0.5),
            resistance: if g > 0.0 { 1.0 / g } else { f64::INFINITY },
            half_width: (x1 - x0) * 0.5,
            half_height: (y1 - y0) * 0.5,
        })
        .collect()
}

/// Map a solution `old_u` (indexed by `old`'s nodes) onto `new`'s node
/// numbering. Untouched nodes are matched by exact position (physical layer +
/// component + quantized (x, y)); newly refined nodes (absent from `old`) get
/// the nearest old node's value via a coarse spatial grid. Used to warm-start
/// the linear solve between adaptive refinement iterations.
pub fn map_solution(
    old: &SheetNetwork,
    old_u: &[f64],
    new: &SheetNetwork,
    cell: f64,
) -> Vec<f64> {
    let q = |v: f64| (v / 1e-12).round() as i64;
    let mut exact: BTreeMap<(String, usize, i64, i64), f64> = BTreeMap::new();
    for i in 0..old.x.len() {
        exact.insert(
            (
                old.physical[i].clone(),
                old.component[i],
                q(old.x[i]),
                q(old.y[i]),
            ),
            old_u[i],
        );
    }
    let mut out = vec![0.0f64; new.x.len()];
    let mut missing: Vec<usize> = Vec::new();
    for i in 0..new.x.len() {
        let key = (
            new.physical[i].clone(),
            new.component[i],
            q(new.x[i]),
            q(new.y[i]),
        );
        match exact.get(&key) {
            Some(&v) => out[i] = v,
            None => missing.push(i),
        }
    }
    if !missing.is_empty() {
        let c = cell.max(1e-9);
        let gq = |v: f64| (v / c).floor() as i64;
        let mut buckets: BTreeMap<(String, usize, i64, i64), Vec<(f64, f64, f64)>> =
            BTreeMap::new();
        for i in 0..old.x.len() {
            buckets
                .entry((
                    old.physical[i].clone(),
                    old.component[i],
                    gq(old.x[i]),
                    gq(old.y[i]),
                ))
                .or_default()
                .push((old.x[i], old.y[i], old_u[i]));
        }
        for &i in &missing {
            let (ph, comp) = (new.physical[i].clone(), new.component[i]);
            let (cx, cy) = (gq(new.x[i]), gq(new.y[i]));
            let mut best: Option<(f64, f64)> = None;
            for dx in -1i64..=1 {
                for dy in -1i64..=1 {
                    if let Some(pts) = buckets.get(&(ph.clone(), comp, cx + dx, cy + dy)) {
                        for &(ox, oy, ov) in pts {
                            let d2 = (ox - new.x[i]).powi(2) + (oy - new.y[i]).powi(2);
                            if best.map(|(bd, _)| d2 < bd).unwrap_or(true) {
                                best = Some((d2, ov));
                            }
                        }
                    }
                }
            }
            if let Some((_, ov)) = best {
                out[i] = ov;
            }
        }
    }
    out
}

#[derive(Debug, Clone, Default)]
pub struct SheetNetwork {
    pub physical: Vec<String>,
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub component: Vec<usize>,
    /// `(node_a, node_b, conductance)`.
    pub edges: Vec<(u32, u32, f64)>,
    /// Grid quad cells as `[n00, n10, n11, n01]` node indices (CCW).
    pub cells: Vec<[u32; 4]>,
    /// Per-cell sheet resistance (Ω/sq), parallel to `cells` (the polygon's
    /// width-dependent value, so field/report output matches the solver).
    pub cell_r_sheet: Vec<f64>,
    /// Per-cell quadtree refinement level (0 = base grid), parallel to `cells`.
    pub cell_level: Vec<u8>,
    /// Per-cell quadtree parent index (`u32::MAX` for a base cell).
    pub cell_parent: Vec<u32>,
    pub poly_nodes: Vec<Vec<u32>>,
    pub physical_nodes: HashMap<String, Vec<u32>>,
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

/// Re-weight every edge by the node-centred finite-volume control volume.
///
/// The 5-point stencil uses the full cell height as the edge face, which
/// over-counts the domain-boundary rows (a full-height bar reads
/// `(r_sheet*L/W)*ny/(ny+1)`). Here each node owns a control volume whose
/// perpendicular extent is `dy` in the interior and `dy/2` on the top/bottom
/// boundary (likewise `dx`/`dx/2` left/right), and each edge's conductance is
/// `face / (r_sheet * distance)` using the averaged control-volume face.
pub fn apply_fv_weights(net: &mut SheetNetwork) {
    use std::collections::BTreeMap;
    let n = net.x.len();
    if n == 0 {
        return;
    }
    let qx = |x: f64| (x / 1e-12).round() as i64;
    let qy = |y: f64| (y / 1e-12).round() as i64;
    // Nodes grouped by (physical layer, column) and (physical layer, row):
    // each layer's sheet is meshed independently, so a node's CV neighbours are
    // on the *same* layer (overlapping layers must not pollute each other).
    let mut cols: BTreeMap<(String, i64), Vec<usize>> = BTreeMap::new();
    let mut rows: BTreeMap<(String, i64), Vec<usize>> = BTreeMap::new();
    for i in 0..n {
        cols.entry((net.physical[i].clone(), qx(net.x[i]))).or_default().push(i);
        rows.entry((net.physical[i].clone(), qy(net.y[i]))).or_default().push(i);
    }
    // Control-volume half-extents: half the gap to the next/prev node along the
    // same column (y) and row (x). Zero on the domain boundary.
    let mut cv_up = vec![0.0f64; n];
    let mut cv_dn = vec![0.0f64; n];
    let mut cv_rt = vec![0.0f64; n];
    let mut cv_lf = vec![0.0f64; n];
    for nodes in cols.values() {
        let mut v = nodes.clone();
        v.sort_by(|&a, &b| net.y[a].partial_cmp(&net.y[b]).unwrap());
        for k in 0..v.len() {
            if k + 1 < v.len() {
                cv_up[v[k]] = (net.y[v[k + 1]] - net.y[v[k]]) * 0.5;
            }
            if k > 0 {
                cv_dn[v[k]] = (net.y[v[k]] - net.y[v[k - 1]]) * 0.5;
            }
        }
    }
    for nodes in rows.values() {
        let mut v = nodes.clone();
        v.sort_by(|&a, &b| net.x[a].partial_cmp(&net.x[b]).unwrap());
        for k in 0..v.len() {
            if k + 1 < v.len() {
                cv_rt[v[k]] = (net.x[v[k + 1]] - net.x[v[k]]) * 0.5;
            }
            if k > 0 {
                cv_lf[v[k]] = (net.x[v[k]] - net.x[v[k - 1]]) * 0.5;
            }
        }
    }
    // Fix hanging-node control volumes: a node at the midpoint of a coarse
    // cell's edge has no neighbour in the coarse direction (so the node-gap CV
    // above records 0 there), but its control volume still extends half the
    // coarse cell's extent into that cell. Detect such nodes and fill the
    // missing extent.
    let mut node_at: BTreeMap<(String, i64, i64), u32> = BTreeMap::new();
    for i in 0..n {
        node_at.insert((net.physical[i].clone(), qx(net.x[i]), qy(net.y[i])), i as u32);
    }
    for (ci, cell) in net.cells.iter().enumerate() {
        let cdx = (net.x[cell[1] as usize] - net.x[cell[0] as usize]).abs();
        let cdy = (net.y[cell[3] as usize] - net.y[cell[0] as usize]).abs();
        let cx = (net.x[cell[0] as usize] + net.x[cell[1] as usize]) / 2.0;
        let cy = (net.y[cell[0] as usize] + net.y[cell[3] as usize]) / 2.0;
        for k in 0..4 {
            let (a, b) = (cell[k], cell[(k + 1) % 4]);
            let mx = (net.x[a as usize] + net.x[b as usize]) / 2.0;
            let my = (net.y[a as usize] + net.y[b as usize]) / 2.0;
            let phys = net.physical[cell[0] as usize].clone();
            let Some(&m) = node_at.get(&(phys, qx(mx), qy(my))) else { continue };
            if m == a || m == b {
                continue;
            }
            // m is a hanging node on this cell's edge (a, b).
            let horizontal = (net.y[a as usize] - net.y[b as usize]).abs() < 1e-15;
            if horizontal {
                // Edge along x: the coarse cell is above or below m.
                if cy < my {
                    cv_dn[m as usize] = cv_dn[m as usize].max(cdy * 0.5);
                } else {
                    cv_up[m as usize] = cv_up[m as usize].max(cdy * 0.5);
                }
            } else if cx < mx {
                // Edge along y: the coarse cell is left of m.
                cv_lf[m as usize] = cv_lf[m as usize].max(cdx * 0.5);
            } else {
                cv_rt[m as usize] = cv_rt[m as usize].max(cdx * 0.5);
            }
        }
        let _ = ci;
    }

    // Per-edge r_sheet from the owning cell.
    let mut edge_rs: BTreeMap<(u32, u32), f64> = BTreeMap::new();
    for (ci, cell) in net.cells.iter().enumerate() {
        let rs = net.cell_r_sheet.get(ci).copied().unwrap_or(1.0);
        for k in 0..4 {
            let (a, b) = (cell[k], cell[(k + 1) % 4]);
            edge_rs.entry((a.min(b), a.max(b))).or_insert(rs);
        }
    }
    for e in &mut net.edges {
        let (a, b) = (e.0 as usize, e.1 as usize);
        let key = (e.0.min(e.1), e.0.max(e.1));
        let Some(&rs) = edge_rs.get(&key) else { continue };
        let horizontal = (net.y[a] - net.y[b]).abs() < 1e-15;
        let dist = if horizontal {
            (net.x[a] - net.x[b]).abs()
        } else {
            (net.y[a] - net.y[b]).abs()
        };
        if dist <= 0.0 {
            continue;
        }
        // Shared control-volume face: the overlap of the two nodes' CVs across
        // the edge (handles coarse/fine junctions correctly).
        let face = if horizontal {
            cv_up[a].min(cv_up[b]) + cv_dn[a].min(cv_dn[b])
        } else {
            cv_rt[a].min(cv_rt[b]) + cv_lf[a].min(cv_lf[b])
        };
        if face <= 0.0 {
            continue;
        }
        e.2 = face / (rs * dist).max(1e-300);
    }
}

impl SheetNetwork {
    /// Build the network. `cell_target` is the maximum cell size (m); large
    /// polygons are capped to `cap` cells per dimension.
    pub fn build(
        polygons: &[SheetPolygon],
        same_layer_connections: &[(usize, usize)],
        vias: &[SheetVia],
        channels: &[SheetChannel],
        cell_target: f64,
        cap: usize,
    ) -> Self {
        Self::build_with_cells(
            polygons,
            same_layer_connections,
            vias,
            channels,
            cell_target,
            cap,
            &std::collections::BTreeMap::new(),
        )
    }

    /// Build with a default cell size plus optional per-physical-layer overrides
    /// (R3D-style layer-based refinement).
    pub fn build_with_cells(
        polygons: &[SheetPolygon],
        same_layer_connections: &[(usize, usize)],
        vias: &[SheetVia],
        channels: &[SheetChannel],
        cell_target: f64,
        cap: usize,
        per_layer_cell: &std::collections::BTreeMap<String, f64>,
    ) -> Self {
        let mut net = SheetNetwork::default();
        // Nodes are welded by (physical, component, x, y): polygons of the same
        // net share nodes, but distinct nets never do (so touching D/S polygons
        // cannot short).
        let mut key_map: HashMap<(String, usize, i64, i64), u32> = HashMap::new();
        let q = |v: f64| (v / 1e-12).round() as i64;

        let add_node = |net: &mut SheetNetwork,
                            key_map: &mut HashMap<(String, usize, i64, i64), u32>,
                            physical: &str,
                            component: usize,
                            x: f64,
                            y: f64|
         -> u32 {
            let key = (physical.to_string(), component, q(x), q(y));
            if let Some(&n) = key_map.get(&key) {
                // Keep the first component assignment for a welded node.
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

        for poly in polygons.iter() {
            let rects = polygon_minus_polygon(&poly.points, &[]);
            let mut pnodes: Vec<u32> = Vec::new();
            // Layer-based refinement: a finer/coarser cell for this layer.
            let cell_target = per_layer_cell
                .get(&poly.physical)
                .copied()
                .filter(|c| *c > 0.0)
                .unwrap_or(cell_target);
            for rect in rects {
                let w = rect.x1 - rect.x0;
                let h = rect.y1 - rect.y0;
                if w <= 0.0 || h <= 0.0 {
                    continue;
                }
                let nx = ((w / cell_target).ceil() as usize).clamp(1, cap);
                let ny = ((h / cell_target).ceil() as usize).clamp(1, cap);
                let dx = w / nx as f64;
                let dy = h / ny as f64;
                // Nodes.
                for j in 0..=ny {
                    for i in 0..=nx {
                        let n = add_node(
                            &mut net,
                            &mut key_map,
                            &poly.physical,
                            poly.component,
                            rect.x0 + dx * i as f64,
                            rect.y0 + dy * j as f64,
                        );
                        pnodes.push(n);
                    }
                }
                // Neighbor resistors. Horizontal: R = Rs * dx/dy; vertical: Rs * dy/dx.
                let idx = |i: usize, j: usize| -> u32 {
                    let key = (
                        poly.physical.clone(),
                        poly.component,
                        q(rect.x0 + dx * i as f64),
                        q(rect.y0 + dy * j as f64),
                    );
                    key_map[&key]
                };
                // Node-based sheet resistance. NOTE: this 5-point stencil uses the
                // full cell height `dy` as the edge face, so a full-height bar
                // reads `(r_sheet*L/W)*ny/(ny+1)` (the boundary rows each carry a
                // full cell). The exact fix is node-centred finite volume: use the
                // control-volume face (dy/2 on the domain-boundary rows) instead
                // of the cell height. That is a reformulation of this loop and of
                // `adaptive::build_network`'s `edge()` closure; it must keep 2-D
                // shapes (necks) correct, so it cannot be a blanket edge halving.
                let rh = poly.r_sheet * dx / dy;
                let rv = poly.r_sheet * dy / dx;
                for j in 0..=ny {
                    for i in 0..nx {
                        net.edges.push((idx(i, j), idx(i + 1, j), 1.0 / rh));
                    }
                }
                for j in 0..ny {
                    for i in 0..=nx {
                        net.edges.push((idx(i, j), idx(i, j + 1), 1.0 / rv));
                    }
                }
                // Record the quad cells (used to redraw the layout as a mosaic).
                for j in 0..ny {
                    for i in 0..nx {
                        net.cells.push([
                            idx(i, j),
                            idx(i + 1, j),
                            idx(i + 1, j + 1),
                            idx(i, j + 1),
                        ]);
                        net.cell_r_sheet.push(poly.r_sheet);
                        net.cell_level.push(0);
                        net.cell_parent.push(u32::MAX);
                    }
                }
            }
            net.poly_nodes.push(pnodes);
        }

        // Same-layer overlap/touch connections: join the nearest nodes.
        for &(a, b) in same_layer_connections {
            if a >= net.poly_nodes.len() || b >= net.poly_nodes.len() {
                continue;
            }
            let pa = &net.poly_nodes[a];
            let pb = &net.poly_nodes[b];
            if pa.is_empty() || pb.is_empty() {
                continue;
            }
            let (ax, ay) = (net.x[pa[0] as usize], net.y[pa[0] as usize]);
            if let Some(na) = nearest(pa, &net.x, &net.y, ax, ay) {
                // nearest node of b to node a
                let (nx, ny) = (net.x[na as usize], net.y[na as usize]);
                if let Some(nb) = nearest(pb, &net.x, &net.y, nx, ny) {
                    let r_sheet = 1.0; // small connection resistance (ohm)
                    net.edges.push((na, nb, 1.0 / r_sheet));
                }
            }
        }

        openrdson_core::log_debug!(
            "same-layer stitch links: {} (1 ohm each)",
            same_layer_connections.len()
        );

        // Vias: connect same-net nodes on the bottom and top layers. A point
        // via (zero footprint) connects the single nearest node on each layer; a
        // via with a footprint spreads its current over every node in the rect
        // [x±half_width] x [y±half_height] (the via is a uniform vertical
        // resistor, so the total conductance 1/R is split evenly across the
        // bottom x top node pairs).
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
            if let Some((d, s)) = net.channel_nodes(ch.drain, ch.source, ch.drain_component, ch.source_component) {
                if ch.resistance > 0.0 && ch.resistance.is_finite() {
                    net.edges.push((d, s, 1.0 / ch.resistance));
                }
            }
        }

        // The adaptive per-polygon grids do not guarantee that every same-net
        // polygon is joined. Merge any graph fragments that belong to the same
        // connectivity net so each net is one connected component.
        apply_fv_weights(&mut net);
        net.merge_nets_by_component();
        openrdson_core::log_info!(
            "sheet network: {} polygons -> {} nodes, {} edges",
            polygons.len(),
            net.x.len(),
            net.edges.len()
        );
        net
    }

    /// Connect graph fragments that share a connectivity component with a small
    /// resistor, so each electrical net is one connected component.
    pub fn merge_nets_by_component(&mut self) {
        // BTreeMap (not HashMap): the "main" fragment must be chosen
        // deterministically, or the merge links (and the solved result) change
        // from run to run.
        use std::collections::BTreeMap;
        let n = self.x.len();
        let mut parent: Vec<u32> = (0..n as u32).collect();
        fn find(parent: &mut [u32], mut x: u32) -> u32 {
            while parent[x as usize] != x {
                parent[x as usize] = parent[parent[x as usize] as usize];
                x = parent[x as usize];
            }
            x
        }
        for &(i, j, _) in &self.edges {
            let (ri, rj) = (find(&mut parent, i), find(&mut parent, j));
            if ri != rj {
                parent[ri as usize] = rj;
            }
        }
        // Group node indices by (connectivity component, graph root).
        let mut groups: BTreeMap<usize, BTreeMap<u32, Vec<u32>>> = BTreeMap::new();
        for i in 0..n {
            let c = self.component[i];
            if c == usize::MAX {
                continue;
            }
            let r = find(&mut parent, i as u32);
            groups
                .entry(c)
                .or_default()
                .entry(r)
                .or_default()
                .push(i as u32);
        }
        let mut extra: Vec<(u32, u32, f64)> = Vec::new();
        for comp_groups in groups.values() {
            let mut it = comp_groups.values();
            let Some(main) = it.next() else { continue };
            for other in it {
                // Nearest node pair between this fragment and the main one.
                let mut best: Option<(u32, u32)> = None;
                let mut bd = f64::INFINITY;
                for &a in other {
                    let (ax, ay) = (self.x[a as usize], self.y[a as usize]);
                    for &b in main {
                        let d = (self.x[b as usize] - ax).powi(2) + (self.y[b as usize] - ay).powi(2);
                        if d < bd {
                            bd = d;
                            best = Some((a, b));
                        }
                    }
                }
                if let Some((a, b)) = best {
                    extra.push((a, b, 1e3)); // ~0 ohm merge link (same net)
                }
            }
        }
        openrdson_core::log_debug!("merge links: {} (1 ohm each)", extra.len());
        self.edges.extend(extra);
    }

    /// Nearest RX nodes for a device's drain/source components.
    /// Nearest node on a connectivity component to `(x, y)`, regardless of
    /// layer. Used to re-project device terminals onto a rebuilt mesh.
    pub fn nearest_on_comp(&self, comp: usize, x: f64, y: f64) -> Option<u32> {
        let mut best: Option<(u32, f64)> = None;
        for n in 0..self.x.len() {
            if self.component[n] != comp {
                continue;
            }
            let d = (self.x[n] - x).powi(2) + (self.y[n] - y).powi(2);
            if best.map(|(_, bd)| d < bd).unwrap_or(true) {
                best = Some((n as u32, d));
            }
        }
        best.map(|(n, _)| n)
    }

    pub fn channel_nodes(
        &self,
        drain: Point,
        source: Point,
        drain_component: usize,
        source_component: usize,
    ) -> Option<(u32, u32)> {
        let nodes = self.physical_nodes.get("RX")?;
        let pick = |comp: usize, px: f64, py: f64| -> Option<u32> {
            let cand: Vec<u32> = nodes
                .iter()
                .copied()
                .filter(|&n| self.component[n as usize] == comp)
                .collect();
            nearest(&cand, &self.x, &self.y, px, py)
        };
        Some((
            pick(drain_component, drain.x, drain.y)?,
            pick(source_component, source.x, source.y)?,
        ))
    }

    /// Solve for node potentials with `a_nodes` at 1 V, `b_nodes` at 0 V, other
    /// nets grounded, and `extra` device edges injected.
    pub fn solve_potentials(
        &self,
        a_nodes: &[u32],
        b_nodes: &[u32],
        extra: &[(u32, u32, f64)],
        tol: f64,
        max_iter: usize,
    ) -> Result<Vec<f64>, String> {
        let n = self.x.len();
        let mut parent: Vec<u32> = (0..n as u32).collect();
        fn find(parent: &mut [u32], mut x: u32) -> u32 {
            while parent[x as usize] != x {
                parent[x as usize] = parent[parent[x as usize] as usize];
                x = parent[x as usize];
            }
            x
        }
        for &(i, j, _) in &self.edges {
            let (ri, rj) = (find(&mut parent, i), find(&mut parent, j));
            if ri != rj {
                parent[ri as usize] = rj;
            }
        }
        let ra = a_nodes.first().map(|&x| find(&mut parent, x));
        let rb = b_nodes.first().map(|&x| find(&mut parent, x));

        let mut a = self.to_sparse();
        for &(i, j, g) in extra {
            a.add_edge(i, j, g);
        }
        a.finalize();
        let mut fixed = BTreeMap::new();
        for &nd in a_nodes {
            fixed.insert(nd, 1.0);
        }
        for &nd in b_nodes {
            fixed.insert(nd, 0.0);
        }
        for i in 0..n {
            let r = find(&mut parent, i as u32);
            if Some(r) != ra && Some(r) != rb {
                fixed.insert(i as u32, 0.0);
            }
        }
        solve_dirichlet(&a, &fixed, tol, max_iter)
    }

    /// Solve with arbitrary Dirichlet potentials at the given nodes.
    ///
    /// Every connected component that contains a fixed node is held at the
    /// given potentials; components with no fixed node are grounded so the
    /// system stays non-singular. This is the general form of
    /// [`SheetNetwork::solve_potentials`] used by the multi-terminal bias.
    pub fn solve_fixed(
        &self,
        fixed: &[(u32, f64)],
        extra: &[(u32, u32, f64)],
        tol: f64,
        max_iter: usize,
    ) -> Result<Vec<f64>, String> {
        self.solve_fixed_warm(fixed, extra, None, tol, max_iter)
    }

    /// Like [`Self::solve_fixed`], but seeds the linear solve with the initial
    /// guess `x0` (length = node count). The direct Cholesky path (feature
    /// `faer-sheet`) ignores `x0`; the CG path uses it as its starting point.
    pub fn solve_fixed_warm(
        &self,
        fixed: &[(u32, f64)],
        extra: &[(u32, u32, f64)],
        x0: Option<&[f64]>,
        tol: f64,
        max_iter: usize,
    ) -> Result<Vec<f64>, String> {
        let n = self.x.len();
        let mut parent: Vec<u32> = (0..n as u32).collect();
        fn find(parent: &mut [u32], mut x: u32) -> u32 {
            while parent[x as usize] != x {
                parent[x as usize] = parent[parent[x as usize] as usize];
                x = parent[x as usize];
            }
            x
        }
        for &(i, j, _) in &self.edges {
            let (ri, rj) = (find(&mut parent, i), find(&mut parent, j));
            if ri != rj {
                parent[ri as usize] = rj;
            }
        }
        let mut a = self.to_sparse();
        for &(i, j, g) in extra {
            a.add_edge(i, j, g);
        }
        a.finalize();
        let mut bc = BTreeMap::new();
        for &(nd, v) in fixed {
            bc.insert(nd, v);
        }
        let mut has_fixed: std::collections::BTreeSet<u32> = std::collections::BTreeSet::new();
        for &(nd, _) in fixed {
            has_fixed.insert(find(&mut parent, nd));
        }
        for i in 0..n {
            let r = find(&mut parent, i as u32);
            if !has_fixed.contains(&r) {
                bc.insert(i as u32, 0.0);
            }
        }
        let x0_vec = x0.map(|x| x.to_vec()).unwrap_or_else(|| vec![0.0f64; n]);
        // 2D sheet systems favour a direct sparse Cholesky (O(n log n) fill);
        // the feature is off by default so the tool stays dependency-free.
        #[cfg(feature = "faer-sheet")]
        {
            openrdson_solver::solve_dirichlet_direct(&a, &bc, tol, max_iter)
        }
        #[cfg(not(feature = "faer-sheet"))]
        {
            solve_dirichlet_warm(&a, &bc, &x0_vec, tol, max_iter)
        }
    }

    /// Resistance between node sets with extra device edges injected.
    pub fn resistance_with_edges(
        &self,
        a_nodes: &[u32],
        b_nodes: &[u32],
        extra: &[(u32, u32, f64)],
        tol: f64,
        max_iter: usize,
    ) -> Result<f64, String> {
        let u = self.solve_potentials(a_nodes, b_nodes, extra, tol, max_iter)?;
        // Dissipated power summed over every edge (including injected channels).
        let mut p = 0.0;
        for &(i, j, g) in self.edges.iter().chain(extra.iter()) {
            let d = u[i as usize] - u[j as usize];
            p += g * d * d;
        }
        Ok(if p.abs() < 1e-300 { f64::INFINITY } else { 1.0 / p })
    }

    pub fn node_count(&self) -> usize {
        self.x.len()
    }

    pub fn cell_count(&self) -> usize {
        self.cells.len()
    }

    /// Per-layer mesh counts keyed by physical layer: `(nodes, same-layer edges,
    /// cells)`. Cross-layer (via) and channel edges are not attributed to a
    /// layer, so `edges` counts only edges whose endpoints share the layer.
    pub fn layer_stats(&self) -> BTreeMap<String, (usize, usize, usize)> {
        let mut out: BTreeMap<String, (usize, usize, usize)> = BTreeMap::new();
        for (layer, nodes) in &self.physical_nodes {
            out.entry(layer.clone()).or_default().0 = nodes.len();
        }
        for &(a, b, _) in &self.edges {
            let la = &self.physical[a as usize];
            if la == &self.physical[b as usize] {
                out.entry(la.clone()).or_default().1 += 1;
            }
        }
        for cell in &self.cells {
            out.entry(self.physical[cell[0] as usize].clone())
                .or_default()
                .2 += 1;
        }
        out
    }

    /// Physical layer, corner coordinates (meters), and cell extents for a cell.
    #[allow(clippy::type_complexity)]
    pub fn cell_info(&self, cell: [u32; 4]) -> (String, [(f64, f64); 4], f64, f64) {
        let c = [
            (self.x[cell[0] as usize], self.y[cell[0] as usize]),
            (self.x[cell[1] as usize], self.y[cell[1] as usize]),
            (self.x[cell[2] as usize], self.y[cell[2] as usize]),
            (self.x[cell[3] as usize], self.y[cell[3] as usize]),
        ];
        let dx = ((c[1].0 - c[0].0).abs() + (c[2].0 - c[3].0).abs()) / 2.0;
        let dy = ((c[3].1 - c[0].1).abs() + (c[2].1 - c[1].1).abs()) / 2.0;
        (
            self.physical[cell[0] as usize].clone(),
            c,
            dx.max(1e-12),
            dy.max(1e-12),
        )
    }

    pub fn edge_count(&self) -> usize {
        self.edges.len()
    }

    /// Nearest node on `physical` to a point.
    pub fn nearest_on(&self, physical: &str, x: f64, y: f64) -> Option<u32> {
        let nodes = self.physical_nodes.get(physical)?;
        nearest(nodes, &self.x, &self.y, x, y)
    }

    pub fn to_sparse(&self) -> Sparse {
        let mut a = Sparse::new(self.x.len());
        for &(i, j, g) in &self.edges {
            a.add_edge(i, j, g);
        }
        a.finalize();
        a
    }

    /// Two-terminal resistance between node sets.
    pub fn resistance(
        &self,
        a_nodes: &[u32],
        b_nodes: &[u32],
        tol: f64,
        max_iter: usize,
    ) -> Result<f64, String> {
        let a = self.to_sparse();
        let mut fixed = BTreeMap::new();
        for &n in a_nodes {
            fixed.insert(n, 1.0);
        }
        for &n in b_nodes {
            fixed.insert(n, 0.0);
        }
        let u = solve_dirichlet(&a, &fixed, tol, max_iter)?;
        let mut p = 0.0;
        for &(i, j, g) in &self.edges {
            let d = u[i as usize] - u[j as usize];
            p += g * d * d;
        }
        Ok(if p.abs() < 1e-300 {
            f64::INFINITY
        } else {
            1.0 / p
        })
    }
}

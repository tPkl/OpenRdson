//! Finite-element meshing for OpenRDSon.
//!
//! v0 implements the "2.5D first" strategy:
//! - axis-aligned box regions become **structured hexahedral grids**;
//! - general prism regions (non-rectangular footprints) are ear-clipped into
//!   **triangular prisms**.
//!
//! This covers flat metal/poly/diffusion layers efficiently and exactly.
//! Cross-region conforming interfaces, vias/contacts and conformal dielectrics
//! are later milestones.

use openrdson_core::geometry::{polygon_signed_area, Point, Solid, SolidModel, SolidRegion};
use openrdson_core::mesh::{Element, ElementType, Mesh, Point3};

#[derive(Debug, Clone, Copy)]
pub struct MeshConfig {
    /// Target lateral element size (m).
    pub lateral_cell: f64,
    /// Target vertical element size (m).
    pub z_cell: f64,
}

impl Default for MeshConfig {
    fn default() -> Self {
        Self {
            lateral_cell: 1e-6,
            z_cell: 0.1e-6,
        }
    }
}

/// Per-layer mesh resolution override.
#[derive(Debug, Clone, Copy)]
pub struct LayerMesh {
    pub lateral_cell: f64,
    pub z_cell: f64,
}

/// R3D-style meshing options: a default resolution plus per-physical-layer
/// overrides (e.g. finer cells on Poly1/gate, coarser on thick top metal).
#[derive(Debug, Clone, Default)]
pub struct MeshOptions {
    pub default: MeshConfig,
    pub per_layer: std::collections::BTreeMap<String, LayerMesh>,
}

impl MeshOptions {
    pub fn with_layer(mut self, layer: &str, lateral_cell: f64, z_cell: f64) -> Self {
        self.per_layer.insert(
            layer.to_string(),
            LayerMesh {
                lateral_cell,
                z_cell,
            },
        );
        self
    }

    /// Resolution to use for a region, based on its physical layer.
    pub fn for_region(&self, region: &SolidRegion) -> MeshConfig {
        region
            .layer
            .as_ref()
            .and_then(|l| self.per_layer.get(l))
            .map(|m| MeshConfig {
                lateral_cell: m.lateral_cell,
                z_cell: m.z_cell,
            })
            .unwrap_or(self.default)
    }
}

/// Mesh a solid model using per-layer resolution overrides.
pub fn mesh_model_with_options(model: &SolidModel, opts: &MeshOptions) -> (Mesh, Vec<String>) {
    let mut mesh = Mesh::default();
    let mut diagnostics = Vec::new();
    for (i, region) in model.regions.iter().enumerate() {
        let cfg = opts.for_region(region);
        match mesh_region(region, &cfg) {
            Ok(m) => merge_mesh(&mut mesh, &m),
            Err(e) => diagnostics.push(format!(
                "region {i} ({}) failed: {e}",
                model.material_name(region.material_id).unwrap_or("?")
            )),
        }
    }
    (mesh, diagnostics)
}

/// Mesh a single tagged solid region.
pub fn mesh_region(region: &SolidRegion, cfg: &MeshConfig) -> Result<Mesh, String> {
    match &region.solid {
        Solid::Box { .. } => mesh_box(region, cfg),
        Solid::Prism { .. } => mesh_prism(region, cfg),
    }
}

/// Mesh every region of a solid model, accumulating per-region diagnostics.
pub fn mesh_model(model: &SolidModel, cfg: &MeshConfig) -> (Mesh, Vec<String>) {
    openrdson_core::log_info!(
        "meshing {} regions (per-region 2.5D, lateral {:.3} um)",
        model.regions.len(),
        cfg.lateral_cell * 1e6
    );
    let mut mesh = Mesh::default();
    let mut diagnostics = Vec::new();
    for (i, region) in model.regions.iter().enumerate() {
        match mesh_region(region, cfg) {
            Ok(m) => merge_mesh(&mut mesh, &m),
            Err(e) => diagnostics.push(format!(
                "region {i} ({}) failed: {e}",
                model.material_name(region.material_id).unwrap_or("?")
            )),
        }
    }
    (mesh, diagnostics)
}

/// Merge coincident nodes (within `tol`) so adjacent regions share interface
/// nodes and the assembled system is conforming.
pub fn weld_nodes(mesh: &Mesh, tol: f64) -> Mesh {
    use std::collections::HashMap;
    let q = |v: f64| (v / tol).round() as i64;
    let mut map: HashMap<(i64, i64, i64), u32> = HashMap::new();
    let mut new_nodes: Vec<Point3> = Vec::new();
    let mut remap = vec![0u32; mesh.nodes.len()];
    for (i, p) in mesh.nodes.iter().enumerate() {
        let key = (q(p.x), q(p.y), q(p.z));
        let idx = *map.entry(key).or_insert_with(|| {
            new_nodes.push(*p);
            (new_nodes.len() - 1) as u32
        });
        remap[i] = idx;
    }
    let mut elements = Vec::with_capacity(mesh.elements.len());
    for e in &mesh.elements {
        let mut e = e.clone();
        for n in &mut e.nodes {
            *n = remap[*n as usize];
        }
        elements.push(e);
    }
    Mesh {
        nodes: new_nodes,
        elements,
        sets: Default::default(),
    }
}

fn merge_mesh(dst: &mut Mesh, src: &Mesh) {
    let offset = dst.nodes.len() as u32;
    dst.nodes.extend_from_slice(&src.nodes);
    for e in &src.elements {
        let mut e = e.clone();
        for n in &mut e.nodes {
            *n += offset;
        }
        dst.elements.push(e);
    }
    for (k, v) in &src.sets {
        dst.sets
            .entry(k.clone())
            .or_default()
            .extend(v.iter().map(|x| x + offset));
    }
}

fn linspace(min: f64, max: f64, n: usize) -> Vec<f64> {
    let mut v = Vec::with_capacity(n + 1);
    if n == 0 {
        v.push(min);
        return v;
    }
    let step = (max - min) / n as f64;
    for i in 0..=n {
        v.push(min + step * i as f64);
    }
    if let Some(last) = v.last_mut() {
        *last = max;
    }
    v
}

/// Number of cells needed to cover `extent` with cells of size `cell`.
fn cells_for(extent: f64, cell: f64) -> usize {
    let q = extent / cell;
    // Shrink slightly so an exact multiple does not round up due to FP error.
    ((q * (1.0 - 1e-12)).ceil() as usize).max(1)
}

fn mesh_box(region: &SolidRegion, cfg: &MeshConfig) -> Result<Mesh, String> {
    let b = region.solid.footprint_bbox();
    let z0 = region.solid.z_bottom();
    let z1 = region.solid.z_top();
    if b.width() <= 0.0 || b.height() <= 0.0 || z1 <= z0 {
        return Err("degenerate box region".into());
    }
    if cfg.lateral_cell <= 0.0 || cfg.z_cell <= 0.0 {
        return Err("mesh cell sizes must be positive".into());
    }

    let nx = cells_for(b.width(), cfg.lateral_cell);
    let ny = cells_for(b.height(), cfg.lateral_cell);
    let nz = cells_for(z1 - z0, cfg.z_cell);

    let xs = linspace(b.min_x, b.max_x, nx);
    let ys = linspace(b.min_y, b.max_y, ny);
    let zs = linspace(z0, z1, nz);

    let mut nodes = Vec::with_capacity((nx + 1) * (ny + 1) * (nz + 1));
    for k in 0..=nz {
        for j in 0..=ny {
            for i in 0..=nx {
                nodes.push(Point3::new(xs[i], ys[j], zs[k]));
            }
        }
    }

    let stride_y = nx + 1;
    let stride_z = (nx + 1) * (ny + 1);
    let idx = |i: usize, j: usize, k: usize| -> u32 { (i + stride_y * j + stride_z * k) as u32 };

    let mut elements = Vec::with_capacity(nx * ny * nz);
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let n = vec![
                    idx(i, j, k),
                    idx(i + 1, j, k),
                    idx(i + 1, j + 1, k),
                    idx(i, j + 1, k),
                    idx(i, j, k + 1),
                    idx(i + 1, j, k + 1),
                    idx(i + 1, j + 1, k + 1),
                    idx(i, j + 1, k + 1),
                ];
                elements.push(Element {
                    kind: ElementType::Hex,
                    nodes: n,
                    material_id: region.material_id,
                    net: region.net.clone(),
                    layer: region.layer.clone(),
                    device_ref: region.device_ref.clone(),
                    source_polygon: region.source_polygon,
                });
            }
        }
    }

    Ok(Mesh {
        nodes,
        elements,
        sets: Default::default(),
    })
}

/// Remove consecutive duplicate vertices and the closing duplicate of a ring.
fn dedup_ring(points: &[Point]) -> Vec<Point> {
    let mut out: Vec<Point> = Vec::with_capacity(points.len());
    for p in points {
        if let Some(last) = out.last() {
            if (last.x - p.x).abs() < f64::EPSILON && (last.y - p.y).abs() < f64::EPSILON {
                continue;
            }
        }
        out.push(*p);
    }
    if out.len() > 1 {
        let first = out[0];
        let last = *out.last().unwrap();
        if (first.x - last.x).abs() < f64::EPSILON && (first.y - last.y).abs() < f64::EPSILON {
            out.pop();
        }
    }
    out
}

fn point_in_triangle(p: Point, a: Point, b: Point, c: Point) -> bool {
    let sign = |p1: Point, p2: Point, p3: Point| {
        (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y)
    };
    let d1 = sign(p, a, b);
    let d2 = sign(p, b, c);
    let d3 = sign(p, c, a);
    let has_neg = d1 < 0.0 || d2 < 0.0 || d3 < 0.0;
    let has_pos = d1 > 0.0 || d2 > 0.0 || d3 > 0.0;
    !(has_neg && has_pos)
}

/// Ear-clip a simple polygon (CCW or CW) into triangles.
///
/// Returns triangles as triples of indices into the *deduplicated* ring, which
/// is also returned so callers can build nodes consistently.
pub fn triangulate_polygon(points: &[Point]) -> Result<(Vec<Point>, Vec<[usize; 3]>), String> {
    let pts = dedup_ring(points);
    let n = pts.len();
    if n < 3 {
        return Err("polygon has fewer than 3 distinct vertices".into());
    }
    let mut idx: Vec<usize> = (0..n).collect();
    if polygon_signed_area(&pts) < 0.0 {
        idx.reverse();
    }

    let mut tris: Vec<[usize; 3]> = Vec::with_capacity(n - 2);
    let mut guard = 0usize;
    while idx.len() > 3 {
        let m = idx.len();
        let mut found = false;
        for i in 0..m {
            let ia = idx[(i + m - 1) % m];
            let ib = idx[i];
            let ic = idx[(i + 1) % m];
            let (a, b, c) = (pts[ia], pts[ib], pts[ic]);
            // Convex corner (CCW): cross > 0.
            let cross = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
            if cross <= 0.0 {
                continue;
            }
            let mut empty = true;
            for &j in &idx {
                if j == ia || j == ib || j == ic {
                    continue;
                }
                if point_in_triangle(pts[j], a, b, c) {
                    empty = false;
                    break;
                }
            }
            if empty {
                tris.push([ia, ib, ic]);
                idx.remove(i);
                found = true;
                break;
            }
        }
        if !found {
            return Err("ear clipping stalled (non-simple or degenerate polygon?)".into());
        }
        guard += 1;
        if guard > n * n + 8 {
            return Err("ear clipping exceeded iteration budget".into());
        }
    }
    tris.push([idx[0], idx[1], idx[2]]);
    Ok((pts, tris))
}

/// True if every edge of the ring is axis-aligned.
fn ring_is_rectilinear(ring: &[Point]) -> bool {
    let n = ring.len();
    if n < 3 {
        return false;
    }
    let scale = ring
        .iter()
        .fold(0.0f64, |m, p| m.max(p.x.abs()).max(p.y.abs()))
        .max(1.0);
    let eps = 1e-9 * scale;
    for i in 0..n {
        let a = ring[i];
        let b = ring[(i + 1) % n];
        if (a.x - b.x).abs() > eps && (a.y - b.y).abs() > eps {
            return false;
        }
    }
    true
}

/// Exact scanline (trapezoidal) decomposition of a ring, meshed as triangular
/// prisms.
///
/// Between consecutive distinct vertex `y` values the polygon is bounded on the
/// left and right by single linear edges, so each interior interval is an exact
/// trapezoid. Crossings are paired by the even-odd rule, which makes this
/// correct for simple rings **and** self-touching / keyhole rings (holes and
/// zero-width slits included). The summed trapezoid area equals the polygon's
/// material area exactly.
fn mesh_scanline(region: &SolidRegion, ring: &[Point], z_bottom: f64, z_top: f64) -> Result<Mesh, String> {
    let n = ring.len();
    if n < 3 {
        return Err("scanline: fewer than 3 vertices".into());
    }
    let mut ys: Vec<f64> = ring.iter().map(|p| p.y).collect();
    ys.sort_by(|a, b| a.partial_cmp(b).unwrap());
    ys.dedup_by(|a, b| (*a - *b).abs() < 1e-15);
    if ys.len() < 2 {
        return Err("scanline: degenerate y extent".into());
    }

    let eval = |e: usize, y: f64| -> f64 {
        let a = ring[e];
        let b = ring[(e + 1) % n];
        if (b.y - a.y).abs() < 1e-300 {
            a.x
        } else {
            a.x + (y - a.y) / (b.y - a.y) * (b.x - a.x)
        }
    };

    let mut nodes: Vec<Point3> = Vec::new();
    let mut elements: Vec<Element> = Vec::new();
    let mut push_quad = |a: Point, b: Point, c: Point, d: Point| {
        // Quad (a,b,c,d) CCW in XY, extruded to z_bottom/z_top as two triangular
        // prisms. Skip degenerate (near-zero area) quads.
        let area = 0.5
            * ((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)
                + (c.x - a.x) * (d.y - a.y) - (c.y - a.y) * (d.x - a.x));
        if area.abs() < 1e-30 {
            return;
        }
        let base = nodes.len() as u32;
        let pts = [a, b, c, d];
        for p in pts {
            nodes.push(Point3::new(p.x, p.y, z_bottom));
        }
        for p in pts {
            nodes.push(Point3::new(p.x, p.y, z_top));
        }
        // Triangles (0,1,2) and (0,2,3).
        for tri in [[0u32, 1, 2], [0, 2, 3]] {
            elements.push(Element {
                kind: ElementType::Prism,
                nodes: vec![
                    base + tri[0],
                    base + tri[1],
                    base + tri[2],
                    base + tri[0] + 4,
                    base + tri[1] + 4,
                    base + tri[2] + 4,
                ],
                material_id: region.material_id,
                net: region.net.clone(),
                layer: region.layer.clone(),
                device_ref: region.device_ref.clone(),
                source_polygon: region.source_polygon,
            });
        }
    };

    for k in 0..ys.len() - 1 {
        let (y0, y1) = (ys[k], ys[k + 1]);
        if y1 <= y0 {
            continue;
        }
        let ymid = 0.5 * (y0 + y1);
        // Collect crossings (x, edge index) of edges spanning the mid-line.
        let mut crossings: Vec<(f64, usize)> = Vec::new();
        for e in 0..n {
            let a = ring[e];
            let b = ring[(e + 1) % n];
            let spans = (a.y <= ymid && b.y > ymid) || (b.y <= ymid && a.y > ymid);
            if spans {
                crossings.push((eval(e, ymid), e));
            }
        }
        crossings.sort_by(|p, q| p.0.partial_cmp(&q.0).unwrap());
        if crossings.len() % 2 != 0 {
            return Err("scanline: odd crossing count (open or self-intersecting ring)".into());
        }
        for pair in crossings.chunks_exact(2) {
            let e_l = pair[0].1;
            let e_r = pair[1].1;
            let (xl0, xl1) = (eval(e_l, y0), eval(e_l, y1));
            let (xr0, xr1) = (eval(e_r, y0), eval(e_r, y1));
            let a = Point::new(xl0, y0);
            let b = Point::new(xr0, y0);
            let c = Point::new(xr1, y1);
            let d = Point::new(xl1, y1);
            push_quad(a, b, c, d);
        }
    }

    if elements.is_empty() {
        return Err("scanline: no interior trapezoids found".into());
    }
    Ok(Mesh {
        nodes,
        elements,
        sets: Default::default(),
    })
}

fn mesh_prism(region: &SolidRegion, _cfg: &MeshConfig) -> Result<Mesh, String> {
    let (footprint, z_bottom, z_top) = match &region.solid {
        Solid::Prism {
            footprint,
            z_bottom,
            z_top,
        } => (footprint, *z_bottom, *z_top),
        _ => return Err("mesh_prism called on a non-prism".into()),
    };
    if z_top <= z_bottom {
        return Err("degenerate prism (z_top <= z_bottom)".into());
    }
    let ring = dedup_ring(footprint);
    // Simple polygons get an ear-clipped triangular-prism mesh; rectilinear or
    // self-touching/keyhole rings use the exact scanline decomposition.
    if ring_is_rectilinear(&ring) {
        return mesh_scanline(region, &ring, z_bottom, z_top);
    }
    match triangulate_polygon(&ring) {
        Ok((ring, tris)) => mesh_tris(region, &ring, &tris, z_bottom, z_top),
        Err(_) => mesh_scanline(region, &ring, z_bottom, z_top),
    }
}

fn mesh_tris(
    region: &SolidRegion,
    ring: &[Point],
    tris: &[[usize; 3]],
    z_bottom: f64,
    z_top: f64,
) -> Result<Mesh, String> {
    // Nodes: ring points at z_bottom (0..n) and z_top (n..2n).
    let n = ring.len();
    let mut nodes = Vec::with_capacity(2 * n);
    for p in ring {
        nodes.push(Point3::new(p.x, p.y, z_bottom));
    }
    for p in ring {
        nodes.push(Point3::new(p.x, p.y, z_top));
    }

    let mut elements = Vec::with_capacity(tris.len());
    for t in tris {
        elements.push(Element {
            kind: ElementType::Prism,
            nodes: vec![
                t[0] as u32,
                t[1] as u32,
                t[2] as u32,
                (t[0] + n) as u32,
                (t[1] + n) as u32,
                (t[2] + n) as u32,
            ],
            material_id: region.material_id,
            net: region.net.clone(),
            layer: region.layer.clone(),
            device_ref: region.device_ref.clone(),
            source_polygon: region.source_polygon,
        });
    }

    Ok(Mesh {
        nodes,
        elements,
        sets: Default::default(),
    })
}

fn triangle_area(a: Point3, b: Point3, c: Point3) -> f64 {
    let ux = b.x - a.x;
    let uy = b.y - a.y;
    let uz = b.z - a.z;
    let vx = c.x - a.x;
    let vy = c.y - a.y;
    let vz = c.z - a.z;
    let cx = uy * vz - uz * vy;
    let cy = uz * vx - ux * vz;
    let cz = ux * vy - uy * vx;
    0.5 * (cx * cx + cy * cy + cz * cz).sqrt()
}

/// Volume of a mesh (m³).
///
/// Exact for the axis-aligned hexes and z-extruded triangular prisms produced
/// here; other element types fall back to their node bounding box.
pub fn mesh_volume(mesh: &Mesh) -> f64 {
    let mut total = 0.0;
    for e in &mesh.elements {
        match e.kind {
            ElementType::Prism if e.nodes.len() >= 6 => {
                let p = |i: usize| mesh.nodes[e.nodes[i] as usize];
                let area = triangle_area(p(0), p(1), p(2));
                let h = ((p(3).z - p(0).z).abs() + (p(4).z - p(1).z).abs() + (p(5).z - p(2).z).abs())
                    / 3.0;
                total += area * h;
            }
            _ => {
                let mut min = Point3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
                let mut max = Point3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
                for &ni in &e.nodes {
                    let p = mesh.nodes[ni as usize];
                    min.x = min.x.min(p.x);
                    min.y = min.y.min(p.y);
                    min.z = min.z.min(p.z);
                    max.x = max.x.max(p.x);
                    max.y = max.y.max(p.y);
                    max.z = max.z.max(p.z);
                }
                total += (max.x - min.x) * (max.y - min.y) * (max.z - min.z);
            }
        }
    }
    total
}

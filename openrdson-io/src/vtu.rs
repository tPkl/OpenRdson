//! Minimal, dependency-free VTK XML (`.vtu`) writer for scientific viewers
//! such as ParaView.
//!
//! GDS is a geometry/mask format: each polygon carries exactly one fill colour,
//! so a scalar field can only ever be shown as discrete bands. A `.vtu` carries
//! per-vertex (and per-cell) arrays on an actual mesh, letting the viewer
//! Gouraud-shade a continuous, interpolated colormap and offer a dropdown to
//! switch between quantities.
//!
//! The output is the ASCII form of the VTK XML `UnstructuredGrid` format
//! (version 0.1). No external crate is required.

use std::io::{self, Write};
use std::path::Path;

use openrdson_core::mesh::{ElementType, Mesh};

/// VTK cell type identifiers (the subset this project emits).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VtkCellType {
    Vertex,
    Line,
    Triangle,
    Quad,
    Tet,
    Hex,
    Prism,
    Pyramid,
}

impl VtkCellType {
    /// The integer id VTK uses in the `types` array.
    pub fn vtk_id(self) -> u8 {
        match self {
            VtkCellType::Vertex => 1,
            VtkCellType::Line => 3,
            VtkCellType::Triangle => 5,
            VtkCellType::Quad => 9,
            VtkCellType::Tet => 10,
            VtkCellType::Hex => 12,
            VtkCellType::Prism => 13,
            VtkCellType::Pyramid => 14,
        }
    }
}

impl From<ElementType> for VtkCellType {
    fn from(e: ElementType) -> Self {
        match e {
            ElementType::Tet => VtkCellType::Tet,
            ElementType::Hex => VtkCellType::Hex,
            ElementType::Prism => VtkCellType::Prism,
            ElementType::Pyramid => VtkCellType::Pyramid,
        }
    }
}

/// A VTK unstructured grid with named point and cell data arrays.
///
/// Each `point_data` array must have one value per point; each `cell_data`
/// array one value per cell. `write_vtu` returns an error if that does not hold.
#[derive(Debug, Clone, Default)]
pub struct VtuMesh {
    pub points: Vec<[f64; 3]>,
    pub cells: Vec<(VtkCellType, Vec<u32>)>,
    /// `(name, values)` shown as selectable scalars on the mesh vertices.
    pub point_data: Vec<(String, Vec<f64>)>,
    /// `(name, values)` shown as selectable scalars on the cells.
    pub cell_data: Vec<(String, Vec<f64>)>,
}

impl VtuMesh {
    /// Build a new mesh keeping only the cells whose `keep[cell]` is true,
    /// remapping and compacting the point set (and carrying the data arrays).
    ///
    /// Used to split a combined export into one file per physical layer.
    pub fn subset_cells(&self, keep: &[bool]) -> VtuMesh {
        let mut remap = vec![usize::MAX; self.points.len()];
        let mut old_of_new: Vec<usize> = Vec::new();
        let mut points: Vec<[f64; 3]> = Vec::new();
        let mut cells: Vec<(VtkCellType, Vec<u32>)> = Vec::new();
        for (ci, (t, nodes)) in self.cells.iter().enumerate() {
            if !keep.get(ci).copied().unwrap_or(false) {
                continue;
            }
            let mut nn = Vec::with_capacity(nodes.len());
            for &n in nodes {
                let r = &mut remap[n as usize];
                if *r == usize::MAX {
                    *r = old_of_new.len();
                    old_of_new.push(n as usize);
                    points.push(self.points[n as usize]);
                }
                nn.push(*r as u32);
            }
            cells.push((*t, nn));
        }
        let point_data = self
            .point_data
            .iter()
            .map(|(name, vals)| {
                (
                    name.clone(),
                    old_of_new.iter().map(|&o| vals[o]).collect::<Vec<f64>>(),
                )
            })
            .collect();
        let cell_data = self
            .cell_data
            .iter()
            .map(|(name, vals)| {
                (
                    name.clone(),
                    keep.iter()
                        .enumerate()
                        .filter(|(_, k)| **k)
                        .map(|(i, _)| vals[i])
                        .collect::<Vec<f64>>(),
                )
            })
            .collect();
        VtuMesh {
            points,
            cells,
            point_data,
            cell_data,
        }
    }
}

/// Hexahedron faces (VTK node order: 0–3 bottom, 4–7 top).
const HEX_FACES: [&[usize]; 6] = [
    &[0, 1, 2, 3],
    &[4, 5, 6, 7],
    &[0, 1, 5, 4],
    &[1, 2, 6, 5],
    &[2, 3, 7, 6],
    &[3, 0, 4, 7],
];

/// Wedge/prism faces (VTK node order: 0–2 bottom triangle, 3–5 top triangle).
const PRISM_FACES: [&[usize]; 5] = [
    &[0, 1, 2],
    &[3, 4, 5],
    &[0, 1, 4, 3],
    &[1, 2, 5, 4],
    &[2, 0, 3, 5],
];

/// Gradient of a nodal field and the element volume, on a convex polyhedron,
/// by the divergence theorem: `∇u = (1/V) Σ_faces u_face · A_face`. Face normals
/// are flipped to point away from the element centroid, so any consistent node
/// ordering works.
fn poly_grad_vol(pts: &[[f64; 3]], u: &[f64], faces: &[&[usize]]) -> ([f64; 3], f64) {
    let n = pts.len() as f64;
    let mut cen = [0.0f64; 3];
    for p in pts {
        cen[0] += p[0];
        cen[1] += p[1];
        cen[2] += p[2];
    }
    for c in &mut cen {
        *c /= n;
    }

    let mut g = [0.0f64; 3];
    let mut ca = 0.0f64;
    for face in faces {
        let m = face.len();
        let mut a = [0.0f64; 3];
        let mut usum = 0.0f64;
        let mut fc = [0.0f64; 3];
        for i in 0..m {
            let p = pts[face[i]];
            let q = pts[face[(i + 1) % m]];
            a[0] += p[1] * q[2] - p[2] * q[1];
            a[1] += p[2] * q[0] - p[0] * q[2];
            a[2] += p[0] * q[1] - p[1] * q[0];
            usum += u[face[i]];
            fc[0] += p[0];
            fc[1] += p[1];
            fc[2] += p[2];
        }
        for k in 0..3 {
            a[k] *= 0.5;
            fc[k] /= m as f64;
        }
        // Make the face normal point outward from the element centroid.
        let outward =
            (fc[0] - cen[0]) * a[0] + (fc[1] - cen[1]) * a[1] + (fc[2] - cen[2]) * a[2];
        if outward < 0.0 {
            for k in 0..3 {
                a[k] = -a[k];
            }
        }
        let uf = usum / m as f64;
        for k in 0..3 {
            g[k] += uf * a[k];
        }
        ca += fc[0] * a[0] + fc[1] * a[1] + fc[2] * a[2];
    }

    let vol = ca / 3.0;
    if vol.abs() < 1e-300 {
        return ([0.0; 3], 0.0);
    }
    ([g[0] / vol, g[1] / vol, g[2] / vol], vol)
}

/// Gradient of a nodal field on a convex polyhedron (see [`poly_grad_vol`]).
fn poly_gradient(pts: &[[f64; 3]], u: &[f64], faces: &[&[usize]]) -> [f64; 3] {
    poly_grad_vol(pts, u, faces).0
}

/// Per-element current-density magnitude `|J| = σ|∇u|` on a hex/prism mesh.
///
/// `potential` is one value per mesh node, `sigma` one conductivity per
/// material id. Elements of other types get `0`. This is the field a viewer
/// colours the 3D mesh by, alongside the nodal `potential`.
pub fn cell_current_density(mesh: &Mesh, potential: &[f64], sigma: &[f64]) -> Vec<f64> {
    mesh.elements
        .iter()
        .map(|e| {
            let faces: &[&[usize]] = match (e.kind, e.nodes.len()) {
                (ElementType::Hex, 8) => HEX_FACES.as_slice(),
                (ElementType::Prism, 6) => PRISM_FACES.as_slice(),
                _ => return 0.0,
            };
            let pts: Vec<[f64; 3]> = e
                .nodes
                .iter()
                .map(|&n| {
                    let p = mesh.nodes[n as usize];
                    [p.x, p.y, p.z]
                })
                .collect();
            let un: Vec<f64> = e.nodes.iter().map(|&n| potential[n as usize]).collect();
            let grad = poly_gradient(&pts, &un, faces);
            let s = sigma.get(e.material_id as usize).copied().unwrap_or(0.0);
            s * (grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2]).sqrt()
        })
        .collect()
}

/// Per-element Joule power (W) and effective current (A) from a nodal potential
/// field.
///
/// `power = σ|∇u|²·V_e`; `current = power / ΔV_e` where `ΔV_e` is the potential
/// drop across the element (max − min node value). Elements of other types get
/// `0`. Requires `sigma` per material id and `potential` per node.
pub fn cell_power_current(mesh: &Mesh, potential: &[f64], sigma: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let mut power = Vec::with_capacity(mesh.elements.len());
    let mut current = Vec::with_capacity(mesh.elements.len());
    for e in &mesh.elements {
        let faces: &[&[usize]] = match (e.kind, e.nodes.len()) {
            (ElementType::Hex, 8) => HEX_FACES.as_slice(),
            (ElementType::Prism, 6) => PRISM_FACES.as_slice(),
            _ => {
                power.push(0.0);
                current.push(0.0);
                continue;
            }
        };
        let pts: Vec<[f64; 3]> = e
            .nodes
            .iter()
            .map(|&n| {
                let p = mesh.nodes[n as usize];
                [p.x, p.y, p.z]
            })
            .collect();
        let un: Vec<f64> = e.nodes.iter().map(|&n| potential[n as usize]).collect();
        let (grad, vol) = poly_grad_vol(&pts, &un, faces);
        let s = sigma.get(e.material_id as usize).copied().unwrap_or(0.0);
        let g2 = grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2];
        let p = s * g2 * vol.abs();
        let hi = un.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let lo = un.iter().cloned().fold(f64::INFINITY, f64::min);
        let dv = (hi - lo).abs();
        let i = if dv > 1e-15 { p / dv } else { 0.0 };
        power.push(p);
        current.push(i);
    }
    (power, current)
}

fn escape(name: &str) -> String {
    name.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

fn write_f64_array<W: Write>(w: &mut W, values: &[f64]) -> io::Result<()> {
    for (i, v) in values.iter().enumerate() {
        if i % 12 == 0 {
            write!(w, "\n         ")?;
        }
        write!(w, " {v:.9e}")?;
    }
    Ok(())
}

/// Write a [`VtuMesh`] to `path` as an ASCII VTK XML unstructured grid.
pub fn write_vtu<P: AsRef<Path>>(path: P, mesh: &VtuMesh) -> io::Result<()> {
    if let Some((name, v)) = mesh.point_data.iter().find(|(_, v)| v.len() != mesh.points.len()) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "vtu: point data '{name}' has {} values but there are {} points",
                v.len(),
                mesh.points.len()
            ),
        ));
    }
    if let Some((name, v)) = mesh.cell_data.iter().find(|(_, v)| v.len() != mesh.cells.len()) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "vtu: cell data '{name}' has {} values but there are {} cells",
                v.len(),
                mesh.cells.len()
            ),
        ));
    }

    let mut f = io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(f, "<?xml version=\"1.0\"?>")?;
    writeln!(
        f,
        "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"
    )?;
    writeln!(f, "  <UnstructuredGrid>")?;
    writeln!(
        f,
        "    <Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">",
        mesh.points.len(),
        mesh.cells.len()
    )?;

    // Points.
    writeln!(f, "      <Points>")?;
    write!(
        f,
        "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"
    )?;
    for (i, p) in mesh.points.iter().enumerate() {
        if i % 4 == 0 {
            write!(f, "\n         ")?;
        }
        write!(f, " {:.9e} {:.9e} {:.9e}", p[0], p[1], p[2])?;
    }
    writeln!(f, "\n        </DataArray>")?;
    writeln!(f, "      </Points>")?;

    // Cells: connectivity, offsets, types.
    writeln!(f, "      <Cells>")?;
    write!(
        f,
        "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">"
    )?;
    for (i, (_, nodes)) in mesh.cells.iter().enumerate() {
        if i % 12 == 0 {
            write!(f, "\n         ")?;
        }
        for &n in nodes {
            write!(f, " {n}")?;
        }
    }
    writeln!(f, "\n        </DataArray>")?;

    write!(
        f,
        "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">"
    )?;
    let mut offset = 0usize;
    for (i, (_, nodes)) in mesh.cells.iter().enumerate() {
        offset += nodes.len();
        if i % 12 == 0 {
            write!(f, "\n         ")?;
        }
        write!(f, " {offset}")?;
    }
    writeln!(f, "\n        </DataArray>")?;

    write!(
        f,
        "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">"
    )?;
    for (i, (t, _)) in mesh.cells.iter().enumerate() {
        if i % 24 == 0 {
            write!(f, "\n         ")?;
        }
        write!(f, " {}", t.vtk_id())?;
    }
    writeln!(f, "\n        </DataArray>")?;
    writeln!(f, "      </Cells>")?;

    // Point data.
    if !mesh.point_data.is_empty() {
        writeln!(f, "      <PointData Scalars=\"potential_V\">")?;
        for (name, values) in &mesh.point_data {
            write!(
                f,
                "        <DataArray type=\"Float64\" Name=\"{}\" format=\"ascii\">",
                escape(name)
            )?;
            write_f64_array(&mut f, values)?;
            writeln!(f, "\n        </DataArray>")?;
        }
        writeln!(f, "      </PointData>")?;
    }

    // Cell data.
    if !mesh.cell_data.is_empty() {
        writeln!(f, "      <CellData>")?;
        for (name, values) in &mesh.cell_data {
            write!(
                f,
                "        <DataArray type=\"Float64\" Name=\"{}\" format=\"ascii\">",
                escape(name)
            )?;
            write_f64_array(&mut f, values)?;
            writeln!(f, "\n        </DataArray>")?;
        }
        writeln!(f, "      </CellData>")?;
    }

    writeln!(f, "    </Piece>")?;
    writeln!(f, "  </UnstructuredGrid>")?;
    writeln!(f, "</VTKFile>")?;
    f.flush()
}

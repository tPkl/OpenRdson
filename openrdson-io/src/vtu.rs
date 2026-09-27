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

use openrdson_core::mesh::ElementType;

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

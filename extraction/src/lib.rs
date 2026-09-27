//! Resistive network extraction for OpenRDSon.
//!
//! The [`sheet`] module implements the full-die 2.5D sheet-resistance network;
//! [`adaptive`] adds quadtree refinement. The free function
//! [`two_terminal_resistance`] is the analytic acceptance helper shared with the
//! validation crate.

pub mod sheet;
pub mod adaptive;
pub use sheet::{group_vias, map_solution, SheetChannel, SheetNetwork, SheetPolygon, SheetVia};

use openrdson_core::mesh::Mesh;
use openrdson_solver::{assemble_conductance_by, effective_resistance_from_power, solve_dirichlet};
use std::collections::BTreeMap;

/// Resistance between two node sets with `ΔV = 1` and all other nodes floating.
pub fn two_terminal_resistance(
    mesh: &Mesh,
    sigma_of: &dyn Fn(u32) -> f64,
    a_nodes: &[u32],
    b_nodes: &[u32],
    tol: f64,
    max_iter: usize,
) -> Result<f64, String> {
    let (a, _diags) = assemble_conductance_by(mesh, sigma_of);
    let mut fixed = BTreeMap::new();
    for &n in a_nodes {
        fixed.insert(n, 1.0);
    }
    for &n in b_nodes {
        fixed.insert(n, 0.0);
    }
    if fixed.is_empty() {
        return Err("no terminal nodes fixed".into());
    }
    let u = solve_dirichlet(&a, &fixed, tol, max_iter)?;
    Ok(effective_resistance_from_power(&a, &u))
}

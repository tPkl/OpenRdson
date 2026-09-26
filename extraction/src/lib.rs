//! Resistive network extraction for OpenRDSon.
//!
//! v0 extracts a **reduced terminal-level resistive network** directly from the
//! FEM solution: terminal pairs are driven with `ΔV = 1` while all other
//! terminals float, and the two-terminal resistance is read from the dissipated
//! power (`R = ΔV² / P`).
//!
//! Note: the trilinear hex FEM stiffness is *not* an M-matrix (some off-diagonal
//! couplings are positive), so it must not be naively converted into a
//! two-terminal resistor per off-diagonal entry. The reduced network produced
//! here is always physically valid.

pub mod sheet;
pub mod adaptive;
pub use sheet::{group_vias, map_solution, SheetChannel, SheetNetwork, SheetPolygon, SheetVia};

use openrdson_core::mesh::Mesh;
use openrdson_solver::{
    assemble_conductance_by, effective_resistance_from_power, solve_dirichlet,
};
use std::collections::BTreeMap;

/// A two-terminal resistor between mesh nodes (used by detailed networks).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Resistor {
    pub a: u32,
    pub b: u32,
    pub resistance: f64,
}

/// A reduced terminal-level parasitic network.
#[derive(Debug, Clone, Default)]
pub struct ParasiticNetwork {
    /// Terminal names, in matrix order.
    pub terminals: Vec<String>,
    /// `rmatrix[i][j]` = open-circuit resistance between terminals `i` and `j`.
    pub rmatrix: Vec<Vec<f64>>,
    /// Optional detailed resistors (populated by later milestones).
    pub resistors: Vec<Resistor>,
}

impl ParasiticNetwork {
    pub fn resistance(&self, a: &str, b: &str) -> Option<f64> {
        let ia = self.terminals.iter().position(|t| t == a)?;
        let ib = self.terminals.iter().position(|t| t == b)?;
        Some(self.rmatrix[ia][ib])
    }
}

/// Map each mesh node to a connectivity component.
///
/// A node is assigned a component only when **every** incident element agrees
/// (nodes on interfaces between two nets stay `None`, so they are not treated as
/// terminal nodes).
pub fn node_components(mesh: &Mesh, component_of: &[usize]) -> Vec<Option<usize>> {
    let n = mesh.nodes.len();
    let mut result: Vec<Option<usize>> = vec![None; n];
    let mut agree = vec![true; n];
    for e in &mesh.elements {
        let Some(pi) = e.source_polygon else { continue };
        let comp = match component_of.get(pi).copied() {
            Some(c) if c != usize::MAX => c,
            _ => continue,
        };
        for &node in &e.nodes {
            let idx = node as usize;
            match result[idx] {
                None => result[idx] = Some(comp),
                Some(c) if c == comp => {}
                Some(_) => agree[idx] = false,
            }
        }
    }
    for i in 0..n {
        if !agree[i] {
            result[i] = None;
        }
    }
    result
}

/// Split mesh nodes that are shared between elements of **different
/// connectivity components**, so distinct nets never share a node.
///
/// A coarse mesh can bridge fine gaps (e.g. a 1 µm source/drain gap at a 5 µm
/// cell size), shorting nets that are electrically separate. Splitting the
/// nodes at component boundaries keeps each net's FEM graph connected internally
/// while isolating it from other nets; the only cross-net coupling is then the
/// channel/device edges injected explicitly.
pub fn split_nodes_by_component(mesh: &Mesh, component_of: &[usize]) -> Mesh {
    use std::collections::{BTreeSet, HashMap};
    let n = mesh.nodes.len();
    let comp_of_elem = |e: &openrdson_core::mesh::Element| -> Option<usize> {
        let pi = e.source_polygon?;
        match component_of.get(pi).copied() {
            Some(c) if c != usize::MAX => Some(c),
            _ => None,
        }
    };

    let mut node_comps: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); n];
    for e in &mesh.elements {
        if let Some(c) = comp_of_elem(e) {
            for &node in &e.nodes {
                node_comps[node as usize].insert(c);
            }
        }
    }

    let mut new_nodes = mesh.nodes.clone();
    let mut dup: HashMap<(usize, usize), u32> = HashMap::new();
    let mut elements = Vec::with_capacity(mesh.elements.len());
    for e in &mesh.elements {
        let comp = comp_of_elem(e);
        let mut e = e.clone();
        for node in &mut e.nodes {
            let idx = *node as usize;
            if node_comps[idx].len() > 1 {
                if let Some(c) = comp {
                    let new_idx = *dup.entry((idx, c)).or_insert_with(|| {
                        new_nodes.push(mesh.nodes[idx]);
                        (new_nodes.len() - 1) as u32
                    });
                    *node = new_idx;
                }
            }
        }
        elements.push(e);
    }
    Mesh {
        nodes: new_nodes,
        elements,
        sets: Default::default(),
    }
}

/// Nodes belonging to a given connectivity component.
pub fn terminal_nodes(node_component: &[Option<usize>], component: usize) -> Vec<u32> {
    node_component
        .iter()
        .enumerate()
        .filter(|(_, c)| **c == Some(component))
        .map(|(i, _)| i as u32)
        .collect()
}

/// Resistance with `ΔV = 1` across `fixed` nodes, plus extra device edges
/// (e.g. channel conductances) injected into the network.
pub fn resistance_with_edges(
    mesh: &Mesh,
    sigma_of: &dyn Fn(u32) -> f64,
    fixed: &[(u32, f64)],
    extra_edges: &[(u32, u32, f64)],
    tol: f64,
    max_iter: usize,
) -> Result<f64, String> {
    let (mut a, _diags) = assemble_conductance_by(mesh, sigma_of);
    for &(i, j, g) in extra_edges {
        a.add_edge(i, j, g);
    }
    a.finalize();
    let mut fixed_map = BTreeMap::new();
    for &(n, v) in fixed {
        fixed_map.insert(n, v);
    }
    if fixed_map.is_empty() {
        return Err("no fixed nodes".into());
    }
    let u = solve_dirichlet(&a, &fixed_map, tol, max_iter)?;
    Ok(effective_resistance_from_power(&a, &u))
}

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

/// Extract the reduced terminal resistance matrix.
///
/// For each ordered terminal pair `(i, j)`, terminal `i` is driven to 1 V,
/// terminal `j` to 0 V and every other terminal is left floating; no current can
/// enter a floating terminal, so `R_ij = ΔV² / P = 1 / P`. The matrix is
/// symmetric with a zero diagonal.
pub fn terminal_network(
    mesh: &Mesh,
    sigma_of: &dyn Fn(u32) -> f64,
    terminals: &[(String, Vec<u32>)],
    tol: f64,
    max_iter: usize,
) -> Result<ParasiticNetwork, String> {
    let (a, _diags) = assemble_conductance_by(mesh, sigma_of);
    let t = terminals.len();
    let mut rmatrix = vec![vec![0.0f64; t]; t];
    for i in 0..t {
        for j in (i + 1)..t {
            let mut fixed = BTreeMap::new();
            for &n in &terminals[i].1 {
                fixed.insert(n, 1.0);
            }
            for &n in &terminals[j].1 {
                fixed.insert(n, 0.0);
            }
            let u = solve_dirichlet(&a, &fixed, tol, max_iter)?;
            let r = effective_resistance_from_power(&a, &u);
            rmatrix[i][j] = r;
            rmatrix[j][i] = r;
        }
    }
    Ok(ParasiticNetwork {
        terminals: terminals.iter().map(|(n, _)| n.clone()).collect(),
        rmatrix,
        resistors: Vec::new(),
    })
}

//! Enforces the built-in analytic acceptance checks (`R = ρL/A` on a bar and
//! series-material addition) under `cargo test` / CI.
//!
//! These are the generic, database-independent correctness anchors for the
//! solver + extraction spine; they were previously only reachable by calling
//! `run_acceptance()` manually.

#[test]
fn acceptance_checks_pass() {
    let checks = openrdson_validation::run_acceptance();
    assert!(!checks.is_empty(), "no acceptance checks were produced");
    let failed: Vec<&str> = checks
        .iter()
        .filter(|c| !c.passed)
        .map(|c| c.name.as_str())
        .collect();
    assert!(failed.is_empty(), "acceptance checks failed: {failed:?}");
}

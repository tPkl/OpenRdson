//! Parse diagnostics for imported files.
//!
//! Every text parser returns a `Vec<Diag>` alongside its data so callers can
//! report syntax/semantic problems with the source line, and the `read_*` file
//! wrappers log them as warnings. This keeps import problems visible to the
//! user instead of silently dropping malformed input.

use std::fmt;

/// A problem found while parsing an imported file.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Diag {
    /// 1-based source line (0 = whole file / not line-specific).
    pub line: usize,
    pub message: String,
}

impl Diag {
    /// A diagnostic attached to a 1-based source line.
    pub fn at(line: usize, message: impl Into<String>) -> Self {
        Self {
            line,
            message: message.into(),
        }
    }

    /// A file-level diagnostic (no specific line).
    pub fn file(message: impl Into<String>) -> Self {
        Self::at(0, message)
    }
}

impl fmt::Display for Diag {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.line == 0 {
            write!(f, "{}", self.message)
        } else {
            write!(f, "line {}: {}", self.line, self.message)
        }
    }
}

/// Log every diagnostic as a warning, prefixed with `source` (usually the file
/// name). Returns the number logged, so a caller can summarise if it wants.
///
/// The same `(source, line, message)` is only reported once per process: files
/// are read by several pipeline stages, and repeating the same import problem
/// adds noise without new information.
pub fn log_diags(source: &str, diags: &[Diag]) -> usize {
    use std::collections::HashSet;
    use std::hash::{Hash, Hasher};
    use std::sync::{Mutex, OnceLock};

    static SEEN: OnceLock<Mutex<HashSet<u64>>> = OnceLock::new();
    let seen = SEEN.get_or_init(|| Mutex::new(HashSet::new()));

    let mut logged = 0;
    for d in diags {
        let mut h = std::collections::hash_map::DefaultHasher::new();
        source.hash(&mut h);
        d.hash(&mut h);
        let key = h.finish();
        let first = seen.lock().map(|mut s| s.insert(key)).unwrap_or(true);
        if first {
            crate::log_warn!("{source}: {d}");
            logged += 1;
        }
    }
    logged
}

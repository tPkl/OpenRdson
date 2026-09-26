//! Lightweight, dependency-free logging with verbosity levels.
//!
//! Default level is `INFO`, so pipeline progress is visible out of the box;
//! `--verbose` raises it to `DEBUG` and `--trace` to `TRACE`, `--quiet` lowers
//! it to `WARN`. Messages go to stderr so stdout stays clean for machine output.

use std::sync::atomic::{AtomicU8, Ordering};
use std::sync::{Mutex, OnceLock};

pub const ERROR: u8 = 0;
pub const WARN: u8 = 1;
pub const INFO: u8 = 2;
pub const DEBUG: u8 = 3;
pub const TRACE: u8 = 4;

static LEVEL: AtomicU8 = AtomicU8::new(INFO);

/// Current log section tag (e.g. `sheet-rds`), rendered as `[INFO] [sheet-rds]`.
static SECTION: OnceLock<Mutex<Option<String>>> = OnceLock::new();

fn section() -> &'static Mutex<Option<String>> {
    SECTION.get_or_init(|| Mutex::new(None))
}

/// Set (or clear, with `None`) the current section tag. Every log line emitted
/// while a section is set carries a `[section]` prefix.
pub fn set_section(name: Option<&str>) {
    if let Ok(mut g) = section().lock() {
        *g = name.map(str::to_string);
    }
}

/// RAII guard that sets a log section and clears it on drop. Bind it with a
/// `let _section = SectionGuard::new("name");` in a command handler so that
/// every log line it emits is tagged with that section.
pub struct SectionGuard;

impl SectionGuard {
    pub fn new(name: &str) -> Self {
        set_section(Some(name));
        SectionGuard
    }
}

impl Drop for SectionGuard {
    fn drop(&mut self) {
        set_section(None);
    }
}

pub fn set_level(level: u8) {
    LEVEL.store(level.min(TRACE), Ordering::Relaxed);
}

pub fn level() -> u8 {
    LEVEL.load(Ordering::Relaxed)
}

pub fn enabled(level: u8) -> bool {
    level <= self::level()
}

pub fn tag(level: u8) -> &'static str {
    match level {
        0 => "ERROR",
        1 => "WARN",
        2 => "INFO",
        3 => "DEBUG",
        _ => "TRACE",
    }
}

pub fn log(level: u8, args: std::fmt::Arguments<'_>) {
    if enabled(level) {
        let sec = section().lock().ok().and_then(|g| g.clone());
        match sec {
            Some(s) => eprintln!("[{}] [{}] {}", tag(level), s, args),
            None => eprintln!("[{}] {}", tag(level), args),
        }
    }
}

#[macro_export]
macro_rules! log_error {
    ($($arg:tt)*) => { $crate::log::log($crate::log::ERROR, format_args!($($arg)*)) };
}
#[macro_export]
macro_rules! log_warn {
    ($($arg:tt)*) => { $crate::log::log($crate::log::WARN, format_args!($($arg)*)) };
}
#[macro_export]
macro_rules! log_info {
    ($($arg:tt)*) => { $crate::log::log($crate::log::INFO, format_args!($($arg)*)) };
}
#[macro_export]
macro_rules! log_debug {
    ($($arg:tt)*) => { $crate::log::log($crate::log::DEBUG, format_args!($($arg)*)) };
}
#[macro_export]
macro_rules! log_trace {
    ($($arg:tt)*) => { $crate::log::log($crate::log::TRACE, format_args!($($arg)*)) };
}

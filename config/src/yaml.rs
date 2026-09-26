//! A minimal, dependency-free YAML subset parser.
//!
//! Supports what a config file needs:
//!   - nested block mappings (`key: value`, indentation-driven),
//!   - block sequences (`- item`),
//!   - scalars: null / bool / int / float / string,
//!   - single- and double-quoted strings,
//!   - full-line and trailing `#` comments.
//!
//! It deliberately does NOT support anchors, aliases, flow collections
//! (`[a, b]`, `{a: b}`), multi-line scalars, or tags.

use std::collections::BTreeMap;

/// A parsed YAML node.
#[derive(Debug, Clone, PartialEq)]
pub enum Value {
    Null,
    Bool(bool),
    Int(i64),
    Float(f64),
    Str(String),
    List(Vec<Value>),
    Map(BTreeMap<String, Value>),
}

impl Value {
    /// Look up a key in a mapping (returns `None` for non-mappings).
    pub fn get(&self, key: &str) -> Option<&Value> {
        match self {
            Value::Map(m) => m.get(key),
            _ => None,
        }
    }

    pub fn as_str(&self) -> Option<&str> {
        match self {
            Value::Str(s) => Some(s),
            _ => None,
        }
    }

    pub fn as_f64(&self) -> Option<f64> {
        match self {
            Value::Int(i) => Some(*i as f64),
            Value::Float(f) => Some(*f),
            Value::Str(s) => s.parse().ok(),
            _ => None,
        }
    }

    pub fn as_i64(&self) -> Option<i64> {
        match self {
            Value::Int(i) => Some(*i),
            Value::Float(f) => Some(*f as i64),
            Value::Str(s) => s.parse().ok(),
            _ => None,
        }
    }

    pub fn as_bool(&self) -> Option<bool> {
        match self {
            Value::Bool(b) => Some(*b),
            Value::Str(s) => match s.to_ascii_lowercase().as_str() {
                "true" | "yes" | "on" => Some(true),
                "false" | "no" | "off" => Some(false),
                _ => None,
            },
            _ => None,
        }
    }

    pub fn as_list(&self) -> Option<&[Value]> {
        match self {
            Value::List(v) => Some(v),
            _ => None,
        }
    }
}

struct Line {
    indent: usize,
    text: String,
}

/// Parse a YAML document (subset) into a [`Value`].
pub fn parse(text: &str) -> Result<Value, String> {
    let mut lines: Vec<Line> = Vec::new();
    for raw in text.lines() {
        let stripped = strip_comment(raw);
        let trimmed = stripped.trim_end();
        if trimmed.trim().is_empty() {
            continue;
        }
        let indent = trimmed.len() - trimmed.trim_start().len();
        lines.push(Line {
            indent,
            text: trimmed.trim_start().to_string(),
        });
    }
    if lines.is_empty() {
        return Ok(Value::Map(BTreeMap::new()));
    }
    let (value, _) = parse_block(&lines, 0, lines[0].indent)?;
    Ok(value)
}

fn parse_block(lines: &[Line], mut i: usize, indent: usize) -> Result<(Value, usize), String> {
    if is_seq_line(&lines[i].text) {
        let mut items = Vec::new();
        while i < lines.len() && lines[i].indent == indent && is_seq_line(&lines[i].text) {
            let text = lines[i].text.clone();
            let rest = text[1..].trim_start().to_string();
            // Column where the item's content begins (after "- ").
            let content_indent = indent + 1 + (text[1..].len() - rest.len());
            let mut sub: Vec<Line> = Vec::new();
            if !rest.is_empty() {
                sub.push(Line {
                    indent: content_indent,
                    text: rest,
                });
            }
            i += 1;
            while i < lines.len() && lines[i].indent > indent {
                sub.push(Line {
                    indent: lines[i].indent,
                    text: lines[i].text.clone(),
                });
                i += 1;
            }
            if sub.is_empty() {
                items.push(Value::Null);
            } else {
                let (v, _) = parse_block(&sub, 0, sub[0].indent)?;
                items.push(v);
            }
        }
        Ok((Value::List(items), i))
    } else if find_colon(&lines[i].text).is_none() {
        // Bare scalar (single line).
        Ok((scalar(&lines[i].text), i + 1))
    } else {
        let mut map = BTreeMap::new();
        while i < lines.len() && lines[i].indent == indent && !is_seq_line(&lines[i].text) {
            let text = lines[i].text.clone();
            let Some(colon) = find_colon(&text) else {
                return Err(format!("expected `key:` but found `{text}`"));
            };
            let key = text[..colon].trim().trim_matches('"').trim_matches('\'').to_string();
            let rest = text[colon + 1..].trim().to_string();
            i += 1;
            if rest.is_empty() {
                if i < lines.len() && lines[i].indent > indent {
                    let (v, ni) = parse_block(lines, i, lines[i].indent)?;
                    map.insert(key, v);
                    i = ni;
                } else {
                    map.insert(key, Value::Null);
                }
            } else {
                map.insert(key, scalar(&rest));
            }
        }
        Ok((Value::Map(map), i))
    }
}

fn is_seq_line(text: &str) -> bool {
    text == "-" || text.starts_with("- ")
}

/// Find the first `:` that is not inside quotes.
fn find_colon(text: &str) -> Option<usize> {
    let mut in_s = false;
    let mut in_d = false;
    for (i, c) in text.char_indices() {
        match c {
            '\'' if !in_d => in_s = !in_s,
            '"' if !in_s => in_d = !in_d,
            ':' if !in_s && !in_d => return Some(i),
            _ => {}
        }
    }
    None
}

/// Remove a trailing `#` comment (outside quotes).
fn strip_comment(s: &str) -> &str {
    let bytes = s.as_bytes();
    let mut in_s = false;
    let mut in_d = false;
    let mut i = 0;
    while i < bytes.len() {
        let c = bytes[i] as char;
        match c {
            '\'' if !in_d => in_s = !in_s,
            '"' if !in_s => in_d = !in_d,
            '#' if !in_s && !in_d => {
                if i == 0 || (bytes[i - 1] as char).is_whitespace() {
                    return &s[..i];
                }
            }
            _ => {}
        }
        i += 1;
    }
    s
}

fn scalar(s: &str) -> Value {
    let s = s.trim();
    if s.is_empty() {
        return Value::Null;
    }
    if s.len() >= 2
        && ((s.starts_with('"') && s.ends_with('"')) || (s.starts_with('\'') && s.ends_with('\'')))
    {
        return Value::Str(s[1..s.len() - 1].to_string());
    }
    match s {
        "null" | "~" => Value::Null,
        "true" | "yes" | "on" => Value::Bool(true),
        "false" | "no" | "off" => Value::Bool(false),
        _ => {
            if let Ok(i) = s.parse::<i64>() {
                Value::Int(i)
            } else if let Ok(f) = s.parse::<f64>() {
                Value::Float(f)
            } else {
                Value::Str(s.to_string())
            }
        }
    }
}

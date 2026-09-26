//! Rectilinear polygon boolean operations via scanline interval arithmetic.
//!
//! These operate on even-odd interval sets per horizontal band. They are exact
//! for rectilinear (axis-aligned) rings, including keyhole rings with holes, and
//! give a staircase approximation for rings with diagonal edges.

use crate::geometry::{Bbox, Point};

/// An axis-aligned rectangle.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Rect {
    pub x0: f64,
    pub y0: f64,
    pub x1: f64,
    pub y1: f64,
}

impl Rect {
    pub fn area(&self) -> f64 {
        (self.x1 - self.x0).max(0.0) * (self.y1 - self.y0).max(0.0)
    }

    /// The rectangle as a closed 5-point ring.
    pub fn to_ring(self) -> Vec<Point> {
        vec![
            Point::new(self.x0, self.y0),
            Point::new(self.x1, self.y0),
            Point::new(self.x1, self.y1),
            Point::new(self.x0, self.y1),
            Point::new(self.x0, self.y0),
        ]
    }

    pub fn bbox(&self) -> Bbox {
        Bbox {
            min_x: self.x0,
            min_y: self.y0,
            max_x: self.x1,
            max_y: self.y1,
        }
    }
}

/// Even-odd x-intervals of `ring` on the horizontal line `y`.
pub fn x_intervals_at(ring: &[Point], y: f64) -> Vec<(f64, f64)> {
    let n = ring.len();
    let mut xs: Vec<f64> = Vec::new();
    for i in 0..n {
        let a = ring[i];
        let b = ring[(i + 1) % n];
        if (a.y <= y && b.y > y) || (b.y <= y && a.y > y) {
            let x = a.x + (y - a.y) / (b.y - a.y) * (b.x - a.x);
            xs.push(x);
        }
    }
    xs.sort_by(|p, q| p.partial_cmp(q).unwrap());
    xs.chunks_exact(2).map(|c| (c[0], c[1])).collect()
}

fn bands(a: &[Point], b: &[Point]) -> Vec<f64> {
    let mut ys: Vec<f64> = a.iter().chain(b.iter()).map(|p| p.y).collect();
    ys.sort_by(|p, q| p.partial_cmp(q).unwrap());
    ys.dedup_by(|p, q| (*p - *q).abs() < 1e-15);
    ys
}

/// `a \ b` as a set of rectangles.
pub fn polygon_minus_polygon(a: &[Point], b: &[Point]) -> Vec<Rect> {
    let mut out = Vec::new();
    if a.len() < 3 {
        return out;
    }
    let ys = bands(a, b);
    for k in 0..ys.len().saturating_sub(1) {
        let (y0, y1) = (ys[k], ys[k + 1]);
        if y1 <= y0 {
            continue;
        }
        let ym = 0.5 * (y0 + y1);
        let ia = x_intervals_at(a, ym);
        let ib = if b.len() >= 3 {
            x_intervals_at(b, ym)
        } else {
            Vec::new()
        };
        for &(a0, a1) in &ia {
            let mut segs = vec![(a0, a1)];
            for &(b0, b1) in &ib {
                let mut next = Vec::new();
                for &(s0, s1) in &segs {
                    if b1 <= s0 || b0 >= s1 {
                        next.push((s0, s1));
                        continue;
                    }
                    if b0 > s0 {
                        next.push((s0, b0));
                    }
                    if b1 < s1 {
                        next.push((b1, s1));
                    }
                }
                segs = next;
            }
            for (s0, s1) in segs {
                if s1 - s0 > 1e-15 {
                    out.push(Rect {
                        x0: s0,
                        y0,
                        x1: s1,
                        y1,
                    });
                }
            }
        }
    }
    out
}

/// `a ∩ b` as a set of rectangles.
pub fn polygon_intersection_rects(a: &[Point], b: &[Point]) -> Vec<Rect> {
    let mut out = Vec::new();
    if a.len() < 3 || b.len() < 3 {
        return out;
    }
    let ys = bands(a, b);
    for k in 0..ys.len().saturating_sub(1) {
        let (y0, y1) = (ys[k], ys[k + 1]);
        if y1 <= y0 {
            continue;
        }
        let ym = 0.5 * (y0 + y1);
        for &(a0, a1) in &x_intervals_at(a, ym) {
            for &(b0, b1) in &x_intervals_at(b, ym) {
                let x0 = a0.max(b0);
                let x1 = a1.min(b1);
                if x1 - x0 > 1e-15 {
                    out.push(Rect {
                        x0,
                        y0,
                        x1,
                        y1,
                    });
                }
            }
        }
    }
    out
}

/// Total area of a rectangle set.
pub fn rects_area(rects: &[Rect]) -> f64 {
    rects.iter().map(|r| r.area()).sum()
}

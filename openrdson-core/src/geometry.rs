//! Basic geometric primitives shared across the pipeline.

/// A 2D point in **meters** (SI) unless otherwise documented.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point {
    pub x: f64,
    pub y: f64,
}

impl Point {
    pub const fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }
}

/// A GDSII layer/datatype pair. `Ord` for stable, deterministic map iteration.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct LayerKey {
    pub layer: i16,
    pub datatype: i16,
}

impl LayerKey {
    pub const fn new(layer: i16, datatype: i16) -> Self {
        Self { layer, datatype }
    }
}

impl std::fmt::Display for LayerKey {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}:{}", self.layer, self.datatype)
    }
}

/// An axis-aligned bounding box.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Bbox {
    pub min_x: f64,
    pub min_y: f64,
    pub max_x: f64,
    pub max_y: f64,
}

impl Default for Bbox {
    fn default() -> Self {
        Self::empty()
    }
}

impl Bbox {
    /// An "inverted" box that absorbs the first point extended into it.
    pub fn empty() -> Self {
        Self {
            min_x: f64::INFINITY,
            min_y: f64::INFINITY,
            max_x: f64::NEG_INFINITY,
            max_y: f64::NEG_INFINITY,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.min_x > self.max_x || self.min_y > self.max_y
    }

    pub fn from_points(points: &[Point]) -> Self {
        let mut b = Self::empty();
        for p in points {
            b.extend_point(*p);
        }
        b
    }

    pub fn extend_point(&mut self, p: Point) {
        if p.x < self.min_x {
            self.min_x = p.x;
        }
        if p.y < self.min_y {
            self.min_y = p.y;
        }
        if p.x > self.max_x {
            self.max_x = p.x;
        }
        if p.y > self.max_y {
            self.max_y = p.y;
        }
    }

    pub fn extend(&mut self, other: &Bbox) {
        if other.is_empty() {
            return;
        }
        self.extend_point(Point::new(other.min_x, other.min_y));
        self.extend_point(Point::new(other.max_x, other.max_y));
    }

    pub fn union(&self, other: &Bbox) -> Bbox {
        let mut b = *self;
        b.extend(other);
        b
    }

    pub fn width(&self) -> f64 {
        if self.is_empty() {
            0.0
        } else {
            self.max_x - self.min_x
        }
    }

    pub fn height(&self) -> f64 {
        if self.is_empty() {
            0.0
        } else {
            self.max_y - self.min_y
        }
    }

    pub fn center(&self) -> Point {
        Point::new(
            (self.min_x + self.max_x) / 2.0,
            (self.min_y + self.max_y) / 2.0,
        )
    }
}

/// A 3D solid primitive. v0 supports the two shapes the 2.5D path needs.
#[derive(Debug, Clone, PartialEq)]
pub enum Solid {
    /// Axis-aligned box: footprint `[min_x..max_x] x [min_y..max_y]`, extruded.
    Box {
        min_x: f64,
        min_y: f64,
        max_x: f64,
        max_y: f64,
        z_bottom: f64,
        z_top: f64,
    },
    /// A general prism: arbitrary footprint extruded between two z-planes.
    Prism {
        footprint: Vec<Point>,
        z_bottom: f64,
        z_top: f64,
    },
}

impl Solid {
    pub fn volume(&self) -> f64 {
        match self {
            Solid::Box {
                min_x,
                min_y,
                max_x,
                max_y,
                z_bottom,
                z_top,
            } => ((max_x - min_x) * (max_y - min_y) * (z_top - z_bottom)).abs(),
            Solid::Prism {
                footprint,
                z_bottom,
                z_top,
            } => polygon_area(footprint) * (z_top - z_bottom).abs(),
        }
    }

    pub fn z_bottom(&self) -> f64 {
        match self {
            Solid::Box { z_bottom, .. } => *z_bottom,
            Solid::Prism { z_bottom, .. } => *z_bottom,
        }
    }

    pub fn z_top(&self) -> f64 {
        match self {
            Solid::Box { z_top, .. } => *z_top,
            Solid::Prism { z_top, .. } => *z_top,
        }
    }

    pub fn footprint_bbox(&self) -> Bbox {
        match self {
            Solid::Box {
                min_x,
                min_y,
                max_x,
                max_y,
                ..
            } => Bbox {
                min_x: *min_x,
                min_y: *min_y,
                max_x: *max_x,
                max_y: *max_y,
            },
            Solid::Prism { footprint, .. } => Bbox::from_points(footprint),
        }
    }
}

/// A tagged 3D region ready for meshing (the `SolidModel` element).
#[derive(Debug, Clone, PartialEq)]
pub struct SolidRegion {
    pub material_id: u32,
    pub net: Option<String>,
    pub layer: Option<String>,
    pub device_ref: Option<String>,
    /// Index of the source layout polygon, used to map mesh elements back to
    /// connectivity components.
    pub source_polygon: Option<usize>,
    pub solid: Solid,
}

impl SolidRegion {
    pub fn volume(&self) -> f64 {
        self.solid.volume()
    }
}

/// The tagged 3D solid model handed from geometry to meshing (the `SolidModel`
/// contract in `goal.md` §5).
#[derive(Debug, Clone, Default)]
pub struct SolidModel {
    pub regions: Vec<SolidRegion>,
    /// `material_names[id - 1]` is the physical layer name for `material_id = id`.
    /// `material_id = 0` means "unassigned".
    pub material_names: Vec<String>,
}

impl SolidModel {
    pub fn total_volume(&self) -> f64 {
        self.regions.iter().map(|r| r.volume()).sum()
    }

    pub fn material_name(&self, id: u32) -> Option<&str> {
        if id == 0 {
            None
        } else {
            self.material_names.get((id - 1) as usize).map(|s| s.as_str())
        }
    }
}

/// True if the polygon is an axis-aligned rectangle.
///
/// Requires every vertex to lie on a corner of the bounding box AND all four
/// corners to be present (so a right triangle spanning three corners is
/// rejected).
pub fn is_axis_aligned_rectangle(points: &[Point]) -> bool {
    let b = Bbox::from_points(points);
    if b.width() <= 0.0 || b.height() <= 0.0 {
        return false;
    }
    let eps = 1e-9 * b.width().max(b.height());
    let corners = [
        (b.min_x, b.min_y),
        (b.max_x, b.min_y),
        (b.max_x, b.max_y),
        (b.min_x, b.max_y),
    ];
    let mut present = [false; 4];
    for p in points {
        let mut on_corner = false;
        for (ci, (cx, cy)) in corners.iter().enumerate() {
            if (p.x - cx).abs() <= eps && (p.y - cy).abs() <= eps {
                present[ci] = true;
                on_corner = true;
            }
        }
        if !on_corner {
            return false;
        }
    }
    present.iter().all(|&b| b)
}

/// Signed area of a polygon via the shoelace formula. Positive = CCW.
pub fn polygon_signed_area(points: &[Point]) -> f64 {
    if points.len() < 3 {
        return 0.0;
    }
    let mut acc = 0.0;
    for i in 0..points.len() {
        let a = points[i];
        let b = points[(i + 1) % points.len()];
        acc += a.x * b.y - b.x * a.y;
    }
    acc / 2.0
}

pub fn polygon_area(points: &[Point]) -> f64 {
    polygon_signed_area(points).abs()
}

/// Whether the polygon winds counter-clockwise (positive signed area).
pub fn is_ccw(points: &[Point]) -> bool {
    polygon_signed_area(points) > 0.0
}

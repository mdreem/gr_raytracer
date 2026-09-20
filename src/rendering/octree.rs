use crate::rendering::star_catalog::Star;
use nalgebra::Vector3;
use nalgebra::{DimMax, DimMin};

pub struct Octree {
    root: Node,
}

impl Octree {
    pub fn new(stars: Vec<Star>, depth: usize) -> Self {
        let mut root = Node {
            bounds: AABB::new(Vector3::new(-1.0, -1.0, -1.0), Vector3::new(1.0, 1.0, 1.0)),
            children: None,
            stars: None,
        };

        for star in stars {
            root.add_star(star, depth);
        }

        Self { root }
    }
}

struct Node {
    bounds: AABB,
    pub children: Option<[Box<Node>; 8]>,
    pub stars: Option<Vec<Star>>,
}

impl Node {
    fn contains_point(&self, point: &Vector3<f64>) -> bool {
        point.x >= self.bounds.min.x
            && point.x < self.bounds.max.x
            && point.y >= self.bounds.min.y
            && point.y < self.bounds.max.y
            && point.z >= self.bounds.min.z
            && point.z < self.bounds.max.z
    }

    fn get_child_index(&self, point: &Vector3<f64>) -> usize {
        let mid_x = self.bounds.center.x;
        let mid_y = self.bounds.center.y;
        let mid_z = self.bounds.center.z;

        let mut index = 0;
        if point.x >= mid_x {
            index |= 1;
        }
        if point.y >= mid_y {
            index |= 2;
        }
        if point.z >= mid_z {
            index |= 4;
        }
        index
    }

    fn get_child(&self, point: &Vector3<f64>) -> Option<&Node> {
        if !self.contains_point(point) {
            return None;
        }
        if let Some(children) = &self.children {
            Some(children[self.get_child_index(point)].as_ref())
        } else {
            None
        }
    }

    fn child_bounds(&self, i: usize) -> AABB {
        let c = self.bounds.center;
        let min = Vector3::new(
            if i & 1 == 0 { self.bounds.min.x } else { c.x },
            if i & 2 == 0 { self.bounds.min.y } else { c.y },
            if i & 4 == 0 { self.bounds.min.z } else { c.z },
        );
        let max = Vector3::new(
            if i & 1 == 0 { c.x } else { self.bounds.max.x },
            if i & 2 == 0 { c.y } else { self.bounds.max.y },
            if i & 4 == 0 { c.z } else { self.bounds.max.z },
        );
        AABB::new(min, max)
    }

    fn subdivide(&mut self) {
        let children = std::array::from_fn(|i| {
            Box::new(Node {
                bounds: self.child_bounds(i),
                children: None,
                stars: None,
            })
        });
        self.children = Some(children);
    }

    fn add_star(&mut self, star: Star, max_depth: usize) {
        if max_depth == 0 {
            if let Some(stars) = &mut self.stars {
                stars.push(star);
            } else {
                self.stars = Some(vec![star]);
            }
            return;
        }

        if self.children.is_none() {
            self.subdivide();
        }

        let index = self.get_child_index(&star.direction);
        if let Some(children) = &mut self.children {
            children[index].add_star(star, max_depth - 1);
        }
    }
}

struct AABB {
    min: Vector3<f64>,
    max: Vector3<f64>,
    pub center: Vector3<f64>,
    pub extents: Vector3<f64>,
}

impl AABB {
    fn new(min: Vector3<f64>, max: Vector3<f64>) -> Self {
        debug_assert!(
            min.x <= max.x && min.y <= max.y && min.z <= max.z,
            "AABB min must be <= max"
        );
        let center = (min + max) / 2.0;
        let extents = (max - min) / 2.0;

        Self {
            min,
            max,
            center,
            extents,
        }
    }
}

struct Triangle {
    vertices: [Vector3<f64>; 3],
}

// https://gdbooks.gitbooks.io/3dcollisions/content/Chapter4/aabb-triangle.html
fn intersects_aabb(aabb: &AABB, triangle: &Triangle) -> bool {
    let v0 = triangle.vertices[0] - aabb.center;
    let v1 = triangle.vertices[1] - aabb.center;
    let v2 = triangle.vertices[2] - aabb.center;
    let e = aabb.extents;

    let f0 = v1 - v0;
    let f1 = v2 - v1;
    let f2 = v0 - v2;

    let u0 = Vector3::new(1.0, 0.0, 0.0);
    let u1 = Vector3::new(0.0, 1.0, 0.0);
    let u2 = Vector3::new(0.0, 0.0, 1.0);

    let axis_u0_f0 = u0.cross(&f0);
    let axis_u0_f1 = u0.cross(&f1);
    let axis_u0_f2 = u0.cross(&f2);

    let axis_u1_f0 = u1.cross(&f0);
    let axis_u1_f1 = u1.cross(&f1);
    let axis_u1_f2 = u1.cross(&f2);

    let axis_u2_f0 = u2.cross(&f0);
    let axis_u2_f1 = u2.cross(&f1);
    let axis_u2_f2 = u2.cross(&f2);

    let triangle_normal = f0.cross(&f1);
    let axes = [
        axis_u0_f0,
        axis_u0_f1,
        axis_u0_f2,
        axis_u1_f0,
        axis_u1_f1,
        axis_u1_f2,
        axis_u2_f0,
        axis_u2_f1,
        axis_u2_f2,
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.0, 1.0, 0.0),
        Vector3::new(0.0, 0.0, 1.0),
        triangle_normal,
    ];

    // SAT tests for the edges of the triangle and the axes of the AABB
    for axis in axes.iter() {
        let p0 = v0.dot(axis);
        let p1 = v1.dot(axis);
        let p2 = v2.dot(axis);

        //  e.x * axis.dot(u0).abs() =  e.x * axis.x.abs()
        let r = e.x * axis.x.abs() + e.y * axis.y.abs() + e.z * axis.z.abs();

        let p_min = p0.min(p1).min(p2);
        let p_max = p0.max(p1).max(p2);

        if (-p_max).max(p_min) > r {
            return false;
        }
    }
    true
}

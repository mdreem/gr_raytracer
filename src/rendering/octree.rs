use nalgebra::{DimMax, DimMin};
use noise::Vector3;

pub struct Octree {}

struct Node {
    bound_x: [f64; 2],
    bound_y: [f64; 2],
    bound_z: [f64; 2],
    pub children: Option<[Box<Node>; 8]>,
}

impl Node {
    fn contains_point(&self, point: &Vector3<f64>) -> bool {
        point.x >= self.bound_x[0]
            && point.x < self.bound_x[1]
            && point.y >= self.bound_y[0]
            && point.y < self.bound_y[1]
            && point.z >= self.bound_z[0]
            && point.z < self.bound_z[1]
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

    let axis_u0_f0 = u0.cross(f0);
    let axis_u0_f1 = u0.cross(f1);
    let axis_u0_f2 = u0.cross(f2);

    let axis_u1_f0 = u1.cross(f0);
    let axis_u1_f1 = u1.cross(f1);
    let axis_u1_f2 = u1.cross(f2);

    let axis_u2_f0 = u2.cross(f0);
    let axis_u2_f1 = u2.cross(f1);
    let axis_u2_f2 = u2.cross(f2);

    let triangle_normal = f0.cross(f1);
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

use crate::rendering::star_catalog::Star;
use nalgebra::Vector3;
use nalgebra::{DimMax, DimMin};

pub struct Octree {
    root: Node,
}

impl Octree {
    pub fn new(stars: Vec<Star>) -> Self {
        let mut root = Node {
            bounds: AABB::new(Vector3::new(-1.0, -1.0, -1.0), Vector3::new(1.0, 1.0, 1.0)),
            children: None,
            stars: None,
        };

        for star in stars {
            root.add_star(star);
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
    const LEAF_CAPACITY: usize = 10;
    // The minimum half-extent of a node's bounding box before it stops subdividing.
    // Ensures there is not infinite subdivision of the octree when stars are very close together.
    const MIN_HALF_EXTENT: f64 = 1e-5;

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

    fn add_star(&mut self, star: Star) {
        if self.children.is_some() {
            let index = self.get_child_index(&star.direction);
            self.children.as_mut().unwrap()[index].add_star(star);
            return;
        }

        self.stars.get_or_insert_with(Vec::new).push(star);

        if self.stars.as_ref().map_or(0, Vec::len) > Self::LEAF_CAPACITY
            && self.bounds.extents.x > Self::MIN_HALF_EXTENT
        {
            self.subdivide();
            let stars = self.stars.take().unwrap();
            for s in stars {
                let index = self.get_child_index(&s.direction);
                self.children.as_mut().unwrap()[index].add_star(s);
            }
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rendering::color::CIETristimulus;

    /// A star at `direction`; every other field is a placeholder, since the
    /// octree only ever looks at `direction`.
    fn star_at(direction: Vector3<f64>) -> Star {
        Star {
            source_id: 0,
            ra_deg: 0.0,
            dec_deg: 0.0,
            g_mag: 0.0,
            bp_mag: 0.0,
            rp_mag: 0.0,
            bp_rp: 0.0,
            direction,
            temperature: 5000.0,
            emission_xyz: CIETristimulus::new(0.0, 0.0, 0.0, 1.0),
        }
    }

    fn unit_box() -> AABB {
        AABB::new(Vector3::new(-1.0, -1.0, -1.0), Vector3::new(1.0, 1.0, 1.0))
    }

    fn leaf(bounds: AABB) -> Node {
        Node { bounds, children: None, stars: None }
    }

    fn tri(a: Vector3<f64>, b: Vector3<f64>, c: Vector3<f64>) -> Triangle {
        Triangle { vertices: [a, b, c] }
    }

    /// Total stars stored anywhere in the subtree.
    fn count_stars(node: &Node) -> usize {
        let here = node.stars.as_ref().map_or(0, Vec::len);
        let below = node
            .children
            .as_ref()
            .map_or(0, |c| c.iter().map(|n| count_stars(n)).sum());
        here + below
    }

    // ---- AABB -------------------------------------------------------------

    #[test]
    fn aabb_center_and_extents() {
        let b = AABB::new(Vector3::new(0.0, 0.0, 0.0), Vector3::new(2.0, 4.0, 6.0));
        assert_eq!(b.center, Vector3::new(1.0, 2.0, 3.0));
        assert_eq!(b.extents, Vector3::new(1.0, 2.0, 3.0));
        // extents == max - center
        assert_eq!(b.extents, b.max - b.center);
    }

    // ---- Node geometry ----------------------------------------------------

    #[test]
    fn contains_point_is_min_inclusive_max_exclusive() {
        let n = leaf(unit_box());
        assert!(n.contains_point(&Vector3::new(0.0, 0.0, 0.0)));
        assert!(n.contains_point(&Vector3::new(-1.0, -1.0, -1.0))); // min inclusive
        assert!(!n.contains_point(&Vector3::new(1.0, 0.0, 0.0))); // max exclusive
        assert!(!n.contains_point(&Vector3::new(2.0, 0.0, 0.0)));
    }

    #[test]
    fn child_index_maps_octant_bits() {
        let n = leaf(unit_box());
        assert_eq!(n.get_child_index(&Vector3::new(-0.5, -0.5, -0.5)), 0);
        assert_eq!(n.get_child_index(&Vector3::new(0.5, -0.5, -0.5)), 1); // +x
        assert_eq!(n.get_child_index(&Vector3::new(-0.5, 0.5, -0.5)), 2); // +y
        assert_eq!(n.get_child_index(&Vector3::new(-0.5, -0.5, 0.5)), 4); // +z
        assert_eq!(n.get_child_index(&Vector3::new(0.5, 0.5, 0.5)), 7); // +x+y+z
    }

    #[test]
    fn child_bounds_partition_the_parent() {
        let n = leaf(unit_box());
        let c0 = n.child_bounds(0);
        assert_eq!(c0.min, Vector3::new(-1.0, -1.0, -1.0));
        assert_eq!(c0.max, Vector3::new(0.0, 0.0, 0.0));

        let c7 = n.child_bounds(7);
        assert_eq!(c7.min, Vector3::new(0.0, 0.0, 0.0));
        assert_eq!(c7.max, Vector3::new(1.0, 1.0, 1.0));

        // A child's index is exactly the octant its own center falls in.
        for i in 0..8 {
            assert_eq!(n.get_child_index(&n.child_bounds(i).center), i);
        }
    }

    // ---- SAT triangle/box test -------------------------------------------

    #[test]
    fn triangle_overlapping_box_intersects() {
        let t = tri(
            Vector3::new(-0.5, -0.5, 0.0),
            Vector3::new(0.5, -0.5, 0.0),
            Vector3::new(0.0, 0.5, 0.0),
        );
        assert!(intersects_aabb(&unit_box(), &t));
    }

    #[test]
    fn triangle_fully_inside_intersects() {
        let t = tri(
            Vector3::new(-0.2, -0.2, 0.0),
            Vector3::new(0.2, -0.2, 0.0),
            Vector3::new(0.0, 0.2, 0.1),
        );
        assert!(intersects_aabb(&unit_box(), &t));
    }

    #[test]
    fn triangle_above_box_is_separated() {
        let t = tri(
            Vector3::new(0.0, 0.0, 2.0),
            Vector3::new(1.0, 0.0, 2.0),
            Vector3::new(0.0, 1.0, 2.0),
        );
        assert!(!intersects_aabb(&unit_box(), &t));
    }

    #[test]
    fn triangle_below_box_is_separated() {
        // Regression for the negative-side case: the original test used
        // max(p_max, p_min) and missed a triangle below the box.
        let t = tri(
            Vector3::new(0.0, 0.0, -2.0),
            Vector3::new(1.0, 0.0, -2.0),
            Vector3::new(0.0, 1.0, -2.0),
        );
        assert!(!intersects_aabb(&unit_box(), &t));
    }

    #[test]
    fn triangle_off_to_the_side_is_separated() {
        let t = tri(
            Vector3::new(5.0, 0.0, 0.0),
            Vector3::new(6.0, 0.0, 0.0),
            Vector3::new(5.0, 1.0, 0.0),
        );
        assert!(!intersects_aabb(&unit_box(), &t));
    }

    #[test]
    fn triangle_parallel_touching_top_face_intersects() {
        // Triangle lies in the plane z = 1 (the box's top face). Touching counts
        // as intersecting under SAT (the gap test is strict).
        let t = tri(
            Vector3::new(-0.5, -0.5, 1.0),
            Vector3::new(0.5, -0.5, 1.0),
            Vector3::new(0.0, 0.5, 1.0),
        );
        assert!(intersects_aabb(&unit_box(), &t));
    }

    #[test]
    fn triangle_straddling_a_face_intersects() {
        // One vertex inside the box, the rest outside past +x.
        let t = tri(
            Vector3::new(0.5, 0.0, 0.0),
            Vector3::new(2.0, 0.5, 0.0),
            Vector3::new(2.0, -0.5, 0.0),
        );
        assert!(intersects_aabb(&unit_box(), &t));
    }

    #[test]
    fn triangle_separated_by_edge_cross_axis() {
        // A diagonal box, triangle tucked past a corner so only an edge-edge
        // cross-product axis separates them (not a face normal).
        let t = tri(
            Vector3::new(2.0, 2.0, 0.0),
            Vector3::new(1.4, 2.0, 0.0),
            Vector3::new(2.0, 1.4, 0.0),
        );
        assert!(!intersects_aabb(&unit_box(), &t));
    }

    // ---- Adaptive insertion ----------------------------------------------

    #[test]
    fn at_capacity_stays_a_single_leaf() {
        let stars: Vec<Star> = (0..Node::LEAF_CAPACITY)
            .map(|i| star_at(Vector3::new(0.9 - 0.05 * i as f64, 0.1, 0.1)))
            .collect();
        let tree = Octree::new(stars);
        assert!(tree.root.children.is_none(), "should not subdivide at capacity");
        assert_eq!(count_stars(&tree.root), Node::LEAF_CAPACITY);
    }

    #[test]
    fn over_capacity_subdivides_and_keeps_all_stars() {
        // One star per octant plus extras: exceeds capacity and spreads out.
        let dirs = [
            Vector3::new(1.0, 1.0, 1.0),
            Vector3::new(-1.0, 1.0, 1.0),
            Vector3::new(1.0, -1.0, 1.0),
            Vector3::new(1.0, 1.0, -1.0),
            Vector3::new(-1.0, -1.0, 1.0),
            Vector3::new(-1.0, 1.0, -1.0),
            Vector3::new(1.0, -1.0, -1.0),
            Vector3::new(-1.0, -1.0, -1.0),
            Vector3::new(0.5, 0.5, 0.5),
            Vector3::new(-0.5, -0.5, -0.5),
            Vector3::new(0.3, -0.3, 0.3),
        ];
        assert!(dirs.len() > Node::LEAF_CAPACITY);
        let stars: Vec<Star> = dirs.iter().map(|d| star_at(d.normalize())).collect();
        let n = stars.len();

        let tree = Octree::new(stars);
        assert!(tree.root.children.is_some(), "should subdivide past capacity");
        assert_eq!(count_stars(&tree.root), n, "no stars lost during subdivision");
    }

    #[test]
    fn coincident_directions_terminate_and_keep_all_stars() {
        // Far more identical directions than capacity. Without the MIN_HALF_EXTENT
        // guard this would recurse forever (they never separate); it must instead
        // stop splitting and pile them all into one deep leaf.
        let dir = Vector3::new(0.3, 0.4, 0.5).normalize();
        let stars: Vec<Star> = (0..100).map(|_| star_at(dir)).collect();

        let tree = Octree::new(stars); // must not stack-overflow
        assert_eq!(count_stars(&tree.root), 100);
    }
}

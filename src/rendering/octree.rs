use crate::rendering::star_catalog::Star;
use nalgebra::Vector3;

/// Inward-facing normals of a spherical triangle's three edge planes (each
/// through the origin), oriented by the opposite vertex so "inside" is the
/// positive side. Three points always share a hemisphere, so unlike a quad this
/// is always a well-formed convex region.
fn triangle_planes(t: &[Vector3<f64>; 3]) -> [Vector3<f64>; 3] {
    let mut n = [t[0].cross(&t[1]), t[1].cross(&t[2]), t[2].cross(&t[0])];
    let opposite = [t[2], t[0], t[1]];
    for (ni, opp) in n.iter_mut().zip(opposite) {
        if ni.dot(&opp) < 0.0 {
            *ni = -*ni;
        }
    }
    n
}

fn point_in_triangle(point: &Vector3<f64>, planes: &[Vector3<f64>; 3]) -> bool {
    planes.iter().all(|n| n.dot(point) >= 0.0)
}

/// How an AABB sits relative to the frustum.
enum Class {
    /// No overlap; the whole subtree can be pruned.
    Outside,
    /// Fully contained; every star beneath can be taken without further tests.
    Inside,
    /// Crosses the boundary; recurse or test individual stars.
    Straddle,
}

pub struct Octree {
    root: Node,
    len: usize,
}

impl Octree {
    pub fn new(stars: Vec<Star>) -> Self {
        let mut root = Node {
            bounds: AABB::new(Vector3::new(-1.0, -1.0, -1.0), Vector3::new(1.0, 1.0, 1.0)),
            children: None,
            stars: None,
        };

        let len = stars.len();
        for star in stars {
            root.add_star(star);
        }

        Self { root, len }
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Append every star whose direction lies inside the spherical triangle.
    pub fn gather_triangle(&self, triangle: &[Vector3<f64>; 3], out: &mut Vec<Star>) {
        let planes = triangle_planes(triangle);
        self.root.gather(&planes, out);
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

    /// Classify this node's box against the triangle's side planes. Uses the same
    /// projection-radius trick as the SAT test: s is the box centre's signed
    /// distance to a plane and r is the box's half-width along the plane
    /// normal, so `[s - r, s + r]` is the box's shadow on that normal.
    fn classify(&self, planes: &[Vector3<f64>; 3]) -> Class {
        let c = self.bounds.center;
        let e = self.bounds.extents;
        let mut all_inside = true;
        for n in planes {
            let s = n.dot(&c);
            let r = e.x * n.x.abs() + e.y * n.y.abs() + e.z * n.z.abs();
            if s + r < 0.0 {
                return Class::Outside; // whole box on the wrong side of this plane
            }
            if s - r < 0.0 {
                all_inside = false; // box crosses this plane
            }
        }
        if all_inside {
            Class::Inside
        } else {
            Class::Straddle
        }
    }

    fn collect_all(&self, out: &mut Vec<Star>) {
        if let Some(stars) = &self.stars {
            out.extend(stars.iter().cloned());
        }
        if let Some(children) = &self.children {
            for child in children.iter() {
                child.collect_all(out);
            }
        }
    }

    fn gather(&self, planes: &[Vector3<f64>; 3], out: &mut Vec<Star>) {
        match self.classify(planes) {
            Class::Outside => {}
            Class::Inside => self.collect_all(out),
            Class::Straddle => {
                if let Some(children) = &self.children {
                    for child in children.iter() {
                        child.gather(planes, out);
                    }
                } else if let Some(stars) = &self.stars {
                    // Boundary leaf: test each star directly.
                    out.extend(
                        stars
                            .iter()
                            .filter(|s| point_in_triangle(&s.direction, planes))
                            .cloned(),
                    );
                }
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
        Node {
            bounds,
            children: None,
            stars: None,
        }
    }

    /// A wide spherical triangle around +z.
    fn plus_z_triangle() -> [Vector3<f64>; 3] {
        [
            Vector3::new(1.0, 1.0, 1.0),
            Vector3::new(-1.0, 1.0, 1.0),
            Vector3::new(0.0, -1.0, 1.0),
        ]
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

    // ---- Triangle query ---------------------------------------------------

    #[test]
    fn point_in_triangle_axis_and_behind() {
        let planes = triangle_planes(&plus_z_triangle());
        assert!(point_in_triangle(&Vector3::new(0.0, 0.0, 1.0), &planes)); // on axis
        assert!(point_in_triangle(&Vector3::new(0.1, 0.1, 1.0), &planes)); // inside
        assert!(!point_in_triangle(&Vector3::new(0.0, 0.0, -1.0), &planes)); // behind apex
        assert!(!point_in_triangle(&Vector3::new(2.0, 0.0, 0.1), &planes)); // off to +x
    }

    #[test]
    fn classify_box_inside_outside_straddle() {
        let planes = triangle_planes(&plus_z_triangle());

        // Small box well inside the triangle, above the apex on the +z axis.
        let inside = AABB::new(Vector3::new(-0.1, -0.1, 0.5), Vector3::new(0.1, 0.1, 0.7));
        assert!(matches!(leaf(inside).classify(&planes), Class::Inside));

        // Box entirely behind the apex.
        let outside = AABB::new(Vector3::new(-0.1, -0.1, -0.7), Vector3::new(0.1, 0.1, -0.5));
        assert!(matches!(leaf(outside).classify(&planes), Class::Outside));

        // The root box spans everything, so it crosses the boundary.
        assert!(matches!(
            leaf(unit_box()).classify(&planes),
            Class::Straddle
        ));
    }

    #[test]
    fn triangle_query_returns_only_inside_stars() {
        let inside = [
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(0.2, 0.1, 1.0),
            Vector3::new(-0.2, 0.1, 1.0),
        ];
        let outside = [
            Vector3::new(0.0, 0.0, -1.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, -1.0, 0.0),
            Vector3::new(-1.0, 0.0, -0.5),
        ];
        let stars: Vec<Star> = inside
            .iter()
            .chain(outside.iter())
            .map(|d| star_at(d.normalize()))
            .collect();

        let tree = Octree::new(stars);
        let mut found = Vec::new();
        tree.gather_triangle(&plus_z_triangle(), &mut found);

        assert_eq!(found.len(), inside.len());
        let planes = triangle_planes(&plus_z_triangle());
        assert!(found.iter().all(|s| point_in_triangle(&s.direction, &planes)));
    }

    #[test]
    fn triangle_query_across_subdivided_tree_matches_bruteforce() {
        // Enough spread stars to force subdivision, exercising the recursive
        // Inside/Straddle branches; the gather must equal the exact filter.
        let dirs: Vec<Vector3<f64>> = (0..200)
            .map(|i| {
                let t = i as f64;
                Vector3::new((t * 1.1).sin(), (t * 1.7).cos(), (t * 0.9).sin() + 0.3).normalize()
            })
            .collect();
        let stars: Vec<Star> = dirs.iter().map(|d| star_at(*d)).collect();
        let tree = Octree::new(stars);
        assert!(tree.root.children.is_some(), "test should subdivide");

        let planes = triangle_planes(&plus_z_triangle());
        let expected = dirs.iter().filter(|d| point_in_triangle(d, &planes)).count();
        let mut found = Vec::new();
        tree.gather_triangle(&plus_z_triangle(), &mut found);

        assert!(expected > 0 && expected < dirs.len(), "triangle should split the set");
        assert_eq!(found.len(), expected);
        assert!(found.iter().all(|s| point_in_triangle(&s.direction, &planes)));
    }

    // ---- Adaptive insertion ----------------------------------------------

    #[test]
    fn at_capacity_stays_a_single_leaf() {
        let stars: Vec<Star> = (0..Node::LEAF_CAPACITY)
            .map(|i| star_at(Vector3::new(0.9 - 0.05 * i as f64, 0.1, 0.1)))
            .collect();
        let tree = Octree::new(stars);
        assert!(
            tree.root.children.is_none(),
            "should not subdivide at capacity"
        );
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
        assert!(
            tree.root.children.is_some(),
            "should subdivide past capacity"
        );
        assert_eq!(
            count_stars(&tree.root),
            n,
            "no stars lost during subdivision"
        );
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

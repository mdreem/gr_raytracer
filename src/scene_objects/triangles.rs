use crate::geometry::coordinates::atlas::{Chart, GeometryChart};
use crate::geometry::coordinates::euclidean::CartesianChart;
use crate::geometry::geometry::Geometry;
use crate::geometry::point::Point;
use crate::rendering::integrator::Step;
use crate::rendering::radiance::Radiance;
use crate::rendering::raytracer::RaytracerError;
use crate::rendering::texture::UVCoordinates;
use crate::scene_objects::hittable::{ColorComputationData, Hittable, Intersection};

pub struct Triangle<C: Chart> {
    v1: Point<C>,
    v2: Point<C>,
    v3: Point<C>,
}

impl<C: Chart> Triangle<C> {
    pub fn new(v1: Point<C>, v2: Point<C>, v3: Point<C>) -> Triangle<C> {
        // TODO: Assert all points are in the same chart, or put them into the same chart.
        Self { v1, v2, v3 }
    }
}

// https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
impl<C: GeometryChart> Hittable<C> for Triangle<C> {
    fn intersects(
        &self,
        y_start: &Step<C>,
        y_end: &Step<C>,
        geometry: &dyn Geometry<Chart = C>,
    ) -> Option<Intersection<C>> {
        let ray_origin = y_start.x.vector;
        let ray_vector = y_end.x.vector - y_start.x.vector;

        let edge1 = self.v2.vector - self.v1.vector;
        let edge2 = self.v3.vector - self.v1.vector;

        let normal = edge1.cross(&edge2); // No need to normalize

        if normal.dot(&ray_vector) > 0.0 {
            return None;
        }

        let ray_cross_e2 = ray_vector.cross(&edge2);
        let det = edge1.dot(&ray_cross_e2);

        if det.abs() < f64::EPSILON {
            // Ray is parallel to triangle
            return None;
        }

        let inv_det = 1.0 / det;
        let s = ray_origin - self.v1.vector;
        let u = inv_det * s.dot(&ray_cross_e2);

        if u < -f64::EPSILON || u - 1.0 > f64::EPSILON {
            // Ray passes outside edge2's bounds
            return None;
        }

        let s_cross_e1 = -ray_vector.cross(&edge1);
        let v = inv_det * ray_vector.dot(&ray_cross_e2);

        if (v < -f64::EPSILON) || (u + v - 1.1 > f64::EPSILON) {
            // Ray passes outside edge1's bounds
            return None;
        }

        // The ray line intersects with the triangle.
        // We compute t to find where on the ray the intersection is.
        let t = inv_det * edge2.dot(&s_cross_e1);

        if t > f64::EPSILON && t < 1.0 + f64::EPSILON {
            let res = ray_origin + t * ray_vector;
            return Some(Intersection {
                uv: UVCoordinates { u: 1.0, v: 1.0 },
                intersection_point: Point::new(0.0, 0.0, 0.0, 0.0, y_start.x.chart),
                t: 1.0,
            });
        }

        None
    }

    fn color_at_uv(
        &self,
        color_computation_data: &ColorComputationData<C>,
        geometry: &dyn Geometry<Chart = C>,
    ) -> Result<Radiance, RaytracerError> {
        todo!()
    }

    fn energy_of_emitter(
        &self,
        geometry: &dyn Geometry<Chart = C>,
        step: &Step<C>,
    ) -> Result<f64, RaytracerError> {
        todo!()
    }

    fn temperature_of_emitter(
        &self,
        point: &Point<C>,
        geometry: &dyn Geometry<Chart = C>,
    ) -> Result<f64, RaytracerError> {
        todo!()
    }
}

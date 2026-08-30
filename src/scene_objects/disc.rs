use crate::geometry::geometry::Geometry;
use crate::geometry::point::Point;
use crate::rendering::color::CIETristimulus;
use crate::rendering::integrator::Step;
use crate::rendering::raytracer::RaytracerError;
use crate::rendering::temperature::TemperatureComputer;
use crate::rendering::texture::{TextureMapHandle, UVCoordinates};
use crate::scene_objects::hittable::{ColorComputationData, Hittable, Intersection};
use crate::scene_objects::objects::SceneObject;
use nalgebra::Vector3;
use roots::{SimpleConvergency, find_root_brent};

pub struct Disc {
    center_disk_inner_radius: f64,
    center_disk_outer_radius: f64,
    texture_mapper: TextureMapHandle,
    temperature_computer: Box<dyn TemperatureComputer>,
}

impl Disc {
    pub fn new(
        center_disk_inner_radius: f64,
        center_disk_outer_radius: f64,
        texture_mapper: TextureMapHandle,
        temperature_computer: Box<dyn TemperatureComputer>,
    ) -> Self {
        Self {
            center_disk_inner_radius,
            center_disk_outer_radius,
            texture_mapper,
            temperature_computer,
        }
    }

    fn create_intersection(
        &self,
        geometry: &dyn Geometry,
        center: Vector3<f64>,
        y_start_spatial: Vector3<f64>,
        direction: Vector3<f64>,
        t: f64,
    ) -> Option<Intersection> {
        let intersection_point = y_start_spatial + t * direction;
        let intersection_point_p = Point::new_cartesian(
            0.0,
            intersection_point[0],
            intersection_point[1],
            intersection_point[2],
        );
        // The disc's inner/outer radius are metric radial coordinates (the
        // convention the ISCO, orbit velocities, and temperature profile use).
        // In Kerr charts the equatorial Cartesian radius obeys
        // x^2 + y^2 = r^2 + a^2, so the Cartesian norm must not be compared
        // against them directly.
        let r = geometry.get_radial_coordinate(&intersection_point_p);

        if r >= self.center_disk_inner_radius && r <= self.center_disk_outer_radius {
            let vector_in_plane = intersection_point - center;

            let phi = vector_in_plane[1].atan2(vector_in_plane[0]); // phi in x-y plane.
            let r_normalized = (r - self.center_disk_inner_radius)
                / (self.center_disk_outer_radius - self.center_disk_inner_radius);

            let u = 0.5 + 0.5 * r_normalized * phi.cos();
            let v = 0.5 + 0.5 * r_normalized * phi.sin();

            return Some(Intersection {
                uv: UVCoordinates { u, v },
                intersection_point: intersection_point_p,
                t,
            });
        } else {
            return None;
        }
    }
}

impl Hittable for Disc {
    // TODO: explicitly construct the ray. Follow the integration. Some intervals seem to be skipped
    // here. See with current test setup. Intersection should be at t=7.63. With z=-2.442748091.
    // The intersection should be with an interval crossing y=0. But it seems to happen near 0 with
    // both coordinates.
    fn intersects(
        &self,
        y_start: &Step,
        y_end: &Step,
        geometry: &dyn Geometry,
    ) -> Option<Intersection> {
        // z x y
        let normal = Vector3::new(0.0, 0.0, 1.0);
        let center = Vector3::new(0.0, 0.0, 0.0);
        let y_start_spatial = y_start.x.get_spatial_vector_cartesian();
        let y_end_spatial = y_end.x.get_spatial_vector_cartesian();
        let direction = y_end_spatial - y_start_spatial;

        // Note: most of the following assumes that normal will be in the z direction.

        // Start and end point along the normal direction.
        let point_start = (center - y_start_spatial).dot(&normal);
        let point_end = direction.dot(&normal);

        // The segment crosses the plane.
        let z0 = y_start_spatial[2];
        let z1 = y_end_spatial[2];
        if z0.signum() != z1.signum() {
            // TODO: p2 can be 0 if parallel -> handle
            let t = point_start / point_end; // plane intersection parameter.

            if !(0.0..=1.0).contains(&t) {
                return None;
            }

            return self.create_intersection(geometry, center, y_start_spatial, direction, t);
        }

        // If the segment does not cross the plane, check if the cubic Hermite interpolation of
        // the z-coordinate intersects the plane.
        let d_lambda = y_end.t - y_start.t;
        let p0 = y_start_spatial[2];
        let p1 = y_end_spatial[2];
        let m0 = y_start.p.get_z_cartesian(&y_start.x) * d_lambda;
        let m1 = y_end.p.get_z_cartesian(&y_end.x) * d_lambda;

        // p(t) = (2 t^3 - 3 t^2 + 1) p0 + (t^3 - 2 t^2 + t) m0 + (-2 t^3 + 3 t^2) p1 + (t^3 - t^2) m1
        //      = p0 + m0 t + (-2 m0 - m1 - 3 p0 + 3 p1) t^2 + (m0 + m1 + 2 p0 - 2 p1) t^3
        // p'(t) = (6 t^2 - 6 t) p0 + (3 t^2 - 4 t + 1) m0 + (-6 t^2 + 6 t) p1 + (3 t^2 - 2 t) m1
        //       = m0 + (-4 m0 - 2 m1 - 6 p0 + 6 p1) t + (3 m0 + 3 m1 + 6 p0 - 6 p1) t^2

        let c0 = p0;
        let c1 = m0;
        let c2 = -3.0 * p0 - 2.0 * m0 + 3.0 * p1 - m1;
        let c3 = 2.0 * p0 + m0 - 2.0 * p1 + m1;

        let a = 3.0 * c3;
        let b = 2.0 * c2;
        let c = c1;
        let discriminant = b * b - 4.0 * a * c;

        // If the discriminant is negative, there are no real roots, and thus no intersection with the plane.
        // We already checked that the segment does not cross the plane.
        if discriminant < 0.0 {
            return None;
        }
        // Safe way to compute the roots of a quadratic equation. This will not blow up all values
        // if a is very small. See https://people.csail.mit.edu/bkph/articles/Quadratics.pdf
        let q = -0.5 * (b + discriminant.sqrt().copysign(b));
        let int_p1 = q / a;
        let int_p2 = c / q;

        let mut points = vec![0.0, 1.0];
        for tp in [int_p1, int_p2] {
            if tp.is_finite() && tp > 0.0 && tp < 1.0 {
                points.push(tp);
            }
        }
        points.sort_by(f64::total_cmp);

        for points in points.windows(2) {
            let p1 = points[0];
            let p2 = points[1];

            let f = |s: f64| c0 + c1 * s + c2 * s * s + c3 * s * s * s;
            let mut conv = SimpleConvergency {
                eps: 1e-12,
                max_iter: 100,
            };
            let root = find_root_brent(p1, p2, &f, &mut conv);
            let t = match root {
                Ok(t) => t,
                Err(_) => continue,
            };

            let intersection =
                self.create_intersection(geometry, center, y_start_spatial, direction, t);

            if intersection.is_some() {
                return intersection;
            }
        }
        None
    }

    fn color_at_uv(
        &self,
        color_computation_data: &ColorComputationData,
        _geometry: &dyn Geometry,
    ) -> Result<CIETristimulus, RaytracerError> {
        self.texture_mapper.color_at_uv(
            &color_computation_data.uv,
            &color_computation_data.temperature_data,
        )
    }

    fn energy_of_emitter(
        &self,
        geometry: &dyn Geometry,
        step: &Step,
    ) -> Result<f64, RaytracerError> {
        let position = step.x;
        let velocity = geometry.get_circular_orbit_velocity_at(&position)?;
        let momentum = step.p;
        Ok(geometry.inner_product(&position, &velocity, &momentum))
    }

    fn temperature_of_emitter(
        &self,
        point: &Point,
        geometry: &dyn Geometry,
    ) -> Result<f64, RaytracerError> {
        let r = geometry.get_radial_coordinate(point);
        let temperature = self.temperature_computer.compute_temperature(r)?;
        Ok(temperature)
    }
}

impl SceneObject for Disc {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::euclidean::EuclideanSpace;
    use crate::geometry::four_vector::FourVector;
    use crate::geometry::kerr::Kerr;
    use crate::rendering::color::Color;
    use crate::rendering::temperature::ConstantTemperatureComputer;
    use crate::rendering::texture::CheckerMapper;
    use std::sync::Arc;

    fn create_disc(inner: f64, outer: f64) -> Disc {
        Disc::new(
            inner,
            outer,
            Arc::new(CheckerMapper::new(
                0.0,
                5.0,
                5.0,
                Color::new(100, 0, 0, 255),
                Color::new(0, 100, 0, 255),
            )),
            Box::new(ConstantTemperatureComputer::new(1000.0)),
        )
    }

    fn crossing_at(x: f64) -> (Step, Step) {
        let start = Point::new_cartesian(0.0, x, 0.0, 0.1);
        let end = Point::new_cartesian(0.0, x, 0.0, -0.1);
        // Straight-line tangent, so the Cubic Hermite reduces to the linear crossing this fixture describes.
        let dir = end.get_spatial_vector_cartesian() - start.get_spatial_vector_cartesian();
        let p = FourVector::new_cartesian(1.0, dir[0], dir[1], dir[2]);
        (
            Step {
                x: start,
                p,
                t: 0.0,
                step: 0,
            },
            Step {
                x: end,
                p,
                t: 1.0,
                step: 1,
            },
        )
    }

    #[test]
    fn test_disc_boundary_uses_metric_radial_coordinate_in_kerr() {
        // Kerr equatorial relation: x^2 + y^2 = r^2 + a^2. A plane crossing
        // at Cartesian radius R corresponds to Boyer-Lindquist
        // r = sqrt(R^2 - a^2), so points slightly outside the inner edge in
        // Cartesian norm can lie INSIDE it in the metric radial coordinate.
        let geometry = Kerr::new(1.0, 0.499, 1e-4);
        let disc = create_disc(0.795, 8.0);

        // Cartesian R = 0.8 -> BL r = sqrt(0.64 - 0.499^2) ~ 0.625 < 0.795:
        // no gas there.
        let (start, end) = crossing_at(0.8);
        assert!(disc.intersects(&start, &end, &geometry).is_none());

        // Cartesian R = 0.95 -> BL r ~ 0.808 > 0.795: inside the disc.
        let (start, end) = crossing_at(0.95);
        assert!(disc.intersects(&start, &end, &geometry).is_some());

        // Outer edge: Cartesian R = 8.01 -> BL r ~ 7.994 < 8.0: still gas.
        let (start, end) = crossing_at(8.01);
        assert!(disc.intersects(&start, &end, &geometry).is_some());
    }

    #[test]
    fn test_disc_boundary_unchanged_in_flat_space() {
        let geometry = EuclideanSpace::new();
        let disc = create_disc(3.05, 8.0);

        let (start, end) = crossing_at(3.0);
        assert!(disc.intersects(&start, &end, &geometry).is_none());
        let (start, end) = crossing_at(3.1);
        assert!(disc.intersects(&start, &end, &geometry).is_some());
        let (start, end) = crossing_at(8.1);
        assert!(disc.intersects(&start, &end, &geometry).is_none());
    }

    #[test]
    fn test_disc_uv_normalization_uses_metric_radius() {
        // At the disc's inner edge the normalized UV radius must be ~0
        // (texture inner boundary), evaluated in the metric radial
        // coordinate, not the Cartesian norm.
        let geometry = Kerr::new(1.0, 0.499, 1e-4);
        let disc = create_disc(0.795, 8.0);

        // Cartesian R such that BL r == inner edge: R = sqrt(r^2 + a^2).
        let r_inner_cartesian = (0.795f64 * 0.795 + 0.499 * 0.499).sqrt();
        let (start, end) = crossing_at(r_inner_cartesian + 1e-6);
        let intersection = disc
            .intersects(&start, &end, &geometry)
            .expect("should hit the inner edge");
        // u = 0.5 + 0.5 * r_normalized * cos(phi); at phi = 0 and
        // r_normalized ~ 0 this is ~0.5.
        assert!((intersection.uv.u - 0.5).abs() < 1e-3);
        assert!((intersection.uv.v - 0.5).abs() < 1e-3);
    }
}

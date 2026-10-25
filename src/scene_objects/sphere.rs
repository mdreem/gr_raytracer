use crate::geometry::coordinates::atlas::GeometryChart;
use crate::geometry::coordinates::euclidean::CartesianChart;
use crate::geometry::geometry::Geometry;
use crate::geometry::point::Point;
use crate::geometry::spherical_coordinates_helper::cartesian_to_spherical;
use crate::rendering::integration::integrator::Step;
use crate::rendering::light::radiance::Radiance;
use crate::rendering::error::RaytracerError;
use crate::rendering::texture::{TextureMapHandle, UVCoordinates};
use crate::scene_objects::hittable::{ColorComputationData, Hittable, Intersection};
use nalgebra::Vector3;
use std::f64::consts::PI;

pub struct Sphere {
    radius: f64,
    texture_mapper: TextureMapHandle,
    position: Point<CartesianChart>,
    temperature: f64,
}

impl Sphere {
    pub fn new(
        radius: f64,
        texture_mapper: TextureMapHandle,
        position: Point<CartesianChart>,
        temperature: f64,
    ) -> Self {
        Self {
            radius,
            texture_mapper,
            position,
            temperature,
        }
    }
}

fn solve_for_t(y_start_spatial: Vector3<f64>, direction: Vector3<f64>, r: f64) -> Option<f64> {
    let a = direction.dot(&direction);
    let b = 2.0 * y_start_spatial.dot(&direction);
    let c = y_start_spatial.dot(&y_start_spatial) - r * r;

    let discriminant = b * b - 4.0 * a * c;

    if discriminant < 0.0 {
        None
    } else {
        let sqrt_disc = discriminant.sqrt();
        // t_near <= t_far since a = |direction|^2 > 0. Return the nearest
        // surface crossing inside the segment: t_near for an entering or
        // tunnelling ray, t_far when the segment starts inside the sphere.
        let t_near = (-b - sqrt_disc) / (2.0 * a);
        let t_far = (-b + sqrt_disc) / (2.0 * a);
        if (0.0..=1.0).contains(&t_near) {
            Some(t_near)
        } else if (0.0..=1.0).contains(&t_far) {
            Some(t_far)
        } else {
            None
        }
    }
}

impl<C: GeometryChart> Hittable<C> for Sphere {
    // y_start and y_end have to be Cartesian.
    fn intersects(
        &self,
        y_start: &Step<C>,
        y_end: &Step<C>,
        _geometry: &dyn Geometry<Chart = C>,
    ) -> Option<Intersection<CartesianChart>> {
        let y_start_cartesian = y_start.x.to_cartesian();
        let y_end_cartesian = y_end.x.to_cartesian();

        let neg_position = -self.position;
        let y_start_shifted = y_start_cartesian + neg_position;
        let y_end_shifted = y_end_cartesian + neg_position;

        let y_start_spatial = y_start_shifted.get_spatial_vector_cartesian();
        let y_end_spatial = y_end_shifted.get_spatial_vector_cartesian();
        let direction = y_end_spatial - y_start_spatial;

        // Solve the ray-sphere quadratic over the whole segment rather than
        // gating on the endpoints' radii: a step that passes fully through the
        // sphere leaves both endpoints outside, so an endpoint-inside/outside
        // test would tunnel straight through it. The quadratic classifies every
        // case (no hit -> negative discriminant or roots outside [0,1]).
        let t = solve_for_t(y_start_spatial, direction, self.radius)?;

        // UV mapping needs the hit point in the sphere's own local frame
        // (so e.g. the texture's north pole is always the sphere's own
        // north pole, regardless of where the sphere sits in the scene).
        let point_on_sphere_spatial = y_start_spatial + t * direction;
        let point_on_sphere_local = cartesian_to_spherical(&Point::new_cartesian(
            0.0,
            point_on_sphere_spatial[0],
            point_on_sphere_spatial[1],
            point_on_sphere_spatial[2],
        ));

        let theta = point_on_sphere_local[2];
        let phi = point_on_sphere_local[3];
        let u = (PI + phi) / (2.0 * PI);
        let v = theta / PI;

        // `Intersection::intersection_point` itself must be world-space:
        // it's used (via `energy_of_emitter`) to evaluate the geometry's
        // fields (e.g. the stationary observer's velocity) at the
        // emitter's actual location relative to the black hole, not
        // relative to the sphere's own center.
        let position_spatial = self.position.get_spatial_vector_cartesian();
        let point_on_sphere_world = Point::new_cartesian(
            0.0,
            point_on_sphere_spatial[0] + position_spatial[0],
            point_on_sphere_spatial[1] + position_spatial[1],
            point_on_sphere_spatial[2] + position_spatial[2],
        );

        Some(Intersection {
            uv: UVCoordinates { u: 1.0 - u, v },
            intersection_point: point_on_sphere_world,
            t,
        })
    }

    fn color_at_uv(
        &self,
        color_computation_data: &ColorComputationData<C>,
        _geometry: &dyn Geometry<Chart = C>,
    ) -> Result<Radiance, RaytracerError> {
        self.texture_mapper
            .color_at_uv(
                &color_computation_data.uv,
                &color_computation_data.temperature_data,
            )
            .map(Radiance::from_straight)
    }

    fn energy_of_emitter(
        &self,
        geometry: &dyn Geometry<Chart = C>,
        step: &Step<C>,
    ) -> Result<f64, RaytracerError> {
        let position = step.x;
        let velocity = geometry.get_stationary_velocity_at(&position);
        let momentum = step.p;
        Ok(geometry.inner_product(&position, &velocity, &momentum))
    }

    fn temperature_of_emitter(
        &self,
        _point: &Point<C>,
        _geometry: &dyn Geometry<Chart = C>,
    ) -> Result<f64, RaytracerError> {
        Ok(self.temperature)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::spacetimes::euclidean::EuclideanSpace;
    use crate::geometry::four_vector::FourVector;
    use crate::geometry::point::Point;
    use crate::rendering::light::color::Color;
    use crate::rendering::texture::CheckerMapper;
    use crate::scene_objects::hittable::Hittable;
    use approx::assert_abs_diff_eq;
    use std::sync::Arc;

    // Sphere intersection only reads the endpoint positions; the momentum is
    // a placeholder so the tests can build `Step`s.
    fn as_step(p: Point<CartesianChart>) -> Step<CartesianChart> {
        Step {
            x: p,
            p: FourVector::new(1.0, 0.0, 0.0, 0.0, p),
            t: 0.0,
            step: 0,
        }
    }

    fn create_sphere_at(x: f64, y: f64, z: f64) -> Sphere {
        Sphere::new(
            1.0,
            Arc::new(CheckerMapper::new(
                3.0,
                5.0,
                5.0,
                Color::new(100, 0, 0, 255),
                Color::new(0, 100, 0, 255),
            )),
            Point::new_cartesian(0.0, x, y, z),
            0.0,
        )
    }

    #[test]
    fn test_sphere_intersection_center_sphere() {
        let sphere = create_sphere_at(0.0, 0.0, 0.0);
        let y_start = Point::new_cartesian(0.0, 1.1, 0.0, 0.0);
        let y_end = Point::new_cartesian(0.0, 0.9, 0.0, 0.0);

        assert!(
            sphere
                .intersects(&as_step(y_start), &as_step(y_end), &EuclideanSpace::new())
                .is_some()
        );
    }

    #[test]
    fn test_sphere_intersection_center_sphere_no_intersection() {
        let sphere = create_sphere_at(0.0, 0.0, 0.0);
        let y_start = Point::new_cartesian(0.0, 1.1, 0.0, 0.0);
        let y_end = Point::new_cartesian(0.0, 1.01, 0.0, 0.0);

        assert!(
            sphere
                .intersects(&as_step(y_start), &as_step(y_end), &EuclideanSpace::new())
                .is_none()
        );
    }

    #[test]
    fn test_sphere_intersection_moved_sphere() {
        let sphere = create_sphere_at(5.0, 0.0, 0.0);
        let y_start = Point::new_cartesian(0.0, 6.1, 0.0, 0.0);
        let y_end = Point::new_cartesian(0.0, 5.9, 0.0, 0.0);

        assert!(
            sphere
                .intersects(&as_step(y_start), &as_step(y_end), &EuclideanSpace::new())
                .is_some()
        );
    }

    #[test]
    fn test_sphere_intersection_moved_sphere_misses() {
        let sphere = create_sphere_at(5.0, 0.0, 0.0);
        let y_start = Point::new_cartesian(0.0, 6.1, 0.0, 0.0);
        let y_end = Point::new_cartesian(0.0, 6.01, 0.0, 0.0);

        assert!(
            sphere
                .intersects(&as_step(y_start), &as_step(y_end), &EuclideanSpace::new())
                .is_none()
        );
    }

    #[test]
    fn test_sphere_intersection_point_is_world_space_not_sphere_local() {
        // `intersection_point` is used downstream (via `energy_of_emitter`)
        // to evaluate the black hole's fields at the emitter's true
        // location, so it must be in world coordinates. Regression test for
        // a bug where it was returned in the sphere's own local frame
        // (i.e. relative to the sphere's center) instead.
        let sphere = create_sphere_at(0.0, 0.0, 20.0);
        let y_start = Point::new_cartesian(0.0, 0.0, 0.0, 22.0);
        let y_end = Point::new_cartesian(0.0, 0.0, 0.0, 19.5);

        let intersection = sphere
            .intersects(&as_step(y_start), &as_step(y_end), &EuclideanSpace::new())
            .expect("ray should hit sphere");
        let hit = intersection
            .intersection_point
            .get_spatial_vector_cartesian();

        // World-space hit point is near z = 21 (top of the sphere, facing
        // the ray's origin). A sphere-local point would incorrectly read
        // z ≈ 1 (relative to the sphere's own center at z = 20).
        assert_abs_diff_eq!(hit[2], 21.0, epsilon = 1e-9);
    }
}

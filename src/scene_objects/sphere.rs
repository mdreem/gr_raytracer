use crate::geometry::point::{CoordinateSystem, Point};
use crate::geometry::spherical_coordinates_helper::cartesian_to_spherical;
use crate::rendering::color::CIETristimulus;
use crate::rendering::texture::{TextureMap, TextureMapHandle, UVCoordinates};
use crate::scene_objects::hittable::{Hittable, Intersection};
use crate::scene_objects::objects::SceneObject;
use nalgebra::Vector3;
use std::f64::consts::PI;

pub struct Sphere {
    radius: f64,
    texture_mapper: TextureMapHandle,
    position: Point,
}

impl Sphere {
    pub fn new(radius: f64, texture_mapper: TextureMapHandle, position: Point) -> Self {
        Self {
            radius,
            texture_mapper,
            position,
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
        let t1 = (-b + sqrt_disc) / (2.0 * a);
        let t2 = (-b - sqrt_disc) / (2.0 * a);
        if (0.0..=1.0).contains(&t1) {
            Some(t1)
        } else if (0.0..=1.0).contains(&t2) {
            Some(t2)
        } else {
            None
        }
    }
}

impl Hittable for Sphere {
    // y_start and y_end have to be Cartesian.
    fn intersects(&self, y_start: &Point, y_end: &Point) -> Option<Intersection> {
        let neg_position = -self.position;
        let y_start_shifted = y_start + &neg_position;
        let y_end_shifted = y_end + &neg_position;
        let r_start = y_start_shifted.radial_distance_spatial_part_squared();
        let r_end = y_end_shifted.radial_distance_spatial_part_squared();

        // Checks if the line element intersects the surface of the sphere.
        if (r_start >= self.radius.powi(2) && r_end <= self.radius.powi(2))
            || (r_start <= self.radius.powi(2) && r_end >= self.radius.powi(2))
        {
            let y_start_spatial = y_start_shifted.get_spatial_vector();
            let y_end_spatial = y_end_shifted.get_spatial_vector();
            let direction = y_end_spatial - y_start_spatial;

            let t = match solve_for_t(y_start_spatial, direction, self.radius) {
                None => {
                    return None;
                }
                Some(t) => t,
            };

            let point_on_sphere_spatial = y_start_spatial + t * direction;
            let point_on_sphere = cartesian_to_spherical(&Point::new(
                0.0,
                point_on_sphere_spatial[0],
                point_on_sphere_spatial[1],
                point_on_sphere_spatial[2],
                CoordinateSystem::Cartesian,
            ));

            let theta = point_on_sphere[2];
            let phi = point_on_sphere[3];
            let u = (PI + phi) / (2.0 * PI);
            let v = theta / PI;

            return Some(Intersection {
                uv: UVCoordinates { u, v: 1.0 - v },
                intersection_point: point_on_sphere,
            });
        }

        None
    }

    fn color_at_uv(&self, uv: UVCoordinates, redshift: f64) -> CIETristimulus {
        self.texture_mapper.color_at_uv(uv, redshift)
    }
}

impl SceneObject for Sphere {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::point::Point;
    use crate::rendering::color::Color;
    use crate::rendering::texture::CheckerMapper;
    use crate::scene_objects::hittable::Hittable;
    use std::sync::Arc;

    fn create_sphere_at(x: f64, y: f64, z: f64) -> Sphere {
        Sphere::new(
            1.0,
            Arc::new(CheckerMapper::new(
                5.0,
                5.0,
                Color::new(100, 0, 0, 255),
                Color::new(0, 100, 0, 255),
            )),
            Point::new_cartesian(0.0, x, y, z),
        )
    }

    #[test]
    fn test_sphere_intersection_center_sphere() {
        let sphere = create_sphere_at(0.0, 0.0, 0.0);
        let y_start = Point::new_cartesian(0.0, 1.1, 0.0, 0.0);
        let y_end = Point::new_cartesian(0.0, 0.9, 0.0, 0.0);

        assert!(sphere.intersects(&y_start, &y_end).is_some());
    }

    #[test]
    fn test_sphere_intersection_center_sphere_no_intersection() {
        let sphere = create_sphere_at(0.0, 0.0, 0.0);
        let y_start = Point::new_cartesian(0.0, 1.1, 0.0, 0.0);
        let y_end = Point::new_cartesian(0.0, 1.01, 0.0, 0.0);

        assert!(sphere.intersects(&y_start, &y_end).is_none());
    }

    #[test]
    fn test_sphere_intersection_moved_sphere() {
        let sphere = create_sphere_at(5.0, 0.0, 0.0);
        let y_start = Point::new_cartesian(0.0, 6.1, 0.0, 0.0);
        let y_end = Point::new_cartesian(0.0, 5.9, 0.0, 0.0);

        assert!(sphere.intersects(&y_start, &y_end).is_some());
    }

    #[test]
    fn test_sphere_intersection_moved_sphere_misses() {
        let sphere = create_sphere_at(5.0, 0.0, 0.0);
        let y_start = Point::new_cartesian(0.0, 6.1, 0.0, 0.0);
        let y_end = Point::new_cartesian(0.0, 6.01, 0.0, 0.0);

        assert!(sphere.intersects(&y_start, &y_end).is_none());
    }
}

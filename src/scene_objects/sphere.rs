use crate::geometry::point::{CoordinateSystem, Point};
use crate::geometry::spherical_coordinates_helper::cartesian_to_spherical;
use crate::rendering::color::Color;
use crate::rendering::texture::{TextureMap, UVCoordinates};
use crate::scene_objects::objects::SceneObject;
use nalgebra::Vector3;
use std::f64::consts::PI;

pub struct Sphere<T: TextureMap> {
    radius: f64,
    texture_mapper: T,
}

impl<T: TextureMap> Sphere<T> {
    pub fn new(radius: f64, texture_mapper: T) -> Self {
        Self {
            radius,
            texture_mapper,
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

impl<T: TextureMap> crate::scene_objects::hittable::Hittable for Sphere<T> {
    fn intersects(&self, y_start: &Point, y_end: &Point) -> Option<UVCoordinates> {
        let r_start = y_start.radial_distance_spatial_part_squared();
        let r_end = y_end.radial_distance_spatial_part_squared();

        if (r_start >= self.radius.powi(2) && r_end <= self.radius.powi(2))
            || (r_start <= self.radius.powi(2) && r_end >= self.radius.powi(2))
        {
            let y_start_spatial = y_start.get_spatial_vector();
            let y_end_spatial = y_end.get_spatial_vector();
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

            return Some(UVCoordinates { u, v: 1.0 - v });
        }

        None
    }
}

impl<T: TextureMap> TextureMap for Sphere<T> {
    fn color_at_uv(&self, uv: UVCoordinates) -> Color {
        self.texture_mapper.color_at_uv(uv)
    }
}

impl<T: TextureMap> SceneObject for Sphere<T> {}

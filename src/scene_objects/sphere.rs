use crate::color::Color;
use crate::four_vector::FourVector;
use crate::scene_objects::objects::SceneObject;
use crate::texture::{TextureMap, UVCoordinates};
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

impl<T: TextureMap> crate::hittable::Hittable for Sphere<T> {
    fn intersects(&self, y_start: &FourVector, y_end: &FourVector) -> Option<UVCoordinates> {
        let r_start = y_start.radial_distance_spatial_part_squared();
        let r_end = y_end.radial_distance_spatial_part_squared();

        if (r_start >= self.radius.powi(2) && r_end <= self.radius.powi(2))
            || (r_start <= self.radius.powi(2) && r_end >= self.radius.powi(2))
        {
            let point_on_sphere = y_start.get_as_spherical(); // approximate y_start als intersection point.

            let theta = point_on_sphere[1];
            let phi = point_on_sphere[2];
            let u = (PI + phi) / (2.0 * PI);
            let v = theta / PI;

            return Some(UVCoordinates { u, v });
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

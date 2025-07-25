use crate::color::Color;
use crate::geometry::four_vector::FourVector;
use crate::scene_objects::hittable::Hittable;
use crate::scene_objects::objects::SceneObject;
use crate::texture::{TextureMap, UVCoordinates};
use nalgebra::Vector3;
use std::f64::consts::PI;

pub struct Disc<T: TextureMap> {
    center_disk_inner_radius: f64,
    center_disk_outer_radius: f64,
    texture_mapper: T,
}

impl<T: TextureMap> Disc<T> {
    pub fn new(
        center_disk_inner_radius: f64,
        center_disk_outer_radius: f64,
        texture_mapper: T,
    ) -> Self {
        Self {
            center_disk_inner_radius,
            center_disk_outer_radius,
            texture_mapper,
        }
    }
}

impl<T: TextureMap> Hittable for Disc<T> {
    // TODO: explicitly construct the ray. Follow the integration. Some intervals seem to be skipped
    // here. See with current test setup. Intersection should be at t=7.63. With z=-2.442748091.
    // The intersection should be with an interval crossing y=0. But it seems to happen near 0 with
    // both coordinates.
    fn intersects(&self, y_start: &FourVector, y_end: &FourVector) -> Option<UVCoordinates> {
        // z x y
        let normal = Vector3::new(0.0, 1.0, 0.0);
        let center = Vector3::new(0.0, 0.0, 0.0);
        let y_start_spatial = y_start.get_spatial_vector();
        let y_end_spatial = y_end.get_spatial_vector();
        let direction = y_end_spatial - y_start_spatial;

        let p1 = (y_start_spatial - center).dot(&normal);
        let p2 = direction.dot(&normal);

        // TODO: p2 can be 0 if parallel -> handle
        let t = p1 / p2; // plane intersection parameter.

        if !(0.0..=1.0).contains(&t) {
            return None;
        }

        let intersection_point = y_start_spatial + t * direction;
        let rr = intersection_point.norm_squared();

        if rr >= self.center_disk_inner_radius * self.center_disk_inner_radius
            && rr <= self.center_disk_outer_radius * self.center_disk_outer_radius
        {
            let vector_in_plane = intersection_point - center;
            // TODO: properly implement finding intersection and compute the values accordingly.

            let phi = vector_in_plane[2].atan2(vector_in_plane[0]); // phi in x-z plane.
            let u = (PI + phi) / (2.0 * PI);
            let v = (rr.sqrt() - self.center_disk_inner_radius)
                / (self.center_disk_outer_radius - self.center_disk_inner_radius);

            Some(UVCoordinates { u, v })
        } else {
            None
        }
    }
}

impl<T: TextureMap> TextureMap for Disc<T> {
    fn color_at_uv(&self, uv: UVCoordinates) -> Color {
        self.texture_mapper.color_at_uv(uv)
    }
}

impl<T: TextureMap> SceneObject for Disc<T> {}

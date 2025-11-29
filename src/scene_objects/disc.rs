use crate::geometry::geometry::Geometry;
use crate::geometry::point::Point;
use crate::rendering::color::CIETristimulus;
use crate::rendering::integrator::Step;
use crate::rendering::texture::{TemperatureData, TextureMapHandle, UVCoordinates};
use crate::scene_objects::hittable::{Hittable, Intersection};
use crate::scene_objects::objects::SceneObject;
use nalgebra::Vector3;

pub struct Disc {
    center_disk_inner_radius: f64,
    center_disk_outer_radius: f64,
    texture_mapper: TextureMapHandle,
    temperature: f64,
}

impl Disc {
    pub fn new(
        center_disk_inner_radius: f64,
        center_disk_outer_radius: f64,
        texture_mapper: TextureMapHandle,
        temperature: f64,
    ) -> Self {
        Self {
            center_disk_inner_radius,
            center_disk_outer_radius,
            texture_mapper,
            temperature,
        }
    }
}

impl Hittable for Disc {
    // TODO: explicitly construct the ray. Follow the integration. Some intervals seem to be skipped
    // here. See with current test setup. Intersection should be at t=7.63. With z=-2.442748091.
    // The intersection should be with an interval crossing y=0. But it seems to happen near 0 with
    // both coordinates.
    fn intersects(&self, y_start: &Point, y_end: &Point) -> Option<Intersection> {
        // z x y
        let normal = Vector3::new(0.0, 0.0, 1.0);
        let center = Vector3::new(0.0, 0.0, 0.0);
        let y_start_spatial = y_start.get_spatial_vector_cartesian();
        let y_end_spatial = y_end.get_spatial_vector_cartesian();
        let direction = y_end_spatial - y_start_spatial;

        let p1 = (center - y_start_spatial).dot(&normal);
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

            let phi = vector_in_plane[1].atan2(vector_in_plane[0]); // phi in x-y plane.
            let r = (rr.sqrt() - self.center_disk_inner_radius)
                / (self.center_disk_outer_radius - self.center_disk_inner_radius);

            let u = 0.5 + 0.5 * r * phi.cos();
            let v = 0.5 + 0.5 * r * phi.sin();

            Some(Intersection {
                uv: UVCoordinates { u, v },
                intersection_point: Point::new_cartesian(
                    0.0,
                    intersection_point[0],
                    intersection_point[1],
                    intersection_point[2],
                ),
            })
        } else {
            None
        }
    }

    fn color_at_uv(&self, uv: UVCoordinates, temperature_data: TemperatureData) -> CIETristimulus {
        self.texture_mapper.color_at_uv(uv, temperature_data)
    }

    fn energy_of_emitter(&self, geometry: &dyn Geometry, step: &Step) -> f64 {
        let position = step.x;
        let velocity = geometry.get_circular_orbit_velocity_at(&position);
        let momentum = step.p;
        geometry.inner_product(&position, &velocity, &momentum)
    }

    fn temperature_of_emitter(&self, _geometry: &dyn Geometry, step: &Step) -> f64 {
        // TODO: implement.
        let position = step.x;
        let r = position.get_as_spherical()[0];

        self.temperature * 0.7 * (1.0 / r).powf(0.75)
    }
}

impl SceneObject for Disc {}

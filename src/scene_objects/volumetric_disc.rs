use crate::geometry::geometry::Geometry;
use crate::geometry::point::Point;
use crate::rendering::color::CIETristimulus;
use crate::rendering::integrator::Step;
use crate::rendering::raytracer::RaytracerError;
use crate::rendering::temperature::TemperatureComputer;
use crate::rendering::texture::{TemperatureData, TextureMapHandle, UVCoordinates};
use crate::scene_objects::hittable::{ColorComputationData, Hittable, Intersection};
use crate::scene_objects::objects::SceneObject;

pub struct VolumetricDisc {
    center_disk_inner_radius: f64,
    center_disk_outer_radius: f64,
    texture_mapper: TextureMapHandle,
    temperature_computer: Box<dyn TemperatureComputer>,
}

impl VolumetricDisc {
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
}

impl Hittable for VolumetricDisc {
    // TODO: ensure sphere is not skipped due to large step sizes.
    fn intersects(&self, y_start: &Point, y_end: &Point) -> Option<Intersection> {
        let y_start_spatial = y_start.get_spatial_vector_cartesian();
        let y_end_spatial = y_end.get_spatial_vector_cartesian();
        let direction = y_end_spatial - y_start_spatial;

        let n_start = y_start_spatial.norm_squared().sqrt();
        let n_end = y_end_spatial.norm_squared().sqrt();

        // If close to center, no intersection.
        if n_end < self.center_disk_inner_radius || n_start < self.center_disk_inner_radius {
            return None;
        }
        // If both points are outside outer radius, no intersection.
        if n_start > self.center_disk_outer_radius && n_end > self.center_disk_outer_radius {
            return None;
        }

        // Check if the segment crosses the outer radius.
        let t = if (n_start < self.center_disk_outer_radius
            && n_end > self.center_disk_outer_radius)
            || n_start > self.center_disk_outer_radius && n_end < self.center_disk_outer_radius
        {
            (self.center_disk_outer_radius - n_start) / (n_end - n_start)
        } else {
            0.0
        };

        let intersection_point = y_start_spatial + t * direction;

        Some(Intersection {
            uv: UVCoordinates { u: 0.0, v: 0.0 },
            intersection_point: Point::new_cartesian(
                0.0,
                intersection_point[0],
                intersection_point[1],
                intersection_point[2],
            ),
        })
    }

    fn color_at_uv(&self, color_computation_data: &ColorComputationData) -> CIETristimulus {
        let disc_thickness = (color_computation_data.intersection_point[3].abs() / 0.5).exp();
        if color_computation_data.intersection_point[3].abs() > 0.5 {
            CIETristimulus::new(0.0, 0.0, 0.0, 0.0)
        } else {
            CIETristimulus::new(100.0, 0.0, disc_thickness, 1.0)
        }
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

    fn temperature_of_emitter(&self, point: &Point) -> Result<f64, RaytracerError> {
        let r = point.get_as_spherical()[0];
        let temperature = self.temperature_computer.compute_temperature(r)?;
        Ok(temperature)
    }
}

impl SceneObject for VolumetricDisc {}

use crate::geometry::geometry::Geometry;
use crate::geometry::point::Point;
use crate::rendering::color::CIETristimulus;
use crate::rendering::integrator::Step;
use crate::rendering::raytracer::RaytracerError;
use crate::rendering::texture::{TemperatureData, UVCoordinates};

pub struct Intersection {
    pub uv: UVCoordinates,
    pub intersection_point: Point,
}

pub struct ColorComputationData {
    pub uv: UVCoordinates,
    pub temperature_data: TemperatureData,
    pub intersection_point: Point,
}

pub trait Hittable: Sync {
    fn intersects(&self, y_start: &Point, y_end: &Point) -> Option<Intersection>;
    fn color_at_uv(&self, color_computation_data: &ColorComputationData) -> CIETristimulus;
    fn energy_of_emitter(
        &self,
        geometry: &dyn Geometry,
        step: &Step,
    ) -> Result<f64, RaytracerError>;
    fn temperature_of_emitter(&self, point: &Point) -> Result<f64, RaytracerError>;
}

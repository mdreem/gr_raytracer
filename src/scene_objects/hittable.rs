use crate::geometry::four_vector::FourVector;
use crate::geometry::geometry::Geometry;
use crate::geometry::point::Point;
use crate::rendering::color::CIETristimulus;
use crate::rendering::integrator::Step;
use crate::rendering::raytracer::RaytracerError;
use crate::rendering::redshift::RayFrequencyData;
use crate::rendering::texture::{TemperatureData, UVCoordinates};

pub struct Intersection {
    pub uv: UVCoordinates,
    pub intersection_point: Point,
    pub direction: FourVector,
    pub t: f64,
}

pub struct ColorComputationData {
    pub uv: UVCoordinates,
    pub temperature_data: TemperatureData,
    pub intersection_point: Point,
    pub direction: FourVector,
    /// Conserved per-ray frequency data; volumetric emitters use it to derive
    /// a per-sample redshift instead of the single per-intersection value in
    /// `temperature_data.redshift`.
    pub frequency: RayFrequencyData,
}

pub trait Hittable: Sync {
    fn intersects(&self, y_start: &Point, y_end: &Point) -> Option<Intersection>;
    fn color_at_uv(
        &self,
        color_computation_data: &ColorComputationData,
        geometry: &dyn Geometry,
    ) -> CIETristimulus;
    fn energy_of_emitter(
        &self,
        geometry: &dyn Geometry,
        step: &Step,
    ) -> Result<f64, RaytracerError>;
    fn temperature_of_emitter(
        &self,
        point: &Point,
        geometry: &dyn Geometry,
    ) -> Result<f64, RaytracerError>;
}

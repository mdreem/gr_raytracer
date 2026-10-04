use crate::geometry::coordinates::atlas::{CartesianConvertible, Chart};
use crate::geometry::coordinates::euclidean::CartesianChart;
use crate::geometry::geometry::Geometry;
use crate::geometry::point::Point;
use crate::rendering::integrator::Step;
use crate::rendering::radiance::Radiance;
use crate::rendering::raytracer::RaytracerError;
use crate::rendering::redshift::RayFrequencyData;
use crate::rendering::texture::{TemperatureData, UVCoordinates};

pub struct Intersection<C: Chart> {
    pub uv: UVCoordinates,
    pub intersection_point: Point<C>,
    pub t: f64,
}

pub struct ColorComputationData<'a, C: Chart + CartesianConvertible> {
    pub uv: UVCoordinates,
    pub temperature_data: TemperatureData,
    /// Conserved per-ray frequency data; volumetric emitters use it to derive
    /// a per-sample redshift instead of the single per-intersection value in
    /// `temperature_data.redshift`.
    pub frequency: RayFrequencyData,
    pub remaining_steps: &'a [Step<C>],
}

pub trait Hittable<C: CartesianConvertible>: Sync {
    fn intersects(
        &self,
        y_start: &Step<C>,
        y_end: &Step<C>,
        geometry: &dyn Geometry<Chart = C>,
    ) -> Option<Intersection<CartesianChart>>; // TODO: Should be the natural chart in the future
    /// Integrated light and background transmittance. Surface implementations
    /// convert straight texture colors here; volumes return their marched light.
    fn color_at_uv(
        &self,
        color_computation_data: &ColorComputationData<C>,
        geometry: &dyn Geometry<Chart = C>,
    ) -> Result<Radiance, RaytracerError>;
    fn energy_of_emitter(
        &self,
        geometry: &dyn Geometry<Chart = C>,
        step: &Step<C>,
    ) -> Result<f64, RaytracerError>;
    fn temperature_of_emitter(
        &self,
        point: &Point<C>,
        geometry: &dyn Geometry<Chart = C>,
    ) -> Result<f64, RaytracerError>;
}

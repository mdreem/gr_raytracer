use crate::geometry::geometry::Geometry;
use crate::geometry::point::Point;
use crate::rendering::color::CIETristimulus;
use crate::rendering::integrator::Step;
use crate::rendering::raytracer::RaytracerError;
use crate::rendering::redshift::RayFrequencyData;
use crate::rendering::texture::{TemperatureData, UVCoordinates};
use nalgebra::Vector3;

/// One step window of a traced ray, as the Cartesian spatial endpoints every
/// scene object works in.
///
/// The integrator runs in the geometry's native chart, so each object used to
/// redo the same spherical/Boyer-Lindquist -> Cartesian conversion (a sin/cos
/// pair per endpoint) for every object on every step of every ray. Handing
/// objects an already-converted segment, built from the step's cached
/// Cartesian position, makes that a single conversion per step.
#[derive(Debug, Clone, Copy)]
pub struct Segment {
    pub start_cartesian: Vector3<f64>,
    pub end_cartesian: Vector3<f64>,
}

impl Segment {
    /// Build a segment from two steps, reusing their cached Cartesian
    /// positions.
    pub fn from_steps(start: &Step, end: &Step) -> Segment {
        Segment {
            start_cartesian: start.spatial_cartesian(),
            end_cartesian: end.spatial_cartesian(),
        }
    }

    /// Build a segment from two bare points, converting them here. Used where
    /// no `Step` is at hand.
    #[allow(dead_code)] // Used by the scene objects' tests.
    pub fn from_points(start: &Point, end: &Point) -> Segment {
        Segment {
            start_cartesian: start.get_spatial_vector_cartesian(),
            end_cartesian: end.get_spatial_vector_cartesian(),
        }
    }
}

pub struct Intersection {
    pub uv: UVCoordinates,
    pub intersection_point: Point,
    pub t: f64,
}

pub struct ColorComputationData<'a> {
    pub uv: UVCoordinates,
    pub temperature_data: TemperatureData,
    /// Conserved per-ray frequency data; volumetric emitters use it to derive
    /// a per-sample redshift instead of the single per-intersection value in
    /// `temperature_data.redshift`.
    pub frequency: RayFrequencyData,
    pub remaining_steps: &'a [Step],
}

pub trait Hittable: Sync {
    fn intersects(&self, segment: &Segment, geometry: &dyn Geometry) -> Option<Intersection>;
    fn color_at_uv(
        &self,
        color_computation_data: &ColorComputationData,
        geometry: &dyn Geometry,
    ) -> Result<CIETristimulus, RaytracerError>;
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

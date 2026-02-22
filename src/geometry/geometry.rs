use crate::cli::cli::GlobalOpts;
use crate::configuration::RenderConfig;
use crate::geometry::four_vector::FourVector;
use crate::geometry::point::{CoordinateSystem, Point};
use crate::geometry::tetrad::Tetrad;
use crate::rendering::ray::Ray;
use crate::rendering::raytracer::RaytracerError;
use crate::rendering::runge_kutta::OdeFunction;
use crate::rendering::scene::EquationOfMotionState;
use crate::rendering::temperature::TemperatureComputer;
use nalgebra::{Const, Matrix4};
use std::io::Write;

pub trait GeodesicSolver: OdeFunction<Const<8>> + HasCoordinateSystem {
    fn geodesic(&self, t: f64, y: &EquationOfMotionState) -> EquationOfMotionState;
    fn create_initial_state(&self, ray: &Ray) -> EquationOfMotionState {
        EquationOfMotionState::from_column_slice(&[
            ray.position[0],
            ray.position[1],
            ray.position[2],
            ray.position[3],
            ray.momentum[0],
            ray.momentum[1],
            ray.momentum[2],
            ray.momentum[3],
        ])
    }
    fn momentum_from_state(&self, y: &EquationOfMotionState) -> FourVector {
        FourVector::new(y[4], y[5], y[6], y[7], self.coordinate_system())
    }
}

pub trait HasCoordinateSystem {
    fn coordinate_system(&self) -> CoordinateSystem;
}

pub trait InnerProduct {
    fn inner_product(&self, position: &Point, v: &FourVector, w: &FourVector) -> f64;
}

pub trait Signature {
    fn signature(&self) -> [f64; 4];
}

/// Support quantities that can be derived from the geometry.
/// This includes things like computing various velocities at a given position or ways to
/// compute temperatures.
pub trait SupportQuantities {
    fn get_stationary_velocity_at(&self, position: &Point) -> FourVector;
    fn get_circular_orbit_velocity_at(
        &self,
        position: &Point,
    ) -> Result<FourVector, RaytracerError>;
    fn get_temperature_computer(
        &self,
        temperature: f64,
        inner_radius: f64,
        outer_radius: f64,
    ) -> Result<Box<dyn TemperatureComputer>, RaytracerError>;
}

#[derive(Debug, Clone, Default)]
#[allow(dead_code)]
pub struct ConstantsOfMotion {
    constants: Vec<(&'static str, f64)>,
}

#[allow(dead_code)]
impl ConstantsOfMotion {
    pub fn push(&mut self, name: &'static str, value: f64) {
        self.constants.push((name, value));
    }

    pub fn as_slice(&self) -> &[(&'static str, f64)] {
        &self.constants
    }
}

pub trait Geometry:
    InnerProduct + HasCoordinateSystem + Signature + SupportQuantities + Sync
{
    fn get_tetrad_at(&self, position: &Point) -> Tetrad;
    fn lorentz_transformation(&self, position: &Point, velocity: &FourVector) -> Matrix4<f64>;
    fn inside_horizon(&self, position: &Point) -> bool;
    fn closed_orbit(&self, position: &Point, step_index: usize, max_steps: usize) -> bool;
    fn get_geodesic_solver(&self, ray: &Ray) -> Box<dyn GeodesicSolver>;
    fn get_radial_coordinate(&self, position: &Point) -> f64;
    /// Return the invariants along a geodesic, e.g. E and L_z.
    #[allow(dead_code)]
    fn get_constants_of_motion(&self, position: &Point, momentum: &FourVector)
    -> ConstantsOfMotion;
}

pub trait RenderableGeometry: Geometry {
    fn render(
        &self,
        opts: GlobalOpts,
        config: RenderConfig,
        camera_position: Point,
        filename: String,
        from_row: Option<u32>,
        from_col: Option<u32>,
        to_row: Option<u32>,
        to_col: Option<u32>,
    ) -> Result<(), RaytracerError>;

    fn render_ray(
        &self,
        row: i64,
        col: i64,
        opts: GlobalOpts,
        config: RenderConfig,
        camera_position: Point,
        write: &mut dyn Write,
    ) -> Result<(), RaytracerError>;

    fn render_ray_at(
        &self,
        position: Point,
        direction: FourVector,
        opts: GlobalOpts,
        write: &mut dyn Write,
    ) -> Result<(), RaytracerError>;
}

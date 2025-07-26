use crate::geometry::four_vector::{CoordinateSystem, FourVector};
use crate::rendering::runge_kutta::OdeFunction;
use crate::rendering::scene::EquationOfMotionState;
use nalgebra::{Const, Matrix4, Vector4};

#[derive(Debug, Clone)]
pub struct Tetrad {
    position: Vector4<f64>,
    pub t: FourVector,
    pub x: FourVector, // vertical wrt. the camera.
    pub y: FourVector, // horizontal wrt. the camera.
    pub z: FourVector, // away from the camera.
}

impl Tetrad {
    pub fn new(
        position: Vector4<f64>,
        t: FourVector,
        x: FourVector,
        y: FourVector,
        z: FourVector,
    ) -> Self {
        Tetrad {
            position,
            t,
            x,
            y,
            z,
        }
    }
}

pub trait GeodesicSolver {
    fn geodesic(&self, t: f64, y: &EquationOfMotionState) -> EquationOfMotionState;
}

pub trait HasCoordinateSystem {
    fn coordinate_system(&self) -> CoordinateSystem;
}

pub trait InnerProduct {
    fn inner_product(&self, position: &Vector4<f64>, v: &FourVector, w: &FourVector) -> f64;
}

pub trait Geometry:
    GeodesicSolver + InnerProduct + HasCoordinateSystem + Clone + Sync + OdeFunction<Const<8>>
{
    fn get_tetrad_at(&self, position: &Vector4<f64>) -> Tetrad;
    fn lorentz_transformation(
        &self,
        position: &Vector4<f64>,
        velocity: &FourVector,
    ) -> Matrix4<f64>;
    fn get_stationary_velocity_at(&self, position: &Vector4<f64>) -> FourVector;
    fn inside_horizon(&self, position: &Vector4<f64>) -> bool;
}

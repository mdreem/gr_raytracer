use crate::four_vector::{CoordinateSystem, FourVector};
use crate::runge_kutta::OdeFunction;
use crate::scene::EquationOfMotionState;
use nalgebra::{Const, Matrix4, Vector4};

#[derive(Debug)]
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

pub trait Geometry: Clone + Sync + OdeFunction<Const<8>> + HasCoordinateSystem {
    fn geodesic(&self, t: f64, y: &EquationOfMotionState) -> EquationOfMotionState;
    fn get_tetrad_at(&self, position: &Vector4<f64>) -> Tetrad;
    fn lorentz_transformation(
        &self,
        position: &Vector4<f64>,
        velocity: &FourVector,
    ) -> Matrix4<f64>;
}

pub trait HasCoordinateSystem {
    fn coordinate_system(&self) -> CoordinateSystem;
}

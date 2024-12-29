use crate::four_vector::FourVector;
use crate::runge_kutta::OdeFunction;
use crate::scene::EquationOfMotionState;
use nalgebra::{Const, Vector4};

pub struct Tetrad {
    position: Vector4<f64>,
    pub x: FourVector, // vertical wrt. the camera.
    pub y: FourVector, // horizontal wrt. the camera.
    pub z: FourVector, // away from the camera.
}

impl Tetrad {
    pub fn new(position: Vector4<f64>, x: FourVector, y: FourVector, z: FourVector) -> Self {
        Tetrad { position, x, y, z }
    }
}

pub trait Geometry: Clone + Sync + OdeFunction<Const<8>> {
    fn geodesic(&self, t: f64, y: &EquationOfMotionState) -> EquationOfMotionState;
    fn get_tetrad_at(&self, position: &Vector4<f64>) -> Tetrad;
}

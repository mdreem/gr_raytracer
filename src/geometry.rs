use crate::four_vector::FourVector;
use crate::runge_kutta::OdeFunction;
use crate::scene::EquationOfMotionState;
use nalgebra::{Const, OVector, Vector4};

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

#[derive(Clone)]
pub struct EuclideanSpace {}

impl EuclideanSpace {
    pub fn new() -> Self {
        EuclideanSpace {}
    }
}

impl OdeFunction<Const<8>> for EuclideanSpace {
    // TODO: maybe just have geodesic being used in solver. This doesn't need to be that generic here.
    fn apply(&self, t: f64, y: &OVector<f64, Const<8>>) -> OVector<f64, Const<8>> {
        self.geodesic(t, y)
    }
}

impl Geometry for EuclideanSpace {
    fn geodesic(&self, _: f64, y: &EquationOfMotionState) -> EquationOfMotionState {
        let y_new =
            EquationOfMotionState::from_column_slice(&[y[4], y[5], y[6], y[7], 0.0, 0.0, 0.0, 0.0]);
        y_new
    }

    // TODO: take into account Lorentz transformations.
    // TODO: take into account rotations.
    fn get_tetrad_at(&self, position: &Vector4<f64>) -> Tetrad {
        Tetrad::new(
            position.clone(),
            FourVector::new(0.0, 1.0, 0.0, 0.0),
            FourVector::new(0.0, 0.0, 1.0, 0.0),
            FourVector::new(0.0, 0.0, 0.0, 1.0),
        )
    }
}

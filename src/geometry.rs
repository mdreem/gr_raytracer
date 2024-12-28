use crate::runge_kutta::OdeFunction;
use crate::scene::EquationOfMotionState;
use nalgebra::{Const, OVector};

pub trait Geometry: Sync + OdeFunction<Const<8>> {
    fn geodesic(&self, t: f64, y: &EquationOfMotionState) -> EquationOfMotionState;
}

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
}

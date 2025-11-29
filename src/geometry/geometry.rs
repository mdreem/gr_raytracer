use crate::geometry::four_vector::FourVector;
use crate::geometry::point::{CoordinateSystem, Point};
use crate::geometry::tetrad::Tetrad;
use crate::rendering::ray::Ray;
use crate::rendering::runge_kutta::OdeFunction;
use crate::rendering::scene::EquationOfMotionState;
use nalgebra::{Const, Matrix4};

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

pub trait SupportQuantities {
    fn conserved_energy(&self, position: &Point, momentum: &FourVector) -> f64;
    fn conserved_angular_momentum(&self, position: &Point) -> f64;
    fn angular_velocity(&self, position: &Point) -> f64;
    fn get_stationary_velocity_at(&self, position: &Point) -> FourVector;
    fn get_circular_orbit_velocity_at(&self, position: &Point) -> FourVector;
}

pub trait Geometry:
    InnerProduct + HasCoordinateSystem + Signature + SupportQuantities + Sync
{
    fn get_tetrad_at(&self, position: &Point) -> Tetrad;
    fn lorentz_transformation(&self, position: &Point, velocity: &FourVector) -> Matrix4<f64>;
    fn inside_horizon(&self, position: &Point) -> bool;
    fn closed_orbit(&self, position: &Point, step_index: usize, max_steps: usize) -> bool;
    fn get_geodesic_solver(&self, ray: &Ray) -> Box<dyn GeodesicSolver>;
}

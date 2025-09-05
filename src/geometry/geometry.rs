use crate::geometry::four_vector::FourVector;
use crate::geometry::point::{CoordinateSystem, Point};
use crate::rendering::ray::Ray;
use crate::rendering::runge_kutta::OdeFunction;
use crate::rendering::scene::EquationOfMotionState;
use nalgebra::{Const, Matrix4};
use std::fmt::{Display, Formatter};

#[derive(Debug, Clone)]
pub struct Tetrad {
    position: Point,
    pub t: FourVector,
    pub x: FourVector, // vertical wrt. the camera.
    pub y: FourVector, // horizontal wrt. the camera.
    pub z: FourVector, // away from the camera.
}

impl Tetrad {
    pub fn new(
        position: Point,
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

impl Display for Tetrad {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "Tetrad\n")?;
        write!(f, "  t: {:?}\n", self.t)?;
        write!(f, "  x: {:?}\n", self.x)?;
        write!(f, "  y: {:?}\n", self.y)?;
        write!(f, "  z: {:?}\n", self.z)?;
        Ok(())
    }
}

pub trait GeodesicSolver: OdeFunction<Const<8>> {
    fn geodesic(&self, t: f64, y: &EquationOfMotionState) -> EquationOfMotionState;
}

pub trait HasCoordinateSystem {
    fn coordinate_system(&self) -> CoordinateSystem;
}

pub trait InnerProduct {
    fn inner_product(&self, position: &Point, v: &FourVector, w: &FourVector) -> f64;
}

pub trait Geometry: InnerProduct + HasCoordinateSystem + Clone + Sync {
    fn get_tetrad_at(&self, position: &Point) -> Tetrad;
    fn lorentz_transformation(&self, position: &Point, velocity: &FourVector) -> Matrix4<f64>;
    fn get_stationary_velocity_at(&self, position: &Point) -> FourVector;
    fn inside_horizon(&self, position: &Point) -> bool;
    fn get_geodesic_solver(&self, ray: &Ray) -> Box<dyn GeodesicSolver>;
}

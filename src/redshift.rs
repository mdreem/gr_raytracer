use crate::geometry::four_vector::FourVector;
use crate::geometry::geometry::Geometry;
use crate::ray::Ray;
use crate::scene::EquationOfMotionState;
use nalgebra::Vector4;

pub struct RedshiftComputer<'a, G: Geometry> {
    geometry: &'a G,
}

impl<'a, G: Geometry> RedshiftComputer<'a, G> {
    pub fn new(geometry: &'a G) -> Self {
        Self { geometry }
    }

    pub fn compute_redshift(&self, y: EquationOfMotionState, observer_energy: f64) -> f64 {
        let emitter_energy = self.energy_of_stationary_emitter(y);
        emitter_energy / observer_energy
    }

    pub fn get_observer_energy(&self, ray: &Ray, velocity: &FourVector) -> f64 {
        self.geometry
            .inner_product(&ray.position, velocity, &ray.momentum)
    }

    fn energy_of_stationary_emitter(&self, y: EquationOfMotionState) -> f64 {
        let position = Vector4::new(y[0], y[1], y[2], y[3]);
        let velocity = self.geometry.get_stationary_velocity_at(&position);
        let momentum = FourVector::new(y[4], y[5], y[6], y[7], self.geometry.coordinate_system());
        self.geometry.inner_product(&position, &velocity, &momentum)
    }
}

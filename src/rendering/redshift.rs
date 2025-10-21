use crate::geometry::four_vector::FourVector;
use crate::geometry::geometry::Geometry;
use crate::rendering::integrator::Step;
use crate::rendering::ray::Ray;

pub struct RedshiftComputer<'a, G: Geometry> {
    geometry: &'a G,
}

impl<'a, G: Geometry> RedshiftComputer<'a, G> {
    pub fn new(geometry: &'a G) -> Self {
        Self { geometry }
    }

    pub fn compute_redshift(&self, step: &Step, observer_energy: f64) -> f64 {
        let emitter_energy = self.energy_of_stationary_emitter(step);
        observer_energy / emitter_energy
    }

    pub fn get_observer_energy(&self, ray: &Ray, velocity: &FourVector) -> f64 {
        self.geometry
            .inner_product(&ray.position, velocity, &ray.momentum)
    }

    fn energy_of_stationary_emitter(&self, step: &Step) -> f64 {
        let position = step.x;
        let velocity = self.geometry.get_stationary_velocity_at(&position);
        let momentum = step.p;
        self.geometry.inner_product(&position, &velocity, &momentum)
    }
}

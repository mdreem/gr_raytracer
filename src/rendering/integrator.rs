use crate::geometry::geometry::Geometry;
use crate::geometry::point::Point;
use crate::rendering::integrator::StopReason::{CelestialSphereReached, HorizonReached};
use crate::rendering::ray::{IntegratedRay, Ray};
use crate::rendering::raytracer::RaytracerError;
use crate::rendering::runge_kutta::rkf45;
use crate::rendering::scene::{EquationOfMotionState, get_position};

#[derive(Debug)]
pub struct Step {
    pub y: EquationOfMotionState,
    pub t: f64,
    pub step: usize,
}

#[derive(Debug, PartialEq)]
pub enum StopReason {
    HorizonReached,
    CelestialSphereReached,
}

pub struct Integrator<'a, G: Geometry> {
    integration_configuration: IntegrationConfiguration,
    geometry: &'a G,
}

impl<'a, G: Geometry> Integrator<'a, G> {
    pub fn new(
        geometry: &'a G,
        integration_configuration: IntegrationConfiguration,
    ) -> Integrator<'a, G> {
        Integrator {
            integration_configuration,
            geometry,
        }
    }
}

pub struct IntegrationConfiguration {
    max_steps: usize,
    max_radius_sq: f64,
    step_size: f64,
    epsilon: f64,
}

impl IntegrationConfiguration {
    pub fn new(
        max_steps: usize,
        max_radius: f64,
        step_size: f64,
        epsilon: f64,
    ) -> IntegrationConfiguration {
        IntegrationConfiguration {
            max_steps,
            max_radius_sq: max_radius * max_radius,
            step_size,
            epsilon,
        }
    }
}

#[derive(Debug)]
pub enum IntegrationError {
    MaxStepsReached,
}

impl<G: Geometry> Integrator<'_, G> {
    pub fn integrate(
        &self,
        ray: &Ray,
    ) -> Result<(IntegratedRay, Option<StopReason>), RaytracerError> {
        let mut t = 0.0;
        let direction = ray.momentum.get_as_vector();
        let mut y = EquationOfMotionState::from_column_slice(&[
            ray.position[0],
            ray.position[1],
            ray.position[2],
            ray.position[3],
            direction[0],
            direction[1],
            direction[2],
            direction[3],
        ]);

        let mut result: Vec<Step> = Vec::with_capacity(self.integration_configuration.max_steps);
        result.push(Step { y, t, step: 0 });

        let mut h = self.integration_configuration.step_size;
        for i in 1..self.integration_configuration.max_steps {
            let last_y = y;
            (y, h) = rkf45(
                &y,
                t,
                h,
                self.integration_configuration.epsilon,
                self.geometry,
            )?;
            t += h;

            result.push(Step { y, t, step: i });
            match self.should_stop(&last_y, &y) {
                None => {}
                Some(r) => return Ok((IntegratedRay::new(result), Some(r))),
            }
        }

        Ok((IntegratedRay::new(result), None))
    }

    fn should_stop(
        &self,
        _last_y: &EquationOfMotionState,
        cur_y: &EquationOfMotionState,
    ) -> Option<StopReason> {
        if self.geometry.inside_horizon(&Point::new(
            cur_y[0],
            cur_y[1],
            cur_y[2],
            cur_y[3],
            self.geometry.coordinate_system(),
        )) {
            return Some(HorizonReached);
        }

        // iterate until the celestial plane distance has been reached.
        if get_position(cur_y, self.geometry.coordinate_system())
            .radial_distance_spatial_part_squared()
            > self.integration_configuration.max_radius_sq
        {
            return Some(CelestialSphereReached);
        }

        None
    }
}

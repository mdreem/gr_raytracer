use crate::geometry::four_vector::FourVector;
use crate::geometry::geometry::Geometry;
use crate::geometry::point::Point;
use crate::rendering::integrator::StopReason::{
    CelestialSphereReached, CoordinateIsNan, HorizonReached,
};
use crate::rendering::ray::{IntegratedRay, Ray};
use crate::rendering::raytracer::RaytracerError;
use crate::rendering::runge_kutta::rkf45;
use crate::rendering::scene::{EquationOfMotionState, get_position};
use log::debug;

#[derive(Debug)]
pub struct Step {
    pub x: Point,
    pub p: FourVector,
    pub t: f64,
    pub step: usize,
}

#[derive(Debug, PartialEq)]
pub enum StopReason {
    HorizonReached,
    CelestialSphereReached,
    CoordinateIsNan,
    ClosedOrbitDetected,
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

#[derive(Debug, thiserror::Error)]
pub enum IntegrationError {
    #[error("Max steps reached")]
    MaxStepsReached,
    #[error("No steps produced")]
    NoStepsProduced,
}

impl<G: Geometry> Integrator<'_, G> {
    pub fn integrate(
        &self,
        ray: &Ray,
    ) -> Result<(IntegratedRay, Option<StopReason>), RaytracerError> {
        let mut t = 0.0;
        let geodesic_solver = self.geometry.get_geodesic_solver(ray);
        let mut y = geodesic_solver.create_initial_state(ray);

        let mut result: Vec<Step> = Vec::with_capacity(self.integration_configuration.max_steps);
        let x = Point::new(y[0], y[1], y[2], y[3], self.geometry.coordinate_system());
        let p = geodesic_solver.momentum_from_state(&y);
        result.push(Step { x, p, t, step: 0 });

        #[cfg(debug_assertions)]
        let initial_constants = self.geometry.get_constants_of_motion(&x, &p);

        #[cfg(debug_assertions)]
        let mut max_k_dot_k_drift = 0.0f64;
        #[cfg(debug_assertions)]
        let mut max_constant_drifts = vec![0.0f64; initial_constants.as_slice().len()];

        let mut h = self.integration_configuration.step_size;
        for i in 1..self.integration_configuration.max_steps {
            let last_y = y;
            let (y_new, h_taken, h_next) = rkf45(
                &y,
                t,
                h,
                self.integration_configuration.epsilon,
                &*geodesic_solver,
            )?;
            y = y_new;
            // Advance the affine parameter by the step actually taken, and carry
            // the controller's suggestion into the next iteration.
            t += h_taken;
            h = h_next;

            let x = Point::new(y[0], y[1], y[2], y[3], self.geometry.coordinate_system());
            let p = geodesic_solver.momentum_from_state(&y);
            result.push(Step { x, p, t, step: i });

            #[cfg(debug_assertions)]
            {
                let k_dot_k = self.geometry.inner_product(&x, &p, &p);
                // For a photon, k.k should always be exactly 0.0
                let k_drift = k_dot_k.abs();
                if k_drift > max_k_dot_k_drift {
                    max_k_dot_k_drift = k_drift;
                }

                let current_constants = self.geometry.get_constants_of_motion(&x, &p);

                for (idx, ((_, initial_val), (_, current_val))) in initial_constants
                    .as_slice()
                    .iter()
                    .zip(current_constants.as_slice().iter())
                    .enumerate()
                {
                    // Relative drift, falling back to absolute if initial value is 0
                    let drift = if initial_val.abs() > 1e-12 {
                        (current_val - initial_val).abs() / initial_val.abs()
                    } else {
                        (current_val - initial_val).abs()
                    };
                    if drift > max_constant_drifts[idx] {
                        max_constant_drifts[idx] = drift;
                    }
                }
            }

            match self.should_stop(&last_y, &y, i) {
                None => {}
                Some(r) => {
                    self.report_drifts(
                        #[cfg(debug_assertions)]
                        max_k_dot_k_drift,
                        #[cfg(debug_assertions)]
                        &initial_constants,
                        #[cfg(debug_assertions)]
                        &max_constant_drifts,
                    );
                    return Ok((IntegratedRay::new(result), Some(r)));
                }
            }
        }

        self.report_drifts(
            #[cfg(debug_assertions)]
            max_k_dot_k_drift,
            #[cfg(debug_assertions)]
            &initial_constants,
            #[cfg(debug_assertions)]
            &max_constant_drifts,
        );

        Ok((IntegratedRay::new(result), None))
    }

    fn report_drifts(
        &self,
        #[cfg(debug_assertions)] max_k_dot_k_drift: f64,
        #[cfg(debug_assertions)] initial_constants: &crate::geometry::geometry::ConstantsOfMotion,
        #[cfg(debug_assertions)] max_constant_drifts: &[f64],
    ) {
        #[cfg(debug_assertions)]
        {
            let tolerance = 1e-4;
            if max_k_dot_k_drift > tolerance {
                log::warn!(
                    "Large invariant drift detected! k.k drift: {:.3e}",
                    max_k_dot_k_drift
                );
            }
            for (idx, (name, _)) in initial_constants.as_slice().iter().enumerate() {
                if max_constant_drifts[idx] > tolerance {
                    log::warn!(
                        "Large invariant drift detected! {} drift: {:.3e}",
                        name,
                        max_constant_drifts[idx]
                    );
                }
            }
        }
    }

    fn should_stop(
        &self,
        _last_y: &EquationOfMotionState,
        cur_y: &EquationOfMotionState,
        step_index: usize,
    ) -> Option<StopReason> {
        // Position must be finite to evaluate anything below; if it isn't,
        // there's no reliable state to check horizon/escape against.
        if !cur_y.iter().take(4).all(|v| v.is_finite()) {
            debug!("last_y: {:?}", _last_y);
            debug!(
                "spherical last_y: {:?}",
                get_position(_last_y, self.geometry.coordinate_system()).get_as_spherical()
            );
            return Some(CoordinateIsNan);
        }

        if self.geometry.inside_horizon(&Point::new(
            cur_y[0],
            cur_y[1],
            cur_y[2],
            cur_y[3],
            self.geometry.coordinate_system(),
        )) {
            return Some(HorizonReached);
        }

        if self.geometry.closed_orbit(
            &Point::new(
                cur_y[0],
                cur_y[1],
                cur_y[2],
                cur_y[3],
                self.geometry.coordinate_system(),
            ),
            step_index,
            self.integration_configuration.max_steps,
        ) {
            return Some(StopReason::ClosedOrbitDetected);
        }

        // iterate until the celestial plane distance has been reached.
        if get_position(cur_y, self.geometry.coordinate_system())
            .radial_distance_spatial_part_squared()
            > self.integration_configuration.max_radius_sq
        {
            return Some(CelestialSphereReached);
        }

        // Momentum can blow up at coordinate singularities (e.g. the BL polar
        // axis) even when the position stays finite for one more step. Check
        // this only after the position-based conditions above, so a ray that
        // has already escaped past max_radius or crossed the horizon isn't
        // misclassified as a failure just because its momentum diverged
        // afterward (e.g. from grazing the polar axis far from the origin).
        if !cur_y.iter().skip(4).all(|v| v.is_finite()) {
            debug!("last_y: {:?}", _last_y);
            debug!(
                "spherical last_y: {:?}",
                get_position(_last_y, self.geometry.coordinate_system()).get_as_spherical()
            );
            return Some(CoordinateIsNan);
        }

        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::schwarzschild::Schwarzschild;

    #[test]
    fn test_should_stop_prefers_celestial_sphere_over_nonfinite_momentum() {
        // A ray whose position has already escaped past max_radius but whose
        // momentum has since diverged (e.g. from grazing a coordinate
        // singularity far from the origin) must be classified as having
        // reached the celestial sphere, not as a NaN failure.
        let geometry = Schwarzschild::new(2.0, 1e-5);
        let integration_configuration = IntegrationConfiguration::new(100, 100.0, 0.01, 1e-5);
        let integrator = Integrator::new(&geometry, integration_configuration);

        let last_y =
            EquationOfMotionState::from_column_slice(&[0.0, 200.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0]);
        let cur_y = EquationOfMotionState::from_column_slice(&[
            0.0,
            200.0, // position: well past max_radius = 100
            1.0,
            0.0,
            1.0,
            1.0,
            0.0,
            f64::INFINITY, // momentum component diverged
        ]);

        assert_eq!(
            integrator.should_stop(&last_y, &cur_y, 1),
            Some(StopReason::CelestialSphereReached)
        );
    }

    #[test]
    fn test_should_stop_detects_nonfinite_momentum_when_not_escaped() {
        // If the ray hasn't escaped/crossed the horizon, non-finite momentum
        // must still be caught.
        let geometry = Schwarzschild::new(2.0, 1e-5);
        let integration_configuration = IntegrationConfiguration::new(100, 100.0, 0.01, 1e-5);
        let integrator = Integrator::new(&geometry, integration_configuration);

        let last_y =
            EquationOfMotionState::from_column_slice(&[0.0, 10.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0]);
        let cur_y = EquationOfMotionState::from_column_slice(&[
            0.0,
            10.0, // position: well within max_radius = 100
            1.0,
            0.0,
            1.0,
            1.0,
            0.0,
            f64::NAN,
        ]);

        assert_eq!(
            integrator.should_stop(&last_y, &cur_y, 1),
            Some(StopReason::CoordinateIsNan)
        );
    }
}

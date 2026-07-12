use crate::geometry::four_vector::FourVector;
use crate::geometry::geometry::Geometry;
use crate::rendering::integrator::Step;
use crate::rendering::ray::Ray;

/// Per-ray frequency data: the observer energy plus the photon's conserved
/// covariant components along the Killing vectors. All three are
/// chart-independent scalars, computed once per ray at the camera (where the
/// exact, un-interpolated momentum is available) and valid anywhere along the
/// geodesic by conservation. Together with a local emitter's Killing
/// coefficients they give the redshift at any sample point without parallel
/// transport:
///
///     g = observer_energy / (u^t * p_t + u^phi * p_phi)
#[derive(Debug, Clone, Copy)]
pub struct RayFrequencyData {
    pub observer_energy: f64,
    pub p_t: f64,
    pub p_phi: f64,
}

pub struct RedshiftComputer<'a, G: Geometry> {
    geometry: &'a G,
}

impl<'a, G: Geometry> RedshiftComputer<'a, G> {
    pub fn new(geometry: &'a G) -> Self {
        Self { geometry }
    }

    pub fn compute_redshift(&self, step: &Step, observer_energy: f64) -> f64 {
        let emitter_energy = self.energy_of_stationary_emitter(step);
        self.compute_redshift_from_energies(emitter_energy, observer_energy)
    }

    pub fn compute_redshift_from_energies(&self, emitter_energy: f64, observer_energy: f64) -> f64 {
        self.to_physical_energy(observer_energy) / self.to_physical_energy(emitter_energy)
    }

    pub fn get_observer_energy(&self, ray: &Ray, velocity: &FourVector) -> f64 {
        self.geometry
            .inner_product(&ray.position, velocity, &ray.momentum)
    }

    pub fn get_ray_frequency_data(&self, ray: &Ray, velocity: &FourVector) -> RayFrequencyData {
        let observer_energy = self.get_observer_energy(ray, velocity);
        let e_t = FourVector::new(1.0, 0.0, 0.0, 0.0, ray.momentum.coordinate_system);
        let p_t = self
            .geometry
            .inner_product(&ray.position, &e_t, &ray.momentum);
        let axial = self.geometry.axial_killing_vector(&ray.position);
        let p_phi = self
            .geometry
            .inner_product(&ray.position, &axial, &ray.momentum);
        RayFrequencyData {
            observer_energy,
            p_t,
            p_phi,
        }
    }

    fn energy_of_stationary_emitter(&self, step: &Step) -> f64 {
        let position = step.x;
        let velocity = self.geometry.get_stationary_velocity_at(&position);
        let momentum = step.p;
        self.geometry.inner_product(&position, &velocity, &momentum)
    }

    fn to_physical_energy(&self, inner_product_energy: f64) -> f64 {
        let energy = self.geometry.signature()[0] * inner_product_energy;
        debug_assert!(energy.is_finite());
        energy
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::euclidean::EuclideanSpace;
    use crate::geometry::four_vector::FourVector;
    use crate::geometry::geometry::{InnerProduct, SupportQuantities};
    use crate::geometry::point::Point;
    use crate::geometry::schwarzschild::Schwarzschild;
    use crate::rendering::camera::Camera;
    use crate::rendering::integrator::{IntegrationConfiguration, Integrator};
    use approx::assert_abs_diff_eq;
    use std::f64::consts::PI;

    const EMITTER_SPEED: f64 = 0.5;

    fn gamma(v: f64) -> f64 {
        (1.0 - v * v).sqrt().recip()
    }

    // ------------------------------------------------------------------
    // Pure RedshiftComputer tests: hand-built photon momenta, no Camera
    // and no Integrator involved.
    // ------------------------------------------------------------------

    /// Flat space, hand-built traced momentum. The observer sits at +x, the
    /// emitter towards -x, so the traced (past-directed) photon momentum is
    /// N - e_t with N = -x_hat: p = (-1, -1, 0, 0).
    fn doppler_redshift_for(emitter_velocity: &FourVector) -> f64 {
        let geometry = EuclideanSpace::new();
        let computer = RedshiftComputer::new(&geometry);
        let position = Point::new_cartesian(0.0, 0.0, 0.0, 0.0);
        let momentum = FourVector::new_cartesian(-1.0, -1.0, 0.0, 0.0);

        let observer = FourVector::new_cartesian(1.0, 0.0, 0.0, 0.0);
        let observer_energy = geometry.inner_product(&position, &observer, &momentum);
        let emitter_energy = geometry.inner_product(&position, emitter_velocity, &momentum);
        computer.compute_redshift_from_energies(emitter_energy, observer_energy)
    }

    #[test]
    fn test_doppler_formulas_with_hand_built_momentum() {
        let v = EMITTER_SPEED;
        // Approaching the observer at +x: nu_obs/nu_em = 1 / (gamma (1 - v)).
        let approaching = FourVector::new_cartesian(gamma(v), gamma(v) * v, 0.0, 0.0);
        assert_abs_diff_eq!(
            doppler_redshift_for(&approaching),
            (gamma(v) * (1.0 - v)).recip(),
            epsilon = 1e-12
        );
        // Receding: nu_obs/nu_em = 1 / (gamma (1 + v)).
        let receding = FourVector::new_cartesian(gamma(v), -gamma(v) * v, 0.0, 0.0);
        assert_abs_diff_eq!(
            doppler_redshift_for(&receding),
            (gamma(v) * (1.0 + v)).recip(),
            epsilon = 1e-12
        );
        // Transverse: only time dilation, nu_obs/nu_em = 1 / gamma.
        let transverse = FourVector::new_cartesian(gamma(v), 0.0, gamma(v) * v, 0.0);
        assert_abs_diff_eq!(
            doppler_redshift_for(&transverse),
            gamma(v).recip(),
            epsilon = 1e-12
        );
        // At rest: no shift.
        let at_rest = FourVector::new_cartesian(1.0, 0.0, 0.0, 0.0);
        assert_abs_diff_eq!(doppler_redshift_for(&at_rest), 1.0, epsilon = 1e-12);
    }

    /// The redshift is a ratio of two inner products with the same transported
    /// momentum, so it must not depend on the overall sign (time orientation)
    /// of that momentum.
    #[test]
    fn test_redshift_is_invariant_under_momentum_negation() {
        let geometry = EuclideanSpace::new();
        let computer = RedshiftComputer::new(&geometry);
        let position = Point::new_cartesian(0.0, 0.0, 0.0, 0.0);
        let observer = FourVector::new_cartesian(1.0, 0.0, 0.0, 0.0);
        let v = EMITTER_SPEED;
        let emitter = FourVector::new_cartesian(gamma(v), gamma(v) * v, 0.0, 0.0);

        let momentum = FourVector::new_cartesian(-1.0, -1.0, 0.0, 0.0);
        let negated = -momentum;

        let g = computer.compute_redshift_from_energies(
            geometry.inner_product(&position, &emitter, &momentum),
            geometry.inner_product(&position, &observer, &momentum),
        );
        let g_negated = computer.compute_redshift_from_energies(
            geometry.inner_product(&position, &emitter, &negated),
            geometry.inner_product(&position, &observer, &negated),
        );
        assert_abs_diff_eq!(g, g_negated, epsilon = 1e-15);
    }

    /// The code's disc redshift must reproduce Luminet (1979), eq. for 1+z of
    /// a circular-orbit emitter in Schwarzschild seen from infinity:
    ///   1 + z = (1 - 3M/r)^{-1/2} (1 + Omega * L/E)
    /// where E = a p^t and L = -r^2 p^phi are the photon's conserved energy
    /// and angular momentum (signature +---). This holds for *any* null
    /// momentum at the emitter, so no integration is needed: it checks the
    /// same inner products `Disc::energy_of_emitter` uses.
    #[test]
    fn test_disc_redshift_matches_luminet_closed_form() {
        let geometry = Schwarzschild::new(1.0, 1e-4);
        let m = 0.5; // r_s = 1

        for r in [2.0_f64, 3.0, 5.0, 10.0] {
            let a = 1.0 - 1.0 / r;
            let omega = (m / (r * r * r)).sqrt();
            let u_t = (1.0 - 3.0 * m / r).sqrt().recip();

            for phi in [0.3_f64, 2.0, 4.5] {
                let position = Point::new_spherical(0.0, r, PI / 2.0, phi);
                let u_emitter = geometry.get_circular_orbit_velocity_at(&position).unwrap();

                // A few null momenta with different directions, including
                // prograde and retrograde azimuthal components. Past-directed
                // (p^t < 0), matching the camera convention.
                for (p_r, p_phi_hat) in [(-0.7_f64, 0.4_f64), (-0.2, -0.9), (0.5, 0.6), (0.0, 1.0)]
                {
                    let p_phi = p_phi_hat / r; // scale to comparable magnitude
                    let spatial = p_r * p_r / a + r * r * p_phi * p_phi;
                    let p_t = -(spatial / a).sqrt();
                    let momentum = FourVector::new_spherical(p_t, p_r, 0.0, p_phi);
                    assert_abs_diff_eq!(
                        geometry.inner_product(&position, &momentum, &momentum),
                        0.0,
                        epsilon = 1e-12
                    );

                    // Conserved quantities of the photon.
                    let e_conserved = a * p_t;
                    let l_conserved = -r * r * p_phi;

                    // Code path: the emitter energy exactly as
                    // Disc::energy_of_emitter computes it.
                    let emitter_energy =
                        geometry.inner_product(&position, &u_emitter, &momentum);

                    // Closed form for the same contraction.
                    assert_abs_diff_eq!(
                        emitter_energy,
                        u_t * (e_conserved + omega * l_conserved),
                        epsilon = 1e-10
                    );

                    // Luminet's 1+z against the code's g = nu_obs/nu_em for an
                    // observer at infinity (u_obs . p = E_conserved there).
                    let g_code = e_conserved / emitter_energy;
                    let g_luminet = (1.0 - 3.0 * m / r).sqrt()
                        / (1.0 + omega * l_conserved / e_conserved);
                    assert_abs_diff_eq!(g_code, g_luminet, epsilon = 1e-10);
                }
            }
        }
    }

    /// Gravitational redshift without any integration: for a radial null
    /// geodesic in Schwarzschild, p_t is conserved, so the momentum at any
    /// radius is (p_t / a, p_t, 0, 0) up to the radial sign. A stationary
    /// emitter deeper in the potential must be redshifted by
    /// sqrt(a_em / a_cam).
    #[test]
    fn test_gravitational_redshift_analytic_without_integration() {
        let geometry = Schwarzschild::new(1.0, 1e-4);
        let computer = RedshiftComputer::new(&geometry);

        let r_camera: f64 = 10.0;
        let r_emitter: f64 = 3.0;
        let a_camera = 1.0 - 1.0 / r_camera;
        let a_emitter = 1.0 - 1.0 / r_emitter;

        // Past-directed trace marching inward: p_t < 0, p^r < 0.
        let p_t = -1.0;
        let camera_position = Point::new_spherical(0.0, r_camera, PI / 2.0, 0.0);
        let emitter_position = Point::new_spherical(0.0, r_emitter, PI / 2.0, 0.0);
        let p_camera = FourVector::new_spherical(p_t / a_camera, p_t, 0.0, 0.0);
        let p_emitter = FourVector::new_spherical(p_t / a_emitter, p_t, 0.0, 0.0);

        // Both ends must be null.
        assert_abs_diff_eq!(
            geometry.inner_product(&camera_position, &p_camera, &p_camera),
            0.0,
            epsilon = 1e-12
        );
        assert_abs_diff_eq!(
            geometry.inner_product(&emitter_position, &p_emitter, &p_emitter),
            0.0,
            epsilon = 1e-12
        );

        let u_camera = geometry.get_stationary_velocity_at(&camera_position);
        let observer_energy = geometry.inner_product(&camera_position, &u_camera, &p_camera);

        let step = Step {
            x: emitter_position,
            p: p_emitter,
            t: 1.0,
            step: 1,
        };
        let redshift = computer.compute_redshift(&step, observer_energy);
        assert_abs_diff_eq!(redshift, (a_emitter / a_camera).sqrt(), epsilon = 1e-12);
    }

    /// Checks 2 and 3 of docs/redshift-doppler-sign-review.md: flat space,
    /// camera at x = +10 looking at the origin (center pixel), emitter on the
    /// line of sight at x = +5 with four-velocity `emitter_velocity`. Parallel
    /// transport is trivial in flat space, so the transported ray momentum at
    /// the emitter is the camera-frame momentum itself. The returned value is
    /// the code's `redshift` variable, i.e. nu_obs / nu_em = 1 / (1 + z).
    fn flat_space_redshift_for(emitter_velocity: FourVector) -> f64 {
        let geometry = EuclideanSpace::new();
        let camera = Camera::new(
            Point::new_cartesian(0.0, 10.0, 0.0, 0.0),
            FourVector::new_cartesian(1.0, 0.0, 0.0, 0.0),
            PI / 2.0,
            11,
            11,
            0.0,
            0.0,
            0.0,
            &geometry,
        )
        .unwrap();

        let ray = camera.get_ray_for(5, 5);
        // The traced ray must march from the camera toward the emitter,
        // i.e. in the -x direction, whatever the time orientation is.
        assert!(ray.momentum[1] < 0.0);

        let computer = RedshiftComputer::new(&geometry);
        let observer_energy = computer.get_observer_energy(&ray, &camera.velocity);

        let emitter_position = Point::new_cartesian(0.0, 5.0, 0.0, 0.0);
        let emitter_energy =
            geometry.inner_product(&emitter_position, &emitter_velocity, &ray.momentum);
        computer.compute_redshift_from_energies(emitter_energy, observer_energy)
    }

    #[test]
    fn test_emitter_at_rest_in_flat_space_has_no_shift() {
        let redshift = flat_space_redshift_for(FourVector::new_cartesian(1.0, 0.0, 0.0, 0.0));
        assert_abs_diff_eq!(redshift, 1.0, epsilon = 1e-12);
    }

    #[test]
    fn test_emitter_moving_toward_camera_is_blueshifted() {
        let v = EMITTER_SPEED;
        // The camera sits at x = +10, the emitter at x = +5, so moving toward
        // the camera means a positive x velocity.
        let u = FourVector::new_cartesian(gamma(v), gamma(v) * v, 0.0, 0.0);
        let redshift = flat_space_redshift_for(u);
        // Longitudinal Doppler: nu_obs / nu_em = 1 / (gamma (1 - v)) > 1.
        assert_abs_diff_eq!(redshift, (gamma(v) * (1.0 - v)).recip(), epsilon = 1e-12);
    }

    #[test]
    fn test_emitter_moving_away_from_camera_is_redshifted() {
        let v = EMITTER_SPEED;
        let u = FourVector::new_cartesian(gamma(v), -gamma(v) * v, 0.0, 0.0);
        let redshift = flat_space_redshift_for(u);
        // Longitudinal Doppler: nu_obs / nu_em = 1 / (gamma (1 + v)) < 1.
        assert_abs_diff_eq!(redshift, (gamma(v) * (1.0 + v)).recip(), epsilon = 1e-12);
    }

    #[test]
    fn test_emitter_moving_transverse_shows_transverse_doppler_redshift() {
        let v = EMITTER_SPEED;
        // Motion perpendicular to the line of sight: only time dilation remains.
        let u = FourVector::new_cartesian(gamma(v), 0.0, gamma(v) * v, 0.0);
        let redshift = flat_space_redshift_for(u);
        assert_abs_diff_eq!(redshift, gamma(v).recip(), epsilon = 1e-12);
    }

    /// Pure gravitational redshift must be independent of the traced ray's
    /// time orientation: a stationary emitter and a stationary camera see
    /// nu_obs / nu_em = sqrt(a_emitter / a_camera).
    #[test]
    fn test_stationary_emitter_gravitational_redshift_matches_analytic() {
        let geometry = Schwarzschild::new(1.0, 1e-4);
        let r_camera: f64 = 10.0;
        let a_camera = 1.0 - 1.0 / r_camera;

        let position = Point::new_spherical(0.0, r_camera, PI / 2.0, 0.0);
        let velocity = FourVector::new_spherical(a_camera.sqrt().recip(), 0.0, 0.0, 0.0);
        let camera = Camera::new(
            position,
            velocity,
            PI / 2.0,
            11,
            11,
            0.0,
            0.0,
            0.0,
            &geometry,
        )
        .unwrap();
        let ray = camera.get_ray_for(5, 5);

        let integrator = Integrator::new(
            &geometry,
            IntegrationConfiguration::new(10000, 100.0, 0.001, 1e-10),
        );
        let (integrated, _) = integrator.integrate(&ray).unwrap();

        let computer = RedshiftComputer::new(&geometry);
        let observer_energy = computer.get_observer_energy(&ray, &camera.velocity);

        // Check the redshift against the analytic value at several points
        // along the trajectory (skipping the camera step itself).
        let mut checked = 0;
        for step in integrated.steps.iter().skip(50).step_by(200) {
            let r = step.x.get_as_spherical()[0];
            if r <= 1.0 + 1e-3 {
                continue;
            }
            let a_emitter = 1.0 - 1.0 / r;
            let redshift = computer.compute_redshift(step, observer_energy);
            assert_abs_diff_eq!(
                redshift,
                (a_emitter / a_camera).sqrt(),
                epsilon = 1e-6
            );
            checked += 1;
        }
        assert!(checked > 0, "no steps were checked");
    }
}

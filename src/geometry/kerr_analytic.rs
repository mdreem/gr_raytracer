use crate::geometry::four_vector::FourVector;
use crate::geometry::geometry::{
    GeodesicSolver, Geometry, HasCoordinateSystem, InnerProduct, Signature,
};
use crate::geometry::point::CoordinateSystem::Spherical;
use crate::geometry::point::{CoordinateSystem, Point};
use crate::geometry::tetrad::Tetrad;
use crate::rendering::ray::Ray;
use crate::rendering::runge_kutta::OdeFunction;
use crate::rendering::scene::EquationOfMotionState;
use log::{debug, trace};
use nalgebra::{Const, Matrix4, OVector};

/// Kerr black hole geometry using analytic geodesic equations in Boyer-Lindquist coordinates.
///
/// This implementation uses first-order geodesic equations with conserved quantities
/// (energy E, angular momentum L_z, and Carter constant Q) to compute photon trajectories.
#[derive(Clone, Debug)]
pub struct KerrAnalytic {
    /// Schwarzschild radius (2M in geometric units)
    radius: f64,
    /// Spin parameter
    a: f64,
    /// Small epsilon to detect horizon crossing
    horizon_epsilon: f64,
}

struct KerrAnalyticSolver {
    radius: f64,
    a: f64,
    /// Conserved energy E = -p_t
    energy: f64,
    /// Conserved angular momentum L_z = p_φ
    angular_momentum: f64,
    /// Carter constant Q
    carter_constant: f64,
}

impl KerrAnalytic {
    pub fn new(radius: f64, a: f64, horizon_epsilon: f64) -> Self {
        KerrAnalytic {
            radius,
            a,
            horizon_epsilon,
        }
    }

    /// Compute Σ = r² + a² cos²θ
    fn sigma(&self, r: f64, theta: f64) -> f64 {
        r * r + self.a * self.a * theta.cos() * theta.cos()
    }

    /// Compute Δ = r² - 2Mr + a²
    fn delta(&self, r: f64) -> f64 {
        r * r - self.radius * r + self.a * self.a
    }
}

impl OdeFunction<Const<8>> for KerrAnalyticSolver {
    fn apply(&self, t: f64, y: &OVector<f64, Const<8>>) -> OVector<f64, Const<8>> {
        self.geodesic(t, y)
    }
}

impl HasCoordinateSystem for KerrAnalyticSolver {
    fn coordinate_system(&self) -> CoordinateSystem {
        Spherical
    }
}

impl KerrAnalyticSolver {
    /// Compute Σ = r² + a² cos²θ
    fn sigma(&self, r: f64, theta: f64) -> f64 {
        r * r + self.a * self.a * theta.cos() * theta.cos()
    }

    /// Compute Δ = r² - 2Mr + a²
    fn delta(&self, r: f64) -> f64 {
        r * r - self.radius * r + self.a * self.a
    }

    /// Compute R = [(r² + a²)E - aL_z]² - Δ[(L_z - aE)² + Q]
    /// For photon geodesics (mass = 0)
    fn r_potential(&self, r: f64) -> f64 {
        let delta = self.delta(r);
        let term1 = (r * r + self.a * self.a) * self.energy - self.a * self.angular_momentum;
        let term2 = self.angular_momentum - self.a * self.energy;
        term1 * term1 - delta * (term2 * term2 + self.carter_constant)
    }

    /// Compute Θ = Q - cos²θ[a²(1 - E²) + L²/sin²θ]
    /// For photon geodesics (mass = 0, so 1 - E² → -E²)
    fn theta_potential(&self, theta: f64) -> f64 {
        let cos_theta = theta.cos();
        let sin_theta = theta.sin();
        if sin_theta.abs() < 1e-10 {
            return self.carter_constant;
        }
        let term = -self.a * self.a * self.energy * self.energy
            + self.angular_momentum * self.angular_momentum / (sin_theta * sin_theta);
        self.carter_constant - cos_theta * cos_theta * term
    }
}

impl GeodesicSolver for KerrAnalyticSolver {
    /// First-order geodesic equations in Boyer-Lindquist coordinates using conserved quantities.
    ///
    /// The equations are:
    /// Σ dt/dλ = (r² + a²)[(r² + a²)E - aL_z]/Δ - a[aE sin²θ - L_z]
    /// Σ dr/dλ = ±√R
    /// Σ dθ/dλ = ±√Θ
    /// Σ dφ/dλ = a[(r² + a²)E - aL_z]/Δ - [aE - L_z/sin²θ]
    ///
    /// where:
    /// - Σ = r² + a² cos²θ
    /// - Δ = r² - 2Mr + a²
    /// - R = [(r² + a²)E - aL_z]² - Δ[(L_z - aE)² + Q]
    /// - Θ = Q - cos²θ[a²(1 - E²) + L²/sin²θ]
    fn geodesic(&self, _: f64, y: &EquationOfMotionState) -> EquationOfMotionState {
        let _t = y[0];
        let r = y[1];
        let theta = y[2];
        let _phi = y[3];

        // Current velocities (stored but not directly used in first-order form)
        let v_r = y[5];
        let v_theta = y[6];

        let sigma = self.sigma(r, theta);
        let delta = self.delta(r);
        let sin_theta = theta.sin();
        let cos_theta = theta.cos();

        // Compute R and Θ potentials
        let r_pot = self.r_potential(r);
        let theta_pot = self.theta_potential(theta);

        trace!("r = {}, theta = {}, sigma = {}, delta = {}", r, theta, sigma, delta);
        trace!("R = {}, Θ = {}", r_pot, theta_pot);

        // Avoid singularities
        if delta.abs() < 1e-10 || sigma.abs() < 1e-10 {
            trace!("Near singularity: delta = {}, sigma = {}", delta, sigma);
            return EquationOfMotionState::from_column_slice(&[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        }

        // Compute velocities from first-order equations
        let r2_plus_a2 = r * r + self.a * self.a;
        let term1 = (r2_plus_a2 * self.energy - self.a * self.angular_momentum) / delta;
        let term2 = self.a * self.energy * sin_theta * sin_theta - self.angular_momentum;

        // dt/dλ = (1/Σ)[(r² + a²)[(r² + a²)E - aL_z]/Δ - a[aE sin²θ - L_z]]
        let v_t = (r2_plus_a2 * term1 - self.a * term2) / sigma;

        // dr/dλ = (1/Σ)√R (use sign from current velocity)
        let dr_dlambda = if r_pot >= 0.0 {
            r_pot.sqrt() * v_r.signum()
        } else {
            trace!("Warning: R < 0, setting dr/dλ = 0");
            0.0
        };

        // dθ/dλ = (1/Σ)√Θ (use sign from current velocity)
        let dtheta_dlambda = if theta_pot >= 0.0 {
            theta_pot.sqrt() * v_theta.signum()
        } else {
            trace!("Warning: Θ < 0, setting dθ/dλ = 0");
            0.0
        };

        // dφ/dλ = (1/Σ)[a[(r² + a²)E - aL_z]/Δ - [aE - L_z/sin²θ]]
        let v_phi = if sin_theta.abs() > 1e-10 {
            (self.a * term1 - (self.a * self.energy - self.angular_momentum / (sin_theta * sin_theta))) / sigma
        } else {
            0.0
        };

        trace!("Velocities: v_t = {}, v_r = {}, v_theta = {}, v_phi = {}", v_t, dr_dlambda, dtheta_dlambda, v_phi);

        // In first-order form, "accelerations" are zero (conserved quantities ensure this)
        EquationOfMotionState::from_column_slice(&[
            v_t,
            dr_dlambda,
            dtheta_dlambda,
            v_phi,
            0.0,
            0.0,
            0.0,
            0.0,
        ])
    }

    fn create_initial_state(&self, ray: &Ray) -> EquationOfMotionState {
        // For first-order form, we need to store the direction of motion
        // We'll use the sign of the velocity components from the initial momentum
        let r = ray.position[1];
        let theta = ray.position[2];

        let sigma = self.sigma(r, theta);
        let delta = self.delta(r);

        // Initial velocity signs from momentum
        let v_r_sign = ray.momentum.vector[1].signum();
        let v_theta_sign = ray.momentum.vector[2].signum();

        // Compute initial radial and theta potentials
        let r_pot = self.r_potential(r);
        let theta_pot = self.theta_potential(theta);

        let v_r = if r_pot >= 0.0 {
            r_pot.sqrt() * v_r_sign / sigma
        } else {
            0.0
        };

        let v_theta = if theta_pot >= 0.0 {
            theta_pot.sqrt() * v_theta_sign / sigma
        } else {
            0.0
        };

        EquationOfMotionState::from_column_slice(&[
            ray.position[0],
            ray.position[1],
            ray.position[2],
            ray.position[3],
            0.0, // These will be computed by geodesic()
            v_r,
            v_theta,
            0.0,
        ])
    }
}

impl HasCoordinateSystem for KerrAnalytic {
    fn coordinate_system(&self) -> CoordinateSystem {
        Spherical
    }
}

impl InnerProduct for KerrAnalytic {
    fn inner_product(&self, position: &Point, v: &FourVector, w: &FourVector) -> f64 {
        debug_assert_eq!(position.coordinate_system, Spherical);
        let r = position[1];
        let theta = position[2];

        let sigma = self.sigma(r, theta);
        let delta = self.delta(r);
        let sin_theta = theta.sin();
        let sin2_theta = sin_theta * sin_theta;
        let a2 = self.a * self.a;
        let r2 = r * r;

        // Metric components in Boyer-Lindquist coordinates
        let g_tt = -(1.0 - self.radius * r / sigma);
        let g_tphi = -self.radius * r * self.a * sin2_theta / sigma;
        let g_rr = sigma / delta;
        let g_thetatheta = sigma;
        let g_phiphi = (r2 + a2 + self.radius * r * a2 * sin2_theta / sigma) * sin2_theta;

        g_tt * v.vector[0] * w.vector[0]
            + 2.0 * g_tphi * v.vector[0] * w.vector[3]
            + g_rr * v.vector[1] * w.vector[1]
            + g_thetatheta * v.vector[2] * w.vector[2]
            + g_phiphi * v.vector[3] * w.vector[3]
    }
}

impl Signature for KerrAnalytic {
    fn signature(&self) -> [f64; 4] {
        [1.0, -1.0, -1.0, -1.0]
    }
}

impl Geometry for KerrAnalytic {
    /// Computes orthonormal tetrad for an observer at position.
    /// This is a stationary observer (ZAMO - Zero Angular Momentum Observer).
    fn get_tetrad_at(&self, position: &Point) -> Tetrad {
        assert_eq!(position.coordinate_system, Spherical);
        let r = position[1];
        let theta = position[2];

        let sigma = self.sigma(r, theta);
        let delta = self.delta(r);
        let sin_theta = theta.sin();
        let a2 = self.a * self.a;
        let r2 = r * r;

        // Metric components
        let g_tt = -(1.0 - self.radius * r / sigma);
        let g_tphi = -self.radius * r * self.a * sin_theta * sin_theta / sigma;
        let g_phiphi = (r2 + a2 + self.radius * r * a2 * sin_theta * sin_theta / sigma) * sin_theta * sin_theta;

        // ZAMO 4-velocity (zero angular momentum)
        let omega = -g_tphi / g_phiphi;
        let ut = 1.0 / (-g_tt - 2.0 * g_tphi * omega - g_phiphi * omega * omega).sqrt();
        let uphi = omega * ut;

        // Construct orthonormal basis
        let e_t = FourVector::new_spherical(ut, 0.0, 0.0, uphi);

        // Spatial basis vectors (normalized)
        let e_r = FourVector::new_spherical(0.0, (delta / sigma).sqrt(), 0.0, 0.0);
        let e_theta = FourVector::new_spherical(0.0, 0.0, 1.0 / sigma.sqrt(), 0.0);
        let e_phi = FourVector::new_spherical(0.0, 0.0, 0.0, 1.0 / (g_phiphi.sqrt()));

        Tetrad::new(*position, e_t, e_phi, e_theta, e_r)
    }

    fn lorentz_transformation(&self, position: &Point, velocity: &FourVector) -> Matrix4<f64> {
        // For now, use a simplified version similar to Schwarzschild
        assert_eq!(position.coordinate_system, Spherical);
        let mut matrix = Matrix4::zeros();
        let tetrad_t = self.get_tetrad_at(position).t;

        let r = position[1];
        let theta = position[2];
        let sigma = self.sigma(r, theta);
        let delta = self.delta(r);
        let sin_theta = theta.sin();
        let sin2_theta = sin_theta * sin_theta;
        let a2 = self.a * self.a;
        let r2 = r * r;

        // Metric diagonal-ish components
        let g_tt = -(1.0 - self.radius * r / sigma);
        let g_tphi = -self.radius * r * self.a * sin2_theta / sigma;
        let g_rr = sigma / delta;
        let g_thetatheta = sigma;
        let g_phiphi = (r2 + a2 + self.radius * r * a2 * sin2_theta / sigma) * sin2_theta;

        debug!(
            "scalar prod tetrad_t: {:?}",
            self.inner_product(position, &tetrad_t, &tetrad_t)
        );
        debug!(
            "scalar prod velocity: {:?}",
            self.inner_product(position, velocity, velocity)
        );

        // Compute gamma (Lorentz factor)
        let mut gamma = g_tt * velocity.vector[0] * tetrad_t.vector[0];
        gamma += g_tphi * (velocity.vector[0] * tetrad_t.vector[3] + velocity.vector[3] * tetrad_t.vector[0]);
        gamma += g_rr * velocity.vector[1] * tetrad_t.vector[1];
        gamma += g_thetatheta * velocity.vector[2] * tetrad_t.vector[2];
        gamma += g_phiphi * velocity.vector[3] * tetrad_t.vector[3];

        // Build transformation matrix
        for mu in 0..4 {
            for nu in 0..4 {
                let mut res = 0.0;
                if mu == nu {
                    res = 1.0;
                }

                let a = 1.0 / (1.0 + gamma);
                let b = tetrad_t.vector[mu] + velocity.vector[mu];

                // Compute metric-lowered component
                let mut c = 0.0;
                c += g_tt * (tetrad_t.vector[0] + velocity.vector[0]);
                if nu == 3 || nu == 0 {
                    c += g_tphi * (tetrad_t.vector[3] + velocity.vector[3]);
                }
                if nu == 1 {
                    c = g_rr * (tetrad_t.vector[1] + velocity.vector[1]);
                } else if nu == 2 {
                    c = g_thetatheta * (tetrad_t.vector[2] + velocity.vector[2]);
                } else if nu == 3 {
                    c += g_phiphi * (tetrad_t.vector[3] + velocity.vector[3]);
                }

                res -= a * b * c;

                // Add term with metric
                let mut metric_tetrad_nu = 0.0;
                metric_tetrad_nu += g_tt * tetrad_t.vector[0];
                if nu == 3 {
                    metric_tetrad_nu += g_tphi * tetrad_t.vector[3];
                }
                if nu == 1 {
                    metric_tetrad_nu = g_rr * tetrad_t.vector[1];
                } else if nu == 2 {
                    metric_tetrad_nu = g_thetatheta * tetrad_t.vector[2];
                } else if nu == 3 {
                    metric_tetrad_nu += g_phiphi * tetrad_t.vector[3];
                }

                res += 2.0 * metric_tetrad_nu * velocity.vector[mu];

                matrix[(mu, nu)] = res;
            }
        }
        matrix
    }

    fn get_stationary_velocity_at(&self, position: &Point) -> FourVector {
        // ZAMO velocity
        let tetrad = self.get_tetrad_at(position);
        tetrad.t
    }

    fn inside_horizon(&self, position: &Point) -> bool {
        if self.a > self.radius {
            return false;
        }
        let r = position[1];
        // Outer horizon radius
        let r_plus = 0.5 * self.radius
            + ((0.5 * self.radius) * (0.5 * self.radius) - self.a * self.a).sqrt();
        r <= r_plus + self.horizon_epsilon
    }

    fn closed_orbit(&self, _position: &Point, _step_index: usize, _max_steps: usize) -> bool {
        false
    }

    fn get_geodesic_solver(&self, ray: &Ray) -> Box<dyn GeodesicSolver> {
        // Compute conserved quantities from the initial ray
        let r = ray.position[1];
        let theta = ray.position[2];

        let sigma = self.sigma(r, theta);
        let delta = self.delta(r);
        let sin_theta = theta.sin();
        let sin2_theta = sin_theta * sin_theta;
        let a2 = self.a * self.a;
        let r2 = r * r;

        // Metric components
        let g_tt = -(1.0 - self.radius * r / sigma);
        let g_tphi = -self.radius * r * self.a * sin2_theta / sigma;
        let g_rr = sigma / delta;
        let g_thetatheta = sigma;
        let g_phiphi = (r2 + a2 + self.radius * r * a2 * sin2_theta / sigma) * sin2_theta;

        // Lower the momentum to get covariant components
        let p_t = g_tt * ray.momentum.vector[0] + g_tphi * ray.momentum.vector[3];
        let p_r = g_rr * ray.momentum.vector[1];
        let p_theta = g_thetatheta * ray.momentum.vector[2];
        let p_phi = g_tphi * ray.momentum.vector[0] + g_phiphi * ray.momentum.vector[3];

        // Conserved quantities
        let energy = -p_t;
        let angular_momentum = p_phi;

        // Carter constant from theta component and other conserved quantities
        let carter_constant = p_theta * p_theta
            + theta.cos() * theta.cos()
                * (a2 * energy * energy - angular_momentum * angular_momentum / sin2_theta);

        debug!("Conserved quantities: E = {}, L_z = {}, Q = {}", energy, angular_momentum, carter_constant);

        Box::new(KerrAnalyticSolver {
            radius: self.radius,
            a: self.a,
            energy,
            angular_momentum,
            carter_constant,
        })
    }
}

#[cfg(test)]
mod tests {
    use crate::geometry::four_vector::FourVector;
    use crate::geometry::geometry::{Geometry, InnerProduct};
    use crate::geometry::kerr_analytic::KerrAnalytic;
    use crate::geometry::point::Point;
    use crate::geometry::spherical_coordinates_helper::cartesian_to_spherical;
    use crate::rendering::camera::Camera;
    use approx::assert_abs_diff_eq;
    use std::f64::consts::PI;

    const NO_ANGULAR_MOMENTUM: f64 = 0.0;
    const EPSILON: f64 = 1e-4;

    #[test]
    fn test_tetrad_orthonormal() {
        let position = cartesian_to_spherical(&Point::new_cartesian(2.0, 3.0, 4.0, 5.0));
        let geometry = KerrAnalytic::new(2.0, NO_ANGULAR_MOMENTUM, EPSILON);

        let tetrad = geometry.get_tetrad_at(&position);

        let k = tetrad.t + (-tetrad.z);
        let s = geometry.inner_product(&position, &k, &k);
        assert_abs_diff_eq!(s, 0.0, epsilon = 1e-8);

        assert_abs_diff_eq!(geometry.inner_product(&position, &tetrad.t, &tetrad.t), 1.0);
        assert_abs_diff_eq!(
            geometry.inner_product(&position, &tetrad.x, &tetrad.x),
            -1.0
        );
        assert_abs_diff_eq!(
            geometry.inner_product(&position, &tetrad.y, &tetrad.y),
            -1.0
        );
        assert_abs_diff_eq!(
            geometry.inner_product(&position, &tetrad.z, &tetrad.z),
            -1.0
        );

        assert_abs_diff_eq!(geometry.inner_product(&position, &tetrad.t, &tetrad.x), 0.0, epsilon = 1e-8);
        assert_abs_diff_eq!(geometry.inner_product(&position, &tetrad.t, &tetrad.y), 0.0, epsilon = 1e-8);
        assert_abs_diff_eq!(geometry.inner_product(&position, &tetrad.t, &tetrad.z), 0.0, epsilon = 1e-8);

        assert_abs_diff_eq!(geometry.inner_product(&position, &tetrad.x, &tetrad.y), 0.0, epsilon = 1e-8);
        assert_abs_diff_eq!(geometry.inner_product(&position, &tetrad.x, &tetrad.z), 0.0, epsilon = 1e-8);

        assert_abs_diff_eq!(geometry.inner_product(&position, &tetrad.y, &tetrad.z), 0.0, epsilon = 1e-8);
    }

    #[test]
    fn test_kerr_analytic_ray() {
        let position = cartesian_to_spherical(&Point::new_cartesian(2.0, 3.0, 4.0, 5.0));
        let radius = 2.0;
        let geometry = KerrAnalytic::new(radius, NO_ANGULAR_MOMENTUM, EPSILON);

        let r = position[1];
        let a = 1.0 - radius / r;
        let velocity = FourVector::new_spherical(1.0 / a, -(radius / r).sqrt(), 0.0, 0.0);

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

        let ray = camera.get_ray_for(1, 6);
        assert_abs_diff_eq!(
            geometry.inner_product(&position, &ray.momentum, &ray.momentum),
            0.0,
            epsilon = 1e-8
        );
    }
}

//! Boyer-Lindquist Kerr geometry using separated geodesic equations (complete integrability).
//!
//! Unlike the Kerr-Schild Cartesian implementation in `geometry/kerr.rs`, this module works
//! in Boyer-Lindquist coordinates (t, r, θ, φ) and exploits the Carter constant to decouple
//! the r and θ equations of motion.

use nalgebra::{Const, Matrix4, OVector};

use crate::geometry::geometry::{
    ConstantsOfMotion, GeodesicSolver, Geometry, HasCoordinateSystem, InnerProduct, Signature,
    SupportQuantities,
};
use crate::geometry::four_vector::FourVector;
use crate::geometry::point::{CoordinateSystem, Point};
use crate::geometry::tetrad::Tetrad;
use crate::rendering::ray::Ray;
use crate::rendering::raytracer::RaytracerError;
use crate::rendering::runge_kutta::OdeFunction;
use crate::rendering::scene::EquationOfMotionState;
use crate::rendering::temperature::{KerrTemperatureComputer, TemperatureComputer};

fn sigma(r: f64, a: f64, theta: f64) -> f64 {
    r * r + a * a * theta.cos().powi(2)
}

fn delta(r: f64, r_s: f64, a: f64) -> f64 {
    r * r - r_s * r + a * a
}

/// R(r) = [(r² + a²)E - aL_z]² - Δ[(L_z - aE)² + Q]
fn potential_r(r: f64, r_s: f64, a: f64, e: f64, l_z: f64, q: f64) -> f64 {
    let del = delta(r, r_s, a);
    let p_r = (r * r + a * a) * e - a * l_z;
    p_r * p_r - del * ((l_z - a * e).powi(2) + q)
}

/// R'(r) = 4rE[(r² + a²)E - aL_z] - (2r - r_s)[(L_z - aE)² + Q]
fn potential_r_derivative(r: f64, r_s: f64, a: f64, e: f64, l_z: f64, q: f64) -> f64 {
    let p_r = (r * r + a * a) * e - a * l_z;
    let carter_term = (l_z - a * e).powi(2) + q;
    4.0 * r * e * p_r - (2.0 * r - r_s) * carter_term
}

/// Θ(θ) = Q + a²E²cos²θ - L_z²cos²θ/sin²θ
fn potential_theta(theta: f64, a: f64, e: f64, l_z: f64, q: f64) -> f64 {
    let cos_t = theta.cos();
    let sin_t = theta.sin();
    q + a * a * e * e * cos_t * cos_t - l_z * l_z * cos_t * cos_t / (sin_t * sin_t)
}

/// Θ'(θ) = -2a²E²cosθsinθ + 2L_z²cosθ/sin³θ
fn potential_theta_derivative(theta: f64, a: f64, e: f64, l_z: f64, q: f64) -> f64 {
    let cos_t = theta.cos();
    let sin_t = theta.sin();
    -2.0 * a * a * e * e * cos_t * sin_t + 2.0 * l_z * l_z * cos_t / (sin_t.powi(3))
}

struct KerrBLSolver {
    radius: f64,
    a: f64,
    e: f64,
    l_z: f64,
    q: f64,
}

impl HasCoordinateSystem for KerrBLSolver {
    fn coordinate_system(&self) -> CoordinateSystem {
        CoordinateSystem::BoyerLindquist { a: self.a }
    }
}

impl OdeFunction<Const<8>> for KerrBLSolver {
    fn apply(&self, t: f64, y: &OVector<f64, Const<8>>) -> OVector<f64, Const<8>> {
        self.geodesic(t, y)
    }
}

impl GeodesicSolver for KerrBLSolver {
    fn geodesic(&self, _t: f64, y: &EquationOfMotionState) -> EquationOfMotionState {
        let r = y[1];
        let theta = y[2];
        let v_r = y[4];
        let v_theta = y[5];

        let del = delta(r, self.radius, self.a);
        let p_r = (r * r + self.a * self.a) * self.e - self.a * self.l_z;

        // dt/dλ = (r²+a²)/Δ * P_r + a(L_z - aE sin²θ)
        let sin_t = theta.sin();
        let sin2 = sin_t * sin_t;
        let dt = (r * r + self.a * self.a) / del * p_r + self.a * (self.l_z - self.a * self.e * sin2);

        // dφ/dλ = a/Δ * P_r + L_z/sin²θ - aE
        let dphi = self.a / del * p_r + self.l_z / sin2 - self.a * self.e;

        // d²r/dλ² = R'(r)/2
        let dv_r = potential_r_derivative(r, self.radius, self.a, self.e, self.l_z, self.q) / 2.0;

        // d²θ/dλ² = Θ'(θ)/2
        let dv_theta = potential_theta_derivative(theta, self.a, self.e, self.l_z, self.q) / 2.0;

        EquationOfMotionState::from_column_slice(&[dt, v_r, v_theta, dphi, dv_r, dv_theta, 0.0, 0.0])
    }

    fn create_initial_state(&self, _ray: &Ray) -> EquationOfMotionState {
        // Will be implemented in Task 5
        todo!()
    }

    fn momentum_from_state(&self, _y: &EquationOfMotionState) -> FourVector {
        // Will be implemented in Task 5
        todo!()
    }
}

/// Covariant BL Kerr metric g_μν. Signature (−,+,+,+).
fn metric_bl(r_s: f64, a: f64, r: f64, theta: f64) -> Matrix4<f64> {
    let sig = sigma(r, a, theta);
    let sin_t = theta.sin();
    let sin2 = sin_t * sin_t;

    let g_tt = -(1.0 - r_s * r / sig);
    let g_rr = sig / delta(r, r_s, a);
    let g_thth = sig;
    let g_phph = (r * r + a * a + a * a * r_s * r * sin2 / sig) * sin2;
    let g_tph = -a * r_s * r * sin2 / sig;

    let mut g = Matrix4::zeros();
    g[(0, 0)] = g_tt;
    g[(1, 1)] = g_rr;
    g[(2, 2)] = g_thth;
    g[(3, 3)] = g_phph;
    g[(0, 3)] = g_tph;
    g[(3, 0)] = g_tph;
    g
}

/// Contravariant BL Kerr metric g^μν.
fn metric_bl_contravariant(r_s: f64, a: f64, r: f64, theta: f64) -> Matrix4<f64> {
    let sig = sigma(r, a, theta);
    let del = delta(r, r_s, a);
    let sin_t = theta.sin();
    let sin2 = sin_t * sin_t;
    let r2 = r * r;
    let a2 = a * a;
    let big_a = (r2 + a2).powi(2) - del * a2 * sin2;

    let mut g = Matrix4::zeros();
    g[(0, 0)] = -big_a / (sig * del);
    g[(1, 1)] = del / sig;
    g[(2, 2)] = 1.0 / sig;
    g[(3, 3)] = (del - a2 * sin2) / (sig * del * sin2);
    g[(0, 3)] = -a * r_s * r / (sig * del);
    g[(3, 0)] = g[(0, 3)];
    g
}

#[derive(Clone, Debug)]
pub struct KerrBL {
    pub radius: f64,          // r_s = 2M (Schwarzschild radius)
    pub a: f64,               // spin parameter
    pub horizon_epsilon: f64,
}

impl KerrBL {
    pub fn new(radius: f64, a: f64, horizon_epsilon: f64) -> Self {
        KerrBL {
            radius,
            a,
            horizon_epsilon,
        }
    }
}

impl HasCoordinateSystem for KerrBL {
    fn coordinate_system(&self) -> CoordinateSystem {
        CoordinateSystem::BoyerLindquist { a: self.a }
    }
}

impl Signature for KerrBL {
    fn signature(&self) -> [f64; 4] {
        [-1.0, 1.0, 1.0, 1.0]
    }
}

impl InnerProduct for KerrBL {
    fn inner_product(&self, position: &Point, v: &FourVector, w: &FourVector) -> f64 {
        debug_assert!(
            matches!(position.coordinate_system, CoordinateSystem::BoyerLindquist { .. }),
            "inner_product expects BL coordinates"
        );
        let r = position[1];
        let theta = position[2];
        let g = metric_bl(self.radius, self.a, r, theta);

        let mut result = 0.0;
        for mu in 0..4 {
            for nu in 0..4 {
                result += g[(mu, nu)] * v.vector[mu] * w.vector[nu];
            }
        }
        result
    }
}

impl SupportQuantities for KerrBL {
    fn get_stationary_velocity_at(&self, _position: &Point) -> FourVector {
        todo!()
    }

    fn get_circular_orbit_velocity_at(
        &self,
        _position: &Point,
    ) -> Result<FourVector, RaytracerError> {
        todo!()
    }

    fn get_temperature_computer(
        &self,
        temperature: f64,
        _inner_radius: f64,
        outer_radius: f64,
    ) -> Result<Box<dyn TemperatureComputer>, RaytracerError> {
        Ok(Box::new(KerrTemperatureComputer::new(
            temperature,
            outer_radius,
            self.a,
            self.radius,
        )?))
    }
}

impl Geometry for KerrBL {
    fn get_tetrad_at(&self, _position: &Point) -> Tetrad {
        todo!()
    }

    fn lorentz_transformation(&self, _position: &Point, _velocity: &FourVector) -> Matrix4<f64> {
        todo!()
    }

    fn inside_horizon(&self, position: &Point) -> bool {
        // r_plus = M + sqrt(M^2 - a^2) where M = radius/2 (since radius = r_s = 2M)
        let m = self.radius / 2.0;
        if self.a > m {
            return false;
        }
        let r_plus = m + (m * m - self.a * self.a).sqrt();
        position[1] <= r_plus + self.horizon_epsilon
    }

    fn closed_orbit(&self, _position: &Point, _step_index: usize, _max_steps: usize) -> bool {
        false
    }

    fn get_geodesic_solver(&self, _ray: &Ray) -> Box<dyn GeodesicSolver> {
        todo!()
    }

    fn get_radial_coordinate(&self, position: &Point) -> f64 {
        position[1]
    }

    fn get_constants_of_motion(
        &self,
        _position: &Point,
        _momentum: &FourVector,
    ) -> ConstantsOfMotion {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_bl_metric_inverse() {
        let r_s = 1.0;
        let a = 0.5;
        let test_points = [(5.0_f64, 1.2_f64), (3.0, 0.8), (10.0, 2.5)];
        for (r, theta) in test_points {
            let g = metric_bl(r_s, a, r, theta);
            let g_inv = g.try_inverse().expect("Metric should be invertible");
            let g_contra = metric_bl_contravariant(r_s, a, r, theta);
            for i in 0..4 {
                for j in 0..4 {
                    assert_abs_diff_eq!(g_contra[(i, j)], g_inv[(i, j)], epsilon = 1e-10);
                }
            }
        }
    }

    #[test]
    fn test_bl_metric_schwarzschild_limit() {
        // When a=0, BL Kerr reduces to Schwarzschild metric
        let r_s = 2.0;
        let a = 0.0;
        let r = 5.0;
        let theta = 1.2_f64;
        let g = metric_bl(r_s, a, r, theta);
        let a_factor = 1.0 - r_s / r;
        assert_abs_diff_eq!(g[(0, 0)], -a_factor, epsilon = 1e-12);
        assert_abs_diff_eq!(g[(1, 1)], 1.0 / a_factor, epsilon = 1e-12);
        assert_abs_diff_eq!(g[(2, 2)], r * r, epsilon = 1e-12);
        assert_abs_diff_eq!(g[(3, 3)], r * r * theta.sin().powi(2), epsilon = 1e-12);
        assert_abs_diff_eq!(g[(0, 3)], 0.0, epsilon = 1e-12);
    }

    #[test]
    fn test_metric_is_symmetric() {
        let g = metric_bl(1.0, 0.4, 4.0, 1.1);
        for i in 0..4 {
            for j in 0..4 {
                assert_abs_diff_eq!(g[(i, j)], g[(j, i)], epsilon = 1e-15);
            }
        }
    }

    #[test]
    fn test_inner_product_null_vector() {
        let kerr = KerrBL::new(1.0, 0.5, 1e-4);
        let r = 5.0;
        let theta = std::f64::consts::FRAC_PI_2;
        let position = Point::new(0.0, r, theta, 0.0, CoordinateSystem::BoyerLindquist { a: 0.5 });

        // Construct a null vector from metric: g_μν k^μ k^ν = 0
        // For a radial null ray: k^t and k^r with g_tt (k^t)^2 + g_rr (k^r)^2 = 0
        let g = metric_bl(1.0, 0.5, r, theta);
        let kt = 1.0;
        let kr = (-g[(0, 0)] / g[(1, 1)]).sqrt() * kt;
        let k = FourVector::new_boyer_lindquist(0.5, kt, kr, 0.0, 0.0);
        let ip = kerr.inner_product(&position, &k, &k);
        assert_abs_diff_eq!(ip, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_inner_product_with_cross_term() {
        // Test that the g_tφ cross-term is correctly included.
        // For a vector with both k^t and k^φ nonzero, the inner product
        // must match the manual computation using the full metric.
        let r_s = 1.0;
        let a = 0.5;
        let r = 5.0;
        let theta = std::f64::consts::FRAC_PI_2;
        let kerr = KerrBL::new(r_s, a, 1e-4);
        let position = Point::new(0.0, r, theta, 0.0, CoordinateSystem::BoyerLindquist { a });

        let g = metric_bl(r_s, a, r, theta);
        // Arbitrary vector with all components nonzero
        let kt = 1.0_f64;
        let kr = 0.3_f64;
        let kth = 0.1_f64;
        let kph = 0.5_f64;
        let k = FourVector::new_boyer_lindquist(a, kt, kr, kth, kph);

        // Manual computation using the metric
        let components = [kt, kr, kth, kph];
        let mut expected = 0.0_f64;
        for mu in 0..4 {
            for nu in 0..4 {
                expected += g[(mu, nu)] * components[mu] * components[nu];
            }
        }

        let ip = kerr.inner_product(&position, &k, &k);
        assert_abs_diff_eq!(ip, expected, epsilon = 1e-12);
    }

    #[test]
    fn test_inside_horizon() {
        let kerr = KerrBL::new(1.0, 0.3, 1e-4);
        let m = 0.5_f64; // M = r_s/2
        let r_plus = m + (m * m - 0.3_f64 * 0.3_f64).sqrt();

        let inside = Point::new(0.0, r_plus - 0.01, 1.0, 0.0, CoordinateSystem::BoyerLindquist { a: 0.3 });
        let outside = Point::new(0.0, r_plus + 0.1, 1.0, 0.0, CoordinateSystem::BoyerLindquist { a: 0.3 });
        assert!(kerr.inside_horizon(&inside));
        assert!(!kerr.inside_horizon(&outside));
    }

    #[test]
    fn test_get_radial_coordinate() {
        let kerr = KerrBL::new(1.0, 0.5, 1e-4);
        let position = Point::new(0.0, 7.5, 1.2, 0.8, CoordinateSystem::BoyerLindquist { a: 0.5 });
        assert_abs_diff_eq!(kerr.get_radial_coordinate(&position), 7.5);
    }

    #[test]
    fn test_potential_r_non_negative_allowed_region() {
        // R(r) must be non-negative for physically allowed radial motion
        let r_s = 1.0;
        let a = 0.5;
        let r = 5.0;
        let _theta = std::f64::consts::FRAC_PI_2;
        let e = 1.0;
        let l_z = 3.0;
        // Q for equatorial orbit (theta=PI/2, p_theta=0): Q = 0
        let q = 0.0;

        let r_potential = potential_r(r, r_s, a, e, l_z, q);
        assert!(r_potential >= 0.0, "R(r) must be non-negative for allowed motion");
    }

    #[test]
    fn test_potential_derivatives_numerical() {
        // Compare analytical derivatives against numerical finite differences
        let r_s = 1.0;
        let a = 0.5;
        let e = 1.0;
        let l_z = 3.0;
        let q = 1.0;
        let r = 5.0;
        let h = 1e-7;

        let dr_analytical = potential_r_derivative(r, r_s, a, e, l_z, q);
        let dr_numerical = (potential_r(r + h, r_s, a, e, l_z, q)
            - potential_r(r - h, r_s, a, e, l_z, q))
            / (2.0 * h);
        assert_abs_diff_eq!(dr_analytical, dr_numerical, epsilon = 1e-4);

        let theta = 1.2;
        let dth_analytical = potential_theta_derivative(theta, a, e, l_z, q);
        let dth_numerical = (potential_theta(theta + h, a, e, l_z, q)
            - potential_theta(theta - h, a, e, l_z, q))
            / (2.0 * h);
        assert_abs_diff_eq!(dth_analytical, dth_numerical, epsilon = 1e-4);
    }

    #[test]
    fn test_geodesic_rhs_structure() {
        // Verify the ODE RHS structure:
        // ẏ[0] = dt/dλ, ẏ[1] = y[4], ẏ[2] = y[5], ẏ[3] = dφ/dλ,
        // ẏ[4] = R'(r)/2, ẏ[5] = Θ'(θ)/2, ẏ[6] = 0, ẏ[7] = 0
        let r_s = 1.0;
        let a = 0.5;
        let e = 1.0;
        let l_z = 3.0;
        let q = 1.0;
        let solver = KerrBLSolver { radius: r_s, a, e, l_z, q };

        let r = 5.0;
        let theta = 1.2;
        let v_r = 0.1;
        let v_theta = -0.05;
        let y = EquationOfMotionState::from_column_slice(&[0.0, r, theta, 0.0, v_r, v_theta, 0.0, 0.0]);

        let rhs = solver.geodesic(0.0, &y);

        // ẏ[1] = v_r = y[4]
        assert_abs_diff_eq!(rhs[1], v_r, epsilon = 1e-12);
        // ẏ[2] = v_theta = y[5]
        assert_abs_diff_eq!(rhs[2], v_theta, epsilon = 1e-12);
        // ẏ[4] = R'(r)/2
        assert_abs_diff_eq!(rhs[4], potential_r_derivative(r, r_s, a, e, l_z, q) / 2.0, epsilon = 1e-12);
        // ẏ[5] = Θ'(θ)/2
        assert_abs_diff_eq!(rhs[5], potential_theta_derivative(theta, a, e, l_z, q) / 2.0, epsilon = 1e-12);
        // ẏ[6] = 0, ẏ[7] = 0
        assert_abs_diff_eq!(rhs[6], 0.0, epsilon = 1e-12);
        assert_abs_diff_eq!(rhs[7], 0.0, epsilon = 1e-12);
    }
}

//! Boyer-Lindquist Kerr geometry using separated geodesic equations (complete integrability).
//!
//! Unlike the Kerr-Schild Cartesian implementation in `geometry/kerr.rs`, this module works
//! in Boyer-Lindquist coordinates (t, r, θ, φ) and exploits the Carter constant to decouple
//! the r and θ equations of motion.

use nalgebra::Matrix4;

use crate::geometry::geometry::{
    ConstantsOfMotion, GeodesicSolver, Geometry, HasCoordinateSystem, InnerProduct, Signature,
    SupportQuantities,
};
use crate::geometry::four_vector::FourVector;
use crate::geometry::point::{CoordinateSystem, Point};
use crate::geometry::tetrad::Tetrad;
use crate::rendering::ray::Ray;
use crate::rendering::raytracer::RaytracerError;
use crate::rendering::temperature::{KerrTemperatureComputer, TemperatureComputer};

fn sigma(r: f64, a: f64, theta: f64) -> f64 {
    r * r + a * a * theta.cos().powi(2)
}

fn delta(r: f64, r_s: f64, a: f64) -> f64 {
    r * r - r_s * r + a * a
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
    use crate::geometry::point::CoordinateSystem;
    use crate::geometry::geometry::InnerProduct;

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
}

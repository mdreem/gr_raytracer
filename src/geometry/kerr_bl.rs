//! Boyer-Lindquist Kerr geometry using separated geodesic equations (complete integrability).
//!
//! Unlike the Kerr-Schild Cartesian implementation in `geometry/kerr.rs`, this module works
//! in Boyer-Lindquist coordinates (t, r, θ, φ) and exploits the Carter constant to decouple
//! the r and θ equations of motion.

use nalgebra::Matrix4;

// Metric helper functions are consumed by the Geometry trait impls added in later tasks.
// The #[allow(dead_code)] is temporary until those impls are complete.
#[allow(dead_code)]
fn sigma(r: f64, a: f64, theta: f64) -> f64 {
    r * r + a * a * theta.cos().powi(2)
}

#[allow(dead_code)]
fn delta(r: f64, r_s: f64, a: f64) -> f64 {
    r * r - r_s * r + a * a
}

/// Covariant BL Kerr metric g_μν. Signature (−,+,+,+).
#[allow(dead_code)]
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
#[allow(dead_code)]
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

#[allow(dead_code)]
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
}

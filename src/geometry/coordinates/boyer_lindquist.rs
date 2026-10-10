use crate::geometry::coordinates::atlas::{Atlas, Chart};
use crate::geometry::four_vector::FourVector;
use crate::geometry::point::Point;
use nalgebra::{Matrix4, Vector4};
use std::f64::consts::{FRAC_PI_2, PI};

const BL_TO_KS_SIN_THETA: f64 = 0.1;
const KS_TO_BL_SIN_THETA: f64 = 0.5;
const KS_TO_BL_MIN_R_FACTOR: f64 = 3.0;

/// Boyer-Lindquist coordinates (t, r, θ, φ) for the Kerr metric.
///
/// The Cartesian embedding used is the null-tetrad convention consistent with
/// the Kerr-Schild form used in `geometry/kerr.rs`:
///   x = (r cosφ − a sinφ) sinθ
///   y = (r sinφ + a cosφ) sinθ
///   z = r cosθ
///
/// Note: this differs from the oblate-spheroidal convention
/// (x = √(r²+a²) sinθ cosφ) used in some references.
///
/// # Note on equality
/// `PartialEq` is derived, which compares `a` and `r_s` with exact `f64` equality.
/// This is correct as long as they are always initialised from the same
/// values (e.g. `kerr_bl.a`) rather than an independent floating-point
/// computation.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct BoyerLindquistData {
    pub a: f64,
    pub r_s: f64,
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct KerrSchildData {
    pub a: f64,
    pub r_s: f64,
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum KerrChart {
    BoyerLindquist(BoyerLindquistData),
    KerrSchild(KerrSchildData),
}

impl BoyerLindquistData {
    /// sin theta of a point given in BL coordinates (t, r, theta, phi).
    pub fn sin_theta(&self, v: &Vector4<f64>) -> f64 {
        v[2].sin()
    }
}

impl KerrSchildData {
    /// BL r of a point given in Kerr-Schild coordinates (t, x, y, z).
    pub fn bl_radius(&self, v: &Vector4<f64>) -> f64 {
        let (x, y, z) = (v[1], v[2], v[3]);
        let rho2 = x * x + y * y + z * z;
        let a2 = self.a * self.a;
        (0.5 * (rho2 - a2 + ((rho2 - a2).powi(2) + 4.0 * a2 * z * z).sqrt())).sqrt()
    }

    /// sin theta_BL of a point given in Kerr-Schild coordinates (t, x, y, z): cos theta = z / r.
    pub fn sin_theta(&self, v: &Vector4<f64>) -> f64 {
        let r = self.bl_radius(v);
        if r == 0.0 {
            return 0.0;
        }
        let cos_theta = (v[3] / r).clamp(-1.0, 1.0);
        (1.0 - cos_theta * cos_theta).sqrt()
    }
}

impl KerrChart {
    pub fn sin_theta(point: &Point<KerrChart>) -> f64 {
        match point.chart {
            KerrChart::BoyerLindquist(d) => d.sin_theta(&point.vector),
            KerrChart::KerrSchild(d) => d.sin_theta(&point.vector),
        }
    }
}

impl Chart for KerrChart {
    fn tetrad_handedness_fallback(&self) -> f64 {
        match self {
            KerrChart::BoyerLindquist(_) => -1.0,
            KerrChart::KerrSchild(_) => 1.0,
        }
    }
}

fn r_plus(a: f64, r_s: f64) -> f64 {
    let m = 0.5 * r_s;
    let disc = (m * m - a * a).max(0.0);
    m + disc.sqrt()
}

fn delta(r: f64, r_s: f64, a: f64) -> f64 {
    r * r - r_s * r + a * a
}

/// Antiderivatives (∫ dr/Δ, ∫ r dr/Δ) with Δ = r² − r_s r + a², integration constants zero.
/// They give the Kerr-Schild twist: φ_KS − φ_BL = a ∫ dr/Δ, t_KS − t_BL = r_s ∫ r dr/Δ.
/// Singular at the horizons r±, so transitions must happen well away from them.
fn twist_integrals(r: f64, a: f64, r_s: f64) -> (f64, f64) {
    let m = 0.5 * r_s;
    let disc = m * m - a * a;
    if disc > 0.0 {
        let s = disc.sqrt();
        let (r_p, r_m) = (m + s, m - s);
        let (l_p, l_m) = ((r - r_p).abs().ln(), (r - r_m).abs().ln());
        (
            (l_p - l_m) / (r_p - r_m),
            (r_p * l_p - r_m * l_m) / (r_p - r_m),
        )
    } else if disc < 0.0 {
        let b = (-disc).sqrt();
        let at = ((r - m) / b).atan();
        // The constant makes the azimuthal twist vanish at infinity, as in the other cases.
        (
            (at - FRAC_PI_2) / b,
            0.5 * delta(r, r_s, a).ln() + m / b * at,
        )
    } else {
        let u = r - m;
        (-1.0 / u, u.abs().ln() - m / u)
    }
}

/// Jacobian ∂(x^Cartesian)/∂(x^BL) as a 4×4 matrix.
/// Row μ, column ν = ∂x^μ_Cart / ∂x^ν_BL
pub(crate) fn jacobian_bl_to_cartesian(
    r_s: f64,
    a: f64,
    r: f64,
    theta: f64,
    phi: f64,
) -> Matrix4<f64> {
    let (st, ct) = (theta.sin(), theta.cos());
    let (sp, cp) = (phi.sin(), phi.cos());
    let del = delta(r, r_s, a);

    // The (ingoing) Kerr-Schild chart's time and azimuth are twisted against
    // Boyer-Lindquist by dt_KS = dt_BL + (r_s r / Δ) dr and
    // dφ_KS = dφ_BL + (a / Δ) dr, so the r column carries the corresponding
    // twist terms; with them J^T g_KS J = g_BL holds to machine precision
    // (see test_bl_jacobian_transforms_metric_exactly).
    let dx_dphi = (-r * sp - a * cp) * st;
    let dy_dphi = (r * cp - a * sp) * st;

    #[rustfmt::skip]
    let data = [
        1.0, r_s * r / del,                    0.0,                    0.0,
        0.0, st * cp + (a / del) * dx_dphi,    (r * cp - a * sp) * ct, dx_dphi,
        0.0, st * sp + (a / del) * dy_dphi,    (r * sp + a * cp) * ct, dy_dphi,
        0.0, ct,                               -r * st,                0.0,
    ];

    Matrix4::from_row_slice(&data)
}

/// BL (t, r, θ, φ_BL) → Kerr-Schild (t_KS, x, y, z), including the t and φ twist.
pub(crate) fn bl_to_ks(a: f64, r_s: f64, v: &Vector4<f64>) -> Vector4<f64> {
    let (t, r, theta, phi) = (v[0], v[1], v[2], v[3]);
    let (i_phi, i_t) = twist_integrals(r, a, r_s);
    let phi_ks = phi + a * i_phi;
    let (st, ct) = theta.sin_cos();
    let (sp, cp) = phi_ks.sin_cos();
    Vector4::new(
        t + r_s * i_t,
        (r * cp - a * sp) * st,
        (r * sp + a * cp) * st,
        r * ct,
    )
}

/// Wrap an azimuth into (−π, π].
fn wrap_phi(phi: f64) -> f64 {
    let wrapped = (phi + PI).rem_euclid(2.0 * PI) - PI;
    if wrapped == -PI { PI } else { wrapped }
}

/// Kerr-Schild (t_KS, x, y, z) → BL (t, r, θ, φ_BL); also returns φ_KS for the Jacobian.
pub(crate) fn ks_to_bl(a: f64, r_s: f64, v: &Vector4<f64>) -> (Vector4<f64>, f64) {
    let (x, y, z) = (v[1], v[2], v[3]);
    let r = KerrSchildData { a, r_s }.bl_radius(v);
    let theta = if r == 0.0 {
        0.0
    } else {
        (z / r).clamp(-1.0, 1.0).acos()
    };
    let phi_ks = (r * y - a * x).atan2(r * x + a * y);
    let (i_phi, i_t) = twist_integrals(r, a, r_s);
    (
        Vector4::new(v[0] - r_s * i_t, r, theta, wrap_phi(phi_ks - a * i_phi)),
        phi_ks,
    )
}

impl Atlas for KerrChart {
    fn best_chart(&self, position: &Point<Self>) -> Self {
        match position.chart {
            KerrChart::BoyerLindquist(d) if d.sin_theta(&position.vector) < BL_TO_KS_SIN_THETA => {
                KerrChart::KerrSchild(KerrSchildData { a: d.a, r_s: d.r_s })
            }
            KerrChart::KerrSchild(d)
                if d.sin_theta(&position.vector) > KS_TO_BL_SIN_THETA
                    && d.bl_radius(&position.vector)
                        > KS_TO_BL_MIN_R_FACTOR * r_plus(d.a, d.r_s) =>
            {
                KerrChart::BoyerLindquist(BoyerLindquistData { a: d.a, r_s: d.r_s })
            }
            current => current,
        }
    }

    fn transition_point(&self, position: &Point<Self>, to: Self) -> Point<Self> {
        match (position.chart, to) {
            (from, to) if from == to => *position,
            (KerrChart::BoyerLindquist(d), KerrChart::KerrSchild(_)) => {
                Point::new_from_vector(bl_to_ks(d.a, d.r_s, &position.vector), to)
            }
            (KerrChart::KerrSchild(d), KerrChart::BoyerLindquist(_)) => {
                Point::new_from_vector(ks_to_bl(d.a, d.r_s, &position.vector).0, to)
            }
            (from, to) => unreachable!(
                "transition between Kerr charts with different parameters: {:?} -> {:?}",
                from, to
            ),
        }
    }

    fn transition_four_vector(&self, vector: &FourVector<Self>, to: Self) -> FourVector<Self> {
        let at = self.transition_point(&vector.at, to);
        let components = match (vector.at.chart, to) {
            (from, to) if from == to => vector.vector,
            (KerrChart::BoyerLindquist(d), KerrChart::KerrSchild(_)) => {
                let (r, theta) = (vector.at[1], vector.at[2]);
                let phi_ks = (at[2] * r - d.a * at[1]).atan2(at[1] * r + d.a * at[2]);
                jacobian_bl_to_cartesian(d.r_s, d.a, r, theta, phi_ks) * vector.vector
            }
            (KerrChart::KerrSchild(d), KerrChart::BoyerLindquist(_)) => {
                let (bl, phi_ks) = ks_to_bl(d.a, d.r_s, &vector.at.vector);
                jacobian_bl_to_cartesian(d.r_s, d.a, bl[1], bl[2], phi_ks)
                    .try_inverse()
                    .expect("BL Jacobian should be invertible away from the horizons and the axis")
                    * vector.vector
            }
            (from, to) => unreachable!(
                "transition between Kerr charts with different parameters: {:?} -> {:?}",
                from, to
            ),
        };
        FourVector::new(
            components[0],
            components[1],
            components[2],
            components[3],
            at,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::kerr::test_helpers::metric_cartesian;
    use crate::geometry::kerr_bl::metric_bl;
    use approx::assert_abs_diff_eq;

    const R_S: f64 = 1.0;
    // Sub-extremal, extremal (a = M = r_s / 2) and over-extremal spins.
    const SPINS: [f64; 3] = [0.3, 0.5, 0.7];
    const POINTS: [(f64, f64, f64); 4] = [
        (5.0, 1.2, 0.7),
        (8.0, 2.1, -1.3),
        (2.5, 0.4, 2.9),
        (40.0, 1.5707963267948966, 0.1),
    ];

    fn bl(a: f64) -> KerrChart {
        KerrChart::BoyerLindquist(BoyerLindquistData { a, r_s: R_S })
    }

    fn ks(a: f64) -> KerrChart {
        KerrChart::KerrSchild(KerrSchildData { a, r_s: R_S })
    }

    #[test]
    fn twist_integrals_are_antiderivatives_of_one_and_r_over_delta() {
        let h = 1e-6;
        for a in SPINS {
            for r in [2.5, 5.0, 8.0, 40.0] {
                let (p0, p1) = twist_integrals(r + h, a, R_S);
                let (m0, m1) = twist_integrals(r - h, a, R_S);
                let del = delta(r, R_S, a);
                assert_abs_diff_eq!((p0 - m0) / (2.0 * h), 1.0 / del, epsilon = 1e-8);
                assert_abs_diff_eq!((p1 - m1) / (2.0 * h), r / del, epsilon = 1e-8);
            }
        }
    }

    #[test]
    fn point_round_trip_bl_ks_bl() {
        for a in SPINS {
            for (r, theta, phi) in POINTS {
                let p = Point::new(3.0, r, theta, phi, bl(a));
                let back = p
                    .chart
                    .transition_point(&p.chart.transition_point(&p, ks(a)), bl(a));
                assert_eq!(back.chart, bl(a));
                assert_abs_diff_eq!(back.vector, p.vector, epsilon = 1e-10);
            }
        }
    }

    #[test]
    fn jacobian_is_the_derivative_of_the_point_map() {
        let h = 1e-6;
        for a in SPINS {
            for (r, theta, phi) in POINTS {
                let v = Vector4::new(3.0, r, theta, phi);
                let ks_point = bl_to_ks(a, R_S, &v);
                let phi_ks =
                    (ks_point[2] * r - a * ks_point[1]).atan2(ks_point[1] * r + a * ks_point[2]);
                let j = jacobian_bl_to_cartesian(R_S, a, r, theta, phi_ks);
                for nu in 0..4 {
                    let mut e = Vector4::zeros();
                    e[nu] = h;
                    let column =
                        (bl_to_ks(a, R_S, &(v + e)) - bl_to_ks(a, R_S, &(v - e))) / (2.0 * h);
                    for mu in 0..4 {
                        assert_abs_diff_eq!(j[(mu, nu)], column[mu], epsilon = 1e-6);
                    }
                }
            }
        }
    }

    #[test]
    fn vector_round_trip_bl_ks_bl() {
        for a in SPINS {
            for (r, theta, phi) in POINTS {
                let p = Point::new(0.0, r, theta, phi, bl(a));
                let v = FourVector::new(1.3, -0.4, 0.02, 0.05, p);
                let there = p.chart.transition_four_vector(&v, ks(a));
                let back = there.at.chart.transition_four_vector(&there, bl(a));
                assert_abs_diff_eq!(back.vector, v.vector, epsilon = 1e-9);
                assert_abs_diff_eq!(back.at.vector, p.vector, epsilon = 1e-10);
            }
        }
    }

    #[test]
    fn vector_transition_preserves_norm_energy_and_angular_momentum() {
        for a in SPINS {
            for (r, theta, phi) in POINTS {
                let p = Point::new(0.0, r, theta, phi, bl(a));
                let v = FourVector::new(1.3, -0.4, 0.02, 0.05, p);
                let w = p.chart.transition_four_vector(&v, ks(a));
                let (x, y, z) = (w.at[1], w.at[2], w.at[3]);

                let g_bl = metric_bl(R_S, a, r, theta);
                let g_ks = metric_cartesian(R_S, a, x, y, z);
                let p_bl = g_bl * v.vector;
                let p_ks = g_ks * w.vector;

                assert_abs_diff_eq!(v.vector.dot(&p_bl), w.vector.dot(&p_ks), epsilon = 1e-9);
                assert_abs_diff_eq!(p_bl[0], p_ks[0], epsilon = 1e-9);
                assert_abs_diff_eq!(p_bl[3], -y * p_ks[1] + x * p_ks[2], epsilon = 1e-9);
            }
        }
    }

    #[test]
    fn best_chart_switches_near_the_axis_and_back_with_hysteresis() {
        let a = 0.3;
        let near_axis = Point::new(0.0, 10.0, 0.05, 0.3, bl(a));
        assert_eq!(near_axis.chart.best_chart(&near_axis), ks(a));
        let off_axis = Point::new(0.0, 10.0, 0.3, 0.3, bl(a));
        assert_eq!(off_axis.chart.best_chart(&off_axis), bl(a));

        let in_band = off_axis.chart.transition_point(&off_axis, ks(a));
        assert_eq!(in_band.chart.best_chart(&in_band), ks(a));
        let far_from_axis = Point::new(0.0, 10.0, 1.2, 0.3, bl(a));
        let far_ks = far_from_axis.chart.transition_point(&far_from_axis, ks(a));
        assert_eq!(far_ks.chart.best_chart(&far_ks), bl(a));
        let near_horizon = Point::new(0.0, 1.5, 1.2, 0.3, bl(a));
        let near_horizon_ks = near_horizon.chart.transition_point(&near_horizon, ks(a));
        assert_eq!(near_horizon_ks.chart.best_chart(&near_horizon_ks), ks(a));
    }
}

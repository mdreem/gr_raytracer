//! Exact capture classification of photons in Kerr (and Schwarzschild).
//!
//! In Mino time the radial motion of a photon decouples from θ:
//! (dr/dσ)² = R(r) with R(r) = (E(r² + a²) − a L_z)² − Δ (Q + (L_z − aE)²) and
//! Δ = r² − r_s r + a². Divided by E² that is the monic quartic
//! r⁴ + (a² − λ² − η) r² + r_s (η + (λ − a)²) r − a² η in λ = L_z/E, η = Q/E².
//! The photon turns around at the simple roots of R and nowhere else, so
//! whether it ends at the horizon or at infinity follows from where those
//! roots lie relative to its start, without integrating anything. That makes
//! it an exact reference for the integrator's shadow (see `Geometry::analytic_capture`).

use crate::geometry::atlases::kerr::{KerrChart, jacobian_bl_to_cartesian};
use crate::geometry::coordinates::kerr_schild::KerrSchildData;
use crate::geometry::four_vector::FourVector;
use crate::geometry::spacetimes::kerr::metric as metric_ks;
use crate::geometry::spacetimes::kerr_bl::metric_bl;
use nalgebra::Matrix4;

/// Eigenvalues with an imaginary part below this (relative to `r_s + |r|`)
/// count as real roots of R.
const REAL_ROOT_TOLERANCE: f64 = 1e-7;

/// Two roots of R closer than this (relative to `r_s + |r|`) are treated as a
/// double root: the photon sits on (or within ~1e-12 of) a spherical photon
/// orbit, where the roots are too ill-conditioned to decide which way it goes.
const DOUBLE_ROOT_TOLERANCE: f64 = 1e-6;

/// Below this |E| (relative to the momentum's size) λ and η are meaningless.
const MIN_ENERGY: f64 = 1e-12;

/// Outer horizon radius r_+ = M + √(M² − a²), or None for |a| > M.
pub fn outer_horizon(r_s: f64, a: f64) -> Option<f64> {
    let m = 0.5 * r_s;
    if a.abs() > m {
        return None;
    }
    // Guard against FP rounding at extremal Kerr (a == M).
    Some(m + (m * m - a * a).max(0.0).sqrt())
}

/// The radial potential R(r) / E² in terms of λ = L_z/E and η = Q/E².
pub fn radial_potential(r: f64, r_s: f64, a: f64, lambda: f64, eta: f64) -> f64 {
    let p = r * r + a * a - a * lambda;
    p * p - (r * r - r_s * r + a * a) * (eta + (lambda - a).powi(2))
}

/// Real roots of R(r) / E², ascending. Each comes with a flag that is set
/// when another root (real or complex) lies within DOUBLE_ROOT_TOLERANCE.
fn radial_roots(r_s: f64, a: f64, lambda: f64, eta: f64) -> Vec<(f64, bool)> {
    let c2 = a * a - lambda * lambda - eta;
    let c1 = r_s * (eta + (lambda - a).powi(2));
    let c0 = -a * a * eta;
    // Companion matrix of r⁴ + c2 r² + c1 r + c0.
    #[rustfmt::skip]
    let companion = Matrix4::new(
        0.0, 0.0, 0.0, -c0,
        1.0, 0.0, 0.0, -c1,
        0.0, 1.0, 0.0, -c2,
        0.0, 0.0, 1.0, 0.0,
    );
    let eigenvalues = companion.complex_eigenvalues();
    let scale = |r: f64| r_s.max(f64::MIN_POSITIVE) + r.abs();
    let polynomial = |r: f64| ((r * r + c2) * r + c1) * r + c0;
    let derivative = |r: f64| (4.0 * r * r + 2.0 * c2) * r + c1;

    let mut roots: Vec<(f64, bool)> = eigenvalues
        .iter()
        .enumerate()
        .filter(|(_, z)| z.im.abs() < REAL_ROOT_TOLERANCE * scale(z.re))
        .map(|(i, z)| {
            // Polish with a few Newton steps; skipped where R' vanishes.
            let mut r = z.re;
            for _ in 0..3 {
                let d = derivative(r);
                if d == 0.0 {
                    break;
                }
                let next = r - polynomial(r) / d;
                if !next.is_finite() || (next - z.re).abs() > DOUBLE_ROOT_TOLERANCE * scale(r) {
                    break;
                }
                r = next;
            }
            let double = eigenvalues
                .iter()
                .enumerate()
                .any(|(j, w)| j != i && (w - z).norm() < DOUBLE_ROOT_TOLERANCE * scale(z.re));
            (r, double)
        })
        .collect();
    roots.sort_by(|p, q| p.0.total_cmp(&q.0));
    roots
}

/// Whether a photon at Boyer-Lindquist radius `r` with λ = L_z/E and η = Q/E²,
/// moving outward (`outgoing`) or inward in r, ends at the outer horizon
/// (`Some(true)`) or at infinity (`Some(false)`). `None` where that is not
/// decided: no horizon (|a| > M), a photon on or within ~1e-12 of a spherical
/// photon orbit, one bound between two turning points, or a start in a
/// forbidden region (R(r) < 0, inconsistent input).
pub fn photon_captured(
    r_s: f64,
    a: f64,
    lambda: f64,
    eta: f64,
    r: f64,
    outgoing: bool,
) -> Option<bool> {
    let r_plus = outer_horizon(r_s, a)?;
    if !(lambda.is_finite() && eta.is_finite() && r.is_finite()) {
        return None;
    }
    if r <= r_plus {
        return Some(true);
    }
    // R is ~r⁴ for large r, so compare against that scale.
    let magnitude = (r * r + lambda * lambda + eta.abs()).powi(2);
    if radial_potential(r, r_s, a, lambda, eta) < -1e-9 * magnitude {
        return None;
    }
    let roots = radial_roots(r_s, a, lambda, eta);

    // Index of the first turning point and the direction after it.
    let (turn, captured_after_turn) = if outgoing {
        match roots.iter().position(|&(root, _)| root > r) {
            None => return Some(false),
            Some(i) => (i, true),
        }
    } else {
        match roots
            .iter()
            .rposition(|&(root, _)| root < r && root > r_plus)
        {
            None => return Some(true),
            Some(i) => (i, false),
        }
    };
    if roots[turn].1 {
        return None;
    }
    // After turning, the photon moves the other way until the next root or the end.
    let next = if captured_after_turn {
        turn.checked_sub(1).filter(|&i| roots[i].0 > r_plus)
    } else {
        Some(turn + 1).filter(|&i| i < roots.len())
    };
    match next {
        None => Some(captured_after_turn),
        // A second turning point: bound between them (or degenerate).
        Some(_) => None,
    }
}

/// λ, η, the BL radius and the sign of dr/dσ for a photon momentum in either
/// Kerr chart. The curve is classified in the direction it is traced; R is
/// even in (E, L_z), so a past-directed momentum (the renderer traces rays
/// backwards from the camera) gives the same λ and η as its reverse.
fn kerr_invariants(
    r_s: f64,
    a: f64,
    momentum: &FourVector<KerrChart>,
) -> Option<(f64, f64, f64, bool)> {
    let at = &momentum.at;
    // (r, θ, E, L_z, p_θ, a quantity with the sign of dr/dσ)
    let (r, theta, e, l_z, p_theta, dr) = match at.chart {
        KerrChart::BoyerLindquist(_) => {
            let (r, theta) = (at[1], at[2]);
            let p = metric_bl(r_s, a, r, theta) * momentum.vector;
            (r, theta, -p[0], p[3], p[2], momentum[1])
        }
        KerrChart::KerrSchild(_) => {
            let (x, y, z) = (at[1], at[2], at[3]);
            let r = KerrSchildData { a, r_s }.bl_radius(&at.vector);
            if r <= 0.0 {
                return None;
            }
            let theta = (z / r).clamp(-1.0, 1.0).acos();
            let p = metric_ks(r_s, a, x, y, z) * momentum.vector;
            // p_θ is the covariant momentum contracted with ∂x^μ_KS/∂θ_BL,
            // which stays regular on the axis (unlike the inverse Jacobian).
            let phi_ks = (r * y - a * x).atan2(r * x + a * y);
            let d_theta = jacobian_bl_to_cartesian(r_s, a, r, theta, phi_ks)
                .column(2)
                .into_owned();
            // dr from differentiating r⁴ − (ρ² − a²) r² − a² z² = 0; the
            // omitted factor 1 / (r √((ρ² − a²)² + 4a²z²)) is positive.
            let dr = r * r * (x * momentum[1] + y * momentum[2] + z * momentum[3])
                + a * a * z * momentum[3];
            // L_z = p_φ with ∂_φ = x ∂_y − y ∂_x.
            (r, theta, -p[0], -y * p[1] + x * p[2], p.dot(&d_theta), dr)
        }
    };
    // Also rejects a NaN E.
    if e.abs().partial_cmp(&(MIN_ENERGY * momentum.vector.amax()))
        != Some(std::cmp::Ordering::Greater)
    {
        return None;
    }
    let (cos_t, sin2) = (theta.cos(), theta.sin().powi(2));
    // Carter constant (photon form); the floor keeps L_z²/sin²θ finite on the
    // axis, where L_z vanishes anyway.
    let q = p_theta * p_theta + cos_t * cos_t * (l_z * l_z / sin2.max(1e-28) - a * a * e * e);
    Some((l_z / e, q / (e * e), r, dr > 0.0))
}

/// `photon_captured` for a photon momentum in either Kerr chart.
pub fn kerr_capture(r_s: f64, a: f64, momentum: &FourVector<KerrChart>) -> Option<bool> {
    outer_horizon(r_s, a)?;
    let (lambda, eta, r, outgoing) = kerr_invariants(r_s, a, momentum)?;
    photon_captured(r_s, a, lambda, eta, r, outgoing)
}

#[cfg(test)]
/// Bardeen's critical impact parameters (λ, η) of the spherical photon
/// orbit at radius `r` (Bardeen 1973; here with M = r_s / 2).
pub(crate) fn critical_parameters(r_s: f64, a: f64, r: f64) -> (f64, f64) {
    let m = 0.5 * r_s;
    let lambda = -(r.powi(3) - 3.0 * m * r * r + a * a * r + a * a * m) / (a * (r - m));
    let eta = r.powi(3) * (4.0 * a * a * m - r * (r - 3.0 * m).powi(2)) / (a * a * (r - m).powi(2));
    (lambda, eta)
}

#[cfg(test)]
/// Prograde and retrograde equatorial photon orbit radii.
pub(crate) fn photon_orbit_range(r_s: f64, a: f64) -> (f64, f64) {
    let m = 0.5 * r_s;
    let orbit = |sign: f64| 2.0 * m * (1.0 + (2.0 / 3.0 * (sign * (a / m).abs()).acos()).cos());
    (orbit(-1.0), orbit(1.0))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::atlases::Atlas;
    use crate::geometry::coordinates::boyer_lindquist::BoyerLindquistData;
    use crate::geometry::point::Point;

    #[test]
    fn bardeen_orbits_are_double_roots_of_the_radial_potential() {
        // Checks the formulas the other tests lean on: R and R' vanish on the orbit.
        let (r_s, a) = (1.0, 0.499);
        let (r_min, r_max) = photon_orbit_range(r_s, a);
        for k in 1..10 {
            let r = r_min + (r_max - r_min) * k as f64 / 10.0;
            let (lambda, eta) = critical_parameters(r_s, a, r);
            let h = 1e-6;
            let d = (radial_potential(r + h, r_s, a, lambda, eta)
                - radial_potential(r - h, r_s, a, lambda, eta))
                / (2.0 * h);
            assert!(
                radial_potential(r, r_s, a, lambda, eta).abs() < 1e-9,
                "R({r})"
            );
            assert!(d.abs() < 1e-6, "R'({r}) = {d}");
        }
    }

    #[test]
    fn rays_just_inside_the_critical_curve_are_captured() {
        // For a = 0.499 r_s (a/M = 0.998) and every spherical photon orbit,
        // move η by ±1e-6 relative to b² = η + λ² off the critical curve. A
        // ray from far away with the smaller η has no turning point and
        // falls in; with the larger one it turns just outside the orbit.
        let (r_s, a) = (1.0, 0.499);
        let (r_min, r_max) = photon_orbit_range(r_s, a);
        for k in 1..40 {
            let r_orbit = r_min + (r_max - r_min) * k as f64 / 40.0;
            let (lambda, eta) = critical_parameters(r_s, a, r_orbit);
            let shift = 1e-6 * (eta + lambda * lambda);
            for r0 in [20.0, 1e3, 1e5] {
                assert_eq!(
                    photon_captured(r_s, a, lambda, eta - shift, r0, false),
                    Some(true),
                    "r_orbit {r_orbit}, r0 {r0}"
                );
                assert_eq!(
                    photon_captured(r_s, a, lambda, eta + shift, r0, false),
                    Some(false),
                    "r_orbit {r_orbit}, r0 {r0}"
                );
                // Outgoing from far away, it escapes either way.
                assert_eq!(
                    photon_captured(r_s, a, lambda, eta - shift, r0, true),
                    Some(false)
                );
            }
            // Started between the horizon and the orbit, the outcome flips:
            // outgoing, the larger η turns inside the orbit and falls back.
            let r0 = 0.5 * (outer_horizon(r_s, a).unwrap() + r_orbit);
            if radial_potential(r0, r_s, a, lambda, eta + shift) > 0.0 {
                assert_eq!(
                    photon_captured(r_s, a, lambda, eta + shift, r0, true),
                    Some(true),
                    "r_orbit {r_orbit}"
                );
                assert_eq!(
                    photon_captured(r_s, a, lambda, eta - shift, r0, true),
                    Some(false),
                    "r_orbit {r_orbit}"
                );
            }
            // On the orbit itself the classifier declines.
            assert_eq!(photon_captured(r_s, a, lambda, eta, 1e3, false), None);
        }
    }

    #[test]
    fn schwarzschild_limit_matches_the_critical_impact_parameter() {
        // a = 0: captured iff b < 3√3 M from far away, whatever the split
        // of b² into λ² and η.
        let r_s = 2.0;
        let b_c = 1.5 * 3.0_f64.sqrt() * r_s;
        for split in [0.0, 0.3, 1.0] {
            for (b, captured) in [(b_c * (1.0 - 1e-6), true), (b_c * (1.0 + 1e-6), false)] {
                let lambda = (split * b * b).sqrt();
                let eta = b * b - lambda * lambda;
                assert_eq!(
                    photon_captured(r_s, 0.0, lambda, eta, 1e4, false),
                    Some(captured),
                    "b {b}, split {split}"
                );
            }
        }
        // Radial infall and a ray sent outward from inside the photon sphere.
        assert_eq!(photon_captured(r_s, 0.0, 0.0, 0.0, 50.0, false), Some(true));
        assert_eq!(photon_captured(r_s, 0.0, 0.0, 0.0, 2.5, true), Some(false));
        assert_eq!(
            photon_captured(r_s, 0.0, 0.0, 1.1 * b_c * b_c, 2.5, true),
            Some(true)
        );
    }

    #[test]
    fn no_horizon_no_answer() {
        assert_eq!(photon_captured(1.0, 0.7, 1.0, 5.0, 100.0, false), None);
        assert_eq!(outer_horizon(1.0, 0.7), None);
    }

    #[test]
    fn invariants_agree_between_the_kerr_charts() {
        // The same momentum in BL and in Kerr-Schild gives the same λ, η and
        // radial direction, and the same as computed from the BL definitions.
        let (r_s, a) = (1.0, 0.499);
        let bl = KerrChart::BoyerLindquist(BoyerLindquistData { a, r_s });
        let ks = KerrChart::KerrSchild(KerrSchildData { a, r_s });
        for (r, theta, phi, p) in [
            (17.0, 1.48, 0.3, [-1.0, 0.9, 0.004, -0.01]),
            (6.0, 2.2, -1.0, [-1.2, -0.7, 0.05, 0.02]),
            (40.0, 0.4, 2.0, [1.0, -0.95, 0.002, 0.003]),
        ] {
            let x = Point::new(0.0, r, theta, phi, bl);
            let g = metric_bl(r_s, a, r, theta);
            // Complete p^θ so the vector is null.
            let mut v = nalgebra::Vector4::from(p);
            let rest = v.dot(&(g * v)) - g[(2, 2)] * v[2] * v[2];
            v[2] = v[2].signum() * (-rest / g[(2, 2)]).sqrt();
            let momentum = FourVector::new(v[0], v[1], v[2], v[3], x);
            let in_ks = bl.transition_four_vector(&momentum, ks);

            let (l1, e1, r1, o1) = kerr_invariants(r_s, a, &momentum).unwrap();
            let (l2, e2, r2, o2) = kerr_invariants(r_s, a, &in_ks).unwrap();
            assert!((l1 - l2).abs() < 1e-8 * (1.0 + l1.abs()), "{l1} vs {l2}");
            assert!((e1 - e2).abs() < 1e-8 * (1.0 + e1.abs()), "{e1} vs {e2}");
            assert!((r1 - r2).abs() < 1e-9 * r1);
            assert_eq!(o1, o2);
            assert_eq!(o1, v[1] > 0.0);
            // η ≥ 0 for these rays and Θ(θ) = (p_θ / E)² ≥ 0 at the start.
            let theta_potential = e1 + a * a * theta.cos().powi(2)
                - l1 * l1 * theta.cos().powi(2) / theta.sin().powi(2);
            assert!(theta_potential >= -1e-9, "{theta_potential}");
            assert!(radial_potential(r, r_s, a, l1, e1) >= 0.0);
        }
    }
}

//! Integrator tests against exact answers.

use crate::geometry::atlases::Atlas;
use crate::geometry::atlases::kerr::KerrChart;
use crate::geometry::capture::{critical_parameters, photon_captured, photon_orbit_range};
use crate::geometry::coordinates::boyer_lindquist::BoyerLindquistData;
use crate::geometry::coordinates::kerr_schild::KerrSchildData;
use crate::geometry::four_vector::FourVector;
use crate::geometry::geometry::{Geometry, InnerProduct};
use crate::geometry::point::Point;
use crate::geometry::spacetimes::kerr::Kerr;
use crate::geometry::spacetimes::kerr_bl::{KerrBL, metric_bl};
use crate::rendering::integration::integrator::{IntegrationConfiguration, Integrator, StopReason};
use crate::rendering::ray::Ray;
use nalgebra::Vector4;

/// A photon at BL (r0, θ0, φ = 0) with E = 1, λ, η, moving inward and with
/// p_θ of sign `theta_sign`. None where Θ(θ0) < 0.
fn kerr_photon(
    r_s: f64,
    a: f64,
    lambda: f64,
    eta: f64,
    r0: f64,
    theta0: f64,
    theta_sign: f64,
) -> Option<FourVector<KerrChart>> {
    let (c, s) = (theta0.cos(), theta0.sin());
    let theta_potential = eta + a * a * c * c - lambda * lambda * c * c / (s * s);
    let p = r0 * r0 + a * a - a * lambda;
    let delta = r0 * r0 - r_s * r0 + a * a;
    let radial_potential = p * p - delta * (eta + (lambda - a).powi(2));
    if theta_potential < 0.0 || radial_potential < 0.0 {
        return None;
    }
    let p_cov = Vector4::new(
        -1.0,
        -radial_potential.sqrt() / delta,
        theta_sign * theta_potential.sqrt(),
        lambda,
    );
    let g_inv = metric_bl(r_s, a, r0, theta0).try_inverse()?;
    let v = g_inv * p_cov;
    let chart = KerrChart::BoyerLindquist(BoyerLindquistData { a, r_s });
    Some(FourVector::new(
        v[0],
        v[1],
        v[2],
        v[3],
        Point::new(0.0, r0, theta0, 0.0, chart),
    ))
}

fn captured<G: Geometry>(geometry: &G, momentum: FourVector<G::Chart>, epsilon: f64) -> bool {
    let configuration = IntegrationConfiguration::new(200_000, 200.0, 0.01, epsilon);
    let integrator = Integrator::new(geometry, configuration);
    let ray = Ray::new(0, 0, momentum.at, momentum);
    let (_, stop) = integrator.integrate(&ray).expect("integration");
    stop != Some(StopReason::CelestialSphereReached)
}

/// Of the rays just inside (`1 - delta`) and outside (`1 + delta`) the
/// critical curve of a = 0.499 r_s, scaled from Bardeen's (λ, η) for
/// `orbits` spherical photon orbits and started at r = 20 at the polar
/// angles `thetas`, how many the integrator classifies differently from the
/// exact solution.
fn kerr_shadow_edge_misses<G: Geometry<Chart = KerrChart>>(
    geometry: &G,
    chart: KerrChart,
    orbits: usize,
    thetas: &[f64],
    delta: f64,
    epsilon: f64,
) -> (usize, usize) {
    let (r_s, a) = (1.0, 0.499);
    let (r_min, r_max) = photon_orbit_range(r_s, a);
    let (mut rays, mut misses) = (0, 0);
    for k in 0..orbits {
        let r_orbit = r_min + (r_max - r_min) * (k as f64 + 0.5) / orbits as f64;
        let (lambda, eta) = critical_parameters(r_s, a, r_orbit);
        for (scale, expected) in [(1.0 - delta, true), (1.0 + delta, false)] {
            // Scaling (λ, η) by (s, s²) scales the image point on the sky of a
            // distant observer by s, and the shadow is convex around the origin.
            let (l, e) = (lambda * scale, eta * scale * scale);
            assert_eq!(photon_captured(r_s, a, l, e, 20.0, false), Some(expected));
            for &theta0 in thetas {
                let Some(m) = kerr_photon(r_s, a, l, e, 20.0, theta0, 1.0) else {
                    continue;
                };
                let m = m.at.chart.transition_four_vector(&m, chart);
                assert_eq!(geometry.analytic_capture(&m.at, &m), Some(expected));
                rays += 1;
                if captured(geometry, m, epsilon) != expected {
                    misses += 1;
                }
            }
        }
    }
    (rays, misses)
}

/// Near the shadow's edge the integrator must reach the exact verdict at the
/// default tolerance, in both Kerr charts. delta = 1e-3 is about one pixel at
/// 4320 rows in the README framing (camera at r = 17). Measured misses with 16
/// orbits and both inclinations (62 rays), BL / KS: at delta 1e-3, 22 / 0 at
/// epsilon 1e-5 and 0 / 0 at 1e-9; at delta 1e-7, 31 / 24 at 1e-5 and 31 / 0
/// at 1e-9. So the BL solver in Mino time is the less accurate one this close
/// to the photon orbit: at 1e-9 it gets one pixel right, but lets rays from
/// within ~1e-5 inside the curve escape, where KS still gets 1e-7 right. The
/// test uses fewer rays, KS fewest: debug builds integrate it slowly.
#[test]
fn kerr_shadow_edge_matches_the_exact_solution_at_default_epsilon() {
    let (r_s, a) = (1.0, 0.499);
    let bl = KerrBL::new(r_s, a, 1e-4);
    let ks = Kerr::new(r_s, a, 1e-4);
    let bl_chart = bl.boyer_lindquist_chart();
    let ks_chart = KerrChart::KerrSchild(KerrSchildData { a, r_s });
    let both = [std::f64::consts::FRAC_PI_2, 1.2];
    let equator = [std::f64::consts::FRAC_PI_2];
    let default = crate::cli::cli::DEFAULT_EPSILON;

    let (rays, misses) = kerr_shadow_edge_misses(&bl, bl_chart, 8, &both, 1e-3, default);
    assert!(rays >= 30, "{rays}");
    assert_eq!(misses, 0, "BL, delta 1e-3");
    let (_, misses) = kerr_shadow_edge_misses(&ks, ks_chart, 4, &equator, 1e-7, default);
    assert_eq!(misses, 0, "KS, delta 1e-7");

    // The test can fail: at the old default both charts miss some.
    let (_, misses) = kerr_shadow_edge_misses(&bl, bl_chart, 8, &both, 1e-3, 1e-5);
    assert!(misses >= 5, "BL at 1e-5: {misses}");
    let (_, misses) = kerr_shadow_edge_misses(&ks, ks_chart, 4, &equator, 1e-7, 1e-5);
    assert!(misses >= 2, "KS at 1e-5: {misses}");
}

/// ∫_a^b f by adaptive Simpson to absolute tolerance `tol`.
fn adaptive_simpson(f: &dyn Fn(f64) -> f64, a: f64, b: f64, tol: f64) -> f64 {
    fn simpson(f: &dyn Fn(f64) -> f64, a: f64, fa: f64, b: f64, fb: f64) -> (f64, f64, f64) {
        let m = 0.5 * (a + b);
        let fm = f(m);
        (m, fm, (b - a) / 6.0 * (fa + 4.0 * fm + fb))
    }
    #[allow(clippy::too_many_arguments)]
    fn recurse(
        f: &dyn Fn(f64) -> f64,
        a: f64,
        fa: f64,
        b: f64,
        fb: f64,
        m: f64,
        fm: f64,
        whole: f64,
        tol: f64,
        depth: u32,
    ) -> f64 {
        let (lm, flm, left) = simpson(f, a, fa, m, fm);
        let (rm, frm, right) = simpson(f, m, fm, b, fb);
        let delta = left + right - whole;
        if depth == 0 || delta.abs() <= 15.0 * tol {
            return left + right + delta / 15.0;
        }
        recurse(f, a, fa, m, fm, lm, flm, left, 0.5 * tol, depth - 1)
            + recurse(f, m, fm, b, fb, rm, frm, right, 0.5 * tol, depth - 1)
    }
    let (fa, fb) = (f(a), f(b));
    let (m, fm, whole) = simpson(f, a, fa, b, fb);
    recurse(f, a, fa, b, fb, m, fm, whole, tol, 50)
}

/// Exact azimuth a Schwarzschild photon with impact parameter b sweeps
/// between radius r_a, periapsis and radius r_b: the integral of
/// du / √(1/b² − u² + r_s u³) from 1/r_a up to u0 = 1/periapsis and back to
/// 1/r_b. With u = u0 (1 − t²) the endpoint singularity at u0 factors out:
/// 1/b² − u² + r_s u³ = u0 t² g(u), g(u) = u + u0 − r_s (u² + u u0 + u0²),
/// so each leg is 2√u0 ∫ dt / √g over t ∈ [0, √(1 − u_a/u0)].
fn schwarzschild_swept_angle(r_s: f64, b: f64, r_a: f64, r_b: f64) -> f64 {
    // Periapsis: the largest root of r³ − b² r + b² r_s.
    let periapsis =
        2.0 * b / 3.0_f64.sqrt() * ((-1.5 * 3.0_f64.sqrt() * r_s / b).acos() / 3.0).cos();
    let u0 = 1.0 / periapsis;
    let g = |t: f64| {
        let u = u0 * (1.0 - t * t);
        1.0 / (u + u0 - r_s * (u * u + u * u0 + u0 * u0)).sqrt()
    };
    let leg =
        |r: f64| 2.0 * u0.sqrt() * adaptive_simpson(&g, 0.0, (1.0 - 1.0 / (r * u0)).sqrt(), 1e-14);
    leg(r_a) + leg(r_b)
}

/// Integrate an equatorial Schwarzschild photon with impact parameter b from
/// r = `r_start` inwards until it is back past r_start; returns the azimuth
/// it swept and its final radius.
fn integrated_swept_angle(r_s: f64, b: f64, r_start: f64, epsilon: f64) -> (f64, f64) {
    use crate::geometry::spacetimes::schwarzschild::Schwarzschild;
    let geometry = Schwarzschild::new(r_s, 1e-4);
    let x = Point::new_spherical(0.0, r_start, std::f64::consts::FRAC_PI_2, 0.0);
    let f = 1.0 - r_s / r_start;
    // E = 1, L = b: null by construction.
    let p_r = -(1.0 - f * b * b / (r_start * r_start)).sqrt();
    let momentum = FourVector::new(1.0 / f, p_r, 0.0, b / (r_start * r_start), x);
    let configuration = IntegrationConfiguration::new(1_000_000, r_start, 0.01, epsilon);
    let (ray, stop) = Integrator::new(&geometry, configuration)
        .integrate(&Ray::new(0, 0, x, momentum))
        .expect("integration");
    assert_eq!(stop, Some(StopReason::CelestialSphereReached), "b = {b}");
    let mut swept = 0.0;
    for pair in ray.steps.windows(2) {
        let (p, q) = (
            pair[0].x.get_spatial_vector_cartesian(),
            pair[1].x.get_spatial_vector_cartesian(),
        );
        swept += (p.x * q.y - p.y * q.x).atan2(p.dot(&q));
    }
    let r_end = ray.steps.last().expect("steps").x.get_as_spherical()[0];
    (swept, r_end)
}

/// The azimuth an equatorial Schwarzschild ray sweeps from r = 1000 past the
/// hole and back out matches the exact integral, for impact parameters from
/// 0.1% above the critical 3√3 M out to the weak field. Measured errors at
/// epsilon 1e-9: 1.6e-7 at 1.001 b_c, 1.6e-8 at 1.01 b_c, under 1e-9 from
/// b = 3 r_s on. At 1e-5 they are 7e-4, 6e-5, 8e-7, 7e-7, 1.6e-8 for
/// b = 1.001 b_c .. 8 r_s, all above the tolerances below.
#[test]
fn schwarzschild_deflection_matches_the_exact_integral() {
    let r_s = 1.0;
    let b_c = 1.5 * 3.0_f64.sqrt() * r_s;

    // The reference itself: in the weak field the total deflection is
    // 2 r_s / b + (15π/16) (r_s / b)² + O(b⁻³), less 2 asin(b / R) for
    // starting and ending at a finite R instead of infinity.
    let (b, far) = (1e4, 1e12);
    let deflection = schwarzschild_swept_angle(r_s, b, far, far) - std::f64::consts::PI;
    let expected = 2.0 * r_s / b + 15.0 * std::f64::consts::PI / 16.0 * (r_s / b).powi(2)
        - 2.0 * (b / far).asin();
    assert!(
        (deflection - expected).abs() < 1e-11,
        "{deflection} vs {expected}"
    );
    // And one strong-field value: 1.01 b_c sweeps 7.36600885 rad between r = 1000 and r = 1000.
    let tolerances = [
        (b_c * 1.001, 1e-6),
        (b_c * 1.01, 1e-7),
        (3.0, 1e-8),
        (4.0, 1e-8),
        (8.0, 1e-8),
        (30.0, 1e-8),
    ];
    for (b, tolerance) in tolerances {
        let (swept, r_end) =
            integrated_swept_angle(r_s, b, 1000.0, crate::cli::cli::DEFAULT_EPSILON);
        let exact = schwarzschild_swept_angle(r_s, b, 1000.0, r_end);
        assert!(
            (swept - exact).abs() < tolerance,
            "b {b}: {swept} vs {exact}"
        );
    }

    // The test can fail: at the old default every b up to 8 r_s is off by more.
    for &(b, tolerance) in &tolerances[..5] {
        let (swept, r_end) = integrated_swept_angle(r_s, b, 1000.0, 1e-5);
        let exact = schwarzschild_swept_angle(r_s, b, 1000.0, r_end);
        assert!(
            (swept - exact).abs() > tolerance,
            "b {b} at 1e-5: {swept} vs {exact}"
        );
    }
}

/// Schwarzschild's exact classification, off the equator and in a rotated
/// chart, against the integrator: b = (1 ± 1e-3) 3√3 M from r = 50.
#[test]
fn schwarzschild_capture_matches_the_integrator_off_the_critical_impact_parameter() {
    use crate::geometry::atlases::spherical::SphericalChart;
    use crate::geometry::coordinates::spherical::RotatedSpherical;
    use crate::geometry::spacetimes::schwarzschild::Schwarzschild;
    let r_s = 2.0;
    let geometry = Schwarzschild::new(r_s, 1e-4);
    let b_c = 1.5 * 3.0_f64.sqrt() * r_s;
    let chart = SphericalChart::Rotated(RotatedSpherical { rotation: 90.0 });
    let (r, theta) = (50.0, 1.0_f64);
    let x = Point::new(0.0, r, theta, 0.4, chart);
    for (b, expected) in [(b_c * (1.0 - 1e-3), true), (b_c * (1.0 + 1e-3), false)] {
        let f = 1.0 - r_s / r;
        // Split L = b between the θ and φ directions.
        let (l_theta, l_phi) = (0.6 * b, 0.8 * b);
        let p_r = -(1.0 - f * b * b / (r * r)).sqrt();
        let momentum = FourVector::new(
            1.0 / f,
            p_r,
            l_theta / (r * r),
            l_phi / (r * r * theta.sin()),
            x,
        );
        assert!(geometry.inner_product(&x, &momentum, &momentum).abs() < 1e-12);
        assert_eq!(
            geometry.analytic_capture(&x, &momentum),
            Some(expected),
            "b {b}"
        );
        let reversed = FourVector::new(-momentum[0], -momentum[1], -momentum[2], -momentum[3], x);
        // Traced outward from r = 50 it escapes, whatever b.
        assert_eq!(
            geometry.analytic_capture(&x, &reversed),
            Some(false),
            "b {b}"
        );
        assert_eq!(
            captured(&geometry, momentum, crate::cli::cli::DEFAULT_EPSILON),
            expected,
            "b {b}"
        );
    }
}

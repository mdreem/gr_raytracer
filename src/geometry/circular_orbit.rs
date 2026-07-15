//! Equatorial circular-orbit quantities shared by all black-hole geometries
//! and the disc temperature model.
//!
//! Everything here is a scalar function of (r_s, a, r) in Boyer-Lindquist
//! terms; the geometries act as chart adapters (they extract r from a Point in
//! their own coordinates and, where a coordinate four-velocity is needed,
//! assemble it from the Killing decomposition returned here).
//!
//! Conventions: geometric units with the Schwarzschild radius r_s = 2M, spin
//! parameter `a` signed (positive = prograde orbits for positive Omega).
//! Off-equatorial callers use these values as the standard near-equatorial
//! approximation, matching the previous per-geometry implementations.

use crate::rendering::raytracer::RaytracerError;
use log::debug;

/// A circular-orbit four-velocity decomposed on the Killing basis:
/// u = u_t * d_t + u_phi * d_phi.
///
/// The coefficients are chart-independent scalars. Together with the ray's
/// conserved covariant components (p_t, p_phi) they give the emitter energy
/// at any point without parallel transport:
///
///     u . p = u_t * p_t + u_phi * p_phi
///
/// Deliberately NOT an assembled coordinate four-vector: recovering the
/// coefficients from one is chart-dependent (in Kerr-Schild Cartesian the
/// vector is (u^t, -y u^phi, x u^phi, 0)), while the coefficients pair
/// directly with the conserved quantities.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OrbitKillingDecomposition {
    /// Coefficient of the time Killing vector d_t.
    pub u_t: f64,
    /// Coefficient of the axial Killing vector d_phi (= Omega * u_t).
    pub u_phi: f64,
}

/// BL metric components (g_tt, g_tphi, g_phiphi) at (r, theta).
fn metric_components_at(r_s: f64, a: f64, r: f64, theta: f64) -> (f64, f64, f64) {
    let sig = r * r + a * a * theta.cos().powi(2);
    let sin2 = theta.sin().powi(2);
    let g_tt = -(1.0 - r_s * r / sig);
    let g_tphi = -a * r_s * r * sin2 / sig;
    let g_phiphi = (r * r + a * a + a * a * r_s * r * sin2 / sig) * sin2;
    (g_tt, g_tphi, g_phiphi)
}

/// Equatorial BL metric components (theta = pi/2).
fn metric_components(r_s: f64, a: f64, r: f64) -> (f64, f64, f64) {
    metric_components_at(r_s, a, r, std::f64::consts::FRAC_PI_2)
}

/// Killing coefficients (u^t, u^phi) of the ZAMO (zero angular momentum
/// observer / locally non-rotating frame) at (r, theta). Like the circular
/// orbit, the ZAMO is a pure Killing combination u = u^t d_t + u^phi d_phi
/// with u^phi = omega_zamo * u^t and omega_zamo = -g_tphi / g_phiphi, so the
/// same chart-independent pairing with conserved (p_t, p_phi) applies.
/// Exists everywhere outside the horizon, including inside the ergosphere.
pub fn zamo_killing_coefficients(
    r_s: f64,
    a: f64,
    r: f64,
    theta: f64,
) -> OrbitKillingDecomposition {
    let (g_tt, g_tphi, g_phiphi) = metric_components_at(r_s, a, r, theta);
    let omega = -g_tphi / g_phiphi;
    let u_t = (-1.0 / (g_tt + 2.0 * g_tphi * omega + g_phiphi * omega * omega)).sqrt();
    OrbitKillingDecomposition {
        u_t,
        u_phi: omega * u_t,
    }
}

/// Coordinate angular velocity Omega = dphi/dt of a prograde circular orbit.
/// https://arxiv.org/abs/1104.5499 equation (36).
pub fn angular_velocity(r_s: f64, a: f64, r: f64) -> f64 {
    let m = 0.5 * r_s;
    let sqrt_m = m.sqrt();
    sqrt_m / (r.powf(1.5) + a * sqrt_m)
}

/// Killing coefficients (u^t, u^phi) of the circular orbit at radius r.
/// Errors where no timelike circular orbit exists.
pub fn killing_coefficients(
    r_s: f64,
    a: f64,
    r: f64,
) -> Result<OrbitKillingDecomposition, RaytracerError> {
    let omega = angular_velocity(r_s, a, r);
    let (g_tt, g_tphi, g_phiphi) = metric_components(r_s, a, r);

    let ut_pre = g_tt + 2.0 * omega * g_tphi + omega * omega * g_phiphi;
    if ut_pre >= 0.0 {
        // debug! not error!: called per raymarch sample, and radii inside the
        // photon orbit legitimately have no timelike orbit.
        debug!(
            "No timelike circular orbit at r = {} (ut_pre = {})",
            r, ut_pre
        );
        return Err(RaytracerError::NoCircularOrbitPossible);
    }

    let u_t = (-ut_pre).sqrt().recip();
    Ok(OrbitKillingDecomposition {
        u_t,
        u_phi: omega * u_t,
    })
}

/// Conserved specific energy E = -u_t(covariant) of the circular orbit.
pub fn conserved_energy(r_s: f64, a: f64, r: f64) -> Result<f64, RaytracerError> {
    let omega = angular_velocity(r_s, a, r);
    let (g_tt, g_tphi, _) = metric_components(r_s, a, r);
    let c = killing_coefficients(r_s, a, r)?;
    Ok(-(g_tt + g_tphi * omega) * c.u_t)
}

/// Conserved specific angular momentum L = u_phi(covariant) of the circular orbit.
pub fn conserved_angular_momentum(r_s: f64, a: f64, r: f64) -> Result<f64, RaytracerError> {
    let omega = angular_velocity(r_s, a, r);
    let (_, g_tphi, g_phiphi) = metric_components(r_s, a, r);
    let c = killing_coefficients(r_s, a, r)?;
    Ok((g_tphi + g_phiphi * omega) * c.u_t)
}

/// Prograde ISCO radius (innermost stable circular orbit).
pub fn r_isco(r_s: f64, a: f64) -> f64 {
    let a_s = 2.0 * a / r_s;

    let z1 = 1.0
        + (1.0 - a_s * a_s).powf(1.0 / 3.0)
            * ((1.0 + a_s).powf(1.0 / 3.0) + (1.0 - a_s).powf(1.0 / 3.0));
    let z2 = (3.0 * a_s * a_s + z1 * z1).sqrt();

    (3.0 + z2 - ((3.0 - z1) * (3.0 + z1 + 2.0 * z2)).sqrt()) * r_s / 2.0
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::four_vector::FourVector;
    use crate::geometry::geometry::{Geometry, SupportQuantities};
    use crate::geometry::kerr::Kerr;
    use crate::geometry::kerr_bl::KerrBL;
    use crate::geometry::point::{CoordinateSystem, Point};
    use crate::geometry::schwarzschild::Schwarzschild;
    use approx::assert_abs_diff_eq;
    use std::f64::consts::PI;

    #[test]
    fn test_r_isco_known_values() {
        // Schwarzschild: ISCO = 6M = 3 r_s.
        assert_abs_diff_eq!(r_isco(1.0, 0.0), 3.0, epsilon = 1e-12);
        // Near-extremal prograde: ISCO approaches M = 0.5 r_s.
        assert!(r_isco(1.0, 0.499) < 0.63);
        assert!(r_isco(1.0, 0.499) > 0.5);
    }

    #[test]
    fn test_schwarzschild_limit_matches_closed_forms() {
        let (r_s, r) = (1.0_f64, 5.0_f64);
        let m = 0.5 * r_s;
        assert_abs_diff_eq!(
            angular_velocity(r_s, 0.0, r),
            (m / (r * r * r)).sqrt(),
            epsilon = 1e-14
        );
        let c = killing_coefficients(r_s, 0.0, r).unwrap();
        assert_abs_diff_eq!(c.u_t, (1.0 - 3.0 * m / r).sqrt().recip(), epsilon = 1e-14);
    }

    #[test]
    fn test_no_timelike_orbit_inside_photon_sphere() {
        // Schwarzschild photon sphere at 1.5 r_s: no timelike circular orbit at/below.
        assert!(killing_coefficients(1.0, 0.0, 1.4).is_err());
        assert!(killing_coefficients(1.0, 0.0, 1.6).is_ok());
    }

    /// ZAMO properties in both Kerr charts: normalized, zero angular
    /// momentum (the defining property), identical Killing coefficients
    /// across charts, and reduction to the static observer at a = 0.
    #[test]
    fn test_zamo_properties_across_charts() {
        let (r_s, a) = (1.0, 0.499);
        let kerr = Kerr::new(r_s, a, 1e-4);
        let kerr_bl = KerrBL::new(r_s, a, 1e-4);

        // Same physical point: equatorial, BL r = 5, phi = 0
        // (x = st(r cp - a sp) etc. with theta = pi/2).
        let pos_bl = Point::new(
            0.0,
            5.0,
            PI / 2.0,
            0.0,
            CoordinateSystem::BoyerLindquist { a },
        );
        let pos_ks = Point::new_cartesian(0.0, 5.0, a, 0.0);

        for (zamo, axial, pos, geom_ip) in [
            (
                kerr.get_zamo_velocity_at(&pos_ks),
                kerr.axial_killing_vector(&pos_ks),
                &pos_ks,
                &kerr as &dyn Geometry,
            ),
            (
                kerr_bl.get_zamo_velocity_at(&pos_bl),
                kerr_bl.axial_killing_vector(&pos_bl),
                &pos_bl,
                &kerr_bl as &dyn Geometry,
            ),
        ] {
            assert_abs_diff_eq!(
                geom_ip.inner_product(pos, &zamo, &zamo),
                -1.0,
                epsilon = 1e-9
            );
            assert_abs_diff_eq!(
                geom_ip.inner_product(pos, &zamo, &axial),
                0.0,
                epsilon = 1e-9
            );
        }

        // Chart-invariant coefficients.
        let c = zamo_killing_coefficients(r_s, a, 5.0, PI / 2.0);
        let zamo_bl = kerr_bl.get_zamo_velocity_at(&pos_bl);
        assert_abs_diff_eq!(zamo_bl[0], c.u_t, epsilon = 1e-12);
        assert_abs_diff_eq!(zamo_bl[3], c.u_phi, epsilon = 1e-12);

        // a = 0: ZAMO reduces to the static observer.
        let c0 = zamo_killing_coefficients(1.0, 0.0, 5.0, 1.1);
        assert_abs_diff_eq!(c0.u_phi, 0.0, epsilon = 1e-15);
        let schwarzschild = Schwarzschild::new(1.0, 1e-4);
        let pos_s = Point::new_spherical(0.0, 5.0, 1.1, 0.3);
        let zamo_s = schwarzschild.get_zamo_velocity_at(&pos_s);
        let stat_s = schwarzschild.get_stationary_velocity_at(&pos_s);
        for i in 0..4 {
            assert_abs_diff_eq!(zamo_s[i], stat_s[i], epsilon = 1e-15);
        }
    }

    /// The assembled coordinate four-velocity must equal
    /// u_t * d_t + u_phi * axial_killing_vector, be normalized (u.u = +/-1
    /// depending on signature), and satisfy the emitter-energy identity
    /// u.p = u_t * p_t + u_phi * p_phi for arbitrary vectors p, in every
    /// geometry/chart.
    #[test]
    fn test_killing_decomposition_identities_across_geometries() {
        fn check<G: Geometry>(
            geometry: &G,
            position: &Point,
            probe_momenta: &[FourVector],
            time_sign: f64, // signature()[0]: u.u == time_sign
        ) {
            let u_vec = geometry.get_circular_orbit_velocity_at(position).unwrap();
            let c = geometry
                .circular_orbit_killing_coefficients(position)
                .unwrap();
            let cs = u_vec.coordinate_system;
            let e_t = FourVector::new(1.0, 0.0, 0.0, 0.0, cs);
            let axial = geometry.axial_killing_vector(position);

            // Assembly identity.
            let assembled = c.u_t * e_t + c.u_phi * axial;
            for i in 0..4 {
                assert_abs_diff_eq!(u_vec[i], assembled[i], epsilon = 1e-12);
            }

            // Normalization.
            assert_abs_diff_eq!(
                geometry.inner_product(position, &u_vec, &u_vec),
                time_sign,
                epsilon = 1e-10
            );

            // Emitter-energy identity for arbitrary probe vectors.
            for p in probe_momenta {
                let p_t = geometry.inner_product(position, &e_t, p);
                let p_phi = geometry.inner_product(position, &axial, p);
                let direct = geometry.inner_product(position, &u_vec, p);
                assert_abs_diff_eq!(direct, c.u_t * p_t + c.u_phi * p_phi, epsilon = 1e-10);
            }
        }

        // Schwarzschild (spherical chart, +--- signature).
        let schwarzschild = Schwarzschild::new(1.0, 1e-4);
        let pos_s = Point::new_spherical(0.0, 6.0, PI / 2.0, 1.1);
        let probes_s = [
            FourVector::new_spherical(1.0, -0.4, 0.02, 0.05),
            FourVector::new_spherical(-1.3, 0.2, 0.0, -0.08),
        ];
        check(&schwarzschild, &pos_s, &probes_s, 1.0);

        // Kerr (Kerr-Schild Cartesian chart, -+++ signature).
        let kerr = Kerr::new(1.0, 0.499, 1e-4);
        let pos_k = Point::new_cartesian(0.0, 3.0, -4.0, 0.0);
        let probes_k = [
            FourVector::new_cartesian(1.0, 0.3, -0.2, 0.1),
            FourVector::new_cartesian(-0.8, -0.5, 0.4, 0.0),
        ];
        check(&kerr, &pos_k, &probes_k, -1.0);

        // KerrBL (Boyer-Lindquist chart, -+++ signature).
        let a = 0.3;
        let kerr_bl = KerrBL::new(1.0, a, 1e-4);
        let pos_bl = Point::new(
            0.0,
            5.0,
            PI / 2.0,
            0.7,
            CoordinateSystem::BoyerLindquist { a },
        );
        let probes_bl = [
            FourVector::new_boyer_lindquist(a, 1.0, -0.4, 0.02, 0.05),
            FourVector::new_boyer_lindquist(a, -1.2, 0.1, -0.03, 0.09),
        ];
        check(&kerr_bl, &pos_bl, &probes_bl, -1.0);
    }
}

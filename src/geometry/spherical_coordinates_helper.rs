use crate::geometry::point::{CoordinateSystem, Point};
use log::trace;
use std::sync::OnceLock;

/// Fixed rotation (degrees, about the x-axis) inserted between physical
/// Cartesian space and the integration's spherical frame, so the coordinate
/// pole (theta = 0/pi, where 1/sin(theta) blows up) can be steered off the
/// field of view. 0 (the default) is the identity and leaves every render
/// byte-identical; set POLE_ROT_DEG to move the pole. Schwarzschild is
/// spherically symmetric, so this cannot change the physics, only where the
/// coordinate singularity sits.
pub fn frame_rotation_deg() -> f64 {
    static DEG: OnceLock<f64> = OnceLock::new();
    *DEG.get_or_init(|| {
        std::env::var("POLE_ROT_DEG")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(0.0)
    })
}

/// Rotate the spatial vector `(x, y, z)` about the x-axis by `deg` degrees.
pub fn rotate_about_x(x: f64, y: f64, z: f64, deg: f64) -> (f64, f64, f64) {
    if deg == 0.0 {
        return (x, y, z);
    }
    let (s, c) = deg.to_radians().sin_cos();
    (x, y * c - z * s, y * s + z * c)
}

// The order of the components is: (r, theta, phi)
pub fn cartesian_to_spherical(cartesian: &Point) -> Point {
    let t = cartesian[0];
    // Physical -> coordinate frame: rotate so the pole moves off the field.
    let (x, y, z) = rotate_about_x(
        cartesian[1],
        cartesian[2],
        cartesian[3],
        frame_rotation_deg(),
    );

    let r = (x * x + y * y + z * z).sqrt();
    if r == 0.0 {
        return Point::new(t, 0.0, 0.0, 0.0, CoordinateSystem::Spherical);
    }

    let theta = (z / r).acos();
    let phi = y.atan2(x);

    trace!("Converting Cartesian to Spherical:");
    trace!("  Cartesian: t={}, x={}, y={}, z={}", t, x, y, z);
    trace!(
        "  Spherical: t={}, r={}, theta={}, phi={}",
        t, r, theta, phi
    );
    Point::new(t, r, theta, phi, CoordinateSystem::Spherical)
}

pub fn spherical_to_cartesian(spherical: &Point) -> Point {
    let t = spherical[0];
    let r = spherical[1];
    let theta = spherical[2];
    let phi = spherical[3];

    let x = r * theta.sin() * phi.cos();
    let y = r * theta.sin() * phi.sin();
    let z = r * theta.cos();

    // Coordinate -> physical frame: undo the pole rotation.
    let (x, y, z) = rotate_about_x(x, y, z, -frame_rotation_deg());
    Point::new(t, x, y, z, CoordinateSystem::Cartesian)
}

/// Convert a Cartesian point to Boyer-Lindquist (t, r, θ, φ) coordinates for
/// spin parameter `a`, using the Kerr-Schild embedding documented on
/// `CoordinateSystem::BoyerLindquist`.
pub fn cartesian_to_boyer_lindquist(a: f64, cartesian: &Point) -> Point {
    let t = cartesian[0];
    let x = cartesian[1];
    let y = cartesian[2];
    let z = cartesian[3];

    let rho_sqr = x * x + y * y + z * z;
    let r_sqr = 0.5 * (rho_sqr - a * a + ((rho_sqr - a * a).powi(2) + 4.0 * a * a * z * z).sqrt());
    let r = r_sqr.sqrt();
    let theta = if r == 0.0 {
        0.0
    } else {
        (z / r).clamp(-1.0, 1.0).acos()
    };
    let phi = (r * y - a * x).atan2(r * x + a * y);

    Point::new(t, r, theta, phi, CoordinateSystem::BoyerLindquist { a })
}

#[cfg(test)]
mod tests {
    use crate::geometry::point::{CoordinateSystem, Point};
    use crate::geometry::spherical_coordinates_helper::{
        cartesian_to_spherical, spherical_to_cartesian,
    };
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_cartesian_to_spherical() {
        let cartesian = Point::new(0.0, 1.0, 2.0, 3.0, CoordinateSystem::Cartesian);
        let spherical = cartesian_to_spherical(&cartesian);
        let back_to_cartesian = spherical_to_cartesian(&spherical);

        assert_abs_diff_eq!(
            cartesian.get_as_vector(),
            back_to_cartesian.get_as_vector(),
            epsilon = 1e-15
        );
    }
}

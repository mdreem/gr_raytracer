use crate::geometry::coordinates;
use crate::geometry::point::Point;
use log::trace;

/// Rotate the spatial vector `(x, y, z)` about the x-axis by `deg` degrees.
pub fn rotate_about_x(x: f64, y: f64, z: f64, deg: f64) -> (f64, f64, f64) {
    if deg == 0.0 {
        return (x, y, z);
    }
    let (s, c) = deg.to_radians().sin_cos();
    (x, y * c - z * s, y * s + z * c)
}

// The order of the components is: (r, theta, phi)
pub fn cartesian_to_spherical(
    cartesian: &Point<coordinates::euclidean::CartesianChart>,
) -> Point<coordinates::spherical::SphericalChart> {
    let t = cartesian[0];
    let (x, y, z) = (cartesian[1], cartesian[2], cartesian[3]);

    let r = (x * x + y * y + z * z).sqrt();
    if r == 0.0 {
        return Point::new_spherical(t, 0.0, 0.0, 0.0);
    }

    let theta = (z / r).acos();
    let phi = y.atan2(x);

    trace!("Converting Cartesian to Spherical:");
    trace!("  Cartesian: t={}, x={}, y={}, z={}", t, x, y, z);
    trace!(
        "  Spherical: t={}, r={}, theta={}, phi={}",
        t, r, theta, phi
    );
    Point::new_spherical(t, r, theta, phi)
}

/// Physical Cartesian position of a point given in a (possibly rotated) spherical chart.
pub fn spherical_to_cartesian(
    spherical: &Point<coordinates::spherical::SphericalChart>,
) -> Point<coordinates::euclidean::CartesianChart> {
    let t = spherical[0];
    let r = spherical[1];
    let theta = spherical[2];
    let phi = spherical[3];

    let x = r * theta.sin() * phi.cos();
    let y = r * theta.sin() * phi.sin();
    let z = r * theta.cos();

    // Coordinate -> physical frame: undo the chart's rotation.
    let (x, y, z) = rotate_about_x(x, y, z, -spherical.chart.rotation());
    Point::new_cartesian(t, x, y, z)
}

#[cfg(test)]
mod tests {
    use crate::geometry::point::Point;
    use crate::geometry::spherical_coordinates_helper::{
        cartesian_to_spherical, spherical_to_cartesian,
    };
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_cartesian_to_spherical() {
        let cartesian = Point::new_cartesian(0.0, 1.0, 2.0, 3.0);
        let spherical = cartesian_to_spherical(&cartesian);
        let back_to_cartesian = spherical_to_cartesian(&spherical);

        assert_abs_diff_eq!(
            cartesian.get_as_vector(),
            back_to_cartesian.get_as_vector(),
            epsilon = 1e-15
        );
    }
}

use crate::geometry::point::{CoordinateSystem, Point};

// The order of the components is: (r, theta, phi)
pub fn cartesian_to_spherical(cartesian: &Point) -> Point {
    let t = cartesian[0];
    let x = cartesian[1];
    let y = cartesian[2];
    let z = cartesian[3];

    let r = (x * x + y * y + z * z).sqrt();
    if r == 0.0 {
        return Point::new(t, 0.0, 0.0, 0.0, CoordinateSystem::Spherical);
    }

    let theta = (y / r).acos();
    let phi = z.atan2(x);

    Point::new(t, r, theta, phi, CoordinateSystem::Spherical)
}

pub fn spherical_to_cartesian(spherical: &Point) -> Point {
    let t = spherical[0];
    let r = spherical[1];
    let theta = spherical[2];
    let phi = spherical[3];

    let x = r * theta.sin() * phi.cos();
    let y = r * theta.cos();
    let z = r * theta.sin() * phi.sin();

    Point::new(t, x, y, z, CoordinateSystem::Cartesian)
}

#[cfg(test)]
mod tests {
    use crate::geometry::point::{CoordinateSystem, Point};
    use crate::geometry::spherical_coordinates_helper::{
        cartesian_to_spherical, spherical_to_cartesian,
    };
    use approx::assert_abs_diff_eq;
    use nalgebra::Vector4;

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

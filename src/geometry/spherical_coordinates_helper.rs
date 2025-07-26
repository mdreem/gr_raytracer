use nalgebra::Vector4;

// The order of the components is: (r, theta, phi)
pub fn cartesian_to_spherical(cartesian: &Vector4<f64>) -> Vector4<f64> {
    let t = cartesian[0];
    let x = cartesian[1];
    let y = cartesian[2];
    let z = cartesian[3];

    let r = (x * x + y * y + z * z).sqrt();
    if r == 0.0 {
        return Vector4::new(t, 0.0, 0.0, 0.0);
    }

    let theta = (y / r).acos();
    let phi = z.atan2(x);

    Vector4::new(t, r, theta, phi)
}

pub fn spherical_to_cartesian(spherical: &Vector4<f64>) -> Vector4<f64> {
    let t = spherical[0];
    let r = spherical[1];
    let theta = spherical[2];
    let phi = spherical[3];

    let x = r * theta.sin() * phi.cos();
    let y = r * theta.cos();
    let z = r * theta.sin() * phi.sin();

    Vector4::new(t, x, y, z)
}

#[cfg(test)]
mod tests {
    use crate::geometry::spherical_coordinates_helper::{
        cartesian_to_spherical, spherical_to_cartesian,
    };
    use approx::assert_abs_diff_eq;
    use nalgebra::Vector4;

    #[test]
    fn test_cartesian_to_spherical() {
        let cartesian = Vector4::new(0.0, 1.0, 2.0, 3.0);
        let spherical = cartesian_to_spherical(&cartesian);
        let back_to_cartesian = spherical_to_cartesian(&spherical);

        assert_abs_diff_eq!(cartesian, back_to_cartesian, epsilon = 1e-15);
    }
}

use nalgebra::{Vector3, Vector4};

// The order of the components is: (r, theta, phi)
pub fn cartesian_to_spherical(cartesian: &Vector4<f64>) -> Vector3<f64> {
    let x = cartesian[1];
    let y = cartesian[2];
    let z = cartesian[3];

    let r = (x * x + y * y + z * z).sqrt();
    if r == 0.0 {
        return Vector3::new(0.0, 0.0, 0.0);
    }

    let theta = (y / r).acos();
    let phi = z.atan2(x);

    Vector3::new(r, theta, phi)
}

pub fn spherical_to_cartesian(spherical: &Vector4<f64>) -> Vector3<f64> {
    let _t = spherical[0];
    let r = spherical[1];
    let theta = spherical[2];
    let phi = spherical[3];

    let x = r * theta.sin() * phi.cos();
    let y = r * theta.cos();
    let z = r * theta.sin() * phi.sin();

    Vector3::new(x, y, z)
}

#[cfg(test)]
mod tests {
    use crate::spherical_coordinates_helper::{cartesian_to_spherical, spherical_to_cartesian};
    use approx::assert_abs_diff_eq;
    use nalgebra::{Vector3, Vector4};

    #[test]
    fn test_cartesian_to_spherical() {
        let cartesian = Vector3::new(1.0, 2.0, 3.0);
        let spherical =
            cartesian_to_spherical(&Vector4::new(1.0, cartesian[0], cartesian[1], cartesian[2]));
        let back_to_cartesian =
            spherical_to_cartesian(&Vector4::new(1.0, spherical[0], spherical[1], spherical[2]));

        assert_abs_diff_eq!(cartesian, back_to_cartesian, epsilon = 1e-15);
    }
}

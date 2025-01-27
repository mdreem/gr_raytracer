use crate::four_vector::CoordinateSystem::Spherical;
use crate::four_vector::{CoordinateSystem, FourVector};
use crate::geometry::{Geometry, HasCoordinateSystem, Tetrad};
use crate::runge_kutta::OdeFunction;
use crate::scene::EquationOfMotionState;
use nalgebra::{Const, Matrix4, OVector, Vector4};

#[derive(Clone, Debug)]
pub struct Schwarzschild {
    radius: f64,
}

impl Schwarzschild {
    pub fn new(radius: f64) -> Self {
        Schwarzschild { radius }
    }
}

impl OdeFunction<Const<8>> for Schwarzschild {
    // TODO: maybe just have geodesic being used in solver. This doesn't need to be that generic here.
    fn apply(&self, t: f64, y: &OVector<f64, Const<8>>) -> OVector<f64, Const<8>> {
        self.geodesic(t, y)
    }
}

impl HasCoordinateSystem for Schwarzschild {
    fn coordinate_system(&self) -> CoordinateSystem {
        Spherical
    }
}

// All coordinates here are spherical coordinates.
impl Geometry for Schwarzschild {
    fn geodesic(&self, _: f64, y: &EquationOfMotionState) -> EquationOfMotionState {
        let _t = y[0];
        let r = y[1];
        let theta = y[2];
        let _phi = y[3];

        let v_t = y[4];
        let v_r = y[5];
        let v_theta = y[6];
        let v_phi = y[7];

        let a = 1.0 - self.radius / r;
        let a_prime = self.radius / (r * r);
        let aprime_over_a = a_prime / a;

        // acceleration
        let a_t = -(aprime_over_a) * v_t * v_r;
        let a_r = -0.5 * a * a_prime * v_t * v_t
            + 0.5 * (aprime_over_a) * v_r * v_r
            + a * r * (v_theta * v_theta + v_phi * v_phi * theta.sin() * theta.sin());
        let a_theta = -(2.0 / r) * v_r * v_theta + theta.sin() * theta.cos() * v_phi * v_phi;
        let a_phi = -(2.0 / r) * v_phi * v_r - 2.0 * theta.cos() / theta.sin() * v_theta * v_phi;

        let y_new = EquationOfMotionState::from_column_slice(&[
            v_t, v_r, v_theta, v_phi, a_t, a_r, a_theta, a_phi,
        ]);
        y_new
    }

    // TODO: take into account rotations.
    fn get_tetrad_at(&self, position: &Vector4<f64>) -> Tetrad {
        let r = position[1];
        let theta = position[2];
        let _phi = position[3];

        let rr0 = self.radius / r;
        let a = 1.0 - rr0;

        Tetrad::new(
            position.clone(),
            FourVector::new_spherical(1.0 / a, -rr0.sqrt(), 0.0, 0.0),
            FourVector::new_spherical(0.0, 0.0, 0.0, 1.0 / (r * theta.sin())), // Phi
            -FourVector::new_spherical(0.0, 0.0, 1.0 / r, 0.0),                // Theta
            -FourVector::new_spherical(-rr0.sqrt() / a, 1.0, 0.0, 0.0),        // R
        )
    }

    fn lorentz_transformation(
        &self,
        position: &Vector4<f64>,
        velocity: &FourVector,
    ) -> Matrix4<f64> {
        let mut matrix = Matrix4::zeros();
        let tetrad_t = self.get_tetrad_at(&position).t;

        let r = position[1];
        let theta = position[2];
        let _phi = position[3];

        let a = 1.0 - self.radius / r;

        let metric_diag = Vector4::new(a, -1.0 / a, -r * r, -r * r * theta.sin() * theta.sin());

        let mut gamma = 0.0;
        for i in 0..4 {
            gamma += metric_diag[i] * velocity.vector[i] * tetrad_t.vector[i];
        }
        for mu in 0..4 {
            for nu in 0..4 {
                let mut res = 0.0;
                if mu == nu {
                    res = 1.0;
                }

                let a = 1.0 / (1.0 + gamma);
                let b = tetrad_t.vector[mu] + velocity.vector[mu];
                let c = metric_diag[nu] * (tetrad_t.vector[nu] + velocity.vector[nu]);
                res -= a * b * c;

                res += 2.0 * metric_diag[nu] * tetrad_t.vector[nu] * velocity.vector[mu];

                matrix[(mu, nu)] = res;
            }
        }
        matrix
    }

    fn mul(&self, position: &Vector4<f64>, v: &FourVector, w: &FourVector) -> f64 {
        let r = position[1];
        let theta = position[2];
        let _phi = position[3];

        let a = 1.0 - self.radius / r;

        a * v.vector[0] * w.vector[0]
            - v.vector[1] * w.vector[1] / a
            - r * r * v.vector[2] * w.vector[2]
            - r * r * theta.sin() * theta.sin() * v.vector[3] * w.vector[3]
    }
}

#[cfg(test)]
mod tests {
    use crate::camera::Camera;
    use crate::debug::save_rays_to_file;
    use crate::four_vector::FourVector;
    use crate::geometry::Geometry;
    use crate::schwarzschild::Schwarzschild;
    use crate::spherical_coordinates_helper::cartesian_to_spherical;
    use approx::assert_abs_diff_eq;
    use nalgebra::Vector4;
    use std::f64::consts::PI;

    #[test]
    fn test_tetrad_orthonormal() {
        let position = cartesian_to_spherical(&Vector4::new(2.0, 3.0, 4.0, 5.0));
        let geometry = Schwarzschild::new(2.0);

        let tetrad = geometry.get_tetrad_at(&position);

        let k = tetrad.t + (-tetrad.z);
        let s = geometry.mul(&position, &k, &k);
        assert_abs_diff_eq!(s, 0.0);

        assert_abs_diff_eq!(geometry.mul(&position, &tetrad.t, &tetrad.t), 1.0);
        assert_abs_diff_eq!(geometry.mul(&position, &tetrad.x, &tetrad.x), -1.0);
        assert_abs_diff_eq!(geometry.mul(&position, &tetrad.y, &tetrad.y), -1.0);
        assert_abs_diff_eq!(geometry.mul(&position, &tetrad.z, &tetrad.z), -1.0);

        assert_abs_diff_eq!(geometry.mul(&position, &tetrad.t, &tetrad.x), 0.0);
        assert_abs_diff_eq!(geometry.mul(&position, &tetrad.t, &tetrad.y), 0.0);
        assert_abs_diff_eq!(geometry.mul(&position, &tetrad.t, &tetrad.z), 0.0);

        assert_abs_diff_eq!(geometry.mul(&position, &tetrad.x, &tetrad.y), 0.0);
        assert_abs_diff_eq!(geometry.mul(&position, &tetrad.x, &tetrad.z), 0.0);

        assert_abs_diff_eq!(geometry.mul(&position, &tetrad.y, &tetrad.z), 0.0);
    }

    #[test]
    fn test_lorentz_transformed_tetrad_orthonormal() {
        let position = cartesian_to_spherical(&Vector4::new(2.0, 3.0, 4.0, 5.0));
        let radius = 2.0;
        let r = position[1];
        let a = 1.0 - radius / position[1];

        let geometry = Schwarzschild::new(radius);

        let velocity = FourVector::new_spherical(1.0 / a, -(radius / r).sqrt(), 0.0, 0.0); // we have a freely falling observer here.

        let original_tetrad = geometry.get_tetrad_at(&position);
        let tetrad = crate::camera::lorentz_transform_tetrad(
            &geometry,
            &original_tetrad,
            &position,
            &velocity,
        );

        let k = tetrad.t + (-tetrad.z);
        let s = geometry.mul(&position, &k, &k);
        assert_abs_diff_eq!(s, 0.0);

        assert_abs_diff_eq!(
            tetrad.t.get_as_vector(),
            velocity.get_as_vector(),
            epsilon = 1e-6
        );

        assert_abs_diff_eq!(geometry.mul(&position, &tetrad.t, &tetrad.x), 0.0);
        assert_abs_diff_eq!(geometry.mul(&position, &tetrad.t, &tetrad.y), 0.0);
        assert_abs_diff_eq!(geometry.mul(&position, &tetrad.t, &tetrad.z), 0.0);

        assert_abs_diff_eq!(geometry.mul(&position, &tetrad.x, &tetrad.y), 0.0);
        assert_abs_diff_eq!(geometry.mul(&position, &tetrad.x, &tetrad.z), 0.0);

        assert_abs_diff_eq!(geometry.mul(&position, &tetrad.y, &tetrad.z), 0.0);
    }

    #[ignore]
    #[test]
    fn save_camera_rays() {
        let rows = 30;
        let cols = 30;

        let position = cartesian_to_spherical(&Vector4::new(0.0, 0.0, 0.0, -10.0));
        let radius = 2.0;
        let geometry = Schwarzschild::new(radius);
        let r = position[1];
        let a = 1.0 - radius / r;
        let velocity = FourVector::new_spherical(1.0 / a, -(radius / r).sqrt(), 0.0, 0.0); // we have a freely falling observer here.
        let camera = Camera::new(
            position,
            velocity,
            PI / 2.0,
            rows,
            cols,
            Schwarzschild::new(2.0),
        );

        save_rays_to_file(rows, cols, &position, geometry, camera);
    }

    #[test]
    fn test_schwarzschild_ray() {
        let position = cartesian_to_spherical(&Vector4::new(2.0, 3.0, 4.0, 5.0));
        let radius = 2.0;
        let geometry = Schwarzschild::new(radius);
        let r = position[1];
        let a = 1.0 - radius / r;
        let velocity = FourVector::new_spherical(1.0 / a, -(radius / r).sqrt(), 0.0, 0.0); // we have a freely falling observer here.
        let camera = Camera::new(
            position,
            velocity,
            PI / 2.0,
            11,
            11,
            Schwarzschild::new(2.0),
        );

        let ray = camera.get_ray_for(6, 6);
        // This is a space-like vector and therefore should be -1.
        assert_abs_diff_eq!(
            geometry.mul(&position, &ray.direction, &ray.direction),
            -1.0
        );
    }
}

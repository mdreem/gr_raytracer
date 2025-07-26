use crate::geometry::four_vector::CoordinateSystem::Spherical;
use crate::geometry::four_vector::{CoordinateSystem, FourVector};
use crate::geometry::geometry::{
    GeodesicSolver, Geometry, HasCoordinateSystem, InnerProduct, Tetrad,
};
use crate::rendering::runge_kutta::OdeFunction;
use crate::rendering::scene::EquationOfMotionState;
use nalgebra::{Const, Matrix4, OVector, Vector4};

#[derive(Clone)]
pub struct EuclideanSpaceSpherical {}

impl EuclideanSpaceSpherical {
    pub fn new() -> Self {
        EuclideanSpaceSpherical {}
    }
}

impl OdeFunction<Const<8>> for EuclideanSpaceSpherical {
    // TODO: maybe just have geodesic being used in solver. This doesn't need to be that generic here.
    fn apply(&self, t: f64, y: &OVector<f64, Const<8>>) -> OVector<f64, Const<8>> {
        self.geodesic(t, y)
    }
}

impl GeodesicSolver for EuclideanSpaceSpherical {
    fn geodesic(&self, _: f64, y: &EquationOfMotionState) -> EquationOfMotionState {
        let _t = y[0];
        let r = y[1];
        let theta = y[2];
        let _phi = y[3];

        let v_t = y[4];
        let v_r = y[5];
        let v_theta = y[6];
        let v_phi = y[7];

        // acceleration
        let a_t = 0.0;
        let a_r = r * (v_theta * v_theta + v_phi * v_phi * theta.sin() * theta.sin());
        let a_theta = -(2.0 / r) * v_r * v_theta + theta.sin() * theta.cos() * v_phi * v_phi;
        let a_phi = -(2.0 / r) * v_phi * v_r - 2.0 * theta.cos() / theta.sin() * v_theta * v_phi;

        // y'
        let y_new = EquationOfMotionState::from_column_slice(&[
            v_t, v_r, v_theta, v_phi, a_t, a_r, a_theta, a_phi,
        ]);
        y_new
    }
}

impl HasCoordinateSystem for EuclideanSpaceSpherical {
    fn coordinate_system(&self) -> CoordinateSystem {
        Spherical
    }
}

impl InnerProduct for EuclideanSpaceSpherical {
    fn inner_product(&self, position: &Vector4<f64>, v: &FourVector, w: &FourVector) -> f64 {
        let r = position[1];
        let theta = position[2];
        let _phi = position[3];

        1.0 * v.vector[0] * w.vector[0]
            - v.vector[1] * w.vector[1]
            - r * r * v.vector[2] * w.vector[2]
            - r * r * theta.sin() * theta.sin() * v.vector[3] * w.vector[3]
    }
}

impl Geometry for EuclideanSpaceSpherical {
    // TODO: take into account rotations.
    fn get_tetrad_at(&self, position: &Vector4<f64>) -> Tetrad {
        let r = position[1];
        let theta = position[2];
        let _phi = position[3];

        Tetrad::new(
            *position,
            FourVector::new_spherical(1.0, 0.0, 0.0, 0.0),
            FourVector::new_spherical(0.0, 0.0, 0.0, 1.0 / (r * theta.sin())), // Phi
            -FourVector::new_spherical(0.0, 0.0, 1.0 / r, 0.0),                // Theta
            -FourVector::new_spherical(0.0, 1.0, 0.0, 0.0),                    // R
        )
    }

    fn lorentz_transformation(
        &self,
        _position: &Vector4<f64>,
        _velocity: &FourVector,
    ) -> Matrix4<f64> {
        let mut matrix = Matrix4::zeros();
        matrix[(0, 0)] = 1.0;
        matrix[(1, 1)] = 1.0;
        matrix[(2, 2)] = 1.0;
        matrix[(3, 3)] = 1.0;

        matrix
    }

    fn get_stationary_velocity_at(&self, _position: &Vector4<f64>) -> FourVector {
        FourVector::new_spherical(1.0, 0.0, 0.0, 0.0)
    }

    fn inside_horizon(&self, _position: &Vector4<f64>) -> bool {
        false
    }
}

#[cfg(test)]
mod tests {
    use crate::geometry::euclidean_spherical::EuclideanSpaceSpherical;
    use crate::geometry::four_vector::FourVector;
    use crate::geometry::spherical_coordinates_helper::cartesian_to_spherical;
    use crate::rendering::camera::Camera;
    use crate::rendering::debug::save_rays_to_file;
    use nalgebra::Vector4;
    use std::f64::consts::PI;

    #[ignore]
    #[test]
    fn save_camera_rays() {
        let rows = 30;
        let cols = 30;

        let position = cartesian_to_spherical(&Vector4::new(0.0, 0.0, 0.0, -10.0));
        let velocity = FourVector::new_spherical(1.0, 0.0, 0.0, 0.0);
        let camera = Camera::new(
            position,
            velocity,
            PI / 2.0,
            rows,
            cols,
            &EuclideanSpaceSpherical::new(),
        );

        save_rays_to_file(
            rows,
            cols,
            &position,
            EuclideanSpaceSpherical::new(),
            camera,
        );
    }
}

use crate::four_vector::CoordinateSystem::Cartesian;
use crate::four_vector::{CoordinateSystem, FourVector};
use crate::geometry::{Geometry, HasCoordinateSystem, Tetrad};
use crate::runge_kutta::OdeFunction;
use crate::scene::EquationOfMotionState;
use nalgebra::{Const, Matrix4, OVector, Vector4};

#[derive(Clone)]
pub struct EuclideanSpace {}

impl EuclideanSpace {
    pub fn new() -> Self {
        EuclideanSpace {}
    }
}

impl OdeFunction<Const<8>> for EuclideanSpace {
    // TODO: maybe just have geodesic being used in solver. This doesn't need to be that generic here.
    fn apply(&self, t: f64, y: &OVector<f64, Const<8>>) -> OVector<f64, Const<8>> {
        self.geodesic(t, y)
    }
}

impl HasCoordinateSystem for EuclideanSpace {
    fn coordinate_system(&self) -> CoordinateSystem {
        Cartesian
    }
}

impl Geometry for EuclideanSpace {
    fn geodesic(&self, _: f64, y: &EquationOfMotionState) -> EquationOfMotionState {
        let y_new =
            EquationOfMotionState::from_column_slice(&[y[4], y[5], y[6], y[7], 0.0, 0.0, 0.0, 0.0]);
        y_new
    }

    // TODO: take into account rotations.
    fn get_tetrad_at(&self, position: &Vector4<f64>) -> Tetrad {
        Tetrad::new(
            position.clone(),
            FourVector::new_cartesian(1.0, 0.0, 0.0, 0.0),
            FourVector::new_cartesian(0.0, 1.0, 0.0, 0.0),
            FourVector::new_cartesian(0.0, 0.0, 1.0, 0.0),
            FourVector::new_cartesian(0.0, 0.0, 0.0, 1.0),
        )
    }

    fn lorentz_transformation(
        &self,
        position: &Vector4<f64>,
        t_velocity: &FourVector,
    ) -> Matrix4<f64> {
        let t_tetrad = self.get_tetrad_at(position).t.get_as_vector();

        let velocity = t_velocity.get_as_vector(); // TODO: add indexing to FourVector.
        let gamma = t_tetrad[0] * velocity[0]
            - t_tetrad[1] * velocity[1]
            - t_tetrad[2] * velocity[2]
            - t_tetrad[3] * velocity[3];

        let mut matrix = Matrix4::zeros();

        for i in 0..4 {
            for j in 0..4 {
                let mut res = 0.0;
                if i == j {
                    res = 1.0;
                }

                let g;
                if j == 0 {
                    g = 1.0;
                } else {
                    g = -1.0;
                }

                let a = -1.0 / (1.0 + gamma);
                let b = t_tetrad[i] + velocity[i];
                let c = g * (t_tetrad[j] + velocity[j]);
                res += a * b * c;

                res += 2.0 * g * t_tetrad[i] * velocity[j];

                matrix[(i, j)] = res;
            }
        }
        matrix
    }

    fn mul(&self, _position: &Vector4<f64>, v: &FourVector, w: &FourVector) -> f64 {
        1.0 * v.vector[0] * w.vector[0]
            + (-1.0) * v.vector[1] * w.vector[1]
            + (-1.0) * v.vector[2] * w.vector[2]
            + (-1.0) * v.vector[3] * w.vector[3]
    }

    fn get_stationary_velocity_at(&self, _position: &Vector4<f64>) -> FourVector {
        FourVector::new_cartesian(1.0, 0.0, 0.0, 0.0)
    }
}

#[cfg(test)]
mod tests {
    use crate::camera::Camera;
    use crate::debug::save_rays_to_file;
    use crate::euclidean::EuclideanSpace;
    use crate::four_vector::FourVector;
    use nalgebra::Vector4;
    use std::f64::consts::PI;

    #[ignore]
    #[test]
    fn save_camera_rays() {
        let rows = 30;
        let cols = 30;

        let position = Vector4::new(0.0, 0.0, 0.0, -10.0);
        let geometry = EuclideanSpace::new();
        let velocity = FourVector::new_cartesian(1.0, 0.0, 0.0, 0.0);
        let camera = Camera::new(position, velocity, PI / 2.0, rows, cols, geometry.clone());
        save_rays_to_file(rows, cols, &position, geometry, camera);
    }
}

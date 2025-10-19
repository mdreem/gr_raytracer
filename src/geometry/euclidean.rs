use crate::geometry::four_vector::FourVector;
use crate::geometry::geometry::{
    GeodesicSolver, Geometry, HasCoordinateSystem, InnerProduct, Tetrad,
};
use crate::geometry::point::CoordinateSystem::Cartesian;
use crate::geometry::point::{CoordinateSystem, Point};
use crate::rendering::runge_kutta::OdeFunction;
use crate::rendering::scene::EquationOfMotionState;
use nalgebra::{Const, Matrix4, OVector};

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
impl GeodesicSolver for EuclideanSpace {
    fn geodesic(&self, _: f64, y: &EquationOfMotionState) -> EquationOfMotionState {
        EquationOfMotionState::from_column_slice(&[y[4], y[5], y[6], y[7], 0.0, 0.0, 0.0, 0.0])
    }
}

impl HasCoordinateSystem for EuclideanSpace {
    fn coordinate_system(&self) -> CoordinateSystem {
        Cartesian
    }
}

impl InnerProduct for EuclideanSpace {
    fn inner_product(&self, _position: &Point, v: &FourVector, w: &FourVector) -> f64 {
        1.0 * v.vector[0] * w.vector[0]
            + (-1.0) * v.vector[1] * w.vector[1]
            + (-1.0) * v.vector[2] * w.vector[2]
            + (-1.0) * v.vector[3] * w.vector[3]
    }
}

impl Geometry for EuclideanSpace {
    /// Returns a local tetrad (orthonormal basis) at the given position.
    ///
    /// The tetrad is constructed from the spherical coordinate unit vectors as follows:
    /// - x-axis: eφ (unit vector in the increasing φ direction)
    /// - y-axis: -eθ (negative unit vector in the θ direction)
    /// - z-axis: -er (negative unit vector in the radial direction)
    ///
    /// This mapping (x = eφ, y = −eθ, z = −er) should match the one used in the Schwarzschild geometry.
    /// The time-like basis vector is e_t.
    fn get_tetrad_at(&self, position: &Point) -> Tetrad {
        let point_in_spherical = position.get_as_spherical();

        let _r = point_in_spherical.x;
        let theta = point_in_spherical.y;
        let phi = point_in_spherical.z;

        let e_t = FourVector::new_cartesian(1.0, 0.0, 0.0, 0.0);
        let e_r = FourVector::new_cartesian(
            0.0,
            theta.sin() * phi.cos(),
            theta.sin() * phi.sin(),
            theta.cos(),
        );
        let e_theta = FourVector::new_cartesian(
            0.0,
            theta.cos() * phi.cos(),
            theta.cos() * phi.sin(),
            -theta.sin(),
        );
        let e_phi = FourVector::new_cartesian(0.0, -phi.sin(), phi.cos(), 0.0);

        Tetrad::new(position.clone(), e_t, e_phi, -e_theta, -e_r)
    }

    fn lorentz_transformation(&self, position: &Point, t_velocity: &FourVector) -> Matrix4<f64> {
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

                let g = if j == 0 { 1.0 } else { -1.0 };

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

    fn get_stationary_velocity_at(&self, _position: &Point) -> FourVector {
        FourVector::new_cartesian(1.0, 0.0, 0.0, 0.0)
    }

    fn inside_horizon(&self, _position: &Point) -> bool {
        false
    }
}

#[cfg(test)]
mod tests {
    use crate::geometry::euclidean::EuclideanSpace;
    use crate::geometry::four_vector::FourVector;
    use crate::geometry::point::Point;
    use crate::rendering::camera::Camera;
    use crate::rendering::debug::save_rays_to_file;
    use std::f64::consts::PI;

    #[ignore]
    #[test]
    fn save_camera_rays() {
        let rows = 30;
        let cols = 30;

        let position = Point::new_cartesian(0.0, 0.0, 0.0, -10.0);
        let geometry = EuclideanSpace::new();
        let velocity = FourVector::new_cartesian(1.0, 0.0, 0.0, 0.0);
        let camera = Camera::new(
            position,
            velocity,
            PI / 2.0,
            rows,
            cols,
            0.0,
            0.0,
            0.0,
            &geometry,
        );
        save_rays_to_file(rows, cols, &position, geometry, camera);
    }
}

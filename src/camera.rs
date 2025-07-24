use crate::four_vector::FourVector;
use crate::geometry::{Geometry, Tetrad};
use nalgebra::Vector4;

#[derive(Debug)]
pub struct Ray {
    pub position: Vector4<f64>,
    pub momentum: FourVector,
    pub row: i64,
    pub col: i64,
}
impl Ray {
    pub fn new(row: i64, col: i64, position: Vector4<f64>, momentum: FourVector) -> Self {
        Self {
            row,
            col,
            position,
            momentum,
        }
    }
}

#[derive(Debug, Clone)]
pub struct Camera {
    alpha: f64,
    pub rows: i64,
    pub columns: i64,
    position: Vector4<f64>,
    pub velocity: FourVector,
    tetrad: Tetrad,
}

pub fn lorentz_transform_tetrad<G: Geometry>(
    geometry: &G,
    tetrad: &Tetrad,
    position: &Vector4<f64>,
    velocity: &FourVector,
) -> Tetrad {
    let lorentz = geometry.lorentz_transformation(position, velocity);

    println!("lorentz transformation: {:?}", lorentz);

    let t_vec = lorentz * tetrad.t.get_as_vector();
    let x_vec = lorentz * tetrad.x.get_as_vector();
    let y_vec = lorentz * tetrad.y.get_as_vector();
    let z_vec = lorentz * tetrad.z.get_as_vector();

    Tetrad::new(
        *position,
        FourVector::new(
            t_vec[0],
            t_vec[1],
            t_vec[2],
            t_vec[3],
            tetrad.t.coordinate_system,
        ),
        FourVector::new(
            x_vec[0],
            x_vec[1],
            x_vec[2],
            x_vec[3],
            tetrad.x.coordinate_system,
        ),
        FourVector::new(
            y_vec[0],
            y_vec[1],
            y_vec[2],
            y_vec[3],
            tetrad.y.coordinate_system,
        ),
        FourVector::new(
            z_vec[0],
            z_vec[1],
            z_vec[2],
            z_vec[3],
            tetrad.z.coordinate_system,
        ),
    )
}

impl Camera {
    // Position is given in cartesian coordinates.
    pub fn new<G: Geometry>(
        position: Vector4<f64>,
        velocity: FourVector,
        alpha: f64,
        rows: i64,
        columns: i64,
        geometry: &G,
    ) -> Camera {
        let original_tetrad = geometry.get_tetrad_at(&position);
        println!("original_tetrad: {:?}", original_tetrad);
        let tetrad = lorentz_transform_tetrad(geometry, &original_tetrad, &position, &velocity);
        println!("tetrad: {:?}", tetrad);
        Self {
            position,
            velocity,
            alpha,
            rows,
            columns,
            tetrad,
        }
    }

    // row, column range from 1..R, 1..C
    fn get_direction_for(&self, row: i64, column: i64) -> FourVector {
        let i_prime = (2.0 * f64::tan(self.alpha / 2.0) / (self.columns as f64))
            * (column as f64 - (self.columns as f64 + 1.0) / 2.0);
        let j_prime = (2.0 * f64::tan(self.alpha / 2.0) / (self.rows as f64))
            * (row as f64 - (self.rows as f64 + 1.0) / 2.0);

        let w = self.tetrad.z + i_prime * self.tetrad.x + j_prime * self.tetrad.y;
        let w_squared = -1.0 - i_prime * i_prime - j_prime * j_prime;

        -self.tetrad.z + 2.0 * w / (-w_squared)
    }

    // row, column range from 1..R, 1..C
    pub fn get_ray_for(&self, row: i64, column: i64) -> Ray {
        let direction = self.get_direction_for(row, column);
        let momentum = direction + self.tetrad.t; // Add T-component of the tetrad to get the momentum.
        Ray::new(row, column, self.position, momentum)
    }
}

#[cfg(test)]
mod tests {
    use crate::camera::Camera;
    use crate::euclidean::EuclideanSpace;

    use crate::four_vector::FourVector;
    use crate::geometry::{Geometry, InnerProduct};
    use approx::assert_abs_diff_eq;
    use nalgebra::Vector4;
    use std::f64::consts::PI;

    #[test]
    fn test_get_direction_for() {
        let camera = Camera::new(
            Vector4::new(0.0, 0.0, 1.0, 0.0),
            FourVector::new_cartesian(1.0, 0.0, 0.0, 0.0),
            PI / 2.0,
            11,
            11,
            &EuclideanSpace::new(),
        );
        let geometry = EuclideanSpace::new();
        let position = Vector4::new(0.0, 0.0, 1.0, 0.0);

        let top_left_corner = camera.get_direction_for(1, 1);
        let top_right_corner = camera.get_direction_for(1, 11);
        let middle = camera.get_direction_for(6, 6);
        let bottom_left_corner = camera.get_direction_for(11, 1);
        let bottom_right_corner = camera.get_direction_for(11, 11);

        let corner = -0.6853582554517135;
        let corner_z = -0.24610591900311507;
        assert_abs_diff_eq!(
            top_left_corner.get_as_vector(),
            Vector4::new(0.0, corner, corner, corner_z)
        );
        let top_left_corner_scalar =
            geometry.inner_product(&position, &top_left_corner, &top_left_corner);
        assert_abs_diff_eq!(top_left_corner_scalar, -1.0);

        assert_abs_diff_eq!(
            top_right_corner.get_as_vector(),
            Vector4::new(0.0, -corner, corner, corner_z)
        );
        let top_right_corner_scalar =
            geometry.inner_product(&position, &top_right_corner, &top_right_corner);
        assert_abs_diff_eq!(top_right_corner_scalar, -1.0);

        assert_abs_diff_eq!(middle.get_as_vector(), Vector4::new(0.0, 0.0, 0.0, 1.0));
        let middle_scalar = geometry.inner_product(&position, &middle, &middle);
        assert_abs_diff_eq!(middle_scalar, -1.0);

        assert_abs_diff_eq!(
            bottom_left_corner.get_as_vector(),
            Vector4::new(0.0, corner, -corner, corner_z)
        );
        let bottom_left_corner_scalar =
            geometry.inner_product(&position, &bottom_left_corner, &bottom_left_corner);
        assert_abs_diff_eq!(bottom_left_corner_scalar, -1.0);

        assert_abs_diff_eq!(
            bottom_right_corner.get_as_vector(),
            Vector4::new(0.0, -corner, -corner, corner_z)
        );
        let bottom_right_corner_scalar =
            geometry.inner_product(&position, &bottom_right_corner, &bottom_right_corner);
        assert_abs_diff_eq!(bottom_right_corner_scalar, -1.0);
    }
}

use crate::four_vector::FourVector;
use crate::geometry::{Geometry, Tetrad};
use nalgebra::Vector4;

#[derive(Debug)]
pub struct Ray {
    pub position: Vector4<f64>,
    pub direction: FourVector,
}
impl Ray {
    pub fn new(position: Vector4<f64>, direction: FourVector) -> Self {
        Self {
            position,
            direction,
        }
    }
}

pub struct Camera<G: Geometry> {
    alpha: f64,
    rows: i64,
    columns: i64,
    position: Vector4<f64>,
    geometry: G,
    tetrad: Tetrad,
}

impl<G: Geometry> Camera<G> {
    pub fn new(
        position: Vector4<f64>,
        alpha: f64,
        rows: i64,
        columns: i64,
        geometry: G,
    ) -> Camera<G> {
        let tetrad = geometry.get_tetrad_at(&position);
        Self {
            position,
            alpha,
            rows,
            columns,
            geometry,
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
        Ray::new(self.position, direction)
    }
}

#[cfg(test)]
mod tests {
    use crate::camera::Camera;
    use crate::geometry::EuclideanSpace;
    use approx::assert_abs_diff_eq;
    use nalgebra::Vector4;

    #[test]
    fn test_get_direction_for() {
        let camera = Camera::new(
            Vector4::new(0.0, 0.0, 1.0, 0.0),
            std::f64::consts::PI / 2.0,
            11,
            11,
            EuclideanSpace::new(),
        );

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
        assert_abs_diff_eq!(top_left_corner * top_left_corner, -1.0);

        assert_abs_diff_eq!(
            top_right_corner.get_as_vector(),
            Vector4::new(0.0, -corner, corner, corner_z)
        );
        assert_abs_diff_eq!(top_right_corner * top_right_corner, -1.0);

        assert_abs_diff_eq!(middle.get_as_vector(), Vector4::new(0.0, 0.0, 0.0, 1.0));
        assert_abs_diff_eq!(middle * middle, -1.0);

        assert_abs_diff_eq!(
            bottom_left_corner.get_as_vector(),
            Vector4::new(0.0, corner, -corner, corner_z)
        );
        assert_abs_diff_eq!(bottom_left_corner * bottom_left_corner, -1.0);

        assert_abs_diff_eq!(
            bottom_right_corner.get_as_vector(),
            Vector4::new(0.0, -corner, -corner, corner_z)
        );
        assert_abs_diff_eq!(bottom_right_corner * bottom_right_corner, -1.0);
    }
}

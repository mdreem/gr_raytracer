use crate::geometry::four_vector::FourVector;
use crate::geometry::geometry::{Geometry, Tetrad};
use crate::geometry::point::Point;
use crate::rendering::ray::Ray;
use log::debug;

#[derive(Debug, Clone)]
pub struct Camera {
    alpha: f64,
    pub rows: i64,
    pub columns: i64,
    position: Point,
    pub velocity: FourVector,
    tetrad: Tetrad,
}

pub fn lorentz_transform_tetrad<G: Geometry>(
    geometry: &G,
    tetrad: &Tetrad,
    position: &Point,
    velocity: &FourVector,
) -> Tetrad {
    let lorentz = geometry.lorentz_transformation(position, velocity);

    debug!("lorentz transformation: {:?}", lorentz);

    let t_vec = lorentz * tetrad.t.get_as_vector();
    let x_vec = lorentz * tetrad.x.get_as_vector();
    let y_vec = lorentz * tetrad.y.get_as_vector();
    let z_vec = lorentz * tetrad.z.get_as_vector();

    Tetrad::new(
        position.clone(),
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
        position: Point,
        velocity: FourVector,
        alpha: f64,
        rows: i64,
        columns: i64,
        geometry: &G,
    ) -> Camera {
        let original_tetrad = geometry.get_tetrad_at(&position);
        debug!("position: {:?}", position);
        debug!("original_tetrad: {}", original_tetrad);
        debug!("inner product checks:");
        debug!(
            "  inner product t.t: {}",
            geometry.inner_product(&position, &original_tetrad.t, &original_tetrad.t)
        );
        debug!(
            "  inner product x.x: {}",
            geometry.inner_product(&position, &original_tetrad.x, &original_tetrad.x)
        );
        debug!(
            "  inner product y.y: {}",
            geometry.inner_product(&position, &original_tetrad.y, &original_tetrad.y)
        );
        debug!(
            "  inner product z.z: {}",
            geometry.inner_product(&position, &original_tetrad.z, &original_tetrad.z)
        );
        debug!("");
        debug!(
            "  inner product t.x: {}",
            geometry.inner_product(&position, &original_tetrad.t, &original_tetrad.x)
        );
        debug!(
            "  inner product t.y: {}",
            geometry.inner_product(&position, &original_tetrad.t, &original_tetrad.y)
        );
        debug!(
            "  inner product t.z: {}",
            geometry.inner_product(&position, &original_tetrad.t, &original_tetrad.z)
        );
        debug!("");
        debug!(
            "  inner product x.y: {}",
            geometry.inner_product(&position, &original_tetrad.x, &original_tetrad.y)
        );
        debug!(
            "  inner product x.z: {}",
            geometry.inner_product(&position, &original_tetrad.x, &original_tetrad.z)
        );
        debug!("");
        debug!(
            "  inner product y.z: {}",
            geometry.inner_product(&position, &original_tetrad.y, &original_tetrad.z)
        );
        let tetrad = lorentz_transform_tetrad(geometry, &original_tetrad, &position, &velocity);
        debug!("tetrad: {}", tetrad);
        Self {
            position,
            velocity,
            alpha,
            rows,
            columns,
            tetrad,
        }
    }

    // row, column range from 1..R, 1..C in https://arxiv.org/abs/1511.06025, but here we work
    // with 0-based indices. This needs to be accounted for below.
    fn get_direction_for(&self, row: i64, column: i64) -> FourVector {
        let shifted_column = (column + 1) as f64; // Convert to 1-based index.
        let shifted_row = (row + 1) as f64; // Convert to 1-based index.
        let i_prime = (2.0 * f64::tan(self.alpha / 2.0) / (self.columns as f64))
            * (shifted_column - (self.columns as f64 + 1.0) / 2.0);
        let j_prime = (2.0 * f64::tan(self.alpha / 2.0) / (self.rows as f64))
            * (shifted_row as f64 - (self.rows as f64 + 1.0) / 2.0);

        let w = self.tetrad.z + i_prime * self.tetrad.x + j_prime * self.tetrad.y;
        let w_squared = -1.0 - i_prime * i_prime - j_prime * j_prime;

        -self.tetrad.z + 2.0 * w / (-w_squared)
    }

    pub fn get_ray_for(&self, row: i64, column: i64) -> Ray {
        let direction = self.get_direction_for(row, column);
        let momentum = direction + self.tetrad.t; // Add T-component of the tetrad to get the momentum.
        Ray::new(
            row,
            column,
            self.rows,
            self.columns,
            self.position,
            momentum,
        )
    }
}

#[cfg(test)]
mod tests {
    use crate::geometry::euclidean::EuclideanSpace;
    use crate::geometry::four_vector::FourVector;
    use crate::geometry::geometry::InnerProduct;
    use crate::geometry::point::{CoordinateSystem, Point};
    use crate::rendering::camera::Camera;
    use approx::assert_abs_diff_eq;
    use std::f64::consts::PI;

    #[test]
    fn test_get_direction_for() {
        let position = Point::new(0.0, 0.0, 0.0, 1.0, CoordinateSystem::Cartesian);
        let camera = Camera::new(
            position.clone(),
            FourVector::new_cartesian(1.0, 0.0, 0.0, 0.0),
            PI / 2.0,
            11,
            11,
            &EuclideanSpace::new(),
        );
        let geometry = EuclideanSpace::new();

        let top_left_corner = camera.get_direction_for(0, 0);
        let top_right_corner = camera.get_direction_for(0, 10);
        let middle = camera.get_direction_for(5, 5);
        let bottom_left_corner = camera.get_direction_for(10, 0);
        let bottom_right_corner = camera.get_direction_for(10, 10);

        let corner = -0.6853582554517135;
        let corner_z = 0.24610591900311507;
        assert_abs_diff_eq!(
            top_left_corner.get_as_vector(),
            FourVector::new_cartesian(0.0, corner_z, corner, corner).get_as_vector()
        );
        let top_left_corner_scalar =
            geometry.inner_product(&position, &top_left_corner, &top_left_corner);
        assert_abs_diff_eq!(top_left_corner_scalar, -1.0);

        assert_abs_diff_eq!(
            top_right_corner.get_as_vector(),
            FourVector::new_cartesian(0.0, corner_z, -corner, corner).get_as_vector()
        );
        let top_right_corner_scalar =
            geometry.inner_product(&position, &top_right_corner, &top_right_corner);
        assert_abs_diff_eq!(top_right_corner_scalar, -1.0);

        assert_abs_diff_eq!(
            middle.get_as_vector(),
            FourVector::new_cartesian(0.0, -1.0, 0.0, 0.0).get_as_vector()
        );
        let middle_scalar = geometry.inner_product(&position, &middle, &middle);
        assert_abs_diff_eq!(middle_scalar, -1.0);

        assert_abs_diff_eq!(
            bottom_left_corner.get_as_vector(),
            FourVector::new_cartesian(0.0, corner_z, corner, -corner).get_as_vector()
        );
        let bottom_left_corner_scalar =
            geometry.inner_product(&position, &bottom_left_corner, &bottom_left_corner);
        assert_abs_diff_eq!(bottom_left_corner_scalar, -1.0);

        assert_abs_diff_eq!(
            bottom_right_corner.get_as_vector(),
            FourVector::new_cartesian(0.0, corner_z, -corner, -corner).get_as_vector()
        );
        let bottom_right_corner_scalar =
            geometry.inner_product(&position, &bottom_right_corner, &bottom_right_corner);
        assert_abs_diff_eq!(bottom_right_corner_scalar, -1.0);
    }
}

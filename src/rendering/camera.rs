use crate::geometry::four_vector::FourVector;
use crate::geometry::geometry::Geometry;
use crate::geometry::point::Point;
use crate::geometry::tetrad::{Tetrad, TetradValidator};
use crate::rendering::ray::Ray;
use log::{debug, trace};

#[derive(Debug, PartialEq, thiserror::Error)]
pub enum CameraError {
    #[error("Tetrad is not orthonormal")]
    TetradNotOrthonormal,
}

#[derive(Debug, Clone)]
pub struct Camera {
    alpha: f64,
    pub rows: i64,
    pub columns: i64,
    pub position: Point,
    pub velocity: FourVector,
    tetrad: Tetrad,
    spatial_signature: f64,
}

pub fn lorentz_transform_tetrad<G: Geometry>(
    geometry: &G,
    tetrad: &Tetrad,
    position: &Point,
    velocity: &FourVector,
) -> Tetrad {
    let lorentz = geometry.lorentz_transformation(position, velocity);

    debug!("lorentz transformation: {:?}", lorentz);

    // TODO: Move matrix multiplication to FourVector.
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

fn rotate(v1: FourVector, v2: FourVector, angle: f64) -> (FourVector, FourVector) {
    let v1_rotated = angle.cos() * v1 + angle.sin() * v2;
    let v2_rotated = -angle.sin() * v1 + angle.cos() * v2;

    (v1_rotated, v2_rotated)
}

impl Camera {
    pub fn new<G: Geometry>(
        position: Point,
        velocity: FourVector,
        alpha: f64,
        rows: i64,
        columns: i64,
        phi: f64,
        theta: f64,
        psi: f64,
        geometry: &G,
    ) -> Result<Camera, CameraError> {
        let original_tetrad = geometry.get_tetrad_at(&position);

        let tetrad_validator = TetradValidator::new(geometry);
        tetrad_validator.validate(&original_tetrad)?;

        let (a_prime, b_prime) = rotate(original_tetrad.x, original_tetrad.y, phi);
        let (z, a_two_prime) = rotate(original_tetrad.z, a_prime, theta);
        let (x, y) = rotate(a_two_prime, b_prime, psi);

        let rotated_tetrad = Tetrad::new(position.clone(), original_tetrad.t, x, y, z);
        debug!("rotated tetrad: {}", rotated_tetrad);

        let tetrad = lorentz_transform_tetrad(geometry, &rotated_tetrad, &position, &velocity);
        debug!("lorentz transformed tetrad: {}", tetrad);
        tetrad_validator.validate(&tetrad)?;

        let signature = geometry.signature();
        debug_assert!(
            (signature[1] - signature[2]).abs() < 1e-12
                && (signature[2] - signature[3]).abs() < 1e-12,
            "Camera projection expects equal spatial signature components"
        );

        Ok(Self {
            position,
            velocity,
            alpha,
            rows,
            columns,
            tetrad,
            spatial_signature: signature[3],
        })
    }

    /// Generates ray direction using pinhole camera projection.
    ///
    /// Algorithm (from https://arxiv.org/abs/1511.06025):
    /// 1. Map pixel (row, col) to image plane coordinates (i', j')
    /// 2. Construct spacelike vector: w = e_z + i' e_x + j' e_y
    /// 3. Project to null direction: k = -e_z + 2w / (-w·w)
    ///
    /// This ensures k·k = 0 (null geodesic) in the observer's tetrad frame.
    ///
    /// Coordinate convention:
    /// - e_x: horizontal (right) in image
    /// - e_y: vertical (up) in image
    /// - e_z: camera forward direction
    ///
    /// row, column range from 1..R, 1..C in https://arxiv.org/abs/1511.06025, but here we work
    /// with 0-based indices. This needs to be accounted for below.
    fn get_direction_for(&self, row: i64, column: i64) -> FourVector {
        let shifted_column = (column + 1) as f64; // Convert to 1-based index.
        let shifted_row = (row + 1) as f64; // Convert to 1-based index.
        let tan_half_alpha = f64::tan(self.alpha / 2.0);
        // To ensure square pixels, we use the same angular scale for both directions.
        // alpha is treated as the vertical field of view (FOV).
        let i_prime = (2.0 * tan_half_alpha / (self.rows as f64))
            * (shifted_column - (self.columns as f64 + 1.0) / 2.0);
        let j_prime = (2.0 * tan_half_alpha / (self.rows as f64))
            * (shifted_row as f64 - (self.rows as f64 + 1.0) / 2.0);

        let w = self.tetrad.z + i_prime * self.tetrad.x + j_prime * self.tetrad.y;
        // Keep the paper-style notation w_squared = w·w. Its sign depends on the
        // metric signature, so normalize with spatial_signature to keep the same formula.
        let w_squared = self.spatial_signature * (1.0 + i_prime * i_prime + j_prime * j_prime);

        -self.tetrad.z + 2.0 * w / (self.spatial_signature * w_squared)
    }

    pub fn get_ray_for(&self, row: i64, column: i64) -> Ray {
        let direction = self.get_direction_for(row, column);
        trace!("direction ({}|{}): {:?}", row, column, direction);
        let momentum = direction + self.tetrad.t; // Add T-component of the tetrad to get the momentum.
        Ray::new(row, column, self.position, momentum)
    }
}

#[cfg(test)]
mod tests {
    use crate::geometry::euclidean::EuclideanSpace;
    use crate::geometry::euclidean_spherical::EuclideanSpaceSpherical;
    use crate::geometry::four_vector::FourVector;
    use crate::geometry::geometry::InnerProduct;
    use crate::geometry::point::{CoordinateSystem, Point};
    use crate::geometry::schwarzschild::Schwarzschild;
    use crate::geometry::spherical_coordinates_helper::{
        cartesian_to_spherical, spherical_to_cartesian,
    };
    use crate::rendering::camera::Camera;
    use approx::assert_abs_diff_eq;
    use std::f64::consts::PI;

    #[test]
    fn test_get_direction_for() {
        let position = Point::new(0.0, 1.0, 0.0, 0.0, CoordinateSystem::Cartesian);
        let camera = Camera::new(
            position.clone(),
            FourVector::new_cartesian(1.0, 0.0, 0.0, 0.0),
            PI / 2.0,
            11,
            11,
            0.0,
            0.0,
            0.0,
            &EuclideanSpace::new(),
        )
        .unwrap();
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

    #[test]
    fn test_get_ray_for_different_geometries_euclidean() {
        let position = Point::new_cartesian(0.0, 0.0, 1.0, 0.0);
        let velocity = FourVector::new_cartesian(1.0, 0.0, 0.0, 0.0);

        let geometry_euclidean = EuclideanSpace::new();
        let camera_euclidean = Camera::new(
            position.clone(),
            velocity.clone(),
            PI / 2.0,
            100,
            100,
            0.0,
            PI / 2.0,
            PI / 2.0,
            &geometry_euclidean,
        )
        .unwrap();

        let geometry_euclidean_spherical = EuclideanSpaceSpherical::new();
        let camera_euclidean_spherical = Camera::new(
            cartesian_to_spherical(&position.clone()),
            velocity.clone(),
            PI / 2.0,
            100,
            100,
            0.0,
            PI / 2.0,
            PI / 2.0,
            &geometry_euclidean_spherical,
        )
        .unwrap();

        for _row in 0..100 {
            for _col in 0..100 {
                let ray_euclidean = camera_euclidean.get_ray_for(_row, _col);
                let ray_euclidean_spherical = camera_euclidean_spherical.get_ray_for(_row, _col);
                assert_abs_diff_eq!(
                    ray_euclidean.position.vector,
                    spherical_to_cartesian(&ray_euclidean_spherical.position).vector,
                    epsilon = 1e-10
                );
            }
        }
    }

    #[test]
    fn test_get_ray_for_different_geometries_schwarzschild() {
        let position = Point::new_cartesian(0.0, 0.0, 1.0, 0.0);
        let velocity = FourVector::new_cartesian(1.0, 0.0, 0.0, 0.0);

        let geometry_euclidean = EuclideanSpace::new();
        let camera_euclidean = Camera::new(
            position.clone(),
            velocity.clone(),
            PI / 2.0,
            100,
            100,
            0.0,
            PI / 2.0,
            PI / 2.0,
            &geometry_euclidean,
        )
        .unwrap();

        let geometry_euclidean_schwarzschild = Schwarzschild::new(0.0, 0.0);
        let camera_euclidean_schwarzschild = Camera::new(
            cartesian_to_spherical(&position.clone()),
            velocity.clone(),
            PI / 2.0,
            100,
            100,
            0.0,
            PI / 2.0,
            PI / 2.0,
            &geometry_euclidean_schwarzschild,
        )
        .unwrap();

        for _row in 0..100 {
            for _col in 0..100 {
                let ray_euclidean = camera_euclidean.get_ray_for(_row, _col);
                let ray_euclidean_schwarzschild =
                    camera_euclidean_schwarzschild.get_ray_for(_row, _col);
                assert_abs_diff_eq!(
                    ray_euclidean.position.vector,
                    spherical_to_cartesian(&ray_euclidean_schwarzschild.position).vector,
                    epsilon = 1e-10
                );
            }
        }
    }
}

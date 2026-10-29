pub mod debug;
use crate::geometry::coordinates::conversion::CartesianConvertible;
use crate::geometry::coordinates::chart::Chart;
use crate::geometry::four_vector::FourVector;
use crate::geometry::geometry::Geometry;
use crate::geometry::point::Point;
use crate::geometry::tetrad::{Tetrad, TetradValidator};
use crate::rendering::ray::Ray;
use log::{debug, trace};
use nalgebra::Vector3;

#[derive(Debug, PartialEq, thiserror::Error)]
pub enum CameraError {
    #[error("Tetrad is not orthonormal")]
    TetradNotOrthonormal,
}

#[derive(Debug, Clone)]
pub struct Camera<C: Chart> {
    alpha: f64,
    pub rows: i64,
    pub columns: i64,
    pub position: Point<C>,
    pub velocity: FourVector<C>,
    tetrad: Tetrad<C>,
    spatial_signature: f64,
    spatial_handedness: f64,
}

pub fn lorentz_transform_tetrad<G: Geometry<Chart = C>, C: Chart>(
    geometry: &G,
    tetrad: &Tetrad<C>,
    position: &Point<C>,
    velocity: &FourVector<C>,
) -> Tetrad<C> {
    let lorentz = geometry.lorentz_transformation(position, velocity);

    debug!("lorentz transformation: {:?}", lorentz);

    // TODO: Move matrix multiplication to FourVector.
    let t_vec = lorentz * tetrad.t.vector;
    let x_vec = lorentz * tetrad.x.vector;
    let y_vec = lorentz * tetrad.y.vector;
    let z_vec = lorentz * tetrad.z.vector;

    Tetrad::new(
        *position,
        FourVector::new(t_vec[0], t_vec[1], t_vec[2], t_vec[3], tetrad.t.at),
        FourVector::new(x_vec[0], x_vec[1], x_vec[2], x_vec[3], tetrad.x.at),
        FourVector::new(y_vec[0], y_vec[1], y_vec[2], y_vec[3], tetrad.y.at),
        FourVector::new(z_vec[0], z_vec[1], z_vec[2], z_vec[3], tetrad.z.at),
    )
}

// TODO: This likely assumes cartesian vectors. Confirm.
fn rotate<C: Chart>(
    v1: FourVector<C>,
    v2: FourVector<C>,
    angle: f64,
) -> (FourVector<C>, FourVector<C>) {
    let v1_rotated = angle.cos() * v1 + angle.sin() * v2;
    let v2_rotated = -angle.sin() * v1 + angle.cos() * v2;

    (v1_rotated, v2_rotated)
}

fn spatial_basis_vector_cartesian<C: CartesianConvertible>(
    _position: &Point<C>,
    vector: &FourVector<C>,
) -> Vector3<f64> {
    let v = vector.to_cartesian().vector;
    Vector3::new(v[1], v[2], v[3])
}

fn spatial_handedness<C: CartesianConvertible>(position: &Point<C>, tetrad: &Tetrad<C>) -> f64 {
    let x = spatial_basis_vector_cartesian(position, &tetrad.x);
    let y = spatial_basis_vector_cartesian(position, &tetrad.y);
    let z = spatial_basis_vector_cartesian(position, &tetrad.z);

    let triple_product = x.dot(&y.cross(&z));
    if !triple_product.is_finite() || triple_product.abs() <= 1e-12 {
        return tetrad.x.at.chart.tetrad_handedness_fallback();
    }
    if triple_product >= 0.0 { 1.0 } else { -1.0 }
}

impl<C: CartesianConvertible> Camera<C> {
    pub fn new<G: Geometry<Chart = C>>(
        position: Point<C>,
        velocity: FourVector<C>,
        alpha: f64,
        rows: i64,
        columns: i64,
        phi: f64,
        theta: f64,
        psi: f64,
        geometry: &G,
    ) -> Result<Camera<C>, CameraError> {
        let original_tetrad = geometry.get_tetrad_at(&position);

        let tetrad_validator = TetradValidator::new(geometry);
        tetrad_validator.validate(&original_tetrad)?;

        let (a_prime, b_prime) = rotate(original_tetrad.x, original_tetrad.y, phi);
        let (z, a_two_prime) = rotate(original_tetrad.z, a_prime, theta);
        let (x, y) = rotate(a_two_prime, b_prime, psi);

        let rotated_tetrad = Tetrad::new(position, original_tetrad.t, x, y, z);
        debug!("rotated tetrad: {}", rotated_tetrad);

        let tetrad = lorentz_transform_tetrad(geometry, &rotated_tetrad, &position, &velocity);
        debug!("lorentz transformed tetrad: {}", tetrad);
        tetrad_validator.validate(&tetrad)?;
        let spatial_handedness = spatial_handedness(&position, &tetrad);

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
            spatial_handedness,
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
    fn get_direction_for(&self, row: f64, column: f64) -> FourVector<C> {
        let shifted_column = column + 1.0; // Convert to 1-based index.
        let shifted_row = row + 1.0; // Convert to 1-based index.
        let tan_half_alpha = f64::tan(self.alpha / 2.0);
        // To ensure square pixels, we use the same angular scale for both directions.
        // alpha is treated as the vertical field of view (FOV).
        let i_prime = self.spatial_handedness
            * (2.0 * tan_half_alpha / (self.rows as f64))
            * (shifted_column - (self.columns as f64 + 1.0) / 2.0);
        let j_prime = (2.0 * tan_half_alpha / (self.rows as f64))
            * (shifted_row - (self.rows as f64 + 1.0) / 2.0);

        let w = self.tetrad.z + i_prime * self.tetrad.x + j_prime * self.tetrad.y;
        // Keep the paper-style notation w_squared = w·w. Its sign depends on the
        // metric signature, so normalize with spatial_signature to keep the same formula.
        let w_squared = self.spatial_signature * (1.0 + i_prime * i_prime + j_prime * j_prime);

        -self.tetrad.z + 2.0 * w / (self.spatial_signature * w_squared)
    }

    /// Angle one pixel spans at the image centre (alpha is the vertical field
    /// of view). An approximation away from the centre, good enough as a
    /// "different pixel" threshold.
    pub fn pixel_angle(&self) -> f64 {
        self.alpha / self.rows as f64
    }

    pub fn get_ray_for(&self, row: i64, column: i64) -> Ray<C> {
        let direction = self.get_direction_for(row as f64, column as f64);
        trace!("direction ({}|{}): {:?}", row, column, direction);
        // Trace the physical photon that arrives at the camera from
        // `direction`: its traced momentum N - e_t is past-directed, so the
        // integration marches backwards in time, from the camera out towards
        // the emitter. The future-directed outbound ray (N + e_t) is not the
        // same geodesic: it flips the Doppler term of the redshift for moving
        // emitters and follows a path with reversed frame dragging in Kerr.
        let momentum = direction + (-self.tetrad.t);
        Ray::new(row, column, self.position, momentum)
    }

    pub fn get_ray_for_offset(&self, row: i64, column: i64, dx: f64, dy: f64) -> Ray<C> {
        // The original pixel is in the middle of the sweeping interval from 0.0...1.0, so we need to subtract 0.5 to get the correct offset.
        let direction =
            self.get_direction_for((row as f64) + (dy - 0.5), (column as f64) + (dx - 0.5));
        trace!("direction ({}|{}): {:?}", row, column, direction);
        let momentum = direction + (-self.tetrad.t);
        Ray::new(row, column, self.position, momentum)
    }
}

#[cfg(test)]
mod tests {
    use crate::geometry::coordinates::conversion::CartesianConvertible;
    use crate::geometry::atlases::kerr::KerrChart;
    use crate::geometry::spacetimes::euclidean::EuclideanSpace;
    use crate::geometry::spacetimes::euclidean_spherical::EuclideanSpaceSpherical;
    use crate::geometry::four_vector::FourVector;
    use crate::geometry::geometry::{InnerProduct, Signature};
    use crate::geometry::spacetimes::kerr::Kerr;
    use crate::geometry::point::Point;
    use crate::geometry::spacetimes::schwarzschild::Schwarzschild;
    use crate::geometry::spherical_coordinates_helper::{
        cartesian_to_spherical, spherical_to_cartesian,
    };
    use crate::rendering::camera::Camera;
    use approx::assert_abs_diff_eq;
    use std::f64::consts::PI;

    #[test]
    fn test_get_direction_for() {
        let position = Point::new_cartesian(0.0, 1.0, 0.0, 0.0);
        let camera = Camera::new(
            position.clone(),
            FourVector::new(1.0, 0.0, 0.0, 0.0, position),
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

        let top_left_corner = camera.get_direction_for(0.0, 0.0);
        let top_right_corner = camera.get_direction_for(0.0, 10.0);
        let middle = camera.get_direction_for(5.0, 5.0);
        let bottom_left_corner = camera.get_direction_for(10.0, 0.0);
        let bottom_right_corner = camera.get_direction_for(10.0, 10.0);

        let corner = -0.6853582554517135;
        let corner_z = 0.24610591900311507;
        assert_abs_diff_eq!(
            top_left_corner.vector,
            FourVector::new(0.0, corner_z, -corner, corner, position).vector
        );
        let top_left_corner_scalar =
            geometry.inner_product(&position, &top_left_corner, &top_left_corner);
        assert_abs_diff_eq!(top_left_corner_scalar, -1.0);

        assert_abs_diff_eq!(
            top_right_corner.vector,
            FourVector::new(0.0, corner_z, corner, corner, position).vector
        );
        let top_right_corner_scalar =
            geometry.inner_product(&position, &top_right_corner, &top_right_corner);
        assert_abs_diff_eq!(top_right_corner_scalar, -1.0);

        assert_abs_diff_eq!(
            middle.vector,
            FourVector::new(0.0, -1.0, 0.0, 0.0, position).vector
        );
        let middle_scalar = geometry.inner_product(&position, &middle, &middle);
        assert_abs_diff_eq!(middle_scalar, -1.0);

        assert_abs_diff_eq!(
            bottom_left_corner.vector,
            FourVector::new(0.0, corner_z, -corner, -corner, position).vector
        );
        let bottom_left_corner_scalar =
            geometry.inner_product(&position, &bottom_left_corner, &bottom_left_corner);
        assert_abs_diff_eq!(bottom_left_corner_scalar, -1.0);

        assert_abs_diff_eq!(
            bottom_right_corner.vector,
            FourVector::new(0.0, corner_z, corner, -corner, position).vector
        );
        let bottom_right_corner_scalar =
            geometry.inner_product(&position, &bottom_right_corner, &bottom_right_corner);
        assert_abs_diff_eq!(bottom_right_corner_scalar, -1.0);
    }

    #[test]
    fn centered_offset_produces_the_base_pixel_ray() {
        let geometry = EuclideanSpace::new();
        let camera = Camera::new(
            Point::new_cartesian(0.0, 1.0, 0.0, 0.0),
            FourVector::new(1.0, 0.0, 0.0, 0.0, Point::new_cartesian(0.0, 1.0, 0.0, 0.0)),
            PI / 2.0,
            11,
            11,
            0.0,
            0.0,
            0.0,
            &geometry,
        )
        .unwrap();

        let base = camera.get_ray_for(3, 7);
        let centered = camera.get_ray_for_offset(3, 7, 0.5, 0.5);

        assert_eq!(centered.row, base.row);
        assert_eq!(centered.col, base.col);
        assert_abs_diff_eq!(centered.position.vector, base.position.vector);
        assert_abs_diff_eq!(centered.momentum.vector, base.momentum.vector);
    }

    #[test]
    fn test_get_ray_for_different_geometries_euclidean() {
        let position = Point::new_cartesian(0.0, 0.0, 1.0, 0.0);
        let velocity = FourVector::new(1.0, 0.0, 0.0, 0.0, position);

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
        let position_spherical = cartesian_to_spherical(&position);
        let camera_euclidean_spherical = Camera::new(
            position_spherical,
            FourVector::new(1.0, 0.0, 0.0, 0.0, position_spherical),
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
        let velocity = FourVector::new(1.0, 0.0, 0.0, 0.0, position);

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
        let position_spherical = cartesian_to_spherical(&position);
        let camera_euclidean_schwarzschild = Camera::new(
            position_spherical,
            FourVector::new(1.0, 0.0, 0.0, 0.0, position_spherical),
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

    // The camera traces the physical photon that arrives at it, so the ray
    // momentum must be past-directed: the physical energy measured by the
    // camera observer, signature[0] * <u, p>, is negative.
    #[test]
    fn test_camera_ray_is_past_directed_for_plus_minus_minus_minus_signature() {
        let position = cartesian_to_spherical(&Point::new_cartesian(0.0, 10.0, 0.0, 0.0));
        let geometry = Schwarzschild::new(0.0, 0.0);
        let velocity = FourVector::new(1.0, 0.0, 0.0, 0.0, position);
        let camera = Camera::new(
            position,
            velocity,
            PI / 2.0,
            11,
            11,
            0.0,
            0.0,
            0.0,
            &geometry,
        )
        .unwrap();

        let ray = camera.get_ray_for(5, 5);
        let orientation = geometry.signature()[0]
            * geometry.inner_product(&ray.position, &camera.velocity, &ray.momentum);
        assert!(orientation < 0.0);
    }

    #[test]
    fn test_camera_ray_is_past_directed_for_minus_plus_plus_plus_signature() {
        let geometry = Kerr::new(0.0, 0.0, 0.0);
        let position = KerrChart::point_from_cartesian(
            geometry.kerr_schild_chart(),
            &Point::new_cartesian(0.0, 10.0, 0.0, 0.0),
        );
        let velocity = FourVector::new(1.0, 0.0, 0.0, 0.0, position);
        let camera = Camera::new(
            position,
            velocity,
            PI / 2.0,
            11,
            11,
            0.0,
            0.0,
            0.0,
            &geometry,
        )
        .unwrap();

        let ray = camera.get_ray_for(5, 5);
        let orientation = geometry.signature()[0]
            * geometry.inner_product(&ray.position, &camera.velocity, &ray.momentum);
        assert!(orientation < 0.0);
    }

    /// Handedness of the geometry's spatial tetrad, measured in the Cartesian view. The point
    /// must be away from places where the view degenerates, so the measurement is meaningful.
    fn measured_handedness<G: crate::geometry::geometry::Geometry>(
        geometry: &G,
        position: &Point<G::Chart>,
    ) -> f64 {
        let tetrad = crate::geometry::geometry::Geometry::get_tetrad_at(geometry, position);
        let x = super::spatial_basis_vector_cartesian(position, &tetrad.x);
        let y = super::spatial_basis_vector_cartesian(position, &tetrad.y);
        let z = super::spatial_basis_vector_cartesian(position, &tetrad.z);
        let triple = x.dot(&y.cross(&z));
        assert!(
            triple.is_finite() && triple.abs() > 1e-6,
            "degenerate tetrad at {position:?}: triple product {triple}"
        );
        triple.signum()
    }

    /// The camera falls back to `tetrad_handedness_fallback` where it cannot measure the
    /// handedness. That is only correct if the fallback equals the handedness measured
    /// everywhere else; otherwise the image would be mirrored at those points.
    #[test]
    fn handedness_fallback_matches_the_measured_tetrad_handedness() {
        use crate::geometry::coordinates::chart::Chart;
        use crate::geometry::coordinates::spherical::RotatedSpherical;
        use crate::geometry::atlases::spherical::SphericalChart;
        use crate::geometry::spacetimes::kerr_bl::KerrBL;

        let spherical_points = [(12.0, 1.1, 0.4), (7.0, 2.5, -2.0), (30.0, 0.3, 3.0)];

        let flat = EuclideanSpace::new();
        for (x, y, z) in [(1.0, 2.0, 3.0), (-4.0, 0.5, -1.0)] {
            let p = Point::new_cartesian(0.0, x, y, z);
            assert_eq!(
                measured_handedness(&flat, &p),
                p.chart.tetrad_handedness_fallback()
            );
        }

        let flat_spherical = EuclideanSpaceSpherical::new();
        let schwarzschild = Schwarzschild::new(1.0, 1e-4);
        for rotation in [0.0, 90.0, 180.0, 270.0] {
            let chart = SphericalChart::Rotated(RotatedSpherical { rotation });
            for (r, theta, phi) in spherical_points {
                let p = Point::new(0.0, r, theta, phi, chart);
                let expected = chart.tetrad_handedness_fallback();
                assert_eq!(
                    measured_handedness(&flat_spherical, &p),
                    expected,
                    "flat, {p:?}"
                );
                assert_eq!(
                    measured_handedness(&schwarzschild, &p),
                    expected,
                    "Schwarzschild, {p:?}"
                );
            }
        }

        let kerr = Kerr::new(1.0, 0.6, 1e-4);
        for (x, y, z) in [(8.0, 3.0, 2.0), (-5.0, 6.0, -4.0)] {
            let p = Point::new(0.0, x, y, z, kerr.kerr_schild_chart());
            assert_eq!(
                measured_handedness(&kerr, &p),
                p.chart.tetrad_handedness_fallback(),
                "Kerr, {p:?}"
            );
        }

        let kerr_bl = KerrBL::new(1.0, 0.6, 1e-4);
        for (r, theta, phi) in spherical_points {
            let p = Point::new(0.0, r, theta, phi, kerr_bl.boyer_lindquist_chart());
            assert_eq!(
                measured_handedness(&kerr_bl, &p),
                p.chart.tetrad_handedness_fallback(),
                "KerrBL, {p:?}"
            );
        }
    }
}

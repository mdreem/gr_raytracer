use crate::geometry::coordinates;
use crate::geometry::four_vector::FourVector;
use crate::geometry::point::Point;
use nalgebra::Vector3;
use std::fmt::Debug;

pub trait Chart: Copy + PartialEq + Debug + Send + Sync + 'static {
    // TODO: Align all tetrads on one handedness convention (the (phi, theta, r) tetrads are
    // left-handed, the Kerr-Schild Gram-Schmidt one right-handed), then delete this and
    // `spatial_handedness` in camera.rs.
    /// Handedness of the spatial tetrad that geometries build in this chart; the camera
    /// falls back to it where the Cartesian triple product degenerates.
    fn tetrad_handedness_fallback(&self) -> f64;
}

/// The chart enum of a geometry is its atlas: it lists every chart, picks the best one for a
/// point and transforms points and vectors between them. Transitions are exact coordinate
/// changes of the same spacetime; the source chart is the one the point or vector carries.
#[allow(dead_code)] // Used once the integrator switches charts.
pub trait Atlas: Chart {
    fn best_chart(&self, position: &Point<Self>) -> Self;
    fn transition_point(&self, position: &Point<Self>, to: Self) -> Point<Self>;
    fn transition_four_vector(&self, vector: &FourVector<Self>, to: Self) -> FourVector<Self>;
}

pub trait SphericalConvertible: Chart {
    // The order of the components is: (r, theta, phi)
    fn point_as_spherical(point: &Point<Self>) -> Vector3<f64>; // TODO: Consider using Point here too.
    fn radial_distance_spatial_part_squared(point: &Point<Self>) -> f64;
}

pub trait CartesianConvertible: Chart {
    fn point_to_cartesian(point: &Point<Self>) -> Point<coordinates::euclidean::CartesianChart>;
    fn point_from_cartesian(
        chart: Self,
        point: &Point<coordinates::euclidean::CartesianChart>,
    ) -> Point<Self>;
    fn vector_to_cartesian(
        vector: &FourVector<Self>,
    ) -> FourVector<coordinates::euclidean::CartesianChart>;
}

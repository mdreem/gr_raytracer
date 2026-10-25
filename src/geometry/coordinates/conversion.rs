//! Renderer-facing Cartesian and spherical views of chart coordinates.
//! These views are distinct from exact atlas transitions within a spacetime.

use super::chart::Chart;
use crate::geometry::coordinates;
use crate::geometry::four_vector::FourVector;
use crate::geometry::point::Point;
use nalgebra::Vector3;

pub trait SphericalConvertible: Chart {
    // The order of the components is: (r, theta, phi)
    fn point_as_spherical(point: &Point<Self>) -> Vector3<f64>; // TODO: Consider using Point here too.
    fn radial_distance_spatial_part_squared(point: &Point<Self>) -> f64;
}

pub trait CartesianConvertible: Chart {
    fn point_to_cartesian(point: &Point<Self>) -> Point<coordinates::cartesian::CartesianChart>;
    fn point_from_cartesian(
        chart: Self,
        point: &Point<coordinates::cartesian::CartesianChart>,
    ) -> Point<Self>;
    fn vector_to_cartesian(
        vector: &FourVector<Self>,
    ) -> FourVector<coordinates::cartesian::CartesianChart>;
}

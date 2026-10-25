//! Stereographic coordinate views and the inversion Jacobian.

use crate::geometry::atlases::Atlas;
use crate::geometry::atlases::stereographic::StereographicChart;
use crate::geometry::atlases::stereographic::StereographicChart::WithNorthPole;
use crate::geometry::coordinates::cartesian::CartesianChart;
use crate::geometry::coordinates::chart::Chart;
use crate::geometry::coordinates::conversion::{CartesianConvertible, SphericalConvertible};
use crate::geometry::four_vector::FourVector;
use crate::geometry::point::Point;
use crate::geometry::spherical_coordinates_helper::cartesian_to_spherical;
use nalgebra::{Matrix3, Vector3};

pub(in crate::geometry) fn compute_u_squared<C: Chart>(point: &Point<C>) -> f64 {
    point.vector[1].powi(2) + point.vector[2].powi(2) + point.vector[3].powi(2)
}

/// The Cartesian view is the north chart, in units of the radius. The south pole of the sphere
/// (the origin of the south chart) lies at infinity there.
impl CartesianConvertible for StereographicChart {
    fn point_to_cartesian(point: &Point<Self>) -> Point<CartesianChart> {
        let north = point.chart.transition_point(point, WithNorthPole);
        Point::new_cartesian(north[0], north[1], north[2], north[3])
    }

    fn point_from_cartesian(chart: Self, point: &Point<CartesianChart>) -> Point<Self> {
        let north = Point::new(point[0], point[1], point[2], point[3], WithNorthPole);
        WithNorthPole.transition_point(&north, chart)
    }

    fn vector_to_cartesian(vector: &FourVector<Self>) -> FourVector<CartesianChart> {
        let north = vector
            .at
            .chart
            .transition_four_vector(vector, WithNorthPole);
        FourVector::new(
            north[0],
            north[1],
            north[2],
            north[3],
            Self::point_to_cartesian(&vector.at),
        )
    }
}

impl SphericalConvertible for StereographicChart {
    fn point_as_spherical(point: &Point<Self>) -> Vector3<f64> {
        let v = cartesian_to_spherical(&point.to_cartesian());
        Vector3::new(v[1], v[2], v[3])
    }

    /// A closed space has no celestial sphere, so nothing ever counts as escaped.
    fn radial_distance_spatial_part_squared(_point: &Point<Self>) -> f64 {
        0.0
    }
}

pub(in crate::geometry) fn jacobian(position: &Point<StereographicChart>) -> Matrix3<f64> {
    let (v1, v2, v3) = (position.vector[1], position.vector[2], position.vector[3]);

    let l_sq = v1.powi(2) + v2.powi(2) + v3.powi(2);
    let l_quad = l_sq * l_sq;

    #[rustfmt::skip]
    let data = [
        l_sq - 2.0*v1*v1, - 2.0*v1*v2, - 2.0*v1*v3,
         -2.0*v2*v1, l_sq - 2.0*v2*v2, - 2.0*v2*v3,
        -2.0*v3*v1, -2.0*v3*v2, l_sq - 2.0*v3*v3,
    ];
    Matrix3::from_row_slice(&data) / l_quad
}

use crate::geometry::coordinates;
use crate::geometry::coordinates::atlas::Atlas;
use crate::geometry::coordinates::atlas::{CartesianConvertible, Chart, SphericalConvertible};
use crate::geometry::coordinates::boyer_lindquist::{
    BoyerLindquistData, KerrSchildData, bl_to_ks, ks_to_bl,
};
use crate::geometry::coordinates::spherical::RotatedSpherical;
use crate::geometry::four_vector::FourVector;
use crate::geometry::spherical_coordinates_helper::{
    cartesian_to_spherical, rotate_about_x, spherical_to_cartesian,
};
use nalgebra::{Vector3, Vector4};
use std::cmp::PartialEq;
use std::f64::consts::PI;
use std::ops::{Add, Index, Neg};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Point<C: Chart> {
    pub chart: C,
    pub vector: Vector4<f64>,
}

impl Neg for Point<coordinates::euclidean::CartesianChart> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Point {
            chart: self.chart,
            vector: self.vector.neg(),
        }
    }
}

impl<C: Chart> Add for Point<C> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Point {
            chart: self.chart,
            vector: self.vector + rhs.vector,
        }
    }
}

impl<C: Chart + Copy> Add for &Point<C> {
    type Output = Point<C>;

    fn add(self, rhs: Self) -> Self::Output {
        Point {
            chart: self.chart,
            vector: self.vector + rhs.vector,
        }
    }
}

/// Wrap any angle to the half-open interval [0, Pi)
fn wrap_theta(theta: f64) -> f64 {
    theta.rem_euclid(PI)
}

/// Wrap any angle to the half-open interval (-Pi/2, Pi/2]
fn wrap_phi(phi: f64) -> f64 {
    (phi + PI).rem_euclid(2.0 * PI) - PI
}

impl Point<coordinates::euclidean::CartesianChart> {
    pub fn new_cartesian(
        x0: f64,
        x1: f64,
        x2: f64,
        x3: f64,
    ) -> Point<coordinates::euclidean::CartesianChart> {
        Point {
            chart: coordinates::euclidean::CartesianChart::Cartesian,
            vector: Vector4::new(x0, x1, x2, x3),
        }
    }
}

impl Point<coordinates::spherical::SphericalChart> {
    pub fn new_spherical(
        t: f64,
        r: f64,
        theta: f64,
        phi: f64,
    ) -> Point<coordinates::spherical::SphericalChart> {
        Point {
            chart: coordinates::spherical::SphericalChart::Rotated(RotatedSpherical {
                rotation: 0.0,
            }),
            vector: Vector4::new(t, r, theta, phi),
        }
    }
}

impl Point<coordinates::boyer_lindquist::KerrChart> {
    pub fn new_boyer_lindquist(
        a: f64,
        r_s: f64,
        t: f64,
        r: f64,
        theta: f64,
        phi: f64,
    ) -> Point<coordinates::boyer_lindquist::KerrChart> {
        Point {
            chart: coordinates::boyer_lindquist::KerrChart::BoyerLindquist(BoyerLindquistData {
                a,
                r_s,
            }),
            vector: Vector4::new(t, r, theta, phi),
        }
    }
}

impl<C: Chart> Point<C> {
    pub fn new(x0: f64, x1: f64, x2: f64, x3: f64, chart: C) -> Point<C> {
        Point {
            chart,
            vector: Vector4::new(x0, x1, x2, x3),
        }
    }

    pub fn new_from_vector(vector: Vector4<f64>, chart: C) -> Point<C> {
        Point { chart, vector }
    }

    #[allow(dead_code)] // For testing
    pub fn get_as_vector(self) -> Vector4<f64> {
        self.vector
    }
}

impl<C: SphericalConvertible> Point<C> {
    // The order of the components is: (r, theta, phi)
    pub fn get_as_spherical(&self) -> Vector3<f64> {
        C::point_as_spherical(self)
    }

    pub fn radial_distance_spatial_part_squared(&self) -> f64 {
        C::radial_distance_spatial_part_squared(self)
    }
}

impl<C: CartesianConvertible> Point<C> {
    pub fn to_cartesian(&self) -> Point<coordinates::euclidean::CartesianChart> {
        C::point_to_cartesian(self)
    }

    pub fn get_as_cartesian(&self) -> Vector3<f64> {
        let v = self.to_cartesian();
        Vector3::new(v[1], v[2], v[3])
    }

    pub fn get_spatial_vector_cartesian(&self) -> Vector3<f64> {
        let v = self.to_cartesian();
        Vector3::new(v[1], v[2], v[3])
    }
}

impl SphericalConvertible for coordinates::euclidean::CartesianChart {
    fn point_as_spherical(point: &Point<Self>) -> Vector3<f64> {
        let v = cartesian_to_spherical(point);
        Vector3::new(v[1], v[2], v[3])
    }

    fn radial_distance_spatial_part_squared(point: &Point<Self>) -> f64 {
        let v = point.vector;
        v[1] * v[1] + v[2] * v[2] + v[3] * v[3]
    }
}

impl SphericalConvertible for coordinates::spherical::SphericalChart {
    /// Angles in the physical frame, also for rotated charts.
    fn point_as_spherical(point: &Point<Self>) -> Vector3<f64> {
        if point.chart.rotation() != 0.0 {
            let v = cartesian_to_spherical(&spherical_to_cartesian(point));
            return Vector3::new(v[1], v[2], v[3]);
        }
        Vector3::new(
            point.vector[1],
            wrap_theta(point.vector[2]),
            wrap_phi(point.vector[3]),
        )
    }

    fn radial_distance_spatial_part_squared(point: &Point<Self>) -> f64 {
        let v = point.vector;
        let r = v[1]; // v[1] represents the radial coordinate r.
        r * r
    }
}

fn kerr_schild_as_cartesian(
    point: &Point<coordinates::boyer_lindquist::KerrChart>,
) -> Point<coordinates::euclidean::CartesianChart> {
    Point::new_from_vector(
        point.vector,
        coordinates::euclidean::CartesianChart::Cartesian,
    )
}

// TODO: Rethink once the atlas exists: untwisted BL embedding vs. Kerr-Schild; to_cartesian must agree across KerrChart variants.
impl SphericalConvertible for coordinates::boyer_lindquist::KerrChart {
    fn point_as_spherical(point: &Point<Self>) -> Vector3<f64> {
        match point.chart {
            coordinates::boyer_lindquist::KerrChart::BoyerLindquist(_) => Vector3::new(
                point.vector[1],
                wrap_theta(point.vector[2]),
                wrap_phi(point.vector[3]),
            ),
            coordinates::boyer_lindquist::KerrChart::KerrSchild(_) => {
                coordinates::euclidean::CartesianChart::point_as_spherical(
                    &kerr_schild_as_cartesian(point),
                )
            }
        }
    }

    fn radial_distance_spatial_part_squared(point: &Point<Self>) -> f64 {
        match point.chart {
            coordinates::boyer_lindquist::KerrChart::BoyerLindquist(_) => {
                point.vector[1] * point.vector[1]
            }
            coordinates::boyer_lindquist::KerrChart::KerrSchild(_) => {
                coordinates::euclidean::CartesianChart::radial_distance_spatial_part_squared(
                    &kerr_schild_as_cartesian(point),
                )
            }
        }
    }
}

impl CartesianConvertible for coordinates::euclidean::CartesianChart {
    fn point_to_cartesian(point: &Point<Self>) -> Point<coordinates::euclidean::CartesianChart> {
        *point
    }

    fn point_from_cartesian(
        _chart: Self,
        point: &Point<coordinates::euclidean::CartesianChart>,
    ) -> Point<Self> {
        *point
    }

    fn vector_to_cartesian(
        vector: &FourVector<Self>,
    ) -> FourVector<coordinates::euclidean::CartesianChart> {
        *vector
    }
}

impl CartesianConvertible for coordinates::spherical::SphericalChart {
    fn point_to_cartesian(point: &Point<Self>) -> Point<coordinates::euclidean::CartesianChart> {
        spherical_to_cartesian(point)
    }

    fn point_from_cartesian(
        chart: Self,
        point: &Point<coordinates::euclidean::CartesianChart>,
    ) -> Point<Self> {
        // Physical -> coordinate frame: apply the chart's rotation.
        let (x, y, z) = rotate_about_x(point[1], point[2], point[3], chart.rotation());
        let rotated = Point::new_cartesian(point[0], x, y, z);
        Point::new_from_vector(cartesian_to_spherical(&rotated).vector, chart)
    }

    fn vector_to_cartesian(
        vector: &FourVector<Self>,
    ) -> FourVector<coordinates::euclidean::CartesianChart> {
        let at = vector.at.to_cartesian();
        let (r, theta, phi) = (vector.at[1], vector.at[2], vector.at[3]);
        let (st, ct) = theta.sin_cos();
        let (sp, cp) = phi.sin_cos();
        let (dr, dtheta, dphi) = (vector[1], vector[2], vector[3]);

        let x = st * cp * dr + r * ct * cp * dtheta - r * st * sp * dphi;
        let y = st * sp * dr + r * ct * sp * dtheta + r * st * cp * dphi;
        let z = ct * dr - r * st * dtheta;

        // Coordinate -> physical frame: undo the chart's rotation, as `spherical_to_cartesian` does for points.
        let (x, y, z) = rotate_about_x(x, y, z, -vector.at.chart.rotation());
        FourVector::new(vector[0], x, y, z, at)
    }
}

/// The Cartesian view of both Kerr charts is the Kerr-Schild one: Boyer-Lindquist points and
/// vectors go through the same (twisted) map as the atlas transition, so the view does not
/// depend on the chart a step happens to be in.
impl CartesianConvertible for coordinates::boyer_lindquist::KerrChart {
    fn point_to_cartesian(point: &Point<Self>) -> Point<coordinates::euclidean::CartesianChart> {
        match point.chart {
            coordinates::boyer_lindquist::KerrChart::BoyerLindquist(BoyerLindquistData {
                a,
                r_s,
            }) => Point::new_from_vector(
                bl_to_ks(a, r_s, &point.vector),
                coordinates::euclidean::CartesianChart::Cartesian,
            ),
            coordinates::boyer_lindquist::KerrChart::KerrSchild(_) => {
                kerr_schild_as_cartesian(point)
            }
        }
    }

    fn point_from_cartesian(
        chart: Self,
        point: &Point<coordinates::euclidean::CartesianChart>,
    ) -> Point<Self> {
        match chart {
            coordinates::boyer_lindquist::KerrChart::BoyerLindquist(BoyerLindquistData {
                a,
                r_s,
            }) => Point::new_from_vector(ks_to_bl(a, r_s, &point.vector).0, chart),
            coordinates::boyer_lindquist::KerrChart::KerrSchild(_) => {
                Point::new_from_vector(point.vector, chart)
            }
        }
    }

    fn vector_to_cartesian(
        vector: &FourVector<Self>,
    ) -> FourVector<coordinates::euclidean::CartesianChart> {
        let in_kerr_schild = match vector.at.chart {
            coordinates::boyer_lindquist::KerrChart::BoyerLindquist(BoyerLindquistData {
                a,
                r_s,
            }) => vector.at.chart.transition_four_vector(
                vector,
                coordinates::boyer_lindquist::KerrChart::KerrSchild(KerrSchildData { a, r_s }),
            ),
            coordinates::boyer_lindquist::KerrChart::KerrSchild(_) => *vector,
        };
        FourVector::new(
            in_kerr_schild[0],
            in_kerr_schild[1],
            in_kerr_schild[2],
            in_kerr_schild[3],
            kerr_schild_as_cartesian(&in_kerr_schild.at),
        )
    }
}

impl<C: Chart> Index<usize> for Point<C> {
    type Output = f64;

    fn index(&self, index: usize) -> &Self::Output {
        &self.vector[index]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_rotated_spherical_chart_round_trips_through_cartesian() {
        use crate::geometry::coordinates::atlas::Atlas;
        use crate::geometry::coordinates::spherical::{RotatedSpherical, SphericalChart};
        use crate::geometry::four_vector::FourVector;

        let physical = Point::new_cartesian(1.0, 3.0, -2.0, 4.0);
        let unrotated = SphericalChart::Rotated(RotatedSpherical { rotation: 0.0 });
        let rotated = SphericalChart::Rotated(RotatedSpherical { rotation: 37.0 });

        let p = SphericalChart::point_from_cartesian(rotated, &physical);
        assert_eq!(p.chart, rotated);
        assert_abs_diff_eq!(p.to_cartesian().vector, physical.vector, epsilon = 1e-12);

        // The same physical vector, expressed in both charts, maps to the same Cartesian vector.
        let q = SphericalChart::point_from_cartesian(unrotated, &physical);
        let v_q = FourVector::new(1.0, 0.3, -0.05, 0.02, q);
        let v_p = q.chart.transition_four_vector(&v_q, rotated);
        assert_abs_diff_eq!(v_p.at.vector, p.vector, epsilon = 1e-12);
        assert_abs_diff_eq!(
            v_p.to_cartesian().vector,
            v_q.to_cartesian().vector,
            epsilon = 1e-12
        );
    }

    #[test]
    fn test_rotated_spherical_chart_reports_physical_angles() {
        use crate::geometry::coordinates::spherical::{RotatedSpherical, SphericalChart};

        let physical = Point::new_cartesian(0.0, 3.0, -2.0, 4.0);
        let expected = cartesian_to_spherical(&physical);
        for rotation in [0.0, 90.0, 180.0, 270.0] {
            let chart = SphericalChart::Rotated(RotatedSpherical { rotation });
            let p = SphericalChart::point_from_cartesian(chart, &physical);
            let spherical = p.get_as_spherical();
            assert_abs_diff_eq!(spherical[0], expected[1], epsilon = 1e-12);
            assert_abs_diff_eq!(spherical[1], expected[2], epsilon = 1e-12);
            assert_abs_diff_eq!(spherical[2], expected[3], epsilon = 1e-12);
        }
    }

    #[test]
    fn test_boyer_lindquist_to_cartesian_a_zero() {
        // When a=0, BL space matches standard spherical (the time keeps the Kerr-Schild twist).
        let bl = Point::new_boyer_lindquist(0.0, 1.0, 0.0, 5.0, 1.2, 0.8);
        let sph = Point::new_spherical(0.0, 5.0, 1.2, 0.8);
        assert_abs_diff_eq!(
            bl.get_spatial_vector_cartesian(),
            sph.get_spatial_vector_cartesian(),
            epsilon = 1e-12
        );
    }

    /// The Cartesian view of a BL point or vector is that of the same event in the Kerr-Schild
    /// chart, so it does not jump when a ray switches charts, and it round-trips.
    #[test]
    fn test_boyer_lindquist_cartesian_view_agrees_with_kerr_schild() {
        use crate::geometry::coordinates::boyer_lindquist::KerrChart;

        for a in [0.5, -0.3, 0.7] {
            let bl_chart = KerrChart::BoyerLindquist(BoyerLindquistData { a, r_s: 1.0 });
            let ks_chart = KerrChart::KerrSchild(KerrSchildData { a, r_s: 1.0 });
            let bl = Point::new(0.4, 5.0, 1.2, 0.8, bl_chart);
            let ks = bl_chart.transition_point(&bl, ks_chart);
            assert_abs_diff_eq!(
                bl.to_cartesian().vector,
                ks.to_cartesian().vector,
                epsilon = 1e-12
            );

            let back = KerrChart::point_from_cartesian(bl_chart, &bl.to_cartesian());
            assert_abs_diff_eq!(back.vector, bl.vector, epsilon = 1e-10);

            let v = FourVector::new(1.0, -0.3, 0.02, 0.05, bl);
            let v_ks = bl_chart.transition_four_vector(&v, ks_chart);
            assert_abs_diff_eq!(
                v.to_cartesian().vector,
                v_ks.to_cartesian().vector,
                epsilon = 1e-12
            );
        }
    }

    #[test]
    fn test_boyer_lindquist_radial_distance() {
        let bl = Point::new_boyer_lindquist(0.5, 1.0, 0.0, 7.0, 1.0, 2.0);
        assert_abs_diff_eq!(
            bl.radial_distance_spatial_part_squared(),
            49.0,
            epsilon = 1e-12
        );
    }
}

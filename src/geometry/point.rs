use crate::geometry::coordinates;
use crate::geometry::coordinates::atlas::{CartesianConvertible, Chart, SphericalConvertible};
use crate::geometry::coordinates::boyer_lindquist::BoyerLindquistData;
use crate::geometry::coordinates::spherical::RotatedSpherical;
use crate::geometry::four_vector::FourVector;
use crate::geometry::spherical_coordinates_helper::{
    cartesian_to_boyer_lindquist, cartesian_to_spherical, frame_rotation_deg, rotate_about_x,
    spherical_to_cartesian,
};
use nalgebra::{Matrix4, Vector3, Vector4};
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
    fn point_as_spherical(point: &Point<Self>) -> Vector3<f64> {
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

    // TODO: Take into account the chart's rotation (cartesian_to_spherical applies the global frame rotation).
    fn point_from_cartesian(
        chart: Self,
        point: &Point<coordinates::euclidean::CartesianChart>,
    ) -> Point<Self> {
        Point::new_from_vector(cartesian_to_spherical(point).vector, chart)
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

        // Coordinate -> physical frame: undo the pole rotation, as `spherical_to_cartesian` does for points.
        let (x, y, z) = rotate_about_x(x, y, z, -frame_rotation_deg());
        FourVector::new(vector[0], x, y, z, at)
    }
}

// TODO: Rethink once the atlas exists: untwisted BL embedding vs. Kerr-Schild; to_cartesian must agree across KerrChart variants.
impl CartesianConvertible for coordinates::boyer_lindquist::KerrChart {
    fn point_to_cartesian(point: &Point<Self>) -> Point<coordinates::euclidean::CartesianChart> {
        let a = match point.chart {
            coordinates::boyer_lindquist::KerrChart::BoyerLindquist(BoyerLindquistData {
                a,
                ..
            }) => a,
            coordinates::boyer_lindquist::KerrChart::KerrSchild(_) => {
                return Point::new_cartesian(
                    point.vector[0],
                    point.vector[1],
                    point.vector[2],
                    point.vector[3],
                );
            }
        };

        let t = point.vector[0];
        let r = point.vector[1];
        let theta = point.vector[2];
        let phi = point.vector[3];
        let x = (r * phi.cos() - a * phi.sin()) * theta.sin();
        let y = (r * phi.sin() + a * phi.cos()) * theta.sin();
        let z = r * theta.cos();
        Point::new_cartesian(t, x, y, z)
    }

    fn point_from_cartesian(
        chart: Self,
        point: &Point<coordinates::euclidean::CartesianChart>,
    ) -> Point<Self> {
        match chart {
            coordinates::boyer_lindquist::KerrChart::BoyerLindquist(BoyerLindquistData {
                a,
                r_s,
            }) => cartesian_to_boyer_lindquist(a, r_s, point),
            coordinates::boyer_lindquist::KerrChart::KerrSchild(_) => {
                Point::new_from_vector(point.vector, chart)
            }
        }
    }

    fn vector_to_cartesian(
        vector: &FourVector<Self>,
    ) -> FourVector<coordinates::euclidean::CartesianChart> {
        let at = vector.at.to_cartesian();
        let a = match vector.at.chart {
            coordinates::boyer_lindquist::KerrChart::BoyerLindquist(BoyerLindquistData {
                a,
                ..
            }) => a,
            coordinates::boyer_lindquist::KerrChart::KerrSchild(_) => {
                return FourVector::new(vector[0], vector[1], vector[2], vector[3], at);
            }
        };

        // Derivative of the embedding in `point_to_cartesian` (no t/φ twist), so the
        // components match `at`.
        let (r, theta, phi) = (vector.at[1], vector.at[2], vector.at[3]);
        let (st, ct) = theta.sin_cos();
        let (sp, cp) = phi.sin_cos();

        #[rustfmt::skip]
        let j = Matrix4::from_row_slice(&[
            1.0, 0.0,     0.0,                    0.0,
            0.0, st * cp, (r * cp - a * sp) * ct, (-r * sp - a * cp) * st,
            0.0, st * sp, (r * sp + a * cp) * ct, (r * cp - a * sp) * st,
            0.0, ct,      -r * st,                0.0,
        ]);

        let v = j * vector.vector;
        FourVector::new(v[0], v[1], v[2], v[3], at)
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
    fn test_boyer_lindquist_to_cartesian_a_zero() {
        // When a=0, BL should match standard spherical
        let bl = Point::new_boyer_lindquist(0.0, 1.0, 0.0, 5.0, 1.2, 0.8);
        let sph = Point::new_spherical(0.0, 5.0, 1.2, 0.8);
        let bl_cart = bl.to_cartesian();
        let sph_cart = sph.to_cartesian();
        assert_abs_diff_eq!(bl_cart.vector, sph_cart.vector, epsilon = 1e-12);
    }

    #[test]
    fn test_boyer_lindquist_to_cartesian_nonzero_a() {
        let a = 0.5_f64;
        let r = 5.0_f64;
        let theta = 1.2_f64;
        let phi = 0.8_f64;
        let bl = Point::new_boyer_lindquist(a, 1.0, 0.0, r, theta, phi);
        let cart = bl.to_cartesian();
        // Expected values computed independently: x = (r cosφ - a sinφ) sinθ, etc.
        // cos(0.8) ≈ 0.6967067093471654, sin(0.8) ≈ 0.7173560908995228
        // sin(1.2) ≈ 0.9320390309173473, cos(1.2) ≈ 0.3623577544766736
        assert_abs_diff_eq!(cart[1], 2.91248746519832302226_f64, epsilon = 1e-10);
        assert_abs_diff_eq!(cart[2], 3.66769851865865170737_f64, epsilon = 1e-10);
        assert_abs_diff_eq!(cart[3], 1.81178877238336810684_f64, epsilon = 1e-10);
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

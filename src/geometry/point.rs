use crate::geometry::point::CoordinateSystem::{BoyerLindquist, Cartesian, Spherical};
use crate::geometry::spherical_coordinates_helper::{
    cartesian_to_spherical, spherical_to_cartesian,
};
use nalgebra::{Vector3, Vector4};
use std::cmp::PartialEq;
use std::f64::consts::PI;
use std::ops::{Add, Index, Neg};

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CoordinateSystem {
    Cartesian,
    Spherical,
    /// Boyer-Lindquist coordinates (t, r, θ, φ) for the Kerr metric.
    ///
    /// The Cartesian embedding used is the null-tetrad convention consistent with
    /// the Kerr-Schild form used in `geometry/kerr.rs`:
    ///   x = (r cosφ − a sinφ) sinθ
    ///   y = (r sinφ + a cosφ) sinθ
    ///   z = r cosθ
    ///
    /// Note: this differs from the oblate-spheroidal convention
    /// (x = √(r²+a²) sinθ cosφ) used in some references.
    BoyerLindquist { a: f64 },
}

#[derive(Clone, Copy, Debug)]
pub struct Point {
    pub coordinate_system: CoordinateSystem,
    pub vector: Vector4<f64>,
}

impl Neg for Point {
    type Output = Self;

    fn neg(self) -> Self::Output {
        debug_assert_eq!(self.coordinate_system, Cartesian); // Negation is only works like this for Cartesian coordinates
        Point {
            coordinate_system: self.coordinate_system,
            vector: self.vector.neg(),
        }
    }
}

impl Add for Point {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        debug_assert_eq!(self.coordinate_system, rhs.coordinate_system);
        Point {
            coordinate_system: self.coordinate_system,
            vector: self.vector + rhs.vector,
        }
    }
}

impl Add for &Point {
    type Output = Point;

    fn add(self, rhs: Self) -> Self::Output {
        debug_assert_eq!(self.coordinate_system, rhs.coordinate_system);
        Point {
            coordinate_system: self.coordinate_system,
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

impl Point {
    pub fn new(x0: f64, x1: f64, x2: f64, x3: f64, coordinate_system: CoordinateSystem) -> Point {
        Point {
            coordinate_system,
            vector: Vector4::new(x0, x1, x2, x3),
        }
    }

    pub fn new_from_vector(vector: Vector4<f64>, coordinate_system: CoordinateSystem) -> Point {
        Point {
            coordinate_system,
            vector,
        }
    }

    pub fn new_cartesian(x0: f64, x1: f64, x2: f64, x3: f64) -> Point {
        Point {
            coordinate_system: Cartesian,
            vector: Vector4::new(x0, x1, x2, x3),
        }
    }

    pub fn new_spherical(t: f64, r: f64, theta: f64, phi: f64) -> Point {
        Point {
            coordinate_system: Spherical,
            vector: Vector4::new(t, r, theta, phi),
        }
    }

    pub fn new_boyer_lindquist(a: f64, t: f64, r: f64, theta: f64, phi: f64) -> Point {
        Point {
            coordinate_system: BoyerLindquist { a },
            vector: Vector4::new(t, r, theta, phi),
        }
    }

    // TODO: deduplicate
    pub fn get_spatial_vector_cartesian(self) -> Vector3<f64> {
        match self.coordinate_system {
            Cartesian => Vector3::new(self.vector[1], self.vector[2], self.vector[3]),
            Spherical => {
                let v = spherical_to_cartesian(&self);
                Vector3::new(v[1], v[2], v[3])
            }
            BoyerLindquist { .. } => {
                let v = self.to_cartesian();
                Vector3::new(v[1], v[2], v[3])
            }
        }
    }

    pub fn to_cartesian(self) -> Point {
        match self.coordinate_system {
            Cartesian => self,
            Spherical => spherical_to_cartesian(&self),
            BoyerLindquist { a } => {
                let t = self.vector[0];
                let r = self.vector[1];
                let theta = self.vector[2];
                let phi = self.vector[3];
                let x = (r * phi.cos() - a * phi.sin()) * theta.sin();
                let y = (r * phi.sin() + a * phi.cos()) * theta.sin();
                let z = r * theta.cos();
                Point::new_cartesian(t, x, y, z)
            }
        }
    }

    // The order of the components is: (r, theta, phi)
    pub fn get_as_spherical(self) -> Vector3<f64> {
        match self.coordinate_system {
            Cartesian => {
                let v = cartesian_to_spherical(&self);
                Vector3::new(v[1], v[2], v[3])
            }
            Spherical | BoyerLindquist { .. } => Vector3::new(
                self.vector[1],
                wrap_theta(self.vector[2]),
                wrap_phi(self.vector[3]),
            ),
        }
    }

    pub fn radial_distance_spatial_part_squared(&self) -> f64 {
        let v = self.vector;
        match self.coordinate_system {
            Cartesian => v[1] * v[1] + v[2] * v[2] + v[3] * v[3],
            Spherical | BoyerLindquist { .. } => {
                // In spherical/BL coordinates, the radial distance is just r^2.
                let r = v[1]; // v[1] represents the radial coordinate r.
                r * r
            }
        }
    }

    #[allow(dead_code)] // For testing
    pub fn get_as_vector(self) -> Vector4<f64> {
        self.vector
    }
}

impl Index<usize> for Point {
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
        let bl = Point::new(0.0, 5.0, 1.2, 0.8, CoordinateSystem::BoyerLindquist { a: 0.0 });
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
        let bl = Point::new(0.0, r, theta, phi, CoordinateSystem::BoyerLindquist { a });
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
        let bl = Point::new(0.0, 7.0, 1.0, 2.0, CoordinateSystem::BoyerLindquist { a: 0.5 });
        assert_abs_diff_eq!(bl.radial_distance_spatial_part_squared(), 49.0, epsilon = 1e-12);
    }
}

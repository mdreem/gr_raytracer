use crate::geometry::point::CoordinateSystem::{Cartesian, Spherical};
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

    // TODO: deduplicate
    pub fn get_spatial_vector_cartesian(self) -> Vector3<f64> {
        match self.coordinate_system {
            Cartesian => Vector3::new(self.vector[1], self.vector[2], self.vector[3]),
            Spherical => {
                let v = spherical_to_cartesian(&self);
                Vector3::new(v[1], v[2], v[3])
            }
        }
    }

    pub fn to_cartesian(self) -> Point {
        match self.coordinate_system {
            Cartesian => self,
            Spherical => spherical_to_cartesian(&self),
        }
    }

    // The order of the components is: (r, theta, phi)
    pub fn get_as_spherical(self) -> Vector3<f64> {
        match self.coordinate_system {
            Cartesian => {
                let v = cartesian_to_spherical(&self);
                Vector3::new(v[1], v[2], v[3])
            }
            Spherical => Vector3::new(
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
            Spherical => {
                // In spherical coordinates, the radial distance is just r^2.
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

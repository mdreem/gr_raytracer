use crate::geometry::spherical_coordinates_helper::cartesian_to_spherical;
use nalgebra::{Vector3, Vector4};
use std::cmp::PartialEq;
use std::ops::Index;

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
            coordinate_system: CoordinateSystem::Cartesian,
            vector: Vector4::new(x0, x1, x2, x3),
        }
    }

    pub fn new_spherical(t: f64, r: f64, theta: f64, phi: f64) -> Point {
        Point {
            coordinate_system: CoordinateSystem::Spherical,
            vector: Vector4::new(t, r, theta, phi),
        }
    }

    pub fn get_spatial_vector(self) -> Vector3<f64> {
        Vector3::new(self.vector[1], self.vector[2], self.vector[3])
    }

    // The order of the components is: (r, theta, phi)
    pub fn get_as_spherical(self) -> Vector3<f64> {
        match self.coordinate_system {
            CoordinateSystem::Cartesian => {
                let v = cartesian_to_spherical(&self);
                Vector3::new(v[1], v[2], v[3])
            }
            CoordinateSystem::Spherical => {
                Vector3::new(self.vector[1], self.vector[2], self.vector[3])
            }
        }
    }

    pub fn radial_distance_spatial_part_squared(&self) -> f64 {
        let v = self.vector;
        match self.coordinate_system {
            CoordinateSystem::Cartesian => v[1] * v[1] + v[2] * v[2] + v[3] * v[3],
            CoordinateSystem::Spherical => {
                // In spherical coordinates, the radial distance is just r^2.
                let r = v[1]; // v[1] represents the radial coordinate r.
                r * r
            }
        }
    }

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

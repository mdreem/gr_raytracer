use crate::spherical_coordinates_helper::cartesian_to_spherical;
use nalgebra::{Vector3, Vector4};
use std::ops::{Add, Div, Mul, Neg};

#[derive(Clone, Copy, Debug)]
pub enum CoordinateSystem {
    Cartesian,
    Spherical,
}

#[derive(Clone, Copy, Debug)]
pub struct FourVector {
    pub coordinate_system: CoordinateSystem,
    pub vector: Vector4<f64>,
}

impl Neg for FourVector {
    type Output = Self;

    fn neg(self) -> Self::Output {
        FourVector {
            coordinate_system: self.coordinate_system,
            vector: self.vector.neg(),
        }
    }
}

impl Add for FourVector {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        FourVector {
            coordinate_system: self.coordinate_system,
            vector: self.vector + rhs.vector,
        }
    }
}

impl Mul<FourVector> for f64 {
    type Output = FourVector;

    fn mul(self, f: FourVector) -> FourVector {
        FourVector {
            coordinate_system: f.coordinate_system,
            vector: self * f.vector,
        }
    }
}

impl Mul<f64> for FourVector {
    type Output = Self;

    fn mul(self, f: f64) -> Self {
        FourVector {
            coordinate_system: self.coordinate_system,
            vector: f * self.vector,
        }
    }
}

impl Div<f64> for FourVector {
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        FourVector {
            coordinate_system: self.coordinate_system,
            vector: self.vector / rhs,
        }
    }
}

impl FourVector {
    pub fn new(
        x0: f64,
        x1: f64,
        x2: f64,
        x3: f64,
        coordinate_system: CoordinateSystem,
    ) -> FourVector {
        FourVector {
            coordinate_system,
            vector: Vector4::new(x0, x1, x2, x3),
        }
    }

    pub fn new_cartesian(x0: f64, x1: f64, x2: f64, x3: f64) -> FourVector {
        FourVector {
            coordinate_system: CoordinateSystem::Cartesian,
            vector: Vector4::new(x0, x1, x2, x3),
        }
    }

    pub fn new_spherical(t: f64, r: f64, theta: f64, phi: f64) -> FourVector {
        FourVector {
            coordinate_system: CoordinateSystem::Spherical,
            vector: Vector4::new(t, r, theta, phi),
        }
    }

    pub fn get_as_vector(self) -> Vector4<f64> {
        self.vector
    }

    pub fn get_spatial_vector(self) -> Vector3<f64> {
        Vector3::new(self.vector[1], self.vector[2], self.vector[3])
    }

    // The order of the components is: (r, theta, phi)
    pub fn get_as_spherical(self) -> Vector3<f64> {
        match self.coordinate_system {
            CoordinateSystem::Cartesian => {
                let v = cartesian_to_spherical(&self.vector);
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
}

#[cfg(test)]
mod tests {
    use crate::euclidean::EuclideanSpace;
    use crate::four_vector::FourVector;
    use crate::geometry::{Geometry, InnerProduct};
    use approx::assert_abs_diff_eq;
    use nalgebra::Vector4;

    #[test]
    fn test_multiplication_self() {
        let geometry = EuclideanSpace::new();
        let v1 = FourVector::new_cartesian(1.0, 2.0, 3.0, 4.0);
        assert_abs_diff_eq!(
            geometry.inner_product(&Vector4::new(0.0, 0.0, 0.0, 0.0), &v1, &v1),
            -28.0
        );
    }

    #[test]
    fn test_multiplication_different() {
        let geometry = EuclideanSpace::new();
        let v1 = FourVector::new_cartesian(1.0, 2.0, 3.0, 4.0);
        let v2 = FourVector::new_cartesian(5.0, 6.0, 7.0, 8.0);

        assert_abs_diff_eq!(
            geometry.inner_product(&Vector4::new(0.0, 0.0, 0.0, 0.0), &v1, &v2),
            -60.0
        );
    }
}

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
    coordinate_system: CoordinateSystem,
    vector: Vector4<f64>,
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

impl Mul<FourVector> for FourVector {
    type Output = f64;

    // Scalar multiplication with metric. For now with a flat metric and signature (+---)
    fn mul(self, v: Self) -> f64 {
        1.0 * self.vector[0] * v.vector[0]
            + (-1.0) * self.vector[1] * v.vector[1]
            + (-1.0) * self.vector[2] * v.vector[2]
            + (-1.0) * self.vector[3] * v.vector[3]
    }
}

impl FourVector {
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
            CoordinateSystem::Cartesian => cartesian_to_spherical(&self.vector),
            CoordinateSystem::Spherical => {
                Vector3::new(self.vector[1], self.vector[2], self.vector[3])
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::four_vector::FourVector;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_multiplication_self() {
        let v1 = FourVector::new_cartesian(1.0, 2.0, 3.0, 4.0);
        assert_abs_diff_eq!(v1 * v1, -28.0);
    }

    #[test]
    fn test_multiplication_different() {
        let v1 = FourVector::new_cartesian(1.0, 2.0, 3.0, 4.0);
        let v2 = FourVector::new_cartesian(5.0, 6.0, 7.0, 8.0);

        assert_abs_diff_eq!(v1 * v2, -60.0);
    }
}

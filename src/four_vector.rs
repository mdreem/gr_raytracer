use nalgebra::Vector4;
use std::ops::{Add, Div, Mul, Neg};

#[derive(Clone, Copy, Debug)]
pub struct FourVector {
    vector: Vector4<f64>,
}

impl Neg for FourVector {
    type Output = Self;

    fn neg(self) -> Self::Output {
        FourVector {
            vector: self.vector.neg(),
        }
    }
}

impl Add for FourVector {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        FourVector {
            vector: self.vector + rhs.vector,
        }
    }
}

impl Mul<FourVector> for f64 {
    type Output = FourVector;

    fn mul(self, f: FourVector) -> FourVector {
        FourVector {
            vector: self * f.vector,
        }
    }
}

impl Mul<f64> for FourVector {
    type Output = Self;

    fn mul(self, f: f64) -> Self {
        FourVector {
            vector: f * self.vector,
        }
    }
}

impl Div<f64> for FourVector {
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        FourVector {
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
    pub fn new(x0: f64, x1: f64, x2: f64, x3: f64) -> FourVector {
        FourVector {
            vector: Vector4::new(x0, x1, x2, x3),
        }
    }

    pub fn get_as_vector(self) -> Vector4<f64> {
        self.vector
    }
}

#[cfg(test)]
mod tests {
    use crate::four_vector::FourVector;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_multiplication_self() {
        let v1 = FourVector::new(1.0, 2.0, 3.0, 4.0);
        assert_abs_diff_eq!(v1 * v1, -28.0);
    }

    #[test]
    fn test_multiplication_different() {
        let v1 = FourVector::new(1.0, 2.0, 3.0, 4.0);
        let v2 = FourVector::new(5.0, 6.0, 7.0, 8.0);

        assert_abs_diff_eq!(v1 * v2, -60.0);
    }
}

use crate::geometry::point::{CoordinateSystem, Point};
use nalgebra::{Vector3, Vector4};
use std::ops::{Add, Div, Index, Mul, Neg};

#[derive(Clone, Copy, Debug, PartialEq)]
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
        debug_assert_eq!(self.coordinate_system, rhs.coordinate_system);
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

    pub fn new_boyer_lindquist(a: f64, t: f64, r: f64, theta: f64, phi: f64) -> FourVector {
        FourVector {
            coordinate_system: CoordinateSystem::BoyerLindquist { a },
            vector: Vector4::new(t, r, theta, phi),
        }
    }

    pub fn get_as_vector(self) -> Vector4<f64> {
        self.vector
    }

    pub fn get_spatial_vector(self) -> Vector3<f64> {
        Vector3::new(self.vector[1], self.vector[2], self.vector[3])
    }

    pub fn get_z_cartesian(self, at: &Point) -> f64 {
        debug_assert_eq!(self.coordinate_system, at.coordinate_system);
        match self.coordinate_system {
            CoordinateSystem::Cartesian => self.vector[3],
            CoordinateSystem::Spherical | CoordinateSystem::BoyerLindquist { .. } => {
                let r = at.vector[1];
                let theta = at.vector[2];
                let (st, ct) = (theta.sin(), theta.cos());
                ct * self.vector[1] - r * st * self.vector[2]
            }
        }
    }
}

impl Index<usize> for FourVector {
    type Output = f64;

    fn index(&self, index: usize) -> &Self::Output {
        &self.vector[index]
    }
}

#[cfg(test)]
mod tests {
    use crate::geometry::euclidean::EuclideanSpace;
    use crate::geometry::four_vector::FourVector;
    use crate::geometry::geometry::InnerProduct;
    use crate::geometry::point::Point;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_multiplication_self() {
        let geometry = EuclideanSpace::new();
        let v1 = FourVector::new_cartesian(1.0, 2.0, 3.0, 4.0);
        assert_abs_diff_eq!(
            geometry.inner_product(&Point::new_cartesian(0.0, 0.0, 0.0, 0.0), &v1, &v1),
            -28.0
        );
    }

    #[test]
    fn test_multiplication_different() {
        let geometry = EuclideanSpace::new();
        let v1 = FourVector::new_cartesian(1.0, 2.0, 3.0, 4.0);
        let v2 = FourVector::new_cartesian(5.0, 6.0, 7.0, 8.0);

        assert_abs_diff_eq!(
            geometry.inner_product(&Point::new_cartesian(0.0, 0.0, 0.0, 0.0), &v1, &v2),
            -60.0
        );
    }

    /// Finite-difference dz_cartesian/dλ by advancing the *native* coordinates
    /// linearly along the velocity: this is the directional derivative of the
    /// coordinate map z(native), which is exactly what `get_z_cartesian`
    /// computes analytically via the Jacobian z-row.
    fn finite_difference_z(at: &Point, v: &FourVector) -> f64 {
        let eps = 1e-6;
        let cs = at.coordinate_system;
        let plus = Point::new(
            at[0] + eps * v[0],
            at[1] + eps * v[1],
            at[2] + eps * v[2],
            at[3] + eps * v[3],
            cs,
        );
        let minus = Point::new(
            at[0] - eps * v[0],
            at[1] - eps * v[1],
            at[2] - eps * v[2],
            at[3] - eps * v[3],
            cs,
        );
        (plus.get_spatial_vector_cartesian()[2] - minus.get_spatial_vector_cartesian()[2])
            / (2.0 * eps)
    }

    #[test]
    fn test_get_z_cartesian_cartesian_is_the_z_component() {
        // In the Cartesian chart the Jacobian z-row is [0,0,0,1], so ż is just
        // the velocity's z-component and the base point is irrelevant.
        let v = FourVector::new_cartesian(1.0, 2.0, 3.0, 4.0);
        let at = Point::new_cartesian(0.0, 5.0, 6.0, 7.0);
        assert_abs_diff_eq!(v.get_z_cartesian(&at), 4.0);
    }

    #[test]
    fn test_get_z_cartesian_spherical_matches_finite_difference() {
        let at = Point::new_spherical(0.0, 5.0, 1.0, 0.7);
        let v = FourVector::new_spherical(1.0, 0.3, 0.2, 0.1);
        assert_abs_diff_eq!(
            v.get_z_cartesian(&at),
            finite_difference_z(&at, &v),
            epsilon = 1e-6
        );
    }

    #[test]
    fn test_get_z_cartesian_boyer_lindquist_matches_finite_difference() {
        // z = r cosθ in Boyer-Lindquist as well (a enters only x and y), so
        // the same z-row applies and must match the Cartesian-z derivative.
        let at = Point::new_boyer_lindquist(0.8, 0.0, 5.0, 1.0, 0.7);
        let v = FourVector::new_boyer_lindquist(0.8, 1.0, 0.3, 0.2, 0.1);
        assert_abs_diff_eq!(
            v.get_z_cartesian(&at),
            finite_difference_z(&at, &v),
            epsilon = 1e-6
        );
    }

    #[test]
    fn test_get_z_cartesian_is_independent_of_spin() {
        // The z-row [0, cosθ, -r sinθ, 0] carries no `a`, so identical (r, θ)
        // and (v_r, v_θ) give identical ż at any spin.
        let at_a = Point::new_boyer_lindquist(0.8, 0.0, 5.0, 1.0, 0.7);
        let v_a = FourVector::new_boyer_lindquist(0.8, 1.0, 0.3, 0.2, 0.1);
        let at_b = Point::new_boyer_lindquist(0.2, 0.0, 5.0, 1.0, 0.7);
        let v_b = FourVector::new_boyer_lindquist(0.2, 1.0, 0.3, 0.2, 0.1);
        assert_abs_diff_eq!(
            v_a.get_z_cartesian(&at_a),
            v_b.get_z_cartesian(&at_b),
            epsilon = 1e-15
        );
    }

    #[test]
    fn test_get_z_cartesian_equatorial_closed_form() {
        use std::f64::consts::FRAC_PI_2;
        // At the equator θ = π/2: purely radial motion stays in the plane
        // (ż = 0); purely polar motion descends at rate -r·θ̇.
        let at = Point::new_spherical(0.0, 5.0, FRAC_PI_2, 0.0);
        let radial = FourVector::new_spherical(1.0, 2.0, 0.0, 0.0);
        assert_abs_diff_eq!(radial.get_z_cartesian(&at), 0.0, epsilon = 1e-12);
        let polar = FourVector::new_spherical(1.0, 0.0, 0.4, 0.0);
        // ż = -r sinθ v_θ = -5 · 1 · 0.4 = -2.0
        assert_abs_diff_eq!(polar.get_z_cartesian(&at), -2.0, epsilon = 1e-12);
    }
}

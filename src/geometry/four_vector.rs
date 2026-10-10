use crate::geometry::coordinates::atlas::{CartesianConvertible, Chart};
use crate::geometry::coordinates::euclidean::CartesianChart;
use crate::geometry::point::Point;
use nalgebra::{Matrix4, Vector3, Vector4};
use std::marker::PhantomData;
use std::ops::{Add, Div, Index, Mul, Neg};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Up;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Down;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Tensor1<C: Chart, I> {
    pub at: Point<C>,
    pub vector: Vector4<f64>,
    _index: PhantomData<I>,
}

pub type FourVector<C> = Tensor1<C, Up>;
pub type Covector<C> = Tensor1<C, Down>;

impl<C: Chart, I> Neg for Tensor1<C, I> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Tensor1 {
            at: self.at,
            vector: self.vector.neg(),
            _index: PhantomData,
        }
    }
}

impl<C: Chart, I> Add for Tensor1<C, I> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Tensor1 {
            at: self.at,
            vector: self.vector + rhs.vector,
            _index: PhantomData,
        }
    }
}

impl<C: Chart, I> Mul<Tensor1<C, I>> for f64 {
    type Output = Tensor1<C, I>;

    fn mul(self, f: Tensor1<C, I>) -> Tensor1<C, I> {
        Tensor1 {
            at: f.at,
            vector: self * f.vector,
            _index: PhantomData,
        }
    }
}

impl<C: Chart, I> Mul<f64> for Tensor1<C, I> {
    type Output = Self;

    fn mul(self, f: f64) -> Self {
        Tensor1 {
            at: self.at,
            vector: f * self.vector,
            _index: PhantomData,
        }
    }
}

impl<C: Chart, I> Div<f64> for Tensor1<C, I> {
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        Tensor1 {
            at: self.at,
            vector: self.vector / rhs,
            _index: PhantomData,
        }
    }
}

impl<C: Chart, I> Index<usize> for Tensor1<C, I> {
    type Output = f64;

    fn index(&self, index: usize) -> &Self::Output {
        &self.vector[index]
    }
}

impl<C: Chart, I> Tensor1<C, I> {
    pub fn new(x0: f64, x1: f64, x2: f64, x3: f64, at: Point<C>) -> Tensor1<C, I> {
        Tensor1 {
            at,
            vector: Vector4::new(x0, x1, x2, x3),
            _index: PhantomData,
        }
    }

    pub fn get_spatial_vector(&self) -> Vector3<f64> {
        Vector3::new(self.vector[1], self.vector[2], self.vector[3])
    }
}

impl<C: Chart> FourVector<C> {
    pub fn lower(self, g: &Matrix4<f64>) -> Covector<C> {
        Tensor1 {
            at: self.at,
            vector: g * self.vector,
            _index: PhantomData,
        }
    }
}

impl<C: Chart> Covector<C> {
    pub fn raise(self, g_inv: &Matrix4<f64>) -> FourVector<C> {
        Tensor1 {
            at: self.at,
            vector: g_inv * self.vector,
            _index: PhantomData,
        }
    }

    pub fn contract(&self, v: &FourVector<C>) -> f64 {
        self.vector.dot(&v.vector)
    }
}

impl<C: CartesianConvertible> FourVector<C> {
    pub fn to_cartesian(&self) -> FourVector<CartesianChart> {
        C::vector_to_cartesian(self)
    }
}

#[cfg(test)]
mod tests {
    use crate::geometry::coordinates::atlas::CartesianConvertible;
    use crate::geometry::euclidean::EuclideanSpace;
    use crate::geometry::four_vector::FourVector;
    use crate::geometry::geometry::InnerProduct;
    use crate::geometry::point::Point;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_multiplication_self() {
        let geometry = EuclideanSpace::new();
        let v1 = FourVector::new(1.0, 2.0, 3.0, 4.0, Point::new_cartesian(0.0, 0.0, 0.0, 0.0));
        assert_abs_diff_eq!(
            geometry.inner_product(&Point::new_cartesian(0.0, 0.0, 0.0, 0.0), &v1, &v1),
            -28.0
        );
    }

    #[test]
    fn test_multiplication_different() {
        let geometry = EuclideanSpace::new();
        let v1 = FourVector::new(1.0, 2.0, 3.0, 4.0, Point::new_cartesian(0.0, 0.0, 0.0, 0.0));
        let v2 = FourVector::new(5.0, 6.0, 7.0, 8.0, Point::new_cartesian(0.0, 0.0, 0.0, 0.0));

        assert_abs_diff_eq!(
            geometry.inner_product(&Point::new_cartesian(0.0, 0.0, 0.0, 0.0), &v1, &v2),
            -60.0
        );
    }

    /// Finite-difference dz_cartesian/dλ by advancing the *native* coordinates
    /// linearly along the velocity: this is the directional derivative of the
    /// coordinate map z(native), which is exactly what `get_z_cartesian`
    /// computes analytically via the Jacobian z-row.
    fn finite_difference_z<C: CartesianConvertible>(at: &Point<C>, v: &FourVector<C>) -> f64 {
        let eps = 1e-6;
        let cs = at.chart;
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
        let at = Point::new_cartesian(0.0, 5.0, 6.0, 7.0);
        let v = FourVector::new(1.0, 2.0, 3.0, 4.0, at);
        assert_abs_diff_eq!(v.to_cartesian().vector[3], 4.0);
    }

    #[test]
    fn test_get_z_cartesian_spherical_matches_finite_difference() {
        let at = Point::new_spherical(0.0, 5.0, 1.0, 0.7);
        let v = FourVector::new(1.0, 0.3, 0.2, 0.1, at);
        assert_abs_diff_eq!(
            v.to_cartesian().vector[3],
            finite_difference_z(&at, &v),
            epsilon = 1e-6
        );
    }

    #[test]
    fn test_get_z_cartesian_boyer_lindquist_matches_finite_difference() {
        // z = r cosθ in Boyer-Lindquist as well (a enters only x and y), so
        // the same z-row applies and must match the Cartesian-z derivative.
        let at = Point::new_boyer_lindquist(0.8, 1.0, 0.0, 5.0, 1.0, 0.7);
        let v = FourVector::new(1.0, 0.3, 0.2, 0.1, at);
        assert_abs_diff_eq!(
            v.to_cartesian().vector[3],
            finite_difference_z(&at, &v),
            epsilon = 1e-6
        );
    }

    #[test]
    fn test_get_z_cartesian_is_independent_of_spin() {
        // The z-row [0, cosθ, -r sinθ, 0] carries no `a`, so identical (r, θ)
        // and (v_r, v_θ) give identical ż at any spin.
        let at_a = Point::new_boyer_lindquist(0.8, 1.0, 0.0, 5.0, 1.0, 0.7);
        let v_a = FourVector::new(1.0, 0.3, 0.2, 0.1, at_a);
        let at_b = Point::new_boyer_lindquist(0.2, 1.0, 0.0, 5.0, 1.0, 0.7);
        let v_b = FourVector::new(1.0, 0.3, 0.2, 0.1, at_b);
        assert_abs_diff_eq!(
            v_a.to_cartesian().vector[3],
            v_b.to_cartesian().vector[3],
            epsilon = 1e-15
        );
    }

    #[test]
    fn test_get_z_cartesian_equatorial_closed_form() {
        use std::f64::consts::FRAC_PI_2;
        // At the equator θ = π/2: purely radial motion stays in the plane
        // (ż = 0); purely polar motion descends at rate -r·θ̇.
        let at = Point::new_spherical(0.0, 5.0, FRAC_PI_2, 0.0);
        let radial = FourVector::new(1.0, 2.0, 0.0, 0.0, at);
        assert_abs_diff_eq!(radial.to_cartesian().vector[3], 0.0, epsilon = 1e-12);
        let polar = FourVector::new(1.0, 0.0, 0.4, 0.0, at);
        // ż = -r sinθ v_θ = -5 · 1 · 0.4 = -2.0
        assert_abs_diff_eq!(polar.to_cartesian().vector[3], -2.0, epsilon = 1e-12);
    }
}

use crate::four_vector::CoordinateSystem::Spherical;
use crate::four_vector::{CoordinateSystem, FourVector};
use crate::geometry::{Geometry, HasCoordinateSystem, Tetrad};
use crate::runge_kutta::OdeFunction;
use crate::scene::EquationOfMotionState;
use nalgebra::{Const, Matrix4, OVector, Vector4};

#[derive(Clone)]
pub struct Schwarzschild {
    radius: f64,
}

impl Schwarzschild {
    pub fn new(radius: f64) -> Self {
        Schwarzschild { radius }
    }
}

impl OdeFunction<Const<8>> for Schwarzschild {
    // TODO: maybe just have geodesic being used in solver. This doesn't need to be that generic here.
    fn apply(&self, t: f64, y: &OVector<f64, Const<8>>) -> OVector<f64, Const<8>> {
        self.geodesic(t, y)
    }
}

impl HasCoordinateSystem for Schwarzschild {
    fn coordinate_system(&self) -> CoordinateSystem {
        Spherical
    }
}

// All coordinates here are spherical coordinates.
impl Geometry for Schwarzschild {
    fn geodesic(&self, _: f64, y: &EquationOfMotionState) -> EquationOfMotionState {
        let _t = y[0];
        let r = y[1];
        let theta = y[2];
        let _phi = y[3];

        let v_t = y[4];
        let v_r = y[5];
        let v_theta = y[6];
        let v_phi = y[7];

        let a = 1.0 - self.radius / r;
        let a_prime = self.radius / r * r;
        let a_div_aprime = self.radius / (r * r - self.radius * r);

        // acceleration
        let a_t = -(a_div_aprime) * v_t * v_r;
        let a_r = -0.5 * a * a_prime * v_t * v_t
            + 0.5 * (a_div_aprime) * v_r * v_r
            + a * r * (v_theta * v_theta + v_phi * v_phi * theta.sin() * theta.sin());
        let a_theta = -(2.0 / r) * v_r * v_theta + theta.sin() * theta.cos() * v_phi * v_phi;
        let a_phi = -(2.0 / r) * v_phi * v_r - 2.0 * theta.cos() / theta.sin() * v_theta * v_phi;

        let y_new = EquationOfMotionState::from_column_slice(&[
            v_t, v_r, v_theta, v_phi, a_t, a_r, a_theta, a_phi,
        ]);
        y_new
    }

    // TODO: take into account Lorentz transformations.
    // TODO: take into account rotations.
    fn get_tetrad_at(&self, position: &Vector4<f64>) -> Tetrad {
        let r = position[1];
        let theta = position[2];
        let _phi = position[3];

        let rr0 = self.radius / r;

        Tetrad::new(
            position.clone(),
            FourVector::new_spherical(1.0 / (1.0 - rr0), -rr0.sqrt(), 0.0, 0.0),
            FourVector::new_spherical(0.0, 0.0, 0.0, 1.0 / (r * theta.sin())), // Phi
            -FourVector::new_spherical(0.0, 0.0, 1.0 / r, 0.0),                // Theta
            -FourVector::new_spherical(-rr0.sqrt() / (1.0 - rr0), 1.0, 0.0, 0.0), // R
        )
    }

    fn lorentz_transformation(
        &self,
        _position: &Vector4<f64>,
        _velocity: &FourVector,
    ) -> Matrix4<f64> {
        let mut matrix = Matrix4::zeros();
        matrix[(0, 0)] = 1.0;
        matrix[(1, 1)] = 1.0;
        matrix[(2, 2)] = 1.0;
        matrix[(3, 3)] = 1.0;

        matrix
    }
}

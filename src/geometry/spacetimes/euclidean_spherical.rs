use crate::geometry::circular_orbit::OrbitKillingDecomposition;
use crate::geometry::coordinates::spherical::SphericalChart;
use crate::geometry::four_vector::FourVector;
use crate::geometry::geometry::{
    GeodesicSolver, Geometry, HasChart, InnerProduct, Signature, SupportQuantities,
};
use crate::geometry::point::Point;
use crate::geometry::tetrad::Tetrad;
use crate::rendering::error::RaytracerError;
use crate::rendering::integration::runge_kutta::OdeFunction;
use crate::rendering::integration::state::EquationOfMotionState;
use crate::rendering::light::temperature::{ConstantTemperatureComputer, TemperatureComputer};
use nalgebra::{Const, Matrix4, OVector};

pub struct EuclideanSpaceSphericalSolver {
    chart: SphericalChart,
}

#[derive(Clone)]
pub struct EuclideanSpaceSpherical {}

impl Default for EuclideanSpaceSpherical {
    fn default() -> Self {
        Self::new()
    }
}

impl EuclideanSpaceSpherical {
    pub fn new() -> Self {
        EuclideanSpaceSpherical {}
    }
}

impl OdeFunction<Const<8>> for EuclideanSpaceSphericalSolver {
    // TODO: maybe just have geodesic being used in solver. This doesn't need to be that generic here.
    fn apply(&self, t: f64, y: &OVector<f64, Const<8>>) -> OVector<f64, Const<8>> {
        self.geodesic(t, y)
    }
}

impl GeodesicSolver for EuclideanSpaceSphericalSolver {
    type Chart = SphericalChart;

    fn chart(&self) -> Self::Chart {
        self.chart
    }

    fn geodesic(&self, _: f64, y: &EquationOfMotionState) -> EquationOfMotionState {
        let _t = y[0];
        let r = y[1];
        let theta = y[2];
        let _phi = y[3];

        let v_t = y[4];
        let v_r = y[5];
        let v_theta = y[6];
        let v_phi = y[7];

        // acceleration
        let a_t = 0.0;
        let a_r = r * (v_theta * v_theta + v_phi * v_phi * theta.sin() * theta.sin());
        let a_theta = -(2.0 / r) * v_r * v_theta + theta.sin() * theta.cos() * v_phi * v_phi;
        let a_phi = -(2.0 / r) * v_phi * v_r - 2.0 * theta.cos() / theta.sin() * v_theta * v_phi;

        // y'

        EquationOfMotionState::from_column_slice(&[
            v_t, v_r, v_theta, v_phi, a_t, a_r, a_theta, a_phi,
        ])
    }
}

impl InnerProduct for EuclideanSpaceSpherical {
    fn inner_product(
        &self,
        position: &Point<SphericalChart>,
        v: &FourVector<SphericalChart>,
        w: &FourVector<SphericalChart>,
    ) -> f64 {
        let r = position[1];
        let theta = position[2];
        let _phi = position[3];

        1.0 * v.vector[0] * w.vector[0]
            - v.vector[1] * w.vector[1]
            - r * r * v.vector[2] * w.vector[2]
            - r * r * theta.sin() * theta.sin() * v.vector[3] * w.vector[3]
    }
}

impl Signature for EuclideanSpaceSpherical {
    fn signature(&self) -> [f64; 4] {
        [1.0, -1.0, -1.0, -1.0]
    }
}

impl HasChart for EuclideanSpaceSpherical {
    type Chart = SphericalChart;
}

impl Geometry for EuclideanSpaceSpherical {
    // TODO: take into account rotations.
    fn get_tetrad_at(&self, position: &Point<SphericalChart>) -> Tetrad<Self::Chart> {
        let r = position[1];
        let theta = position[2];
        let _phi = position[3];

        Tetrad::new(
            *position,
            FourVector::new(1.0, 0.0, 0.0, 0.0, *position),
            FourVector::new(0.0, 0.0, 0.0, 1.0 / (r * theta.sin()), *position), // Phi
            -FourVector::new(0.0, 0.0, 1.0 / r, 0.0, *position),                // Theta
            -FourVector::new(0.0, 1.0, 0.0, 0.0, *position),                    // R
        )
    }

    fn lorentz_transformation(
        &self,
        _position: &Point<SphericalChart>,
        _velocity: &FourVector<SphericalChart>,
    ) -> Matrix4<f64> {
        let mut matrix = Matrix4::zeros();
        matrix[(0, 0)] = 1.0;
        matrix[(1, 1)] = 1.0;
        matrix[(2, 2)] = 1.0;
        matrix[(3, 3)] = 1.0;

        matrix
    }

    fn inside_horizon(&self, _position: &Point<SphericalChart>) -> bool {
        false
    }

    fn closed_orbit(
        &self,
        _position: &Point<SphericalChart>,
        _step_index: usize,
        _max_steps: usize,
    ) -> bool {
        false
    }

    fn get_geodesic_solver(
        &self,
        momentum: &FourVector<Self::Chart>,
    ) -> Box<dyn GeodesicSolver<Chart = Self::Chart>> {
        Box::new(EuclideanSpaceSphericalSolver {
            chart: momentum.at.chart,
        })
    }

    fn get_radial_coordinate(&self, position: &Point<SphericalChart>) -> f64 {
        position[1]
    }

    fn get_constants_of_motion(
        &self,
        position: &Point<SphericalChart>,
        momentum: &FourVector<SphericalChart>,
    ) -> crate::geometry::geometry::ConstantsOfMotion {
        let mut constants = crate::geometry::geometry::ConstantsOfMotion::default();

        let r = position[1];
        let theta = position[2];

        // E = p_t = v^t
        let e = momentum.vector[0];
        constants.push("E", e);

        // L_z = p_phi = g_phiphi v^phi = - r^2 sin^2(theta) v^phi
        let l_z = -r * r * theta.sin() * theta.sin() * momentum.vector[3];
        constants.push("L_z", l_z);

        constants
    }
}

impl SupportQuantities for EuclideanSpaceSpherical {
    fn get_stationary_velocity_at(
        &self,
        _position: &Point<SphericalChart>,
    ) -> FourVector<SphericalChart> {
        FourVector::new(1.0, 0.0, 0.0, 0.0, *_position)
    }

    fn get_zamo_velocity_at(&self, position: &Point<SphericalChart>) -> FourVector<SphericalChart> {
        // No frame dragging: the ZAMO is the static observer.
        self.get_stationary_velocity_at(position)
    }

    fn get_circular_orbit_velocity_at(
        &self,
        _position: &Point<SphericalChart>,
    ) -> Result<FourVector<SphericalChart>, RaytracerError> {
        Ok(FourVector::new(1.0, 0.0, 0.0, 0.0, *_position))
        // TODO: implement proper circular orbit velocity in Euclidean space.
    }

    fn axial_killing_vector(&self, position: &Point<SphericalChart>) -> FourVector<SphericalChart> {
        SphericalChart::physical_axial_vector(position)
    }

    fn circular_orbit_killing_coefficients(
        &self,
        _position: &Point<SphericalChart>,
    ) -> Result<OrbitKillingDecomposition, RaytracerError> {
        // Emitters are at rest in flat space (matches
        // get_circular_orbit_velocity_at above).
        Ok(OrbitKillingDecomposition {
            u_t: 1.0,
            u_phi: 0.0,
        })
    }

    fn get_temperature_computer(
        &self,
        temperature: f64,
        _inner_radius: f64,
        _outer_radius: f64,
    ) -> Result<Box<dyn TemperatureComputer>, RaytracerError> {
        Ok(Box::new(ConstantTemperatureComputer::new(temperature)))
    }
}

#[cfg(test)]
mod tests {
    use crate::geometry::spacetimes::euclidean_spherical::EuclideanSpaceSpherical;
    use crate::geometry::four_vector::FourVector;
    use crate::geometry::point::Point;
    use crate::geometry::spherical_coordinates_helper::cartesian_to_spherical;
    use crate::rendering::camera::Camera;
    use crate::rendering::camera::debug::save_rays_to_file;
    use std::f64::consts::PI;

    #[ignore]
    #[test]
    fn save_camera_rays() {
        let rows = 30;
        let cols = 30;

        let position = cartesian_to_spherical(&Point::new_cartesian(0.0, 0.0, 0.0, -10.0));
        let velocity = FourVector::new(1.0, 0.0, 0.0, 0.0, position);
        let camera = Camera::new(
            position,
            velocity,
            PI / 2.0,
            rows,
            cols,
            0.0,
            0.0,
            0.0,
            &EuclideanSpaceSpherical::new(),
        )
        .unwrap();

        save_rays_to_file(
            rows,
            cols,
            &position,
            EuclideanSpaceSpherical::new(),
            camera,
        );
    }
}

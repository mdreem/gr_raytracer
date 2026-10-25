use crate::geometry::circular_orbit::OrbitKillingDecomposition;
use crate::geometry::coordinates::stereographic::StereographicChart;
use crate::geometry::four_vector::FourVector;
use crate::geometry::geometry::{
    ConstantsOfMotion, GeodesicSolver, Geometry, HasChart, InnerProduct, Signature,
    SupportQuantities,
};
use crate::geometry::point::Point;
use crate::geometry::tetrad::Tetrad;
use crate::rendering::error::RaytracerError;
use crate::rendering::integration::runge_kutta::OdeFunction;
use crate::rendering::integration::state::EquationOfMotionState;
use crate::rendering::light::temperature::{ConstantTemperatureComputer, TemperatureComputer};
use nalgebra::{Const, Matrix4, OVector, Vector3};
use std::f64::consts::PI;

/// The static universe R x S^3 in dimensionless stereographic coordinates u,
/// ds^2 = dt^2 - R^2 Omega^2 |du|^2 with Omega = 2 / (1 + |u|^2).
///
/// The radius only scales the spatial metric by a constant, so it drops out of the geodesic
/// equation and appears only where lengths are measured: inner product, tetrad, distances.
#[derive(Clone, Debug)]
pub struct HyperSphere {
    pub radius: f64,
}

impl HyperSphere {
    pub fn new(radius: f64) -> Self {
        Self { radius }
    }
}

/// Conformal factor Omega of the unit 3-sphere in stereographic coordinates.
fn conformal_factor(position: &Point<StereographicChart>) -> f64 {
    let u_sq = position[1].powi(2) + position[2].powi(2) + position[3].powi(2);
    2.0 / (1.0 + u_sq)
}

struct HyperSphereSolver {
    chart: StereographicChart,
}

impl OdeFunction<Const<8>> for HyperSphereSolver {
    // TODO: maybe just have geodesic being used in solver. This doesn't need to be that generic here.
    fn apply(&self, t: f64, y: &OVector<f64, Const<8>>) -> OVector<f64, Const<8>> {
        self.geodesic(t, y)
    }
}

impl GeodesicSolver for HyperSphereSolver {
    type Chart = StereographicChart;

    fn chart(&self) -> Self::Chart {
        self.chart
    }

    fn geodesic(&self, _: f64, y_state: &EquationOfMotionState) -> EquationOfMotionState {
        let _t = y_state[0];
        let x = y_state[1];
        let y = y_state[2];
        let z = y_state[3];

        let p_t = y_state[4];
        let p_x = y_state[5];
        let p_y = y_state[6];
        let p_z = y_state[7];

        let x_vec = Vector3::new(x, y, z);
        let v_vec = Vector3::new(p_x, p_y, p_z);

        let x_v = x_vec.dot(&v_vec);
        let x_sq = x_vec.dot(&x_vec);
        let v_sq = v_vec.dot(&v_vec);
        let denom = 1.0 + x_sq;

        let new_v_vec = (4.0 * x_v * v_vec - 2.0 * v_sq * x_vec) / denom;

        let dt = p_t;
        let dx = p_x;
        let dy = p_y;
        let dz = p_z;

        let a_t = 0.0;
        let a_x = new_v_vec[0];
        let a_y = new_v_vec[1];
        let a_z = new_v_vec[2];

        EquationOfMotionState::from_column_slice(&[dt, dx, dy, dz, a_t, a_x, a_y, a_z])
    }
}

impl HasChart for HyperSphere {
    type Chart = StereographicChart;
}

impl InnerProduct for HyperSphere {
    fn inner_product(
        &self,
        position: &Point<StereographicChart>,
        v: &FourVector<StereographicChart>,
        w: &FourVector<StereographicChart>,
    ) -> f64 {
        let scale = (self.radius * conformal_factor(position)).powi(2);
        v[0] * w[0] - scale * (v[1] * w[1] + v[2] * w[2] + v[3] * w[3])
    }
}

impl Signature for HyperSphere {
    fn signature(&self) -> [f64; 4] {
        [1.0, -1.0, -1.0, -1.0]
    }
}

impl Geometry for HyperSphere {
    /// The coordinate directions normalized by the metric. The south chart is the inversion of
    /// the north one, so the same construction there has the opposite handedness in the
    /// Cartesian view (see `tetrad_handedness_fallback`).
    fn get_tetrad_at(&self, position: &Point<StereographicChart>) -> Tetrad<StereographicChart> {
        let e = 1.0 / (self.radius * conformal_factor(position));
        Tetrad::new(
            *position,
            FourVector::new(1.0, 0.0, 0.0, 0.0, *position),
            FourVector::new(0.0, e, 0.0, 0.0, *position),
            FourVector::new(0.0, 0.0, e, 0.0, *position),
            FourVector::new(0.0, 0.0, 0.0, e, *position),
        )
    }

    /// Only the static observer is supported.
    fn lorentz_transformation(
        &self,
        _position: &Point<StereographicChart>,
        _velocity: &FourVector<StereographicChart>,
    ) -> Matrix4<f64> {
        Matrix4::identity()
    }

    fn inside_horizon(&self, _position: &Point<StereographicChart>) -> bool {
        false
    }

    /// Light travels along great circles, so a ray that hits nothing within one full circuit
    /// never will. Camera rays have dt/dlambda = -1 and unit proper spatial speed, so -t is the
    /// distance traveled.
    fn closed_orbit(
        &self,
        position: &Point<StereographicChart>,
        step_index: usize,
        max_steps: usize,
    ) -> bool {
        -position[0] >= 2.0 * PI * self.radius || step_index == max_steps - 1
    }

    fn get_geodesic_solver(
        &self,
        momentum: &FourVector<Self::Chart>,
    ) -> Box<dyn GeodesicSolver<Chart = Self::Chart>> {
        Box::new(HyperSphereSolver {
            chart: momentum.at.chart,
        })
    }

    /// Geodesic distance from the north pole, the origin of the north chart:
    /// chi = 2 atan |u| there, and pi - 2 atan |u| from the south chart.
    fn get_radial_coordinate(&self, position: &Point<StereographicChart>) -> f64 {
        let u = (position[1].powi(2) + position[2].powi(2) + position[3].powi(2)).sqrt();
        let chi = 2.0 * u.atan();
        match position.chart {
            StereographicChart::WithNorthPole => self.radius * chi,
            StereographicChart::WithSouthPole => self.radius * (PI - chi),
        }
    }

    fn get_constants_of_motion(
        &self,
        _position: &Point<StereographicChart>,
        momentum: &FourVector<StereographicChart>,
    ) -> ConstantsOfMotion {
        let mut constants = ConstantsOfMotion::default();
        constants.push("E", momentum[0]);
        constants
    }
}

impl SupportQuantities for HyperSphere {
    fn get_stationary_velocity_at(
        &self,
        position: &Point<StereographicChart>,
    ) -> FourVector<StereographicChart> {
        FourVector::new(1.0, 0.0, 0.0, 0.0, *position)
    }

    fn get_zamo_velocity_at(
        &self,
        position: &Point<StereographicChart>,
    ) -> FourVector<StereographicChart> {
        // No frame dragging: the ZAMO is the static observer.
        self.get_stationary_velocity_at(position)
    }

    fn get_circular_orbit_velocity_at(
        &self,
        position: &Point<StereographicChart>,
    ) -> Result<FourVector<StereographicChart>, RaytracerError> {
        Ok(self.get_stationary_velocity_at(position))
    }

    /// Rotation about the u_3 axis. Rotations commute with the inversion, so it has the same
    /// form in both charts.
    fn axial_killing_vector(
        &self,
        position: &Point<StereographicChart>,
    ) -> FourVector<StereographicChart> {
        FourVector::new(0.0, -position[2], position[1], 0.0, *position)
    }

    fn circular_orbit_killing_coefficients(
        &self,
        _position: &Point<StereographicChart>,
    ) -> Result<OrbitKillingDecomposition, RaytracerError> {
        // Emitters are at rest (matches get_circular_orbit_velocity_at above).
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
    use super::*;
    use crate::geometry::coordinates::atlas::Atlas;
    use crate::rendering::integration::integrator::{IntegrationConfiguration, Integrator};
    use crate::rendering::ray::Ray;
    use nalgebra::Vector4;

    const RADIUS: f64 = 2.0;
    const NORTH: StereographicChart = StereographicChart::WithNorthPole;
    const SOUTH: StereographicChart = StereographicChart::WithSouthPole;

    /// Unit embedding in R^4, (X_0, X_spatial), with the north pole at X_0 = 1.
    fn embedding(p: &Point<StereographicChart>) -> Vector4<f64> {
        let u = Vector3::new(p[1], p[2], p[3]);
        let d = 1.0 + u.norm_squared();
        let x0 = (1.0 - u.norm_squared()) / d;
        let x0 = if p.chart == NORTH { x0 } else { -x0 };
        Vector4::new(x0, 2.0 * u[0] / d, 2.0 * u[1] / d, 2.0 * u[2] / d)
    }

    /// A camera-like ray in the north chart: p^t = -1 and unit proper spatial speed along
    /// `direction`, together with the great circle it has to follow.
    fn ray(
        u: Vector3<f64>,
        direction: Vector3<f64>,
    ) -> (Ray<StereographicChart>, Vector4<f64>, Vector4<f64>) {
        let geometry = HyperSphere::new(RADIUS);
        let p = Point::new(0.0, u[0], u[1], u[2], NORTH);
        let v = direction.normalize() / (RADIUS * conformal_factor(&p));
        let momentum = FourVector::new(-1.0, v[0], v[1], v[2], p);
        assert!(geometry.inner_product(&p, &momentum, &momentum).abs() < 1e-12);

        let h = 1e-6;
        let ahead = Point::new(
            0.0,
            u[0] + h * v[0],
            u[1] + h * v[1],
            u[2] + h * v[2],
            NORTH,
        );
        let behind = Point::new(
            0.0,
            u[0] - h * v[0],
            u[1] - h * v[1],
            u[2] - h * v[2],
            NORTH,
        );
        let tangent = (embedding(&ahead) - embedding(&behind)).normalize();
        (Ray::new(0, 0, p, momentum), embedding(&p), tangent)
    }

    #[test]
    fn tetrad_is_orthonormal_in_both_charts() {
        let geometry = HyperSphere::new(RADIUS);
        for p in [
            Point::new(0.0, 0.3, 0.1, -0.2, NORTH),
            Point::new(0.0, -1.5, 0.4, 0.9, SOUTH),
        ] {
            let tetrad = geometry.get_tetrad_at(&p);
            let e = [tetrad.t, tetrad.x, tetrad.y, tetrad.z];
            for i in 0..4 {
                for j in 0..4 {
                    let expected = match (i, j) {
                        (0, 0) => 1.0,
                        _ if i == j => -1.0,
                        _ => 0.0,
                    };
                    let g = geometry.inner_product(&p, &e[i], &e[j]);
                    assert!((g - expected).abs() < 1e-12, "g({i}, {j}) = {g}");
                }
            }
        }
    }

    /// Omega^2 |v|^2 is the squared proper speed, conserved along an affinely parametrised
    /// geodesic.
    #[test]
    fn proper_speed_is_conserved_by_the_equation_of_motion() {
        let (ray, _, _) = ray(Vector3::new(0.3, 0.1, -0.2), Vector3::new(0.2, 0.5, 0.1));
        let solver = HyperSphere::new(RADIUS).get_geodesic_solver(&ray.momentum);
        let y = solver.create_initial_state(&ray.momentum);
        let dy = solver.geodesic(0.0, &y);
        // d/dlambda (Omega^2 |v|^2) = 2 Omega^2 (v . a + |v|^2 v . grad ln Omega)
        let u = Vector3::new(y[1], y[2], y[3]);
        let v = Vector3::new(y[5], y[6], y[7]);
        let a = Vector3::new(dy[5], dy[6], dy[7]);
        let grad_ln_omega = -2.0 * u / (1.0 + u.norm_squared());
        let derivative = v.dot(&a) + v.norm_squared() * v.dot(&grad_ln_omega);
        assert!(derivative.abs() < 1e-14, "{derivative}");
    }

    /// A full circuit around a great circle that leaves the north chart: the integrator switches
    /// to the south chart and back, and the ray stays on the exact great circle and null.
    #[test]
    fn rays_follow_great_circles_across_chart_switches() {
        let geometry = HyperSphere::new(RADIUS);
        let (ray, start, tangent) = ray(Vector3::new(0.3, 0.1, -0.2), Vector3::new(0.2, 0.5, 0.1));
        let integrator = Integrator::new(
            &geometry,
            IntegrationConfiguration::new(100_000, 1e9, 0.01, 1e-10),
        );
        let (steps, _) = integrator.integrate(&ray).expect("integration");

        let switches = steps
            .steps
            .windows(2)
            .filter(|pair| pair[0].x.chart != pair[1].x.chart)
            .count();
        assert_eq!(switches, 2);
        assert!(steps.steps.iter().any(|step| step.x.chart == SOUTH));

        let mut max_error = 0.0f64;
        for step in &steps.steps {
            let angle = -step.x[0] / RADIUS;
            let expected = angle.cos() * start + angle.sin() * tangent;
            max_error = max_error.max((embedding(&step.x) - expected).norm());
            let k = geometry.inner_product(&step.x, &step.p, &step.p);
            assert!(k.abs() < 1e-6, "null condition drifted: {k}");
        }
        assert!(-steps.steps.last().unwrap().x[0] >= 2.0 * PI * RADIUS);
        assert!(
            max_error < 1e-6,
            "deviation from the great circle: {max_error}"
        );
    }

    #[test]
    fn radial_coordinate_is_the_distance_from_the_north_pole_in_both_charts() {
        let geometry = HyperSphere::new(RADIUS);
        let p = Point::new(0.0, 1.3, -0.4, 2.1, NORTH);
        let q = NORTH.transition_point(&p, SOUTH);
        let chi = embedding(&p)[0].acos();
        assert!((geometry.get_radial_coordinate(&p) - RADIUS * chi).abs() < 1e-12);
        assert!((geometry.get_radial_coordinate(&q) - RADIUS * chi).abs() < 1e-12);
    }
}

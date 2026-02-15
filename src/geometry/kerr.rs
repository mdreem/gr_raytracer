use crate::geometry::four_vector::FourVector;
use crate::geometry::geometry::{
    GeodesicSolver, Geometry, HasCoordinateSystem, InnerProduct, Signature, SupportQuantities,
};
use crate::geometry::gram_schmidt::gram_schmidt;
use crate::geometry::point::CoordinateSystem::Cartesian;
use crate::geometry::point::{CoordinateSystem, Point};
use crate::geometry::tetrad::Tetrad;
use crate::rendering::ray::Ray;
use crate::rendering::raytracer::RaytracerError;
use crate::rendering::runge_kutta::OdeFunction;
use crate::rendering::scene::EquationOfMotionState;
use crate::rendering::temperature::{KerrTemperatureComputer, TemperatureComputer};
use log::{debug, error, trace};
use nalgebra::{Const, Matrix4, OVector, Vector3, Vector4};

#[derive(Clone, Debug)]
pub struct Kerr {
    pub radius: f64,
    pub a: f64,
    horizon_epsilon: f64,
}

struct KerrSolver {
    radius: f64,
    a: f64,
}

fn compute_r_sqr(a: f64, x: f64, y: f64, z: f64) -> f64 {
    let rho_sqr = x * x + y * y + z * z;
    0.5 * (rho_sqr - a * a + ((rho_sqr - a * a) * (rho_sqr - a * a) + 4.0 * a * a * z * z).sqrt())
}

fn k_vector(a: f64, x: f64, y: f64, z: f64) -> Vector4<f64> {
    let r_sqr = compute_r_sqr(a, x, y, z);
    let r = r_sqr.sqrt();

    let k_0 = 1.0;
    let k_x = (r * x + a * y) / (r_sqr + a * a);
    let k_y = (r * y - a * x) / (r_sqr + a * a);
    let k_z = z / r;

    Vector4::new(k_0, k_x, k_y, k_z)
}

// Ignore negative r for now.
fn metric(radius: f64, a: f64, x: f64, y: f64, z: f64) -> Matrix4<f64> {
    let r_sqr = compute_r_sqr(a, x, y, z);
    let r = r_sqr.sqrt();
    let f = (r * r * r * radius) / (r * r * r * r + a * a * z * z);

    let k = k_vector(a, x, y, z);
    let k_0 = k[0];
    let k_x = k[1];
    let k_y = k[2];
    let k_z = k[3];

    let mut metric = Matrix4::zeros();
    metric[(0, 0)] = k_0 * k_0 * f - 1.0;
    metric[(0, 1)] = k_0 * k_x * f;
    metric[(0, 2)] = k_0 * k_y * f;
    metric[(0, 3)] = k_0 * k_z * f;

    metric[(1, 0)] = metric[(0, 1)];
    metric[(1, 1)] = k_x * k_x * f + 1.0;
    metric[(1, 2)] = k_x * k_y * f;
    metric[(1, 3)] = k_x * k_z * f;

    metric[(2, 0)] = metric[(0, 2)];
    metric[(2, 1)] = metric[(1, 2)];
    metric[(2, 2)] = k_y * k_y * f + 1.0;
    metric[(2, 3)] = k_y * k_z * f;

    metric[(3, 0)] = metric[(0, 3)];
    metric[(3, 1)] = metric[(1, 3)];
    metric[(3, 2)] = metric[(2, 3)];
    metric[(3, 3)] = k_z * k_z * f + 1.0;

    trace!("metric: {:?}", metric);

    metric
}

// For Kerr-Schild form g_{mu nu} = eta_{mu nu} + f k_mu k_nu, the inverse is
// g^{mu nu} = eta^{mu nu} - f k^mu k^nu.
fn metric_contravariant(radius: f64, a: f64, x: f64, y: f64, z: f64) -> Matrix4<f64> {
    let r_sqr = compute_r_sqr(a, x, y, z);
    let r = r_sqr.sqrt();
    let f = (r * r * r * radius) / (r * r * r * r + a * a * z * z);

    let k = k_vector(a, x, y, z);
    let ku = [-k[0], k[1], k[2], k[3]];

    let mut metric = Matrix4::zeros();
    metric[(0, 0)] = -1.0;
    metric[(1, 1)] = 1.0;
    metric[(2, 2)] = 1.0;
    metric[(3, 3)] = 1.0;

    for i in 0..4 {
        for j in 0..4 {
            metric[(i, j)] -= f * ku[i] * ku[j];
        }
    }

    metric
}

impl Kerr {
    pub fn new(radius: f64, a: f64, horizon_epsilon: f64) -> Self {
        Kerr {
            radius,
            a,
            horizon_epsilon,
        }
    }

    fn jacobian_spherical_to_cartesian(&self, position: &Point) -> Matrix4<f64> {
        let spherical = position.get_as_spherical();

        let r = spherical[0];
        let theta = spherical[1];
        let phi = spherical[2];

        let (st, ct) = (theta.sin(), theta.cos());
        let (sp, cp) = (phi.sin(), phi.cos());

        let data = vec![
            // row 0: t = t
            1.0,
            0.0,
            0.0,
            0.0,
            // row 1: x(r,theta,phi)
            0.0,
            st * cp,
            r * ct * cp,
            -r * st * sp,
            // row 2: y(r,theta,phi)
            0.0,
            st * sp,
            r * ct * sp,
            r * st * cp,
            // row 3: z(r,theta,phi)
            0.0,
            ct,
            -r * st,
            0.0,
        ];

        Matrix4::from_row_slice(&data)
    }

    fn ut_contra(&self, position: &Point) -> Result<f64, RaytracerError> {
        let r = position.get_as_spherical()[0];
        let a = self.a;
        let r_s = self.radius;
        let omega = self.angular_velocity(r);

        let g_tt = -(1.0 - r_s / r);
        let g_tphi = -a * r_s / r;
        let g_phiphi = r * r + a * a + a * a * r_s / r;

        let ut_pre = g_tt + 2.0 * omega * g_tphi + omega * omega * g_phiphi;

        if ut_pre >= 0.0 {
            error!(
                "No timelike circular orbit at r = {} (ut_pre = {})",
                r, ut_pre
            );
            return Err(RaytracerError::NoCircularOrbitPossible);
        }

        Ok((-ut_pre).sqrt().recip())
    }

    // https://arxiv.org/abs/1104.5499 equation (36)
    fn angular_velocity(&self, r: f64) -> f64 {
        let a = self.a;
        let r_s = self.radius;
        let m = 0.5 * r_s;
        let sqrt_m = m.sqrt();

        sqrt_m / (r.powf(1.5) + a * sqrt_m)
    }
}

impl OdeFunction<Const<8>> for KerrSolver {
    // TODO: maybe just have geodesic being used in solver. This doesn't need to be that generic here.
    fn apply(&self, t: f64, y: &OVector<f64, Const<8>>) -> OVector<f64, Const<8>> {
        self.geodesic(t, y)
    }
}

impl KerrSolver {
    /// Derivative of the contravariant metric with respect to the coordinate at `index`.
    fn d_matrix_contravariant(
        &self,
        index: usize,
        x: f64,
        y: f64,
        z: f64,
        contravariant_metric: &Matrix4<f64>,
    ) -> Matrix4<f64> {
        let d_contra = self.d_matrix_covariant(index, x, y, z);
        -(contravariant_metric * d_contra * contravariant_metric)
    }

    /// Numerical derivative of the covariant metric with respect to the coordinate at `index`.
    fn d_matrix_covariant(&self, index: usize, x: f64, y: f64, z: f64) -> Matrix4<f64> {
        if index == 0 {
            return Matrix4::zeros();
        }
        let base = 1e-10; // epsilon = 1e-12
        let h = base
            * match index {
                1 => x.abs().max(1.0),
                2 => y.abs().max(1.0),
                3 => z.abs().max(1.0),
                _ => 1.0,
            };

        let (dx, dy, dz) = match index {
            1 => (h, 0.0, 0.0),
            2 => (0.0, h, 0.0),
            3 => (0.0, 0.0, h),
            _ => (0.0, 0.0, 0.0),
        };

        let m_plus = metric(self.radius, self.a, x + dx, y + dy, z + dz);
        let m_minus = metric(self.radius, self.a, x - dx, y - dy, z - dz);

        (m_plus - m_minus) / (2.0 * h)
    }
}

impl HasCoordinateSystem for KerrSolver {
    fn coordinate_system(&self) -> CoordinateSystem {
        Cartesian
    }
}

impl GeodesicSolver for KerrSolver {
    /// Hamiltonian geodesics
    /// H(x, p) = 0.5 * g^{\mu\nu} p_\mu p_\nu
    /// d x^\mu / d\lambda = \partial H / \partial p_\mu = g^{\mu\nu} p_\nu
    /// d p_\mu / d\lambda = - \partial H / \partial x^\mu = -0.5 * \partial g^{\alpha\beta} / \partial x^\mu * p_\alpha p_\beta
    fn geodesic(&self, _: f64, y_state: &EquationOfMotionState) -> EquationOfMotionState {
        let _t = y_state[0];
        let x = y_state[1];
        let y = y_state[2];
        let z = y_state[3];

        let p_t = y_state[4];
        let p_x = y_state[5];
        let p_y = y_state[6];
        let p_z = y_state[7];

        trace!("y_state = {:?}", y_state);

        let p = OVector::<f64, Const<4>>::from_column_slice(&[p_t, p_x, p_y, p_z]);

        let contravariant_metric = metric_contravariant(self.radius, self.a, x, y, z);
        trace!("contravariant_metric = {:?}", contravariant_metric);

        debug_assert!(!p_t.is_nan());

        // Compute dot{x}^\mu.
        let xdot = contravariant_metric * p;
        let dt = xdot[0];
        let dx = xdot[1];
        let dy = xdot[2];
        let dz = xdot[3];

        // Compute accelerations.
        let d_gcontra_dx = self.d_matrix_contravariant(1, x, y, z, &contravariant_metric);
        let d_gcontra_dy = self.d_matrix_contravariant(2, x, y, z, &contravariant_metric);
        let d_gcontra_dz = self.d_matrix_contravariant(3, x, y, z, &contravariant_metric);

        let a_t = 0.0;
        let a_x = -0.5 * (p.transpose() * d_gcontra_dx * p)[(0, 0)];
        let a_y = -0.5 * (p.transpose() * d_gcontra_dy * p)[(0, 0)];
        let a_z = -0.5 * (p.transpose() * d_gcontra_dz * p)[(0, 0)];

        let hamiltonian = 0.5 * (p.transpose() * contravariant_metric * p)[(0, 0)];
        trace!("Hamiltonian H = {}", hamiltonian);

        EquationOfMotionState::from_column_slice(&[dt, dx, dy, dz, a_t, a_x, a_y, a_z])
    }

    fn create_initial_state(&self, ray: &Ray) -> EquationOfMotionState {
        let (x, y, z) = (ray.position[1], ray.position[2], ray.position[3]);
        let covariant_metric = metric(self.radius, self.a, x, y, z);
        trace!("covariant_metric = {:?}", covariant_metric);

        let momentum_covariant = covariant_metric * ray.momentum.vector;

        EquationOfMotionState::from_column_slice(&[
            ray.position[0],
            ray.position[1],
            ray.position[2],
            ray.position[3],
            momentum_covariant[0],
            momentum_covariant[1],
            momentum_covariant[2],
            momentum_covariant[3],
        ])
    }

    fn momentum_from_state(&self, y: &EquationOfMotionState) -> FourVector {
        let contravariant_metric = metric_contravariant(self.radius, self.a, y[1], y[2], y[3]);
        let covariant = Vector4::new(y[4], y[5], y[6], y[7]);
        let contravariant = contravariant_metric * covariant;

        FourVector::new_cartesian(
            contravariant[0],
            contravariant[1],
            contravariant[2],
            contravariant[3],
        )
    }
}

impl HasCoordinateSystem for Kerr {
    fn coordinate_system(&self) -> CoordinateSystem {
        Cartesian
    }
}

impl InnerProduct for Kerr {
    fn inner_product(&self, position: &Point, v: &FourVector, w: &FourVector) -> f64 {
        assert_eq!(v.coordinate_system, w.coordinate_system);
        let metric = metric(self.radius, self.a, position[1], position[2], position[3]);
        (v.vector.transpose() * metric * w.vector).x
    }
}

impl Signature for Kerr {
    fn signature(&self) -> [f64; 4] {
        [-1.0, 1.0, 1.0, 1.0]
    }
}

impl Geometry for Kerr {
    fn get_tetrad_at(&self, position: &Point) -> Tetrad {
        let x = position[1];
        let y = position[2];
        let z = position[3];

        let r_sqr = compute_r_sqr(self.a, x, y, z);
        let r = r_sqr.sqrt();
        let f = (r * r * r * self.radius) / (r * r * r * r + self.a * self.a * z * z);

        let k = k_vector(self.a, x, y, z);
        let k_x = k[1];
        let k_y = k[2];
        let k_z = k[3];

        let alpha = 1.0 / (1.0 + f).sqrt();
        let bfac = f / (1.0 + f);
        let beta = Vector3::new(bfac * k_x, bfac * k_y, bfac * k_z);

        debug!("alpha = {:?}, beta = {:?}", alpha, beta);
        let e_t = FourVector::new_cartesian(
            1.0 / alpha,
            -beta[0] / alpha,
            -beta[1] / alpha,
            -beta[2] / alpha,
        );
        let e_1 = FourVector::new_cartesian(0.0, 1.0, 0.0, 0.0);
        let e_2 = FourVector::new_cartesian(0.0, 0.0, 1.0, 0.0);
        let e_3 = FourVector::new_cartesian(0.0, 0.0, 0.0, 1.0);

        let basis = gram_schmidt(self, position, &vec![e_t, e_1, e_2, e_3]);

        debug!("{:?}", basis);

        Tetrad::new(position.clone(), basis[0], basis[1], basis[2], basis[3])
    }

    fn lorentz_transformation(&self, position: &Point, velocity: &FourVector) -> Matrix4<f64> {
        debug!(
            "Lorentz transformation at position {:?} with velocity {:?}",
            position, velocity
        );
        let mut matrix = Matrix4::zeros();
        let metric = metric(self.radius, self.a, position[1], position[2], position[3]);
        debug!("metric: {:?}", metric);
        let tetrad_t = self.get_tetrad_at(position).t;

        let gamma = -(tetrad_t.vector.transpose() * metric * velocity.vector)[(0, 0)];
        debug!("gamma: {}", gamma);

        debug!(
            "scalar prod tetrad_t: {:?}",
            self.inner_product(position, &tetrad_t, &tetrad_t)
        );
        debug!(
            "scalar prod velocity: {:?}",
            self.inner_product(position, velocity, velocity)
        );

        let uv = tetrad_t.vector + velocity.vector;
        let uv_lower = metric * uv;

        debug!("uv = {:?}", uv);
        debug!("uv_lower = {:?}", uv_lower);

        for mu in 0..4 {
            for nu in 0..4 {
                let mut res = 0.0;
                if mu == nu {
                    res = 1.0;
                }

                let a = 1.0 / (1.0 + gamma);
                let b = uv[mu];
                let c = uv_lower[nu];
                res += a * b * c;

                res -= 2.0 * (metric * tetrad_t.vector)[nu] * velocity.vector[mu];

                matrix[(mu, nu)] = res;
            }
        }
        debug!("lorentz transformation = {:?}", matrix);
        matrix
    }

    fn inside_horizon(&self, position: &Point) -> bool {
        if self.a > self.radius {
            return false;
        }
        let (x, y, z) = (position[1], position[2], position[3]);
        let r = compute_r_sqr(self.a, x, y, z).sqrt();
        let rp = 0.5 * self.radius
            + ((0.5 * self.radius) * (0.5 * self.radius) - self.a * self.a).sqrt();
        r <= rp + self.horizon_epsilon
    }

    fn closed_orbit(&self, position: &Point, step_index: usize, max_steps: usize) -> bool {
        let r = position.get_as_spherical()[0];
        if step_index == max_steps - 1 && r <= self.radius {
            return true;
        }
        false
    }

    fn get_geodesic_solver(&self, _ray: &Ray) -> Box<dyn GeodesicSolver> {
        Box::new(KerrSolver {
            radius: self.radius,
            a: self.a,
        })
    }
}

impl SupportQuantities for Kerr {
    fn get_stationary_velocity_at(&self, position: &Point) -> FourVector {
        let (x, y, z) = (position[1], position[2], position[3]);
        let r_sqr = compute_r_sqr(self.a, x, y, z);
        let r = r_sqr.sqrt();
        let f = (r * r * r * self.radius) / (r * r * r * r + self.a * self.a * z * z);
        FourVector::new_cartesian((1.0 - f).sqrt().recip(), 0.0, 0.0, 0.0)
    }

    // See https://arxiv.org/abs/1104.5499.
    fn get_circular_orbit_velocity_at(
        &self,
        position: &Point,
    ) -> Result<FourVector, RaytracerError> {
        let r = compute_r_sqr(self.a, position[1], position[2], position[3]).sqrt();
        let omega = self.angular_velocity(r);
        let ut = self.ut_contra(&position)?;
        let uphi = omega * ut;

        let velocity_spherical = FourVector::new_spherical(ut, 0.0, 0.0, uphi);

        let jacobian = self.jacobian_spherical_to_cartesian(position);
        let velocity_cartesian = jacobian * velocity_spherical.vector;

        Ok(FourVector::new_cartesian(
            velocity_cartesian[0],
            velocity_cartesian[1],
            velocity_cartesian[2],
            velocity_cartesian[3],
        ))
    }

    fn get_temperature_computer(
        &self,
        temperature: f64,
        _inner_radius: f64,
        outer_radius: f64,
    ) -> Result<Box<dyn TemperatureComputer>, RaytracerError> {
        Ok(Box::new(KerrTemperatureComputer::new(
            temperature,
            outer_radius,
            self.a,
            self.radius,
        )?))
    }
}

#[cfg(test)]
mod tests {
    use crate::geometry::four_vector::FourVector;
    use crate::geometry::geometry::{Geometry, InnerProduct, SupportQuantities};
    use crate::geometry::kerr::Kerr;
    use crate::geometry::point::Point;
    use crate::geometry::spherical_coordinates_helper::cartesian_to_spherical;
    use crate::rendering::camera::Camera;
    use crate::rendering::debug::save_rays_to_file;
    use crate::rendering::scene;
    use crate::rendering::scene::Scene;
    use approx::assert_abs_diff_eq;
    use std::f64::consts::PI;

    #[test]
    fn test_tetrad_orthonormal() {
        let position = Point::new_cartesian(2.0, 3.0, 4.0, 5.0);
        let geometry = Kerr::new(2.0, NO_ANGULAR_MOMENTUM, 1e-4);

        let tetrad = geometry.get_tetrad_at(&position);

        let k = tetrad.t + (-tetrad.z);
        let s = geometry.inner_product(&position, &k, &k);
        assert_abs_diff_eq!(s, 0.0);

        assert_abs_diff_eq!(
            geometry.inner_product(&position, &tetrad.t, &tetrad.t),
            -1.0
        );
        assert_abs_diff_eq!(geometry.inner_product(&position, &tetrad.x, &tetrad.x), 1.0);
        assert_abs_diff_eq!(geometry.inner_product(&position, &tetrad.y, &tetrad.y), 1.0);
        assert_abs_diff_eq!(geometry.inner_product(&position, &tetrad.z, &tetrad.z), 1.0);

        assert_abs_diff_eq!(geometry.inner_product(&position, &tetrad.t, &tetrad.x), 0.0);
        assert_abs_diff_eq!(geometry.inner_product(&position, &tetrad.t, &tetrad.y), 0.0);
        assert_abs_diff_eq!(geometry.inner_product(&position, &tetrad.t, &tetrad.z), 0.0);

        assert_abs_diff_eq!(geometry.inner_product(&position, &tetrad.x, &tetrad.y), 0.0);
        assert_abs_diff_eq!(geometry.inner_product(&position, &tetrad.x, &tetrad.z), 0.0);

        assert_abs_diff_eq!(geometry.inner_product(&position, &tetrad.y, &tetrad.z), 0.0);
    }

    const KERR_RADIUS_EPSILON: f64 = 1e-4;
    const NO_ANGULAR_MOMENTUM: f64 = 0.0;

    #[test]
    fn test_lorentz_transformed_tetrad_orthonormal() {
        let position = cartesian_to_spherical(&Point::new_cartesian(2.0, 3.0, 4.0, 5.0));
        let radius = 2.0;
        let r = (position[1] * position[1] + position[2] * position[2] + position[3] * position[3])
            .sqrt();
        let a = 1.0 - radius / r;

        let geometry = Kerr::new(radius, NO_ANGULAR_MOMENTUM, KERR_RADIUS_EPSILON);

        let velocity = FourVector::new_cartesian(1.0 / a.sqrt(), 0.0, 0.0, 0.0);

        let original_tetrad = geometry.get_tetrad_at(&position);
        let tetrad = crate::rendering::camera::lorentz_transform_tetrad(
            &geometry,
            &original_tetrad,
            &position,
            &velocity,
        );

        let k = tetrad.t + (-tetrad.z);
        let s = geometry.inner_product(&position, &k, &k);
        assert_abs_diff_eq!(s, 0.0, epsilon = 1e-8);

        assert_abs_diff_eq!(
            tetrad.t.get_as_vector(),
            velocity.get_as_vector(),
            epsilon = 1e-6
        );

        assert_abs_diff_eq!(geometry.inner_product(&position, &tetrad.t, &tetrad.x), 0.0);
        assert_abs_diff_eq!(geometry.inner_product(&position, &tetrad.t, &tetrad.y), 0.0);
        assert_abs_diff_eq!(geometry.inner_product(&position, &tetrad.t, &tetrad.z), 0.0);

        assert_abs_diff_eq!(geometry.inner_product(&position, &tetrad.x, &tetrad.y), 0.0);
        assert_abs_diff_eq!(geometry.inner_product(&position, &tetrad.x, &tetrad.z), 0.0);

        assert_abs_diff_eq!(geometry.inner_product(&position, &tetrad.y, &tetrad.z), 0.0);
    }

    #[ignore]
    #[test]
    fn save_camera_rays() {
        let rows = 30;
        let cols = 30;

        let position = Point::new_cartesian(0.0, 0.0, 0.0, -10.0);
        let radius = 2.0;
        let geometry = Kerr::new(radius, NO_ANGULAR_MOMENTUM, 1e-4);
        let velocity = FourVector::new_cartesian(1.0, 0.0, 0.0, 0.0);
        let camera = Camera::new(
            position,
            velocity,
            PI / 2.0,
            rows,
            cols,
            0.0,
            PI / 2.0,
            PI / 2.0,
            &Kerr::new(2.0, NO_ANGULAR_MOMENTUM, 1e-4),
        )
        .unwrap();

        save_rays_to_file(rows, cols, &position, geometry, camera);
    }

    #[test]
    fn test_kerr_ray() {
        let position = cartesian_to_spherical(&Point::new_cartesian(2.0, 3.0, 4.0, 5.0));
        let radius = 2.0;
        let geometry = Kerr::new(radius, NO_ANGULAR_MOMENTUM, 1e-4);
        let r = (position[1] * position[1] + position[2] * position[2] + position[3] * position[3])
            .sqrt();
        let a = 1.0 - radius / r;
        let velocity = FourVector::new_cartesian(1.0 / a.sqrt(), 0.0, 0.0, 0.0);
        let camera = Camera::new(
            position,
            velocity,
            PI / 2.0,
            11,
            11,
            0.0,
            PI / 2.0,
            PI / 2.0,
            &Kerr::new(2.0, NO_ANGULAR_MOMENTUM, 1e-4),
        )
        .unwrap();

        let ray = camera.get_ray_for(1, 6);
        assert_abs_diff_eq!(
            geometry.inner_product(&position, &ray.momentum, &ray.momentum),
            0.0,
            epsilon = 1e-8
        );
    }

    fn create_camera(position: Point, radius: f64) -> Camera {
        let r = (position[1] * position[1] + position[2] * position[2] + position[3] * position[3])
            .sqrt();
        let a = 1.0 - radius / r;
        let momentum = FourVector::new_cartesian(1.0 / a.sqrt(), 0.0, 0.0, 0.0);
        let camera = Camera::new(
            position,
            momentum,
            PI / 2.0,
            11,
            11,
            0.0,
            PI / 2.0,
            PI / 2.0,
            &Kerr::new(radius, NO_ANGULAR_MOMENTUM, 1e-4),
        )
        .unwrap();
        camera
    }

    #[test]
    fn test_ray_null_condition_momentum() {
        let position = Point::new_cartesian(2.0, 5.0, 0.0, 0.0);
        let radius = 2.0;
        let geometry = Kerr::new(radius, NO_ANGULAR_MOMENTUM, 1e-4);
        let camera = create_camera(position, radius);

        let m_s = geometry.inner_product(&position, &camera.velocity, &camera.velocity);
        assert_abs_diff_eq!(m_s, -1.0, epsilon = 1e-8);

        for i in 1..11 {
            let ray = camera.get_ray_for(6, i);
            let m_s = geometry.inner_product(&position, &ray.momentum, &ray.momentum);
            assert_abs_diff_eq!(m_s, 0.0, epsilon = 1e-8);
        }
        for i in 1..11 {
            let ray = camera.get_ray_for(i, 6);
            let m_s = geometry.inner_product(&position, &ray.momentum, &ray.momentum);
            assert_abs_diff_eq!(m_s, 0.0, epsilon = 1e-8);
        }
    }

    #[test]
    fn test_trajectories_equal_with_rotated_momentum() {
        let position = Point::new_cartesian(2.0, 1.0, 0.0, 0.0);
        let radius = 0.0;
        let geometry = Kerr::new(radius, NO_ANGULAR_MOMENTUM, 1e-4);
        let camera = create_camera(position, radius);
        let scene: Scene<Kerr> =
            scene::test_scene::create_scene_with_camera(1.0, 2.0, 7.0, &geometry, camera, 1e-5)
                .unwrap();

        let ray_a = scene.camera.get_ray_for(5, 10);
        let ray_b = scene.camera.get_ray_for(0, 5);

        println!("ray_a: {:?}", ray_a);
        println!("ray_b: {:?}", ray_b);

        // ensure rays are rotated by 90 degrees.
        assert_abs_diff_eq!(ray_a.momentum.vector[2], -ray_b.momentum.vector[3]);
        assert_abs_diff_eq!(ray_a.momentum.vector[3], ray_b.momentum.vector[2]);

        let (trajectory_a, _) = scene
            .integrator
            .integrate(&ray_a)
            .expect("unable to integrate ray_a");
        let (trajectory_b, _) = scene
            .integrator
            .integrate(&ray_b)
            .expect("unable to integrate ray_b");
        assert_eq!(trajectory_a.len(), trajectory_b.len());

        for i in 0..trajectory_a.len() {
            let step_a = &trajectory_a[i];
            let step_b = &trajectory_b[i];

            assert_abs_diff_eq!(step_a.x[1], step_b.x[1], epsilon = 1e-5);
            assert_abs_diff_eq!(step_a.x[2], -step_b.x[3], epsilon = 1e-5);
            assert_abs_diff_eq!(step_a.x[3], step_b.x[2], epsilon = 1e-5);
        }
    }

    #[test]
    fn test_circular_orbit_velocity() {
        let radius = 1.0;
        let geometry = Kerr::new(radius, NO_ANGULAR_MOMENTUM, 1e-4);

        let position = Point::new_cartesian(0.0, 0.0, 3.0, 0.0);
        let velocity = geometry.get_circular_orbit_velocity_at(&position).unwrap();

        assert_abs_diff_eq!(velocity[0], 1.414213562373095, epsilon = 1e-8);
        assert_abs_diff_eq!(velocity[1], -0.5773502691896257, epsilon = 1e-8);
        assert_abs_diff_eq!(velocity[2], 0.0, epsilon = 1e-8);
        assert_abs_diff_eq!(velocity[3], 0.0, epsilon = 1e-8);
    }
}

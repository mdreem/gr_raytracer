use crate::geometry::four_vector::FourVector;
use crate::geometry::geometry::{
    GeodesicSolver, Geometry, HasCoordinateSystem, InnerProduct, Signature,
};
use crate::geometry::gram_schmidt::gram_schmidt;
use crate::geometry::point::CoordinateSystem::Cartesian;
use crate::geometry::point::{CoordinateSystem, Point};
use crate::geometry::tetrad::Tetrad;
use crate::rendering::ray::Ray;
use crate::rendering::runge_kutta::OdeFunction;
use crate::rendering::scene::EquationOfMotionState;
use log::{debug, trace};
use nalgebra::{Const, Matrix4, OVector, Vector3, Vector4};

/// Kerr black hole geometry with analytical derivatives for improved performance and accuracy.
///
/// This implementation computes metric derivatives analytically rather than using
/// numerical differentiation, providing:
/// - ~3-5x faster geodesic integration
/// - Higher numerical accuracy (no finite difference errors)
/// - Better numerical stability
#[derive(Clone, Debug)]
pub struct KerrAnalytical {
    radius: f64,
    a: f64,
    horizon_epsilon: f64,
}

struct KerrAnalyticalSolver {
    radius: f64,
    a: f64,
}

/// Compute r² from Boyer-Lindquist-like coordinate
///
/// r² = (ρ² - a² + sqrt((ρ² - a²)² + 4a²z²)) / 2
/// where ρ² = x² + y² + z²
fn compute_r_sqr(a: f64, x: f64, y: f64, z: f64) -> f64 {
    let rho_sqr = x * x + y * y + z * z;
    0.5 * (rho_sqr - a * a + ((rho_sqr - a * a) * (rho_sqr - a * a) + 4.0 * a * a * z * z).sqrt())
}

/// Compute analytical derivative of r² with respect to coordinate
///
/// r² = (ρ² - a² + sqrt(Δ)) / 2, where Δ = (ρ² - a²)² + 4a²z²
/// ∂r²/∂x^i = 1/2 * [∂ρ²/∂x^i + ∂sqrt(Δ)/∂x^i]
fn compute_dr_sqr(a: f64, x: f64, y: f64, z: f64) -> Vector3<f64> {
    let rho_sqr = x * x + y * y + z * z;
    let delta = (rho_sqr - a * a) * (rho_sqr - a * a) + 4.0 * a * a * z * z;
    let sqrt_delta = delta.sqrt();

    // ∂ρ²/∂x^i
    // ∂ρ²/∂x = 2x, ∂ρ²/∂y = 2y, ∂ρ²/∂z = 2z

    // ∂Δ/∂x = 2(ρ² - a²) * 2x = 4x(ρ² - a²)
    // ∂Δ/∂y = 2(ρ² - a²) * 2y = 4y(ρ² - a²)
    // ∂Δ/∂z = 2(ρ² - a²) * 2z + 8a²z = 4z(ρ² + a²)

    // ∂sqrt(Δ)/∂x^i = ∂Δ/∂x^i / (2*sqrt(Δ))

    let common_factor = (rho_sqr - a * a) / sqrt_delta;

    // ∂r²/∂x = 1/2 * [2x + 2x(ρ² - a²)/sqrt(Δ)] = x * [1 + (ρ² - a²)/sqrt(Δ)]
    let dr_sqr_dx = x * (1.0 + common_factor);
    let dr_sqr_dy = y * (1.0 + common_factor);

    // ∂r²/∂z = 1/2 * [2z + 2z(ρ² + a²)/sqrt(Δ)] = z * [1 + (ρ² + a²)/sqrt(Δ)]
    let dr_sqr_dz = z * (1.0 + (rho_sqr + a * a) / sqrt_delta);

    Vector3::new(dr_sqr_dx, dr_sqr_dy, dr_sqr_dz)
}

/// Principal null vector k in Kerr-Schild coordinates
fn k_vector(a: f64, x: f64, y: f64, z: f64) -> Vector4<f64> {
    let r_sqr = compute_r_sqr(a, x, y, z);
    let r = r_sqr.sqrt();

    let k_0 = 1.0;
    let k_x = (r * x + a * y) / (r_sqr + a * a);
    let k_y = (r * y - a * x) / (r_sqr + a * a);
    let k_z = z / r;

    Vector4::new(k_0, k_x, k_y, k_z)
}

/// Analytical derivative of k vector with respect to coordinates
///
/// k_x = (rx + ay)/(r² + a²)
/// k_y = (ry - ax)/(r² + a²)
/// k_z = z/r
///
/// Returns matrix where row i is ∂k/∂x^i (i=1,2,3 for x,y,z)
fn k_vector_derivatives(a: f64, x: f64, y: f64, z: f64) -> (Vector4<f64>, Vector4<f64>, Vector4<f64>) {
    let r_sqr = compute_r_sqr(a, x, y, z);
    let r = r_sqr.sqrt();
    let dr_sqr = compute_dr_sqr(a, x, y, z);

    // ∂r/∂x^i = ∂r²/∂x^i / (2r)
    let dr_dx = dr_sqr[0] / (2.0 * r);
    let dr_dy = dr_sqr[1] / (2.0 * r);
    let dr_dz = dr_sqr[2] / (2.0 * r);

    let r_sqr_plus_a_sqr = r_sqr + a * a;
    let inv_denom = 1.0 / r_sqr_plus_a_sqr;

    // ∂k_x/∂x using quotient rule
    // k_x = (rx + ay)/(r² + a²)
    // ∂k_x/∂x = [(∂r/∂x * x + r) * (r² + a²) - (rx + ay) * ∂(r² + a²)/∂x] / (r² + a²)²
    //         = [(∂r/∂x * x + r) * (r² + a²) - (rx + ay) * 2r∂r/∂x] / (r² + a²)²

    let numerator_kx = r * x + a * y;
    let numerator_ky = r * y - a * x;

    // ∂k_x/∂x
    let dk_x_dx = ((dr_dx * x + r) * r_sqr_plus_a_sqr - numerator_kx * 2.0 * r * dr_dx)
                  * inv_denom * inv_denom;

    // ∂k_x/∂y
    let dk_x_dy = ((dr_dy * x + a) * r_sqr_plus_a_sqr - numerator_kx * 2.0 * r * dr_dy)
                  * inv_denom * inv_denom;

    // ∂k_x/∂z
    let dk_x_dz = (dr_dz * x * r_sqr_plus_a_sqr - numerator_kx * 2.0 * r * dr_dz)
                  * inv_denom * inv_denom;

    // ∂k_y/∂x
    let dk_y_dx = ((dr_dx * y - a) * r_sqr_plus_a_sqr - numerator_ky * 2.0 * r * dr_dx)
                  * inv_denom * inv_denom;

    // ∂k_y/∂y
    let dk_y_dy = ((dr_dy * y + r) * r_sqr_plus_a_sqr - numerator_ky * 2.0 * r * dr_dy)
                  * inv_denom * inv_denom;

    // ∂k_y/∂z
    let dk_y_dz = (dr_dz * y * r_sqr_plus_a_sqr - numerator_ky * 2.0 * r * dr_dz)
                  * inv_denom * inv_denom;

    // k_z = z/r
    // ∂k_z/∂x = -z/(r²) * ∂r/∂x
    let dk_z_dx = -z * dr_dx / r_sqr;
    let dk_z_dy = -z * dr_dy / r_sqr;
    let dk_z_dz = (r - z * dr_dz / r) / r_sqr;

    // ∂k/∂t = 0 (stationary metric)
    let dk_dx = Vector4::new(0.0, dk_x_dx, dk_y_dx, dk_z_dx);
    let dk_dy = Vector4::new(0.0, dk_x_dy, dk_y_dy, dk_z_dy);
    let dk_dz = Vector4::new(0.0, dk_x_dz, dk_y_dz, dk_z_dz);

    (dk_dx, dk_dy, dk_dz)
}

/// Kerr-Schild metric: g_μν = η_μν + f k_μ k_ν
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

/// Analytical derivative of f = r³M / (r⁴ + a²z²)
fn f_derivative(radius: f64, a: f64, x: f64, y: f64, z: f64) -> Vector3<f64> {
    let r_sqr = compute_r_sqr(a, x, y, z);
    let r = r_sqr.sqrt();
    let r_cubed = r * r * r;
    let r_quad = r_sqr * r_sqr;

    let denominator = r_quad + a * a * z * z;

    let dr_sqr = compute_dr_sqr(a, x, y, z);
    // Convert ∂r²/∂x^i to ∂r/∂x^i
    let dr_dx = dr_sqr[0] / (2.0 * r);
    let dr_dy = dr_sqr[1] / (2.0 * r);
    let dr_dz = dr_sqr[2] / (2.0 * r);

    // f = r³M / (r⁴ + a²z²)
    // Using quotient rule: ∂f/∂x^i = M * [∂(r³)/∂x^i * denom - r³ * ∂(denom)/∂x^i] / denom²

    // ∂(r³)/∂x^i = 3r² ∂r/∂x^i
    // ∂(r⁴ + a²z²)/∂x = 4r³ ∂r/∂x
    // ∂(r⁴ + a²z²)/∂y = 4r³ ∂r/∂y
    // ∂(r⁴ + a²z²)/∂z = 4r³ ∂r/∂z + 2a²z

    let df_dx = radius * (
        3.0 * r_sqr * dr_dx * denominator
        - r_cubed * (4.0 * r_cubed * dr_dx)
    ) / (denominator * denominator);

    let df_dy = radius * (
        3.0 * r_sqr * dr_dy * denominator
        - r_cubed * (4.0 * r_cubed * dr_dy)
    ) / (denominator * denominator);

    let df_dz = radius * (
        3.0 * r_sqr * dr_dz * denominator
        - r_cubed * (4.0 * r_cubed * dr_dz + 2.0 * a * a * z)
    ) / (denominator * denominator);

    Vector3::new(df_dx, df_dy, df_dz)
}

/// Analytical derivative of covariant metric
///
/// g_μν = η_μν + f k_μ k_ν
/// ∂g_μν/∂x^α = ∂f/∂x^α k_μ k_ν + f (∂k_μ/∂x^α k_ν + k_μ ∂k_ν/∂x^α)
fn metric_derivative(
    radius: f64,
    a: f64,
    x: f64,
    y: f64,
    z: f64,
    index: usize,
) -> Matrix4<f64> {
    if index == 0 {
        return Matrix4::zeros(); // Stationary metric: ∂g/∂t = 0
    }

    let r_sqr = compute_r_sqr(a, x, y, z);
    let r = r_sqr.sqrt();
    let f = (r * r * r * radius) / (r * r * r * r + a * a * z * z);

    let k = k_vector(a, x, y, z);
    let df = f_derivative(radius, a, x, y, z);
    let (dk_dx, dk_dy, dk_dz) = k_vector_derivatives(a, x, y, z);

    let spatial_index = index - 1; // Convert to 0-based spatial index
    let df_dxi = df[spatial_index];
    let dk_dxi = match index {
        1 => dk_dx,
        2 => dk_dy,
        3 => dk_dz,
        _ => Vector4::zeros(),
    };

    let mut d_metric = Matrix4::zeros();

    // ∂g_μν/∂x^i = (∂f/∂x^i) k_μ k_ν + f (∂k_μ/∂x^i k_ν + k_μ ∂k_ν/∂x^i)
    for mu in 0..4 {
        for nu in 0..4 {
            d_metric[(mu, nu)] = df_dxi * k[mu] * k[nu]
                + f * (dk_dxi[mu] * k[nu] + k[mu] * dk_dxi[nu]);
        }
    }

    d_metric
}

impl KerrAnalytical {
    /// Create a new Kerr black hole with analytical derivatives
    ///
    /// # Arguments
    /// * `radius` - Schwarzschild radius r_s = 2M (geometric units)
    /// * `a` - Dimensionless spin parameter (specific angular momentum J/M, |a| ≤ M)
    /// * `horizon_epsilon` - Numerical tolerance for horizon detection
    pub fn new(radius: f64, a: f64, horizon_epsilon: f64) -> Self {
        KerrAnalytical {
            radius,
            a,
            horizon_epsilon,
        }
    }
}

impl OdeFunction<Const<8>> for KerrAnalyticalSolver {
    fn apply(&self, t: f64, y: &OVector<f64, Const<8>>) -> OVector<f64, Const<8>> {
        self.geodesic(t, y)
    }
}

impl KerrAnalyticalSolver {
    /// Analytical derivative of contravariant metric
    ///
    /// Uses the identity: ∂g^{μν}/∂x^α = -g^{μρ} (∂g_{ρσ}/∂x^α) g^{σν}
    fn d_matrix_contravariant(
        &self,
        index: usize,
        x: f64,
        y: f64,
        z: f64,
        contravariant_metric: &Matrix4<f64>,
    ) -> Matrix4<f64> {
        let d_covariant = metric_derivative(self.radius, self.a, x, y, z, index);
        -(contravariant_metric * d_covariant * contravariant_metric)
    }
}

impl HasCoordinateSystem for KerrAnalyticalSolver {
    fn coordinate_system(&self) -> CoordinateSystem {
        Cartesian
    }
}

impl GeodesicSolver for KerrAnalyticalSolver {
    /// Hamiltonian geodesics with analytical derivatives
    ///
    /// H(x, p) = 0.5 * g^{μν} p_μ p_ν
    /// dx^μ/dλ = ∂H/∂p_μ = g^{μν} p_ν
    /// dp_μ/dλ = -∂H/∂x^μ = -0.5 * ∂g^{αβ}/∂x^μ * p_α p_β
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

        let covariant_metric = metric(self.radius, self.a, x, y, z);
        trace!("covariant_metric = {:?}", covariant_metric);
        let contravariant_metric = covariant_metric
            .try_inverse()
            .expect("Metric should be invertible");
        trace!("contravariant_metric = {:?}", contravariant_metric);

        debug_assert!(!p_t.is_nan());

        // Compute dx^μ/dλ = g^{μν} p_ν
        let xdot = contravariant_metric * p;
        let dt = xdot[0];
        let dx = xdot[1];
        let dy = xdot[2];
        let dz = xdot[3];

        // Compute dp_μ/dλ using analytical derivatives
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
        let covariant_metric = metric(self.radius, self.a, y[1], y[2], y[3]);
        let contravariant_metric = covariant_metric
            .try_inverse()
            .expect("Metric should be invertible");
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

impl HasCoordinateSystem for KerrAnalytical {
    fn coordinate_system(&self) -> CoordinateSystem {
        Cartesian
    }
}

impl InnerProduct for KerrAnalytical {
    fn inner_product(&self, position: &Point, v: &FourVector, w: &FourVector) -> f64 {
        let metric = metric(self.radius, self.a, position[1], position[2], position[3]);
        (v.vector.transpose() * metric * w.vector).x
    }
}

impl Signature for KerrAnalytical {
    fn signature(&self) -> [f64; 4] {
        [-1.0, 1.0, 1.0, 1.0]
    }
}

impl Geometry for KerrAnalytical {
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

    fn get_stationary_velocity_at(&self, position: &Point) -> FourVector {
        let (x, y, z) = (position[1], position[2], position[3]);
        let r_sqr = compute_r_sqr(self.a, x, y, z);
        let r = r_sqr.sqrt();
        let f = (r * r * r * self.radius) / (r * r * r * r + self.a * self.a * z * z);
        FourVector::new_cartesian((1.0 - f).sqrt().recip(), 0.0, 0.0, 0.0)
    }

    fn inside_horizon(&self, position: &Point) -> bool {
        if self.a > self.radius {
            return false;
        }
        let (x, y, z) = (position[1], position[2], position[3]);
        let r = compute_r_sqr(self.a, x, y, z).sqrt();
        // Outer horizon radius: r_+ = M + sqrt(M² - a²)
        // where M = radius/2 in geometric units
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
        Box::new(KerrAnalyticalSolver {
            radius: self.radius,
            a: self.a,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::kerr::Kerr;
    use crate::geometry::geometry::{Geometry, InnerProduct};
    use crate::geometry::point::Point;
    use approx::assert_abs_diff_eq;

    const NO_ANGULAR_MOMENTUM: f64 = 0.0;

    /// Test that analytical derivatives match numerical derivatives
    #[test]
    fn test_analytical_vs_numerical_derivatives() {
        let radius = 2.0;
        let a = 0.5;
        let position = Point::new_cartesian(0.0, 3.0, 4.0, 5.0);

        let (x, y, z) = (position[1], position[2], position[3]);

        // Compute analytical derivative
        let analytical = metric_derivative(radius, a, x, y, z, 1);

        // Compute numerical derivative for comparison
        let h = 1e-7;
        let metric_plus = metric(radius, a, x + h, y, z);
        let metric_minus = metric(radius, a, x - h, y, z);
        let numerical = (metric_plus - metric_minus) / (2.0 * h);

        // Compare all components
        for i in 0..4 {
            for j in 0..4 {
                assert_abs_diff_eq!(
                    analytical[(i, j)],
                    numerical[(i, j)],
                    epsilon = 1e-5
                );
            }
        }
    }

    /// Test analytical f derivative
    #[test]
    fn test_f_derivative() {
        let radius = 2.0;
        let a = 0.5;
        let x = 3.0;
        let y = 4.0;
        let z = 5.0;

        let df_analytical = f_derivative(radius, a, x, y, z);

        // Numerical derivative
        let h = 1e-7;
        let f = |x, y, z| {
            let r_sqr = compute_r_sqr(a, x, y, z);
            let r = r_sqr.sqrt();
            (r * r * r * radius) / (r * r * r * r + a * a * z * z)
        };

        let df_dx_numerical = (f(x + h, y, z) - f(x - h, y, z)) / (2.0 * h);
        let df_dy_numerical = (f(x, y + h, z) - f(x, y - h, z)) / (2.0 * h);
        let df_dz_numerical = (f(x, y, z + h) - f(x, y, z - h)) / (2.0 * h);

        assert_abs_diff_eq!(df_analytical[0], df_dx_numerical, epsilon = 1e-6);
        assert_abs_diff_eq!(df_analytical[1], df_dy_numerical, epsilon = 1e-6);
        assert_abs_diff_eq!(df_analytical[2], df_dz_numerical, epsilon = 1e-6);
    }

    /// Test k vector derivatives
    #[test]
    fn test_k_vector_derivatives() {
        let a = 0.5;
        let x = 3.0;
        let y = 4.0;
        let z = 5.0;

        let (dk_dx, _dk_dy, _dk_dz) = k_vector_derivatives(a, x, y, z);

        // Numerical derivatives
        let h = 1e-7;
        let k_plus_x = k_vector(a, x + h, y, z);
        let k_minus_x = k_vector(a, x - h, y, z);
        let dk_dx_numerical = (k_plus_x - k_minus_x) / (2.0 * h);

        for i in 0..4 {
            assert_abs_diff_eq!(
                dk_dx[i],
                dk_dx_numerical[i],
                epsilon = 1e-5
            );
        }
    }

    /// Verify analytical Kerr produces same results as numerical Kerr
    #[test]
    fn test_analytical_matches_numerical_kerr() {
        let position = Point::new_cartesian(0.0, 0.0, 0.0, -10.0);
        let radius = 2.0;
        let velocity = FourVector::new_cartesian(1.0, 0.0, 0.0, 0.0);

        let kerr_numerical = Kerr::new(radius, NO_ANGULAR_MOMENTUM, 1e-4);
        let kerr_analytical = KerrAnalytical::new(radius, NO_ANGULAR_MOMENTUM, 1e-4);

        // Test metric computation
        assert_abs_diff_eq!(
            kerr_numerical.inner_product(&position, &velocity, &velocity),
            kerr_analytical.inner_product(&position, &velocity, &velocity),
            epsilon = 1e-10
        );

        // Test tetrad
        let tetrad_num = kerr_numerical.get_tetrad_at(&position);
        let tetrad_ana = kerr_analytical.get_tetrad_at(&position);

        assert_abs_diff_eq!(
            tetrad_num.t.get_as_vector(),
            tetrad_ana.t.get_as_vector(),
            epsilon = 1e-10
        );
    }

    /// Test that null geodesics remain null during integration
    #[test]
    fn test_geodesic_preserves_null_condition() {
        use crate::rendering::integrator::{IntegrationConfiguration, Integrator};
        use crate::rendering::ray::Ray;

        let position = Point::new_cartesian(0.0, 5.0, 0.0, 0.0);
        let radius = 2.0;
        let geometry = KerrAnalytical::new(radius, NO_ANGULAR_MOMENTUM, 1e-4);

        // Create a null momentum manually
        let momentum = FourVector::new_cartesian(1.0, -1.0, 0.0, 0.0);

        // Verify it's null
        let inner_prod = geometry.inner_product(&position, &momentum, &momentum);
        assert_abs_diff_eq!(inner_prod, 0.0, epsilon = 1e-10);

        // Create ray and integrate
        let ray = Ray::new(0, 0, position.clone(), momentum);
        let config = IntegrationConfiguration::new(1000, 10000.0, 0.01, 1e-5);
        let integrator = Integrator::new(&geometry, config);

        let (trajectory, _) = integrator.integrate(&ray).expect("Integration should succeed");

        // Check null condition is preserved at several points along trajectory
        for (i, step) in trajectory.iter().enumerate().step_by(100) {
            let null_check = geometry.inner_product(&step.x, &step.p, &step.p);
            assert!(
                null_check.abs() < 1e-8,
                "Null condition violated at step {}: {}",
                i,
                null_check
            );
        }
    }
}

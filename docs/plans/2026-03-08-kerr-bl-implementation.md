# KerrBL Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add a new `KerrBL` geometry that solves Kerr geodesics using separated equations (Carter constant / complete integrability) in Boyer-Lindquist coordinates with Mino time parameterization.

**Architecture:** New standalone geometry in `src/geometry/kerr_bl.rs` with a `KerrBLSolver` that stores constants of motion (E, L_z, Q) computed once per ray. Requires a new `BoyerLindquist` coordinate system variant since BL-to-Cartesian conversion differs from standard spherical when a≠0.

**Tech Stack:** Rust, nalgebra, existing GeodesicSolver/Geometry trait system, RKF45 integrator

---

### Task 1: Add BoyerLindquist Coordinate System

**Files:**
- Modify: `src/geometry/point.rs`
- Modify: `src/geometry/four_vector.rs`
- Modify: `src/rendering/scene.rs`
- Modify: `src/rendering/camera.rs`

The BL-to-Cartesian conversion for Kerr is:
```
x = (r cosφ - a sinφ) sinθ
y = (r sinφ + a cosφ) sinθ
z = r cosθ
```
This differs from standard spherical (x = r sinθ cosφ) when a≠0. We need a new coordinate system variant to get correct intersection testing.

**Step 1: Write failing tests for BL coordinate conversion**

In `src/geometry/point.rs`, add to the existing test module (or create one):

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_boyer_lindquist_to_cartesian_a_zero() {
        // When a=0, BL should match standard spherical
        let bl = Point::new(0.0, 5.0, 1.2, 0.8, CoordinateSystem::BoyerLindquist { a: 0.0 });
        let sph = Point::new_spherical(0.0, 5.0, 1.2, 0.8);
        let bl_cart = bl.to_cartesian();
        let sph_cart = sph.to_cartesian();
        assert_abs_diff_eq!(bl_cart.vector, sph_cart.vector, epsilon = 1e-12);
    }

    #[test]
    fn test_boyer_lindquist_to_cartesian_nonzero_a() {
        let a = 0.5;
        let r = 5.0;
        let theta = 1.2;
        let phi = 0.8;
        let bl = Point::new(0.0, r, theta, phi, CoordinateSystem::BoyerLindquist { a });
        let cart = bl.to_cartesian();
        // x = (r cosφ - a sinφ) sinθ
        let expected_x = (r * phi.cos() - a * phi.sin()) * theta.sin();
        // y = (r sinφ + a cosφ) sinθ
        let expected_y = (r * phi.sin() + a * phi.cos()) * theta.sin();
        // z = r cosθ
        let expected_z = r * theta.cos();
        assert_abs_diff_eq!(cart[1], expected_x, epsilon = 1e-12);
        assert_abs_diff_eq!(cart[2], expected_y, epsilon = 1e-12);
        assert_abs_diff_eq!(cart[3], expected_z, epsilon = 1e-12);
    }

    #[test]
    fn test_boyer_lindquist_radial_distance() {
        let bl = Point::new(0.0, 7.0, 1.0, 2.0, CoordinateSystem::BoyerLindquist { a: 0.5 });
        assert_abs_diff_eq!(bl.radial_distance_spatial_part_squared(), 49.0);
    }
}
```

**Step 2: Run tests to verify they fail**

Run: `cargo test --lib point::tests -- --nocapture 2>&1 | head -30`
Expected: Compilation error — `BoyerLindquist` variant does not exist.

**Step 3: Add BoyerLindquist variant and implement conversions**

In `src/geometry/point.rs`, update the enum:

```rust
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CoordinateSystem {
    Cartesian,
    Spherical,
    BoyerLindquist { a: f64 },
}
```

Update `to_cartesian`:
```rust
pub fn to_cartesian(self) -> Point {
    match self.coordinate_system {
        Cartesian => self,
        Spherical => spherical_to_cartesian(&self),
        BoyerLindquist { a } => {
            let t = self.vector[0];
            let r = self.vector[1];
            let theta = self.vector[2];
            let phi = self.vector[3];
            let x = (r * phi.cos() - a * phi.sin()) * theta.sin();
            let y = (r * phi.sin() + a * phi.cos()) * theta.sin();
            let z = r * theta.cos();
            Point::new_cartesian(t, x, y, z)
        }
    }
}
```

Update `get_spatial_vector_cartesian`:
```rust
pub fn get_spatial_vector_cartesian(self) -> Vector3<f64> {
    match self.coordinate_system {
        Cartesian => Vector3::new(self.vector[1], self.vector[2], self.vector[3]),
        Spherical => {
            let v = spherical_to_cartesian(&self);
            Vector3::new(v[1], v[2], v[3])
        }
        BoyerLindquist { .. } => {
            let v = self.to_cartesian();
            Vector3::new(v[1], v[2], v[3])
        }
    }
}
```

Update `get_as_spherical`:
```rust
pub fn get_as_spherical(self) -> Vector3<f64> {
    match self.coordinate_system {
        Cartesian => {
            let v = cartesian_to_spherical(&self);
            Vector3::new(v[1], v[2], v[3])
        }
        Spherical => Vector3::new(
            self.vector[1],
            wrap_theta(self.vector[2]),
            wrap_phi(self.vector[3]),
        ),
        BoyerLindquist { .. } => Vector3::new(
            self.vector[1],
            wrap_theta(self.vector[2]),
            wrap_phi(self.vector[3]),
        ),
    }
}
```

Update `radial_distance_spatial_part_squared`:
```rust
pub fn radial_distance_spatial_part_squared(&self) -> f64 {
    let v = self.vector;
    match self.coordinate_system {
        Cartesian => v[1] * v[1] + v[2] * v[2] + v[3] * v[3],
        Spherical | BoyerLindquist { .. } => {
            let r = v[1];
            r * r
        }
    }
}
```

Update `Neg` impl debug_assert to also accept BL (or remove the Cartesian-only assertion if it blocks BL).

In `src/rendering/scene.rs`, update `get_position` (line 32):
```rust
pub fn get_position(y: &EquationOfMotionState, coordinate_system: CoordinateSystem) -> Point {
    match coordinate_system {
        CoordinateSystem::Cartesian => Point::new_cartesian(y[0], y[1], y[2], y[3]),
        CoordinateSystem::Spherical => spherical_to_cartesian(&Point::new(
            y[0], y[1], y[2], y[3], CoordinateSystem::Spherical,
        )),
        CoordinateSystem::BoyerLindquist { a } => {
            let bl = Point::new(y[0], y[1], y[2], y[3], CoordinateSystem::BoyerLindquist { a });
            bl.to_cartesian()
        }
    }
}
```

In `src/rendering/camera.rs`, update `spatial_basis_vector_cartesian` (line 83) to add BL case. The Jacobian for BL Kerr is:
```
dx/dr = sinθ cosφ
dx/dθ = (r cosφ - a sinφ) cosθ
dx/dφ = (-r sinφ - a cosφ) sinθ

dy/dr = sinθ sinφ
dy/dθ = (r sinφ + a cosφ) cosθ
dy/dφ = (r cosφ - a sinφ) sinθ

dz/dr = cosθ
dz/dθ = -r sinθ
dz/dφ = 0
```

```rust
CoordinateSystem::BoyerLindquist { a } => {
    let spherical_position = position.get_as_spherical();
    let r = spherical_position[0];
    let theta = spherical_position[1];
    let phi = spherical_position[2];

    let dr = vector[1];
    let dtheta = vector[2];
    let dphi = vector[3];

    let (st, ct) = (theta.sin(), theta.cos());
    let (sp, cp) = (phi.sin(), phi.cos());

    Vector3::new(
        st * cp * dr + (r * cp - a * sp) * ct * dtheta + (-r * sp - a * cp) * st * dphi,
        st * sp * dr + (r * sp + a * cp) * ct * dtheta + (r * cp - a * sp) * st * dphi,
        ct * dr - r * st * dtheta,
    )
}
```

In `src/geometry/four_vector.rs`, add a new constructor:
```rust
pub fn new_boyer_lindquist(a: f64, t: f64, r: f64, theta: f64, phi: f64) -> FourVector {
    FourVector {
        coordinate_system: CoordinateSystem::BoyerLindquist { a },
        vector: Vector4::new(t, r, theta, phi),
    }
}
```

**Step 4: Run tests to verify they pass**

Run: `cargo test --lib 2>&1 | tail -5`
Expected: All existing tests pass, new BL coordinate tests pass.

**Step 5: Commit**

```bash
git add src/geometry/point.rs src/geometry/four_vector.rs src/rendering/scene.rs src/rendering/camera.rs
git commit -m "feat: add BoyerLindquist coordinate system variant"
```

---

### Task 2: Add KerrBL Struct and BL Metric

**Files:**
- Create: `src/geometry/kerr_bl.rs`
- Modify: `src/geometry/mod.rs`

**Step 1: Write failing tests for the BL metric**

Create `src/geometry/kerr_bl.rs` with tests first:

```rust
use nalgebra::Matrix4;

fn sigma(r: f64, a: f64, theta: f64) -> f64 {
    r * r + a * a * theta.cos().powi(2)
}

fn delta(r: f64, r_s: f64, a: f64) -> f64 {
    r * r - r_s * r + a * a
}

/// Covariant BL metric g_μν. Signature (-,+,+,+).
fn metric_bl(r_s: f64, a: f64, r: f64, theta: f64) -> Matrix4<f64> {
    todo!()
}

/// Contravariant BL metric g^μν.
fn metric_bl_contravariant(r_s: f64, a: f64, r: f64, theta: f64) -> Matrix4<f64> {
    todo!()
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_bl_metric_inverse() {
        let r_s = 1.0;
        let a = 0.5;
        let test_points = [(5.0, 1.2), (3.0, 0.8), (10.0, 2.5)];
        for (r, theta) in test_points {
            let g = metric_bl(r_s, a, r, theta);
            let g_inv = g.try_inverse().expect("Metric should be invertible");
            let g_contra = metric_bl_contravariant(r_s, a, r, theta);
            for i in 0..4 {
                for j in 0..4 {
                    assert_abs_diff_eq!(g_contra[(i, j)], g_inv[(i, j)], epsilon = 1e-10);
                }
            }
        }
    }

    #[test]
    fn test_bl_metric_schwarzschild_limit() {
        // When a=0, BL Kerr should reduce to Schwarzschild metric
        let r_s = 2.0;
        let a = 0.0;
        let r = 5.0;
        let theta = 1.2;
        let g = metric_bl(r_s, a, r, theta);
        let a_factor = 1.0 - r_s / r;
        // g_tt = -(1 - r_s/r)
        assert_abs_diff_eq!(g[(0, 0)], -a_factor, epsilon = 1e-12);
        // g_rr = 1/a_factor = 1/(1 - r_s/r)
        assert_abs_diff_eq!(g[(1, 1)], 1.0 / a_factor, epsilon = 1e-12);
        // g_θθ = r²
        assert_abs_diff_eq!(g[(2, 2)], r * r, epsilon = 1e-12);
        // g_φφ = r² sin²θ
        assert_abs_diff_eq!(g[(3, 3)], r * r * theta.sin().powi(2), epsilon = 1e-12);
        // g_tφ = 0
        assert_abs_diff_eq!(g[(0, 3)], 0.0, epsilon = 1e-12);
    }
}
```

Register the module in `src/geometry/mod.rs`:
```rust
pub mod kerr_bl;
```

**Step 2: Run tests to verify they fail**

Run: `cargo test --lib kerr_bl 2>&1 | head -20`
Expected: FAIL — `todo!()` panics.

**Step 3: Implement the metric functions**

Replace the `todo!()` in `metric_bl`:
```rust
fn metric_bl(r_s: f64, a: f64, r: f64, theta: f64) -> Matrix4<f64> {
    let sig = sigma(r, a, theta);
    let del = delta(r, r_s, a);
    let sin_t = theta.sin();
    let sin2 = sin_t * sin_t;

    let g_tt = -(1.0 - r_s * r / sig);
    let g_rr = sig / del;
    let g_thth = sig;
    let g_phph = (r * r + a * a + a * a * r_s * r * sin2 / sig) * sin2;
    let g_tph = -a * r_s * r * sin2 / sig;

    let mut g = Matrix4::zeros();
    g[(0, 0)] = g_tt;
    g[(1, 1)] = g_rr;
    g[(2, 2)] = g_thth;
    g[(3, 3)] = g_phph;
    g[(0, 3)] = g_tph;
    g[(3, 0)] = g_tph;
    g
}
```

Replace the `todo!()` in `metric_bl_contravariant`:
```rust
fn metric_bl_contravariant(r_s: f64, a: f64, r: f64, theta: f64) -> Matrix4<f64> {
    let sig = sigma(r, a, theta);
    let del = delta(r, r_s, a);
    let sin_t = theta.sin();
    let sin2 = sin_t * sin_t;
    let a2 = a * a;
    let r2 = r * r;
    let big_a = (r2 + a2).powi(2) - del * a2 * sin2;

    let mut g = Matrix4::zeros();
    g[(0, 0)] = -big_a / (sig * del);
    g[(1, 1)] = del / sig;
    g[(2, 2)] = 1.0 / sig;
    g[(3, 3)] = (del - a2 * sin2) / (sig * del * sin2);
    g[(0, 3)] = -a * r_s * r / (sig * del);
    g[(3, 0)] = g[(0, 3)];
    g
}
```

Also add the `KerrBL` struct stub:
```rust
#[derive(Clone, Debug)]
pub struct KerrBL {
    pub radius: f64,  // r_s = 2M
    pub a: f64,       // spin parameter
    pub horizon_epsilon: f64,
}

impl KerrBL {
    pub fn new(radius: f64, a: f64, horizon_epsilon: f64) -> Self {
        KerrBL { radius, a, horizon_epsilon }
    }
}
```

**Step 4: Run tests to verify they pass**

Run: `cargo test --lib kerr_bl 2>&1 | tail -5`
Expected: PASS

**Step 5: Commit**

```bash
git add src/geometry/kerr_bl.rs src/geometry/mod.rs
git commit -m "feat: add KerrBL struct and BL metric functions"
```

---

### Task 3: Implement Geometry Trait Basics (InnerProduct, Signature, HasCoordinateSystem, get_radial_coordinate, inside_horizon)

**Files:**
- Modify: `src/geometry/kerr_bl.rs`

**Step 1: Write failing tests**

```rust
#[test]
fn test_inner_product_null_vector() {
    let kerr = KerrBL::new(1.0, 0.5, 1e-4);
    let r = 5.0;
    let theta = std::f64::consts::FRAC_PI_2;
    let position = Point::new(0.0, r, theta, 0.0, CoordinateSystem::BoyerLindquist { a: 0.5 });

    // Construct a null vector from metric: g_μν k^μ k^ν = 0
    // For a radial null ray: k^t and k^r with g_tt (k^t)^2 + g_rr (k^r)^2 = 0
    let g = metric_bl(1.0, 0.5, r, theta);
    let kt = 1.0;
    let kr = (-g[(0, 0)] / g[(1, 1)]).sqrt() * kt;
    let k = FourVector::new_boyer_lindquist(0.5, kt, kr, 0.0, 0.0);
    let ip = kerr.inner_product(&position, &k, &k);
    assert_abs_diff_eq!(ip, 0.0, epsilon = 1e-10);
}

#[test]
fn test_inside_horizon() {
    let kerr = KerrBL::new(1.0, 0.3, 1e-4);
    let m = 0.5; // M = r_s/2
    let r_plus = m + (m * m - 0.3 * 0.3).sqrt();

    let inside = Point::new(0.0, r_plus - 0.01, 1.0, 0.0, CoordinateSystem::BoyerLindquist { a: 0.3 });
    let outside = Point::new(0.0, r_plus + 0.1, 1.0, 0.0, CoordinateSystem::BoyerLindquist { a: 0.3 });
    assert!(kerr.inside_horizon(&inside));
    assert!(!kerr.inside_horizon(&outside));
}

#[test]
fn test_get_radial_coordinate() {
    let kerr = KerrBL::new(1.0, 0.5, 1e-4);
    let position = Point::new(0.0, 7.5, 1.2, 0.8, CoordinateSystem::BoyerLindquist { a: 0.5 });
    assert_abs_diff_eq!(kerr.get_radial_coordinate(&position), 7.5);
}
```

**Step 2: Run tests to verify they fail**

Run: `cargo test --lib kerr_bl 2>&1 | head -20`
Expected: Compilation error — traits not implemented.

**Step 3: Implement the trait impls**

```rust
use crate::geometry::geometry::{
    ConstantsOfMotion, GeodesicSolver, Geometry, HasCoordinateSystem,
    InnerProduct, Signature, SupportQuantities,
};
use crate::geometry::four_vector::FourVector;
use crate::geometry::point::{CoordinateSystem, Point};
use crate::geometry::tetrad::Tetrad;
use crate::rendering::ray::Ray;
use crate::rendering::raytracer::RaytracerError;
use nalgebra::Matrix4;

impl HasCoordinateSystem for KerrBL {
    fn coordinate_system(&self) -> CoordinateSystem {
        CoordinateSystem::BoyerLindquist { a: self.a }
    }
}

impl Signature for KerrBL {
    fn signature(&self) -> [f64; 4] {
        [-1.0, 1.0, 1.0, 1.0]
    }
}

impl InnerProduct for KerrBL {
    fn inner_product(&self, position: &Point, v: &FourVector, w: &FourVector) -> f64 {
        let r = position[1];
        let theta = position[2];
        let g = metric_bl(self.radius, self.a, r, theta);

        let mut result = 0.0;
        for mu in 0..4 {
            for nu in 0..4 {
                result += g[(mu, nu)] * v.vector[mu] * w.vector[nu];
            }
        }
        result
    }
}
```

Add `inside_horizon`, `get_radial_coordinate`, and a stub `Geometry` impl (with `todo!()` for methods not yet implemented like `get_tetrad_at`, `lorentz_transformation`, `get_geodesic_solver`, `get_constants_of_motion`, `closed_orbit`).

**Step 4: Run tests to verify they pass**

Run: `cargo test --lib kerr_bl 2>&1 | tail -5`
Expected: PASS

**Step 5: Commit**

```bash
git add src/geometry/kerr_bl.rs
git commit -m "feat: implement basic Geometry traits for KerrBL"
```

---

### Task 4: Implement KerrBLSolver — Separated Geodesic Equations

**Files:**
- Modify: `src/geometry/kerr_bl.rs`

This is the core physics task. The solver stores E, L_z, Q and implements the separated 2nd-order equations.

**Step 1: Write failing tests for potential functions and their derivatives**

```rust
#[test]
fn test_potential_r_at_initial_conditions() {
    // R(r) should equal (dr/dλ)^2 at the initial point
    // This validates that R(r) is computed consistently with initial velocities
    let r_s = 1.0;
    let a = 0.5;
    let r = 5.0;
    let theta = std::f64::consts::FRAC_PI_2;
    let e = 1.0;
    let l_z = 3.0;
    // Q for equatorial orbit (theta=PI/2, p_theta=0): Q = 0
    let q = 0.0;

    let r_potential = potential_r(r, r_s, a, e, l_z, q);
    assert!(r_potential >= 0.0, "R(r) must be non-negative for allowed motion");
}

#[test]
fn test_potential_derivatives_numerical() {
    // Compare analytical R'(r) against numerical finite difference
    let r_s = 1.0;
    let a = 0.5;
    let e = 1.0;
    let l_z = 3.0;
    let q = 1.0;
    let r = 5.0;
    let h = 1e-7;

    let dr_analytical = potential_r_derivative(r, r_s, a, e, l_z, q);
    let dr_numerical = (potential_r(r + h, r_s, a, e, l_z, q)
        - potential_r(r - h, r_s, a, e, l_z, q))
        / (2.0 * h);
    assert_abs_diff_eq!(dr_analytical, dr_numerical, epsilon = 1e-4);

    let theta = 1.2;
    let dth_analytical = potential_theta_derivative(theta, a, e, l_z, q);
    let dth_numerical = (potential_theta(theta + h, a, e, l_z, q)
        - potential_theta(theta - h, a, e, l_z, q))
        / (2.0 * h);
    assert_abs_diff_eq!(dth_analytical, dth_numerical, epsilon = 1e-4);
}
```

**Step 2: Run tests to verify they fail**

Run: `cargo test --lib kerr_bl::tests::test_potential 2>&1 | head -20`
Expected: Compilation error — functions not defined.

**Step 3: Implement potential functions and derivatives**

```rust
/// R(r) = [(r² + a²)E - aL_z]² - Δ[(L_z - aE)² + Q]
fn potential_r(r: f64, r_s: f64, a: f64, e: f64, l_z: f64, q: f64) -> f64 {
    let del = delta(r, r_s, a);
    let p_r = (r * r + a * a) * e - a * l_z;
    p_r * p_r - del * ((l_z - a * e).powi(2) + q)
}

/// R'(r) = 4rE[(r² + a²)E - aL_z] - (2r - r_s)[(L_z - aE)² + Q]
fn potential_r_derivative(r: f64, r_s: f64, a: f64, e: f64, l_z: f64, q: f64) -> f64 {
    let p_r = (r * r + a * a) * e - a * l_z;
    let carter_term = (l_z - a * e).powi(2) + q;
    4.0 * r * e * p_r - (2.0 * r - r_s) * carter_term
}

/// Θ(θ) = Q + a²E²cos²θ - L_z²cos²θ/sin²θ
fn potential_theta(theta: f64, a: f64, e: f64, l_z: f64, q: f64) -> f64 {
    let cos_t = theta.cos();
    let sin_t = theta.sin();
    q + a * a * e * e * cos_t * cos_t - l_z * l_z * cos_t * cos_t / (sin_t * sin_t)
}

/// Θ'(θ) = -2a²E²cosθsinθ + 2L_z²cosθ/sin³θ
fn potential_theta_derivative(theta: f64, a: f64, e: f64, l_z: f64, q: f64) -> f64 {
    let cos_t = theta.cos();
    let sin_t = theta.sin();
    -2.0 * a * a * e * e * cos_t * sin_t + 2.0 * l_z * l_z * cos_t / (sin_t.powi(3))
}
```

**Step 4: Run tests to verify they pass**

Run: `cargo test --lib kerr_bl::tests::test_potential 2>&1 | tail -5`
Expected: PASS

**Step 5: Commit**

```bash
git add src/geometry/kerr_bl.rs
git commit -m "feat: implement Kerr potential functions R(r), Θ(θ) and derivatives"
```

**Step 6: Write failing test for geodesic RHS**

```rust
#[test]
fn test_geodesic_rhs_structure() {
    // Verify the ODE RHS has the right structure:
    // ẏ[0] = dt/dλ, ẏ[1] = y[4], ẏ[2] = y[5], ẏ[3] = dφ/dλ,
    // ẏ[4] = R'(r)/2, ẏ[5] = Θ'(θ)/2, ẏ[6] = 0, ẏ[7] = 0
    let r_s = 1.0;
    let a = 0.5;
    let e = 1.0;
    let l_z = 3.0;
    let q = 1.0;
    let solver = KerrBLSolver { radius: r_s, a, e, l_z, q };

    let r = 5.0;
    let theta = 1.2;
    let v_r = 0.1;
    let v_theta = -0.05;
    let y = EquationOfMotionState::from_column_slice(&[0.0, r, theta, 0.0, v_r, v_theta, 0.0, 0.0]);

    let rhs = solver.geodesic(0.0, &y);

    // ẏ[1] = v_r = y[4]
    assert_abs_diff_eq!(rhs[1], v_r, epsilon = 1e-12);
    // ẏ[2] = v_theta = y[5]
    assert_abs_diff_eq!(rhs[2], v_theta, epsilon = 1e-12);
    // ẏ[4] = R'(r)/2
    assert_abs_diff_eq!(rhs[4], potential_r_derivative(r, r_s, a, e, l_z, q) / 2.0, epsilon = 1e-12);
    // ẏ[5] = Θ'(θ)/2
    assert_abs_diff_eq!(rhs[5], potential_theta_derivative(theta, a, e, l_z, q) / 2.0, epsilon = 1e-12);
    // ẏ[6] = 0, ẏ[7] = 0
    assert_abs_diff_eq!(rhs[6], 0.0);
    assert_abs_diff_eq!(rhs[7], 0.0);
}
```

**Step 7: Run test to verify it fails**

Run: `cargo test --lib kerr_bl::tests::test_geodesic_rhs 2>&1 | head -20`
Expected: Compilation error — `KerrBLSolver` not defined.

**Step 8: Implement KerrBLSolver and geodesic equations**

```rust
use crate::rendering::runge_kutta::OdeFunction;
use crate::rendering::scene::EquationOfMotionState;
use nalgebra::{Const, OVector};

struct KerrBLSolver {
    radius: f64,
    a: f64,
    e: f64,
    l_z: f64,
    q: f64,
}

impl OdeFunction<Const<8>> for KerrBLSolver {
    fn apply(&self, t: f64, y: &OVector<f64, Const<8>>) -> OVector<f64, Const<8>> {
        self.geodesic(t, y)
    }
}

impl HasCoordinateSystem for KerrBLSolver {
    fn coordinate_system(&self) -> CoordinateSystem {
        CoordinateSystem::BoyerLindquist { a: self.a }
    }
}

impl GeodesicSolver for KerrBLSolver {
    fn geodesic(&self, _: f64, y: &EquationOfMotionState) -> EquationOfMotionState {
        let r = y[1];
        let theta = y[2];
        let v_r = y[4];
        let v_theta = y[5];

        let del = delta(r, self.radius, self.a);
        let p_r = (r * r + self.a * self.a) * self.e - self.a * self.l_z;

        // dt/dλ = (r²+a²)/Δ * P_r + a(L_z - aE sin²θ)
        let sin_t = theta.sin();
        let sin2 = sin_t * sin_t;
        let dt = (r * r + self.a * self.a) / del * p_r + self.a * (self.l_z - self.a * self.e * sin2);

        // dφ/dλ = a/Δ * P_r + L_z/sin²θ - aE
        let dphi = self.a / del * p_r + self.l_z / sin2 - self.a * self.e;

        // d²r/dλ² = R'(r)/2
        let dv_r = potential_r_derivative(r, self.radius, self.a, self.e, self.l_z, self.q) / 2.0;

        // d²θ/dλ² = Θ'(θ)/2
        let dv_theta = potential_theta_derivative(theta, self.a, self.e, self.l_z, self.q) / 2.0;

        EquationOfMotionState::from_column_slice(&[dt, v_r, v_theta, dphi, dv_r, dv_theta, 0.0, 0.0])
    }

    fn create_initial_state(&self, ray: &Ray) -> EquationOfMotionState {
        // Will be implemented in Task 5
        todo!()
    }

    fn momentum_from_state(&self, y: &EquationOfMotionState) -> FourVector {
        // Will be implemented in Task 5
        todo!()
    }
}
```

**Step 9: Run tests to verify they pass**

Run: `cargo test --lib kerr_bl::tests::test_geodesic_rhs 2>&1 | tail -5`
Expected: PASS

**Step 10: Commit**

```bash
git add src/geometry/kerr_bl.rs
git commit -m "feat: implement KerrBLSolver with separated geodesic equations"
```

---

### Task 5: Implement Initial Conditions and Momentum Reconstruction

**Files:**
- Modify: `src/geometry/kerr_bl.rs`

This task implements `create_initial_state`, `momentum_from_state`, and `get_geodesic_solver`.

**Step 1: Write failing test for initial condition conversion**

```rust
#[test]
fn test_initial_conditions_constants_match_kerr() {
    use crate::geometry::kerr::Kerr;
    use crate::geometry::geometry::Geometry;
    use crate::rendering::camera::Camera;

    let radius = 1.0;
    let a = 0.5;
    let position = Point::new_cartesian(0.0, -10.0, 0.0, 2.0);

    let kerr = Kerr::new(radius, a, 1e-4);
    let kerr_bl = KerrBL::new(radius, a, 1e-4);

    // Create a camera ray using the Kerr geometry (Cartesian)
    let velocity = FourVector::new_cartesian(1.0, 0.0, 0.0, 0.0);
    let camera_kerr = Camera::new(
        position, velocity, std::f64::consts::FRAC_PI_2,
        11, 11, 0.0, 0.0, 0.0, &kerr,
    ).unwrap();
    let ray = camera_kerr.get_ray_for(5, 5);

    // Get constants from both geometries
    let kerr_constants = kerr.get_constants_of_motion(&ray.position, &ray.momentum);
    let kerr_e = kerr_constants.as_slice()[0].1;
    let kerr_lz = kerr_constants.as_slice()[1].1;

    // Get constants from KerrBL via get_geodesic_solver
    let solver = kerr_bl.get_geodesic_solver(&ray);
    // The solver should have computed E, L_z internally
    // Verify by creating initial state and checking consistency
    let state = solver.create_initial_state(&ray);
    let p = solver.momentum_from_state(&state);
    let position_bl = Point::new(
        state[0], state[1], state[2], state[3],
        CoordinateSystem::BoyerLindquist { a },
    );
    let bl_constants = kerr_bl.get_constants_of_motion(&position_bl, &p);
    let bl_e = bl_constants.as_slice()[0].1;
    let bl_lz = bl_constants.as_slice()[1].1;

    assert_abs_diff_eq!(kerr_e, bl_e, epsilon = 1e-6);
    assert_abs_diff_eq!(kerr_lz, bl_lz, epsilon = 1e-6);
}

#[test]
fn test_initial_null_condition() {
    let radius = 1.0;
    let a = 0.5;
    let position = Point::new_cartesian(0.0, -10.0, 0.0, 2.0);
    let kerr_bl = KerrBL::new(radius, a, 1e-4);

    // Use the existing Kerr geometry to create a valid null ray
    let kerr = crate::geometry::kerr::Kerr::new(radius, a, 1e-4);
    let velocity = FourVector::new_cartesian(1.0, 0.0, 0.0, 0.0);
    let camera = crate::rendering::camera::Camera::new(
        position, velocity, std::f64::consts::FRAC_PI_2,
        11, 11, 0.0, 0.0, 0.0, &kerr,
    ).unwrap();
    let ray = camera.get_ray_for(5, 5);

    let solver = kerr_bl.get_geodesic_solver(&ray);
    let state = solver.create_initial_state(&ray);
    let p = solver.momentum_from_state(&state);
    let pos = Point::new(
        state[0], state[1], state[2], state[3],
        CoordinateSystem::BoyerLindquist { a },
    );

    let null_check = kerr_bl.inner_product(&pos, &p, &p);
    assert_abs_diff_eq!(null_check, 0.0, epsilon = 1e-8);
}
```

**Step 2: Run tests to verify they fail**

Run: `cargo test --lib kerr_bl::tests::test_initial 2>&1 | head -20`
Expected: FAIL — `todo!()` panics.

**Step 3: Implement the conversions**

Implement `get_geodesic_solver` on `KerrBL`:
```rust
fn get_geodesic_solver(&self, ray: &Ray) -> Box<dyn GeodesicSolver> {
    // Convert Cartesian position to BL
    let (x, y, z) = (ray.position[1], ray.position[2], ray.position[3]);
    let r_sqr = compute_r_sqr(self.a, x, y, z);
    let r = r_sqr.sqrt();
    let theta = if r == 0.0 { 0.0 } else { (z / r).max(-1.0).min(1.0).acos() };
    let phi = (r * y - self.a * x).atan2(r * x + self.a * y);

    // Compute Jacobian ∂x^Cart/∂x^BL and its inverse
    let jacobian = self.jacobian_bl_to_cartesian(r, theta, phi);
    let jacobian_inv = jacobian.try_inverse().expect("Jacobian should be invertible");

    // Convert contravariant momentum: p^BL = J^{-1} p^Cart
    let p_bl_contra = jacobian_inv * ray.momentum.vector;

    // Lower with BL metric: p_μ = g_μν p^ν
    let g = metric_bl(self.radius, self.a, r, theta);
    let p_bl_cov = g * p_bl_contra;

    // Extract constants of motion
    let e = -p_bl_cov[0];
    let l_z = p_bl_cov[3];
    let p_theta = p_bl_cov[2];
    let cos_t = theta.cos();
    let sin_t = theta.sin();
    let q = p_theta * p_theta + cos_t * cos_t * (l_z * l_z / (sin_t * sin_t) - self.a * self.a * e * e);

    Box::new(KerrBLSolver {
        radius: self.radius,
        a: self.a,
        e,
        l_z,
        q,
    })
}
```

Add the Jacobian helper to `KerrBL`:
```rust
impl KerrBL {
    fn jacobian_bl_to_cartesian(&self, r: f64, theta: f64, phi: f64) -> Matrix4<f64> {
        let a = self.a;
        let (st, ct) = (theta.sin(), theta.cos());
        let (sp, cp) = (phi.sin(), phi.cos());

        #[rustfmt::skip]
        let data = [
            1.0, 0.0, 0.0, 0.0,
            0.0, st * cp, (r * cp - a * sp) * ct, (-r * sp - a * cp) * st,
            0.0, st * sp, (r * sp + a * cp) * ct, (r * cp - a * sp) * st,
            0.0, ct, -r * st, 0.0,
        ];

        Matrix4::from_row_slice(&data)
    }
}
```

Reuse `compute_r_sqr` from existing Kerr module or duplicate it:
```rust
fn compute_r_sqr(a: f64, x: f64, y: f64, z: f64) -> f64 {
    let rho_sqr = x * x + y * y + z * z;
    0.5 * (rho_sqr - a * a + ((rho_sqr - a * a).powi(2) + 4.0 * a * a * z * z).sqrt())
}
```

Implement `create_initial_state`:
```rust
fn create_initial_state(&self, ray: &Ray) -> EquationOfMotionState {
    let (x, y, z) = (ray.position[1], ray.position[2], ray.position[3]);
    let r_sqr = compute_r_sqr(self.a, x, y, z);
    let r = r_sqr.sqrt();
    let theta = if r == 0.0 { 0.0 } else { (z / r).max(-1.0).min(1.0).acos() };
    let phi = (r * y - self.a * x).atan2(r * x + self.a * y);
    let t = ray.position[0];

    // Compute initial dr/dλ and dθ/dλ from the potentials
    let r_pot = potential_r(r, self.radius, self.a, self.e, self.l_z, self.q);
    let th_pot = potential_theta(theta, self.a, self.e, self.l_z, self.q);

    // Determine signs from the initial momentum direction
    // Convert momentum to BL to get the signs of p^r and p^θ
    let jacobian = KerrBL::new(self.radius, self.a, 0.0)
        .jacobian_bl_to_cartesian(r, theta, phi);
    let jacobian_inv = jacobian.try_inverse().expect("Jacobian should be invertible");
    let p_bl_contra = jacobian_inv * ray.momentum.vector;

    let sign_r = p_bl_contra[1].signum();
    let sign_theta = p_bl_contra[2].signum();

    let v_r = sign_r * r_pot.max(0.0).sqrt();
    let v_theta = sign_theta * th_pot.max(0.0).sqrt();

    EquationOfMotionState::from_column_slice(&[t, r, theta, phi, v_r, v_theta, 0.0, 0.0])
}
```

Implement `momentum_from_state`:
```rust
fn momentum_from_state(&self, y: &EquationOfMotionState) -> FourVector {
    let r = y[1];
    let theta = y[2];
    let v_r = y[4];
    let v_theta = y[5];

    let del = delta(r, self.radius, self.a);
    let sig = sigma(r, self.a, theta);
    let sin2 = theta.sin().powi(2);
    let p_r_coeff = (r * r + self.a * self.a) * self.e - self.a * self.l_z;

    let dt_dlambda = (r * r + self.a * self.a) / del * p_r_coeff
        + self.a * (self.l_z - self.a * self.e * sin2);
    let dphi_dlambda = self.a / del * p_r_coeff + self.l_z / sin2 - self.a * self.e;

    // Convert Mino-time velocities to affine-parameter momentum: p^μ = (1/Σ) dx^μ/dλ
    FourVector::new_boyer_lindquist(
        self.a,
        dt_dlambda / sig,
        v_r / sig,
        v_theta / sig,
        dphi_dlambda / sig,
    )
}
```

**Step 4: Run tests to verify they pass**

Run: `cargo test --lib kerr_bl 2>&1 | tail -5`
Expected: PASS

**Step 5: Commit**

```bash
git add src/geometry/kerr_bl.rs
git commit -m "feat: implement initial conditions and momentum reconstruction for KerrBL"
```

---

### Task 6: Implement Tetrad, Lorentz Transformation, and SupportQuantities

**Files:**
- Modify: `src/geometry/kerr_bl.rs`

**Step 1: Write failing test for tetrad orthonormality**

```rust
#[test]
fn test_tetrad_orthonormal() {
    let a = 0.5;
    let kerr_bl = KerrBL::new(1.0, a, 1e-4);
    let position = Point::new(0.0, 5.0, 1.2, 0.8, CoordinateSystem::BoyerLindquist { a });
    let tetrad = kerr_bl.get_tetrad_at(&position);

    assert_abs_diff_eq!(kerr_bl.inner_product(&position, &tetrad.t, &tetrad.t), -1.0, epsilon = 1e-10);
    assert_abs_diff_eq!(kerr_bl.inner_product(&position, &tetrad.x, &tetrad.x), 1.0, epsilon = 1e-10);
    assert_abs_diff_eq!(kerr_bl.inner_product(&position, &tetrad.y, &tetrad.y), 1.0, epsilon = 1e-10);
    assert_abs_diff_eq!(kerr_bl.inner_product(&position, &tetrad.z, &tetrad.z), 1.0, epsilon = 1e-10);

    assert_abs_diff_eq!(kerr_bl.inner_product(&position, &tetrad.t, &tetrad.x), 0.0, epsilon = 1e-10);
    assert_abs_diff_eq!(kerr_bl.inner_product(&position, &tetrad.t, &tetrad.y), 0.0, epsilon = 1e-10);
    assert_abs_diff_eq!(kerr_bl.inner_product(&position, &tetrad.t, &tetrad.z), 0.0, epsilon = 1e-10);
    assert_abs_diff_eq!(kerr_bl.inner_product(&position, &tetrad.x, &tetrad.y), 0.0, epsilon = 1e-10);
    assert_abs_diff_eq!(kerr_bl.inner_product(&position, &tetrad.x, &tetrad.z), 0.0, epsilon = 1e-10);
    assert_abs_diff_eq!(kerr_bl.inner_product(&position, &tetrad.y, &tetrad.z), 0.0, epsilon = 1e-10);
}
```

**Step 2: Run to verify it fails**

Run: `cargo test --lib kerr_bl::tests::test_tetrad 2>&1 | head -20`
Expected: FAIL — `todo!()` panics.

**Step 3: Implement the ZAMO tetrad and Gram-Schmidt orthonormalization**

The ZAMO (zero angular momentum observer) tetrad in BL coordinates:
```rust
use crate::geometry::gram_schmidt::gram_schmidt;

fn get_tetrad_at(&self, position: &Point) -> Tetrad {
    let r = position[1];
    let theta = position[2];
    let a = self.a;
    let r_s = self.radius;

    let sig = sigma(r, a, theta);
    let del = delta(r, r_s, a);
    let sin_t = theta.sin();
    let sin2 = sin_t * sin_t;

    // ZAMO angular velocity: ω = -g_tφ/g_φφ
    let g_tph = -a * r_s * r * sin2 / sig;
    let g_phph = (r * r + a * a + a * a * r_s * r * sin2 / sig) * sin2;
    let omega = -g_tph / g_phph;

    // ZAMO four-velocity: u^μ = u^t (1, 0, 0, ω)
    // Normalization: g_μν u^μ u^ν = -1
    let g_tt = -(1.0 - r_s * r / sig);
    let ut_sq = -1.0 / (g_tt + 2.0 * g_tph * omega + g_phph * omega * omega);
    let ut = ut_sq.sqrt();

    let coord_sys = CoordinateSystem::BoyerLindquist { a };
    let e_t = FourVector::new(ut, 0.0, 0.0, ut * omega, coord_sys);
    let e_r = FourVector::new(0.0, 1.0, 0.0, 0.0, coord_sys);
    let e_th = FourVector::new(0.0, 0.0, 1.0, 0.0, coord_sys);
    let e_ph = FourVector::new(0.0, 0.0, 0.0, 1.0, coord_sys);

    let basis = gram_schmidt(self, position, &[e_t, e_r, e_th, e_ph]);
    Tetrad::new(*position, basis[0], basis[1], basis[2], basis[3])
}
```

Implement `lorentz_transformation` using the full BL metric (non-diagonal). Follow the same pattern as the existing Kerr implementation but with the BL metric:
```rust
fn lorentz_transformation(&self, position: &Point, velocity: &FourVector) -> Matrix4<f64> {
    let r = position[1];
    let theta = position[2];
    let g = metric_bl(self.radius, self.a, r, theta);
    let tetrad_t = self.get_tetrad_at(position).t;

    let gamma = -(tetrad_t.vector.transpose() * g * velocity.vector)[(0, 0)];

    let uv = tetrad_t.vector + velocity.vector;
    let uv_lower = g * uv;

    let mut matrix = Matrix4::zeros();
    for mu in 0..4 {
        for nu in 0..4 {
            let mut res = if mu == nu { 1.0 } else { 0.0 };
            res += uv[mu] * uv_lower[nu] / (1.0 + gamma);
            res -= 2.0 * (g * tetrad_t.vector)[nu] * velocity.vector[mu];
            matrix[(mu, nu)] = res;
        }
    }
    matrix
}
```

Implement `SupportQuantities`:
```rust
impl SupportQuantities for KerrBL {
    fn get_stationary_velocity_at(&self, position: &Point) -> FourVector {
        // ZAMO velocity (same as tetrad e_t)
        let r = position[1];
        let theta = position[2];
        let sig = sigma(r, self.a, theta);
        let del = delta(r, self.radius, self.a);
        let sin2 = theta.sin().powi(2);
        let g_tph = -self.a * self.radius * r * sin2 / sig;
        let g_phph = (r * r + self.a * self.a + self.a * self.a * self.radius * r * sin2 / sig) * sin2;
        let g_tt = -(1.0 - self.radius * r / sig);
        let omega = -g_tph / g_phph;
        let ut = (-1.0 / (g_tt + 2.0 * g_tph * omega + g_phph * omega * omega)).sqrt();
        FourVector::new_boyer_lindquist(self.a, ut, 0.0, 0.0, ut * omega)
    }

    fn get_circular_orbit_velocity_at(&self, position: &Point) -> Result<FourVector, RaytracerError> {
        let r = position[1];
        let m = 0.5 * self.radius;
        let omega = m.sqrt() / (r.powf(1.5) + self.a * m.sqrt());

        let theta = position[2];
        let sig = sigma(r, self.a, theta);
        let sin2 = theta.sin().powi(2);
        let g_tt = -(1.0 - self.radius * r / sig);
        let g_tph = -self.a * self.radius * r * sin2 / sig;
        let g_phph = (r * r + self.a * self.a + self.a * self.a * self.radius * r * sin2 / sig) * sin2;

        let ut_pre = g_tt + 2.0 * omega * g_tph + omega * omega * g_phph;
        if ut_pre >= 0.0 {
            return Err(RaytracerError::NoCircularOrbitPossible);
        }
        let ut = (-ut_pre).sqrt().recip();
        let uphi = omega * ut;
        Ok(FourVector::new_boyer_lindquist(self.a, ut, 0.0, 0.0, uphi))
    }

    fn get_temperature_computer(
        &self,
        temperature: f64,
        _inner_radius: f64,
        outer_radius: f64,
    ) -> Result<Box<dyn TemperatureComputer>, RaytracerError> {
        use crate::rendering::temperature::KerrTemperatureComputer;
        Ok(Box::new(KerrTemperatureComputer::new(
            temperature, outer_radius, self.a, self.radius,
        )?))
    }
}
```

Implement `get_constants_of_motion` and `closed_orbit`:
```rust
fn get_constants_of_motion(&self, position: &Point, momentum: &FourVector) -> ConstantsOfMotion {
    let r = position[1];
    let theta = position[2];
    let g = metric_bl(self.radius, self.a, r, theta);
    let p_cov = g * momentum.vector;
    let e = -p_cov[0];
    let l_z = p_cov[3];

    let mut constants = ConstantsOfMotion::default();
    constants.push("E", e);
    constants.push("L_z", l_z);
    constants
}

fn closed_orbit(&self, position: &Point, step_index: usize, max_steps: usize) -> bool {
    let r = position[1];
    if step_index == max_steps - 1 && r <= self.radius {
        return true;
    }
    false
}
```

**Step 4: Run tests to verify they pass**

Run: `cargo test --lib kerr_bl 2>&1 | tail -5`
Expected: PASS

**Step 5: Commit**

```bash
git add src/geometry/kerr_bl.rs
git commit -m "feat: implement tetrad, Lorentz transformation, and SupportQuantities for KerrBL"
```

---

### Task 7: Wire KerrBL into Configuration and CLI

**Files:**
- Modify: `src/configuration.rs`
- Create: `src/cli/kerr_bl.rs`
- Modify: `src/cli/mod.rs`

**Step 1: Add GeometryType::KerrBL variant**

In `src/configuration.rs`, add to the enum:
```rust
KerrBL {
    radius: f64,
    a: f64,
    horizon_epsilon: f64,
},
```

Add match arm in `get_renderable_geometry`:
```rust
GeometryType::KerrBL {
    radius,
    a,
    horizon_epsilon,
} => Box::new(KerrBL::new(*radius, *a, *horizon_epsilon)),
```

Add the import at the top of `configuration.rs`:
```rust
use crate::geometry::kerr_bl::KerrBL;
```

**Step 2: Create CLI module**

Create `src/cli/kerr_bl.rs` following the same pattern as `src/cli/kerr.rs`. The key difference is that `create_scene_internal` needs to compute the ZAMO velocity at the camera position in BL coordinates:

```rust
use crate::cli::cli::GlobalOpts;
use crate::cli::shared::{assert_future_directed, create_scene, integrate_and_save_ray, render};
use crate::configuration::RenderConfig;
use crate::geometry::four_vector::FourVector;
use crate::geometry::geometry::{Geometry, InnerProduct, RenderableGeometry, SupportQuantities};
use crate::geometry::kerr_bl::KerrBL;
use crate::geometry::point::{CoordinateSystem, Point};
use crate::rendering::raytracer;
use crate::rendering::raytracer::RaytracerError;
use crate::rendering::scene::Scene;
use log::info;
use std::io::Write;

fn create_scene_internal<'a>(
    geometry: &'a KerrBL,
    opts: GlobalOpts,
    config: &RenderConfig,
    camera_position: Point,
) -> Result<Scene<'a, KerrBL>, RaytracerError> {
    let momentum = geometry.get_stationary_velocity_at(&camera_position);
    assert_future_directed("KerrBL camera four-velocity", geometry, &camera_position, &momentum)?;
    create_scene(geometry, camera_position, momentum, opts, config.clone())
}

impl RenderableGeometry for KerrBL {
    fn render(
        &self,
        opts: GlobalOpts,
        config: RenderConfig,
        camera_position: Point,
        filename: String,
        from_row: Option<u32>,
        from_col: Option<u32>,
        to_row: Option<u32>,
        to_col: Option<u32>,
    ) -> Result<(), RaytracerError> {
        let camera_position_bl = self.cartesian_to_bl(&camera_position);
        let scene = create_scene_internal(self, opts, &config, camera_position_bl)?;
        render(scene, filename, config.color_normalization, from_row, from_col, to_row, to_col)
    }

    fn render_ray(
        &self,
        row: i64,
        col: i64,
        opts: GlobalOpts,
        config: RenderConfig,
        camera_position: Point,
        write: &mut dyn Write,
    ) -> Result<(), RaytracerError> {
        let camera_position_bl = self.cartesian_to_bl(&camera_position);
        let scene = create_scene_internal(self, opts, &config, camera_position_bl)?;
        let raytracer = raytracer::Raytracer::new(scene, config.color_normalization);
        let (integrated_ray, stop_reason) = raytracer.integrate_ray_at_point(row, col)?;
        info!("Stop reason: {:?}", stop_reason);
        integrated_ray.save(write)?;
        Ok(())
    }

    fn render_ray_at(
        &self,
        position: Point,
        direction: FourVector,
        opts: GlobalOpts,
        write: &mut dyn Write,
    ) -> Result<(), RaytracerError> {
        let position_bl = self.cartesian_to_bl(&position);
        let tetrad = self.get_tetrad_at(&position_bl);
        info!("Tetrad at position {:?}: {}", position_bl, tetrad);

        let space_part =
            tetrad.x * direction[1] + tetrad.y * direction[2] + tetrad.z * direction[3];
        let norm_space_part = self
            .inner_product(&position_bl, &space_part, &space_part)
            .sqrt();

        let momentum = tetrad.t * 1.0
            + tetrad.x * direction[1] / norm_space_part
            + tetrad.y * direction[2] / norm_space_part
            + tetrad.z * direction[3] / norm_space_part;
        assert_future_directed("KerrBL render_ray_at momentum", self, &position_bl, &momentum)?;

        integrate_and_save_ray(self, position_bl, momentum, opts, write)
    }
}
```

Add a `cartesian_to_bl` helper to `KerrBL`:
```rust
pub fn cartesian_to_bl(&self, position: &Point) -> Point {
    let x = position[1];
    let y = position[2];
    let z = position[3];
    let r_sqr = compute_r_sqr(self.a, x, y, z);
    let r = r_sqr.sqrt();
    let theta = if r == 0.0 { 0.0 } else { (z / r).max(-1.0).min(1.0).acos() };
    let phi = (r * y - self.a * x).atan2(r * x + self.a * y);
    Point::new(position[0], r, theta, phi, CoordinateSystem::BoyerLindquist { a: self.a })
}
```

Register in `src/cli/mod.rs`:
```rust
pub mod kerr_bl;
```

**Step 3: Run full build**

Run: `cargo build 2>&1 | tail -10`
Expected: Compiles successfully.

**Step 4: Commit**

```bash
git add src/configuration.rs src/cli/kerr_bl.rs src/cli/mod.rs src/geometry/kerr_bl.rs
git commit -m "feat: wire KerrBL into configuration and CLI"
```

---

### Task 8: Cross-Geometry Validation Tests

**Files:**
- Modify: `src/geometry/kerr_bl.rs` (test module)

These are the critical validation tests that compare KerrBL against the existing Kerr implementation.

**Step 1: Write trajectory comparison test**

```rust
#[test]
fn test_trajectory_agreement_with_kerr() {
    use crate::geometry::kerr::Kerr;
    use crate::rendering::camera::Camera;
    use crate::rendering::scene;
    use crate::rendering::scene::Scene;

    let radius = 1.0;
    let a = 0.3;
    let position_cart = Point::new_cartesian(0.0, -10.0, 0.0, 2.0);

    // Set up Kerr (Cartesian)
    let kerr = Kerr::new(radius, a, 1e-5);
    let velocity_cart = FourVector::new_cartesian(1.0, 0.0, 0.0, 0.0);
    let camera_kerr = Camera::new(
        position_cart, velocity_cart, std::f64::consts::FRAC_PI_2,
        11, 11, 0.0, 0.0, 0.0, &kerr,
    ).unwrap();
    let scene_kerr: Scene<Kerr> =
        scene::test_scene::create_scene_with_camera(1.0, 2.0, 7.0, &kerr, camera_kerr, 1e-6)
            .unwrap();
    let ray_kerr = scene_kerr.camera.get_ray_for(5, 8);

    // Set up KerrBL (Boyer-Lindquist)
    let kerr_bl = KerrBL::new(radius, a, 1e-5);
    let position_bl = kerr_bl.cartesian_to_bl(&position_cart);
    let velocity_bl = kerr_bl.get_stationary_velocity_at(&position_bl);
    let camera_bl = Camera::new(
        position_bl, velocity_bl, std::f64::consts::FRAC_PI_2,
        11, 11, 0.0, 0.0, 0.0, &kerr_bl,
    ).unwrap();
    let scene_bl: Scene<KerrBL> =
        scene::test_scene::create_scene_with_camera(1.0, 2.0, 7.0, &kerr_bl, camera_bl, 1e-6)
            .unwrap();
    let ray_bl = scene_bl.camera.get_ray_for(5, 8);

    // Integrate both
    let (traj_kerr, stop_kerr) = scene_kerr.integrator.integrate(&ray_kerr).unwrap();
    let (traj_bl, stop_bl) = scene_bl.integrator.integrate(&ray_bl).unwrap();

    // Same stop reason
    assert_eq!(stop_kerr, stop_bl);

    // Compare Cartesian positions at similar trajectory points
    // Note: different parameterizations mean different step counts, so compare
    // start and end positions rather than step-by-step
    let first_kerr = &traj_kerr[0];
    let first_bl = &traj_bl[0];
    let first_kerr_cart = first_kerr.x.get_spatial_vector_cartesian();
    let first_bl_cart = first_bl.x.get_spatial_vector_cartesian();
    assert_abs_diff_eq!(first_kerr_cart, first_bl_cart, epsilon = 1e-4);

    let last_kerr = traj_kerr.last().unwrap();
    let last_bl = traj_bl.last().unwrap();
    let last_kerr_cart = last_kerr.x.get_spatial_vector_cartesian();
    let last_bl_cart = last_bl.x.get_spatial_vector_cartesian();
    // Final positions should agree to within a reasonable tolerance
    let distance = (last_kerr_cart - last_bl_cart).norm();
    assert!(
        distance < 1.0,
        "Final positions differ by {} (should be < 1.0)",
        distance
    );
}
```

**Step 2: Write constants of motion conservation test**

```rust
#[test]
fn test_constants_of_motion_conservation() {
    use crate::geometry::kerr::Kerr;
    use crate::rendering::camera::Camera;
    use crate::rendering::scene;

    let radius = 1.0;
    let a = 0.4;
    let kerr_bl = KerrBL::new(radius, a, 1e-5);
    let position_bl = kerr_bl.cartesian_to_bl(&Point::new_cartesian(0.0, -10.0, 0.0, 2.0));
    let velocity_bl = kerr_bl.get_stationary_velocity_at(&position_bl);
    let camera = Camera::new(
        position_bl, velocity_bl, std::f64::consts::FRAC_PI_2,
        11, 11, 0.0, 0.0, 0.0, &kerr_bl,
    ).unwrap();
    let scene: scene::Scene<KerrBL> =
        scene::test_scene::create_scene_with_camera(1.0, 2.0, 7.0, &kerr_bl, camera, 1e-6)
            .unwrap();
    let ray = scene.camera.get_ray_for(3, 7);

    let (trajectory, _) = scene.integrator.integrate(&ray).unwrap();
    assert!(trajectory.len() > 10, "Trajectory too short for meaningful test");

    let initial = &trajectory[0];
    let initial_constants = kerr_bl.get_constants_of_motion(&initial.x, &initial.p);
    let e_init = initial_constants.as_slice()[0].1;
    let lz_init = initial_constants.as_slice()[1].1;

    for step in &trajectory[1..] {
        let constants = kerr_bl.get_constants_of_motion(&step.x, &step.p);
        let e = constants.as_slice()[0].1;
        let lz = constants.as_slice()[1].1;

        let e_drift = if e_init.abs() > 1e-12 {
            (e - e_init).abs() / e_init.abs()
        } else {
            (e - e_init).abs()
        };
        let lz_drift = if lz_init.abs() > 1e-12 {
            (lz - lz_init).abs() / lz_init.abs()
        } else {
            (lz - lz_init).abs()
        };

        assert!(e_drift < 1e-4, "Energy drift {:.3e} at step {}", e_drift, step.step);
        assert!(lz_drift < 1e-4, "L_z drift {:.3e} at step {}", lz_drift, step.step);
    }
}
```

**Step 3: Write null condition preservation test**

```rust
#[test]
fn test_null_condition_preserved() {
    use crate::rendering::camera::Camera;
    use crate::rendering::scene;

    let radius = 1.0;
    let a = 0.4;
    let kerr_bl = KerrBL::new(radius, a, 1e-5);
    let position_bl = kerr_bl.cartesian_to_bl(&Point::new_cartesian(0.0, -10.0, 0.0, 2.0));
    let velocity_bl = kerr_bl.get_stationary_velocity_at(&position_bl);
    let camera = Camera::new(
        position_bl, velocity_bl, std::f64::consts::FRAC_PI_2,
        11, 11, 0.0, 0.0, 0.0, &kerr_bl,
    ).unwrap();
    let scene: scene::Scene<KerrBL> =
        scene::test_scene::create_scene_with_camera(1.0, 2.0, 7.0, &kerr_bl, camera, 1e-6)
            .unwrap();
    let ray = scene.camera.get_ray_for(5, 5);

    let (trajectory, _) = scene.integrator.integrate(&ray).unwrap();

    for step in &trajectory {
        let k_dot_k = kerr_bl.inner_product(&step.x, &step.p, &step.p);
        assert!(
            k_dot_k.abs() < 1e-4,
            "Null condition violated: k.k = {:.3e} at step {}",
            k_dot_k,
            step.step
        );
    }
}
```

**Step 4: Write Schwarzschild limit test**

```rust
#[test]
fn test_schwarzschild_limit() {
    use crate::geometry::schwarzschild::Schwarzschild;
    use crate::rendering::camera::Camera;
    use crate::rendering::scene;
    use crate::geometry::spherical_coordinates_helper::cartesian_to_spherical;

    let radius = 1.0;
    let a = 0.0; // Schwarzschild limit
    let position_cart = Point::new_cartesian(0.0, -10.0, 0.0, 2.0);

    // KerrBL with a=0
    let kerr_bl = KerrBL::new(radius, a, 1e-5);
    let position_bl = kerr_bl.cartesian_to_bl(&position_cart);
    let velocity_bl = kerr_bl.get_stationary_velocity_at(&position_bl);
    let camera_bl = Camera::new(
        position_bl, velocity_bl, std::f64::consts::FRAC_PI_2,
        11, 11, 0.0, 0.0, 0.0, &kerr_bl,
    ).unwrap();
    let scene_bl: scene::Scene<KerrBL> =
        scene::test_scene::create_scene_with_camera(1.0, 2.0, 7.0, &kerr_bl, camera_bl, 1e-6)
            .unwrap();
    let ray_bl = scene_bl.camera.get_ray_for(5, 8);
    let (traj_bl, stop_bl) = scene_bl.integrator.integrate(&ray_bl).unwrap();

    // Schwarzschild
    let schwarzschild = Schwarzschild::new(radius, 1e-5);
    let position_sph = cartesian_to_spherical(&position_cart);
    let r = position_sph[1];
    let a_factor = 1.0 - radius / r;
    let velocity_sph = FourVector::new_spherical(a_factor.sqrt().recip(), 0.0, 0.0, 0.0);
    let camera_sch = Camera::new(
        position_sph, velocity_sph, std::f64::consts::FRAC_PI_2,
        11, 11, 0.0, 0.0, 0.0, &schwarzschild,
    ).unwrap();
    let scene_sch: scene::Scene<Schwarzschild> =
        scene::test_scene::create_scene_with_camera(1.0, 2.0, 7.0, &schwarzschild, camera_sch, 1e-6)
            .unwrap();
    let ray_sch = scene_sch.camera.get_ray_for(5, 8);
    let (traj_sch, stop_sch) = scene_sch.integrator.integrate(&ray_sch).unwrap();

    assert_eq!(stop_bl, stop_sch);

    // Compare final positions in Cartesian
    let last_bl = traj_bl.last().unwrap().x.get_spatial_vector_cartesian();
    let last_sch = traj_sch.last().unwrap().x.get_spatial_vector_cartesian();
    let distance = (last_bl - last_sch).norm();
    assert!(distance < 0.5, "Schwarzschild limit: positions differ by {}", distance);
}
```

**Step 5: Run all tests**

Run: `cargo test --lib kerr_bl 2>&1 | tail -10`
Expected: PASS

**Step 6: Commit**

```bash
git add src/geometry/kerr_bl.rs
git commit -m "test: add cross-geometry validation tests for KerrBL"
```

---

### Task 9: Create TOML Scene Definition and Integration Test

**Files:**
- Create: `scene-definitions/kerr-bl.toml`
- Modify: `src/geometry/kerr_bl.rs` (optional: add render integration test)

**Step 1: Create TOML scene definition**

Copy `scene-definitions/kerr.toml` and change the geometry type:
```toml
[geometry_type]
type = "KerrBL"
radius = 1.0
a = 0.5
horizon_epsilon = 0.001
```
Keep all other settings (objects, textures, etc.) the same as `kerr.toml`.

**Step 2: Verify it parses**

Run: `cargo run --release -- --width=50 --height=50 --max-steps=1000 --camera-position=-10.0,0.0,0.5 --theta=1.57 --psi=-1.57 --phi=0 --config-file scene-definitions/kerr-bl.toml render --filename=/tmp/test-kerr-bl.png 2>&1 | tail -10`
Expected: Renders successfully (small image for speed).

**Step 3: Commit**

```bash
git add scene-definitions/kerr-bl.toml
git commit -m "feat: add KerrBL scene definition"
```

---

### Task 10: Run Full Test Suite and Final Verification

**Step 1: Run all tests**

Run: `cargo test 2>&1 | tail -20`
Expected: All tests pass, no regressions.

**Step 2: Run clippy**

Run: `cargo clippy 2>&1 | tail -20`
Expected: No new warnings.

**Step 3: Final commit if any cleanup needed**

```bash
git add -A
git commit -m "chore: final cleanup for KerrBL implementation"
```

# KerrBL: Completely Integrable Kerr Geodesics in Boyer-Lindquist Coordinates

## Overview

Add a new geometry `KerrBL` that implements Kerr geodesics using the complete integrability of the system. Instead of numerically integrating coupled 2nd-order geodesic equations with finite-differenced metric derivatives (as the current `Kerr` does in Kerr-Schild Cartesian coordinates), `KerrBL` works in Boyer-Lindquist coordinates and exploits the Carter constant to separate the equations of motion.

The separated 2nd-order form (d²r/dλ² = R'(r)/2, d²θ/dλ² = Θ'(θ)/2 in Mino time) naturally handles turning points without sign tracking, uses purely analytical derivatives, and computes the constants of motion (E, L_z, Q) only once per ray.

## Architecture

### New File

`src/geometry/kerr_bl.rs` containing:

- **`KerrBL`** struct — implements `Geometry`, `InnerProduct`, `HasCoordinateSystem`, `Signature`, `SupportQuantities`
- **`KerrBLSolver`** struct — implements `GeodesicSolver`, `OdeFunction<Const<8>>`, `HasCoordinateSystem`; stores constants of motion (E, L_z, Q)

### Wiring

- New variant `GeometryType::KerrBL { radius, a, horizon_epsilon }` in `configuration.rs`
- Register `pub mod kerr_bl;` in `geometry/mod.rs`
- Add TOML/CLI support for the new geometry type

### Coordinate System

`CoordinateSystem::Spherical` — same as Schwarzschild. Scene objects receive positions converted to Cartesian via `to_cartesian()` for intersection testing (already handled by the existing pipeline).

## Physics

### Metric

Kerr in Boyer-Lindquist with signature (-,+,+,+):

- Sigma = r^2 + a^2 cos^2(theta)
- Delta = r^2 - r_s r + a^2, where r_s = 2M
- g_tt = -(1 - r_s r / Sigma)
- g_rr = Sigma / Delta
- g_theta_theta = Sigma
- g_phi_phi = (r^2 + a^2 + a^2 r_s r sin^2(theta) / Sigma) sin^2(theta)
- g_t_phi = -a r_s r sin^2(theta) / Sigma (cross term)

### Constants of Motion

Computed once from the initial ray in `get_geodesic_solver(ray)`:

- E = -p_t (energy, from time-translation symmetry)
- L_z = p_phi (angular momentum, from axial symmetry)
- Q = p_theta^2 + cos^2(theta) (L_z^2 / sin^2(theta) - a^2 E^2) (Carter constant)

### Potentials

- R(r) = [(r^2 + a^2)E - a L_z]^2 - Delta [(L_z - aE)^2 + Q]
- Theta(theta) = Q + a^2 E^2 cos^2(theta) - L_z^2 cos^2(theta) / sin^2(theta)

### Separated Geodesic Equations (Mino time)

Mino time lambda is defined by d(tau) = Sigma d(lambda), where tau is affine parameter.

**2nd-order ODEs** (evolved by RKF45):

- d^2 r / d(lambda)^2 = R'(r) / 2
- d^2 theta / d(lambda)^2 = Theta'(theta) / 2

Where the analytical derivatives are:

- R'(r) = 4rE [(r^2+a^2)E - a L_z] - (2r - r_s) [(L_z - aE)^2 + Q]
- Theta'(theta) = -2 a^2 E^2 cos(theta) sin(theta) + 2 L_z^2 cos(theta) / sin^3(theta)

**Algebraic equations** (computed at each step, not evolved):

- dt/d(lambda) = (r^2+a^2)/Delta * [(r^2+a^2)E - a L_z] + a(L_z - aE sin^2(theta))
- d(phi)/d(lambda) = a/Delta * [(r^2+a^2)E - a L_z] + L_z/sin^2(theta) - aE

Turning points (R=0 or Theta=0) are handled naturally by the 2nd-order form without explicit sign tracking.

### State Vector

8D, fits `OdeFunction<Const<8>>`:

| Index | Meaning       | ODE RHS                              |
|-------|---------------|--------------------------------------|
| y[0]  | t             | T_r(r) + T_theta(theta) (algebraic)  |
| y[1]  | r             | y[4]                                 |
| y[2]  | theta         | y[5]                                 |
| y[3]  | phi           | Phi_r(r) + Phi_theta(theta) (algebraic) |
| y[4]  | dr/d(lambda)  | R'(r) / 2                            |
| y[5]  | d(theta)/d(lambda) | Theta'(theta) / 2               |
| y[6]  | unused        | 0                                    |
| y[7]  | unused        | 0                                    |

### Momentum Reconstruction

In `momentum_from_state`, convert from Mino-time velocities to affine-parameter momentum:

p^mu = (1 / Sigma) * dx^mu / d(lambda)

where dt/d(lambda) and d(phi)/d(lambda) are computed algebraically from E, L_z, r, theta.

## Initial Conditions

`create_initial_state(ray)` in `KerrBLSolver`:

1. **Position conversion** (Cartesian to BL):
   - r from `compute_r_sqr` (same formula as existing Kerr)
   - theta = acos(z / r)
   - phi = atan2(r y - a x, r x + a y)

2. **Momentum conversion** (Cartesian to BL):
   - Compute contravariant BL momentum via inverse Jacobian
   - Lower with BL metric to get covariant components p_t, p_phi

3. **Extract constants**: E = -p_t, L_z = p_phi, Q from p_theta

4. **Set initial velocities**:
   - dr/d(lambda) = +/-sqrt(R(r)) with sign from initial radial momentum direction
   - d(theta)/d(lambda) = +/-sqrt(Theta(theta)) with sign from initial polar momentum direction

5. **Return state**: [t, r, theta, phi, dr/d(lambda), d(theta)/d(lambda), 0, 0]

## Geometry Trait Implementations

### `HasCoordinateSystem`
Returns `CoordinateSystem::Spherical`.

### `Signature`
Returns `[-1.0, 1.0, 1.0, 1.0]` (matches existing Kerr).

### `InnerProduct`
Full BL metric including g_t_phi cross term:
g_tt v^0 w^0 + g_rr v^1 w^1 + g_theta_theta v^2 w^2 + g_phi_phi v^3 w^3 + g_t_phi (v^0 w^3 + v^3 w^0)

### `inside_horizon`
r <= r_+ + epsilon, where r_+ = M + sqrt(M^2 - a^2). Trivial since position[1] is r directly.

### `get_radial_coordinate`
Returns `position[1]` directly (no quadratic solve needed).

### `get_tetrad_at`
Analytical BL tetrad for a ZAMO (zero angular momentum observer / locally non-rotating frame), orthonormalized via the BL metric.

### `lorentz_transformation`
Same structure as existing Kerr but using the non-diagonal BL metric.

### `SupportQuantities`
- `get_stationary_velocity_at`: ZAMO observer velocity in BL
- `get_circular_orbit_velocity_at`: Same angular velocity formula as existing Kerr, velocity expressed in BL coordinates directly
- `get_temperature_computer`: Reuse `KerrTemperatureComputer`

### `get_constants_of_motion`
Returns E, L_z computed from BL position and covariant momentum components.

## Testing

All tests are `#[test]` unit tests in `#[cfg(test)] mod tests` inside `kerr_bl.rs` (and optionally a cross-geometry test module).

### Test 1: Initial condition conversion
- Same camera position and ray for both `Kerr` and `KerrBL`
- Convert to BL, verify E, L_z, Q match what `Kerr::get_constants_of_motion` returns in Cartesian
- Also test against `Schwarzschild` for a=0 (BL conversion simplifies, initial state should match `SchwarzschildSolver::create_initial_state`)

### Test 2: Trajectory agreement against existing Kerr
- Identical physical parameters (r_s, a) and identical initial rays for both `Kerr` and `KerrBL`
- Integrate with both, convert positions to Cartesian, compare
- Tolerance: ~1e-4 (different parameterizations accumulate slightly different errors)
- Test with a=0 (Schwarzschild limit), a=0.3 (moderate), a=0.49 (near-extremal)

### Test 3: Constants of motion conservation
- Integrate a ray with `KerrBL`, record E, L_z, Q at initialization
- At each step, recompute from the current state
- Assert they remain constant within ~1e-8

### Test 4: Null condition preservation
- At each step, compute g_mu_nu p^mu p^nu using the BL metric
- Assert it remains zero within ~1e-8
- Validates the 1/Sigma momentum scaling

### Test 5: Horizon and escape agreement
- For rays that fall in, both `Kerr` and `KerrBL` must agree on `StopReason::HorizonReached`
- For rays that escape, both must agree on `StopReason::CelestialSphereReached`

### Test 6: Schwarzschild limit (a=0)
- With zero spin, `KerrBL` trajectories should match `Schwarzschild` geometry
- Tighter tolerance (~1e-6) since both use BL spherical coordinates

### Test 7: Redshift agreement
- For a ray hitting a disc, compare computed redshift between `Kerr` and `KerrBL`
- Validates `momentum_from_state`, `get_circular_orbit_velocity_at`, and the shading pipeline

### Test 8: Tetrad orthonormality
- Verify the ZAMO tetrad is orthonormal under the BL metric at several positions
- Same pattern as existing `test_tetrad_orthonormal` in `kerr.rs`

### Test 9: Metric inverse consistency
- Verify g^mu_nu computed analytically is the inverse of g_mu_nu
- Same pattern as `test_metric_contravariant_matches_inverse` in `kerr.rs`

## Performance Characteristics

Compared to the current `Kerr` (Kerr-Schild Cartesian, Hamiltonian geodesics):

- **No numerical metric derivatives**: The current Kerr computes 3 finite-differenced 4x4 matrix derivatives per step (each requiring 2 metric evaluations). KerrBL uses analytical R'(r) and Theta'(theta) — a few multiplications.
- **No Cartesian-to-BL conversion per step**: The current Kerr solves a quadratic (`compute_r_sqr`) at every evaluation. KerrBL works natively in BL.
- **Simpler RHS**: 6 meaningful components (2 evolved, 4 algebraic) vs 8 coupled components with matrix multiplications.
- **Constants computed once**: E, L_z, Q are fixed per ray, not implicitly re-derived at every step.

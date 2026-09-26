# Integrator roadmap: stability now, bundles and new spacetimes later

Scope: `src/rendering/runge_kutta.rs`, `src/rendering/integrator.rs`, the
state contract of `GeodesicSolver`, and the stop logic in the geometries.
Companion to `physics-review-2026-09.md` (findings 4.1, 4.8, §11, premortem
12) and `extensions-and-ideas-2026-09.md` (§8 bundles). The last section
answers whether the usual wormhole metrics can be traced in one chart.

## 0. Where the integrator stands

Facts, not judgements:

- RKF45 with an absolute error norm over all eight state components
  (`(CT·k).norm()`), safety 0.9, order-5 exponent, `H_MAX = 1`,
  `H_MIN = 1e-12`, growth cap 4, up to 100 retries per step.
- A step whose error is still above `epsilon` at `H_MIN` is **accepted**
  (`runge_kutta.rs`, the `h_cur <= H_MIN` branch).
- The loop stores every step, then calls `should_stop(last_y, y)` after
  the step; objects are tested on the segment between consecutive states.
- Constants-of-motion drift is tracked only in debug builds.
- Known consequences: prograde-edge error at high spin in the BL backend
  (finding 4.1), tunnel-and-explode on rays that approach a coordinate
  singularity without a guarding stop (§11), and the `H_MAX` comment in the
  code that documents why a larger cap tunnels through small objects.

Each item below is one change, with the reason, the test that proves it,
and what it unlocks.

## 1. Stop when the step collapses (do this first)

Change: when the error is still above tolerance at `H_MIN`, return a new
`StopReason::StepCollapsed` instead of accepting the step. Independently,
`N` consecutive accepted steps with `h_taken <= 10·H_MIN` (N ≈ 20) also
end the ray with the same reason. The stop is chart-independent.

Classification happens in the geometry, not the integrator: if the point
is within a few `horizon_epsilon` of `r₊` (from either side) the ray is
`HorizonReached` and painted as frozen; otherwise it is a
`CoordinateSingularity` and painted by the NaN policy (finding 4.8), never
silently black.

Why: this is the one change that makes an interior camera, over-extremal
Kerr, wormholes and any future metric safe. Today the horizon stop masks
the problem for the exterior camera only.

Test: turn the §11 scratch run into a unit test. A past-directed ray with
inward spatial momentum from `r = 0.65` (a = 0) must end as
`HorizonReached` in under 2 000 steps, with `r` never exceeding `1 + 1e-3`.
Today it ends at `r ≈ 1e24`.

## 2. Error norm with per-component scaling

Change: replace the absolute norm by Hairer's scaled norm

```
sc_i  = atol_i + rtol · max(|y_i|, |y_new_i|)
err   = sqrt( (1/n) · Σ_i (e_i / sc_i)² )      accept if err ≤ 1
```

with `atol` per component class: positions in units of `r_s`, momenta in
units of the ray's energy (seed rays with `p_t = −1` so momenta are O(1)),
deviation vectors (§7) normalised to unit initial size.

Caveat: a plain relative norm was tried during the review as a fix for
finding 4.1 and made things worse. The BL problem is conditioning of the
equations, not the norm. So introduce `--rtol`/`--atol` as new options
whose defaults reproduce today's images, keep `--epsilon` meaning what it
means, and validate with the shadow-edge bisection script.

Optional: a PI step controller (Gustafsson),
`h_new = h · 0.9 · err^(−0.7/5) · err_prev^(0.4/5)`, gives smoother step
sequences and fewer rejections. Low priority.

## 3. Geometry-aware step ceiling

Today `H_MAX = 1` is global, and the code comment records that a larger
cap tunnels through a radius-2 sphere at r = 10. Replace the global cap by
two local bounds applied to the controller's proposal:

```
h ≤ safety · d(x) / |dx/dλ|       d = signed distance to the nearest object (SDF from the object list / octree), safety ≈ 0.7
h ≤ c · L(x)                       L = local geometric scale: r for black holes, r(l) for wormholes; c ≈ 0.5
```

Then `H_MAX` can be large in empty regions (fewer steps at r = 1000) and
thin discs or a wormhole throat cannot be stepped over. The hit test
becomes a sign change of the SDF between steps, refined by §4.

Test: the existing `test_color_of_ray_hits_sphere*` with `H_MAX = 100`
must still hit.

## 4. Dense output and event location

RKF45 already has `y`, `f(y)`, `y_new`, `f(y_new)` per step, which is a
cubic Hermite interpolant for free. Use it to locate, inside a step:

- object crossings (disc, sphere): root of the SDF on the interpolant
  instead of linear interpolation between states;
- the celestial-sphere crossing at the exact `max_radius`;
- equatorial-plane crossings, counted for the image-order channel;
- the horizon band for §1's classification.

Why: disc hit positions and redshifts stop depending on the step size
where the ray crosses, and the volumetric disc (finding 4.6) gets exact
entry/exit points.

Test: straight ray in `EuclideanSpace` hitting a plane at a known point
with `H_MAX = 100`; error below 1e-8.

## 5. First-integral monitoring, optional projection

Promote the debug-only drift accounting to a cheap always-on per-ray
record: max `|H|/|p|²`, `ΔE/E`, `ΔL` (and `ΔQ` for BL). Expose it in the
`render-ray-at` CSV and assert bounds in tests. It is the diagnostic that
found §11 and 4.1.

Optional, only if drift becomes visible: project onto `H = 0` every k
steps by re-solving `p_t` from the null condition. Cheap and safe for the
Hamiltonian backends. Never "project" the BL Mino velocities `v_r`, `v_θ`;
their sign is exactly what finding 4.1 is about.

## 6. One first-order Hamiltonian state for all backends

- BL: adopt the super-Hamiltonian first-order form (DNGR A.15, finding
  4.1) with `b`, `q` as parameters.
- Schwarzschild spherical: keep or switch; either is fine, but the state
  must be `(x, p)` plus optional constants.
- Fold constants of motion into the state with zero derivative when a
  backend uses them, so §7 and §8 are one code path.

The `Point` type stays four-component. This is a solver-side refactor
only; the `Geometry` trait boundary is right where it is.

## 7. Deviation vectors as extra state (ray bundles)

State grows from 8 to `8 + 8·n_dev` with `n_dev = 2`. `rkf45` is already
generic over the dimension. The right-hand side of the deviation part is
the directional derivative of the geodesic right-hand side, obtained by
evaluating the same `geodesic` function on a dual-valued state (§8).

Rules:

- include the deviation components in the scaled norm (§2), normalised;
- near the photon ring they grow exponentially (Lyapunov exponent of the
  review's photon-ring section). That is physical. Exclude them from the
  stall criterion of §1, and renormalise a deviation vector when its
  norm exceeds ~1e6 while keeping a log-scale factor;
- readout at the escape radius: the part of `δp` orthogonal to `p`,
  divided by `|p|`, is the sky angle offset; two of them form the 2×2
  Jacobian; singular values × beam radius are the ellipse, `1/|det|` the
  magnification, the sign the parity;
- diagnostics: `ω(δy₁, δy₂) = δx₁·δp₂ − δx₂·δp₁` constant along the ray,
  and `δH = 0`.

Tests: radial ray in Schwarzschild stays circular; Perlick §4.3 closed
forms for `D₊`, `D₋`; far-field magnification `(u²+2)/(u√(u²+4))`.

## 8. Automatic differentiation of the metric

Make `metric` (and the Christoffels where used) generic over the scalar
type. First-order dual numbers give `∂g` exactly and replace the
36-evaluation stencil in `KerrSolver::d_matrix_covariant`; hyper-dual (or
nested dual) numbers give `∂F/∂y · δy` for §7. `num-dual` provides both.
This is also the change the GPU discussion needs (generic scalar), and it
removes the finite-difference noise that sits at the integrator's
tolerance today.

## 9. Budgets and classification

- Add an affine-parameter (or arc-length) budget next to `max_steps`; the
  two fail differently and should be reported differently.
- Classify trapped rays by winding count, not by "ran out of steps inside
  5 r_s".
- One NaN policy (finding 4.8), applied by the compositor, never a
  hard-coded black.

## 10. Order of work

1. §1 stall stop (small, unlocks everything else).
2. §5 diagnostics (small, needed to judge the rest).
3. §2 scaled norm behind new options.
4. §4 Hermite event location.
5. §3 SDF step ceiling (needs §4).
6. §6 first-order BL state.
7. §8 automatic differentiation.
8. §7 deviation vectors and the star filter that goes with them.

## 11. Tests that should exist afterwards

- `stall_stop_past_horizon` (§1), interior and exterior camera.
- `invariants_harness`: max drifts per ray below tolerance-scaled bounds
  for each backend, at ε = 1e-5 and 1e-9.
- `shadow_edge_regression`: the bisection script's numbers as a test.
- `hermite_event_location` (§4), `sdf_step_ceiling_hits_thin_disc` (§3).
- `bundle_radial_ray_stays_circular`, `bundle_matches_perlick_4_3`,
  `bundle_symplectic_form_constant` (§7).

## 12. Wormholes in one chart

Yes, one chart covers the whole spacetime for the usual wormholes. Use
`(t, l, θ, φ)` with the proper radial coordinate `l ∈ (−∞, ∞)`:

```
ds² = −e^{2Φ(l)} dt² + dl² + r(l)² (dθ² + sin²θ dφ²)
```

- Ellis: `r(l) = √(l² + b₀²)`, `Φ = 0`. Ultrastatic, one parameter.
- DNGR wormhole (James et al., arXiv:1502.03809): `r(l) = ρ` for
  `|l| ≤ a`, and `r(l) = ρ + M (x·arctan x − ½ ln(1 + x²))` with
  `x = 2(|l| − a)/(πM)` for `|l| > a`; `Φ = 0`. Three parameters: throat
  radius, length, lensing width.
- Morris–Thorne with a shape function `b(r)`: the `r` chart is
  double-valued and has `g_rr → ∞` at the throat, a coordinate
  singularity. Convert once to `l` by quadrature, `dl = dr/√(1 − b/r)`,
  tabulate `r(l)`, and trace in the `l` chart.

Why one chart is enough: the metric functions are smooth and
non-degenerate for every `l`, since `r(l) ≥ b₀ > 0` and `e^{2Φ}` is
finite. There is no horizon. The only coordinate singularity is the polar
axis of the spherical chart, exactly as in the Schwarzschild spherical
backend; either accept it, or tilt the chart's axis away from the
camera's view, or use the Hamiltonian form where the axis term vanishes
identically for `p_φ = 0`. A two-chart alternative (a Cartesian chart per
side, glued at the throat) is only worth it if objects must be described
in Cartesian coordinates on both sides. Not recommended.

What the code needs:

- a `CoordinateSystem` variant with a signed radial coordinate, and
  `to_cartesian` that carries the side (`sign(l)`) so objects live on one
  side;
- a Hamiltonian solver in this chart (closed-form `∂g`, or §8);
- `inside_horizon` always false; stop on `|l| > l_max` with the side
  attached to `CelestialSphereReached`, and two sky textures;
- the static tetrad `e_t = e^{−Φ} ∂_t` exists everywhere, so the camera
  can sit anywhere, including at the throat;
- redshift `g = e^{Φ(obs) − Φ(src)}`, which is 1 for Ellis and DNGR;
- the throat is a photon sphere (Ellis: unstable circular orbit at
  `l = 0`), so rays can circle there; the step budget or winding count of
  §9 classifies them, and the Einstein ring and photon rings produce the
  same folds as a black hole, so §7 applies unchanged;
- step ceiling `h ≤ c · r(l)` (§3) so nothing steps across the throat
  region in one go.

Validation: the exact lens map of Perlick §4.3 with `R(r) = r(l)`, which
gives the total bending `Φ(Θ)` for Ellis in closed form, and the Einstein
ring radius; for the DNGR wormhole, the figures in arXiv:1502.03809.

Concepts to look up: proper radial distance as a global chart, the
Morris–Thorne shape and redshift functions, embedding diagrams of throats,
photon spheres of wormholes, and the exact lens map for spherically
symmetric static spacetimes.

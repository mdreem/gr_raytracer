# Physics and correctness review (September 2026)

Scope: every file under `src/geometry` and `src/rendering`, the scene objects,
the CLI wrappers, the scene definitions and the gallery recipes, read in full.
The physics was checked by hand against the standard references and against
the two papers the project builds on: Riazuelo, *Seeing relativity I*
(arXiv:1511.06025) and James, von Tunzelmann, Franklin & Thorne,
*Gravitational Lensing by Spinning Black Holes in Astrophysics, and in the
Movie Interstellar* (CQG 32, 065001; arXiv:1502.03808, the DNGR paper).
On top of the reading, four numerical experiments were run against closed-form
results; their scripts are in `scripts/validation/`. The full test suite passes
(187 tests).

Units in this document follow the code: `r_s = 2M` is the length unit, spin
`a` is in the same unit, so `a = 0.499` with `r_s = 1` is `a/M = 0.998`.

## 1. Verdict in short

The physics that is implemented is, with one important exception, implemented
correctly: the metrics, tetrads, Lorentz boost, camera projection, geodesic
equations, conserved quantities, redshift/Doppler factors, blackbody
transformation, Novikov–Thorne flux shape and the ray-tube star gather all
check out analytically and numerically. The architecture (geometry as a chart
adapter, conserved Killing quantities instead of parallel transport, tubes for
point stars) is the same design DNGR and Riazuelo use and is sound.

The exception is a real bug with visible consequences: the Boyer–Lindquist
backend (`KerrBL`) integrates the radial equation in a form that loses the
radial first integral, and this shifts the *prograde* edge of the shadow and
the prograde photon ring at high spin unless the integrator tolerance is much
tighter than the CLI default. Measured: with the default `--epsilon=1e-5` and a
camera at `r = 18`, the prograde edge at `a/M = 0.998` is 22 % too close to the
hole; with the gallery's `1e-9` it is correct to 2·10⁻⁴. Section 3 has the
evidence, Section 4.1 the root cause and the fix.

Three medium findings follow: the star layer is not radiometrically consistent
with the disc (it multiplies by a coordinate-basis pixel solid angle), the two
Kerr backends disagree on what a "Cartesian" position means by an r-dependent
azimuth twist, and the Novikov–Thorne temperature calibration overshoots by 18 %
at high spin. Everything else is documentation, robustness or scope.

## 2. What was verified as correct

Analytic checks (by hand, against the references named):

- Schwarzschild geodesic equations (`schwarzschild.rs`): all four Christoffel
  combinations correct for signature (+,−,−,−).
- Kerr–Schild metric `g = η + f k k` with `f = r_s r³/(r⁴ + a² z²)` and the
  ingoing null covector `k`; inverse `g⁻¹ = η − f k k` (valid because `k` is
  null). Eulerian-observer tetrad `(1/α, −β/α)` with `α² = 1/(1+f)`,
  `β_i = f k_i/(1+f)`: correct lapse and shift of the KS 3+1 split.
- Boyer–Lindquist metric, ZAMO (`ω = −g_tφ/g_φφ`), circular-orbit `Ω`,
  `u^t`, conserved `E`, `L`, ISCO (Bardeen–Press–Teukolsky, both branches).
- Mino-time potentials `R(r)`, `Θ(θ)` and their derivatives (null forms),
  Carter constant `Q = p_θ² + cos²θ (L²/sin²θ − a²E²)`.
- The BL→KS Jacobian including the `dt = r_s r/Δ dr`, `dφ = a/Δ dr` twist
  (also pinned by an existing exact test).
- Lorentz boost: identical to Riazuelo eq. (6), with the sign flip for the
  (−,+,+,+) geometries handled correctly.
- Camera: `w = e_z + i' e_x + j' e_y`, `n = −e_z + 2w/(w·w)` is exactly
  Riazuelo eqs. (4)–(5). See finding 4.5 for what that implies.
- Redshift `g = (u_obs·p)/(u_em·p)` with the emitter velocity assembled from
  Killing coefficients: chart-invariant, sign-invariant, and the existing test
  reproduces Luminet (1979) for the circular-orbit emitter.
- Blackbody: `I_λ^obs(λ) = g⁵ B_λ(gλ, T) ≡ B_λ(λ, gT)`; star flux `g⁴` times
  lensing magnification; both are the standard bolometric/spectral results.
- Novikov–Thorne flux shape (Page & Thorne 1974) with `√(−g) = r` for the
  (t, r, φ) equatorial metric; the constant prefactor is calibrated away.
- Star-tube approach (four corner rays, subdivision on winding spread,
  magnification `Ω_screen/Ω_sky`): the finite-difference version of what
  Riazuelo does with two nearby directions (his eqs. 45–50) and DNGR does with
  geodesic deviation.
- RKF45 tableau (the Wikipedia "Formula 2" coefficients), local extrapolation
  with the 5th-order solution, standard step controller.

## 3. Numerical validation

### 3.1 Shadow edge versus Bardeen's critical curve

For an equatorial observer, the two points where the shadow boundary crosses
the equatorial line of the image have impact parameters `b = ξ(r_ph±)` with
`ξ(r) = [M(r² − a²) − rΔ]/[a(r − M)]` evaluated at the prograde/retrograde
circular photon orbits (James et al. eqs. A.5–A.6, with `q = 0`). The script
`scripts/validation/shadow_edge_bisect.py` fires equatorial rays from a distant
observer with `render-ray-at`, bisects on the local angle between captured and
escaped rays, converts the angle to `b = L/E` exactly for that observer's
frame, and compares. This exercises the integrator, the tetrad, and the
conserved-quantity seeding end to end.

Camera at `r = 1000` (aberration negligible), 18 bisection steps:

| a (r_s) | a/M | backend | ε | side | b numeric | b analytic | rel. error |
|---|---|---|---|---|---|---|---|
| 0 | 0 | KerrBL | 1e-9 | both | ±2.59775 | ±2.59808 | −1.2e-4 |
| 0 | 0 | Kerr (KS) | 1e-9 | both | ±2.60140 | ±2.59808 | +1.3e-3 (Eulerian-frame aberration, expected) |
| 0.25 | 0.5 | KerrBL | 1e-9 | prograde | +2.04523 | +2.04813 | −1.4e-3 |
| 0.25 | 0.5 | KerrBL | 1e-9 | retrograde | −3.06857 | −3.06908 | −1.7e-4 |
| 0.499 | 0.998 | KerrBL | 1e-9 | prograde | **+0.95437** | +1.05544 | **−9.6e-2** |
| 0.499 | 0.998 | KerrBL | 1e-9 | retrograde | −3.49805 | −3.49833 | −8.1e-5 |
| 0.499 | 0.998 | KerrBL | 1e-7 | prograde | +0.27141 | +1.05544 | −7.4e-1 |
| 0.499 | 0.998 | KerrBL | 1e-5 | prograde | ≈0 (every ray escapes) | +1.05544 | −1.0 |
| 0.499 | 0.998 | Kerr (KS) | 1e-5 | prograde | +1.05815 | +1.05544 | +2.6e-3 |
| 0.499 | 0.998 | Kerr (KS) | 1e-5 | retrograde | −3.50272 | −3.49833 | +1.3e-3 |

Camera at `r = 18` (the gallery's camera distance), `a/M = 0.998`, KerrBL:

| ε | prograde rel. error | retrograde rel. error |
|---|---|---|
| 1e-5 (CLI default) | **−2.2e-1** | +6e-5 |
| 1e-7 | −2.4e-2 | +3e-4 |
| 1e-9 (gallery recipes) | −1.7e-4 | +4e-4 |
| 1e-11 | +3.1e-4 | +4e-4 |

Camera at `r = 18`, ε = 1e-5, KerrBL, lower spins: `a/M = 0.5` prograde
−4.8e-3, retrograde −4.9e-4; `a = 0`: −1.3e-3 both sides. The Schwarzschild
backend (spherical chart, affine parameter) at `r = 18`: +2.1e-5 at both
ε = 1e-5 and 1e-9, i.e. limited by the bisection itself.

Reading: the Kerr–Schild backend is right at all spins with default settings.
The Boyer–Lindquist backend is right for retrograde light and at low spin, and
wrong on the prograde side at high spin unless ε is very small. The error is
a strong function of ε, which points at the integrator formulation rather than
the equations (Section 4.1).

### 3.2 First-integral diagnostic (KerrBL, `a/M = 0.998`, `b = 1.00`)

A ray with `b = 1.00 < b_c = 1.055` must be captured (`R(r) > 0` all the way
to the horizon: `R = 1.1e-3` at `r_+`). Instrumenting the trajectory: at
ε = 1e-9 the ray turns around at `r = 0.5556` with `(dr/dλ')² = 3e-8` where
`R(r) = 2.9e-3`, i.e. the first integral `(dr/dλ')² = R(r)` is violated by
100 %, and `|k·k|` reaches 4.5 near the hole. At ε = 1e-12 the turning point
moves to `r = 0.548`; a `b = 1.03` ray is captured at 1e-12 but escapes at
1e-9.

### 3.3 Cross-chart azimuth twist

Re-running the setup of `test_trajectory_agreement_with_kerr` and comparing
final *directions* instead of positions: the KS and BL escape directions differ
by 1.71°. The predicted twist between the KS azimuth and the BL azimuth used
by `Point::to_cartesian`, `[F(∞) − F(r_cam)]·sin θ` with `F' = a/Δ`, is
1.77° × 0.965 = 1.71°. The discrepancy the test tolerates is exactly this
chart difference (Section 4.3), not different time coordinates.

### 3.4 Novikov–Thorne calibration

`KerrTemperatureComputer::new(T, outer, a, r_s)` calibrates `Ṁ` so that the
maximum over ten coarse samples equals `T`. Scanning the resulting profile
finely:

| a (r_s) | requested peak T | actual peak T | peak radius / r_ISCO |
|---|---|---|---|
| 0.499 | 8000 K | 9419 K (+18 %) | 1.28 |
| 0.25 | 8000 K | 8270 K (+3 %) | 1.56 |
| 0 | 8000 K | 8000 K | 1.59 |

At high spin the true peak sits inside the first sample interval
(`r_ISCO + 0.5·dr` with `dr = (outer − r_ISCO)/10`), so the "temperature"
parameter is not the peak temperature the gallery captions quote.

## 4. Findings

Severity: **High** = wrong pixels in supported configurations; **Medium** =
systematic error or inconsistency that changes what is rendered;
**Low** = documentation, robustness, or corner cases.

### 4.1 High: KerrBL loses the radial first integral (prograde shadow edge and photon ring wrong at high spin unless ε ≤ 1e-9)

Where: `kerr_bl.rs` (`KerrBLSolver::geodesic`, `create_initial_state`),
`runge_kutta.rs` (absolute error control), `cli.rs` (default `epsilon = 1e-5`).

Mechanism. The solver integrates `d²r/dλ'² = R'(r)/2` and
`d²θ/dλ'² = Θ'(θ)/2` in Mino time, with the state carrying `v_r = dr/dλ'`.
That form is well behaved at turning points (no `±√R` sign bookkeeping), but
the physically decisive quantity near the photon shell, `C = v_r² − R(r)`,
is not a state variable and is not preserved by the integrator. The scale
mismatch is enormous: `v_r ≈ r² E` is ~3·10² at a camera at `r = 18` and
~10⁶ at `r = 1000`, so `v_r²` is 10⁵–10¹², while `R(r)` at the prograde
photon orbit of a near-extremal hole is ~10⁻³ (`P = (r²+a²)E − aL` and `Δ`
are both nearly zero there). RKF45's *absolute* tolerance ε bounds the per-step
error in `v_r`, so each step can inject `δC ≈ 2 v_r ε`, and near the prograde
orbit the sensitivity is `δb ≈ δC / (∂R/∂L) ≈ 100·δC` at `a/M = 0.998`
(versus ~0.05·δC on the retrograde side, whose photon orbit sits at `r ≈ 2`
with `R` of order one). That is why only prograde, only high spin, and why the
error tracks ε by orders of magnitude (Section 3.1). Tightening ε is a
workaround, not a fix: at `r_cam = 1000` even ε = 1e-9 leaves 10 %, and an
absolute tolerance on a state that spans twelve decades has no stable meaning.
A mixed absolute/relative error norm was tried and makes it worse (a relative
tolerance on `v_r ~ 10⁶` permits absolute errors of 10⁻³).

Who is affected. Every KerrBL render with the CLI default ε: the
`docs/example-render-commands.md` KerrBL command, `scene-definitions/kerr-bl*.toml`
and `kerr-animation*.toml` run from `scripts/rendering/create_kerr_images.sh`
(no `--epsilon`, camera at `r ≈ 22`), and any user following the README. The
gallery images (`--epsilon=1e-9`, camera `r ≈ 18–22`) are fine at the 2·10⁻⁴
level. The Kerr–Schild backend is unaffected because it integrates the
covariant momentum `p_μ`, which stays O(1) everywhere outside the horizon.

Fix (recommended, in order of preference):

1. Switch the BL solver to the super-Hamiltonian first-order form in the
   affine parameter, with covariant `p_r`, `p_θ` as state (James et al.
   eq. A.15; MTW §33.5): `dr/dζ = Δ p_r/ρ²`, `dθ/dζ = p_θ/ρ²`,
   `dp_r/dζ = ∂_r[−Δp_r²/(2ρ²) − p_θ²/(2ρ²) + (R + ΔΘ)/(2Δρ²)]`, and
   similarly for `p_θ`. `p_r = ±√R/Δ` is O(1) from the camera to the celestial
   sphere, turning points need no sign logic, and the Hamiltonian constraint
   `H = 0` is an O(1) quantity the tolerance actually controls. DNGR uses
   exactly this form and its footnote warns against the `±√R` family. The
   Mino-time separability is still used for `dt/dζ` and `dφ/dζ`.
2. If Mino time is kept: integrate `p_r = v_r/Δ` instead of `v_r` (bounded
   outside the horizon), and add a constraint monitor: abort or refine when
   `|v_r² − R(r)| > κ·max(|R|, v_r²)` near the hole, and log it as a stop
   reason instead of letting the ray decide capture/escape on noise.
3. Regardless: make the shadow-edge bisection a regression test (a `#[ignore]`
   release test is fine), because it is the only check in the suite sensitive
   to this. The existing tests assert `E`, `L_z`, `Q` drift < 1e-4, but `E`
   and `L_z` are inputs to the BL solver and are conserved by construction, and
   the `k·k < 1e-4` test only looks at one benign ray.

Documentation until fixed: state that KerrBL needs `--epsilon ≤ 1e-9` above
`a/M ≈ 0.9`, and set the default ε per geometry.

### 4.2 Medium: the star layer is not radiometrically consistent with the disc, and uses coordinate-basis solid angles

Where: `raytracer.rs` (`compute_color`, `flat_gather`, `gather_curved`,
`magnification`).

Two separate problems share one line: the tube's star light is
`(Ω_screen / Ω_sky) · Σ g⁴ F_*`, with `Ω_screen` computed from the four corner
momenta's coordinate components (`momentum.get_cartesian_vector`).

(a) Radiometry. The disc and textures are rendered as radiance (light per
solid angle); the pixel buffer is a radiance map. A star of flux `F_*` lensed
with magnification `μ = Ω_screen/Ω_sky` contributes flux `g⁴ μ F_*` to the
pixel, so its *radiance* contribution is `g⁴ μ F_*/Ω_screen = g⁴ F_*/Ω_sky`.
The screen solid angle cancels. Multiplying by `Ω_screen` turns the star layer
into flux-per-pixel while the disc stays radiance. With the stereographic camera
the pixel solid angle falls as `4/(1+ρ²)²` from the image centre (ρ = 0) to
the edge of a 90° opening (ρ = 1), so relative to the disc, stars at the edge
of every gallery frame are 4× too faint (9× in the corners). A uniform star
background would render with strong vignetting. `flux_scale` cannot absorb a
position-dependent factor.

(b) Frame. Even if flux-per-pixel were intended, the pixel solid angle must be
the physical one in the camera's orthonormal frame (Riazuelo computes his
amplification "in the observer's frame, expanded in `X, Y, Z`", eqs. 45–50).
Coordinate components of `p^μ` are not physical angles: in Schwarzschild the
radial component is stretched by `1/√(1 − r_s/r)`, in BL by `√(Σ/Δ)`, and the
camera's own boost (aberration) is ignored entirely. At `r = 18` this is a few
percent and direction dependent; for the ZAMO close-vantage renders at
`r ≈ 1.1 r_s` it is factors of several across the field.

Fix: drop `Ω_screen` and return `g⁴ Σ F_* / Ω_sky` as the tube's radiance
(the sky footprint is computed at `r = 15000` where space is flat, so it is
already correct). If a pixel-solid-angle weighting is ever wanted, take it from
the frame components of the ray direction, which `get_direction_for` already
has (`n = −e_z + 2w/|w|²` in the tetrad basis); store them in `Ray`.
Existing tube tests set `Ω_screen = Ω_sky` (flat space, small FOV) and would
not notice either change; add a test with a uniform synthetic catalogue and a
wide FOV asserting constant radiance across the frame.

### 4.3 Medium: the two Kerr backends interpret "Cartesian" through charts that differ by an r-dependent rotation

Where: `point.rs` (`to_cartesian` for `BoyerLindquist`),
`spherical_coordinates_helper.rs` (`cartesian_to_boyer_lindquist`), and the
comment in `kerr_bl.rs::test_trajectory_agreement_with_kerr`.

The Kerr–Schild Cartesian coordinates satisfy `x + iy = (r + ia) e^{iφ_KS} sin θ`
with `φ_KS = φ_BL + F(r)`, `F'(r) = a/Δ`. KerrBL's embedding uses `φ_BL` in
place of `φ_KS`. Because the metric is axisymmetric this never affects a
geodesic, and the BL momentum conversion is right (the Jacobian is evaluated at
the embedding angle, which is `φ_KS`), but every Cartesian *position* KerrBL
reports is rotated by `F(r_cam) − F(r)` relative to the KS backend. Measured:
1.71° at the celestial sphere for `a = 0.3`, `r_cam ≈ 10`, matching the
prediction (Section 3.3). Consequences: the Gaia sky is rotated by
`≈ a/r_cam` (1.6° at `a = 0.499`, `r = 18`) between the two backends; a
checker disc's phase differs by a radius-dependent spiral; a sphere placed at
"(0, 0, 20)" is at a slightly different azimuth in the two backends; the
cross-chart test needs a 1393-unit tolerance and its comment blames "different
time coordinates" and a 7 % difference in `E`, which is wrong (`E = −p·∂_t` is
chart-invariant, and the file's own `test_e_lz_consistency_between_ks_and_bl`
asserts it to 1e-10).

Fix: define the BL↔Cartesian map to *be* the KS chart. `F` has the closed form
`F(r) = a/(r_+ − r_-) · ln|(r − r_+)/(r − r_-)|` for `a < M` (with the
constant chosen so `F(∞) = 0`), so `to_cartesian` uses `φ_BL + F(r)` and
`cartesian_to_boyer_lindquist` subtracts it. Then the two backends agree on
positions to integration accuracy and the trajectory test can be tightened by
two orders of magnitude.

### 4.4 Medium/Low: Novikov–Thorne calibration overshoots at high spin

Where: `temperature.rs` (`NUM_STEPS_FIND_MAXIMUM = 10`). Evidence in
Section 3.4: +18 % peak temperature at `a/M = 0.998`, +3 % at `a/M = 0.5`.
Fix: locate the maximum by a fine scan or golden-section search on the LUT
(the LUT is built anyway; calibrate after building it), or search on a
logarithmic radial grid starting at `r_ISCO`. Also replace the finite-difference
`dΩ/dr` (`h = 1e-10`) with the analytic derivative of
`Ω = √M/(r^{3/2} + a√M)`.

### 4.5 Low: camera is stereographic with α = half opening angle; the doc comment says pinhole and the FOV is hard-coded

`Camera::get_direction_for` implements Riazuelo eqs. (4)–(5) exactly: a
stereographic (conformal) projection, chosen by him so that the shadow stays a
circle on screen. In that formula `α` is the *half* opening angle, so the
hard-coded `PI / 4.0` in `cli/shared.rs::create_scene` is a 90° vertical
opening angle, and with `alpha = π/2` (used in several tests) the corner pixels
look 104° off-axis. The doc comment on `get_direction_for` ("pinhole camera
projection", "k·k = 0") is misleading, the README never states the field of
view, and it cannot be changed from the CLI. Not a bug; document it, expose
`--fov`, and consider a gnomonic option for flat-sky comparisons with other
codes.

### 4.6 Low: volumetric disc temperature uses the Cartesian cylindrical radius

`volumetric_disc.rs::march_constant_step` passes `ρ = |p × axis|` to the
temperature LUT, while `Disc` and the ISCO logic use the BL radial coordinate
(`ρ² = r² + a²` on the equator). At `a = 0.499` the inner edge is evaluated at
`r` off by up to 0.2 `r_s`. Use `geometry.get_radial_coordinate` as the flat
disc does.

### 4.7 Low: angle wrapping and pole crossing

`Point::get_as_spherical` wraps `θ` with `rem_euclid(π)`, which maps
`θ ∈ (π, 2π)` to `θ − π` instead of `2π − θ` with `φ + π`. Only rays that
cross a coordinate pole are affected (they are already numerically fragile
because of `1/sin θ`), so this is latent, but it silently misplaces such rays
on the celestial texture. `wrap_phi`'s comment says `(−π/2, π/2]` while the
code returns `(−π, π]`.

### 4.8 Low: silent failure modes paint black

`scene.rs::color_of_ray` maps `CoordinateIsNan`, step-budget exhaustion and
`ClosedOrbitDetected` all to `Captured`/black. Near the BL polar axis or with
too few steps this produces dark pixels indistinguishable from the shadow, and
the adaptive sampler then treats them as real edges. Consider a distinct stop
reason count in the render summary and, for NaN, a retry with a tighter ε or
a visibly wrong debug colour.

### 4.9 Low: sphere emitters use the static observer

`sphere.rs::energy_of_emitter` uses `get_stationary_velocity_at`, which does
not exist inside the ergosphere (NaN → `UnphysicalRedshift` → the whole pixel
errors). Use the ZAMO, as the camera already can.

### 4.10 Informational

- `kerr.rs` evaluates the metric 36 times per right-hand side for numerical
  derivatives. The KS structure `g = η + f k k` has closed-form derivatives
  (`∂f`, `∂k` via `∂r/∂x = ...` from the `r⁴ − (ρ² − a²) r² − a² z² = 0`
  relation); analytic Christoffels would make the KS backend roughly as fast
  as BL while keeping its accuracy, which given 4.1 is the more valuable
  backend.
- `apply_beaming` with `beaming_exponent = 3` on bitmap discs and the celestial
  texture is an artistic `g³` in addition to nothing physical (bitmaps carry no
  spectrum). Fine, but the scene files and README should say it is artistic;
  the blackbody path is the physical one.
- `closed_orbit`'s "ran out of steps inside 5 r_s" heuristic is reasonable but
  its threshold is not documented in the README's `--max-steps` guidance.
- `H_MAX = 1.0` is a step cap in whatever parameter the solver uses; in Mino
  time it is effectively no cap far from the hole (a Mino step of 1 at
  `r = 100` moves `r` by 10⁴).

## 5. Architecture assessment

The user asked whether the architectural style is right. Mostly yes, and the
places where it is not are the same places the bugs live.

What is good and matches the literature:

- `Geometry` as a chart adapter (metric, tetrad, horizon, radial coordinate,
  Killing vectors) with the integrator, camera and objects written against it.
  DNGR is structured the same way (FIDO frame + boost + conserved quantities).
- The Killing decomposition (`OrbitKillingDecomposition`, `RayFrequencyData`):
  computing `u·p = u^t p_t + u^φ p_φ` from conserved scalars instead of
  parallel-transporting anything is exactly right for stationary axisymmetric
  spacetimes and avoids a whole class of chart bugs. Keep it; extend it to the
  sphere emitter.
- Premultiplied radiance with explicit transmittance (`Radiance::over`) and the
  separation of straight texture colour from integrated light is clean and the
  compositing tests are strong.
- Ray tubes for point stars with winding-based subdivision and the curved
  boundary option: the right idea (DNGR §2.2, Riazuelo §V.B).

What should change:

1. The integrator state is untyped. `EquationOfMotionState` means
   `(x, p_μ)` for KS, `(x, ẋ)` for Schwarzschild, `(x, dr/dλ', dθ/dλ', 0, 0)`
   for BL, and the affine parameter means three different things. That leak
   already forced `disc.rs` to special-case BL with a `param_scale`
   closure, and it is why the tolerance problem in 4.1 went unnoticed: the
   integrator cannot know which components are physical. Give each solver its
   own state type or, at minimum, a `parametrization()` and a
   `constraint_residual(&state) -> f64` method the integrator can monitor.
2. `Ray` should carry the frame components of the pixel direction (the
   `(i', j')` stereographic point or the unit vector in the tetrad basis).
   Every place that currently reaches for `momentum.get_cartesian_vector` to
   get "the direction on screen" is approximating that.
3. "Cartesian" must mean one chart. Make the KS chart the definition and give
   `KerrBL` the exact map (4.3). Then `Point::to_coordinate_system` is a true
   change of chart, not an approximate one.
4. Radiometry should be one currency (radiance) from the camera to the file;
   flux-per-pixel appears only at the sensor model, if ever (4.2).
5. The global `OnceLock` for `--pole-rotation-deg` is process state that a
   library user cannot set twice; it belongs in `Schwarzschild` (or in a
   `Chart` object). Also note the option silently does not rotate the flat
   disc, which the CLI help mentions but the config does not enforce.
6. Physics regression tests should be against closed forms, not against
   previously observed numbers with tolerances chosen to pass. The suite has
   good examples (Luminet, Shakura–Sunyaev slope, ISCO limits, BL Jacobian) and
   bad ones (the 1393-unit cross-chart tolerance, the `> 200 matching points`
   trajectory comparisons). Add: Bardeen shadow edges for both backends and
   several spins, the Schwarzschild deflection angle `4M/b` at large `b`, the
   Einstein-ring radius versus `θ_E`, and a uniform-sky star-layer flatness
   test.

## 6. Comparison with the papers

Riazuelo (1511.06025):

- Camera and tetrad: identical (free-fall tetrad from rest at infinity, his
  eq. 8; boost, his eq. 6; stereographic projection, eqs. 4–5; wave vector
  `u − N`, which the code writes as `N − u`, the same null line).
- Stars: he uses spherical symmetry to precompute the deflection function
  `χ(δ)` once per observer (his §IV), then places each star's images by
  inverting it, including all ghost images (`φ_∞ + 2kπ`), and computes the
  amplification from two nearby directions in the observer frame (§V.B). For
  the Schwarzschild backend this is a 100–1000× speed-up over per-pixel tracing
  and gives exact, flicker-free star positions; it is a natural extension here
  (Section 8).
- Redshift: his eq. (22), `1 + z = E/(k·u_obs)`, is what `RedshiftComputer`
  does.

James et al. (DNGR, 1502.03808):

- Frame: FIDO (= ZAMO) orthonormal basis plus the camera's boost `β` relative
  to it (their A.7–A.9). The code's `CameraVelocityConfig::Zamo` plus
  `lorentz_transform_tetrad` is the same construction; their circular-orbit
  camera `β = ϖ(Ω − ω)/α` is a one-line addition here.
- Geodesics: they integrate the super-Hamiltonian first-order system (A.15)
  and explicitly reject `±√R` forms. The code's second-order Mino form avoids
  the sign problem but not the conditioning problem (4.1); adopting A.15 is
  the fix.
- Ray bundles: DNGR propagates the bundle's cross-section with the geodesic
  deviation equation (Pineault–Roeder form, their A.2), obtaining the ellipse
  `(δ+, δ−, μ)` on the celestial sphere per pixel; the star's intensity comes
  from the beam's solid-angle change (A.3.1), and they smooth with a truncated
  Gaussian over ~2 pixel radii to hold star flicker to 2 %. The code's
  four-corner tube is the finite-difference cousin; it is correct away from
  folds but needs the subdivision/`magnification_cap` machinery exactly where
  the Jacobi-field approach is exact and cheap (Section 8).
- Disc: DNGR's disc for the film deliberately omitted Doppler and gravitational
  colour shifts and intensity changes, and the disc scenes used `a/M = 0.6`;
  the paper's Fig. 15-style "physically correct" images show the approaching
  side blue and bright. This code does the physical version by default, which
  is the right choice for a physics renderer.
- Shadow: the near-extremal flattening of the prograde limb is Bardeen's
  result that the code's gallery text describes correctly, and it is the
  feature finding 4.1 gets wrong under default settings, because the flat
  side *is* the prograde photon orbit hugging the horizon.

## 7. Premortem: "a render was wrong, why?"

Ranked by how likely it is to bite and how hard it is to notice.

1. A high-spin KerrBL render with default ε: the shadow's flat side is
   inside its true position, the prograde photon ring is missing or displaced,
   comparisons with the KS backend disagree on the prograde side only. Detect:
   run `scripts/validation/shadow_edge_bisect.py` for the scene's `a`, `r_cam`,
   ε.
2. A star field with the wrong brightness gradient (edge stars 4× faint
   relative to the disc, close-vantage star fields wrong by factors): 4.2.
   Detect: render a synthetic uniform catalogue; the star layer must be flat.
3. A Gaia constellation mirrored or rotated: the KS/BL twist (4.3) rotates the
   sky by `a/r_cam`; a handedness slip in `spatial_handedness` would mirror it.
   Detect: render a known asterism (Orion) with the hole out of frame and
   compare orientation with a planetarium view from the same RA/Dec.
4. "8000 K disc" that is actually 9400 K at peak (4.4); colour captions in
   the gallery are off by that much at `a/M = 0.998`.
5. Black speckles or arcs that are not shadow: NaN/step-exhausted rays painted
   black (4.8), especially near the BL polar axis or in animations that lower
   `--max-steps`.
6. Two backends that "almost" agree: any future test that compares KS and BL
   positions will keep needing loose tolerances until 4.3 is fixed, hiding real
   regressions.
7. Redshifted colour of the celestial *texture* uses the stationary emitter at
   `r = 15000`, fine; but if `--max-radius` is ever lowered to a few `r_s` for
   speed, the sky picks up a gravitational blueshift and the tube footprints
   are no longer flat-space solid angles. Guard `max_radius ≥ 1000 r_s`.
8. FOV surprises: someone sets `alpha = π/2` expecting 90° and gets a 180°
   stereographic image (4.5).
9. Camera inside the ergosphere with the default `StaticObserver`: fails
   loudly (good); a sphere object inside the ergosphere fails per pixel (4.9).
10. Volumetric disc at high spin: temperature edge and ISCO test at the wrong
    radius (4.6); small but systematic.

## 8. Extensions, with the concepts to look up

Ordered by value for a physics-first renderer.

1. Super-Hamiltonian BL integration (James et al. A.15; MTW §33.5;
   Levin & Perez-Giz 2008 for the same equations in the periodic-table
   context). Fixes 4.1 and keeps BL's speed advantage.
2. Geodesic deviation / Jacobi fields for the ray bundle (Pineault & Roeder
   1977; the Sachs optical scalars; DNGR appendix A.2). One extra 2×2 (or
   4×4) linear ODE per ray gives the exact lensing Jacobian, hence the star
   magnification and shear per pixel, the caustic structure, and it removes
   `magnification_cap`, corner subdivision and most curved-membership logic.
   Concepts: optical scalars (expansion, shear), Jacobi equation, caustics as
   zeros of the Jacobian.
3. Deflection-function star rendering for Schwarzschild (Riazuelo §IV–V):
   precompute `φ_∞(δ)` once per observer, place every star image (all orders)
   by inverting it, amplify by the local Jacobian. Exact positions, no tubes,
   two orders of magnitude faster.
4. Analytic Kerr lensing: Gralla & Lupsasca (2020, "Null geodesics of the Kerr
   exterior") give all null geodesics in closed form with elliptic integrals;
   Bardeen's critical curve `(ξ(r), η(r))` as a first-class test object
   (already used in the validation script). Concepts: Mino time, radial roots,
   the photon-shell parameter `r̃`.
5. Emission-time / light-travel-time consistency: everything here is
   stationary, so a rotating disc texture is a "frozen" pattern. For time
   dependence (hot spots, flares, camera motion), record the coordinate time
   at the hit and sample the source at `t_hit`. Concept: the Cunningham–Bardeen
   time-delay and "slow light" rendering.
6. Polarization: the Walker–Penrose constant transports polarization along
   Kerr null geodesics algebraically (no ODE); with a Novikov–Thorne disc and a
   simple scattering-atmosphere polarization model this gives EHT-style
   polarization maps. Concepts: Walker–Penrose constant, EVPA rotation.
7. A physical disc spectrum instead of a blackbody at `T(r)`: colour
   correction / spectral hardening factor (Shimura & Takahara 1995), limb
   darkening, and optionally the `kerrbb` parametrisation (Li et al. 2005) so
   `temperature` becomes `Ṁ` in physical units.
8. Cunningham transfer functions (Cunningham 1975) for thin-disc images:
   `g` and the redshift extremes per radius as an analytic check of the disc
   Doppler pattern, and as a fast path for equatorial discs.
9. Camera options: gnomonic (pinhole) projection for flat-sky comparisons,
   domemaster (Riazuelo §VII), a geodesic-orbit camera with DNGR's `β`, and
   DNGR-style analytic motion blur (their A.3.2).
10. Interior renders with the KS backend (already possible in principle): the
    Kerr–Schild chart is regular at the horizon; add a stop at the inner
    horizon and Riazuelo-style interior views become available.
11. Performance: analytic KS Christoffels (4.10), and an `atol/rtol` pair per
    state component once the state is typed (5.1).

## Appendix: reproducing the numbers

```sh
cargo build --release
# camera radius, spin a (r_s units), epsilon, backend, bracket (rad)
python3 scripts/validation/shadow_edge_bisect.py target/release/gr_raytracer 1000 0.499 1e-9 KerrBL 0.012
python3 scripts/validation/shadow_edge_bisect.py target/release/gr_raytracer 18   0.499 1e-5 KerrBL 0.6
python3 scripts/validation/shadow_edge_bisect.py target/release/gr_raytracer 18   0.0   1e-5 Schwarzschild 0.6
```

The Kerr–Schild backend can be bisected with the same script only at large
`r` (its `render-ray-at` tetrad is the Eulerian observer, whose velocity
relative to the ZAMO is `≈ r_s/r`, so the `b` conversion is exact only in the
limit); at `r = 1000` that adds ~10⁻³.

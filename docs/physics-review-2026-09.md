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
(187 tests). Both papers were read in full text (via alphaXiv); equation
numbers below refer to the arXiv versions (1511.06025v2, 1502.03808v2).

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
at high spin. A fourth, confined to the two flat-space geometries, was found
and fixed while checking the camera boost: they boosted a moving camera to
the reflected velocity or not at all (4.11). Everything else is documentation,
robustness or scope.

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
- Lorentz boost: identical to Riazuelo eq. (6) in Schwarzschild, Kerr and
  KerrBL, with the sign flip for the (−,+,+,+) geometries handled correctly.
  The two flat-space geometries had it wrong (finding 4.11, fixed on this
  branch).
- Momentum construction and frequency bookkeeping under a boost: the traced
  momentum is `p = N − u` with `N` a unit direction in the boosted tetrad, so
  `u·p = ∓1` by construction and the emitter energy `u_em·p` carries the
  whole shift. A boosted flat-space camera reproduces
  `ν_obs/ν_em = 1/(γ(1 − v cos θ_cam))` and the aberration formula to 10⁻¹²
  (regression test `boosted_camera_matches_special_relativity_in_flat_space`).
  The camera angles are applied to the reference tetrad before the boost,
  which is equivalent to Riazuelo's boost-then-rotate: the configured tilt is
  exactly the tilt in the camera's rest frame (checked numerically in
  Schwarzschild at `r = 4, 18, 50`).
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

### 4.11 Medium (flat space only): both Euclidean geometries boosted the camera wrongly — fixed on this branch

`EuclideanSpace::lorentz_transformation` had the last term of the boost as
`2 T^μ u_ν` instead of `2 u^μ T_ν`, which maps `T` to `2γT − u`: the camera
was boosted to **−v**. A camera configured to move toward a source saw it
redshifted and aberrated the wrong way. `EuclideanSpaceSpherical`'s boost was
the identity matrix, so a moving camera there was rendered as static. Neither
affects the black-hole geometries, which have their own correct
implementations, and every existing flat-space test used a static camera, so
nothing caught it. Both are fixed in this branch (one line and a ~20-line
generic boost respectively) and pinned by a test that checks the boosted
tetrad's time axis equals the velocity and that the rendered frequency ratio
and source directions match special relativity for every probed pixel.

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
- Star photometry (his §V.C, eqs. 52–53): `T_obs = T/(1+z)` and a bolometric
  `(1+z)⁻⁴` flux factor times the amplification `f = Ω/Ω*`, which he notes is
  "the convergence part of the optical scalar equations". This is the code's
  `g⁴` and `g·T` treatment; the only difference is the screen solid angle
  discussed in 4.2. He draws stars as truncated-Gaussian blobs whose size
  grows with brightness (the Akira Fujii diffusion-filter look), which is what
  `scripts/grade.py`'s bloom approximates after the fact.
- Sky texture under redshift (his §III.C): he keeps the hue and scales the
  intensity by an ad-hoc power law of `(1+z)⁻¹` for blueshift and an
  exponential for redshift, because a bitmap has no spectrum. The code's
  `beaming_exponent = 3` on bitmaps is the same kind of artistic choice; the
  README should say so as plainly as he does.
- He assumes photopic (colour) vision at all intensities, as the code does
  through the CIE 1931 curves.

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
- Disc (their §4 and A.6): the film disc is an artist's image on an infinitely
  thin plane with an optical-thickness map, "marginally optically thick",
  *not accreting*, at a uniform 4500 K blackbody; frequency shifts are applied
  as a temperature shift of the blackbody (exactly this code's `g·T`), then
  convolved with film sensitivity curves. Nolan and Franklin dropped the
  Doppler/gravitational colour and brightness shifts and lowered the spin from
  `a/M ≈ 1` to `0.6` because the flattened left shadow edge, the multiple disc
  images along it and the lopsided brightness (their Fig. 15c, "the hole's
  shadow barely discernible") were judged confusing. This code's default is
  the physical version with a Novikov–Thorne temperature profile, which is
  the better choice for a physics renderer. Their numbers for `a/M = 0.6`,
  `r_c ≈ 10M`-class views, disc speeds ≈ 0.55c, net frequency factors ≈ 1.5
  (approaching) and ≈ 0.4 (receding) including a ≈ 20 % gravitational
  redshift, are a cheap sanity check for the disc `g` field here.
- Their bug story (§3.4) is worth knowing: an early DNGR showed a
  "fingerprint-like" pattern of lensed stars inside the secondary critical
  curve on the flattened (prograde) side; Riazuelo's images did not, the
  discrepancy exposed a bug, and the corrected code also matched the SXS
  imaging code. That is the same region finding 4.1 affects, and the same
  method (cross-code comparison on the prograde limb) that would have caught
  it here; the two backends in this repository can play that role for each
  other once 4.3 is fixed.
- Caustics (§3.3–3.4): the primary caustic on the celestial sphere is a small
  astroid, secondary and tertiary caustics wrap the sphere once and six-plus
  times for a camera at `2.6M`, and the number of images of a sky patch
  between critical curves (three, then eight) is independent of the camera's
  velocity because caustics depend on the camera's location only. Multiple
  nested images along the flat edge are therefore real, not artifacts.
- Camera velocity (§3.5): for the same location at `r = 2.6M` the static,
  ZAMO and geodesic-orbit cameras see wildly different shadows through
  aberration alone (the static camera's shadow exceeds half the sky). The
  `CameraVelocityConfig` choice is not cosmetic; the gallery captions should
  state which observer each close vantage uses.
- Integration (A.4–A.5): RKF with "empirically determined tolerances", tighter
  on the ray's position than on the beam's shape; renders near the shadow
  dominate run time. Same integrator family as here, and the same place where
  the tolerance matters most.
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
9b. "Too many" nested disc or star images along the flat edge of a high-spin
   shadow are real (DNGR §3.4: three images between the primary and secondary
   critical curves, eight between the secondary and tertiary); a *missing*
   set of them on that side is the symptom of 4.1.
10. Volumetric disc at high spin: temperature edge and ISCO test at the wrong
    radius (4.6); small but systematic.
11. A fast camera in a black-hole scene: the redshift and star *positions*
    are right, but the star *brightness* field is not, because the tube
    magnification ignores the camera's aberration (4.2). Before this branch,
    a fast camera in flat space was also aberrated the wrong way (4.11).
12. A camera inside the horizon, an over-extremal Kerr scene, or any future
    metric without a guarding stop: rays that asymptote to a coordinate
    singularity (the past horizon in ingoing Kerr–Schild, the ring for
    |a| > M) collapse the step to `H_MIN`, are then accepted regardless of
    error, and either explode to r ≈ 10²⁴ or burn the whole step budget.
    Today the `horizon_epsilon` stop hides this (§11).

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
4. Analytic Kerr lensing: Gralla & Lupsasca, "Null geodesics of the Kerr
   exterior" (arXiv:1910.12881) give all null geodesics in closed form with
   elliptic integrals, and Gralla, Lupsasca & Marrone (arXiv:2008.03879) the
   photon-ring shape; Bardeen's critical curve `(ξ(r), η(r))` as a first-class
   test object (already used in the validation script). Concepts: Mino time,
   radial roots, the photon-shell parameter `r̃`.
4b. Cross-code validation as a habit: GYOTO's validation paper
   (arXiv:1605.04195) lists the standard test problems (shadow size, thin-disc
   spectra, redshift of circular orbits); GYOTO 2.0 (arXiv:2311.18802) and the
   2026 TARTARUS code (arXiv:2609.08989) are open references to compare disc
   images against. DNGR's own bug was found exactly this way.
5. Emission-time / light-travel-time consistency: everything here is
   stationary, so a rotating disc texture is a "frozen" pattern. For time
   dependence (hot spots, flares, camera motion), record the coordinate time
   at the hit and sample the source at `t_hit`. Concept: the Cunningham–Bardeen
   time-delay and "slow light" rendering.
6. Polarization: the Walker–Penrose constant transports polarization along
   Kerr null geodesics algebraically (no ODE); with a Novikov–Thorne disc and a
   simple scattering-atmosphere polarization model this gives EHT-style
   polarization maps. Concepts: Walker–Penrose constant, EVPA rotation;
   GYOTO 2.0 (arXiv:2311.18802) for a full polarized transfer implementation,
   and arXiv:2407.14897 for polarized hot-spot images as a small worked
   example.
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

## 9. The disc in depth: energy, temperature, blackbody, redshift

This section follows one photon from the disc to the pixel and checks each
link against the literature: Shakura & Sunyaev (1973), Novikov & Thorne
(1973), Page & Thorne (1974), Cunningham (1975, 1976), Luminet (1979),
Li, Zimmerman, Narayan & McClintock 2005 (the `kerrbb` paper,
astro-ph/0411583), Gralla, Holz & Wald 2019 (arXiv:1906.00873) and
Johnson et al. 2020 (arXiv:1907.04329).

### 9.1 The chain as implemented, and what is right about it

1. **Energy release** (`temperature.rs`): the local flux is Page & Thorne's
   eq. (15), `F(r) = −(Ṁ/(4π√−g)) · Ω_{,r}/(E − ΩL)² · ∫_{r_ISCO}^{r} (E − ΩL) L_{,r} dr`,
   with `√−g = r` for the (t, r, φ) equatorial metric and the zero-torque
   inner boundary at the ISCO (η = 0 in `kerrbb`'s notation). The code has the
   same integrand, prefactor structure and lower limit; only the constant
   differs, and the constant is calibrated away by `Ṁ`. `E(r)`, `L(r)`, `Ω(r)`
   are the same circular-orbit functions the redshift uses, which is the
   consistency Page & Thorne require. Independent check (Simpson rule,
   analytic `Ω_{,r}`, 2·10⁴ intervals): the code's `T(r)/T(1.5 r_ISCO)`
   agrees with Page–Thorne to 4·10⁻⁴ at `a = 0`, 4·10⁻³ at `a/M = 0.5`,
   and 10⁻³ at `a/M = 0.998` for `r ≥ 1.1 r_ISCO` (see 9.2 for the rim).
   The far-field `T ∝ r^{−3/4}(1 − √(r_ISCO/r))^{1/4}` Shakura–Sunyaev law
   is also pinned by an existing test.
2. **Temperature**: `T_eff = (F/σ)^{1/4}`, local blackbody, both faces emit
   equally, opaque (`alpha = 1`). This is the NT/Page–Thorne assumption and
   Luminet's; `kerrbb` adds a colour correction (9.3).
3. **Emitter frame**: the Keplerian circular geodesic at the BL radius of the
   hit, via Killing coefficients `(u^t, u^φ)`. Correct above the ISCO; the
   Luminet closed form `1 + z = (1 − 3M/r)^{−1/2}(1 + Ω b)` is reproduced by
   an existing test, and `kerrbb`'s appendix F gives the same formula (their
   eq. F1) as the one they had to correct in XSPEC's `GRAD`.
4. **Frequency shift**: `g = (u_cam·p)/(u_disc·p)`, chart-invariant and
   sign-invariant.
5. **Intensity**: `I_λ^obs(λ) = g⁵ B_λ(gλ, T) ≡ B_λ(λ, gT)`. This is the
   Liouville invariance of `I_ν/ν³` written for `I_λ`; integrated over all
   frequencies it gives the bolometric `g⁴` (Luminet's `(1+z)^{−4}`, Gralla
   et al. eqs. 10–11), and it is exactly DNGR's "temperature shift of the
   blackbody" (their A.6). The XYZ integral is taken over 380–830 nm of the
   *shifted* spectrum, which is the right band-limited quantity: a strongly
   redshifted patch goes dark because its light leaves the visible band, not
   only because of `g⁴` (Riazuelo makes the same point in his Appendix B).
6. **Beaming**: none extra; `beaming_exponent = 0.0` in the blackbody scenes.
   Correct, since `g⁵` already contains the Doppler boost.

Verdict: the physics chain is the standard thin-disc one and it is
implemented consistently. Everything below is either a numerical detail or
physics the standard model leaves out.

### 9.2 Disc-specific findings

**(a) Peak-temperature calibration** (finding 4.4): +18 % at `a/M = 0.998`.

**(b) Inner-rim resolution of the temperature LUT.** Near the ISCO
`L_{,r} → 0`, so `F ∝ (r − r_ISCO)²` and `T ∝ (r − r_ISCO)^{1/2}`: a square-root
cusp. The LUT is 1000 points uniform in `r` from `r_ISCO` to `outer_radius`,
interpolated linearly in `T`. Measured ratio code/Page–Thorne:

| a/M | outer radius | 1.01 r_ISCO | 1.02 | 1.03 | 1.05 | 1.10 | ≥ 1.2 |
|---|---|---|---|---|---|---|---|
| 0.998 | 12 (gallery) | 0.71 | 0.99 | 0.98 | 0.995 | 0.999 | 1.000 |
| 0.998 | 40 | 0.33 | 0.48 | 0.61 | 0.85 | 0.98 | 0.999 |

With the gallery's outer radius the error is confined to the innermost 2 %
of the rim, where the disc is nearly dark, so it is invisible; with a wide
disc (outer 40) the whole inner rim is 15–50 % too cool. Fix: tabulate `F`
(or `T⁴`, which is smooth and quadratic at the rim) and take the fourth root
after interpolation, or grid in `s = √(r − r_ISCO)`. One line either way.

**(c) No physical units.** `temperature` is a display target; the disc has no
`M` or `Ṁ`. For a physically parametrised disc use `kerrbb`'s normalisation
(their eqs. 16–17: the spectrum depends on `f_col`, `Ṁ^{1/4} M^{1/2}` and
`M²/D²`), or the Shakura–Sunyaev scaling `T_max ∝ M^{−1/2} Ṁ^{1/4}` (Frank,
King & Raine, eq. 5.43). A scene could then say "10 M☉, 10 % Eddington" and
get kelvin out; the visible-band colour of a stellar-mass disc (10⁷ K) and of
an AGN disc (10⁴–10⁵ K) differ enormously, and the gallery's 4000–20000 K
correspond to the AGN/quasar regime.

**(d) Physics the NT chain omits, in order of visual impact:**

1. *Colour correction / spectral hardening*: the disc atmosphere is
   scattering-dominated, so the emergent spectrum is a diluted blackbody
   `I_ν = f_col^{−4} B_ν(f_col T_eff)` with `f_col ≈ 1.5–1.9` (Shimura &
   Takahara 1995; Davis et al. 2005 favour 1.5–1.6; `kerrbb` uses 1.7). The
   code's colours are `f_col = 1`. For an 8000 K disc, `f_col = 1.7` moves the
   hue from amber to white-blue at the same bolometric flux. One-line change
   in `BlackBodyMapper` plus a scene parameter.
2. *Limb darkening*: `kerrbb` eq. (D20), the Chandrasekhar electron-scattering
   law `I(μ) ∝ 1 + 2.06 μ` with `μ = cos` of the emission angle measured in
   the emitter frame. The code emits isotropically. The emission angle is
   available for free: `μ = (p·n)/(p·u_disc)` with `n` the unit normal in the
   disc frame. Matters most for the grazing views in the vantage series.
3. *Returning radiation* (Cunningham 1976; `kerrbb` §3.1): light from the
   inner disc lensed back onto the disc. `F_in` stays finite at the ISCO
   where `F_0 → 0`, so the dark ISCO rim in the gallery renders is partly an
   artefact of ignoring it; `kerrbb` finds the effect is equivalent to raising
   `Ṁ` by ≈ 1.7 at `a/M = 0.999`. This renderer can compute it exactly by
   tracing rays *from* the disc, which is the same machinery as rendering.
4. *Inner-edge torque* (Agol & Krolik 2000; `kerrbb` eq. 2, parameter η):
   MHD discs have `η ~ 0.2`; it adds an `F ∝ r^{−7/2}` component and brightens
   the inner disc. A scene parameter and one extra term in the flux.
5. *Plunging-region emission*: NT assumes none inside the ISCO and the
   code refuses (`BelowRISCO`). Gralla, Holz & Wald show the size of the
   central dark area is set by the lensed inner edge of the *emission*, not
   by the critical curve: for Schwarzschild it is `b ≈ 2.9M` if emission
   reaches the horizon versus `5.2M` for the "shadow". Chael et al. 2021
   (arXiv:2106.00683) call the former the *inner shadow*. Allowing an
   emissivity profile inside the ISCO (with the infalling, not Keplerian,
   four-velocity) would let the renderer show it.

**(e) Photon-ring windings caption.** The gallery says successive windings
are "≈ 23× thinner". Using Johnson et al. eq. (29), `γ = √(R''(r_γ)/2) · G_θ`,
at the two equatorial limb points of an edge-on observer the exponent per
polar half-oscillation (which is what counts successive images of an
equatorial disc, since each is a new equatorial crossing) comes out `γ = π`
for every spin, so `e^π ≈ 23` holds there. Two caveats for the caption:
`γ` varies around the ring at other angles (their Fig. 6), and the number of
*azimuthal* turns per image is very spin dependent: per azimuthal half-turn
the demagnification at `a/M = 0.998` is only `e^{0.18} ≈ 1.2` on the prograde
limb versus `e^{4.1} ≈ 59` on the retrograde limb. That is DNGR's remark that
frame dragging "moves the critical curves outward from the shadow's flattened
edge", and it is why the prograde side of that image shows the nested disc
images spread out while the retrograde side stacks them. It is also, again,
the region finding 4.1 affects.

**(f) Lensing ring versus photon ring.** Gralla, Holz & Wald's decomposition
(direct image `n = 1`, lensing ring `n = 2` at `5.02M < b < 6.17M` in
Schwarzschild, photon ring `n ≥ 3` within `5.19–5.23M`) explains what the
zoomed gallery frames show: the broad second image of the disc's far side is
the lensing ring and carries a few percent of the flux; everything inside it
is the photon ring and carries `e^{−π}` of that per order. The caption's
"surface brightness is conserved, so each stays at full luminance" is right
for an opaque disc; Luminet's 1979 remark that an opaque disc occults most of
its own secondary image applies here too, and is why the secondary image
appears only as a thin band hugging the shadow.

**(g) Small consistencies.** The volumetric disc evaluates `T` at the
cylindrical radius (4.6) and applies an extra artistic `(T/T_ref)⁴` on top of
the Planck amplitude, which double-counts the temperature dependence of
brightness; document it as artistic or remove it. The sphere object uses a
static emitter (4.9). The celestial `celestial_temperature` blackbody sky is
fine.

### 9.3 Sanity numbers you can check against the renders

- Peak of `T(r)`: `1.59 r_ISCO` at `a = 0` (the classical `(49/36) r_in`
  is for the Newtonian `(1 − √(r_in/r))/r³` law, the relativistic value is
  slightly different), `1.56` at `a/M = 0.5`, `1.28` at `a/M = 0.998`.
- Frequency factors at the disc, `a/M = 0.6`, viewed near the plane: ≈ 1.5
  approaching, ≈ 0.4 receding, including ≈ 20 % gravitational redshift
  (DNGR §4.1.2). An `--epsilon`-independent check of the redshift field.
- `a/M = 0.998` is Thorne's (1974) spin-up limit for a disc-fed hole, so the
  gallery's "essentially extremal" is also the astrophysically maximal case;
  `a/M = 0.9995` (`a = 0.49975`) is beyond it.
- Shadow edges (Section 3.1) are the disc-independent part of the same
  geometry.

### 9.4 Papers for the disc, beyond the two in the README

- Page & Thorne 1974, ApJ 191, 499: the flux formula implemented.
- Thorne 1974, ApJ 191, 507: the 0.998 spin limit.
- Cunningham 1975, ApJ 202, 788: transfer functions and the `g`-distribution
  of a Kerr disc; Cunningham 1976, ApJ 208, 534: returning radiation.
- Luminet 1979, A&A 75, 228 and Luminet's history (arXiv:1902.11196):
  the bolometric image, isoradial curves, the opaque-disc occultation.
- Shimura & Takahara 1995, ApJ 445, 780; Davis et al. 2005: `f_col`.
- Agol & Krolik 2000, ApJ 528, 161: inner torque.
- Li et al. 2005 (`kerrbb`, astro-ph/0411583): the complete modern thin-disc
  ray-tracing spec, with every formula in its appendices.
- Gralla, Holz & Wald 2019 (arXiv:1906.00873): direct / lensing / photon
  ring, inner dark area set by emission not by the critical curve.
- Johnson et al. 2020 (arXiv:1907.04329): Lyapunov exponents and subring
  demagnification around the Kerr ring.
- Chael, Johnson & Lupsasca 2021 (arXiv:2106.00683): the inner shadow.

## 10. Would a GPU help, given the early-stopping logic?

Yes, by an order of magnitude or more for the geodesic pass, and the stop
conditions (horizon, celestial sphere, NaN, trapped orbit, object hit) are not
the obstacle. What decides it:

1. **Divergence from early stopping is manageable.** Rays end at different
   steps, so lockstep threads idle; every GPU geodesic code has this (GRay,
   Odyssey) and still reports 10–100× over CPU. The standard cure is a
   persistent-thread work queue with ray compaction: a finished thread takes
   the next unstarted pixel, so warp cost is bounded by the longest *active*
   ray. Rays near the critical curve take ~100× more steps than background
   rays, so a naive one-pixel-per-thread kernel would idle most of the device.
2. **The adaptive integrator diverges more than the stop tests.** RKF45 with
   per-ray retry loops branches every step. Either make accept/reject
   branch-free (compute both, mask), or use a fixed-step high-order scheme
   with a per-ray step chosen from the local curvature scale, as most GPU
   codes do.
3. **The CPU architecture cannot be ported as is.** `color_of_ray`
   integrates the whole trajectory into a `Vec<Step>` (up to `max_steps` × 64
   bytes per ray) and intersects objects afterwards; the volumetric disc
   marches over `remaining_steps`. A GPU version must fuse integration,
   intersection and shading in one loop (test the disc crossing each step,
   march the volume as it goes, accumulate radiance front to back, stop when
   opaque). That restructuring is the real work, and it would also speed up
   the CPU path and cut its memory traffic.
4. **Precision decides the API.** Finding 4.1 shows near-critical geodesics
   are ill-conditioned even in f64; f32 is not usable near the shadow. That
   rules out wgpu/WGSL (no f64). From Rust: CUDA via `cudarc`/`cust`, or
   Vulkan compute with the `Float64` capability. Consumer GPUs run f64 at
   1/32–1/64 of f32 rate, roughly a 16–32-core CPU under rayon; data-centre
   parts (A100/H100 class) are where the 10–50× is.
5. **Cheaper wins first.** Analytic Kerr–Schild Christoffels (the backend
   currently spends 36 metric evaluations per right-hand side on finite
   differences) give several× on the CPU today, and the BL super-Hamiltonian
   form is both the correctness fix and GPU-friendly (polynomial right-hand
   sides). The star gather and adaptive supersampling stay on the CPU either
   way; they are not the bottleneck.

## 11. Is the Kerr–Schild chart numerically viable inside the horizon?

Question raised after the review: if a camera is placed inside the black
hole, does the Kerr–Schild (KS) chart hold up numerically, or does one need
an atlas with chart transitions along the geodesic?

**Answer: KS is viable for every ray that crosses the future horizon, in
either parameter direction, at the default tolerance. The rays that do not
cross are the ones that asymptote to the *past* horizon, which ingoing KS
does not cover; those need a stop condition, and the current one only works
from outside.** No atlas is needed for interior renders.

### 11.1 Experiment

Rays were seeded from the Eulerian tetrad (`Kerr::get_tetrad_at`, which is
regular inside since the KS time slices stay spacelike there) as
p = ±e_t + n, with n a unit spatial direction in the tetrad. The sign −e_t
is the past-directed ray the camera actually traces; +e_t is its
future-directed twin. Each ray was integrated with the production `rkf45`
and `KerrSolver` right-hand side, with the horizon stop disabled, until it
reached r = 60 ("escaped"), an inner stop just outside r₋ (or r = 0.05 for
a = 0), a NaN, or 400 000 steps. Spins a/r_s ∈ {0, 0.4, 0.499}, camera at
r₀ between the horizons, tolerances ε ∈ {10⁻⁵, 10⁻⁹}. Quality measures: the
Hamiltonian constraint |ΔH|/|p|² (should stay ≈ 0), and drift of E = −p_t
and L_z. The scratch test is kept in the session scratchpad
(`ks_interior_scratch.rs`); it is 120 lines and could become a regression
test.

### 11.2 Results

Past-directed rays (what the camera sees) from r₀ inside, that leave through
r₊:

| a/r_s | direction | ε | steps to r₊ | steps to r = 60 | |ΔH|/|p|² | ΔL_z |
|---|---|---|---|---|---|---|
| 0 | outward | 1e-5 | 5 | 42 | 1.5e-12 | 0 |
| 0 | tangential | 1e-5 | 7 | 104 | 9.8e-6 | 1.3e-6 |
| 0.4 | outward | 1e-5 | 5 | 45 | 8.9e-8 | 5.7e-8 |
| 0.4 | oblique | 1e-5 | 5 | 59 | 1.1e-6 | 6e-10 |
| 0.499 | outward | 1e-5 | 3 | 46 | 1.8e-7 | 7.1e-8 |
| 0.499 | oblique | 1e-9 | 4 | 100 | 4.5e-10 | 7e-13 |

Future-directed rays falling in from r = 10 through r₊ to the inner stop
behave the same way (18–435 steps total, |ΔH|/|p|² ≤ 2e-4 at ε = 1e-5,
≤ 1e-7 at ε = 1e-9). The step size never drops below 1e-5 on any crossing
ray; the horizon is simply not there numerically, which is the whole point
of the chart.

The other family: past-directed rays whose spatial direction points *inward*
in the interior tetrad. Traced to the past they move outward toward r₊ but
never cross it; coordinate time runs to −∞ (t ≈ −40 to −160 at abort) and
the step collapses to `H_MIN`:

| a/r_s | ε | outcome |
|---|---|---|
| 0 | 1e-5 | tunnels through r₊ after ≈1000 steps at H_MIN, then explodes: r → 10²⁴, |ΔH|/|p|² ≈ 10⁷² |
| 0 | 1e-9 | 400 000 steps burnt at h = 10⁻¹² sitting on r = 1.000 |
| 0.4 | 1e-5 | NaN after ≈9000 steps |
| 0.499 | 1e-5 | NaN after ≈72 000 steps |

Exactly the same happens to the ordinary exterior shadow rays if the
`r ≤ r₊ + horizon_epsilon` stop is removed. These rays are physical: they
are the ones that would show the collapsing star (or the white-hole region
in the eternal solution) and they are correctly rendered as "frozen at the
horizon" black. The chart is not at fault; the stop condition is. Ingoing
KS covers the future horizon only, so the past horizon shows up as t → −∞
at r = r₊, and no step-size controller can integrate through it.

### 11.3 Consequences for the code

1. **Interior camera needs a two-sided horizon stop.** Today's
   `inside_horizon` fires for every point with r ≤ r₊ + ε, so an interior
   camera would stop every ray at step 0. Replace it by a *stall* criterion
   that is chart-independent: N consecutive steps at `H_MIN` (or
   Δr < ε over the last N steps while |r − r₊| < ε) ⇒ `HorizonReached`. A
   plain band |r − r₊| < ε is not safe, because crossing rays step through
   r₊ at finite dr/dλ and can land inside the band by chance, giving black
   speckles.
2. **Accepting a step at `H_MIN` regardless of error is dangerous**
   (`runge_kutta.rs:148ff`). It is what turns the stall into a tunnel and
   an explosion to r ≈ 10²⁴ at the default ε. The horizon stop masks this
   today. It is unmasked for any geometry without a guarding stop:
   over-extremal Kerr (`inside_horizon` returns `false` for |a| > M, so
   rays approaching the ring have no stop at all), and any future wormhole
   or space-only metric. Recommendation: treat "error still above ε at
   H_MIN" as a stop reason, not as an accepted step. This belongs in the
   premortem list (§7) as item 12.
3. **BL cannot do this at all.** The ZAMO tetrad has an imaginary lapse
   inside r₊ and the Mino-time solver has Δ → 0 at both horizons. So an
   interior camera is a KS-only feature unless an atlas is added; the
   hybrid "BL outside, KS inside r_switch" is a speed optimisation, not a
   requirement.
4. **Between r₋ and r₊ every past-directed ray moves outward** (r never
   decreased on any of the past-directed runs), so the inner horizon needs
   no stop for a camera in that region. A camera inside r₋ is a different
   problem (Cauchy horizon, mass inflation) and out of scope.

Concepts to look up: ingoing versus outgoing Eddington–Finkelstein
coordinates and which horizon each covers; the Penrose diagram of
Schwarzschild and Kerr (regions I–IV); the Cauchy horizon at r₋.

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

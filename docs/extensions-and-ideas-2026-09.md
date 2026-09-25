# Extensions and ideas (September 2026)

Companion to `physics-review-2026-09.md`. That document holds the findings.
This one collects the design discussion that followed the review: feature
ideas, a space-only differential-geometry mode, rendering without light
sources, which non-Euclidean geometries are tractable, and strange metrics
worth rendering. Units: r_s = 1 unless stated.

## 1. Future features

**Disc physics**

- Colour-correction factor f_col ≈ 1.7 (kerrbb style): T_col = f_col·T_eff, radiance scaled by f_col⁻⁴.
- Chandrasekhar limb darkening as a function of the emission angle in the disc rest frame.
- Returning radiation (disc light re-emitted after falling back on the disc).
- Inner torque parameter instead of the zero-torque ISCO boundary.
- Plunging-region emission and the resulting inner shadow.
- Thick/slim discs and Polish-doughnut tori with self-shadowing.
- Non-uniform temperature LUT with denser sampling near the (r − r_ISCO)^½ rim cusp.

**Light transport**

- Ray-tube or Jacobi-field magnification instead of the current cap heuristic. *(pick)*
- Polarisation transport (Walker–Penrose constant in Kerr, or parallel transport).
- Spectral rendering with a camera response instead of bolometric colour mapping.
- Higher-order images: count equatorial crossings per ray and expose the order as a channel.

**Spacetimes**

- Super-Hamiltonian first-order BL solver (closes the prograde-edge gap at high spin). *(pick)*
- Reissner–Nordström, Kerr–Newman, Johannsen–Psaltis.
- Wormholes: Ellis, Morris–Thorne, the three-parameter DNGR wormhole.
- Naked-singularity Kerr (a > M) and the r < 0 sheet.

**Camera and observer**

- Expose the stereographic half-angle α in the config (currently π/4 hard-coded).
- Pinhole and equirectangular projections.
- Observer on a circular geodesic orbit or in free fall, via the existing tetrad + boost machinery.
- Motion blur / exposure time for time-dependent scenes.

**Stars and sky**

- Star layer radiance g⁴·ΣF/Ω_sky with no coordinate solid-angle factor. *(pick)*
- Gaussian star filtering (DNGR) against flicker under magnification.
- Real catalogue (Hipparcos or Gaia subset) with colour from B−V.
- Doppler colour shift of the sky under camera boosts.

**Numerics and tooling**

- Relative + absolute error norm in RKF45 with per-component scales.
- Dual-number automatic differentiation for metric derivatives (removes the 36-evaluation stencil).
- Shadow-edge bisection (`scripts/validation/shadow_edge_bisect.py`) as a CI regression test.
- Fused integrate-and-shade loop; GPU only after that (review §10).
- Progressive preview, tiled rendering, resume from partial output.

## 2. Why fusing integrate and shade also speeds up the CPU path

Today a ray is integrated to completion, its trajectory stored, and a second
pass walks the stored states to find crossings and shade. Fusing the passes
helps the CPU independently of any GPU:

- **Memory traffic.** The trajectory buffer is written once and read once and does not stay in L1 at thousands of steps per ray. The fused loop keeps one state in registers.
- **Early exit.** A ray hitting an opaque disc at step 300 currently still runs to the step or radius limit. Fusing the hit test stops it at the first hit. For disc-heavy views this is the largest saving.
- **Allocation.** No per-ray `Vec` of states.
- **Enables batching.** A self-contained loop with no growing buffer is what SIMD lanes or GPU threads need; the buffer is what blocks that today.

Concepts: loop fusion, arithmetic intensity, the roofline model.

## 3. A general differential-geometry (space-only) raytracer

**Physics.** An *ultrastatic* spacetime ds² = −dt² + h_ij dx^i dx^j has null
geodesics that project to arc-length geodesics of the 3-metric h with p_t
constant. Any Riemannian 3-manifold can therefore be fed to the existing
Hamiltonian integrator by lifting it with g_tt = −1, g_ti = 0. For a static
but non-ultrastatic metric, the *optical (Fermat) metric* h_opt = h/(−g_tt)
gives a 3-metric whose geodesics are exactly the spatial light-ray tracks.
Same image geometry, different radiometry.

**Code changes**

| Area | Current | Space-only mode |
|---|---|---|
| Geometry trait | 4-metric + solver | `Riemannian3` adaptor lifting any h_ij to an ultrastatic 4-metric |
| Metric derivatives | central differences, 36 evals/RHS | dual-number AD (`num-dual` or a hand-rolled dual type); enables user-supplied metric closures |
| Stop logic | horizon, escape radius, max steps | arc length, winding / domain-crossing count, hit; escape radius optional |
| Point type | 4 components | keep 4; let t carry arc length |
| Identification hook | none | optional `wrap(&mut state)` after each accepted step (deck transformation) |
| Observer frame | tetrad + boost | orthonormal triad of h; boost with v = 0 |

**Slider between flat and Schwarzschild.** In isotropic coordinates the
spatial Schwarzschild metric is conformally flat, h = (1 + M/2ρ)⁴ δ. The
PPN form h = (1 + 2γM/ρ) δ, γ ∈ [0, 1], interpolates from flat space to GR
spatial curvature; at γ = 1 it must reproduce half the light-bending angle
(the Eddington split), a good first test of the adaptor.

**Architecture.** The geometry/integrator trait boundary is in the right
place. The identification hook belongs in the solver, not the geometry:
the integrator re-expresses the state after the wrap while the metric is
unchanged.

Concepts: ultrastatic spacetimes, Fermat's principle in GR (Gibbons–Werner
optical metric), conformal flatness, forward-mode AD with dual numbers,
deck transformations and fundamental domains.

## 4. Rendering without light sources

**Shading**

- *Headlight shading*: brightness max(0, n·(−v)) with both vectors in the local orthonormal frame at the hit.
- *Distance fog by arc length*, not coordinate distance. The most effective depth cue in curved / multiply connected spaces.
- *Surface brightness is conserved* (étendue), so geometric renders have no natural inverse-square falloff. Fog must supply depth.

**Content**

- Lattices of tripods, spheres or wireframe grids filling the space; repetition makes curvature and identification readable.
- Faint translucent fundamental-domain walls to show where wrapping happens.
- Put the camera inside the structure.

**Colour layers computed from the ray**

| Layer | Shows | Cost |
|---|---|---|
| Image order | domain / equatorial crossings before the hit | one counter |
| Arrival time / arc length | isochrones of the light front | one accumulator |
| Magnification | det of the geodesic-deviation Jacobian; caustics as natural "lighting" | Jacobi equation alongside the ray (6 extra components in 3D, 12 in 4D) |
| Holonomy | rotation of a parallel-transported frame camera → hit | transport a frame |
| Curvature at the hit | Ricci / Kretschmann scalar | second derivatives, cheap with AD |

The magnification layer is also the correct replacement for the current
magnification cap in the black-hole renders.

## 5. Thurston geometries and ray marching

The relevant property is *homogeneity*, not flatness. Sphere tracing needs
closed-form geodesics and a lower bound on the distance to the nearest
surface; homogeneity supplies both.

| Geometry | Geodesics | Distance | Marching |
|---|---|---|---|
| E³ | lines | closed form | standard |
| S³ | great circles | closed form | works |
| H³ | closed form (hyperboloid model) | closed form | works |
| S²×R, H²×R | product, closed form | closed form | works |
| Nil | closed form (helices) | needs a solve | numeric geodesics + conservative bound |
| SL(2,R)~ | closed form | partial | works with care near the fibre |
| Sol | no closed form | no closed form | numerical integration only |

Five of eight are a plain yes; Nil and SL(2,R)~ need numerical geodesics
but keep usable bounds; Sol falls back to exactly what this repository does.

**Idea to take from ray marching: an SDF step bound.** The RKF45 step is
chosen from the error estimate alone, so thin objects can be stepped over.
Bound the step by the distance to the nearest surface:

```rust
// after computing the RKF45 proposal h_err:
let d = scene.signed_distance(&state.position);
let speed = state.spatial_speed();      // |dx/dλ| in the 3-metric
let h_geom = safety * d / speed;        // safety ≈ 0.8
let h = h_err.min(h_geom).max(H_MIN);
```

Hit testing becomes a sign check of the SDF, thin discs are never skipped,
and the step count near objects drops because a globally small step is no
longer needed. Works unchanged for the Kerr renders.

Concepts: geometrisation and the eight model geometries, sphere tracing
(Hart 1996), signed distance functions, the hyperboloid model, Heisenberg
group geodesics.

## 6. Other geometries, including the inside of a 3-torus

**Quotients of flat space.** The 3-torus is flat; what is interesting is
topology. A ray leaving through one face re-enters through the opposite one
with the same direction, so the camera sees a periodic lattice of copies of
the scene, including itself from behind at the box length. Implementation:
the identification hook reduces coordinates modulo the box after each step;
the metric is the identity. Cheapest new geometry, and it exercises the
stop logic, image-order counting and fog that the curved ones need. There
are ten closed flat 3-manifolds; the half-turn and quarter-turn spaces show
rotated copies, the non-orientable ones mirrored copies, and Hantzsche–Wendt
has no translational symmetry.

**Quotients of curved spaces.**

- Spherical: lens spaces L(p,q), Poincaré dodecahedral space. Finite volume; the antipodal image of the camera is magnified without bound.
- Hyperbolic: Seifert–Weber space, figure-eight knot complement. Copies multiply exponentially with distance; fog is mandatory.
- Quotients of Nil and Sol: where the current integrator is the right tool.

**Bookkeeping a quotient needs.** Fundamental-domain membership test and
deck transformation; the *injectivity radius* as a hard step bound (the SDF
bound covers it if the walls are in the SDF); an image-order counter.

**Spatial slices of black-hole spacetimes.** The t = const slice of
Schwarzschild (Flamm's paraboloid) is a Riemannian 3-manifold. Its geodesics
are not the light rays (those follow the optical metric), which makes a good
teaching pair.

Concepts: covering spaces, the Bieberbach classification in 3D, lens spaces,
injectivity radius, Flamm's paraboloid.

## 7. Strange curved metrics worth rendering

**Throats and wormholes**

- *Ellis wormhole* (first pick): ds² = −dt² + dl² + (l² + b₀²) dΩ². Ultrastatic, one parameter, no horizon; a second sky through a ring, Einstein ring and photon sphere on the throat. Stop logic: escaped on either side.
- *Morris–Thorne* general form; the DNGR three-parameter wormhole (throat radius, length, lensing width).

**Causality breakers**

- *Gödel universe* (pick): rotating dust with closed timelike curves; a ring beyond which the sky is mirrored and rotated. Homogeneous with exact geodesics, so it tests the Hamiltonian solver in a metric with g_tφ ≠ 0 everywhere.
- *Alcubierre warp bubble*: flat inside and out, thin curved shell; strong blueshift in front, redshift behind. Time-dependent, so the solver needs the ∂_t g terms it currently drops.
- *Taub–NUT*: gravitomagnetic monopole; rays near the axis twist along the Misner string.

**Exotic but astrophysically discussed**

- Kerr with a > M and the r < 0 sheet (naked ring singularity). BL solver handles it once the horizon stop is removed.
- Boson stars, gravastars: horizonless, Schwarzschild exterior; shadow with a bright core.
- Janis–Newman–Winicour naked singularity (reduces to Schwarzschild).
- Johannsen–Psaltis and other deformed Kerr metrics (EHT shadow tests).
- Hairy black holes (Kerr with scalar hair): shadows far smaller than Kerr's for the same mass.
- Superposed binaries (Majumdar–Papapetrou extremal pair): two shadows with eyebrow-shaped higher-order images.

**Cosmological**

- de Sitter / anti-de Sitter: constant curvature; AdS reaches its boundary in finite affine parameter, dS has a cosmological horizon.
- Kottler (Schwarzschild–de Sitter): two horizons.
- FLRW and Milne: expansion redshift without gravity; Milne is a patch of Minkowski and a good sanity test.

**Spatial only (space-only mode)**

- Nil and Sol as pure 3-metrics; Flamm's paraboloid.
- *Gradient-index lenses* (pick): Maxwell fish-eye n = 2/(1 + r²), Luneburg n = √(2 − r²), Eaton n = √(2/r − 1); optical metric n²δ. The fish-eye is the stereographic image of S³ (every ray a circle, perfect antipodal imaging); the Eaton lens is a perfect retroreflector. Easiest striking cases and analytic tests of the optical-metric adaptor.
- Cone spaces (cosmic-string metric): double images without distortion.
- Random metrics h = (1 + εf(x)) δ with band-limited noise f: weak-lensing shimmer and a stress test for adaptive stepping.

**Recommended order.** Ellis wormhole → Eaton or fish-eye lens → Gödel. The
3-torus fits in at any point as a one-evening exercise for the hook.

## 8. Stars, filters and folds

**Sizes are a camera question.** The camera is stereographic with a
hard-coded 45° half angle (90° frame height). Shadow radius is verified
against Bardeen to < 1 %. For a static camera at distance D (units r_s),
Schwarzschild:

| D | shadow diameter | Einstein ring diameter | gap: shadow edge → 2nd ring |
|---|---|---|---|
| 5 | 60 % of frame height | 100 % | 0.5° |
| 10 | 30 % | 64 % | 0.34° |
| 18 | 17 % | 46 % | 0.22° |
| 50 | 6 % | 26 % | 0.08° |
| 1000 | 0.3 % | 5.5 % | 0.003° |

All higher-order images live in the gap column (a few pixels at 1000 px).
Film stills use a narrower lens and a textured sky. Check without reference
pictures: a small bright sphere exactly behind the hole must render as a ring
of the tabulated diameter.

**The Gaussian filter is camera-side.** DNGR §3.3: the beam starts at the
camera as a 2-pixel circle in the image plane and is carried to the sky by
the deviation equation; a star inside the resulting ellipse lights the pixel
with the Gaussian weight of its position in the beam. Riazuelo: lens the
star direction first, then draw a fixed-angular-size blob as seen by the
observer. Both are blur *after* lensing. Blur before lensing = extended
source: stars near the Einstein ring become arcs, a star behind the hole a
ring. With the camera-side PSF a star is always ≈ 2 px; only brightness
follows magnification; arcs arise from overlapping neighbouring footprints.
For this code: 2-px-radius overlapping tubes, map each star back into
image coordinates via the footprint's local linear map, Gaussian weight of
the image-plane offset, normalise weights per star to 1 (flux
conservation, DNGR's 2 % flicker figure). Brightness stays the solid-angle
ratio.

**Why the magnification cap exists.** Pixel brightness from a star is flux
× magnification = image solid angle / footprint solid angle = 1/|det J| of
the lens map. The code measures the footprint from four traced corners. On
a critical curve (Einstein ring, photon rings) det J = 0 and the map
folds: locally u ≈ x², v ≈ y. A pixel straddling the curve maps to a strip
covered twice; its corners land pairwise on the same sky points, the quad
collapses to a line, the area → 0, and the two triangles wind oppositely.
Consequences: spurious magnification blow-up, missed stars (dark gaps in
the ring), sign cancellation in polygon gathers. The folded-quad test,
subdivision to single-sheet pieces, flat-gather fallback, hidden-fold
signed-area check and `star_magnification_cap` are all patches for this
one cause; the cap sets the photon-ring peak brightness by fiat.

The physical divergence is integrable: a point source at caustic offset u
has two images of magnification ∝ u^−1/2 each, so expected flux per pixel
from a random star field is finite and only a star exactly on the caustic
is formally infinite. The collapsed quad reproduces none of this.

**Jacobi-field footprint.** DNGR integrates the geodesic deviation
equation (App. A.2, Pineault–Roeder) for two deviation vectors spanning the
pixel; at the sky they are the Jacobian's columns, the ellipse semi-axes
its singular values, magnification 1/|det J|. Across a critical curve det J
crosses zero linearly, one axis shrinks smoothly to zero, nothing folds or
self-intersects, and sign(det J) is the image parity. The remaining
divergence is the physical one, bounded by the pixel spacing and spread by
the filter; a cap becomes a safety valve. Cost: 16 extra ODE components
plus the Riemann tensor (second metric derivatives: expensive with finite
differences, cheap with dual-number AD). Cheap intermediate: shrink the
four-corner stencil to a small fraction of a pixel around the centre so it
measures J at the centre; needs tight tolerance. Filter and footprint are
independent; do the filter first.

**DNGR equation (A.23) decoded.** With dots = d/dζ along the reference ray:

```
ü = −Ψ (g cos ψ + h sin ψ)     v̈ = −Ψ (g sin ψ − h cos ψ)
g̈ = −Ψ (u cos ψ + v sin ψ)     ḧ = −Ψ (u sin ψ − v cos ψ)
χ̇ = M,   Ψ = |Ψ₀*|,   ψ = arg Ψ₀* − 2χ                    (A.23, A.24)
start: u̇ = 1, rest 0 (A.22); end: δ± = δ′(√(u̇²+v̇²) ± √(ġ²+ḣ²)),
μ = χ + ½ arg[(u̇+iv̇)(ġ+iḣ)]                              (A.25, A.27)
```

Geodesic deviation in the transverse plane, complex form: ξ = u+iv (size),
η = g+ih (shear); bundle edge Y = (ξe^{iσ} + ηe^{−iσ})e^{iχ}, an ellipse
with semi-axes |ξ| ± |η|. The four real equations are ξ̈ = −Ψ₀* η̄,
η̈ = −Ψ₀* ξ̄ with Ψ₀* = C(k, m, k, m) the Weyl scalar on the ray (closed form
in Kerr, A.19–A.21); vacuum ⇒ no Ricci focusing term. χ tracks the twist of
the coordinate transverse basis relative to parallel transport (χ̇ = M,
A.18); −2χ because Ψ₀ has spin weight 2. Same Jacobi equation as the
Hamiltonian variational form: DNGR = 9 transverse components + closed-form
Weyl (Kerr-specific); variational = 16 components, metric-agnostic via AD.
Sky Jacobian singular values = δ±/2, rotation = μ, magnification
δ_cs²/(δ₊δ₋). Checks: principal-null-direction bundle stays circular
(Ψ₀ = 0); far-field magnification (u²+2)/(u√(u²+4)).

**Dictionary DNGR (u,v,g,h,χ) ↔ Jacobi matrix D.** Perlick: parallel
screen basis (E₁,E₂), Y = D θ. DNGR: twisting coordinate basis (a,b), χ =
angle to the parallel basis (χ̇ = M). In the parallel basis Y = ξθ + ηθ̄,
so D is Z ↦ ξZ + ηZ̄:

```
ξ = u+iv, η = g+ih;  D = [[u+g, h−v],[v+h, u−g]]
ξ = ½[(D₁₁+D₂₂) + i(D₂₁−D₁₂)],  η = ½[(D₁₁−D₂₂) + i(D₁₂+D₂₁)]
det D = |ξ|²−|η|²;  singular values |ξ|±|η|;  parity sign(|ξ|−|η|)
major axis: ½ arg(ξη) in (E₁,E₂), + χ in (a,b)
Ÿ = Φ₀₀Y − ψ̄₀Ȳ  ⇔  ξ̈ = Φ₀₀ξ − ψ̄₀η̄, η̈ = Φ₀₀η − ψ̄₀ξ̄
DNGR vacuum: ξ̈ = −Ψe^{iψ}η̄, η̈ = −Ψe^{iψ}ξ̄, ψ = arg Ψ₀* − 2χ  (spin-weight-2 rotation)
D(0)=0, Ḋ(0)=1 ⇔ ξ=η=0, ξ̇=1, η̇=0;  elliptical start Ḋ(0)=diag(1,e) ⇔ ξ̇=(1+e)/2, η̇=(1−e)/2
sky: semi-axes ∝ |ξ̇|±|η̇| = singular values of Ḋ; magnification = 1/|det Ḋ|
```

Traps: Perlick's χ (eq. 19) is the ellipse angle, DNGR's χ the basis twist;
Perlick's D̈ = D·R is the transpose of the vector form. DNGR vs
Pineault–Roeder: PR write (x,y,p,q) for (u,v,g,h); DNGR flipped the null
tetrad vector l from 2^−½(1,−1,0,0) to 2^−½(1,+1,0,0) (ingoing instead of
outgoing rays), turning 1/(p_t̂+p_r̂) into 1/(p_t̂−p_r̂) in Ψ₀* and M; the
evolution equations are identical. For the KS chart, use Perlick's
basis-free form with ψ₀ from the Riemann tensor, or the variational form.

**Route 2: Jacobi equation = linearised Hamiltonian flow.** With
ẋ^μ = g^{μν}p_ν, ṗ_μ = −½∂_μg^{αβ}p_αp_β and a family x(λ,s), p(λ,s),
Y = ∂x/∂s, Π = ∂p/∂s:

```
Ẏ^μ = ∂_ν g^{μα} p_α Y^ν + g^{μν} Π_ν
Π̇_μ = −½ ∂_μ∂_ν g^{αβ} p_α p_β Y^ν − ∂_μ g^{αβ} p_α Π_β      (δy' = (∂F/∂y) δy)
```

Y is Perlick's deviation vector; Π = g·∇_kY + Christoffel term. One more
derivative gives ∇_k∇_kY = R(k,Y)k (textbook derivation: Carroll 3.10, Wald
3.3, MTW 11.3). Numerically: evaluate the RHS on y + ε·δy with a dual ε
(∂g must be exact ⇒ first-order dual for ∂g, hyper-dual for the variational
step). Checks: symplectic form δx₁·δp₂ − δx₂·δp₁ constant along the ray;
δH = 0. References: Grasso–Korzyński–Serbenta arXiv:1811.10284 §II.A (8×8
resolvent W, eqs. 10–21; Jacobi map = 2×2 block, eq. 54; reduced eqs.
56–58); Uzun arXiv:1811.10917 §3 (ż = L_H z, L_H = [[0,1],[R,0]], eqs.
53–62, symplectic ABCD transfer matrix); Hairer–Nørsett–Wanner, Solving
ODEs I, §I.14 (variational equation); Fike & Alonso 2011 (hyper-dual
numbers).

Concepts: point spread function; pullback of a pixel footprint through the
lens map; point-source vs extended-source magnification; partition of unity
for flux-conserving resampling; critical curves, caustics and the fold
catastrophe; integrable singularities; geodesic deviation / Jacobi fields;
image parity; Einstein radius vs lens distance; lensing ring vs photon ring
(Gralla–Holz–Wald); stereographic vs gnomonic projection; Weyl scalars and the Newman–Penrose tetrad; spin-weighted quantities; Sachs optical scalars; Fermi–Walker transport of a screen basis; Pineault–Roeder complex bundle representation.

## 9. Interior camera and the Kerr–Schild chart

See `physics-review-2026-09.md` §11 for the measurements. Summary: rays
crossing the future horizon integrate cleanly (≤ 105 steps to r = 60,
|ΔH|/|p|² ≤ 1e-6 at ε = 1e-5); rays asymptoting to the past horizon stall at
H_MIN and, because the integrator accepts steps at H_MIN regardless of
error, tunnel and explode. Needed: a chart-independent stall stop (N
consecutive H_MIN steps), "error above ε at H_MIN" as a stop reason, no
inner-horizon stop for a camera between r₋ and r₊. Boyer–Lindquist cannot
host an interior camera (imaginary ZAMO lapse, Δ = 0 at both horizons); a
BL-outside/KS-inside atlas is a speed optimisation. If an atlas is added:
transitions between accepted steps only, closed-form azimuth twist and time
shift, covariant momentum via the Jacobian transpose, E/L/Q checked across
the switch, hysteresis wider than H_MAX.

Concepts: atlas and transition functions; cotangent lift of a coordinate
change; ingoing vs outgoing Eddington–Finkelstein / Kerr–Schild and which
horizon each covers; Penrose diagram of Kerr; coordinate vs curvature
singularity; maximal analytic extension; Cauchy horizon.

## 10. GPU notes

Only the per-step inner loop moves (metric, inverse, derivatives, RHS,
RKF45 step and controller, stop conditions, hit tests, redshift and
temperature lookups): ≈ 1500–2000 of 15 000 lines. Host keeps config,
camera/tetrad, star gather and octree, temperature table, supersampling,
compositing, output, tests. Polymorphism survives as generics
(`fn trace<G: Geometry>`, one monomorphised kernel per geometry); scene
objects become a tagged enum with a uniform `match`; scalar type generic
(f64 or dual). Shared crate `no_std`, no heap / `dyn` / recursion /
formatted panics; nalgebra and logging host-side behind `cfg`. Routes:
rust-gpu (SPIR-V, vendor neutral, f64 via Float64 capability, host via
`ash` or `wgpu` passthrough); Rust-CUDA (PTX, NVIDIA); CubeCL (`#[cube]`
DSL, f64 on CUDA/HIP); or foreign-language kernels. Ruled out by f64: WGSL,
browsers, Apple GPUs; consumer NVIDIA f64 is 1/32–1/64 rate. Do fused loop,
AD and the first-order BL solver first.

Concepts: loop fusion, arithmetic intensity, roofline; static vs dynamic
polymorphism, monomorphisation; uniform vs divergent control flow; enum
dispatch; SPIR-V and PTX; rustc codegen backends; Float64 capability;
double-precision throughput per GPU generation; stream compaction.

## 11. Reading

- Coulon, Matsumoto, Segerman, Trettel, *Ray-marching Thurston geometries*, arXiv:2010.15801.
- Hart, Hawksley, Matsumoto, Segerman, *Non-euclidean virtual reality I/II*, arXiv:1702.04004, arXiv:1702.04862.
- James, von Tunzelmann, Franklin, Thorne, *Visualizing Interstellar's Wormhole*, arXiv:1502.03809.
- Gibbons, Werner, *Applications of the Gauss–Bonnet theorem to gravitational lensing*, arXiv:0807.0854.
- Müller, Weiskopf, *Distortion of the stellar sky by a Schwarzschild black hole*; Müller's GeoViS / Motion4D as the closest general-metric raytracers.
- Grave, Buser, *Visiting the Gödel universe*, IEEE TVCG 2008.
- Ray bundles, explained better than DNGR's appendix: Perlick, *Gravitational lensing from a spacetime perspective*, arXiv:1010.3416 (§2.3: D̈ = D·R, D(0)=0, Ḋ(0)=1, tidal matrix from Φ₀₀ and ψ₀, polar decomposition into D₊, D₋, rotation; §4.3 eqs. 96–100: closed-form D± in spherically symmetric static spacetimes, the analytic test for a Schwarzschild bundle).
- Seitz, Schneider, Ehlers, *Light propagation in arbitrary spacetimes and the gravitational lens approximation*, arXiv:astro-ph/9403056 (§2 same equation, |det D| = δA/δΩ; §3 behaviour at conjugate points: det D ∝ ε at a fold, ∝ ε² at a focus).
- Fleury, Larena, Uzan, *Weak gravitational lensing of finite beams*, arXiv:1706.09383 (validity conditions of the infinitesimal-beam description).
- Grasso, Korzyński, Serbenta, *Geometric optics in general relativity using bilocal operators*, arXiv:1811.10284; Serbenta thesis arXiv:2305.18843 (Jacobi map as an 8×8 resolvent of position and momentum deviations = the variational-equation form).
- Pineault, Roeder, *Applications of geometrical optics to the Kerr metric* I and II, ApJ 212, 541 and 213, 548 (1977). Original ξ, η formulation.
- Igehy, *Tracing ray differentials*, SIGGRAPH 1999; Heckbert, *Fundamentals of texture mapping and image warping*, 1989 (ray differentials and EWA filtering, the flat-space graphics counterparts).
- Hart, *Sphere tracing*, The Visual Computer 1996.
- James, von Tunzelmann, Franklin, Thorne, *Gravitational lensing by spinning black holes in astrophysics, and in the movie Interstellar*, arXiv:1502.03808 (§3.3 star filter, App. A.2 ray bundles).
- Riazuelo, *Seeing relativity I*, arXiv:1511.06025.
- Gralla, Holz, Wald, *Black hole shadows, photon rings, and lensing rings*, arXiv:1906.00873.
- Ellis, *Ether flow through a drainhole*, J. Math. Phys. 14, 104 (1973).
- Leonhardt, Philbin, *Geometry and Light: The Science of Invisibility*.

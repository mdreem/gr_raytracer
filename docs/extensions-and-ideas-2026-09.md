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

## 8. Reading

- Coulon, Matsumoto, Segerman, Trettel, *Ray-marching Thurston geometries*, arXiv:2010.15801.
- Hart, Hawksley, Matsumoto, Segerman, *Non-euclidean virtual reality I/II*, arXiv:1702.04004, arXiv:1702.04862.
- James, von Tunzelmann, Franklin, Thorne, *Visualizing Interstellar's Wormhole*, arXiv:1502.03809.
- Gibbons, Werner, *Applications of the Gauss–Bonnet theorem to gravitational lensing*, arXiv:0807.0854.
- Müller, Weiskopf, *Distortion of the stellar sky by a Schwarzschild black hole*; Müller's GeoViS / Motion4D as the closest general-metric raytracers.
- Grave, Buser, *Visiting the Gödel universe*, IEEE TVCG 2008.
- Hart, *Sphere tracing*, The Visual Computer 1996.
- Ellis, *Ether flow through a drainhole*, J. Math. Phys. 14, 104 (1973).
- Leonhardt, Philbin, *Geometry and Light: The Science of Invisibility*.

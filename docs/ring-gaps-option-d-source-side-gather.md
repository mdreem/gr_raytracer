# Option D: Source-side (inverse) star gather

Guarantee that every image of every star is found by inverting the lensing map,
for each star solve for the screen positions whose ray lands on it, instead of
casting screen tubes and hoping their footprints tile the sky and catch the
star.

## Why invert the map

The current pipeline is a forward map: screen -> sky. It gathers a star only if
some tube's sky footprint happens to cover the star's direction `S`. Near the
critical curve the footprints mistile and slivers open, so images are missed
(the ring gaps). No forward-side patch fully removes that reliance, because the
guarantee it needs ("every point of the sky is covered by exactly one tube, even
where the map folds") is exactly what breaks at a caustic.

The inverse map removes the reliance. The lensing map `f: screen -> sky` is
smooth and, away from caustics, locally invertible. A star at sky direction `S`
has a discrete set of preimages on the screen, its images:

- a primary image (weakly deflected, near the unlensed position),
- a secondary image on the far side of the shadow,
- an infinite series of higher-order images spiralling into the photon ring,
  each exponentially fainter and closer to the shadow edge.

If we find those preimages directly, each star contributes its images with no
dependence on tube tiling, so the ring is complete by construction, down to
wherever we choose to truncate the higher-order series.

## The math

Let a screen point be `u = (row, col)`. The renderer already maps `u` to a sky
direction by integrating the geodesic: `f(u) = escape_direction(ray(u))` (this
is `Scene::color_of_ray` / the escape info used at `raytracer.rs:396-419`). For
a star direction `S` we want all `u` with

```
f(u) = S    (as unit vectors, i.e. f(u) x S = 0)
```

This is a 2D root find (two angular residual components, two screen unknowns).
Newton's method needs the 2x2 Jacobian `J = df/du`. Two ways to get it:

- **Finite differences.** Trace `f` at `u`, `u + (h,0)`, `u + (0,h)`; three
  geodesics per Newton step. Simple, robust, the recommended first cut.
- **Ray-adjoint / geodesic deviation.** Integrate the Jacobi equation alongside
  the geodesic to get `df/du` analytically in one pass. Faster per step, much
  more code, and needs the deviation equations for each metric. A later
  optimization, not the starting point.

The image magnification is then `mu = |det J|^-1` up to the screen/sky solid-angle
Jacobian, which is the same quantity the forward path computes as
`ratio = Omega_pixel / Omega_traced`. So the flux splat per image is
`flux_star * mu`, distributed onto the pixels the image covers (a small PSF
splat, see Regularization).

## Finding all the images

Seeding Newton is the crux, because the map is multi-valued. The key
simplification for Schwarzschild: by spherical symmetry, **every image of a star
lies on the single line through the shadow centre and the star's undeflected
screen position.** So seeding is never a 2D search, it is placing a few points
along one radial line at radii the lens equation predicts.

Set up three scene constants once (independent of the star): `O`, the shadow
centre on screen (project the camera->hole direction); `theta_E`, the
Einstein-ring radius in pixels (the radius of the rendered ring); and
`theta_ph`, the shadow-edge / photon-ring radius in pixels. Then per star,
project it to its undeflected screen point `P0`, let `n_hat` be the unit
direction from `O` to `P0`, and `beta = |P0 - O|`. All seeds are `O + t * n_hat`.

1. **Primary and secondary seeds.** The point-mass lens equation
   `beta = theta - theta_E^2 / theta` has two roots,
   `theta_pm = 0.5 * (beta +/- sqrt(beta^2 + 4 theta_E^2))`: the primary at
   `O + theta_+ * n_hat` (outside the ring, on the star's side) and the secondary
   at `O + theta_- * n_hat` (`theta_- < 0`, inside the ring, opposite side), with
   `|theta_+| |theta_-| = theta_E^2`. At perfect alignment (`beta = 0`) both
   collapse to `+/- theta_E`, the full ring. Newton polishes each under the true
   metric. Two seeds, two images, cheap.
2. **Higher-order seeds.** The photon ring stacks images indexed by winding
   number `n = 1, 2, ...`, each at a screen radius offset from the shadow edge
   that shrinks like `exp(-pi n)` for Schwarzschild (the Lyapunov
   exponent of the photon orbit). Seed each `n` just outside the shadow at the
   predicted radius and let Newton converge. Truncate at `n_max` where the
   contribution drops below a flux epsilon (usually `n_max` = 2-3 suffices; the
   series is exponentially convergent).
3. **Continuation.** Neighbouring stars have neighbouring images, so once one
   star's images are found, use them as seeds for the next star (sort stars by
   sky position and walk). This turns most Newton solves into one or two steps.

Robustness details:

- **Caustic crossings.** When a star is near the critical curve, `det J -> 0`,
  Newton is ill-conditioned and the primary/secondary images merge into the
  ring. Detect small `det J` and fall back to a bracketed 1D solve along the
  radial line, or hand that star to the forward path for the merged region.
- **Duplicate / missing images.** Deduplicate converged roots (two seeds landing
  on the same `u`), and verify each root actually satisfies `f(u) = S` to reject
  spurious convergence.
- **Shadow test.** Reject seeds whose ray is captured (no escape); those
  directions have no image.

## Regularization (the caustic singularity)

An idealized point star exactly on a caustic has `det J = 0`, so `mu = infinity`.
This is the same singularity the resolution experiment exposed (peak intensity
103 -> 406 -> 5408 with no convergence). The inverse map does not remove it; it
relocates it into `mu`. The physically correct cure is a **finite source**: give
each star an angular size (or a fixed screen-space PSF) and splat its flux with
that kernel, which integrates the divergent `mu` over a small area to a finite
value. This is natural on the source side: convolve the image with the star's
PSF as you splat. It is the main reason D is attractive beyond gap-filling, it
makes the ring brightness a well-defined, convergent number.

## Cost

Per star: a handful of Newton solves (2 low-order images + a few higher-order),
each a few geodesic traces. Cost scales as `O(n_stars * images_per_star *
newton_steps * traces_per_step)`. Two accelerations make it tractable:

- Only stars whose sky direction lies within a modest angle of the shadow get
  extra (secondary, higher-order) images; distant field stars have one nearly
  undeflected image and can stay on the cheap forward path. So the expensive
  inverse solve runs only for the small near-line-of-sight subset.
- Continuation seeding (above) collapses most solves to one or two steps.

Compared with C, D does not blow up at the caustic the way tube refinement does,
because it never subdivides a diverging region; it places a bounded number of
images and regularizes their brightness with the PSF.

## Integration: a hybrid pipeline

D does not have to replace the forward path wholesale. The natural design:

- **Forward path** renders the background, the shadow, opaque objects, and the
  bulk star field (one image per distant star) as today.
- **Inverse path** runs only for stars near the line of sight to the hole, adds
  their secondary and higher-order images, and splats them with a PSF.
- Composite the two by summing flux (both are emission-only, no opacity
  interaction), matching how `add_star_layer` already composites starlight.

This keeps the common case on the fast, well-tested path and pays the inverse
cost only where it is needed.

## Validation plan

- Single-star aligned: a fully closed, symmetric Einstein ring with no gaps at
  any zoom, and a convergent brightness once the PSF is enabled.
- Single-star offset: exactly two low-order images at the analytic point-mass
  positions, plus a resolvable first higher-order pair, with correct flux ratio.
- Higher-order spacing: measured radii of the `n = 1, 2` images match the
  `exp(-pi n)` Schwarzschild law.
- Flux: total over all images of an off-caustic star matches the analytic
  magnification; on-caustic flux is finite and converges as the PSF narrows to
  the pixel scale.
- Gallery frames: identical bulk field; the near-shadow annulus gains complete
  arcs and faint higher-order echoes.

## Files / scope

Substantial, a new module rather than a patch:

- New `src/rendering/star_inverse.rs`: the Newton solver, seed generators,
  continuation, PSF splat, dedup, caustic fallback.
- `src/rendering/raytracer.rs`: hook the inverse pass after the forward pass and
  composite; reuse `color_of_ray` for `f(u)` and the existing solid-angle code
  for `mu`.
- `src/rendering/scene.rs` / `configuration.rs`: config for `n_max`, flux
  epsilon, PSF width, and the near-line-of-sight angular cutoff.
- Optional later: Jacobi-equation Jacobian per metric for the analytic `J`.

## Limitations

- **Complexity and per-metric work.** Robust multi-image root finding through a
  numerically integrated map is real work: seeding, continuation, caustic
  handling, dedup. The analytic Jacobian, if pursued, must be derived per metric
  (Schwarzschild, Kerr, Ellis, ...).
- **Kerr and non-symmetric caustics.** The clean radial seeding above is
  Schwarzschild's spherical symmetry. Kerr's caustics are non-circular and the
  higher-order structure is richer; seeding and continuation get harder, though
  the approach still holds.
- **Not needed for the current gallery look.** As the experiment shows, the
  visible payoff on the committed frames is small; D earns its cost only if you
  want physically faithful, zoomable rings and higher-order photon-ring
  substructure.

## Verdict

The correct answer if the goal is guaranteed image completeness and convergent
ring brightness, including higher-order photon-ring images and zoom studies. It
sidesteps both the tiling reliance and the caustic-cost blowup, and the PSF
splat gives the singularity a principled cure. But it is a project: a new
solver module with nontrivial numerics, best deployed as a hybrid alongside the
existing forward path rather than as a replacement.

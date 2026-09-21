# Plan 13: Point-spread ("real telescope") star rendering

## Status: not started (design, 2026-09-21)

## Motivation

Catalogue stars are currently gathered as **delta functions binned into tubes**:
`compute_star_collection_data` asks the octree (`gather_triangle`) for the stars
whose escape direction falls inside a tube's two sky triangles, sums each star's
XYZ flux, and multiplies the tube total by the lensing magnification
`ratio = pixel_solid_angle / traced_sky_solid_angle`. A star lands wholly in
whichever tube contains it.

Total flux is conserved, so the physics is not *wrong*, but binning deltas into
pixels has three practical costs:

1. **Resolution-dependent brightness (the exposure footgun).** A pixel-tube's
   solid angle grows as resolution drops, so a coarse tube gathers more stars
   and reads brighter. The same scene needs a different exposure at every
   resolution (≈11x between 480x270 and 1600x900). Total image flux is the same;
   it is just packed into fewer pixels.
2. **Temporal aliasing.** As a camera moves by a sub-pixel amount, a star hops
   from one tube to the next and its pixel flickers. Isolated stars twinkle and
   crawl in fly-throughs — an artifact, not physics.
3. **No point-spread function.** Stars are single hard pixels: no diffraction
   glow, poor anti-aliasing, and near the photon ring the "star trail" is a
   smear of the tube footprint rather than a crisp lensed arc. (The bloom in the
   external tone-map pass is a crude, non-physical stand-in.)

A real telescope images a point source not as a pixel but as a **point-spread
function** (PSF): a diffraction Airy disc set by the aperture, convolved with
atmospheric seeing, well-approximated by a 2D Gaussian of fixed *angular* width.
Adopting that model fixes all three problems at once, and — because a star then
occupies a fixed solid angle carrying a fixed flux — it makes brightness
**resolution-independent by construction**, removing the exposure footgun as a
side effect rather than papering over it.

This extends `plan-12-star-catalog.md` step 5 ("optional PSF: splat as a small
Gaussian"), which named the idea but did not design it.

## Physical model

- A point source of flux `F` is imaged as `F * PSF(x)` on the detector, where the
  PSF is normalised (`∫ PSF = 1`), so flux is conserved. Use a 2D Gaussian of
  angular FWHM `theta_psf` as the default (Airy is a later refinement).
- The PSF lives in the **observer image plane** — the detector — which is a
  uniform grid. Lensing has already done its work *before* the detector: it set
  each star image's apparent position and its magnification `mu`. The telescope
  then blurs that image by a fixed angular PSF. So in screen space the PSF is a
  fixed-*angular*-width Gaussian (a fixed number of pixels for a given FOV and
  resolution), uniform across the frame; it does **not** itself vary with local
  magnification.
- **Multiple images come for free.** A star behind the hole has many apparent
  images (Einstein ring, higher-order images). Backward tracing already produces
  them: every tube whose escaped ray points at that star yields one image, each
  with its own `mu`. The PSF is applied independently at each apparent position,
  so ring arcs and Einstein crosses fall out of the same machinery.
- **Brightness per image is unchanged:** each image contributes `F * mu`,
  exactly as today. The PSF only *redistributes* that energy over a small fixed
  angular neighbourhood; it never changes the per-image total.

### Why a point source is magnified by `mu` (Liouville)

This is not an assumption; it follows from **Liouville's theorem for photons**:
specific intensity `I_nu / nu^3` is conserved along a geodesic, i.e. surface
brightness is invariant under lensing (up to redshift). A star is unresolved, so
its only observable is the integrated flux `F = integral of I dOmega`. The lens
map stretches the apparent solid angle by the Jacobian `mu = dOmega_image /
dOmega_source` (that *is* the definition of magnification, `mu = 1/|det A|`).
Since `I` is conserved,

    F_obs = I * dOmega_image = I * mu * dOmega_source = mu * F_source.

So a point source magnifies in **flux** by `mu`, precisely because surface
brightness is conserved and its apparent solid angle is stretched. The renderer's
`ratio = pixel_solid_angle / traced_sky_solid_angle` is exactly that Jacobian, so
`F * ratio` is the correct point-source photometry. Redshift adds the `g^4`
beaming from the `nu^3` factor, applied separately in `redshifted_star_emission`,
so `mu` and `g^4` are not double-counted.

### Pixel space vs sky space

Keep the two spaces distinct. **Lensing lives in the sky/source map** (the
screen<->sky Jacobian, `mu`, the folds that make multiple images). **The PSF
lives in the detector/pixel plane**: light arrives at the aperture from an
*apparent* direction, and the telescope blurs it there. So lensing acts first
(placing each image's apparent screen position and brightness `F * mu`), and the
PSF acts second, as a fixed, symmetric, uniform blur in pixel space. Because the
detector pixel angular scale is uniform, a fixed *angular* PSF is a fixed number
of *pixels* everywhere in the frame; it does not shrink or grow with local
magnification.

Why this removes the resolution footgun: today a star's energy is dumped into
one pixel, so packing N stars into a coarse pixel sums N deltas. With a PSF, each
star deposits `F * mu` spread over a fixed *solid angle*. Surface brightness
(flux per steradian) is then well-defined and resolution-invariant; a finer grid
splits the same PSF over more pixels but the star looks identical. Exposure
becomes a single value independent of resolution.

## Design

### Sub-pixel position of a gathered star

The octree tells us a star is inside a tube's sky triangle, not *where*. Recover
the sub-pixel screen position by interpolation: the tube maps a screen quad to a
sky quad (its four traced corner directions). Compute the star's barycentric /
bilinear coordinates within the sky triangle, then apply the same coordinates in
the tube's screen quad. That yields the star's apparent screen position to
sub-pixel accuracy — the location its PSF is centred on. (This is exact only to
the tube's linear approximation of the lensing map, which is already the
approximation the magnification `ratio` assumes, so no new error is introduced.)

### Recommended approach: sub-pixel splatting into a shared star buffer

1. Allocate a frame-wide **linear XYZ star accumulation buffer** (separate from
   the foreground), plus an optional per-pixel coverage/weight if needed.
2. For each tube, for each gathered star: compute its screen position (above) and
   its per-image brightness `F * mu`, and **splat a normalised Gaussian kernel**
   of fixed angular width into the shared buffer at that position, scaled by
   `F * mu`. A 5x5 (or radius = 2-3 sigma) kernel suffices for a few-pixel PSF.
   Splatting into a *shared* buffer (not per-tube) is essential: a star's PSF
   straddles tube boundaries, and per-tube accumulation would clip it.
3. After all tubes, composite the star buffer under the foreground, weighted by
   the foreground's surviving transmittance, as `add_star_layer` does now.

This keeps the existing octree gather and magnification untouched; it changes
only *how* a gathered star deposits its energy (a kernel instead of a single
tube add). It is temporally stable (sub-pixel motion moves the kernel smoothly),
resolution-stable (fixed angular kernel), and produces diffraction glow directly.

**Where the ring arcs come from.** Not from a stretched PSF. The Einstein-ring
arc is built from **many tubes near the ring each gathering the same source** at
slightly different apparent positions (the fold / high tangential magnification),
each splatted with a *symmetric* kernel; those overlapping blobs strung along the
ring fill the arc. The kernel stays symmetric in pixel space. At low resolution a
single highly-magnified tube under-represents a sub-pixel tangential stretch; if
that matters, make the kernel anisotropic from the local lensing Jacobian (the
same `mu` machinery gives the stretch direction). That is an optional refinement.

### The octree is unchanged

The octree indexes star **directions** (points) and answers "which tube owns this
image" — a point query. The PSF does not touch it, because *collection* and
*deposition* are decoupled: a star is gathered by its centre's tube, and its
kernel spilling into neighbouring pixels is handled at deposition via the shared
buffer, not by any tube re-gathering it. A star is never stored in multiple
boxes. The only query-side tweak is the half-open boundary convention below,
which lives in the `point_in_triangle` predicate, not in the tree.

### Alternatives considered

- **Supersample + convolve.** Render the star gather at Sx angular resolution,
  convolve with a fixed-angular Gaussian, downsample. Simple and correct, but
  pays a full Sx render cost. Splatting gets the same result far cheaper because
  we already know each star's sub-pixel position analytically.
- **Analytic forward projection.** Project each star's apparent screen position
  by inverting the ray map globally. Backward tracing does not give this
  directly; the per-tube interpolation above is the tractable equivalent.

## Collection and deduplication

A star gathered by several tubes is two very different situations:

- **Multiple images (a feature).** Near the ring the screen<->sky map folds, so
  one star legitimately appears at several apparent positions, each in a
  different tube with its own `mu`. Splat each — that is what draws the arcs.
  These are **not** a partition of a luminosity budget: each image gets the full
  `F * mu_i`, and the total `sum_i F * mu_i` can far exceed `F`. Lensing
  amplifies; it does not divide. Never normalise luminosity across images.
- **Tiling seams (a bug).** The same apparent image is gathered twice at a shared
  boundary: the tube's diagonal (both triangles contain a star on it, today's
  `hits = 2`), or the edge between adjacent un-folded tubes. Measure-zero
  geometrically, but a finite PSF turns it into a bright seam.

Fix the seams with **half-open tiling**: make each boundary owned by exactly one
region (one triangle owns the shared diagonal with `>= 0`, the other excludes it
with `> 0`; pick a consistent owner for tube edges). Then each apparent image is
collected once, and you splat once per gather. Splatting the PSF within one image
distributes `F * mu_i` over neighbouring pixels with a normalised kernel
(weights sum to 1) — that partition conserves the image's flux. So: partition
*within* an image (the PSF), never *across* images (magnification).

## PSF vs bloom

They both spread a bright point into a glow but are opposite things and both are
wanted:

- **PSF** is physical image formation (diffraction + seeing), applied **per
  star** at its true sub-pixel position, **before** tone mapping, in linear
  light. It is **flux-conserving** (kernel integral 1), tight (a few pixels), set
  by the optics, and it fixes aliasing, temporal flicker, and resolution-
  dependent brightness. It is the signal.
- **Bloom** is a post-process over the **whole image**, at/after tone mapping,
  modelling veiling glare in the eye/lens/sensor. It **adds** energy (blurred
  copies screened back), is wide and soft (tens of pixels), and is purely
  aesthetic. It is the flare.

The current pipeline has only bloom, which quietly doubles as a *fake* PSF
(softening single-pixel stars). Bloom cannot fix flicker, resolution-dependent
brightness, or photometry; only a real PSF can. Once the PSF is in, dial bloom's
threshold/radius down so it only glows the genuinely bright stars.

## Flux normalisation (resolution-stable exposure)

Deposit surface brightness, not per-pixel flux: when compositing the star buffer,
divide by the pixel solid angle (or fold a reference solid angle into
`star_flux_scale`) so the stored quantity is flux per steradian. Then one
exposure works at any resolution, and the current "retune exposure per
resolution" step disappears. This is the physically meaningful normalisation the
PSF model enables; without a PSF it is ill-defined for point sources.

## Colour and redshift

Unchanged. Each star's XYZ (with the per-tube frequency-shift retint from
`redshifted_star_emission`, `T_obs = g * T`) is the quantity splatted. The PSF is
achromatic by default; a wavelength-dependent Airy radius is a possible later
refinement.

## Parameters

- `theta_psf`: PSF angular FWHM. Physically the diffraction limit `~1.22 * lambda
  / D` for an aperture `D`, or a chosen "seeing" value; practically a small
  multiple of a pixel (e.g. 1.5-2 px at the reference resolution) tuned for look.
- Kernel radius in sigma (default ~3) and shape (Gaussian now, Airy later).
- A single resolution-independent exposure, once normalisation is in.

## Anti-aliasing and movies

Sub-pixel splatting is already the main AA win. For fly-throughs, jitter the
tube-corner sample offsets per frame (the camera path plus existing supersampling
hooks) so residual sampling noise averages rather than crawls. The fixed-angular
PSF guarantees an isolated star's size and brightness are constant frame to
frame, which is what makes smooth motion possible.

## Validation

- **Flux conservation:** sum of the star buffer over all pixels equals
  `sum(F_i * mu_i)` over gathered images, independent of `theta_psf` and
  resolution.
- **Analytic magnification:** render one star at a known angular offset, sum its
  image brightnesses, and compare to the analytic point-source magnification
  (`mu(u) = (u^2 + 2) / (u * sqrt(u^2 + 4))` for a point mass, `u` in Einstein
  radii; the Schwarzschild analogue exactly). Also check `mu -> 1` far from the
  hole (flat-space flux recovery).
- **Resolution invariance:** the same scene at two resolutions yields matching
  surface brightness (up to sampling), with a *single* exposure.
- **Isolated-star invariance:** one star's integrated brightness is independent
  of its sub-pixel position and of resolution.
- **Temporal smoothness:** a star drifting sub-pixel across frames brightens and
  moves continuously, with no step when it crosses a pixel boundary.
- **Lensing:** a star swept behind the hole forms a continuous Einstein-ring arc,
  not a dotted trail.

## Rollout and risks

Incremental:
1. Sub-pixel position from tube interpolation + fixed screen-space Gaussian splat
   into a shared buffer. (Fixes flicker and hard pixels.)
2. Solid-angle normalisation. (Fixes the exposure footgun.)
3. Refine PSF (Airy, chromatic), and per-frame jitter for movies.

Risks:
- **Cost near the ring:** many gathered images per tube means many splats. Bound
  the kernel radius and consider a brightness threshold below which a star splats
  as a single pixel.
- **Double counting:** the magnification `mu` already carries the point-source
  photometry; the PSF must only redistribute `F * mu`, never rescale it. Keep the
  kernel normalised to 1.
- **Boundary correctness:** splat into the shared frame buffer, never per-tube,
  or PSFs clip at tube edges.
- **Interaction with the diagonal double-hit:** a star on a tube's shared
  diagonal is currently gathered by both triangles (weight 2). With explicit
  positions this should be deduplicated (gather once per apparent image), which
  the splatting rewrite is a natural place to fix.

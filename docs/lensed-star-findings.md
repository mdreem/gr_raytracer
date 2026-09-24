# Lensed point-star rendering: findings and open problems

A working record of the problems found while making the Gaia point-star gather
render the lensed ring correctly. Each entry states the problem, its root cause,
what was measured, the current status (fixed / stopgap / open), and the proper
fix. Companion design docs: `ring-gaps-option-c-fold-aware-refinement.md`,
`ring-gaps-option-d-source-side-gather.md`, `coordinate-atlas.md`,
`plan-13-star-psf.md`.

## The pipeline in one paragraph

Stars are point sources indexed in an octree by their celestial direction. For
each screen pixel the renderer casts a four-corner "tube", integrates the corner
rays as geodesics, and forms a spherical quad from their escape directions. Every
star whose direction falls inside that quad is gathered, and its flux is scaled
by the lensing magnification `ratio = A / B`, where `A` is the pixel's screen
solid angle and `B` is the tube's traced sky footprint solid angle. The star
layer is composited under the foreground. This is a forward map (screen -> sky)
with per-tube point gathering.

## Two critical curves, and why folds happen "far out"

The fold map (render with folded tubes painted) shows two concentric critical
curves, not one:

- an **inner ring** hugging the shadow: the **photon ring** (strongly wound
  rays), and
- a **much larger outer ring**: the **primary Einstein ring**.

Folds and magnification blow-ups happen on both, so a lensed-star artifact can
appear well away from the black disc (on the outer Einstein ring), which looks
"almost in the flat part" but is a genuine strong-lensing critical curve. Any
reasoning that assumes lensing artifacts cluster at the shadow edge is wrong.

## Problem 1: flat-quad membership leaves gaps in the ring

**What.** The gather tests whether a star's direction is inside the tube's
*flat* spherical quad (two triangles through the four traced corners). Near the
ring this mistiles and the lensed ring comes out broken/gappy.

**Root cause.** The tube's edges map to curved *arcs* on the sky, not straight
lines. The flat quad is a linear approximation whose error grows with the map's
curvature. Adjacent tubes share corners (hence chords), so the error along shared
edges is misassignment, not holes; the true holes come from the strong curvature
and folding at the critical curve.

**Status: addressed (curved membership, opt-in).** Refine each edge into a
polyline that follows the arc (`refine_edge`/`tube_boundary`), gather against the
fan-triangulated curved polygon (`gather_polygon`). Enabled by
`--curved-star-membership` or `[star_catalog].curved_star_membership`. On a
sparse real field the visible change is tiny (RMSE 2.5e-4, confined to a thin
outline at the critical curve); it matters for full Einstein rings and zoom.

**Proper fix (not done).** The edge refinement makes each tube's own boundary
accurate, but subdivision-created children still mistile at folds (see Problem
2). True completeness needs either sliver-free curved tiling (winding-number
containment over the true arcs, or dilated child polygons with partition-of-unity
weights) or the inverse/source-side gather in `ring-gaps-option-d`.

## Problem 2: folds at the critical curve, and the star-drop bug

**What.** A tube the critical curve crosses has a *self-intersecting* sky
footprint (the map folds). Its two tiling triangles wind opposite ways. The
naive fan gather on a folded polygon is meaningless, and a star there can be
dropped.

**The specific bug found.** With curved membership, a folded tube subdivides, but
a star can fall into a **sliver between the child polygons**, inside the parent's
flat quad but inside no child. The child that should own it is neither
corner-flagged folded nor self-intersecting, so its gather silently returns zero.
Concretely: a bright star on the primary Einstein ring rendered as `0.232`
(flat) vs `0.0` (curved). This is a *collection* failure, "the star was not
properly addressed when found via the tube."

**Status: stopgapped.** `gather_curved` now: (a) flat-falls-back when a fold
cannot be split (depth cap); (b) detects hidden double-folds the corner parity
test misses via signed-vs-absolute polygon area; and (c) for a folded tube that
does split, takes `max(subtube_sum, flat_gather)` so it never gathers less than
the flat path. Result: zero bright drops, curved gathers `>= flat` everywhere
(ratio 1.0001).

**Why it is only a stopgap.** `max(subtube, flat)` is a floor, not a
reconstruction. It papers over dropped slivers with the flat value rather than
tiling correctly; it can keep a flat over-count; and it cannot both fill a gap
and capture a dropped star in the same tube. The real fix is Problem 1's proper
fix (sliver-free tiling / inverse gather).

**Fold detection subtlety.** The corner-orientation parity test only sees fold
*parity*; an even number of critical-curve crossings (a double fold) cancels and
is missed. Robust detection needs interior sampling (or the recursive split to
expose it). The signed-area check catches some, not all, hidden folds.

## Problem 3: the coordinate polar-axis singularity

**What.** A straight bad meridian runs through the shadow center; for a full
Einstein ring it punches a gap where that meridian crosses the ring.

**Root cause.** The integrator works in spherical / Boyer-Lindquist coordinates.
The escape direction is built from `get_cartesian_vector`, full of `1/sin(theta)`
behaviour, and the geodesic equations carry `dphi/dlambda ~ 1/sin^2(theta)`. Rays
grazing the polar axis (`theta -> 0/pi`) return corrupted directions or NaN
(`StopReason::CoordinateIsNan`). The polar axis projects to that meridian.

**Confirmed** by rotating the coordinate frame: moving the pole moves the gap.

**Status: stopgapped (pole rotation).** `--pole-rotation-deg` (default 0 =
identity, byte-identical) rotates the spherical frame so the pole sits off the
field of view (`spherical_coordinates_helper::set_frame_rotation_deg`,
`rotate_about_x`). It is a true no-op at 0. Only the ray/escape path is rotated;
a flat disc's plane test reads `get_z_cartesian` directly and would be misplaced
under a nonzero rotation, so disc scenes must keep the default.

**Why it cannot fully fix a full ring.** Every camera ray's geodesic plane
contains the line of sight, so for any pole there is always one meridian of rays
whose plane contains it; that line crosses a centered ring at (up to) two points.
A single chart can only *move* the gap, never remove it.

**Proper fix (not done): a coordinate atlas.** Two spherical charts with poles 90
degrees apart, switching mid-integration so no ray is pushed through a pole. This
generalizes to horizon-penetrating charts (Eddington-Finkelstein / Kruskal) for
interior views. Designed in `coordinate-atlas.md`.

## Problem 4: caustic magnification blow-up (resolution-dependent)

**What.** Near a caustic the tube footprint `B` collapses toward zero, so
`ratio = A/B` blows up, producing resolution-dependent fireflies on the ring.

**Measured.** For one star on the Einstein ring the effective magnification was
`~197x` at 1000px and `~77x` at 4000px (peak radiance `0.232` vs `0.091`, a 2.5x
swing).

**Root cause, two parts.**
1. **Coarse over-averaging.** A tube reports the *average* magnification over its
   sky patch, over-weighted by the part nearest the caustic. Finer tubes average
   over a smaller patch and converge toward the star's finite *local*
   magnification. This is a discretization artifact and converges with
   resolution.
2. **Genuine divergence.** A star exactly on the caustic has infinite
   magnification (measure zero) and never converges.

**Status: stopgapped (magnification cap).** `[star_catalog].magnification_cap`
(default none) clamps `ratio` (`Raytracer::magnification`). Measured with cap =
50: it leaves field stars untouched (peak `0.4804` capped == uncapped) and bounds
the per-tube blow-up, but it does **not** fully stabilize the ring across
resolution (caustic-star peak `0.127` @1000px vs `0.063` @2000px, still ~2x).

**Why the cap alone is not enough (key finding).** The cap bounds a *single
tube's* `A/B`, but the ring is built from *many overlapping tubes* (the arc is
gathered by many tubes near the fold, plus genuine multiple images). A coarse
pixel packs more of these capped contributions and sums them, so the per-*pixel*
brightness stays resolution-dependent, this is Problem 6's density/packing
effect reappearing at the ring. The cap kills fireflies; it does not make the
ring's surface brightness resolution-independent. That needs the PSF's
fixed-angular deposition (Problem 6) *in addition to* the cap. So the full recipe
is: cap the per-tube blow-up, then splat over a fixed angular kernel.

**Proper fixes (not done).**
- **Local magnification** at the star's sub-pixel position (interpolate the
  Jacobian from the corner directions) instead of the tube average: converges to
  the true finite value, reduces coarse over-averaging.
- **Finite source disc**: integrate the magnification over the star's angular
  disc. The `1/sqrt(distance)` fold singularity integrates to a finite value over
  a finite disc, so this regularizes the exact-caustic divergence physically. The
  cap then emerges as `~ 1 / source_size`.

## Problem 5: point-source-at-caustic singularity

**What.** A mathematical point source on a critical curve has divergent
magnification; no finite render converges (finer resolution keeps finding a
brighter, thinner ring pixel).

**Status: open.** The magnification cap (Problem 4) bounds it arbitrarily; the
principled cure is finite star size (a source disc, or the detector PSF's angular
width acting as an effective source scale). Distinct from the density issue in
Problem 6.

## Problem 6: resolution-dependent star brightness

**Important clarification.** A *single unlensed* star is NOT resolution
dependent: `ratio = A/B ~ 1` (footprint tracks the pixel), so its pixel value is
its flux at any resolution (measured field-star peak ratio `1.01` across 1000 vs
4000px). Resolution dependence appears only where `A/B != 1`:

1. **Density / packing.** A coarse pixel's large `B` swallows several stars, so
   the *aggregate* per-pixel brightness (and the exposure needed) shifts with
   resolution (~11x between 480x270 and 1600x900 per `plan-13`). Not a single
   star, a sum.
2. **Lensing.** Where `B` behaves differently from `A` (the caustic, Problem 4).

**Status: density is designed-away by PSF (not implemented), lensing is Problem
4.** The point-spread splat in `plan-13-star-psf.md` deposits each star's flux
over a fixed *angular* kernel into a shared buffer, making surface brightness
well-defined and resolution-independent by construction (also fixing temporal
twinkle and adding real diffraction glow). It changes *deposition*, not
collection.

**Scope caveat (important).** The PSF is a *deposition* upgrade and does not fix
Problem 2 (a *collection* drop) nor *correct* the Problem 4 magnification (it
still deposits the per-tube `A/B`). Collection and deposition are orthogonal:
"which tube owns the star" vs "how that star's light lands on the detector."

## Problem 7: the star gather ran single-threaded

**What.** `add_star_layer` was a serial nested loop; the geodesic pass is
parallel but the gather pinned one core. Cheap for the flat path, a hard
bottleneck for curved membership (its per-tube fold work near the ring ran on one
core, minutes at 1000px).

**Status: fixed.** Rewritten as a `par_iter` over pixels (each tube writes its own
pixel, reads the shared buffer). ~16x faster gather; curved 1000px went from
6+ min single-threaded to ~108s.

## Status summary

| Problem | Status | Proper fix |
|---|---|---|
| 1 Flat-quad ring gaps | curved membership (opt-in) | sliver-free tiling / inverse gather |
| 2 Fold star-drop | stopgap (max vs flat) | robust containment (Problem 1) |
| 3 Polar-axis singularity | stopgap (pole rotation) | coordinate atlas (two-patch) |
| 4 Caustic magnification blow-up | stopgap (magnification cap) | local magnification / finite source |
| 5 Point-source caustic singularity | open | finite source size / PSF width |
| 6 Resolution-dependent brightness | density: PSF (planned); lensing: Problem 4 | plan-13 PSF |
| 7 Serial star gather | fixed | - |

## Config / CLI added this session

- `--curved-star-membership` / `[star_catalog].curved_star_membership` (bool)
- `--pole-rotation-deg <deg>` (default 0)
- `[star_catalog].magnification_cap` (optional f64)

## Diagnostics that proved useful

- **Fold map:** paint every folded base tube bright to see the critical-curve
  locus (revealed the two rings).
- **Tube extraction:** dump one tube's corners, curved boundary, subdivision
  children, and gathered star to a file; plot in the tube's own PCA frame
  (microradian). Showed the footprint is a ~259:1 razor sliver with the star
  inside it, i.e. the star is genuinely in the footprint (so the drop was a real
  bug) and the sliver's near-zero area is the magnification blow-up.
- **Resolution ground truth, with a caveat:** comparing region *means* across
  resolutions is invalid for point stars (a star's flux sits in ~one pixel, so a
  4x finer grid spreads a region over 16x more pixels). Compare **peak radiance**
  instead (resolution-stable for the radiance-based gather: field star `1.01`,
  caustic star `2.54`).

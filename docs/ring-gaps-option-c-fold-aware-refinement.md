# Option C: Fold-aware tube refinement

Closing the sampling gaps in the lensed star ring by refining the tubes that
straddle the critical curve, and making the flat-quad star-membership test
trustworthy where it currently is not.

## Background: what the gap is

The star pipeline is a forward map. Each screen pixel casts a small "tube" of
four corner rays; the rays are integrated as geodesics; the four escape
directions form a spherical quad on the sky; every catalogue star whose
direction falls inside that quad is gathered and its flux is magnified by
`ratio = Omega_pixel / Omega_traced` (screen solid angle over traced-footprint
solid angle). See `compute_color` and `compute_star_collection_data` in
`src/rendering/raytracer.rs`.

A star sits at one fixed sky direction `S`. A tube collects it when `S` lands
inside the tube's footprint. The star's image (an Einstein ring for a source
directly behind the hole) is therefore the set of screen tubes whose footprint
covers `S`. For the ring to render as a closed curve, neighbouring footprints
must tile the sky with no slivers. Near the critical curve the magnification is
large, each base tube maps to a big, curved sky patch, and the flat-quad
approximation of adjacent tubes does not abut cleanly. Where a sliver opens over
`S`, no tube collects the star (a dark gap); where two footprints overlap `S`,
two tubes collect it (a too-bright spot).

Two separable defects produce this:

1. **The refinement trigger never fires at the primary ring.** For an
   all-escaped tube the only subdivision trigger is
   `winding_spread > winding_spread_threshold` (default `pi`, see
   `Scene::winding_spread_threshold`). Winding grows steeply only near the
   photon sphere, so the strongly wound higher-order rings do refine, but the
   primary Einstein ring, where winding is small and smooth, never trips the
   test. This is why depth 6 and depth 8 render a pixel-identical ring
   (measured RMSE 2.6e-7).
2. **Flat-quad membership leaves slivers.** Even a refined tube tests
   membership against the flat spherical quad of its four corners
   (`gather_triangle`, two triangles). Where the true footprint is large and
   curved, the quad is a poor stand-in and neighbours mistile.

Defect 1 must be fixed to even reach the gaps; defect 2 governs how cleanly the
residual closes.

## Core idea

Add a geometric refinement trigger that fires exactly where the forward map
stretches or folds, so those tubes subdivide until their footprints are small
enough that the flat-quad membership is accurate and neighbours tile. Concretely
refine a tube when either of these holds (in addition to the existing winding
trigger):

- **Footprint too large.** The largest corner-to-corner sky angle exceeds a
  threshold `max_footprint_angle`. A large footprint is precisely when the
  flat-quad membership is untrustworthy. This is the primary trigger.
- **Footprint folds.** The tube's image is non-convex or its orientation
  inverts across the quad, i.e. the two triangles that tile it (see
  `compute_solid_angle`) have opposite signed area. A sign flip means the
  critical curve passes through the tube's interior: one part of the tube images
  the sky right-side-out, the other mirror-reversed, and a flat quad cannot
  represent that. This is the fold-aware part and the reason for the name.

Both are cheap: the corner escape directions are already computed at
`raytracer.rs:396-419`, and the signed triangle areas are already computed
inside `compute_solid_angle` (`solid_angle_from_vecs` returns a signed angle;
the current code takes `.abs()` of each precisely so a fold does not cancel the
area to zero, which is the same signal we reuse as a trigger).

## Why the fold is the right signal

The critical curve is where the lensing Jacobian determinant changes sign: it
separates sky regions the observer sees right-side-out from regions seen
mirror-imaged (the secondary images). A screen tube that contains a piece of the
critical curve has corners on both sides, so the orientation of its mapped quad
is not well defined: the two tiling triangles wind in opposite senses. Detecting
that sign disagreement flags every tube the ring passes through, independent of
how much the ray wound. That is the coverage the winding trigger misses.

Refining a folded tube splits it until each child lies on one side of the
critical curve (or hits the depth cap). Once a child is unfolded and small, its
flat quad is a faithful footprint and its star membership is reliable.

## Algorithm

In `compute_color`, for the all-escaped branch, after the existing winding-spread
check and after computing the corner directions and the two triangle areas:

```
let a1 = solid_angle_from_vecs(&o_a, &o_b, &o_c)?; // signed
let a2 = solid_angle_from_vecs(&o_b, &o_d, &o_c)?; // signed
let folded = a1.signum() != a2.signum();
let footprint = max_corner_to_corner_angle(&o_a, &o_b, &o_c, &o_d);
if (folded || footprint > self.scene.max_footprint_angle)
    && let Some(color) = self.subdivide(sample_tube, depth)? {
    return Ok(color);
}
```

Everything downstream (the `ratio` magnification and the gather) is unchanged;
we only add reasons to subdivide. `solid_angle_from_vecs` currently returns a
positive angle from `atan2` but the sign of `num = a . (b x c)` is the
orientation; expose that sign (a two-line change) rather than only its
magnitude.

Termination:

- `max_subdivision_depth` caps recursion as today. With the new triggers, depth
  becomes a real convergence knob again: deeper allows thinner, more complete
  ring segments. Keep it modest (8-10); see cost below.
- A tube smaller than the footprint threshold and unfolded stops refining, so
  away from the ring nothing changes.

## Making the residual clean (defect 2)

Refinement alone shrinks the slivers but cannot remove them at the depth cap,
and it cannot represent the true magnification at the exact caustic (see
Limitations). Two optional add-ons close the residual:

- **Conservative membership at the leaves.** At a max-depth leaf that is still
  folded or large, dilate the quad slightly (a fraction of a footprint) so
  neighbours overlap rather than leave a sliver. This closes gaps but
  double-counts on overlaps, biasing the ring bright. To keep flux exact, weight
  each contribution by a partition of unity over the overlap (each star's flux
  split between the tubes that claim it in proportion to depth of containment).
  This is the more invasive half of the option.
- **Star-inside-tube fallback.** When a leaf is flagged folded, additionally
  query the octree for any star within the tube's bounding cap and test it
  against the true (curved) footprint by refining the membership test locally,
  instead of the flat quad. Guarantees no star inside the cap is missed.

The minimal version of Option C is just the two triggers; the add-ons trade code
for a visually seamless ring.

## Cost

Cost is `O(ring-length-in-pixels * 4^depth)`. Away from the ring nothing extra
happens, so for the sparse gallery frames (few near-caustic stars) the overhead
is negligible. For a bright star sitting on the caustic (the single-star test)
the whole ring refines to the cap, which is the case that detonated the earlier
`max_subdivision_depth = 12` run. Mitigations:

- Cap depth at 8-10, not 12.
- Optional per-frame or per-pixel subdivision budget: once a pixel has spent N
  traces, stop refining and accept the current leaves.
- The footprint threshold sets how aggressively tubes split; loosen it until the
  ring is acceptable rather than chasing the caustic.

## What the experiment says it buys

Rendering the single-star aligned case at rising resolution (the brute-force
equivalent of finer tubes) gives mean intensity 0.0743 (400px), 0.0759 (800px),
0.0686 (1600px) while peak intensity climbs 103 -> 406 -> 5408. The total
captured flux only wobbles about 10 percent even for a star placed exactly on
the caustic; the dramatic change is peak brightness of the thinnest ring pixels,
which a bloom or tone map largely clips to white. So Option C mostly buys ring
completeness and correctness under zoom or for bright near-caustic stars, not a
visible change to the current gallery frames.

## Limitations

- **The point-source caustic is singular.** An idealized point star exactly on
  the critical curve has infinite magnification, so no amount of refinement
  converges to a finite ring brightness; finer tubes just find brighter, thinner
  pixels (the diverging peak above). Option C makes the ring *closed and
  convergent away from the exact caustic*, but the true fix for the singularity
  is to give stars a finite angular size or point-spread function, which
  regularizes the magnification. That is orthogonal to C and would pair with it.
- **Still a forward map.** C improves the screen-side tiling but never
  guarantees that a given star's image is found; it relies on some tube's
  footprint covering `S`. Pathological sub-pixel arcs between refined leaves can
  still be missed. Option D removes that reliance entirely.

## Validation plan

- Single-star aligned: ring completeness (fraction of ring azimuth with nonzero
  flux) rises toward 1 as depth increases, and converges at fixed depth.
- Single-star offset: still exactly two images, primary brighter than secondary,
  positions unchanged (no new false images introduced).
- Flux conservation away from the caustic: for a star offset well off the ring,
  total gathered flux matches the analytic point-mass magnification within
  tolerance and is depth-independent.
- Gallery frames: diff against the current renders is confined to a thin annulus
  at the shadow; the broad field is bit-identical.

## Files to touch

- `src/rendering/raytracer.rs`: `compute_color` (add triggers), a
  `max_corner_to_corner_angle` helper, expose the sign from
  `solid_angle_from_vecs`; optionally the conservative-membership leaf path in
  `compute_star_collection_data`.
- `src/rendering/scene.rs`: add `max_footprint_angle` (and any fold tolerance)
  beside `winding_spread_threshold` / `max_subdivision_depth`.
- `src/configuration.rs`: serde field + default for the new threshold, defaulted
  loose so existing scenes are unchanged except at high magnification.

## Verdict

A small, local, physically principled patch. Fixes defect 1 outright and defect
2 to the depth limit, and restores depth as a convergence knob. It does not
resolve the point-source caustic singularity (needs finite star size) and does
not guarantee image completeness (needs Option D). Recommended as the practical
fix if the goal is a faithful, zoomable ring at bounded cost.

# Coordinate atlas (multi-chart integration)

Cover the spacetime with an atlas of overlapping coordinate charts, integrate
each geodesic in whichever chart is regular locally, and switch charts (with a
coordinate + momentum transformation) in the overlap. This is the principled
generalization of the pole-rotation prototype, and it also unlocks
horizon-crossing and interior views.

## Why

The integrator currently works in a single chart per geometry (Schwarzschild in
`Spherical`, Kerr in `BoyerLindquist` / Kerr-Schild). Single charts carry
coordinate singularities that corrupt geodesics passing through them, even
though the spacetime there is perfectly smooth:

- **Polar axis** (`theta = 0/pi`): the spherical metric has `1/sin(theta)`
  terms and `dphi/dlambda` blows up, so rays grazing the axis return corrupted
  escape directions or NaN. This shows up as a straight bad meridian through the
  shadow centre and, for a full Einstein ring, a gap where that meridian crosses
  the ring. The `POLE_ROT_DEG` prototype (see below) relocates it but, with one
  chart, can never remove it: the meridian always crosses a centred ring.
- **Event horizon** (`r = 2M` in Schwarzschild coordinates): a coordinate, not
  physical, singularity. The metric degenerates there, so you cannot integrate
  across it. Fine for the shadow image (horizon-crossing rays are just
  `Captured`), but a hard wall for anything *inside* the hole.

An atlas removes both: no ray is ever pushed through a chart's singular region,
because a neighbouring chart that is regular there takes over.

## The idea

A manifold is covered by overlapping charts; in each overlap a transition map
relates the two coordinate systems. A geodesic is a physical curve, independent
of chart, so its coordinate representation can be swapped mid-integration
wherever a chart goes bad, seamlessly, as long as the transition map is exact.

Three pieces:

1. **A chart** carries its coordinates, its metric and Christoffel symbols, a
   validity region (where it is well-conditioned), and a map to a common
   reference (Cartesian, or physical invariants like `r`). The existing
   `CoordinateSystem` enum (`Cartesian`, `Spherical`, `BoyerLindquist`) plus the
   conversion helpers are the seed of this.
2. **Transition maps.** In an overlap, transform both the position and the
   momentum. Position follows the coordinate map; momentum follows its Jacobian,
   `p'^mu = (partial x'^mu / partial x^nu) p^nu`. For analytic chart pairs
   (Schwarzschild <-> Eddington-Finkelstein <-> Kruskal; rotated-spherical pairs)
   the Jacobian is closed-form; otherwise finite-difference the coordinate map.
3. **A switch criterion** in the integrator loop. Monitor a validity signal
   (`theta` near `0/pi`, `r` near `2M`, or the metric condition number) and, when
   the point enters the overlap near the current chart's bad region, apply the
   transition and continue in the neighbour. The affine parameter is
   chart-independent, so it just carries on; re-checking the null condition
   `g_{mu nu} p^mu p^nu = 0` after a switch is a cheap correctness assert.

Switches happen only at region boundaries, so the runtime cost is negligible.

## What it unifies

The polar-axis fix and the horizon fix are the *same* machinery, different chart
pairs. Building the chart + transition + switch infrastructure once yields:

- **Pole-free full rings** — two spherical charts whose poles are 90 degrees
  apart; each covers the other's bad meridian. (The `POLE_ROT_DEG` prototype is
  the one-chart, relocate-only version of this.)
- **Horizon crossing and interior views** — a horizon-penetrating chart takes
  over near `r = 2M`.
- **Kerr** — already integrates in Kerr-Schild form (horizon-penetrating), so
  Kerr interiors are partly handled; the axis charts still apply.
- **Ellis wormhole** — naturally two charts glued at the throat, already a
  two-patch manifold; the atlas is its native description.

## Chart catalogue (Schwarzschild)

- **Exterior**: `Spherical` (or isotropic) for `r > 2M + epsilon`. Cheap, what we
  have.
- **Rotated spherical**: a second `Spherical` chart rotated so its pole sits
  where the first chart's field of view needs coverage. Transition is a pure
  rotation (the `rotate_about_x` prototype generalized to an arbitrary axis).
- **Eddington-Finkelstein (ingoing)**: the simplest horizon-penetrating chart,
  regular at `r = 2M`, covers exterior + interior down to the singularity. The
  natural first interior chart to add.
- **Painleve-Gullstrand**: also horizon-penetrating, flat spatial slices, the
  intuitive "river" picture.
- **Kruskal-Szekeres**: the maximal extension (exterior, interior, white hole,
  parallel universe). Reach for it only when you actually want the analytically
  extended spacetime; for imaging a real hole it is overkill, EF suffices.

For imaging *from outside*, interior charts matter only when you want to see
inside (an infalling observer's view, the approach to the singularity, interior
structure). Otherwise a horizon-crossing ray is classified `Captured` and
stopped, no interior chart required.

## Architecture in this codebase

- **`Geometry` becomes a chart provider.** Given a point, it reports which chart
  is valid and how to transition to a neighbour, instead of assuming one fixed
  coordinate system. The metric / Christoffel methods dispatch per chart.
- **The integrator loop gains a switch step.** After each step (or when a
  validity signal trips), ask the geometry whether to switch charts; if so,
  transform the ray state and continue. Add a hysteresis band so a ray skimming
  a boundary does not thrash between charts.
- **Transitions reuse the conversion layer.** `spherical_coordinates_helper` and
  `four_vector::get_cartesian_vector` already convert positions and momenta
  through Cartesian; a transition is "to common reference, then into the target
  chart," with the momentum Jacobian applied. The `POLE_ROT_DEG` prototype
  already threads a rotation through exactly these functions.
- **Escape / capture classification stays.** A ray that reaches large `r`
  escapes (its chart's asymptotic region); one that reaches the singularity in
  an interior chart terminates there.

## Relationship to the existing prototypes

- **Pole rotation** (`POLE_ROT_DEG`, `frame_rotation_deg` / `rotate_about_x` in
  `spherical_coordinates_helper`, plus the inverse rotation in
  `four_vector::get_cartesian_vector`): a single rotated chart. `0` is the
  identity (every render byte-identical); nonzero relocates the pole. It is the
  degenerate one-chart case of the atlas: it *moves* the singularity but cannot
  remove it for a full ring. The atlas's two-chart version removes it.
- **Curved-boundary star membership** (the `CURVED_MEMBERSHIP` prototype): fills
  ring gaps caused by flat-quad tiling. Orthogonal to the atlas: curved
  membership fixes the *gather* where the rays are valid; the atlas fixes the
  *rays* so they are valid everywhere. They compose.

## Caveats

- **Overlaps must be real.** Switch *inside* the overlap, never exactly at a
  boundary, or you switch on singular data.
- **Transition maps should be exact.** Analytic Jacobians avoid kinks; a
  finite-difference Jacobian near a coordinate singularity is itself noisy.
- **Hysteresis on the switch.** A band around the switch threshold stops a
  boundary-skimming ray from oscillating between charts.
- **Downstream consumers of raw coordinate components.** Anything that reads a
  coordinate component directly rather than through the conversion layer must be
  made chart-aware. The flat disc already does this (`disc.rs` uses
  `get_z_cartesian` for its plane-crossing test), so a rotated or swapped chart
  would misplace it until it is routed through the same transform, exactly the
  caveat the `POLE_ROT_DEG` prototype hit.
- **Null-condition drift.** Re-normalise / assert `g p p = 0` after a transition
  to catch a wrong Jacobian early.

## Suggested phases

1. Generalize the pole rotation to an arbitrary-axis rotation and add the second
   spherical chart with a rotation transition; switch on `theta` near `0/pi`.
   Delivers provably seam-free full rings.
2. Add ingoing Eddington-Finkelstein as an interior chart with a Schwarzschild
   <-> EF transition; switch on `r` near `2M`. Delivers horizon crossing and
   basic interior views.
3. Make `Geometry` a first-class chart provider and route all raw-component
   consumers (e.g. `disc.rs`) through the transform.
4. Optional: Kruskal for the maximal extension; per-geometry axis charts for
   Kerr and the Ellis wormhole.

## Verdict

The principled home for every coordinate-singularity problem in the renderer:
the polar axis, the horizon, the interior, and the naturally multi-chart Ellis
wormhole all reduce to "pick the regular chart and transition in the overlap."
Larger than a patch (it touches the integrator and the geometry abstraction),
but each phase is independently useful, and phase 1 alone supersedes the
pole-rotation workaround with a provably seam-free result.

# Choosing the sky behind the hole

Every starfield render puts the lensed structure against some patch of the Gaia
DR3 catalogue, and the patch is decided entirely by where the camera sits. This
directory records which patches are worth using, how to aim at them, and the
camera conventions that are easy to get wrong.

The frames here are flat space: `geometry_type.Euclidean`, no hole and no disc
(`sky-survey.toml`). Rays travel straight, so each frame is the raw sky in the
direction the camera faces, which is the patch a lensed scene would spread
around its ring. They render in about 45 s each at 1920x1080.

Every PNG in this directory carries its own parameters in PNG text chunks:

```
python3 -c "from PIL import Image; print(Image.open('images/sky-survey/candidate-A-combined-best.png').text)"
```

## The catalogue scan

Rather than sampling frames by eye, `data/gaia_mag12.parquet` was binned into a
1-degree equal-area map of star counts (the catalogue stores `ra`/`dec`, and
`star_catalog.rs:185` turns those into scene Cartesian directions unchanged, so
sky direction and scene direction are the same vectors). Each candidate view was
scored over the real frame footprint:

- **brightness** = mean log star density inside the frame
- **contrast** = 90th minus 10th percentile of that density, so a bright band
  crossing dark sky scores high and a uniformly dense field scores low

| frame | camera | brightness | contrast | what is in it |
|-------|--------|-----------|----------|----------------|
| 01 | `-17,0,0` | 1.29 | 0.29 | **the vantage every current black-hole render uses**: the band only clips the top-left corner |
| 02 | `17,0,0` | 1.25 | 0.34 | emptiest of the survey |
| 03 | `0,-17,0` | 1.78 | 0.58 | band on a diagonal, dark lane, star cloud low left |
| 04 | `0,17,0` | 1.78 | 0.64 | rich band, near-vertical through the frame |
| 05 | `0,0,17` | 1.69 | 0.69 | band across the lower half, bright knot low left |
| 06 | `0,0,-17` | 1.54 | 0.43 | band horizontal across the middle; the axis renders use this |
| 07 | `-12,-12,0` | 1.30 | 0.32 | sparse |
| 08 | `-12,12,0` | 1.67 | 0.64 | band with heavy dust texture |
| 09 | `12,-12,0` | 1.58 | 0.66 | broad band low left |
| 10 | `12,12,0` | 1.36 | 0.36 | rich part jammed into the upper-right corner |
| 11 | `-12,0,12` | 1.29 | 0.30 | sparse |
| 12 | `-12,0,-12` | 1.76 | 0.63 | band low in frame, upper half dark |

Contact sheet: `survey-12-directions.png`.

## Scan picks

The five highest-scoring directions the scan found, rendered in
`scan-candidates.png` and stored individually:

| pick | camera | brightness | contrast | character |
|------|--------|-----------|----------|-----------|
| A | `4.58,14.09,8.33` | 1.85 | 0.80 | best combined: bright band cut by a dark rift on the diagonal |
| B | `12.56,0,11.45` | 1.79 | **0.81** | highest contrast in the sky: band across the top, dark below |
| C | `3.22,-1.05,16.66` | 1.86 | 0.78 | broad bright band corner to corner |
| D | `-2.14,10.06,13.54` | 1.87 | 0.73 | band steeply diagonal, dust lanes through it |
| E | `4.18,9.40,13.54` | **2.01** | 0.51 | densest field in the catalogue, less internal structure |

For comparison, the vantage in use today scores 1.29 / 0.29, the lowest of
everything measured.

## Adjusted framings

Starting from survey frames and adjusted by hand:

| name | camera | extra | change |
|------|--------|-------|--------|
| `adjusted-04-roll-ccw50` | `0,17,0` | `--psi=-0.873` | band swung ~50 degrees counter-clockwise, from +66 to +14 degrees |
| `adjusted-05-panned-down` | `8.16,0,14.92` | | bright knot moved from y = 0.82 to y = 0.51 |
| `adjusted-10-panned-right` | `6.63,15.62,0` | | rich region moved from x = 0.77 to x = 0.63 |
| `adjusted-12-panned-down` | `-8.87,0,-14.47` | | band moved from y = 0.63 to y = 0.49 |

## Aiming without moving the camera

Camera position decides two unrelated things at once: the inclination to the
hole's spin axis, and which patch of sky ends up behind it. `[star_catalog.rotation]`
separates them by rotating the celestial sphere instead. To put the B patch
behind a camera that stays on the equatorial `-x` axis:

```toml
[star_catalog.rotation]
from = [-0.3782, 0.0, -0.9257]   # what the B vantage sees behind the hole
to = [1.0, 0.0, 0.0]             # the view axis of a camera at -17,0,0
```

The `from` vector for any row in the tables above is just the normalised
negative of its camera position, since these cameras face the origin. The
rotation leaves the roll about the view axis free, so pair it with `--psi` when
a particular orientation matters.

## Camera conventions

- **Aiming.** In flat space the camera faces the origin, so to look at sky
  direction `d`, put the camera at `-17 * d`. The scan tables above already give
  camera positions, not directions.
- **Roll** is `--psi`, and it is clean: negative rolls the content
  counter-clockwise, positive clockwise. `--psi=-0.873` gives ~50 degrees CCW.
- **Panning is not `--theta`.** `--theta` tilts the view axis toward the
  *pre-rolled horizontal image axis* of the geometry's own tetrad, so it rotates
  the frame as well as panning it: from `12,12,0`, both signs moved content left
  and one also swung the band by 40 degrees. Pan by moving the camera on its
  sphere instead.
- **A 10 degree camera move pans the view about 5 degrees**, measured, not the
  10 you would expect. Solve for the angle from that ratio and verify.
- **The poles are degenerate.** At `0,0,±R` the up vector has no preferred
  choice; moving that camera 10 degrees about x re-picks it and rolls the frame
  by ~60 degrees. Moving about y behaves normally.
- **Angles do not transfer between geometries.** These flat-space frames need no
  angles at all. Kerr-Schild (`geometry_type.Kerr`) needs `--theta=1.5708
  --psi=-1.5708` for an equatorial view and `--theta=0` down the axis; the
  KerrBL scenes' `--theta=-3.14159` points away from the hole under the
  Kerr-Schild tetrad and puts it off-frame entirely.

## Grade

All frames here: `--white 0.11042 --bloom 0.2 --tonemap aces`, with the white
point fixed across the set so brightness differences between directions survive
the grade instead of being normalised away. The value is the median of the
per-frame 99th percentiles at 1920x1080; it is resolution-dependent, because a
larger pixel gathers more star flux, so a 640x360 probe needs a much higher
white for the same look.

## Reproducing a frame

```
target/release/gr_raytracer --width=1920 --height=1080 \
  --camera-position=4.58,14.09,8.33 --max-radius=200 --max-steps=20000 \
  --config-file=images/sky-survey/sky-survey.toml render --filename=out.hdr
uv run scripts/grade.py out.hdr out.png --white 0.11042 --bloom 0.2 --tonemap aces
```

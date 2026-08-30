#!/usr/bin/env bash
# Reproduces the Kerr showcase/diagnostic images (images/showcase_*.png):
# spin-vs-no-spin pairs that isolate individual GR effects using synthetic
# textures. Every pair shares its camera exactly, so any difference between
# the two images of a pair is purely the spin (a = 0.499, i.e. a/M = 0.998,
# vs a = 0; scene units use r_s = 1, so M = 0.5).
#
# Sets rendered:
#   1. shadow      - equatorial camera, checker sky, no disc: shadow shape
#                    comparison (near-extremal D-shape vs Schwarzschild circle).
#   2. framedrag   - camera on the spin axis, checker sky: frame dragging
#                    twists the sky's radial spokes; straight at a = 0.
#   3. framedrag_zoom - same camera, 3x section zoom onto the shadow, where
#                    the twist concentrates (it falls off ~1/r^3).
#   4. ringspoke   - flat disc with a rainbow ring-and-spoke bitmap on a
#                    black sky: rings are coded by ABSOLUTE radius (same
#                    color = same Boyer-Lindquist r in both images; each
#                    scene's bitmap is pre-warped for its annulus by
#                    scripts/generate_showcase_textures.py), spokes reveal
#                    azimuthal shear between lensing image orders. The spin
#                    scene's inner edge hugs the prograde ISCO (r = 0.62);
#                    the a = 0 scene starts outside its ISCO (r = 3.05).
#
# Usage:  images/create-showcase-images.sh [set ...]
#   With no arguments renders all sets. Set names: shadow framedrag
#   framedrag_zoom ringspoke.
#   SIZE env var overrides the base resolution (default 900).
#   The Kerr (Kerr-Schild) chart is slow: expect ~15-20 min per image at
#   the default size, ~2 hours for everything.
#
# The scene TOMLs live in scene-definitions/showcase-*.toml; the disc
# bitmaps are regenerated with:
#   uv run --with pillow python scripts/generate_showcase_textures.py
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

SIZE="${SIZE:-900}"
BIN="target/release/gr_raytracer"
SETS=("$@")
[ ${#SETS[@]} -eq 0 ] && SETS=(shadow framedrag framedrag_zoom ringspoke)

if [ ! -x "$BIN" ]; then
    echo "Building release binary..."
    cargo build --release
fi

render() { # scene camera-args... -> filename
    local scene="$1" out="$2"; shift 2
    echo "Rendering $out ..."
    "$BIN" --config-file "scene-definitions/$scene" "$@" render \
        --filename "images/$out"
}

for set in "${SETS[@]}"; do
    case "$set" in
    shadow)
        # Equatorial view; the a = 0.499 shadow flattens on the prograde
        # side and shifts off-center; a = 0 is a centered circle.
        CAM=(--width="$SIZE" --height="$SIZE" \
             --camera-position=-18.0,0.0,0.0 --theta=1.5708 --psi=-1.5708 --phi=0.0)
        render showcase-kerr-checker-sky.toml    showcase_kerr_shadow_spin.png "${CAM[@]}"
        render showcase-kerr-checker-sky-a0.toml showcase_kerr_shadow_a0.png   "${CAM[@]}"
        ;;
    framedrag)
        # Pole-on view down the spin axis; sky spokes twist near the shadow
        # at a = 0.499, ruler-straight at a = 0.
        # NEVER place the camera exactly on the axis (x=y=0): the camera
        # tetrad degenerates there (in the KS chart too, not only BL) and
        # the winding band around the shadow silently disappears; the 0.15
        # offset (0.5 degrees of parallax) is the workaround.
        CAM=(--width="$SIZE" --height="$SIZE" \
             --camera-position=0.15,0.0,-18.0 --theta=0.0 --psi=0.0 --phi=0.0)
        render showcase-kerr-checker-sky.toml    showcase_kerr_framedrag_spin.png "${CAM[@]}"
        render showcase-kerr-checker-sky-a0.toml showcase_kerr_framedrag_a0.png   "${CAM[@]}"
        ;;
    framedrag_zoom)
        # Same pole-on camera at 3x virtual resolution, rendering only the
        # central third (section render): a 3x zoom onto the shadow region.
        Z=$((SIZE * 3)); LO=$((SIZE)); HI=$((SIZE * 2))
        # Same 0.15 off-axis offset as the framedrag set (axis-degenerate
        # camera tetrad).
        CAM=(--width="$Z" --height="$Z" \
             --camera-position=0.15,0.0,-18.0 --theta=0.0 --psi=0.0 --phi=0.0)
        SEC=(--from-row="$LO" --to-row="$HI" --from-col="$LO" --to-col="$HI")
        echo "Rendering showcase_kerr_framedrag_zoom_spin.png ..."
        "$BIN" --config-file scene-definitions/showcase-kerr-checker-sky.toml \
            "${CAM[@]}" render "${SEC[@]}" \
            --filename images/showcase_kerr_framedrag_zoom_spin.png
        echo "Rendering showcase_kerr_framedrag_zoom_a0.png ..."
        "$BIN" --config-file scene-definitions/showcase-kerr-checker-sky-a0.toml \
            "${CAM[@]}" render "${SEC[@]}" \
            --filename images/showcase_kerr_framedrag_zoom_a0.png
        ;;
    ringspoke)
        # Slightly-inclined view of the ring-and-spoke disc; theta = 1.56
        # centers the composition.
        CAM=(--width="$SIZE" --height="$SIZE" \
             --camera-position=-18.0,0.0,2.5 --theta=1.56 --psi=-1.5708 --phi=0.0)
        render showcase-kerr-ringspoke-disc.toml    showcase_kerr_ringspoke_spin.png "${CAM[@]}"
        render showcase-kerr-ringspoke-disc-a0.toml showcase_kerr_ringspoke_a0.png   "${CAM[@]}"
        ;;
    *)
        echo "Unknown set: $set (known: shadow framedrag framedrag_zoom ringspoke)" >&2
        exit 1
        ;;
    esac
done

echo "Done. Outputs in images/showcase_*.png"

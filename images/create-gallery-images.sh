#!/usr/bin/env bash
# Reproduces the Kerr gallery images in images/kerr.md (frozen recipes,
# 2026-08-30, post plan-02/NT-fix renderer). Companion to
# create-main-image.sh (README hero) and create-showcase-images.sh
# (spin-vs-no-spin diagnostics); the spin-sweep animations come from
# scripts/rendering/create_kerr_images.sh.
#
# Usage:  images/create-gallery-images.sh [set ...]
#   Sets: checker volumetric grid_gif   (default: all)
#   All scenes are derived from committed scene-definitions at run time;
#   the flat-disc inner radius tracks 1.02 * r_isco(a).
#
# Not scripted here: kerr_trajectory_near_horizon.png (produced by the
# ray-export + scripts/plotting pipeline; see scripts/Readme.md).
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

BIN="target/release/gr_raytracer"
SETS=("$@")
[ ${#SETS[@]} -eq 0 ] && SETS=(checker volumetric grid_gif)

if [ ! -x "$BIN" ]; then
    echo "Building release binary..."
    cargo build --release
fi

TMPDIR_SCENES="$(mktemp -d -t gallery-scenes-XXXXXX)"
trap 'rm -rf "$TMPDIR_SCENES"' EXIT

isco_inner() { # a (r_s units) -> 1.02 * r_isco
    python3 -c "
import math
a = 2.0 * $1
z1 = 1.0 + (1.0 - a*a)**(1.0/3.0) * ((1.0+a)**(1.0/3.0) + (1.0-a)**(1.0/3.0))
z2 = math.sqrt(3.0*a*a + z1*z1)
print(f'{1.02 * (3.0 + z2 - math.sqrt((3.0-z1)*(3.0+z1+2.0*z2))) * 0.5:.3f}')"
}

for set in "${SETS[@]}"; do
    case "$set" in
    checker)
        # Near-edge-on flat blackbody disc against the labeled checker sky,
        # moderate and near-extremal spin.
        for a in 0.250 0.499; do
            scene="$TMPDIR_SCENES/checker_$a.toml"
            # The checkerboard DISC (disk.png) is the subject, per the
            # original images; inner edge tracks the ISCO.
            sed -e "s/^a = .*/a = $a/" \
                -e "s/^inner_radius = .*/inner_radius = $(isco_inner "$a")/" \
                scene-definitions/kerr-checker-disc.toml > "$scene"
            case "$a" in
                0.250) out=images/render_kerr_checker_texture.png ;;
                *)     out=images/render_kerr_large_a_checker_texture.png ;;
            esac
            echo "Rendering $out ..."
            "$BIN" --width=1000 --height=1000 --max-steps=1000000 \
                --camera-position=-19,0,1.4 --theta=-3.14159 --psi=0 --phi=0 \
                --config-file "$scene" render --filename="$out"
        done
        ;;
    volumetric)
        # The streaky volumetric disc (KerrBL) against the checker sky and
        # against the NGC 6355 star cluster.
        NGC="resources/tmp/hubble_ngc6355_potw2301a.jpg"
        if [ ! -f "$NGC" ]; then
            echo "Downloading NGC 6355 background from Wikimedia Commons..."
            mkdir -p "$(dirname "$NGC")"
            curl -L --fail -o "$NGC" \
                "https://commons.wikimedia.org/wiki/Special:FilePath/NGC6355%20-%20HST%20-%20Potw2301a.jpg"
        fi
        for bg in checker stars; do
            scene="$TMPDIR_SCENES/volumetric_$bg.toml"
            case "$bg" in
                checker) sky="resources/celestial.png"
                         out=images/render_kerr_checker_texture_volumetric.png ;;
                *)       sky="resources/tmp/hubble_ngc6355_potw2301a.jpg"
                         out=images/render_kerr_stars_volumetric.png ;;
            esac
            sed -e "s|resources/celestial.png|$sky|" \
                scene-definitions/kerr-bl-volumetric-streaky.toml > "$scene"
            echo "Rendering $out ..."
            "$BIN" --width=1500 --height=1500 --max-steps=1000000 \
                --camera-position=-20,0,-0.6 --theta=-3.14159 --psi=0 --phi=0 \
                --config-file "$scene" render --filename="$out"
        done
        ;;
    grid_gif)
        # Spin sweep with the coordinate-grid overlay.
        # Cyan scene + black labels (matching the historical grid gif):
        # the measurement grid reads cleanly on a flat background.
        rm -rf kerr_images
        SWEEP_SCENE=scene-definitions/kerr-animation-cyan.toml \
        SWEEP_EXPOSURE=6 SWEEP_CAMERA=-14,0,1.2 SWEEP_LABEL_COLOR=black \
        SWEEP_OUTPUT=images/kerr_animation_grid.gif \
            bash scripts/rendering/create_kerr_images.sh "$BIN" --grid --grid-size 50
        ;;
    *)
        echo "Unknown set: $set (known: checker volumetric grid_gif)" >&2
        exit 1
        ;;
    esac
done

echo "Done."

#!/usr/bin/env bash
# Reproduces the README main image (kerr_black_hole_with_stars.png):
# a near-extremal Kerr black hole (a/M = 0.998) with an ISCO-hugging
# volumetric accretion disc in front of the M25 open-cluster star field.
#
# The star-field background is NOT checked into the repository; this script
# downloads it from Wikimedia Commons on first use
# (https://commons.wikimedia.org/wiki/File:Messier_object_025.jpg).
#
# Usage:  images/create-main-image.sh [output.png]
#   WIDTH/HEIGHT env vars override the resolution (default 1280x720).
#   Rendering takes roughly half an hour at the default resolution.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

BACKGROUND="resources/tmp/Messier_object_025.jpg"
OUTPUT="${1:-kerr_black_hole_with_stars.png}"
WIDTH="${WIDTH:-1280}"
HEIGHT="${HEIGHT:-720}"

if [ ! -f "$BACKGROUND" ]; then
    echo "Downloading M25 star-field background from Wikimedia Commons..."
    mkdir -p "$(dirname "$BACKGROUND")"
    curl -L --fail -o "$BACKGROUND" \
        "https://commons.wikimedia.org/wiki/Special:FilePath/Messier_object_025.jpg"
fi

SCENE="$(mktemp -t kerr-main-image-XXXXXX).toml"
trap 'rm -f "$SCENE"' EXIT

# Frozen recipe (2026-07-12): scene-definitions/kerr-volumetric-streaky.toml
# with the inner edge pulled to just outside the prograde ISCO (0.618 r_s at
# this spin), the temperature raised so the whole disc stays in the visible
# range, the brightness reference chosen so the Reinhard tone mapping does not
# clip away the Perlin texture, and the star-field background. The bright side
# of the disc is the approaching (left) side; renders made before the
# redshift-sign fix show it mirrored.
cat > "$SCENE" <<'EOF'
celestial_temperature = 0.0

[celestial_texture.Bitmap]
beaming_exponent = 3.0
path = "resources/tmp/Messier_object_025.jpg"

[geometry_type.Kerr]
radius = 1.0
a = 0.499
horizon_epsilon = 1e-4

[[objects]]

[objects.VolumetricDisc]
inner_radius = 0.8
outer_radius = 16.0
temperature = 12000.0
num_octaves = 8
max_steps = 50000
step_size = 0.0002
thickness = 0.03
density_multiplier = 500.0
brightness_reference_temperature = 6000.0
absorption = 0.3
scattering = 0.4
noise_scale = [60.0, 4.0, 45.0]
noise_offset = -0.45

[objects.VolumetricDisc.texture.BlackBody]
beaming_exponent = 0.0
EOF

cargo build --release

./target/release/gr_raytracer \
    --width="$WIDTH" --height="$HEIGHT" \
    --camera-position=-18,0,0.8 --theta=1.52 --psi=-1.57 --phi=0 \
    --config-file "$SCENE" \
    render --filename="$OUTPUT"

echo "Wrote $OUTPUT (${WIDTH}x${HEIGHT})"

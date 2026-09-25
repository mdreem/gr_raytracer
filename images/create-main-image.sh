#!/usr/bin/env bash
# Reproduces the README main image (kerr_black_hole_with_stars.png):
# a near-extremal Kerr black hole (a/M = 0.998; here a = 0.499 with r_s = 1) with
# an ISCO-hugging volumetric accretion disc in front of the real Gaia DR3 star
# field, graded with bloom + a luminance-preserving ACES tone map.
#
# The star field is the real Gaia DR3 catalogue (G <= 12). It is NOT checked into
# the repository; this script downloads it to data/gaia_mag12.parquet on first
# use via scripts/gaia/download.py.
#
# Usage:  images/create-main-image.sh [output.png]
#   WIDTH/HEIGHT env vars override the resolution (default 1280x720).
#   TEMPERATURE overrides the peak disc temperature (default 10000.0).
#   WHITE/EXPOSURE/BLOOM override the grade (defaults 1.2 / 1.0 / 0.12): lower
#   WHITE lifts the star field but brightens the disc; BLOOM is the glow strength.
#   Rendering takes roughly half an hour at the default resolution.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

CATALOGUE="data/gaia_mag12.parquet"
OUTPUT="${1:-kerr_black_hole_with_stars.png}"
WIDTH="${WIDTH:-1280}"
HEIGHT="${HEIGHT:-720}"
TEMPERATURE="${TEMPERATURE:-10000.0}"
WHITE="${WHITE:-1.2}"
EXPOSURE="${EXPOSURE:-1.0}"
BLOOM="${BLOOM:-0.12}"

if [ ! -f "$CATALOGUE" ]; then
    echo "Downloading Gaia DR3 (G <= 12) star catalogue to $CATALOGUE ..."
    uv run --group gaia scripts/gaia/download.py download \
        --max-magnitude 12 --output "$CATALOGUE"
fi

SCENE="$(mktemp -t kerr-main-image-XXXXXX).toml"
HDR="$SCENE.hdr"
trap 'rm -f "$SCENE" "$HDR"' EXIT

# Frozen recipe: KerrBL a = 0.499 with the volumetric disc inner edge just outside
# the prograde ISCO. The disc is dimmed (density_multiplier 180, down from the
# earlier 500) and the star flux scaled up (flux_scale 120) so the real Gaia field
# reads alongside the disc rather than being washed out. Peak temperature 10000 K
# keeps a golden-orange palette. The bright side of the disc is the approaching
# (left) side.
#
# The in-engine tone map (Reinhard/GlobalLinear in src/rendering/color.rs) has no
# bloom, so we render a LINEAR .hdr and grade it with scripts/grade.py, which adds
# multi-scale bloom and a luminance-preserving ACES filmic curve.
cat > "$SCENE" <<EOF
celestial_temperature = 0.0

[adaptive_sampling]
enabled = true
samples_per_axis = 2
minimum_luminance = 0.0

[celestial_texture.BlackBody]
beaming_exponent = 0.0

[geometry_type.KerrBL]
radius = 1.0
a = 0.499
horizon_epsilon = 1e-4

[star_catalog]
path = "$CATALOGUE"
flux_scale = 120.0
max_subdivision_depth = 6

[[objects]]

[objects.VolumetricDisc]
inner_radius = 0.8
outer_radius = 16.0
temperature = $TEMPERATURE
num_octaves = 8
max_steps = 50000
step_size = 0.0002
thickness = 0.03
density_multiplier = 180.0
brightness_reference_temperature = 7000.0
absorption = 0.3
scattering = 0.4
noise_scale = [60.0, 4.0, 45.0]
noise_offset = 0.0

[objects.VolumetricDisc.texture.BlackBody]
beaming_exponent = 0.0
EOF

cargo build --release

# Render a linear HDR (writing .hdr bypasses the in-engine tone map), then grade.
./target/release/gr_raytracer \
    --width="$WIDTH" --height="$HEIGHT" \
    --camera-position=-17,0,1.5 --theta=-3.14159 --psi=0.0 --phi=0 \
    --config-file "$SCENE" \
    render --filename="$HDR"

uv run scripts/grade.py "$HDR" "$OUTPUT" \
    --white "$WHITE" --exposure "$EXPOSURE" --bloom "$BLOOM"

echo "Wrote $OUTPUT (${WIDTH}x${HEIGHT})"

#!/usr/bin/env bash
# Reproduces the close-vantage series (frozen 2026-08-30): first-person
# camera placements near the hole, Schwarzschild vs near-extremal Kerr
# (ZAMO observer) as comparison pairs where both exist. Scenes:
# scene-definitions/vantages/ (3500 K disc - close-in observers see
# gravitationally blueshifted gas, so the emitted temperature is kept low
# to stay in warm hues after the shift).
# Kerr camera positions are CARTESIAN; on the equator x^2+y^2 = r_BL^2+a^2,
# so e.g. BL r=0.70 needs cartesian 0.86 at a=0.499.
set -euo pipefail
REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"
BIN="target/release/gr_raytracer"
[ -x "$BIN" ] || cargo build --release
SV=scene-definitions/vantages/schwarzschild-vantage.toml
KV=scene-definitions/vantages/kerr-zamo-vantage.toml
render() { # out scene pos theta exposure
  echo "Rendering $1 ..."
  # Tight integrator tolerance: disc hits are tested on straight chords
  # between adaptive steps, and the grazing plane-dips these close-in
  # cameras produce get skipped at the default epsilon, leaving ragged
  # false-sky bites along the wound disc edges.
  "$BIN" --width=800 --height=800 --max-steps=1000000 --epsilon=1e-9 \
    --camera-position="$3" --theta="$4" --psi=0 --phi=0 --exposure="$5" \
    --config-file "$2" render --filename="images/$1"
}
# Schwarzschild vantages
render vantage_porthole.png            "$SV" "-1.05,0,0"  0        4
render vantage_photonsphere_edge.png   "$SV" "-1.25,0,0"  -1.5708  4
render vantage_grazing_winding.png     "$SV" "-1.55,0,0"  -2.1     4
render vantage_rearview.png            "$SV" "-2.0,0,0.4" -0.5     4
# Kerr companions (ZAMO)
# The photonsphere/grazing cameras have BL radii (0.808, 1.003) inside
# the disc annulus, so at z=0 they would sit exactly ON the razor-thin
# disc surface and which side they see depends on floating-point noise
# in the chart conversion; the 0.02 lift keeps them just above it.
render vantage_kerr_dragged_porthole.png "$KV" "-0.75,0,0"   0       4
render vantage_kerr_photonsphere_edge.png "$KV" "-0.95,0,0.02" -1.75 4
render vantage_kerr_grazing_winding.png "$KV" "-1.12,0,0.02" -2.25   4
render vantage_kerr_rearview.png       "$KV" "-2.55,0,0.4" -0.5     4
# Kerr-only vantages
render vantage_kerr_ergosphere_prograde.png "$KV" "-0.86,0,0" -1.5708 4
# NEVER place the camera exactly on the axis (x=y=0): the BL chart is
# degenerate at sin(theta)=0 and the rays get mis-launched; the 0.05
# offset is the workaround.
render vantage_kerr_polar.png          "$KV" "0.05,0,-3"   -3.14159 1
render vantage_kerr_gap_rim.png        "$KV" "-0.86,0,0"   -2.4    4
render vantage_kerr_gap_outward.png    "$KV" "-0.86,0,0"   0       4
echo "Done."

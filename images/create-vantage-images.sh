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
  "$BIN" --width=800 --height=800 --max-steps=1000000 \
    --camera-position="$3" --theta="$4" --psi=0 --phi=0 --exposure="$5" \
    --config-file "$2" render --filename="images/$1"
}
# Schwarzschild vantages
render vantage_porthole.png            "$SV" "-1.05,0,0"  0        4
render vantage_photonsphere_edge.png   "$SV" "-1.25,0,0"  -1.5708  4
render vantage_grazing_winding.png     "$SV" "-1.55,0,0"  -2.1     4
render vantage_rearview.png            "$SV" "-2.0,0,0.4" -0.5     4
# Kerr companions (ZAMO)
render vantage_kerr_dragged_porthole.png "$KV" "-0.75,0,0"   0       4
render vantage_kerr_photonsphere_edge.png "$KV" "-0.95,0,0"  -1.75   4
render vantage_kerr_grazing_winding.png "$KV" "-1.12,0,0"   -2.25   4
render vantage_kerr_rearview.png       "$KV" "-2.55,0,0.4" -0.5     4
# Kerr-only vantages
render vantage_kerr_ergosphere_prograde.png "$KV" "-0.86,0,0" -1.5708 4
render vantage_kerr_polar.png          "$KV" "0,0,-1.3"    0.3     0.5
render vantage_kerr_gap_rim.png        "$KV" "-0.86,0,0"   -2.4    4
render vantage_kerr_gap_outward.png    "$KV" "-0.86,0,0"   0       4
echo "Done."

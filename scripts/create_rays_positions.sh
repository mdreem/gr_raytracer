#!/bin/bash

# Create rays by varying the vertical position of the ray's origin, starting from -5 and going up to 5 in increments
# of 0.1.

set -euo pipefail -x

geometry=${1:-}
if [ -z "$geometry" ]; then
  echo "Usage: $0 <geometry> [command]"
  echo "  <geometry> must be one of: schwarzschild, kerr"
  exit 2
fi

case "$geometry" in
  schwarzschild|kerr)
    ;;
  *)
    echo "Invalid geometry: $geometry"
    echo "Allowed: schwarzschild, kerr"
    exit 2
    ;;
esac

shift || true

shift_parameter=${1:-0}

command=${2:-"cargo run --release --"}

run_kerr() {
  local i="$1"
  local formatted
  formatted=$(printf "%07.2f" "$i")
  local pos_y
  pos_y=$(echo "$i - $shift_parameter" | bc -l)
  $command --width 401 --height 401 --max-radius=13 --max-steps=1000000 --theta=1.57 --psi=1.57 --phi=0 --config-file scene-definitions/kerr.toml render-ray-at --position=-5,"${pos_y},0" --direction=1.0,0.0,0.0 --filename "rays/ray-${formatted}.csv" &> /tmp/create_rays_log
}

run_schwarzschild() {
  local i="$1"
  local formatted
  local pos_y
  pos_y=$(echo "$i - $shift_parameter" | bc -l)
  formatted=$(printf "%07.2f" "$i")
  $command --width 401 --height 401 --max-radius=13 --max-steps=1000000 --theta=0.0 --psi=1.57 --phi=1.57 --config-file scene-definitions/schwarzschild.toml render-ray-at --position=-5,"${pos_y},0" --direction=-1.0,0.0,0.0 --filename "rays/ray-${formatted}.csv" &> /tmp/create_rays_log
}

mkdir -p rays
for i in $(seq 0 0.1 10); do
  echo "Creating ray $i"
  case "$geometry" in
    schwarzschild)
      run_schwarzschild "$i"
      ;;
    kerr)
      run_kerr "$i"
      ;;
  esac
done

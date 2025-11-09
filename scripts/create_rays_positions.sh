#!/bin/bash

set -euo pipefail

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

command=${1:-"cargo run --release --"}

run_kerr() {
  local formatted
  formatted=$(printf "%07.2f" "$i")
  $command --width 401 --height 401 --max-radius=13 --max-steps=1000000 --theta=1.57 --psi=1.57 --phi=0 --config-file kerr.toml render-ray-at --position=-5,0,${i} --direction=1.0,0.0,0.0 --filename "rays/ray-${formatted}.csv" &> /tmp/create_rays_log
}

run_schwarzschild() {
  local formatted
  formatted=$(printf "%07.2f" "$i")
  $command --width 401 --height 401 --max-radius=13 --max-steps=1000000 --theta=0.0 --psi=1.57 --phi=1.57 --config-file schwarzschild.toml render-ray-at --position=-5,0,${i} --direction=1.0,0.0,0.0 --filename "rays/ray-${formatted}.csv" &> /tmp/create_rays_log
}

mkdir -p rays
for i in $(seq 0 0.1 10); do
  echo "Creating ray $i"
  case "$geometry" in
    schwarzschild)
      run_schwarzschild
      ;;
    kerr)
      run_kerr
      ;;
  esac
done

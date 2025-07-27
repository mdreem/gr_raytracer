#!/bin/bash

set -euo pipefail

command=${1:-"cargo run --release --"}

mkdir -p rays
for i in $(seq 120 10 280); do
  echo "Creating ray $i"
  $command -w401 -h401 --max-radius=10.0 --config-file schwarzschild.toml --camera-position "0.0,0.0,-18.0" render-ray -r"${i}" -c200 --filename "rays/ray-200-${i}.csv"
done

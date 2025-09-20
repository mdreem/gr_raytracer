#!/bin/bash

set -euo pipefail

command=${1:-"cargo run --release --"}

mkdir -p rays
for i in $(seq 0 0.1 10); do
  echo "Creating ray $i"
  $command -w401 -h401 --max-radius=13 --max-steps=30000 --config-file schwarzschild.toml render-ray-at --position=8.0,0.0,${i} --direction=0.0,1.0,0.0 --filename "rays/ray-200-${i}.csv" &> /dev/null
done

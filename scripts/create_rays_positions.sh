#!/bin/bash

set -euo pipefail

command=$1

mkdir -p rays
for i in $(seq 0 0.1 10); do
  echo "Creating ray $i"
  $command -w401 -h401 --max-steps=30000 --config-file schwarzschild.toml render-ray-at --position=0.0,${i},-10.0 --direction=0.0,1.0,0.0 --filename "rays/ray-200-${i}.csv"
done

#!/bin/bash

set -euo pipefail

command=${1:-"cargo run --release --"}

mkdir -p kerr_images
#for i in $(seq 0 0.025 0.5); do
#  echo "Creating kerr black hole with a=$i"
#  cp kerr.toml kerr_images/kerr_a_${i}.toml
#  sed -i '' "s/^a=.*/a=${i}/" kerr_images/kerr_a_${i}.toml
#  echo "Rendering image for a=$i"
#  $command --width=500 --height=500 --max-steps=1000000 --camera-position=-10,0,0.5 --theta=1.62 --psi=1.57 --phi=0 --config-file kerr_images/kerr_a_${i}.toml render --filename "kerr_images/kerr_a_${i}.png"
#done

for i in $(seq 0 0.025 0.5); do
  echo "kerr_images/kerr_a_${i}.png"
done > kerr_images/image_list.txt

magick -delay 20 -loop 0 $(cat kerr_images/image_list.txt) kerr_animation.gif

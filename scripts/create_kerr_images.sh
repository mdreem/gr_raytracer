#!/bin/bash

set -euo pipefail

command=${1:-"cargo run --release --"}

mkdir -p kerr_images
for i in $(seq 0 0.025 0.5); do
  a_value=$(printf "%01.3f" ${i})
  if [[ "$a_value" == "0.500" ]]; then
    a_value="0.499"
  fi
  echo "Creating kerr black hole with a=$a_value"
  cp scripts/resources/kerr.toml "kerr_images/kerr_a_${a_value}.toml"
  sed -i.bak "s/^a=.*/a=${a_value}/" "kerr_images/kerr_a_${a_value}.toml"
  echo "Rendering image for a=$a_value"
  FILE="kerr_images/kerr_a_${a_value}.png"
  FILE_ANNOTATED="kerr_images/kerr_a_${a_value}_annotated.png"
  if [ -f "$FILE" ]; then
    echo "$FILE exists, skipping rendering."
  else
    $command --width=500 --height=500 --max-steps=1000000 --camera-position=-15,0,1.0 --theta=1.62 --psi=1.57 --phi=0 --config-file "kerr_images/kerr_a_${a_value}.toml" render --filename "$FILE"
  fi
  magick "$FILE" -gravity Northwest -pointsize 30 -annotate +20+20 "a = ${a_value}" "$FILE_ANNOTATED"
done

for i in $(seq 0 0.025 0.5); do
  a_value=$(printf "%01.3f" ${i})
  if [[ "$a_value" == "0.500" ]]; then
    a_value="0.499"
  fi
  echo "kerr_images/kerr_a_${a_value}_annotated.png"
done > kerr_images/image_list.txt

magick -delay 20 -loop 0 $(cat kerr_images/image_list.txt) kerr_animation.gif

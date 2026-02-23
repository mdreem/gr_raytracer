#!/bin/bash

set -euo pipefail

default_command="cargo run --release --"
command="$default_command"

render_width=500
render_height=500

grid_enabled=0
grid_width=50
grid_height=50
grid_line_width=1
grid_color="rgba(0,0,0,0.30)"

usage() {
  cat <<EOF
Usage: $(basename "$0") [command] [options]

Options:
  --grid                 Enable a semi-transparent grid overlay.
  --no-grid              Disable grid overlay (default).
  --grid-size N          Grid cell size in pixels for both width and height.
  --grid-width N         Grid cell width in pixels.
  --grid-height N        Grid cell height in pixels.
  -h, --help             Show this help message.

Examples:
  $(basename "$0")
  $(basename "$0") --grid --grid-size 40
  $(basename "$0") "cargo run --release --" --grid --grid-width 40 --grid-height 30
EOF
}

require_positive_int() {
  local flag="$1"
  local value="$2"
  if ! [[ "$value" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: ${flag} requires a positive integer, got: ${value}" >&2
    exit 1
  fi
}

if [[ $# -gt 0 && "$1" != --* && "$1" != "-h" ]]; then
  command="$1"
  shift
fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    --grid)
      grid_enabled=1
      shift
      ;;
    --no-grid)
      grid_enabled=0
      shift
      ;;
    --grid-size)
      if [[ $# -lt 2 ]]; then
        echo "Error: --grid-size requires a value." >&2
        exit 1
      fi
      require_positive_int "--grid-size" "$2"
      grid_width="$2"
      grid_height="$2"
      shift 2
      ;;
    --grid-width)
      if [[ $# -lt 2 ]]; then
        echo "Error: --grid-width requires a value." >&2
        exit 1
      fi
      require_positive_int "--grid-width" "$2"
      grid_width="$2"
      shift 2
      ;;
    --grid-height)
      if [[ $# -lt 2 ]]; then
        echo "Error: --grid-height requires a value." >&2
        exit 1
      fi
      require_positive_int "--grid-height" "$2"
      grid_height="$2"
      shift 2
      ;;
    -h | --help)
      usage
      exit 0
      ;;
    *)
      echo "Error: unknown argument '$1'." >&2
      usage
      exit 1
      ;;
  esac
done

declare -a grid_draw_args=()
if (( grid_enabled )); then
  for ((x=grid_width; x<render_width; x+=grid_width)); do
    grid_draw_args+=(-draw "line ${x},0 ${x},${render_height}")
  done
  for ((y=grid_height; y<render_height; y+=grid_height)); do
    grid_draw_args+=(-draw "line 0,${y} ${render_width},${y}")
  done
fi

mkdir -p kerr_images
for i in $(seq 0 0.010 0.5); do
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
    $command --width="$render_width" --height="$render_height" --max-steps=1000000 --camera-position=-15,0,1.0 --theta=1.62 --psi=1.57 --phi=0 --config-file "kerr_images/kerr_a_${a_value}.toml" render --filename "$FILE"
  fi
  magick "$FILE" -gravity Northwest -pointsize 30 -annotate +20+20 "a = ${a_value}" "$FILE_ANNOTATED"
  if (( grid_enabled )) && (( ${#grid_draw_args[@]} > 0 )); then
    magick "$FILE_ANNOTATED" \
      \( -size "${render_width}x${render_height}" xc:none \
      -stroke "$grid_color" -strokewidth "$grid_line_width" -fill none \
      "${grid_draw_args[@]}" \) \
      -compose over -composite \
      "$FILE_ANNOTATED"
  fi
done

for i in $(seq 0 0.010 0.5); do
  a_value=$(printf "%01.3f" ${i})
  if [[ "$a_value" == "0.500" ]]; then
    a_value="0.499"
  fi
  echo "kerr_images/kerr_a_${a_value}_annotated.png"
done > kerr_images/image_list.txt

magick -delay 20 -loop 0 $(cat kerr_images/image_list.txt) kerr_animation.gif
echo "Created kerr_animation.gif."
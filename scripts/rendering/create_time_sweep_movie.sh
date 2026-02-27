#!/bin/bash

set -euo pipefail

default_command="cargo run --release --"
command="$default_command"

config_file="scene-definitions/schwarzschild-checker.toml"
output_dir="movie_timesteps"
output_file="time_sweep.mp4"
frame_prefix="render"

width=500
height=500
max_steps=1000000
camera_position="-16,0,3.5"
theta="-3.142"
psi="0.0"
phi="0.0"

time_start=0.0
time_end=20.0
time_step=0.2
fps=30

overwrite=0

usage() {
  cat <<EOF
Usage: $(basename "$0") [command] [options]

Renders a frame sequence with increasing --initial-time and creates an MP4.

Options:
  --config-file PATH       Scene TOML to render.
  --output-dir DIR         Directory for frames + resulting movie.
  --output-file NAME       Output movie filename inside output dir.
  --frame-prefix NAME      Frame filename prefix (default: render).
  --width N                Render width.
  --height N               Render height.
  --max-steps N            Integrator max steps.
  --camera-position X,Y,Z  Camera position.
  --theta VALUE            Camera theta.
  --psi VALUE              Camera psi.
  --phi VALUE              Camera phi.
  --time-start VALUE       First initial-time value (inclusive).
  --time-end VALUE         Last initial-time value (inclusive).
  --time-step VALUE        Increment per frame.
  --fps N                  Movie frame rate.
  --overwrite              Re-render existing frames and overwrite movie.
  -h, --help               Show this help.

Examples:
  $(basename "$0")
  $(basename "$0") --time-start 0 --time-end 40 --time-step 0.5 --fps 24
  $(basename "$0") "target/release/gr_raytracer" --output-dir movie_checker
EOF
}

if [[ $# -gt 0 && "$1" != --* && "$1" != "-h" ]]; then
  command="$1"
  shift
fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config-file)
      config_file="$2"
      shift 2
      ;;
    --output-dir)
      output_dir="$2"
      shift 2
      ;;
    --output-file)
      output_file="$2"
      shift 2
      ;;
    --frame-prefix)
      frame_prefix="$2"
      shift 2
      ;;
    --width)
      width="$2"
      shift 2
      ;;
    --height)
      height="$2"
      shift 2
      ;;
    --max-steps)
      max_steps="$2"
      shift 2
      ;;
    --camera-position)
      camera_position="$2"
      shift 2
      ;;
    --theta)
      theta="$2"
      shift 2
      ;;
    --psi)
      psi="$2"
      shift 2
      ;;
    --phi)
      phi="$2"
      shift 2
      ;;
    --time-start)
      time_start="$2"
      shift 2
      ;;
    --time-end)
      time_end="$2"
      shift 2
      ;;
    --time-step)
      time_step="$2"
      shift 2
      ;;
    --fps)
      fps="$2"
      shift 2
      ;;
    --overwrite)
      overwrite=1
      shift
      ;;
    -h | --help)
      usage
      exit 0
      ;;
    *)
      echo "Error: unknown argument '$1'" >&2
      usage
      exit 1
      ;;
  esac
done

if ! command -v ffmpeg >/dev/null 2>&1; then
  echo "Error: ffmpeg is required but not found in PATH." >&2
  exit 1
fi

if [[ ! -f "$config_file" ]]; then
  echo "Error: config file not found: $config_file" >&2
  exit 1
fi

if ! awk "BEGIN { exit !($time_step > 0) }"; then
  echo "Error: --time-step must be > 0, got $time_step" >&2
  exit 1
fi

if ! awk "BEGIN { exit !($time_end >= $time_start) }"; then
  echo "Error: --time-end must be >= --time-start." >&2
  exit 1
fi

mkdir -p "$output_dir"

frame_manifest="${output_dir}/frame_times.csv"
echo "frame,time,filename" > "$frame_manifest"

index=0
for t in $(seq "$time_start" "$time_step" "$time_end"); do
  time_value=$(printf "%.6f" "$t")
  frame_file=$(printf "%s/%s-%04d.png" "$output_dir" "$frame_prefix" "$index")

  if [[ -f "$frame_file" && "$overwrite" -eq 0 ]]; then
    echo "Skipping existing frame: $frame_file"
  else
    echo "Rendering frame $index with initial_time=$time_value"
    $command \
      --width="$width" \
      --height="$height" \
      --max-steps="$max_steps" \
      --camera-position="$camera_position" \
      --theta="$theta" \
      --psi="$psi" \
      --phi="$phi" \
      --initial-time="$time_value" \
      --config-file "$config_file" \
      render --filename "$frame_file"
  fi

  printf "%d,%s,%s\n" "$index" "$time_value" "$frame_file" >> "$frame_manifest"
  index=$((index + 1))
done

if [[ "$index" -eq 0 ]]; then
  echo "No frames were generated." >&2
  exit 1
fi

movie_path="${output_dir}/${output_file}"
if [[ -f "$movie_path" && "$overwrite" -eq 0 ]]; then
  echo "Error: movie already exists at $movie_path (use --overwrite to replace)." >&2
  exit 1
fi

echo "Creating movie at $movie_path"
ffmpeg_overwrite_flag="-n"
if [[ "$overwrite" -eq 1 ]]; then
  ffmpeg_overwrite_flag="-y"
fi

ffmpeg "$ffmpeg_overwrite_flag" \
  -framerate "$fps" \
  -i "${output_dir}/${frame_prefix}-%04d.png" \
  -c:v libx264 \
  -pix_fmt yuv420p \
  "$movie_path"

echo "Done. Frames: $index"
echo "Frame manifest: $frame_manifest"
echo "Movie: $movie_path"

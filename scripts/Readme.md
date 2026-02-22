# Scripts

## How to create a movie

First, use `create_camera_trajectory.py` to create a camera trajectory. This will generate a file with the camera
positions.

These will serve as inputs to `create_images_from_camera_positions.py`, which will generate images based on the camera
positions. Here also the scene definition is references.

Then, use `create_movie.py` to create the movie. This script just uses ffmpeg to create the movie from the images
generated in the previous step.

## Create a trajectory

```sh
cargo run --release -- --width 401 --height 401 --max-radius=13 --max-steps=1000000 --theta=1.57 --psi=1.57 --phi=0 --config-file kerr.toml render-ray-at --position=-5,-5,0 --direction=1.0,0.0,0.0 --filename "output.csv"
```

## Manim animation

The manim animation expects a set of rays to animate in the `rays/` directory. You can generate these rays using the
`create_rays_positions.sh` or `create_rays_from_camera.sh` scripts.

```sh
python -m manim scripts/animate-rays/main.py AnimateRays
```

# General Relativity Raytracer (Rust)

[![CI](https://github.com/mdreem/gr_raytracer/actions/workflows/ci.yaml/badge.svg)](https://github.com/mdreem/gr_raytracer/actions/workflows/ci.yaml)
[![Release](https://github.com/mdreem/gr_raytracer/actions/workflows/release.yaml/badge.svg)](https://github.com/mdreem/gr_raytracer/actions/workflows/release.yaml)
[![Latest Release](https://img.shields.io/github/v/release/mdreem/gr_raytracer)](https://github.com/mdreem/gr_raytracer/releases)
[![License](https://img.shields.io/github/license/mdreem/gr_raytracer)](https://github.com/mdreem/gr_raytracer/blob/main/LICENSE)
[![Downloads](https://img.shields.io/github/downloads/mdreem/gr_raytracer/total)](https://github.com/mdreem/gr_raytracer/releases)

<p align="center">
  <img src="kerr_black_hole_with_stars.png" alt="Kerr black hole raytraced in Rust with stars in the background" title="Kerr black hole with stars in the background">
</p>

A Rust ray tracer for **general relativity** and **black hole visualization**.  
It solves geodesic equations numerically and renders gravitational lensing, redshift, beaming, accretion disks, and
photon trajectories in Schwarzschild and Kerr spacetimes.

Inspired
by [Seeing relativity -- I. Ray tracing in a Schwarzschild metric to explore the maximal analytic extension of the metric and making a proper rendering of the stars](https://arxiv.org/abs/1511.06025)
and [BlackHoleViz_v2](https://github.com/HollowaySean/BlackHoleViz_v2).

See the [Image Gallery](./images/images.md) for more rendered outputs.

## Table of Contents

- [Features](#features)
- [Quick Start](#quick-start)
- [How Rendering Works](#how-rendering-works)
- [Scripts](#scripts)
- [Examples](#examples)
- [References](#references)
- [License](#license)

## Features

- Multi-geometry ray tracing: `Euclidean`, `EuclideanSpherical`, `Schwarzschild`, and `Kerr`.
- Geodesic integration with RKF45.
- Relativistic shading effects: gravitational and Doppler redshift plus relativistic beaming.
- Physically motivated emission: black-body spectrum integration in CIE XYZ.
- Scene primitives: `Sphere`, `Disc`, and Perlin-noise-based `VolumetricDisc`.
- Flexible materials and textures: bitmap, checker, and black-body mappers.
- Image output formats: standard image export plus HDR (`.hdr`).
- Debug and analysis tools: per-pixel ray export and arbitrary ray-at-position export.
- Config-driven scenes through TOML definitions.

Note: `VolumetricDisc` is primarily a visual effect for accretion-disk appearance rather than a strict physical
model.

## Quick Start

### Requirements

- Rust toolchain (`cargo`)

### Build

```sh
cargo build --release
```

### Render a First Image

```sh
cargo run --release -- --width=1500 --height=1500 --camera-position=-16.0,0.0,3.5 --theta=-3.142 --psi=0.0 --phi=0.0 --config-file scene-definitions/schwarzschild.toml render --filename=render.png
```

- `--width` and `--height`: output resolution.
- `--camera-position`: camera location.
- `render`: full image render command.
- `--filename`: output image path.

## How Rendering Works

Predefined scenes live in [`scene-definitions`](./scene-definitions) as TOML files.  
You can swap geometry, textures, and objects by choosing or editing a scene file.

## Scripts

This repository includes helper scripts for ray export and animation (Manim-based).

### Create rays

- `scripts/create_rays_positions.sh`: Generates rays in Schwarzschild spacetime from a position and direction, then
  writes CSV files to `rays/`.
- `scripts/create_rays_from_camera.sh`: Generates Schwarzschild rays starting from a camera and selected pixel, then
  writes CSV files to `rays/`.

### Plot/animate rays

```sh
python -m manim scripts/animate-rays/main.py AnimateRays
```

This command renders an animation from CSV ray data in `rays/`.

## Examples

### Schwarzschild black hole with accretion disk

Checkerboard texturing helps visualize lensing and warped geometry.

<p align="center">
  <img src="./images/render_schwarzschild_checker_texture.png" alt="Schwarzschild black hole with checkerboard accretion disk showing lensing distortions" title="Schwarzschild black hole with accretion disk">
</p>

### Rays in a Schwarzschild metric (video)

https://github.com/user-attachments/assets/c1ce889b-6186-4ce5-b613-fecae3af03ef

### Flyover of a Schwarzschild black hole (video)

https://github.com/user-attachments/assets/914d3134-53db-4084-8a5f-1728d8460594

Background image source: https://commons.wikimedia.org/wiki/File:Messier_object_025.jpg

### Lensing with background object (video)

A Schwarzschild black hole with a spherical object behind it, showing lensing while the camera moves.

https://github.com/user-attachments/assets/6907c6a2-970a-4d19-be60-5e0f6f340709

Background image source: https://commons.wikimedia.org/wiki/File:Messier_object_025.jpg

### Kerr black hole with accretion disk

Example render command:

```sh
gr_raytracer --width=500 --height=500 --max-steps=1000000 --camera-position=-10,0,-0.5 --theta=1.52 --psi=-1.57 --phi=0 --config-file scene-definitions/kerr.toml render
```

Kerr scenes often require a high `--max-steps` value because of complex geodesic behavior near the black hole.

#### Example 1: `r_s = 1.0`, `a = 0.5`

<p align="center">
  <img src="./images/render_kerr_checker_texture.png" alt="Kerr black hole render with spin parameter a equals 0.5" title="Kerr black hole with r_s = 1.0 and a = 0.5">
</p>

##### Near-horizon trajectory

```sh
gr_raytracer --width=501 --height=501 --max-steps=1000000 --camera-position=-5,0,0.5 --theta=1.57 --psi=1.57 --phi=0 --config-file scene-definitions/kerr.toml render-ray --col=195 --row=250
```

<p align="center">
  <img src="./images/kerr_trajectory_near_horizon.png" alt="Photon trajectory near the horizon of a Kerr black hole" title="Trajectory near horizon for Kerr black hole with r_s = 1.0 and a = 0.5">
</p>

#### Example 2: `r_s = 1.0`, `a = 0.51`

<p align="center">
  <img src="./images/render_kerr_large_a_checker_texture.png" alt="Kerr black hole render with higher spin parameter a equals 0.51" title="Kerr black hole with r_s = 1.0 and a = 0.51">
</p>

### Kerr spin animation

Animation of a Kerr black hole with `r_s = 1.0` and spin parameter `a` increasing from `0.0` to `0.5`.

<p align="center">
  <img src="./images/kerr_animation.gif" alt="Animation of Kerr black hole as spin parameter increases from 0 to 0.5" title="Kerr black hole animation with increasing spin">
</p>

## References

- [Seeing relativity -- I. Ray tracing in a Schwarzschild metric to explore the maximal analytic extension of the metric and making a proper rendering of the stars](https://arxiv.org/abs/1511.06025)
- [BlackHoleViz_v2](https://github.com/HollowaySean/BlackHoleViz_v2)
- Novikov, I. D., & Thorne, K. S. (1973). *Astrophysics of black holes*. In C. DeWitt & B. S. DeWitt (Eds.), *Black
  Holes (Les Astres Occlus)*, p.
    343. [Chapter bibliographic entry](https://cir.nii.ac.jp/crid/1370025430666224928), [Book record](https://lccn.loc.gov/73169355)
- https://commons.wikimedia.org/wiki/File:Messier_object_025.jpg

## License

This project is licensed under the terms in [`LICENSE`](./LICENSE).

# General Relativity Raytracer

This is a raytracer for general relativity, which can be used to visualize the effects of gravity on light paths. It is
based on the principles of general relativity and uses numerical methods to solve the geodesic equations.

It is based on the
paper [Seeing relativity -- I. Ray tracing in a Schwarzschild metric to explore the maximal analytic extension of the metric and making a proper rendering of the stars](https://arxiv.org/abs/1511.06025).

# Scripts

There are various scripts. Some of them create importable CSV files, others create images and animations based
on [Manim](https://github.com/3b1b/manim).

## Create rays to be plotted.

- `scripts/create_rays_positions.py`: Creates rays in a Schwarzschild metric based on a given position and direction and
  saves them to a CSV file in the directory `rays/`.
- `scripts/create_rays_from_camera.py`: Creates rays in a Schwarzschild metric using the camera given its position and
  a selected pixel. The data will be saved to a CSV file in the directory `rays/`.

## Plot rays.

Running `python -m manim scripts/animate-rays/main.py AnimateRays` will create an animation of the rays saved in
CSV files in the directory `rays/`.

# Example

Plot of the Schwarzschild solution with a accretion disk using a checkerboard texture to visualize the relations.
![alt text](./images/render_checker_texture.png "Black Hole with accretion disk")

# Scripts

## How to create a movie

First, use `create_camera_trajectory.py` to create a camera trajectory. This will generate a file with the camera
positions.

These will serve as inputs to `create_images_from_camera_positions.py`, which will generate images based on the camera
positions. Here also the scene definition is references.

Then, use `create_movie.py` to create the movie. This script just uses ffmpeg to create the movie from the images
generated in the previous step.

#!/bin/bash
directory="movie"
file_pattern="render-%03d.png"
output_file="output.mp4"
fps=30

ffmpeg -framerate "${fps}" -i "${directory}/${file_pattern}" -c:v libx264 -pix_fmt yuv420p "${directory}/${output_file}"

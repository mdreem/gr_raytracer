#!/usr/bin/env -S uv run python
import numpy as np
from scipy.interpolate import make_splprep
import pandas as pd

OUTPUT_FILE = "camera_positions.csv"

# Alternative trajectory:
# 1. Starts at [15, 0, 0] (equatorial plane)
# 2. Goes below the disc a little bit
# 3. Moves above the disc until it's almost above the "north pole"
# 4. Moves down and back below the accretion disc
# 5. Backs off to the original state
points = np.array(
    [
        [15.0, 0.0, 0.0],    # Start: Edge of disc
        [14.0, 2.0, -3.0],   # Below the disc, slightly "right" (y=2)
        [10.0, -2.0, 4.0],   # Rising up, swinging "left" (y=-2)
        [2.0, 3.0, 14.0],    # High above north pole, swinging "right" (y=3)
        [-5.0, -3.0, 10.0],  # Coming down, swinging "left" (y=-3)
        [-12.0, 2.0, -3.0],  # Below disc again, swinging "right" (y=2)
        [5.0, -5.0, -1.0],   # Swinging back around, "left" (y=-5)
        [15.0, 0.0, 0.0],    # Home
    ]
)

x, y, z = points.T

# Create a periodic spline if the first and last points are the same
# but here we just use s=0 for exact interpolation.
spl, u = make_splprep([x, y, z], s=0)

N = 240 # Increased number of frames for a longer path
u_new = np.linspace(0, 1, N)
x_new, y_new, z_new = spl(u_new)

positions = np.vstack([x_new, y_new, z_new]).T

df = pd.DataFrame(positions, columns=["x", "y", "z"])
df.to_csv(OUTPUT_FILE, index=False)

print(f"Saved {N} positions to {OUTPUT_FILE}")

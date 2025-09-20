#!/usr/bin/env python3
import numpy as np
from scipy.interpolate import make_splprep
import pandas as pd

OUTPUT_FILE = "camera_positions.csv"

points = np.array(
    [
        [15, 0, 0],
        [13, 1, 1.5],
        [9, -1, -0.5],
        [6, 0, -2],
    ]
)

x, y, z = points.T

spl, u = make_splprep([x, y, z], s=0)

N = 120
u_new = np.linspace(0, 1, N)
x_new, y_new, z_new = spl(u_new)

positions = np.vstack([x_new, y_new, z_new]).T

df = pd.DataFrame(positions, columns=["x", "y", "z"])
df.to_csv(OUTPUT_FILE, index=False)

print(f"Saved positions to {OUTPUT_FILE}")

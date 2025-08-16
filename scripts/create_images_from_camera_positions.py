#!/usr/bin/env python3
import pandas as pd
import subprocess

WIDTH = 500
HEIGHT = 500
CONFIG_FILE = "schwarzschild.toml"
POSITIONS_FILE = "camera_positions.csv"

# takes a csv created using the script create_camera_trajectory.py
df = pd.read_csv(POSITIONS_FILE)

for idx, row in df.iterrows():
    x, y, z = row["x"], row["y"], row["z"]
    position_str = f"{x},{y},{z}"
    idx_str = f"{idx:03d}"

    res = subprocess.run(
        [
            "cargo",
            "run",
            "--release",
            "--",
            f"--width={WIDTH}",
            f"--height={HEIGHT}",
            f"--camera-position={position_str}",
            f"--config-file={CONFIG_FILE}",
            "render",
            f"--filename=movie/render-{idx_str}.png",
        ]
    )
    if res.returncode != 0:
        print(f"Error rendering image for position {position_str}")
        break

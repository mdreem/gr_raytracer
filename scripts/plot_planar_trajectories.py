#!/usr/bin/env -S uv run python
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

SCHWARZSCHILD_RADIUS = 1.0


def read_trajectories(file_paths):
    trajectories = []
    for file_path in file_paths:
        if os.path.isfile(file_path) and file_path.endswith(".csv"):
            df = pd.read_csv(file_path)

            if "r" in df.columns and "phi" in df.columns:
                r = df["r"].values
                phi = df["phi"].values  # phi is already in radians

                # Convert to Cartesian coordinates
                x = r * np.cos(phi)
                y = r * np.sin(phi)

                trajectories.append(
                    (x, y, os.path.basename(file_path))
                )  # Store with filename for labeling
            elif "x" in df.columns and "y" in df.columns:
                x = df["x"].values
                y = df["z"].values

                # Store Cartesian coordinates directly
                trajectories.append((x, y, os.path.basename(file_path)))
            else:
                print(f"Skipping {file_path}: Missing required columns 'r' and 'phi'")
        else:
            print(f"Skipping {file_path}: Not a valid CSV file")
    return trajectories


def plot_trajectories(trajectories, disk_radius):
    plt.figure(figsize=(8, 8))

    # Plot trajectories
    for x, y, label in trajectories:
        plt.plot(x, y, label=label)

    # Add central disk
    disk = plt.Circle(
        (0, 0),
        disk_radius,
        color="gray",
        alpha=0.5,
        label=f"Black Hole (r={disk_radius})",
    )
    plt.gca().add_patch(disk)

    plt.xlabel("X")
    plt.ylabel("Y")
    plt.title("Trajectories in Cartesian Coordinates with Central Disk")
    plt.legend()
    plt.axis("equal")  # Keep aspect ratio square
    plt.grid()
    plt.show()


# Main execution
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: uv run python script.py file1.csv file2.csv ...")
        sys.exit(1)

    file_paths = sys.argv[1:]
    trajectories = read_trajectories(file_paths)

    if trajectories:
        plot_trajectories(trajectories, SCHWARZSCHILD_RADIUS)
    else:
        print("No valid trajectories found.")

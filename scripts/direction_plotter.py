import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

csv_file = "rays.csv"
data = pd.read_csv(csv_file)

x1 = data["d_x1"]
x2 = data["d_x2"]
x3 = data["d_x3"]

mx1 = data["m_x1"]
mx2 = data["m_x2"]
mx3 = data["m_x3"]


fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")

ax.scatter(x1, x2, x3, color="r", label="Direction")
ax.scatter(mx1, mx2, mx3, color="g", label="Momentum")

basis_vectors = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]  # i-hat  # j-hat  # k-hat
colors = ["b", "g", "c"]
labels = ["X", "Y", "Z"]

for vec, color, label in zip(basis_vectors, colors, labels):
    ax.quiver(
        0,
        0,
        0,
        vec[0],
        vec[1],
        vec[2],
        color=color,
        arrow_length_ratio=0.2,
        label=label,
    )

ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
ax.legend()

plt.show()

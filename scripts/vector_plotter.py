import numpy as np
import matplotlib.pyplot as plt

# Prepare a 3D figure
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')

# ---------------------------
# 1. Plot a basic coordinate system
# ---------------------------
# You can think of these as the unit vectors i, j, k at the origin (0,0,0).
ax.quiver(0, 0, 0, 1, 0, 0, color='red', arrow_length_ratio=0.1, linewidth=2, label='x-axis')
ax.quiver(0, 0, 0, 0, 1, 0, color='green', arrow_length_ratio=0.1, linewidth=2, label='y-axis')
ax.quiver(0, 0, 0, 0, 0, 1, color='blue', arrow_length_ratio=0.1, linewidth=2, label='z-axis')

# ---------------------------
# 2. Define your own vectors
# ---------------------------
# Suppose we have N vectors. Each vector has:
#   - A "start" (position) in 3D space, (px, py, pz)
#   - A "direction" (vx, vy, vz)
positions = np.array([
    [0, 0, -5],
    [0, 0, -5],
    [0, 0, -5]
])
directions = np.array([
    [-0.003999968000255998, -0.003999968000255998, 0.9999840001279994],
    [-0.6662216323892758, -0.6662216323892758, -0.3351081513081079],
    [0.6679945915286996, 0.6679945915286996, -0.3279732479590546]
])

# Plot each vector at its position, pointing in the given direction
for pos, dir_ in zip(positions, directions):
    ax.quiver(
        pos[0], pos[1], pos[2],  # starting point of the vector
        dir_[0], dir_[1], dir_[2],  # direction of the vector
        color='black',
        arrow_length_ratio=0.2,
        linewidth=1.5
    )

# Plot sphere

r = 2  # radius
phi, theta = np.mgrid[0 : np.pi : 50j, 0 : 2*np.pi : 50j]
x = r * np.sin(phi) * np.cos(theta)
y = r * np.sin(phi) * np.sin(theta)
z = r * np.cos(phi)

# Plot the sphere as a surface
ax.plot_surface(x, y, z, color='cyan', alpha=0.3, linewidth=0)


# ---------------------------
# 3. Make it look nicer
# ---------------------------
max_range = 10
ax.set_xlim([-max_range, max_range])
ax.set_ylim([-max_range, max_range])
ax.set_zlim([-max_range, max_range])
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.legend()

plt.title("3D Vectors in Matplotlib")
plt.show()

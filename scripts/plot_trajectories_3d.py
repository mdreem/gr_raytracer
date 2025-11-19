import argparse
import pandas as pd
import matplotlib.pyplot as plt

def set_axes_equal(ax, xs, ys, zs):
    xr = xs.max() - xs.min()
    yr = ys.max() - ys.min()
    zr = zs.max() - zs.min()
    m = max(xr, yr, zr)
    xm = (xs.max() + xs.min()) / 2
    ym = (ys.max() + ys.min()) / 2
    zm = (zs.max() + zs.min()) / 2
    ax.set_xlim(xm - m/2, xm + m/2)
    ax.set_ylim(ym - m/2, ym + m/2)
    ax.set_zlim(zm - m/2, zm + m/2)

def main():
    p = argparse.ArgumentParser(description="Plot multiple 3D trajectories from CSVs with x,y,z columns.")
    p.add_argument("csv", nargs="+", help="Path(s) to CSV file(s)")
    p.add_argument("--save", default=None, help="Optional path to save the figure (e.g., plot.png)")
    args = p.parse_args()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Lists to compute global axis bounds
    all_x, all_y, all_z = [], [], []

    for path in args.csv:
        df = pd.read_csv(path)[['x', 'y', 'z']].dropna()

        x, y, z = df['x'].values, df['y'].values, df['z'].values

        ax.plot(x, y, z, linewidth=1.5, label=f"{path}")
        ax.scatter(x[0], y[0], z[0], s=30)        # start
        ax.scatter(x[-1], y[-1], z[-1], s=30)    # end

        all_x.append(x)
        all_y.append(y)
        all_z.append(z)

    # Compute global axis bounds
    import numpy as np
    all_x = np.concatenate(all_x)
    all_y = np.concatenate(all_y)
    all_z = np.concatenate(all_z)
    set_axes_equal(ax, all_x, all_y, all_z)

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.legend()
    plt.tight_layout()

    if args.save:
        plt.savefig(args.save, dpi=200, bbox_inches='tight')
    else:
        plt.show()

if __name__ == "__main__":
    main()

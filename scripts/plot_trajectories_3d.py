import argparse
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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
    p = argparse.ArgumentParser(description="Plot a 3D trajectory from a CSV with x,y,z columns.")
    p.add_argument("csv", help="Path to CSV file")
    p.add_argument("--save", default=None, help="Optional path to save the figure (e.g., plot.png)")
    args = p.parse_args()

    df = pd.read_csv(args.csv)[['x', 'y', 'z']].dropna()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(df['x'].values, df['y'].values, df['z'].values, linewidth=1.5)
    ax.scatter(df['x'].iloc[0], df['y'].iloc[0], df['z'].iloc[0], s=30, label='start')
    ax.scatter(df['x'].iloc[-1], df['y'].iloc[-1], df['z'].iloc[-1], s=30, label='end')
    ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z'); ax.legend()
    set_axes_equal(ax, df['x'].values, df['y'].values, df['z'].values)
    plt.tight_layout()

    if args.save:
        plt.savefig(args.save, dpi=200, bbox_inches='tight')
    else:
        plt.show()

if __name__ == "__main__":
    main()

#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy", "pillow"]
# ///
"""Paint an accretion-disc texture for the thin `Disc` object.

A warm, art-directed disc in the style of the one in Interstellar: cream-white
at the inner edge falling to deep orange outside, with filaments sheared by
differential rotation. Alpha carries surface coverage, so the gaps let the far
side of the disc and the sky through instead of the disc reading as a painted
ring.

The `Disc` object maps its annulus to UV as a top-down polar plot: the image
centre is the inner edge, radius 0.5 in UV is the outer edge, and the angle is
phi. The texture is therefore authored directly in polar coordinates.

A bitmap-textured disc ignores the Novikov-Thorne temperature profile and uses
a constant temperature, so the radial brightness falloff has to be painted in
here rather than left to the renderer.

Usage:  uv run scripts/textures/make_interstellar_disc.py [--output PATH]
"""

import argparse

import numpy as np
from PIL import Image

SIZE = 2048
NR, NTH = 1024, 2048  # polar noise grid; the theta axis wraps
SHEAR = 2.6  # how far the inner material has wound ahead of the outer

INNER_RGB = np.array([1.00, 0.97, 0.90], dtype=np.float32)
OUTER_RGB = np.array([0.95, 0.56, 0.20], dtype=np.float32)
HOT_RGB = np.array([1.00, 0.98, 0.94], dtype=np.float32)


def octave(rng, nr, nth):
    """One smooth random field on the (r, theta) grid, periodic in theta."""
    coarse = rng.random((nr, nth)).astype(np.float32)
    # Wrap one column from each end before resampling so the seam is seamless.
    tiled = np.concatenate([coarse[:, -1:], coarse, coarse[:, :1]], axis=1)
    pad = NTH // nth
    wide = Image.fromarray((tiled * 255).astype(np.uint8)).resize(
        (NTH + 2 * pad, NR), Image.BICUBIC
    )
    return np.asarray(wide, dtype=np.float32)[:, pad : pad + NTH] / 255.0


def fbm(rng, levels=7):
    out = np.zeros((NR, NTH), dtype=np.float32)
    amp, total = 1.0, 0.0
    for i in range(levels):
        out += amp * octave(rng, 4 * 2**i, 8 * 2**i)
        total += amp
        amp *= 0.64
    return out / total


def build(seed=7):
    rng = np.random.default_rng(seed)
    filaments = fbm(rng)
    # A radially stretched second field so the streaks read as long arcs
    # rather than as isotropic blobs.
    streaks = fbm(rng, levels=4)
    streaks = np.repeat(streaks[::4], 4, axis=0)[:NR]

    yy, xx = np.mgrid[0:SIZE, 0:SIZE].astype(np.float32)
    centre = (SIZE - 1) / 2.0
    dx, dy = xx - centre, yy - centre
    rho = np.sqrt(dx * dx + dy * dy) / (SIZE / 2.0)  # 0 inner edge, 1 outer
    theta = np.arctan2(dy, dx)

    inside = rho <= 1.0
    r = np.clip(rho, 1e-3, 1.0)

    theta_sheared = theta + SHEAR / (r + 0.55)
    ti = np.mod(theta_sheared / (2 * np.pi) * NTH, NTH)
    ri = np.clip(r * (NR - 1), 0, NR - 1)

    def sample(field):
        r0 = np.floor(ri).astype(np.int32)
        t0 = np.floor(ti).astype(np.int32)
        fr, ft = ri - r0, ti - t0
        r1 = np.minimum(r0 + 1, NR - 1)
        t1 = (t0 + 1) % NTH
        return (
            field[r0, t0] * (1 - fr) * (1 - ft)
            + field[r1, t0] * fr * (1 - ft)
            + field[r0, t1] * (1 - fr) * ft
            + field[r1, t1] * fr * ft
        )

    n = 0.62 * sample(filaments) + 0.38 * sample(streaks)
    n = (n - n[inside].min()) / (n[inside].max() - n[inside].min())
    # Ridged transform: thin bright threads instead of broad blobs.
    threads = 1.0 - np.abs(2.0 * n - 1.0)
    n = np.clip(0.55 * n + 0.45 * threads**1.6, 0.0, 1.0)

    radial = 1.0 / (0.45 + 0.75 * r) ** 0.9
    radial /= radial.max()
    fil = np.clip((n - 0.26) / 0.54, 0.0, 1.0) ** 1.15
    brightness = np.clip(radial * (0.50 + 0.70 * fil), 0.10, 1.0)

    mix = np.clip(r**0.55, 0.0, 1.0)[..., None]
    rgb = INNER_RGB * (1 - mix) + OUTER_RGB * mix
    # Hot filaments wash toward white; the gaps between them stay orange.
    rgb = rgb * (1 - 0.30 * fil[..., None]) + HOT_RGB * (0.30 * fil[..., None])
    rgb = np.clip(rgb * brightness[..., None], 0.0, 1.0)

    edge_in = np.clip(r / 0.06, 0.0, 1.0)
    edge_out = np.clip((1.0 - r) / 0.22, 0.0, 1.0)
    alpha = np.clip(0.30 + 0.95 * fil, 0.0, 1.0) * edge_in * edge_out

    alpha = np.where(inside, alpha, 0.0)
    rgb = np.where(inside[..., None], rgb, 0.0)
    return np.dstack([rgb, alpha[..., None]]), alpha[inside].mean()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="resources/interstellar-disc.png")
    parser.add_argument("--seed", type=int, default=7)
    args = parser.parse_args()

    image, mean_alpha = build(args.seed)
    Image.fromarray((image * 255).astype(np.uint8), mode="RGBA").save(args.output)
    print(f"wrote {args.output} ({SIZE}x{SIZE}, mean alpha {mean_alpha:.3f})")


if __name__ == "__main__":
    main()

# Generate synthetic diagnostic textures for the Kerr showcase scenes.
# Run with: uv run --with pillow python scripts/generate_showcase_textures.py
#
# Ring-and-spoke textures for the flat Disc UV map (polar map centered at
# (0.5, 0.5), u = 0.5 + 0.5*r_norm*cos(phi) where r_norm normalizes the
# disc's own annulus [inner, outer] to [0, 1]). The rings identify which
# disc radius each lensed image order comes from; the spokes reveal
# azimuthal shear between image orders.
#
# The rings are coded by ABSOLUTE radius (bands on a fixed global grid), so
# the same color means the same r across scenes with different inner radii;
# each output bitmap is pre-warped for its scene's annulus. This is what
# makes the a = 0.499 vs a = 0 comparison renders comparable band-by-band.

import colorsys
import math

from PIL import Image

SIZE = 1024
N_SPOKES = 12
SPOKE_HALF_WIDTH_RAD = 0.02

# Global band grid in absolute disc radius (units of r_s): shared by all
# scenes so colors align across renders.
R_MIN = 0.62
R_MAX = 8.0
N_BANDS = 8
BAND_WIDTH = (R_MAX - R_MIN) / N_BANDS
RING_GAP_FRACTION = 0.12  # dark separator inside each band

# (filename, scene annulus inner, outer)
VARIANTS = [
    ("resources/ringspoke.png", 0.62, 8.0),
    ("resources/ringspoke_a0.png", 3.05, 8.0),
]


def ring_color(band: int) -> tuple[int, int, int]:
    # Rainbow from red (innermost global band) to violet (outermost),
    # saturated and bright so colors survive tone mapping.
    hue = 0.83 * band / max(N_BANDS - 1, 1)
    r, g, b = colorsys.hsv_to_rgb(hue, 0.9, 1.0)
    return int(255 * r), int(255 * g), int(255 * b)


def generate(path: str, inner: float, outer: float) -> None:
    img = Image.new("RGBA", (SIZE, SIZE), (0, 0, 0, 0))
    px = img.load()
    center = (SIZE - 1) / 2.0
    for y in range(SIZE):
        for x in range(SIZE):
            dx = (x - center) / (SIZE / 2.0)
            dy = (y - center) / (SIZE / 2.0)
            r_norm = math.hypot(dx, dy)
            if r_norm > 1.0:
                continue
            r_abs = inner + r_norm * (outer - inner)
            phi = math.atan2(dy, dx)
            band_pos = (r_abs - R_MIN) / BAND_WIDTH
            band = min(max(int(band_pos), 0), N_BANDS - 1)
            in_gap = (band_pos - band) > (1.0 - RING_GAP_FRACTION)
            spoke_phase = (phi % (2.0 * math.pi / N_SPOKES)) - (
                math.pi / N_SPOKES
            )
            in_spoke = abs(spoke_phase) < SPOKE_HALF_WIDTH_RAD
            if in_gap or in_spoke:
                px[x, y] = (10, 10, 10, 255)
            else:
                px[x, y] = (*ring_color(band), 255)
    img.save(path)
    print(f"wrote {path} (annulus {inner}-{outer})")


def main() -> None:
    for path, inner, outer in VARIANTS:
        generate(path, inner, outer)


if __name__ == "__main__":
    main()

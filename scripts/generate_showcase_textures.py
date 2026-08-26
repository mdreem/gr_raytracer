# Generate synthetic diagnostic textures for the Kerr showcase scenes.
# Run with: uv run --with pillow python scripts/generate_showcase_textures.py
#
# ringspoke.png: concentric rainbow rings crossed by dark radial spokes,
# laid out for the flat Disc UV map (polar map centered at (0.5, 0.5),
# u = 0.5 + 0.5*r*cos(phi), v = 0.5 + 0.5*r*sin(phi), r normalized so the
# outer edge sits at image radius 0.5). The rings identify which disc
# radius each lensed image order comes from; the spokes reveal azimuthal
# shear between image orders.

import colorsys
import math

from PIL import Image

SIZE = 1024
N_RINGS = 8
N_SPOKES = 12
SPOKE_HALF_WIDTH_RAD = 0.02
RING_GAP_FRACTION = 0.12  # dark separator inside each ring band


def ring_color(band: int) -> tuple[int, int, int]:
    # Rainbow from red (inner) to violet (outer), saturated and bright so
    # colors survive tone mapping.
    hue = 0.83 * band / max(N_RINGS - 1, 1)
    r, g, b = colorsys.hsv_to_rgb(hue, 0.9, 1.0)
    return int(255 * r), int(255 * g), int(255 * b)


def main() -> None:
    img = Image.new("RGBA", (SIZE, SIZE), (0, 0, 0, 0))
    px = img.load()
    center = (SIZE - 1) / 2.0
    for y in range(SIZE):
        for x in range(SIZE):
            dx = (x - center) / (SIZE / 2.0)
            dy = (y - center) / (SIZE / 2.0)
            r = math.hypot(dx, dy)
            if r > 1.0:
                continue
            phi = math.atan2(dy, dx)
            band_pos = r * N_RINGS
            band = min(int(band_pos), N_RINGS - 1)
            in_gap = (band_pos - band) > (1.0 - RING_GAP_FRACTION)
            spoke_phase = (phi % (2.0 * math.pi / N_SPOKES)) - (
                math.pi / N_SPOKES
            )
            in_spoke = abs(spoke_phase) < SPOKE_HALF_WIDTH_RAD
            if in_gap or in_spoke:
                px[x, y] = (10, 10, 10, 255)
            else:
                px[x, y] = (*ring_color(band), 255)
    img.save("resources/ringspoke.png")
    print("wrote resources/ringspoke.png")


if __name__ == "__main__":
    main()

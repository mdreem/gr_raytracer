#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy", "pillow"]
# ///
"""Grade a linear HDR render to an 8-bit PNG: exposure, bloom, ACES, gamma.

The renderer's in-engine tone mapping (src/rendering/color.rs) is Reinhard or
GlobalLinear and has no bloom. This is the external grade used for the star and
accretion-disc gallery images: render to a linear .hdr (with the in-engine tone
map effectively bypassed, i.e. GlobalLinear or a high exposure), then run it
through here for the glow and the filmic curve.

Pipeline:
  1. exposure  : multiply linear radiance by the exposure.
  2. bloom     : threshold the luminance, keep the excess, blur it at several
                 scales (tight radii glow point stars, broad radii glow the
                 disc), and add a fraction back. Radii scale with width so the
                 glow is the same fraction of the frame at any resolution.
  3. ACES      : the Narkowicz filmic curve applied to luminance, with the RGB
                 scaled by Lt/L so chromaticity is preserved (a hot pixel keeps
                 its colour instead of clipping to white).
  4. gamma     : encode with 1/2.2 and quantise to 8-bit.

Usage:
  uv run scripts/grade.py in.hdr out.png [--exposure 86] [--bloom 0.7]
                          [--threshold 0.12] [--no-bloom]

Accepts a Radiance .hdr (what the renderer writes) or a .pfm.
"""
import argparse

import numpy as np
from PIL import Image

# Bloom scales (Gaussian sigma, weight), tuned at 1600 px wide. The tight radii
# glow point-source stars; the broad radii glow the accretion disc.
BLOOM_SCALES = [(1, 0.5), (3, 0.5), (9, 0.6), (24, 0.8), (60, 0.6)]
LUMA = np.array([0.2126, 0.7152, 0.0722])


def _read_header_line(f):
    line = f.readline()
    while line.startswith(b"#"):
        line = f.readline()
    return line


def read_pfm(path):
    with open(path, "rb") as f:
        assert _read_header_line(f).strip() == b"PF", "not a colour PFM"
        w, h = map(int, _read_header_line(f).split())
        scale = float(_read_header_line(f))
        data = np.fromfile(f, dtype="<f4" if scale < 0 else ">f4", count=w * h * 3)
    return data.reshape(h, w, 3)[::-1].astype(np.float64)


def read_hdr(path):
    """Decode a Radiance RGBE (.hdr) image to a linear float array."""
    with open(path, "rb") as f:
        if not f.readline().startswith(b"#?"):
            raise ValueError("not a Radiance HDR file")
        while f.readline().strip():  # skip header until the blank line
            pass
        res = f.readline().split()
        if len(res) != 4 or res[0] != b"-Y" or res[2] != b"+X":
            raise ValueError("unsupported HDR orientation: %r" % res)
        h, w = int(res[1]), int(res[3])
        rgbe = np.zeros((h, w, 4), dtype=np.uint8)
        for y in range(h):
            head = f.read(4)
            if head[0] == 2 and head[1] == 2 and ((head[2] << 8) | head[3]) == w:
                for c in range(4):  # new-format run-length encoding, per channel
                    x = 0
                    while x < w:
                        n = f.read(1)[0]
                        if n > 128:  # a run of (n - 128) copies
                            rgbe[y, x:x + n - 128, c] = f.read(1)[0]
                            x += n - 128
                        else:  # n literal bytes
                            rgbe[y, x:x + n, c] = np.frombuffer(f.read(n), np.uint8)
                            x += n
            else:  # flat scanline; the 4 bytes we read are its first pixel
                flat = head + f.read(4 * w - 4)
                rgbe[y] = np.frombuffer(flat, np.uint8).reshape(w, 4)
    exponent = rgbe[..., 3].astype(np.int32)
    scale = np.where(exponent > 0, np.ldexp(1.0, exponent - (128 + 8)), 0.0)
    return rgbe[..., :3].astype(np.float64) * scale[..., None]


def read_image(path):
    img = read_pfm(path) if path.lower().endswith(".pfm") else read_hdr(path)
    return np.maximum(np.nan_to_num(img, posinf=0.0, neginf=0.0), 0.0)


def gaussian_blur(a, sigma):
    """Gaussian blur via an FFT, zero-padded to twice the size so the circular
    convolution does not wrap a bright source across the opposite edge (which
    otherwise shows up as horizontal/vertical streaks off a very bright core)."""
    h, w = a.shape[:2]
    ph, pw = 2 * h, 2 * w
    fy = np.fft.fftfreq(ph)[:, None]
    fx = np.fft.fftfreq(pw)[None, :]
    kernel = np.exp(-2 * np.pi ** 2 * sigma ** 2 * (fx ** 2 + fy ** 2))
    out = np.empty_like(a)
    for c in range(3):
        padded = np.zeros((ph, pw))
        padded[:h, :w] = a[..., c]
        out[..., c] = np.fft.ifft2(np.fft.fft2(padded) * kernel).real[:h, :w]
    return out


def aces(x):
    """Narkowicz ACES filmic approximation, clamped to [0, 1]."""
    a, b, c, d, e = 2.51, 0.03, 2.43, 0.59, 0.14
    return np.clip((x * (a * x + b)) / (x * (c * x + d) + e), 0.0, 1.0)


# Cap the bloom source at this multiple of white, so a single very bright core
# (a hot disc pixel can be 10^6x the sky) cannot dump unbounded energy into the
# blur and swamp the frame. Highlights still bloom, just not without limit.
BLOOM_SOURCE_CAP = 8.0


def grade(img, exposure, bloom_strength, threshold, white=None, tonemap="aces"):
    # Normalize to a white point BEFORE anything else, so the bloom threshold and
    # the ACES knee are relative to the scene's own brightness, not its absolute
    # radiance scale. A raw disc render peaks in the millions; the gallery HDRs
    # were near unity. Auto white = a high luminance percentile (the disc level,
    # not the single brightest firefly). Exposure is then the user's knob on top.
    lum = (img * LUMA).sum(2)
    if white is None:
        white = float(np.percentile(lum, 99.0))
    if not (white > 0.0):
        white = 1.0
    img = img / white * exposure

    if bloom_strength > 0.0:
        lum = (img * LUMA).sum(2, keepdims=True)
        source = np.minimum(img, BLOOM_SOURCE_CAP)
        bright = source * np.clip((lum - threshold) / threshold, 0.0, 1.0)
        s = img.shape[1] / 1600.0
        bloom = sum(w * gaussian_blur(bright, sigma * s) for sigma, w in BLOOM_SCALES)
        img = img + bloom_strength * bloom

    lum = (img * LUMA).sum(2)
    # Tone-map luminance through the chosen curve, scale RGB by the same factor
    # so chromaticity is untouched (hue-preserving). Reinhard L/(1+L) has a
    # gentler toe and keeps saturation; ACES is the filmic look with a harder
    # toe (crushes shadows) and highlight desaturation from the post-scale clip.
    curve = lum / (1.0 + lum) if tonemap == "reinhard" else aces(lum)
    scale = np.divide(curve, lum, out=np.zeros_like(lum), where=lum > 1e-9)
    rgb = np.clip(img * scale[..., None], 0.0, 1.0) ** (1 / 2.2)
    return (rgb * 255.0).astype(np.uint8)


def main():
    ap = argparse.ArgumentParser(description="Grade a linear HDR to PNG (bloom + ACES).")
    ap.add_argument("input", help="linear .hdr (Radiance) or .pfm")
    ap.add_argument("output", help="output .png")
    ap.add_argument("--exposure", type=float, default=1.0,
                    help="brightness on top of the white-point normalization (1 puts the disc near white)")
    ap.add_argument("--white", type=float, default=None,
                    help="white point in linear units (default: auto, the 99th luminance percentile)")
    ap.add_argument("--bloom", type=float, default=0.7, help="bloom strength (0 disables)")
    ap.add_argument("--threshold", type=float, default=0.12, help="bloom luminance threshold (relative to white)")
    ap.add_argument("--no-bloom", action="store_true", help="skip bloom entirely")
    ap.add_argument("--tonemap", choices=["aces", "reinhard"], default="aces",
                    help="tone curve: aces (filmic) or reinhard (gentler toe, more saturated)")
    args = ap.parse_args()

    img = read_image(args.input)
    strength = 0.0 if args.no_bloom else args.bloom
    out = grade(img, args.exposure, strength, args.threshold, args.white, args.tonemap)
    Image.fromarray(out).save(args.output)
    print(f"wrote {args.output} ({out.shape[1]}x{out.shape[0]})")


if __name__ == "__main__":
    main()

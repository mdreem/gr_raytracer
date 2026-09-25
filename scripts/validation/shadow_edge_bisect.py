#!/usr/bin/env python3
"""Shadow-edge check against Bardeen's analytic Kerr critical curve.

Usage:
    shadow_edge_bisect.py BIN R0 A EPS GEOMETRY HI_RAD

    BIN       path to the release binary (target/release/gr_raytracer)
    R0        camera radius on the -x axis (Cartesian), in units of r_s
    A         spin a in units of r_s (a/M = 2a/r_s); ignored for Schwarzschild
    EPS       integrator --epsilon
    GEOMETRY  KerrBL | Kerr | Schwarzschild
    HI_RAD    upper bisection bracket on the local angle (rad); must escape
              (0.012 is enough at R0 = 1000, 0.6 at R0 = 18)

For each side (prograde, retrograde) equatorial rays are fired with
`render-ray-at`, the capture/escape boundary is bisected on the local angle
psi, psi is converted exactly to the impact parameter b = L/E for the tetrad
that `render-ray-at` uses (ZAMO for KerrBL, free-fall for Schwarzschild), and
b is compared with xi(r_ph) at the equatorial photon orbits (James et al. 2015,
eqs. A.5-A.6 with q = 0). This validates the geodesic integrator, the tetrad
and the conserved-quantity seeding end to end.

For the Kerr (Kerr-Schild) backend the `render-ray-at` tetrad is the Eulerian
observer, whose velocity relative to the ZAMO is ~ r_s/R0, so the conversion
is only exact in the large-R0 limit; use R0 >= 1000 for it.

See docs/physics-review-2026-09.md, section 3.1, for reference results.
"""
import csv, math, os, subprocess, sys, tempfile

BIN, R0, A, EPS, GEOM, HI = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), sys.argv[5], float(sys.argv[6])
RS = 1.0
M = 0.5


def scene():
    return f"""celestial_temperature = 0.0
objects = []
[celestial_texture.Checker]
beaming_exponent = 0.0
width = 10.0
height = 10.0
color1 = [255, 255, 255]
color2 = [0, 0, 0]
[geometry_type.{GEOM}]
radius = {RS}
{"" if GEOM == "Schwarzschild" else f"a = {A}"}
horizon_epsilon = 1e-4
"""


def run_ray(cfg, direction, out):
    cmd = [BIN, "--config-file", cfg, "--max-steps", "400000", "--max-radius", "3000",
           "--epsilon", str(EPS), "--step-size", "0.001",
           "render-ray-at", f"--position=-{R0},0,0",
           "--direction=" + ",".join(f"{d:.17g}" for d in direction), "--filename", out]
    subprocess.run(cmd, check=True, capture_output=True)
    rows = list(csv.DictReader(open(out)))
    last = rows[-1]
    return math.sqrt(float(last["x"]) ** 2 + float(last["y"]) ** 2 + float(last["z"]) ** 2)


def zamo_b(psi):
    if GEOM == "Schwarzschild":
        # free-fall tetrad (from rest at infinity): p = T - cos(psi) R + sin(psi) Phi
        return R0 * math.sin(psi) / (1.0 + math.cos(psi) * math.sqrt(RS / R0))
    r, th, a = R0, math.pi / 2, A
    sig = r * r + a * a * math.cos(th) ** 2
    sin2 = math.sin(th) ** 2
    g_tt = -(1 - RS * r / sig)
    g_tp = -a * RS * r * sin2 / sig
    g_pp = (r * r + a * a + a * a * RS * r * sin2 / sig) * sin2
    omega = -g_tp / g_pp
    alpha = math.sqrt(-(g_tt + omega * g_tp))
    n_phi = math.sin(psi)
    return n_phi * math.sqrt(g_pp) / (alpha + omega * n_phi * math.sqrt(g_pp))


def analytic_edges(a):
    def xi(r):
        delta = r * r - 2 * M * r + a * a
        return (M * (r * r - a * a) - r * delta) / (a * (r - M))
    if abs(a) < 1e-12:
        b = 3 * math.sqrt(3) * M
        return b, -b
    r_pro = 2 * M * (1 + math.cos(2.0 / 3.0 * math.acos(-a / M)))
    r_ret = 2 * M * (1 + math.cos(2.0 / 3.0 * math.acos(a / M)))
    return xi(r_pro), xi(r_ret)


def direction_for(psi):
    """`render-ray-at` direction for local angle psi (positive = prograde)."""
    if GEOM == "Schwarzschild":
        # Cartesian direction; the CLI maps it onto the free-fall tetrad.
        return (-math.cos(psi), -math.sin(psi), 0.0)
    if GEOM == "Kerr":
        # Kerr-Schild Cartesian: camera at (-R0,0,0), inward = +x, phi_hat = -y.
        return (math.cos(psi), -math.sin(psi), 0.0)
    # KerrBL tetrad axes are (phi_hat, theta_hat, r_hat).
    return (math.sin(psi), 0.0, -math.cos(psi))


def bisect(cfg, sign, iters=16):
    lo, hi = 0.0, HI
    tmp = tempfile.NamedTemporaryFile(suffix=".csv", delete=False).name
    r_hi = run_ray(cfg, direction_for(sign * hi), tmp)
    assert r_hi > 2000.0, f"bracket hi={hi} does not escape (r_last={r_hi})"
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        psi = sign * mid
        r_last = run_ray(cfg, direction_for(psi), tmp)
        if r_last < 5.0:
            lo = mid
        elif r_last > 2000.0:
            hi = mid
        else:
            raise RuntimeError(f"ray neither captured nor escaped: r_last={r_last}")
    os.unlink(tmp)
    return sign * 0.5 * (lo + hi)


b_pro, b_ret = analytic_edges(A)
cfg = tempfile.NamedTemporaryFile(suffix=".toml", delete=False, mode="w")
cfg.write(scene())
cfg.close()
for sign, ref, label in [(+1, b_pro, "prograde"), (-1, b_ret, "retrograde")]:
    b = zamo_b(bisect(cfg.name, sign))
    print(f"{os.path.basename(BIN):<10} R0={R0:<6} a={A:<6} {GEOM:<7} eps={EPS:<6} {label:<10} "
          f"b_num={b:+.5f} b_analytic={ref:+.5f} rel.err={(b - ref) / ref:+.2e}")
    sys.stdout.flush()
os.unlink(cfg.name)

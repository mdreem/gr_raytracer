# Example Render Commands

Example commands that point the camera at the black hole from the same physical position
(`x=-10, y=0, z=-0.5`) with 1000×1000 resolution. The camera angle parameters differ per
geometry because each uses a different local tetrad convention.

## Camera orientation notes

| Geometry | Coordinate system | Tetrad.z (forward) | Angles to face BH |
|----------|-------------------|--------------------|-------------------|
| Schwarzschild | Spherical | Radial outward (+r) | `theta=-π` flips to inward |
| Kerr | Cartesian (KS) | Cartesian z-axis | `theta≈π/2, psi≈-π/2` rotates to +x |
| KerrBL | Boyer-Lindquist | Radial outward (+r) | `theta=-π` flips to inward |

## Schwarzschild

```bash
cargo run --release -- \
  --width=1000 --height=1000 \
  --max-steps=1000000 \
  --camera-position=-10,0,-0.5 \
  --theta=-3.14159 --psi=0 --phi=0 \
  --config-file scene-definitions/schwarzschild.toml \
  render --filename render-schwarzschild.png
```

## Kerr (Kerr-Schild Cartesian coordinates)

```bash
cargo run --release -- \
  --width=1000 --height=1000 \
  --max-steps=1000000 \
  --camera-position=-10,0,-0.5 \
  --theta=1.52 --psi=-1.57 --phi=0 \
  --config-file scene-definitions/kerr.toml \
  render --filename render-kerr.png
```

## KerrBL (Boyer-Lindquist coordinates)

```bash
cargo run --release -- \
  --width=1000 --height=1000 \
  --max-steps=1000000 \
  --camera-position=-10,0,-0.5 \
  --theta=-3.14159 --psi=0 --phi=0 \
  --config-file scene-definitions/kerr-bl.toml \
  render --filename render-kerr-bl.png
```

## Notes

- `--max-steps=1000000` is used here for safety, but the default (15000) works fine for all
  three geometries at this camera position. Increase it if rays near the horizon terminate
  early with a `MaxSteps` stop reason.
- KerrBL is significantly faster than Kerr (~8-9× at 500×500) because it uses separated
  equations of motion (Carter constant) instead of a generic ODE integrator.
- The spin parameter `a=0.499` (nearly maximal) is set in the TOML files.
- Rendering is parallelised with Rayon. To control the number of threads, set
  `RAYON_NUM_THREADS` before running:
  ```bash
  RAYON_NUM_THREADS=4 cargo run --release -- ...
  ```

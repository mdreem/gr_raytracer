# Analysis of `README.md` vs Codebase (Correctness-Checked)

This document validates the claims in the current `analysis.md` against the implementation and lists what is still
missing from `README.md`.

## Feature Snapshot

- Multi-geometry raytracing: Euclidean, EuclideanSpherical, Schwarzschild, Kerr.
- Adaptive geodesic integration with RKF45 and error-controlled step sizing.
- Relativistic effects in shading: gravitational/Doppler redshift and beaming.
- Physically motivated emission: black-body spectrum integration in CIE XYZ.
- Rich scene primitives: Sphere, Disc, and Perlin-noise-based VolumetricDisc.
- Flexible materials/textures: Bitmap, Checker, and BlackBody mappers.
- High-fidelity outputs: standard image export plus HDR (`.hdr`) support.
- Debug/inspection tooling: per-pixel ray export and ray-at-position export (geometry-limited).
- Config-driven scenes via TOML with pluggable geometry, textures, and object stacks.

## Corrections to Previous Analysis

1. Integrator claim needs correction:

- Previous text said **RK4**.
- Code uses an **adaptive Runge-Kutta-Fehlberg RKF45 (4/5)** integrator with dynamic step size (
  `src/rendering/runge_kutta.rs`).

2. `render-ray-at` support is partial:

- Previous text implied it as a general CLI feature.
- It works for **Schwarzschild** and **Kerr**.
- For **Euclidean** and **EuclideanSpherical**, it currently panics with “not supported yet” (`src/main.rs`).

3. “Background stars” wording is too specific:

- Implementation is a **generic celestial sphere texture map** (`celestial_texture`), not a dedicated star
  catalog/engine.

4. Example scene files are not all valid as committed:

- `scene-definitions/euclidean.toml` and `scene-definitions/euclidean-spherical.toml` define `objects.Sphere` without
  required `temperature`.
- Running with `scene-definitions/euclidean.toml` currently fails TOML parsing with `missing field temperature`.

## Confirmed Implemented Features (and still under-documented)

1. Geometry support:

- `Euclidean`, `EuclideanSpherical`, `Schwarzschild`, `Kerr` are all first-class geometry types in config (
  `src/configuration.rs`).

2. Objects:

- `Sphere`, `Disc`, and `VolumetricDisc` are implemented and configurable (`src/configuration.rs`, `src/cli/shared.rs`).

3. Textures:

- `Bitmap`, `Checker`, and `BlackBody` texture modes are supported.
- Bitmap texture sampling uses **bilinear interpolation** (`src/rendering/texture.rs`).

4. Physics/appearance pipeline:

- Redshift computation (`src/rendering/redshift.rs`).
- Beaming factor application (`src/rendering/color.rs`, `src/rendering/texture.rs`).
- Black-body spectral integration -> CIE XYZ (`src/rendering/black_body_radiation.rs`).
- sRGB conversion for non-HDR output (`src/rendering/color.rs`).

5. Temperature models:

- Constant temperature model.
- Kerr/Page-Thorne-style radial temperature model with ISCO handling (`src/rendering/temperature.rs`).

6. Output modes:

- Standard images (PNG etc.) and **HDR** (`.hdr`) output are supported (`src/rendering/raytracer.rs`).

7. CLI rendering capabilities:

- `render` (full frame or sub-window via `--from-row/--from-col/--to-row/--to-col`).
- `render-ray` (pixel-index ray trajectory to CSV).
- `render-ray-at` (custom position/direction trajectory to CSV, geometry-limited as noted above).

## Missing from `README.md` (recommended additions)

1. Explicit CLI reference section:

- Document all subcommands and key flags from `--help`:
    - global: `--width`, `--height`, `--step-size`, `--max-steps`, `--max-radius`, `--epsilon`, `--camera-position`,
      `--phi`, `--theta`, `--psi`, `--config-file`
    - `render` windowing options
    - `render-ray` / `render-ray-at` CSV outputs

2. Clarify geometry limitation for `render-ray-at`:

- Mention it currently supports Kerr/Schwarzschild only.

3. Scene TOML schema docs:

- Explain `geometry_type` variants.
- Explain `objects` variants and important fields.
- Explain texture variants (`Bitmap`, `Checker`, `BlackBody`).
- Include `color_normalization` modes (`NoNormalization`, `Chromaticity`, `EqualLuminance`).

4. Volumetric disc configuration docs:

- Document parameters like `thickness`, `density_multiplier`, `absorption`, `scattering`, `num_octaves`, `noise_scale`,
  `noise_offset`, `brightness_reference_temperature`, etc.

5. Integrator description:

- Replace RK4-style wording with adaptive RKF45 and mention tolerance control via `--epsilon`.

6. Coordinate/frame orientation docs:

- Describe camera rotation angles `phi/theta/psi`, and that camera rays are generated using local tetrads and Lorentz
  transformation.

7. Output semantics:

- Non-HDR path applies XYZ normalization (config-dependent), then converts to sRGB bytes.
- HDR path currently writes raw XYZ components into RGB float channels (not XYZ->linear-sRGB transformed).

8. Document scene-definition validity expectations:

- Mention that all object variants have required fields (e.g. `Sphere.temperature` is required).
- Add a note that sample scene files should be kept parse-valid and tested.

9. Document runtime validation for volumetric discs:

- Current code rejects invalid values at runtime (`outer_radius > inner_radius`, `thickness > 0`, `max_steps > 0`,
  `step_size > 0`, `brightness_reference_temperature > 0`, `absorption >= 0`, `scattering >= 0`).
- These constraints are currently undocumented in README.

## References Missing from README

Useful to add in the README “Sources” section (already reflected in code/comments):

1. [Seeing Relativity (Schwarzschild ray tracing)](https://arxiv.org/abs/1511.06025)
2. [Johannsen 2011 (angular velocity used in Kerr-related formulas)](https://arxiv.org/abs/1104.5499)
3. Page & Thorne (1974), thin-disk accretion model background used for temperature profile
4. CIE 1931 color matching functions and sRGB transfer/conversion references

## To Implement

1. HDR output color-space fix:

- In the `.hdr` render path, convert XYZ tristimulus values to **linear sRGB** before writing.
- Use existing `xyz_to_linear_srgb` from `src/rendering/color.rs`.
- Current behavior writes XYZ directly into RGB channels, which is not a proper linear RGB HDR encoding.
- Decide negative-channel handling for Radiance HDR (`.hdr`): clamp to `0.0` or move to EXR for signed linear output.

2. TOML schema documentation generation:

- Add code-level schema derives (e.g. `schemars::JsonSchema`) to config types (`RenderConfig`, `GeometryType`,
  `ObjectsConfig`, `TextureConfig`, etc.).
- Add a schema export command/binary to generate JSON Schema from Rust types.
- Add an optional docs step/script to generate Markdown from schema (e.g. `docs/config-schema.md`).
- Keep examples in README/manual docs for TOML table style where generated schema is not sufficiently TOML-specific.

3. Frame-level postprocessor pipeline (without deep renderer changes):

- Add a single postprocessing stage that consumes final-frame XYZ values and produces output-ready pixels based on
  CLI/config parameters.
- Keep physics/integration unchanged; do not push postprocessing logic into geometry/object/intersection code.
- Suggested first operators: exposure control, tone mapping (Reinhard or ACES fit), optional bloom/glare, and output
  transform.
- Output transform should branch by format:
    - non-HDR: XYZ -> linear RGB -> tone map -> sRGB compand/encode
    - HDR: XYZ -> linear RGB (no compand), with explicit negative-channel policy for Radiance `.hdr`.

4. Keep sample TOML files schema-valid in CI:

- Add a lightweight test/command that deserializes every file under `scene-definitions/*.toml`.
- This would catch drift like missing required fields (e.g. `Sphere.temperature`) before docs/examples go stale.

5. `manim` animation setup documentation:
The Manim animation command is presented without any guidance on installing the required Python dependencies. Since
scripts/animate-rays/main.py imports pandas, numpy, and scipy (in addition to manim), consider documenting how to set up
a venv and install the script requirements (and ensure the listed requirements include scipy) so the command works when
copy/pasted.


# Scripts

This directory contains helper scripts for generating trajectories, rendering batches, creating animations, and plotting
ray data.

Run commands from the repository root.

## Layout

```text
scripts/
  animate-rays/           # Manim scene + config for ray animations
  plotting/               # Trajectory and direction plotting tools
  rays/                   # Ray export scripts
  rendering/              # Camera paths, frame rendering, movie/Kerr batch generation
    resources/            # Rendering script templates/assets
  textures/               # UV/checker texture generation
```

## Requirements

- `cargo` (Rust toolchain)
- `uv` (for Python scripts and Manim)
- `ffmpeg` (for `rendering/create_movie_from_images.sh`)
- `magick` / ImageMagick (for `rendering/create_kerr_images.sh`)

Install Python dependencies:

```sh
uv sync
```

## Common workflows

### Camera fly-through movie

1. Generate camera positions:

```sh
uv run scripts/rendering/create_camera_trajectory.py
```

2. Render one frame per camera position:

```sh
uv run scripts/rendering/create_images_from_camera_positions.py
```

3. Build an MP4 from the rendered frames:

```sh
bash scripts/rendering/create_movie_from_images.sh
```

### Ray export + Manim animation

Generate rays:

```sh
bash scripts/rays/create_rays_positions.sh schwarzschild
# or
bash scripts/rays/create_rays_from_camera.sh
```

Then render the animation:

```sh
uv run manim scripts/animate-rays/main.py AnimateRays
```

## Script reference

### `rendering/create_camera_trajectory.py`

- Creates `camera_positions.csv` with interpolated camera positions from hard-coded control points.
- Edit `points` and `N` in the script to change path shape and number of frames.

```sh
uv run scripts/rendering/create_camera_trajectory.py
```

### `rendering/create_images_from_camera_positions.py`

- Reads `camera_positions.csv`.
- Runs the Rust renderer for each row and writes frames to `movie/render-XXX.png`.
- Adjust `WIDTH`, `HEIGHT`, `CONFIG_FILE`, and `POSITIONS_FILE` in the script as needed.

```sh
uv run scripts/rendering/create_images_from_camera_positions.py
```

### `rendering/create_movie_from_images.sh`

- Uses `ffmpeg` to combine `movie/render-%03d.png` into `movie/output.mp4` (30 FPS by default).

```sh
bash scripts/rendering/create_movie_from_images.sh
```

### `rendering/create_kerr_images.sh`

- Sweeps Kerr spin parameter `a` from `0.000` to `0.499`.
- Creates per-frame scene files in `kerr_images/`, renders PNGs, annotates each frame, then builds
  `kerr_animation.gif`.
- Optional first argument: custom render command (default: `cargo run --release --`).
- Optional grid overlay flags:
  - `--grid` / `--no-grid`
  - `--grid-size N`
  - `--grid-width N --grid-height N`

```sh
bash scripts/rendering/create_kerr_images.sh
# Example with binary and grid overlay:
bash scripts/rendering/create_kerr_images.sh "target/release/gr_raytracer" --grid --grid-size 40
```

### `rays/create_rays_positions.sh`

- Exports ray trajectories to `rays/` by sweeping the initial ray position.
- Usage: `create_rays_positions.sh <geometry> [shift] [command]`
- `<geometry>`: `schwarzschild` or `kerr`
- `shift` defaults to `0`
- `command` defaults to `cargo run --release --`

```sh
bash scripts/rays/create_rays_positions.sh schwarzschild
bash scripts/rays/create_rays_positions.sh kerr 2.0
```

### `rays/create_rays_from_camera.sh`

- Exports rays from a fixed camera position while varying horizontal pixel coordinate.
- Writes CSV files to `rays/`.
- Optional argument: custom render command (default: `cargo run --release --`).

```sh
bash scripts/rays/create_rays_from_camera.sh
```

### `animate-rays/main.py`

- Manim scene that reads `rays/*.csv` and animates trajectories (`AnimateRays` scene).
- Config defaults are in `scripts/animate-rays/manim.cfg`.

```sh
uv run manim scripts/animate-rays/main.py AnimateRays
```

### `plotting/plot_trajectories_3d.py`

- Plots one or more trajectory CSVs with `x,y,z` columns.
- Draws a central sphere (`r=1`) and supports clipping by `--max-radius`.

```sh
uv run scripts/plotting/plot_trajectories_3d.py rays/ray-0000.00.csv
uv run scripts/plotting/plot_trajectories_3d.py rays/ray-0000.00.csv --save trajectory.png --max-radius 15
```

### `plotting/plot_planar_trajectories.py`

- 2D plot for multiple trajectories.
- Accepts CSVs with either:
  - `r,phi` columns (polar data converted to Cartesian), or
  - `x,y,z` Cartesian data (plots `x` against `z`).

```sh
uv run scripts/plotting/plot_planar_trajectories.py rays/ray-0000.00.csv rays/ray-0005.00.csv
```

### `plotting/direction_plotter.py`

- 3D scatter plot of direction (`d_x1,d_x2,d_x3`) and momentum (`m_x1,m_x2,m_x3`) from `rays.csv`.

```sh
uv run scripts/plotting/direction_plotter.py
```

### `textures/create_uv_map.py`

- Generates checkerboard texture atlases in `resources/`:
  - `celestial.png`
  - `sphere.png`
  - `disk.png`

```sh
uv run scripts/textures/create_uv_map.py
```

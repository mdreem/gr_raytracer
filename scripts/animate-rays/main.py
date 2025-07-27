from manim import *
import pandas as pd
from scipy.interpolate import interp1d
import glob
import dataclasses


@dataclasses.dataclass
class Interpolators:
    t_points: np.ndarray
    x_interpolator: interp1d
    y_interpolator: interp1d


class AnimateRays(Scene):
    """Animate rays from a central point"""

    def get_color_linear(self, index, total):
        ratio = index / (total - 1) if total > 1 else 0
        return interpolate_color(BLUE, RED, ratio)

    def read_trajectories(self, file_paths):
        trajectories = []

        csv_files = [f for f in file_paths if f.endswith(".csv")]
        csv_files.sort()

        for i, file_path in enumerate(csv_files):
            print(f"Reading {file_path}...")
            df = pd.read_csv(file_path)

            scale = 1.0/2.0
            df["x"] *= scale
            df["y"] *= scale
            df["z"] *= scale

            trajectories.append((df, self.get_color_linear(i, len(csv_files))))
        return trajectories

    def construct(self):
        self.camera.background_color = DARK_GRAY

        polarplane_pi = PolarPlane(
            azimuth_units="PI radians",
            radius_step=2.0,
            radius_max=10.0,
            size=10,
            background_line_style={"stroke_opacity": 0.4},
        ).add_coordinates()
        self.add(polarplane_pi)

        file_paths = glob.glob("rays/*.csv")
        trajectories = self.read_trajectories(file_paths)

        # Create dots and paths for each trajectory
        dots = []
        paths = []

        for i, (df, color) in enumerate(trajectories):
            dot = Dot(color=color, radius=0.02)
            path = TracedPath(dot.get_center, stroke_color=color, stroke_width=2)

            # Set initial position
            dot.move_to([df.iloc[0]["y"], df.iloc[0]["z"], 0])

            dots.append(dot)
            paths.append(path)
            self.add(path, dot)

        center_disc = Circle(
            radius=1.0/2.0,           # Radius = 1 unit
            color=RED,          # Border color
            fill_color=RED_E, # Fill color
            fill_opacity=0.8,   # Transparency
            stroke_width=2      # Border thickness
        )
        self.add(center_disc)

        interpolators = self.create_interpolation_functions(trajectories)

        # Animate all trajectories
        def update_all_dots(alpha):
            updates = []
            for i, (dot, interpolator) in enumerate(zip(dots, interpolators)):
                sim_time = interpolator.t_points[0] + alpha * (
                    interpolator.t_points[-1] - interpolator.t_points[0]
                )
                x_pos = interpolator.x_interpolator(sim_time)
                y_pos = interpolator.y_interpolator(sim_time)
                updates.append(dot.animate.move_to([x_pos, y_pos, 0]))
            return updates

        # Run animation in steps
        steps = 100
        for i in range(steps):
            alpha = i / (steps - 1)
            animations = update_all_dots(alpha)

            self.play(*animations, run_time=0.08, rate_func=linear)

        self.wait(2)

    def create_interpolation_functions(self, trajectories):
        interpolators = []
        for df, _ in trajectories:
            t_points = df["t"].values
            x_interp = interp1d(t_points, df["y"].values, kind="linear")
            y_interp = interp1d(t_points, df["z"].values, kind="linear")
            interpolators.append(
                Interpolators(
                    t_points=t_points, x_interpolator=x_interp, y_interpolator=y_interp
                )
            )
        return interpolators

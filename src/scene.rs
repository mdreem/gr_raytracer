use crate::camera::Ray;
use crate::four_vector::FourVector;
use crate::runge_kutta::rk4;
use nalgebra::{Const, OVector};

pub struct Scene {
    max_steps: usize,
    step_size: f64,
    center_sphere_radius: f64,
}

#[derive(Debug, PartialEq)]
pub struct Color {
    r: u8,
    g: u8,
    b: u8,
}

impl Color {
    fn new(r: u8, g: u8, b: u8) -> Color {
        Color { r, g, b }
    }

    pub fn get_as_array(&self) -> [u8; 3] {
        [self.r, self.g, self.b]
    }
}

type EquationOfMotionState = OVector<f64, Const<8>>;

fn geodesic(_: f64, y: &EquationOfMotionState) -> EquationOfMotionState {
    let y_new =
        EquationOfMotionState::from_column_slice(&[y[4], y[5], y[6], y[7], 0.0, 0.0, 0.0, 0.0]);
    y_new
}

pub fn get_position(y: &EquationOfMotionState) -> FourVector {
    FourVector::new(y[0], y[1], y[2], y[3])
}

impl Scene {
    pub fn new() -> Self {
        Scene {
            max_steps: 1000,
            step_size: 0.01,
            center_sphere_radius: 2.0,
        }
    }

    fn compute_sphere_intersection(&self, y_start: &FourVector, y_end: &FourVector) -> bool {
        let y_start_position = y_start.get_as_vector();
        let y_end_position = y_end.get_as_vector();
        let r_start = y_start_position[1] * y_start_position[1]
            + y_start_position[2] * y_start_position[2]
            + y_start_position[3] * y_start_position[3];
        let r_end = y_end_position[1] * y_end_position[1]
            + y_end_position[2] * y_end_position[2]
            + y_end_position[3] * y_end_position[3];

        if (r_start >= self.center_sphere_radius.powi(2)
            && r_end <= self.center_sphere_radius.powi(2))
            || (r_start <= self.center_sphere_radius.powi(2)
                && r_end >= self.center_sphere_radius.powi(2))
        {
            return true;
        }

        false
    }

    pub fn color_of_ray(&self, ray: &Ray) -> Color {
        let mut t = 0.0;

        let direction = ray.direction.get_as_vector();
        let mut y = EquationOfMotionState::from_column_slice(&[
            ray.position[0],
            ray.position[1],
            ray.position[2],
            ray.position[3],
            direction[0],
            direction[1],
            direction[2],
            direction[3],
        ]);

        for _i in 0..self.max_steps {
            let last_y = y;
            y = rk4(&y, t, self.step_size, geodesic);

            if self.compute_sphere_intersection(&get_position(&last_y), &get_position(&y)) {
                return Color::new(255, 0, 0);
            }

            t += self.step_size;
        }
        Color::new(0, 0, 0)
    }
}

#[cfg(test)]
mod tests {
    use crate::camera::Camera;
    use crate::scene::{Color, Scene};
    use nalgebra::Vector4;

    #[test]
    fn test_color_of_ray_hits_sphere() {
        let camera = Camera::new(
            Vector4::new(0.0, -10.0, 0.0, 0.0),
            std::f64::consts::PI / 2.0,
            11,
            11,
        );
        let ray = camera.get_ray_for(6, 6);
        let scene = Scene::new();

        let color = scene.color_of_ray(&ray);

        assert_eq!(color, Color::new(255, 0, 0));
    }

    #[test]
    fn test_color_of_ray_misses_sphere() {
        let camera = Camera::new(
            Vector4::new(0.0, -10.0, 0.0, 0.0),
            std::f64::consts::PI / 2.0,
            11,
            11,
        );
        let ray = camera.get_ray_for(0, 0);
        let scene = Scene::new();

        let color = scene.color_of_ray(&ray);

        assert_eq!(color, Color::new(0, 0, 0));
    }
}

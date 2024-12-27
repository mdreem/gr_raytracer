use crate::camera::Ray;
use crate::four_vector::FourVector;
use crate::runge_kutta::rk4;
use nalgebra::{Const, OVector, Vector3};

pub struct Scene {
    max_steps: usize,
    step_size: f64,
    center_sphere_radius: f64,
    center_disk_outer_radius: f64,
    center_disk_inner_radius: f64,
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

fn radial_distance_spacial_part(pos: &FourVector) -> f64 {
    let v = pos.get_as_vector();
    v[1] * v[1] + v[2] * v[2] + v[3] * v[3]
}

impl Scene {
    pub fn new(
        max_steps: usize,
        step_size: f64,
        center_sphere_radius: f64,
        center_disk_inner_radius: f64,
        center_disk_outer_radius: f64,
    ) -> Scene {
        Scene {
            max_steps,
            step_size,
            center_sphere_radius,
            center_disk_outer_radius,
            center_disk_inner_radius,
        }
    }

    // TODO: explicitly construct the ray. Follow the integration. Some intervals seem to be skipped
    // here. See with current test setup. Intersection should be at t=7.63. With z=-2.442748091.
    // The intersection should be with an interval crossing y=0. But it seems to happen near 0 with
    // both coordinates.
    fn intersects_with_disk(&self, y_start: &FourVector, y_end: &FourVector) -> bool {
        // z x y
        let normal = Vector3::new(0.0, 0.0, 1.0);
        let center = Vector3::new(0.0, 0.0, 0.0);
        let y_start_spacial = y_start.get_spacial_vector();
        let y_end_spacial = y_end.get_spacial_vector();
        let direction = y_end_spacial - y_start_spacial;

        let p1 = (y_start_spacial - center).transpose() * normal;
        let p2 = direction.transpose() * normal;

        // TODO: p2 can be 0 if parallel -> handle
        let t = p1[0] / p2[0]; // plane intersection parameter.

        if !(t >= 0.0 && t <= 1.0) {
            return false;
        }

        let intersection_point = y_start_spacial + t * direction;
        let rr = intersection_point[0] * intersection_point[0]
            + intersection_point[1] * intersection_point[1]
            + intersection_point[2] * intersection_point[2];

        rr >= self.center_disk_inner_radius * self.center_disk_inner_radius
            && rr <= self.center_disk_outer_radius * self.center_disk_outer_radius
    }

    fn intersects_with_sphere(&self, y_start: &FourVector, y_end: &FourVector) -> bool {
        let r_start = radial_distance_spacial_part(&y_start);
        let r_end = radial_distance_spacial_part(&y_end);

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

            if self.intersects_with_disk(&get_position(&last_y), &get_position(&y)) {
                // println!("intersection: {:?}", last_y);
                return Color::new(0, 0, 255);
            }

            if self.intersects_with_sphere(&get_position(&last_y), &get_position(&y)) {
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
        let scene = Scene::new(1000, 0.01, 2.0, 0.2, 0.3);

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
        let scene = Scene::new(1000, 0.01, 2.0, 0.2, 0.3);

        let color = scene.color_of_ray(&ray);

        assert_eq!(color, Color::new(0, 0, 0));
    }

    #[test]
    fn test_intersects_with_disk() {
        let camera = Camera::new(
            Vector4::new(0.0, -10.0, 0.0, 1.0),
            std::f64::consts::PI / 4.0,
            101,
            101,
        );
        let ray = camera.get_ray_for(43, 51);
        let scene = Scene::new(1000, 0.01, 1.0, 2.0, 7.0);

        let color = scene.color_of_ray(&ray);
        assert_eq!(color, Color::new(0, 0, 255));
    }
}

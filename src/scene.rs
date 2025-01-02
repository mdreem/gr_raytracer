use crate::camera::Ray;
use crate::four_vector::{CoordinateSystem, FourVector};
use crate::geometry::Geometry;
use crate::runge_kutta::rk4;
use crate::spherical_coordinates_helper::spherical_to_cartesian;
use image::{DynamicImage, GenericImageView, ImageReader};
use nalgebra::{Const, OVector, Vector3, Vector4};
use std::f64::consts::PI;

pub struct Scene<T: TextureMap, G: Geometry> {
    max_steps: usize,
    max_radius_sq: f64,
    step_size: f64,
    center_sphere_radius: f64,
    center_disk_outer_radius: f64,
    center_disk_inner_radius: f64,
    celestial_map: T,
    center_disk_map: T,
    center_sphere_map: T,
    pub geometry: G,
}

#[derive(Debug, PartialEq, Copy, Clone)]
pub struct Color {
    r: u8,
    g: u8,
    b: u8,
}

pub trait TextureMap: Sync {
    fn color_at_uv(&self, u: f64, v: f64) -> Color;
}

pub struct TextureMapper {
    image: DynamicImage,
}

pub struct CheckerMapper {
    width: f64,
    height: f64,
    c1: Color,
    c2: Color,
}

impl Color {
    pub fn new(r: u8, g: u8, b: u8) -> Color {
        Color { r, g, b }
    }

    pub fn get_as_array(&self) -> [u8; 3] {
        [self.r, self.g, self.b]
    }
}

pub type EquationOfMotionState = OVector<f64, Const<8>>;

// TODO: replace this with a more natural function, that does not need to use cartesian coordinates.
pub fn get_position(y: &EquationOfMotionState, coordinate_system: CoordinateSystem) -> FourVector {
    match coordinate_system {
        CoordinateSystem::Cartesian => FourVector::new_cartesian(y[0], y[1], y[2], y[3]),
        CoordinateSystem::Spherical => {
            let t_vec = spherical_to_cartesian(&Vector4::new(y[0], y[1], y[2], y[3]));
            FourVector::new_cartesian(0.0, t_vec[0], t_vec[1], t_vec[2]) // TODO: try to do all of this without this conversion.
        }
    }
}

fn radial_distance_spatial_part_squared(pos: &FourVector) -> f64 {
    let v = pos.get_as_vector();
    v[1] * v[1] + v[2] * v[2] + v[3] * v[3]
}

impl TextureMapper {
    pub fn new(filename: String) -> TextureMapper {
        let image = ImageReader::open(filename)
            .expect("Failed to open image file")
            .decode()
            .expect("Failed to decode image");

        TextureMapper { image }
    }
}

impl TextureMap for TextureMapper {
    fn color_at_uv(&self, u: f64, v: f64) -> Color {
        let (width, height) = self.image.dimensions();
        let pixel = self.image.get_pixel(
            (((width - 1) as f64) * u) as u32,
            (((height - 1) as f64) * v) as u32,
        );
        Color::new(pixel[0], pixel[1], pixel[2])
    }
}

impl CheckerMapper {
    pub fn new(width: f64, height: f64, c1: Color, c2: Color) -> CheckerMapper {
        CheckerMapper {
            width,
            height,
            c1,
            c2,
        }
    }
}

impl TextureMap for CheckerMapper {
    fn color_at_uv(&self, u: f64, v: f64) -> Color {
        let ut = (u * self.width).floor() as usize;
        let vt = (v * self.height).floor() as usize;

        if (ut + vt) % 2 == 0 {
            self.c1
        } else {
            self.c2
        }
    }
}

impl<T: TextureMap, G: Geometry> Scene<T, G> {
    pub fn new(
        max_steps: usize,
        max_radius: f64,
        step_size: f64,
        center_sphere_radius: f64,
        center_disk_inner_radius: f64,
        center_disk_outer_radius: f64,
        celestial_map: T,
        center_disk_map: T,
        center_sphere_map: T,
        geometry: G,
    ) -> Scene<T, G> {
        Scene {
            max_steps,
            max_radius_sq: max_radius * max_radius,
            step_size,
            center_sphere_radius,
            center_disk_outer_radius,
            center_disk_inner_radius,
            celestial_map,
            center_disk_map,
            center_sphere_map,
            geometry,
        }
    }

    // TODO: explicitly construct the ray. Follow the integration. Some intervals seem to be skipped
    // here. See with current test setup. Intersection should be at t=7.63. With z=-2.442748091.
    // The intersection should be with an interval crossing y=0. But it seems to happen near 0 with
    // both coordinates.
    fn intersects_with_disk(&self, y_start: &FourVector, y_end: &FourVector) -> Option<Color> {
        // z x y
        let normal = Vector3::new(0.0, 1.0, 0.0);
        let center = Vector3::new(0.0, 0.0, 0.0);
        let y_start_spatial = y_start.get_spatial_vector();
        let y_end_spatial = y_end.get_spatial_vector();
        let direction = y_end_spatial - y_start_spatial;

        let p1 = (y_start_spatial - center).transpose() * normal;
        let p2 = direction.transpose() * normal;

        // TODO: p2 can be 0 if parallel -> handle
        let t = p1[0] / p2[0]; // plane intersection parameter.

        if !(t >= 0.0 && t <= 1.0) {
            return None;
        }

        let intersection_point = y_start_spatial + t * direction;
        let rr = intersection_point[0] * intersection_point[0]
            + intersection_point[1] * intersection_point[1]
            + intersection_point[2] * intersection_point[2];

        if rr >= self.center_disk_inner_radius * self.center_disk_inner_radius
            && rr <= self.center_disk_outer_radius * self.center_disk_outer_radius
        {
            let vector_in_plane = intersection_point - center;
            // TODO: properly implement finding intersection and compute the values accordingly.

            let phi = vector_in_plane[2].atan2(vector_in_plane[0]); // phi in x-z plane.
            let u = (PI + phi) / (2.0 * PI);
            let v = (rr.sqrt() - self.center_disk_inner_radius)
                / (self.center_disk_outer_radius - self.center_disk_inner_radius);

            let color = self.center_disk_map.color_at_uv(u, v);
            Some(color)
        } else {
            None
        }
    }

    // TODO: handle spherical coordinates here!
    fn intersects_with_sphere(&self, y_start: &FourVector, y_end: &FourVector) -> Option<Color> {
        let r_start = radial_distance_spatial_part_squared(&y_start);
        let r_end = radial_distance_spatial_part_squared(&y_end);

        if (r_start >= self.center_sphere_radius.powi(2)
            && r_end <= self.center_sphere_radius.powi(2))
            || (r_start <= self.center_sphere_radius.powi(2)
                && r_end >= self.center_sphere_radius.powi(2))
        {
            let point_on_sphere = y_start.get_as_spherical(); // approximate y_start als intersection point.

            let u = (PI + point_on_sphere[2]) / (2.0 * PI);
            let v = point_on_sphere[1] / PI;

            let color = self.center_sphere_map.color_at_uv(u, v);
            return Some(color);
        }

        None
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
            y = rk4(&y, t, self.step_size, &self.geometry);

            // Check if there is a big jump. This happens when crossing the horizon and is a
            // heuristic here to mark this ray as entering the black hole.
            // TODO: find a better way.
            let position_jump = get_position(&(last_y - y), self.geometry.coordinate_system());
            if position_jump.get_spatial_vector().norm() > 1.0 {
                return Color::new(0, 0, 0);
            }

            match self.intersects_with_disk(
                &get_position(&last_y, self.geometry.coordinate_system()),
                &get_position(&y, self.geometry.coordinate_system()),
            ) {
                None => {}
                Some(c) => {
                    return c;
                }
            }

            match self.intersects_with_sphere(
                &get_position(&last_y, self.geometry.coordinate_system()),
                &get_position(&y, self.geometry.coordinate_system()),
            ) {
                None => {}
                Some(c) => {
                    return c;
                }
            }

            t += self.step_size;

            // iterate until the celestial plane distance has been reached.
            if radial_distance_spatial_part_squared(&get_position(
                &y,
                self.geometry.coordinate_system(),
            )) > self.max_radius_sq
            {
                let point_on_celestial_sphere =
                    get_position(&y, self.geometry.coordinate_system()).get_as_spherical();
                let u = (PI + point_on_celestial_sphere[2]) / (2.0 * PI);
                let v = point_on_celestial_sphere[1] / PI;
                return self.celestial_map.color_at_uv(u, v);
            }
        }

        Color::new(0, 0, 0)
    }
}

#[cfg(test)]
mod tests {
    use crate::camera::Camera;
    use crate::euclidean::EuclideanSpace;
    use crate::euclidean_spherical::EuclideanSpaceSpherical;
    use crate::geometry::Geometry;
    use crate::scene::{CheckerMapper, Color, Scene};
    use crate::schwarzschild::Schwarzschild;
    use crate::spherical_coordinates_helper::cartesian_to_spherical;
    use nalgebra::Vector4;

    #[test]
    fn test_color_of_ray_hits_sphere() {
        let camera = Camera::new(
            Vector4::new(0.0, 0.0, 0.0, -10.0),
            std::f64::consts::PI / 2.0,
            11,
            11,
            EuclideanSpace::new(),
        );
        let scene = create_scene(2.0, 0.2, 0.3, EuclideanSpace::new());

        let ray = camera.get_ray_for(6, 6);
        let color = scene.color_of_ray(&ray);

        assert_eq!(color, Color::new(100, 0, 0));
    }

    #[test]
    fn test_color_of_ray_hits_sphere_spherical() {
        let spatial_position = cartesian_to_spherical(&Vector4::new(0.0, 0.0, 0.0, -10.0));
        let camera = Camera::new(
            Vector4::new(
                0.0,
                spatial_position[0],
                spatial_position[1],
                spatial_position[2],
            ),
            std::f64::consts::PI / 2.0,
            11,
            11,
            EuclideanSpaceSpherical::new(),
        );
        let scene = create_scene(2.0, 0.2, 0.3, EuclideanSpaceSpherical::new());

        let ray = camera.get_ray_for(6, 6);
        let color = scene.color_of_ray(&ray);

        assert_eq!(color, Color::new(100, 0, 0));
    }

    // TODO: Euclidean -> Schwarzschild
    #[test]
    fn test_color_of_ray_hits_sphere_schwarzschild() {
        let spatial_position = cartesian_to_spherical(&Vector4::new(0.0, 0.0, 0.0, -10.0));
        let camera = Camera::new(
            Vector4::new(
                0.0,
                spatial_position[0],
                spatial_position[1],
                spatial_position[2],
            ),
            std::f64::consts::PI / 2.0,
            11,
            11,
            Schwarzschild::new(1.0),
        );
        let scene = create_scene(2.0, 0.2, 0.3, EuclideanSpaceSpherical::new());

        let ray = camera.get_ray_for(6, 6);
        let color = scene.color_of_ray(&ray);

        assert_eq!(color, Color::new(100, 0, 0));
    }

    #[test]
    fn test_color_of_ray_misses_sphere() {
        let camera = Camera::new(
            Vector4::new(0.0, 0.0, 0.0, -10.0),
            std::f64::consts::PI / 2.0,
            11,
            11,
            EuclideanSpace::new(),
        );
        let scene: Scene<CheckerMapper, EuclideanSpace> =
            create_scene(2.0, 0.2, 0.3, EuclideanSpace::new());

        let ray = camera.get_ray_for(0, 0);
        let color = scene.color_of_ray(&ray);

        assert_eq!(color, Color::new(0, 100, 0));
    }

    #[test]
    fn test_color_of_ray_misses_sphere_schwarzschild() {
        let spatial_position = cartesian_to_spherical(&Vector4::new(0.0, 0.0, 0.0, -10.0));

        let camera = Camera::new(
            Vector4::new(
                0.0,
                spatial_position[0],
                spatial_position[1],
                spatial_position[2],
            ),
            std::f64::consts::PI / 2.0,
            11,
            11,
            Schwarzschild::new(2.0),
        );
        let scene = create_scene(2.0, 0.2, 0.3, Schwarzschild::new(2.0));

        let ray = camera.get_ray_for(0, 0);
        let color = scene.color_of_ray(&ray);

        assert_eq!(color, Color::new(0, 255, 0));
    }

    #[test]
    fn test_color_of_ray_hits_horizon_schwarzschild() {
        let spatial_position = cartesian_to_spherical(&Vector4::new(0.0, 0.0, 0.0, -10.0));

        let camera = Camera::new(
            Vector4::new(
                0.0,
                spatial_position[0],
                spatial_position[1],
                spatial_position[2],
            ),
            std::f64::consts::PI / 2.0,
            11,
            11,
            Schwarzschild::new(2.0),
        );
        let scene = create_scene(2.0, 0.2, 0.3, Schwarzschild::new(2.0));

        let ray = camera.get_ray_for(6, 6);
        let color = scene.color_of_ray(&ray);

        assert_eq!(color, Color::new(0, 0, 0));
    }

    #[test]
    fn test_intersects_with_disk() {
        let camera = Camera::new(
            Vector4::new(0.0, 0.0, 0.8, -7.0),
            std::f64::consts::PI / 4.0,
            101,
            101,
            EuclideanSpace::new(),
        );
        let scene: Scene<CheckerMapper, EuclideanSpace> =
            create_scene(1.0, 2.0, 7.0, EuclideanSpace::new());

        let ray = camera.get_ray_for(0, 51);
        let color = scene.color_of_ray(&ray);

        assert_eq!(color, Color::new(0, 0, 255));
    }

    fn create_scene<G: Geometry>(
        center_sphere_radius: f64,
        center_disk_inner_radius: f64,
        center_disk_outer_radius: f64,
        geometry: G,
    ) -> Scene<CheckerMapper, G> {
        let texture_mapper_celestial =
            CheckerMapper::new(100.0, 100.0, Color::new(0, 255, 0), Color::new(0, 100, 0));
        let texture_mapper_disk =
            CheckerMapper::new(200.0, 10.0, Color::new(0, 0, 255), Color::new(0, 0, 100));
        let texture_mapper_sphere =
            CheckerMapper::new(10.0, 10.0, Color::new(255, 0, 0), Color::new(100, 0, 0));
        let scene = Scene::new(
            10000,
            15.0,
            0.01,
            center_sphere_radius,
            center_disk_inner_radius,
            center_disk_outer_radius,
            texture_mapper_celestial,
            texture_mapper_disk,
            texture_mapper_sphere,
            geometry,
        );
        scene
    }
}

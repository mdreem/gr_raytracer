use crate::camera::{Camera, Ray};
use crate::color::{wavelength_to_rgb, Color};
use crate::four_vector::{CoordinateSystem, FourVector};
use crate::geometry::Geometry;
use crate::runge_kutta::rk4;
use crate::scene::StopReason::{CelestialSphereReached, HorizonReached};
use crate::spherical_coordinates_helper::spherical_to_cartesian;
use image::{DynamicImage, GenericImageView, ImageReader};
use nalgebra::{Const, OVector, Vector3, Vector4};
use std::f64::consts::PI;
use std::fs::File;
use std::io::Write;

pub struct Scene<T: TextureMap, G: Geometry> {
    integration_configuration: IntegrationConfiguration,
    center_sphere_radius: f64,
    center_disk_outer_radius: f64,
    center_disk_inner_radius: f64,
    celestial_map: T,
    center_disk_map: T,
    center_sphere_map: T,
    pub geometry: G,
    pub camera: Camera<G>,
    save_ray_data: bool,
}

#[derive(Debug)]
pub struct Step {
    pub y: EquationOfMotionState,
    pub t: f64,
    pub step: usize,
}

#[derive(Debug, PartialEq)]
pub enum StopReason {
    HorizonReached,
    CelestialSphereReached,
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

pub type EquationOfMotionState = OVector<f64, Const<8>>;

// TODO: replace this with a more natural function, that does not need to use cartesian coordinates.
pub fn get_position(y: &EquationOfMotionState, coordinate_system: CoordinateSystem) -> FourVector {
    match coordinate_system {
        CoordinateSystem::Cartesian => FourVector::new_cartesian(y[0], y[1], y[2], y[3]),
        CoordinateSystem::Spherical => {
            let t_vec = spherical_to_cartesian(&Vector4::new(y[0], y[1], y[2], y[3]));
            FourVector::new_cartesian(t_vec[0], t_vec[1], t_vec[2], t_vec[3]) // TODO: try to do all of this without this conversion.
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

pub struct IntegrationConfiguration {
    max_steps: usize,
    max_radius_sq: f64,
    step_size: f64,
    max_steps_celestial_continuation: usize,
    max_radius_celestial_continuation_sq: f64,
    step_size_celestial_continuation: f64,
}

impl IntegrationConfiguration {
    pub fn new(
        max_steps: usize,
        max_radius: f64,
        step_size: f64,
        max_steps_celestial_continuation: usize,
        max_radius_celestial_continuation: f64,
        step_size_celestial_continuation: f64,
    ) -> IntegrationConfiguration {
        IntegrationConfiguration {
            max_steps,
            max_radius_sq: max_radius * max_radius,
            step_size,
            max_steps_celestial_continuation,
            max_radius_celestial_continuation_sq: max_radius_celestial_continuation
                * max_radius_celestial_continuation,
            step_size_celestial_continuation,
        }
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
        integration_configuration: IntegrationConfiguration,
        center_sphere_radius: f64,
        center_disk_inner_radius: f64,
        center_disk_outer_radius: f64,
        celestial_map: T,
        center_disk_map: T,
        center_sphere_map: T,
        geometry: G,
        camera: Camera<G>,
        save_ray_data: bool,
    ) -> Scene<T, G> {
        Scene {
            integration_configuration,
            center_sphere_radius,
            center_disk_outer_radius,
            center_disk_inner_radius,
            celestial_map,
            center_disk_map,
            center_sphere_map,
            geometry,
            camera,
            save_ray_data,
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

            let theta = point_on_sphere[1];
            let phi = point_on_sphere[2];
            let u = (PI + phi) / (2.0 * PI);
            let v = theta / PI;

            let color = self.center_sphere_map.color_at_uv(u, v);
            return Some(color);
        }

        None
    }

    fn should_stop(
        &self,
        last_y: &EquationOfMotionState,
        cur_y: &EquationOfMotionState,
    ) -> Option<StopReason> {
        let position = get_position(&cur_y, self.geometry.coordinate_system());

        // Check if there is a big jump. This happens when crossing the horizon and is a
        // heuristic here to mark this ray as entering the black hole.
        // TODO: find a better way.
        let position_jump = get_position(&(last_y - cur_y), self.geometry.coordinate_system());
        if position_jump.get_spatial_vector().norm() > 1.0 {
            return Some(HorizonReached);
        }

        if self
            .geometry
            .inside_horizon(&Vector4::new(cur_y[0], cur_y[1], cur_y[2], cur_y[3]))
        {
            return Some(HorizonReached);
        }

        // iterate until the celestial plane distance has been reached.
        if radial_distance_spatial_part_squared(&get_position(
            &cur_y,
            self.geometry.coordinate_system(),
        )) > self.integration_configuration.max_radius_sq
        {
            return Some(CelestialSphereReached);
        }

        None
    }

    pub fn integrate(&self, ray: &Ray) -> (Vec<Step>, Option<StopReason>) {
        let mut t = 0.0;
        let direction = ray.momentum.get_as_vector();
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

        let mut result: Vec<Step> = Vec::new();
        result.push(Step { y, t, step: 0 });

        for i in 1..self.integration_configuration.max_steps {
            let last_y = y;

            y = rk4(
                &y,
                t,
                self.integration_configuration.step_size,
                &self.geometry,
            );
            t += self.integration_configuration.step_size;

            match self.should_stop(&last_y, &y) {
                None => {}
                Some(r) => return (result, Some(r)),
            }

            result.push(Step { y, t, step: i });
        }

        (result, None)
    }

    fn save_steps(&self, steps: &Vec<Step>, filename: String) {
        let mut file = File::create(filename).expect("Unable to create file");

        file.write_all(b"i,t,tau,x,y,z\n")
            .expect("Unable to write file");

        for step in steps {
            let position = get_position(&step.y, self.geometry.coordinate_system()).get_as_vector();

            file.write_all(
                format!(
                    "{},{},{},{},{},{}\n",
                    step.step, step.t, position[0], position[1], position[2], position[3],
                )
                .as_bytes(),
            )
            .expect("Unable to write file");
        }
    }

    pub fn color_of_ray(&self, ray: &Ray) -> (Color, Option<f64>) {
        let (steps, stop_reason) = self.integrate(&ray);
        let mut y = steps[0].y;

        if self.save_ray_data {
            self.save_steps(
                &steps,
                String::from(format!("ray-{}-{}.csv", ray.row, ray.col)),
            );
        }

        let velocity = self.camera.velocity;
        let observer_energy = self.geometry.mul(&ray.position, &velocity, &ray.momentum);

        for step in steps.iter().skip(1) {
            let last_y = y;
            y = step.y;

            match self.intersects_with_disk(
                &get_position(&last_y, self.geometry.coordinate_system()),
                &get_position(&y, self.geometry.coordinate_system()),
            ) {
                None => {}
                Some(c) => {
                    let redshift = self.compute_redshift(y, observer_energy);
                    let tune_redshift = 1.0;
                    let redshift = (redshift - 1.0) * tune_redshift + 1.0;
                    let wavelength = 400.0 * redshift;
                    let c = wavelength_to_rgb(wavelength);
                    return (c, Some(redshift));
                }
            }

            match self.intersects_with_sphere(
                &get_position(&last_y, self.geometry.coordinate_system()),
                &get_position(&y, self.geometry.coordinate_system()),
            ) {
                None => {}
                Some(c) => {
                    let redshift = self.compute_redshift(y, observer_energy);
                    return (c, Some(redshift));
                }
            }
        }

        if let Some(reason) = stop_reason {
            match reason {
                HorizonReached => {
                    return (Color::new(0, 0, 0), None);
                }
                CelestialSphereReached => {
                    let mut y = steps.last().unwrap().y;
                    let mut t = steps.last().unwrap().t;
                    let step_size = self
                        .integration_configuration
                        .step_size_celestial_continuation;

                    // integrate further until we are far out.
                    for _ in 1..self
                        .integration_configuration
                        .max_steps_celestial_continuation
                    {
                        y = rk4(&y, t, step_size, &self.geometry);
                        t += step_size;

                        if radial_distance_spatial_part_squared(&get_position(
                            &y,
                            self.geometry.coordinate_system(),
                        )) > self
                            .integration_configuration
                            .max_radius_celestial_continuation_sq
                        {
                            break;
                        }
                    }

                    let point_on_celestial_sphere =
                        get_position(&y, self.geometry.coordinate_system()).get_as_spherical();
                    let theta = point_on_celestial_sphere[1];
                    let phi = point_on_celestial_sphere[2];
                    let u = (PI + phi) / (2.0 * PI);
                    let v = theta / PI;

                    let redshift = self.compute_redshift(y, observer_energy);
                    return (
                        self.celestial_map.color_at_uv(1.0 - u, 1.0 - v),
                        Some(redshift),
                    );
                }
            }
        }
        (Color::new(0, 0, 0), None)
    }

    fn compute_redshift(&self, y: EquationOfMotionState, observer_energy: f64) -> f64 {
        let emitter_energy = self.energy_of_stationary_emitter(y);
        let shift = emitter_energy / observer_energy;
        shift
    }

    fn energy_of_stationary_emitter(&self, y: EquationOfMotionState) -> f64 {
        let position = Vector4::new(y[0], y[1], y[2], y[3]);
        let velocity = self.geometry.get_stationary_velocity_at(&position);
        let momentum = FourVector::new(y[4], y[5], y[6], y[7], self.geometry.coordinate_system());
        let energy = self.geometry.mul(&position, &velocity, &momentum);
        energy
    }
}

#[cfg(test)]
pub mod test_scene {
    use crate::camera::Camera;
    use crate::color::Color;
    use crate::four_vector::FourVector;
    use crate::geometry::Geometry;
    use crate::scene::{CheckerMapper, IntegrationConfiguration, Scene};
    use nalgebra::Vector4;
    use std::f64::consts::PI;

    pub const CELESTIAL_SPHERE_RADIUS: f64 = 15.0;

    pub fn create_scene<G: Geometry>(
        center_sphere_radius: f64,
        center_disk_inner_radius: f64,
        center_disk_outer_radius: f64,
        geometry: G,
        camera_position: Vector4<f64>,
        camera_velocity: FourVector,
    ) -> Scene<CheckerMapper, G> {
        let camera = Camera::new(
            camera_position,
            camera_velocity,
            PI / 4.0,
            500,
            500,
            geometry.clone(), // TODO see how geometry can be distributed to all needed places.
        );

        create_scene_with_camera(
            center_sphere_radius,
            center_disk_inner_radius,
            center_disk_outer_radius,
            geometry,
            camera,
        )
    }

    pub fn create_scene_with_camera<G: Geometry>(
        center_sphere_radius: f64,
        center_disk_inner_radius: f64,
        center_disk_outer_radius: f64,
        geometry: G,
        camera: Camera<G>,
    ) -> Scene<CheckerMapper, G> {
        let texture_mapper_celestial =
            CheckerMapper::new(100.0, 100.0, Color::new(0, 255, 0), Color::new(0, 100, 0));
        let texture_mapper_disk =
            CheckerMapper::new(200.0, 10.0, Color::new(0, 0, 255), Color::new(0, 0, 100));
        let texture_mapper_sphere =
            CheckerMapper::new(10.0, 10.0, Color::new(255, 0, 0), Color::new(100, 0, 0));

        let integration_configuration = IntegrationConfiguration::new(
            30000,
            CELESTIAL_SPHERE_RADIUS,
            0.001,
            15000,
            10000.0,
            1.0,
        );

        let scene = Scene::new(
            integration_configuration,
            center_sphere_radius,
            center_disk_inner_radius,
            center_disk_outer_radius,
            texture_mapper_celestial,
            texture_mapper_disk,
            texture_mapper_sphere,
            geometry,
            camera,
            false,
        );
        scene
    }
}

#[cfg(test)]
mod tests {
    use crate::camera::Camera;
    use crate::euclidean::EuclideanSpace;
    use crate::euclidean_spherical::EuclideanSpaceSpherical;
    use crate::four_vector::FourVector;
    use crate::scene::test_scene::create_scene_with_camera;
    use crate::scene::{CheckerMapper, Color, Scene};
    use crate::schwarzschild::Schwarzschild;
    use crate::spherical_coordinates_helper::cartesian_to_spherical;
    use approx::assert_abs_diff_eq;
    use nalgebra::{ComplexField, Vector4};
    use std::f64::consts::PI;

    #[test]
    fn test_color_of_ray_hits_sphere() {
        let camera = Camera::new(
            Vector4::new(0.0, 0.0, 0.0, -10.0),
            FourVector::new_cartesian(1.0, 0.0, 0.0, 0.0),
            std::f64::consts::PI / 2.0,
            11,
            11,
            EuclideanSpace::new(),
        );
        let scene = create_scene_with_camera(2.0, 0.2, 0.3, EuclideanSpace::new(), camera);

        let ray = scene.camera.get_ray_for(6, 6);
        let (color, _) = scene.color_of_ray(&ray);

        assert_eq!(color, Color::new(100, 0, 0));
    }

    #[test]
    fn test_color_of_ray_hits_sphere_spherical() {
        let position = cartesian_to_spherical(&Vector4::new(0.0, 0.0, 0.0, -10.0));
        let camera = Camera::new(
            position,
            FourVector::new_spherical(1.0, 0.0, 0.0, 0.0),
            std::f64::consts::PI / 2.0,
            11,
            11,
            EuclideanSpaceSpherical::new(),
        );
        let scene = create_scene_with_camera(2.0, 0.2, 0.3, EuclideanSpaceSpherical::new(), camera);

        let ray = scene.camera.get_ray_for(6, 6);
        let (color, _) = scene.color_of_ray(&ray);

        assert_eq!(color, Color::new(100, 0, 0));
    }

    #[test]
    fn test_color_of_ray_hits_sphere_schwarzschild() {
        let position = Vector4::new(0.0, 10.0, PI / 2.0, 0.0);
        let radius = 1.0;
        let r = position[1];
        let a = 1.0 - radius / r;
        let velocity = FourVector::new_spherical(1.0 / a, -(radius / r).sqrt(), 0.0, 0.0); // we have a freely falling observer here.

        let geometry = Schwarzschild::new(radius);

        let camera = Camera::new(position, velocity, PI / 2.0, 11, 11, geometry.clone());
        let scene = create_scene_with_camera(2.0, 0.2, 0.3, geometry, camera);

        let ray = scene.camera.get_ray_for(6, 6);
        let (color, Some(redshift)) = scene.color_of_ray(&ray) else {
            panic!("No redshift found");
        };

        assert_eq!(color, Color::new(255, 0, 0));
    }

    #[test]
    fn test_color_of_ray_hits_sphere_schwarzschild_stationary_observer() {
        let position = Vector4::new(0.0, 10.0, PI / 2.0, 0.0);
        let radius = 1.0;
        let sphere_radius = 2.0;
        let r = position[1];
        let a = 1.0 - radius / r;
        let velocity = FourVector::new_spherical(a.sqrt().recip(), 0.0, 0.0, 0.0); // we have a freely falling observer here.

        let geometry = Schwarzschild::new(radius);

        let camera = Camera::new(position, velocity, PI / 2.0, 11, 11, geometry.clone());
        let scene = create_scene_with_camera(sphere_radius, 0.2, 0.3, geometry, camera);

        let ray = scene.camera.get_ray_for(6, 6);
        let (color, Some(redshift)) = scene.color_of_ray(&ray) else {
            panic!("No redshift found");
        };

        let a_emitter = 1.0 - radius / sphere_radius;
        let expected_redshift = (a / a_emitter).sqrt();

        assert_abs_diff_eq!(redshift, expected_redshift, epsilon = 1e-3);
        assert_eq!(color, Color::new(255, 0, 0));
    }

    #[test]
    fn test_color_of_ray_misses_sphere() {
        let camera = Camera::new(
            Vector4::new(0.0, 0.0, 0.0, -10.0),
            FourVector::new_spherical(1.0, 0.0, 0.0, 0.0),
            PI / 2.0,
            11,
            11,
            EuclideanSpace::new(),
        );
        let scene: Scene<CheckerMapper, EuclideanSpace> =
            create_scene_with_camera(2.0, 0.2, 0.3, EuclideanSpace::new(), camera);

        let ray = scene.camera.get_ray_for(0, 0);
        let (color, _) = scene.color_of_ray(&ray);

        assert_eq!(color, Color::new(0, 255, 0));
    }

    #[test]
    fn test_color_of_ray_misses_sphere_schwarzschild() {
        let position = cartesian_to_spherical(&Vector4::new(0.0, 0.0, 0.0, -10.0));
        let radius = 2.0;
        let r = position[1];
        let a = 1.0 - radius / r;
        let velocity = FourVector::new_spherical(1.0 / a, -(radius / r).sqrt(), 0.0, 0.0); // we have a freely falling observer here.

        let camera = Camera::new(
            position,
            velocity,
            PI / 2.0,
            11,
            11,
            Schwarzschild::new(radius),
        );
        let scene = create_scene_with_camera(2.0, 0.2, 0.3, Schwarzschild::new(radius), camera);

        let ray = scene.camera.get_ray_for(0, 0);
        let (color, _) = scene.color_of_ray(&ray);

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
            FourVector::new_spherical(1.0, 0.0, 0.0, 0.0),
            PI / 2.0,
            11,
            11,
            Schwarzschild::new(2.0),
        );
        let scene = create_scene_with_camera(2.0, 0.2, 0.3, Schwarzschild::new(2.0), camera);

        let ray = scene.camera.get_ray_for(6, 6);
        let (color, _) = scene.color_of_ray(&ray);

        assert_eq!(color, Color::new(0, 0, 0));
    }

    #[test]
    fn test_intersects_with_disk() {
        let camera = Camera::new(
            Vector4::new(0.0, 0.0, 0.8, -7.0),
            FourVector::new_cartesian(1.0, 0.0, 0.0, 0.0),
            PI / 4.0,
            101,
            101,
            EuclideanSpace::new(),
        );
        let scene: Scene<CheckerMapper, EuclideanSpace> =
            create_scene_with_camera(1.0, 2.0, 7.0, EuclideanSpace::new(), camera);

        let ray = scene.camera.get_ray_for(0, 51);
        let (color, redshift) = scene.color_of_ray(&ray);

        assert_eq!(redshift, Some(1.0));
        assert_eq!(color, Color::new(1, 0, 16));
    }
}

use crate::geometry::geometry::Geometry;
use crate::geometry::point::{CoordinateSystem, Point};
use crate::geometry::spherical_coordinates_helper::spherical_to_cartesian;
use crate::rendering::camera::Camera;
use crate::rendering::color::CIETristimulus;
use crate::rendering::integrator::StopReason::{CelestialSphereReached, HorizonReached};
use crate::rendering::integrator::{IntegrationConfiguration, Integrator, StopReason};
use crate::rendering::ray::{IntegratedRay, Ray};
use crate::rendering::redshift::RedshiftComputer;
use crate::rendering::texture::{TextureData, UVCoordinates};
use crate::scene_objects::objects::Objects;
use nalgebra::{Const, OVector};
use std::f64::consts::PI;
use std::fs::File;

pub struct Scene<'a, G: Geometry> {
    pub integrator: Integrator<'a, G>,
    objects: Objects<'a, G>,
    texture_data: TextureData,
    pub geometry: &'a G,
    pub camera: Camera,
    save_ray_data: bool,
    redshift_computer: RedshiftComputer<'a, G>,
}

pub type EquationOfMotionState = OVector<f64, Const<8>>;

// TODO: replace this with a more natural function, that does not need to use cartesian coordinates.
pub fn get_position(y: &EquationOfMotionState, coordinate_system: CoordinateSystem) -> Point {
    match coordinate_system {
        CoordinateSystem::Cartesian => Point::new_cartesian(y[0], y[1], y[2], y[3]),
        CoordinateSystem::Spherical => {
            let t_vec = spherical_to_cartesian(&Point::new(
                y[0],
                y[1],
                y[2],
                y[3],
                CoordinateSystem::Spherical,
            ));
            Point::new_cartesian(t_vec[0], t_vec[1], t_vec[2], t_vec[3]) // TODO: try to do all of this without this conversion.
        }
    }
}

impl<'a, G: Geometry> Scene<'a, G> {
    pub fn new(
        integration_configuration: IntegrationConfiguration,
        objects: Objects<'a, G>,
        texture_data: TextureData,
        geometry: &'a G,
        camera: Camera,
        save_ray_data: bool,
    ) -> Scene<'a, G> {
        let integrator = Integrator::new(geometry, integration_configuration);

        Scene {
            integrator,
            objects,
            texture_data,
            geometry,
            camera,
            save_ray_data,
            redshift_computer: RedshiftComputer::new(geometry),
        }
    }

    pub fn integrate_ray(&self, ray: &Ray) -> (IntegratedRay, Option<StopReason>) {
        self.integrator.integrate(ray)
    }

    pub fn color_of_ray(&self, ray: &Ray) -> CIETristimulus {
        let (steps, stop_reason) = self.integrate_ray(ray);
        let mut y = steps[0].y;

        if self.save_ray_data {
            let mut file = File::create(format!("ray-{}-{}.csv", ray.row, ray.col))
                .expect("Unable to create file");
            steps.save(&mut file, self.geometry);
        }

        let velocity = self.camera.velocity;
        let observer_energy = self.redshift_computer.get_observer_energy(ray, &velocity);

        let mut intersections = Vec::new();
        for step in steps.iter().skip(1) {
            let last_y = y;
            y = step.y;

            if let Some(intersection_color) =
                self.objects
                    .intersects(&last_y, &y, observer_energy, &self.redshift_computer)
            {
                intersections.push(intersection_color);
            }
        }

        if let Some(reason) = stop_reason {
            match reason {
                HorizonReached => {
                    intersections.push(CIETristimulus::new(0.0, 0.0, 0.0, 1.0));
                }
                CelestialSphereReached => {
                    let uv = self.get_uv_coordinates(&y);
                    let redshift = self.redshift_computer.compute_redshift(&y, observer_energy);
                    intersections.push(
                        self.texture_data
                            .celestial_map
                            .color_at_uv(uv, 2500.0, redshift),
                    );
                }
            };
        } else {
            println!("ERROR: Ray did not hit anything: {:?}", ray);
        }
        let mut result = CIETristimulus::new(0.0, 0.0, 0.0, 1.0);

        for color in intersections.iter().rev() {
            result = result.blend(&color)
        }

        result
    }

    fn get_uv_coordinates(&self, y_far: &EquationOfMotionState) -> UVCoordinates {
        let point_on_celestial_sphere =
            get_position(y_far, self.geometry.coordinate_system()).get_as_spherical();

        let theta = point_on_celestial_sphere[1];
        let phi = point_on_celestial_sphere[2];

        let u = (PI + phi) / (2.0 * PI);
        let v = theta / PI;
        UVCoordinates {
            u: 1.0 - u,
            v: 1.0 - v,
        }
    }
}

#[cfg(test)]
pub mod test_scene {
    use crate::geometry::four_vector::FourVector;
    use crate::geometry::geometry::Geometry;
    use crate::geometry::point::Point;
    use crate::rendering::camera::Camera;
    use crate::rendering::color::Color;
    use crate::rendering::scene::{IntegrationConfiguration, Scene, TextureData};
    use crate::rendering::texture::CheckerMapper;
    use crate::scene_objects;
    use std::f64::consts::PI;
    use std::sync::Arc;

    pub const CELESTIAL_SPHERE_RADIUS: f64 = 10000.0;

    pub fn create_scene<G: Geometry>(
        center_sphere_radius: f64,
        center_disk_inner_radius: f64,
        center_disk_outer_radius: f64,
        geometry: &G,
        camera_position: Point,
        camera_velocity: FourVector,
    ) -> Scene<'_, G> {
        let camera = Camera::new(
            camera_position,
            camera_velocity,
            PI / 4.0,
            500,
            500,
            geometry,
        );

        create_scene_with_camera(
            center_sphere_radius,
            center_disk_inner_radius,
            center_disk_outer_radius,
            geometry,
            camera,
            1e-12,
        )
    }

    pub fn create_scene_with_camera<G: Geometry>(
        center_sphere_radius: f64,
        center_disk_inner_radius: f64,
        center_disk_outer_radius: f64,
        geometry: &G,
        camera: Camera,
        epsilon: f64,
    ) -> Scene<'_, G> {
        let texture_mapper_celestial = Arc::new(CheckerMapper::new(
            100.0,
            100.0,
            Color::new(0, 255, 0, 255),
            Color::new(0, 100, 0, 255),
        ));
        let texture_mapper_disk = Arc::new(CheckerMapper::new(
            200.0,
            10.0,
            Color::new(0, 0, 255, 255),
            Color::new(0, 0, 100, 255),
        ));
        let texture_mapper_sphere = Arc::new(CheckerMapper::new(
            10.0,
            10.0,
            Color::new(255, 0, 0, 255),
            Color::new(100, 0, 0, 255),
        ));

        let integration_configuration =
            IntegrationConfiguration::new(30000, CELESTIAL_SPHERE_RADIUS, 0.001, epsilon);

        let texture_data = TextureData {
            celestial_map: texture_mapper_celestial,
        };
        let mut objects = scene_objects::objects::Objects::new(geometry);
        objects.add_object(Box::new(scene_objects::sphere::Sphere::new(
            center_sphere_radius,
            texture_mapper_sphere,
            Point::new_cartesian(0.0, 0.0, 0.0, 0.0),
        )));
        objects.add_object(Box::new(scene_objects::disc::Disc::new(
            center_disk_inner_radius,
            center_disk_outer_radius,
            texture_mapper_disk,
        )));

        let scene = Scene::new(
            integration_configuration,
            objects,
            texture_data,
            geometry,
            camera,
            false,
        );
        scene
    }
}

#[cfg(test)]
mod tests {
    use crate::geometry::euclidean::EuclideanSpace;
    use crate::geometry::euclidean_spherical::EuclideanSpaceSpherical;
    use crate::geometry::four_vector::FourVector;
    use crate::geometry::point::{CoordinateSystem, Point};
    use crate::geometry::schwarzschild::Schwarzschild;
    use crate::geometry::spherical_coordinates_helper::cartesian_to_spherical;
    use crate::rendering::camera::Camera;
    use crate::rendering::color::CIETristimulus;
    use crate::rendering::scene::test_scene::create_scene_with_camera;
    use crate::rendering::scene::Scene;
    use std::f64::consts::PI;

    #[test]
    fn test_color_of_ray_hits_sphere() {
        let camera = Camera::new(
            Point::new(0.0, 0.0, 0.0, -10.0, CoordinateSystem::Cartesian),
            FourVector::new_cartesian(1.0, 0.0, 0.0, 0.0),
            PI / 2.0,
            11,
            11,
            &EuclideanSpace::new(),
        );
        let space = EuclideanSpace::new();
        let scene = create_scene_with_camera(2.0, 0.2, 0.3, &space, camera, 1e-12);

        let ray = scene.camera.get_ray_for(5, 5);
        let color = scene.color_of_ray(&ray);

        assert_eq!(
            color,
            CIETristimulus::new(
                0.27813592476573396,
                0.28536499372678936,
                0.4364990815074766,
                1.0
            )
        );
    }

    #[test]
    fn test_color_of_ray_hits_sphere_spherical() {
        let position = cartesian_to_spherical(&Point::new(
            0.0,
            0.0,
            0.0,
            -10.0,
            CoordinateSystem::Spherical,
        ));
        let camera = Camera::new(
            position,
            FourVector::new_spherical(1.0, 0.0, 0.0, 0.0),
            PI / 2.0,
            11,
            11,
            &EuclideanSpaceSpherical::new(),
        );
        let space = EuclideanSpaceSpherical::new();
        let scene = create_scene_with_camera(2.0, 0.2, 0.3, &space, camera, 1e-12);

        let ray = scene.camera.get_ray_for(5, 5);
        let color = scene.color_of_ray(&ray);

        assert_eq!(
            color,
            CIETristimulus::new(
                0.27813592476573396,
                0.28536499372678936,
                0.4364990815074766,
                1.0
            )
        );
    }

    #[test]
    fn test_color_of_ray_hits_sphere_schwarzschild() {
        let position = Point::new(0.0, 10.0, PI / 2.0, 0.0, CoordinateSystem::Spherical);
        let radius = 1.0;
        let r = position[1];
        let a = 1.0 - radius / r;
        let velocity = FourVector::new_spherical(1.0 / a, -(radius / r).sqrt(), 0.0, 0.0); // we have a freely falling observer here.

        let geometry = Schwarzschild::new(radius, 1e-4);

        let camera = Camera::new(position, velocity, PI / 2.0, 11, 11, &geometry);
        let scene = create_scene_with_camera(2.0, 0.2, 0.3, &geometry, camera, 1e-12);

        let ray = scene.camera.get_ray_for(5, 5);
        let color = scene.color_of_ray(&ray);
        assert_eq!(
            color,
            CIETristimulus::new(
                0.32882834826291946,
                0.3383856312341892,
                0.33278602050289136,
                1.0
            )
        );
    }

    #[test]
    fn test_color_of_ray_hits_sphere_schwarzschild_stationary_observer() {
        let position = Point::new(0.0, 10.0, PI / 2.0, 0.0, CoordinateSystem::Spherical);
        let radius = 1.0;
        let sphere_radius = 2.0;
        let r = position[1];
        let a = 1.0 - radius / r;
        let velocity = FourVector::new_spherical(a.sqrt().recip(), 0.0, 0.0, 0.0); // we have a freely falling observer here.

        let geometry = Schwarzschild::new(radius, 1e-4);

        let camera = Camera::new(position, velocity, PI / 2.0, 11, 11, &geometry);
        let scene = create_scene_with_camera(sphere_radius, 0.2, 0.3, &geometry, camera, 1e-12);

        let ray = scene.camera.get_ray_for(5, 5);
        let color = scene.color_of_ray(&ray);
        let a_emitter = 1.0 - radius / sphere_radius;
        let expected_redshift = (a / a_emitter).sqrt();

        assert_eq!(
            color,
            CIETristimulus::new(
                0.29677218433467684,
                0.3066953444984932,
                0.3965324711668299,
                1.0
            )
        );
    }

    #[test]
    fn test_color_of_ray_misses_sphere() {
        let camera = Camera::new(
            Point::new(0.0, 0.0, 0.0, -10.0, CoordinateSystem::Cartesian),
            FourVector::new_spherical(1.0, 0.0, 0.0, 0.0),
            PI / 2.0,
            11,
            11,
            &EuclideanSpace::new(),
        );
        let space = EuclideanSpace::new();
        let scene: Scene<EuclideanSpace> =
            create_scene_with_camera(2.0, 0.2, 0.3, &space, camera, 1e-12);

        let ray = scene.camera.get_ray_for(0, 0);
        let color = scene.color_of_ray(&ray);

        assert_eq!(
            color,
            CIETristimulus::new(
                0.04556866876322511,
                0.09113733752645022,
                0.015189552006485689,
                1.0
            )
        );
    }

    #[test]
    fn test_color_of_ray_misses_sphere_schwarzschild() {
        let position = cartesian_to_spherical(&Point::new(
            0.0,
            0.0,
            0.0,
            -10.0,
            CoordinateSystem::Cartesian,
        ));
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
            &Schwarzschild::new(radius, 1e-4),
        );
        let space = Schwarzschild::new(radius, 1e-4);
        let scene = create_scene_with_camera(2.0, 0.2, 0.3, &space, camera, 1e-12);

        let ray = scene.camera.get_ray_for(0, 0);
        let color = scene.color_of_ray(&ray);

        assert_eq!(
            color,
            CIETristimulus::new(0.3575761, 0.7151522, 0.119192, 1.0)
        );
    }

    #[test]
    fn test_color_of_ray_hits_horizon_schwarzschild() {
        let spatial_position = cartesian_to_spherical(&Point::new(
            0.0,
            0.0,
            0.0,
            -10.0,
            CoordinateSystem::Cartesian,
        ));

        let camera = Camera::new(
            Point::new(
                0.0,
                spatial_position[0],
                spatial_position[1],
                spatial_position[2],
                CoordinateSystem::Cartesian,
            ),
            FourVector::new_spherical(1.0, 0.0, 0.0, 0.0),
            PI / 2.0,
            11,
            11,
            &Schwarzschild::new(2.0, 1e-4),
        );
        let space = Schwarzschild::new(2.0, 1e-4);
        let scene = create_scene_with_camera(2.0, 0.2, 0.3, &space, camera, 1e-12);

        let ray = scene.camera.get_ray_for(6, 6);
        let color = scene.color_of_ray(&ray);

        assert_eq!(color, CIETristimulus::new(0.0, 0.0, 0.0, 1.0));
    }

    #[test]
    fn test_intersects_with_disk() {
        let camera = Camera::new(
            Point::new(0.0, 0.0, 0.8, -7.0, CoordinateSystem::Cartesian),
            FourVector::new_cartesian(1.0, 0.0, 0.0, 0.0),
            PI / 4.0,
            101,
            101,
            &EuclideanSpace::new(),
        );
        let space = EuclideanSpace::new();
        let scene: Scene<EuclideanSpace> =
            create_scene_with_camera(1.0, 2.0, 7.0, &space, camera, 1e-12);

        let ray = scene.camera.get_ray_for(0, 51);
        let color = scene.color_of_ray(&ray);

        assert_eq!(
            color,
            CIETristimulus::new(
                0.27813592476573396,
                0.28536499372678936,
                0.4364990815074766,
                1.0
            )
        );
    }
}

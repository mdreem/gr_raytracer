use crate::cli::cli::GlobalOpts;
use crate::configuration::RenderConfig;
use crate::geometry::four_vector::FourVector;
use crate::geometry::geometry::Geometry;
use crate::rendering::camera::Camera;
use crate::rendering::integrator::IntegrationConfiguration;
use crate::rendering::raytracer;
use crate::rendering::scene::Scene;
use crate::rendering::texture::{TextureData, TextureMapper};
use crate::scene_objects::objects::Objects;
use crate::{configuration, scene_objects};
use nalgebra::Vector4;
use std::f64::consts::PI;

pub fn render<G: Geometry>(scene: Scene<TextureMapper, G>, filename: String) {
    let raytracer = raytracer::Raytracer::new(scene);
    raytracer.render(filename);
}

pub fn create_scene<G: Geometry>(
    geometry: &G,
    camera_position: Vector4<f64>,
    camera_momentum: FourVector,
    opts: GlobalOpts,
    config: RenderConfig,
) -> Scene<TextureMapper, G> {
    let integration_configuration = IntegrationConfiguration::new(
        opts.max_steps,
        opts.max_radius,
        opts.step_size,
        opts.max_steps_celestial_continuation,
        opts.max_radius_celestial_continuation,
        opts.step_size_celestial_continuation,
    );

    let texture_mapper_celestial = TextureMapper::new(String::from("./resources/celestial.png"));
    let texture_mapper_disk = TextureMapper::new(String::from("./resources/disk.png"));
    let texture_mapper_sphere = TextureMapper::new(String::from("./resources/sphere.png"));

    let camera = Camera::new(
        camera_position,
        camera_momentum,
        PI / 4.0,
        opts.height,
        opts.width,
        geometry,
    );

    let texture_data = TextureData {
        celestial_map: texture_mapper_celestial,
        center_disk_map: texture_mapper_disk.clone(),
        center_sphere_map: texture_mapper_sphere.clone(),
    };

    let mut objects = Objects::new();
    for object in config.objects {
        match object {
            configuration::ObjectsConfig::Sphere { radius } => {
                println!("Adding sphere with radius: {}", radius);
                let sphere =
                    scene_objects::sphere::Sphere::new(radius, texture_mapper_sphere.clone());
                objects.add_object(Box::new(sphere));
            }
            configuration::ObjectsConfig::Disc {
                inner_radius,
                outer_radius,
            } => {
                println!(
                    "Adding disc with inner radius: {}, outer radius: {}",
                    inner_radius, outer_radius
                );
                let disc = scene_objects::disc::Disc::new(
                    inner_radius,
                    outer_radius,
                    texture_mapper_disk.clone(),
                );
                objects.add_object(Box::new(disc));
            }
        }
    }

    let scene = Scene::new(
        integration_configuration,
        objects,
        texture_data,
        geometry,
        camera,
        false,
        true,
    );
    scene
}

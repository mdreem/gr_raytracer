use crate::cli::cli::GlobalOpts;
use crate::configuration::RenderConfig;
use crate::geometry::four_vector::FourVector;
use crate::geometry::geometry::Geometry;
use crate::geometry::point::Point;
use crate::rendering::camera::Camera;
use crate::rendering::color::Color;
use crate::rendering::integrator::IntegrationConfiguration;
use crate::rendering::raytracer;
use crate::rendering::scene::Scene;
use crate::rendering::texture::{CheckerMapper, TextureData, TextureMapperFactory};
use crate::scene_objects::objects::Objects;
use crate::{configuration, scene_objects};
use std::f64::consts::PI;
use std::sync::Arc;

pub fn render<G: Geometry>(scene: Scene<G>, filename: String) {
    let raytracer = raytracer::Raytracer::new(scene);
    raytracer.render(filename);
}

pub fn create_scene<G: Geometry>(
    geometry: &G,
    camera_position: Point,
    camera_momentum: FourVector,
    opts: GlobalOpts,
    config: RenderConfig,
) -> Scene<'_, G> {
    let integration_configuration = IntegrationConfiguration::new(
        opts.max_steps,
        opts.max_radius,
        opts.step_size,
        opts.epsilon,
    );

    let mut texture_mapper_factory = TextureMapperFactory::new();

    let texture_mapper_celestial = if let Some(texture) = config.celestial_texture {
        texture_mapper_factory.get_texture_mapper(texture)
    } else {
        Arc::new(CheckerMapper::new(
            20.0,
            20.0,
            Color::new(255, 255, 255, 255),
            Color::new(0, 0, 0, 255),
        ))
    };

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
    };

    let mut objects = Objects::new(geometry);
    for object in config.objects {
        match object {
            configuration::ObjectsConfig::Sphere {
                radius,
                position,
                texture,
            } => {
                println!(
                    "Adding sphere with radius: {} at ({},{},{})",
                    radius, position.0, position.1, position.2
                );

                let texture_mapper_sphere = if let Some(texture) = texture {
                    texture_mapper_factory.get_texture_mapper(texture)
                } else {
                    Arc::new(CheckerMapper::new(
                        10.0,
                        10.0,
                        Color::new(200, 0, 0, 255),
                        Color::new(0, 200, 0, 255),
                    ))
                };
                let sphere = scene_objects::sphere::Sphere::new(
                    radius,
                    texture_mapper_sphere,
                    Point::new_cartesian(0.0, position.0, position.1, position.2),
                );
                objects.add_object(Box::new(sphere));
            }
            configuration::ObjectsConfig::Disc {
                inner_radius,
                outer_radius,
                texture,
            } => {
                println!(
                    "Adding disc with inner radius: {}, outer radius: {}",
                    inner_radius, outer_radius
                );
                let texture_mapper_disc = if let Some(texture) = texture {
                    texture_mapper_factory.get_texture_mapper(texture)
                } else {
                    Arc::new(CheckerMapper::new(
                        10.0,
                        10.0,
                        Color::new(200, 200, 0, 255),
                        Color::new(0, 200, 200, 255),
                    ))
                };
                let disc =
                    scene_objects::disc::Disc::new(inner_radius, outer_radius, texture_mapper_disc);
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
    );
    scene
}

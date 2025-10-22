use crate::cli::cli::GlobalOpts;
use crate::configuration::{RenderConfig, TextureConfig};
use crate::geometry::four_vector::FourVector;
use crate::geometry::geometry::Geometry;
use crate::geometry::point::Point;
use crate::rendering::camera::Camera;
use crate::rendering::color::{CIETristimulusNormalization, Color};
use crate::rendering::integrator::IntegrationConfiguration;
use crate::rendering::raytracer;
use crate::rendering::raytracer::RaytracerError;
use crate::rendering::scene::Scene;
use crate::rendering::texture::{
    BlackBodyMapper, CheckerMapper, TextureData, TextureMapHandle, TextureMapperFactory,
};
use crate::scene_objects::objects::Objects;
use crate::{configuration, scene_objects};
use log::debug;
use std::f64::consts::PI;
use std::sync::Arc;

pub fn render<G: Geometry>(
    scene: Scene<G>,
    filename: String,
    color_normalization: CIETristimulusNormalization,
) -> Result<(), RaytracerError> {
    let raytracer = raytracer::Raytracer::new(scene, color_normalization);
    raytracer.render(filename)
}

pub fn create_scene<G: Geometry>(
    geometry: &G,
    camera_position: Point,
    camera_momentum: FourVector,
    opts: GlobalOpts,
    config: RenderConfig,
) -> Result<Scene<'_, G>, RaytracerError> {
    let integration_configuration = IntegrationConfiguration::new(
        opts.max_steps,
        opts.max_radius,
        opts.step_size,
        opts.epsilon,
    );

    let mut texture_mapper_factory = TextureMapperFactory::new();

    let texture_mapper_celestial =
        get_texture_mapper(&mut texture_mapper_factory, config.celestial_texture)?;

    let camera = Camera::new(
        camera_position,
        camera_momentum,
        PI / 4.0,
        opts.height,
        opts.width,
        opts.phi,
        opts.theta,
        opts.psi,
        geometry,
    )
    .map_err(RaytracerError::CameraError)?;

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
                debug!(
                    "Adding sphere with radius: {} at ({},{},{})",
                    radius, position.0, position.1, position.2
                );
                let texture_mapper_sphere =
                    get_texture_mapper(&mut texture_mapper_factory, texture)?;

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
                debug!(
                    "Adding disc with inner radius: {}, outer radius: {}",
                    inner_radius, outer_radius
                );
                let texture_mapper_disc = get_texture_mapper(&mut texture_mapper_factory, texture)?;
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
    Ok(scene)
}

fn get_texture_mapper(
    texture_mapper_factory: &mut TextureMapperFactory,
    texture: TextureConfig,
) -> Result<TextureMapHandle, RaytracerError> {
    let texture_mapper_sphere = match texture {
        TextureConfig::Bitmap {
            path,
            color_normalization,
        } => texture_mapper_factory.get_texture_mapper(path, color_normalization)?,
        TextureConfig::Checker {
            width,
            height,
            color1,
            color2,
            color_normalization,
        } => Arc::new(CheckerMapper::new(
            width,
            height,
            Color::new(color1.0, color1.1, color1.2, 255),
            Color::new(color2.0, color2.1, color2.2, 255),
            color_normalization,
        )),
        TextureConfig::BlackBody {
            temperature,
            color_normalization,
        } => Arc::new(BlackBodyMapper::new(temperature, color_normalization)),
    };
    Ok(texture_mapper_sphere)
}

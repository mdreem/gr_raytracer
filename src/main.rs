use crate::cli::{Action, App, GlobalOpts};
use crate::configuration::RenderConfig;
use clap::Parser;
use nalgebra::Vector4;
use rendering::camera::Camera;
use rendering::scene::Scene;
use std::f64::consts::PI;
use std::fs;
use std::time::Instant;

mod cli;
mod configuration;
mod geometry;
mod rendering;
mod scene_objects;

use crate::geometry::euclidean::EuclideanSpace;
use crate::geometry::euclidean_spherical::EuclideanSpaceSpherical;
use crate::geometry::four_vector::FourVector;
use crate::geometry::geometry::Geometry;
use crate::geometry::schwarzschild::Schwarzschild;
use crate::geometry::spherical_coordinates_helper::cartesian_to_spherical;
use crate::scene_objects::objects::Objects;
use rendering::integrator::IntegrationConfiguration;
use rendering::raytracer;
use rendering::texture::{TextureData, TextureMapper};

fn render_euclidean(
    opts: GlobalOpts,
    config: RenderConfig,
    camera_position: Vector4<f64>,
    filename: String,
) {
    let geometry = EuclideanSpace::new();
    let momentum = FourVector::new_cartesian(1.0, 0.0, 0.0, 0.0);

    render(geometry, camera_position, momentum, opts, config, filename);
}

fn render_euclidean_spherical(
    opts: GlobalOpts,
    config: RenderConfig,
    camera_position: Vector4<f64>,
    filename: String,
) {
    let camera_position = cartesian_to_spherical(&camera_position);
    let momentum = FourVector::new_spherical(1.0, 0.0, 0.0, 0.0);
    let geometry = EuclideanSpaceSpherical::new();

    render(geometry, camera_position, momentum, opts, config, filename);
}

fn render_schwarzschild(
    radius: f64,
    opts: GlobalOpts,
    config: RenderConfig,
    camera_position: Vector4<f64>,
    filename: String,
) {
    let camera_position = cartesian_to_spherical(&camera_position);

    let r = camera_position[1];
    let a = 1.0 - radius / r;
    let momentum = FourVector::new_spherical(1.0 / a.sqrt(), 0.0, 0.0, 0.0);

    let geometry = Schwarzschild::new(radius);

    render(geometry, camera_position, momentum, opts, config, filename);
}

fn render<G: Geometry>(
    geometry: G,
    camera_position: Vector4<f64>,
    camera_momentum: FourVector,
    opts: GlobalOpts,
    config: RenderConfig,
    filename: String,
) {
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
        &geometry,
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
        &geometry,
        camera,
        false,
    );

    let raytracer = raytracer::Raytracer::new(scene);
    raytracer.render(filename);
}

fn main() {
    let args = App::parse();

    let start = Instant::now();
    match args.action {
        Action::Render { filename } => {
            let config_file =
                fs::read_to_string(args.config_file).expect("Unable to open configuration file");

            let config: RenderConfig = toml::from_str(config_file.as_str()).unwrap();

            if args.global_opts.camera_position.len() != 3 {
                panic!("Camera position must be a vector of length 3");
            }
            let position = Vector4::new(
                0.0,
                args.global_opts.camera_position[0],
                args.global_opts.camera_position[1],
                args.global_opts.camera_position[2],
            );

            match config.geometry_type {
                configuration::GeometryType::Euclidean => {
                    println!("Rendering Euclidean geometry");
                    render_euclidean(args.global_opts, config, position, filename);
                }
                configuration::GeometryType::EuclideanSpherical => {
                    println!("Rendering Euclidean spherical geometry");
                    render_euclidean_spherical(args.global_opts, config, position, filename);
                }
                configuration::GeometryType::Schwarzschild { radius } => {
                    println!("Rendering Schwarzschild geometry with radius: {}", radius);
                    render_schwarzschild(radius, args.global_opts, config, position, filename);
                }
            }
        }
    }

    let duration = start.elapsed();
    println!("Elapsed time: {:.2?}", duration);
}

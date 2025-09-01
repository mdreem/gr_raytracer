mod cli;
mod configuration;
mod geometry;
mod rendering;
mod scene_objects;

use crate::cli::cli::{Action, App};
use crate::cli::euclidean::{render_euclidean, render_euclidean_ray};
use crate::cli::euclidean_spherical::{render_euclidean_spherical, render_euclidean_spherical_ray};
use crate::cli::schwarzschild::{
    render_schwarzschild, render_schwarzschild_ray, render_schwarzschild_ray_at,
};
use crate::configuration::{GeometryType, RenderConfig};
use crate::geometry::four_vector::FourVector;
use crate::geometry::point::{CoordinateSystem, Point};
use crate::rendering::raytracer::RaytracerError;
use clap::Parser;
use nalgebra::Vector4;
use std::fs;
use std::fs::File;
use std::time::Instant;

fn main() {
    run().expect("Error running raytracer");
}

fn run() -> Result<(), RaytracerError> {
    let args = App::parse();

    let start = Instant::now();
    match args.action {
        Action::Render { filename } => {
            let config_file = fs::read_to_string(args.config_file)
                .map_err(RaytracerError::ConfigurationFileError)?;

            let config: RenderConfig =
                toml::from_str(config_file.as_str()).map_err(RaytracerError::TomlError)?;

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
                GeometryType::Euclidean => {
                    println!("Rendering Euclidean geometry");
                    render_euclidean(
                        args.global_opts,
                        config,
                        Point::new_from_vector(position, CoordinateSystem::Cartesian),
                        filename,
                    )?;
                }
                GeometryType::EuclideanSpherical => {
                    println!("Rendering Euclidean spherical geometry");
                    render_euclidean_spherical(
                        args.global_opts,
                        config,
                        Point::new_from_vector(position, CoordinateSystem::Cartesian),
                        filename,
                    )?;
                }
                GeometryType::Schwarzschild {
                    radius,
                    horizon_epsilon,
                } => {
                    println!("Rendering Schwarzschild geometry with radius: {}", radius);
                    render_schwarzschild(
                        radius,
                        horizon_epsilon,
                        args.global_opts,
                        config,
                        Point::new_from_vector(position, CoordinateSystem::Cartesian),
                        filename,
                    )?;
                }
            }
        }
        Action::RenderRay { row, col, filename } => {
            let config_file = fs::read_to_string(args.config_file)
                .map_err(RaytracerError::ConfigurationFileError)?;

            let config: RenderConfig =
                toml::from_str(config_file.as_str()).map_err(RaytracerError::TomlError)?;

            if args.global_opts.camera_position.len() != 3 {
                panic!("Camera position must be a vector of length 3");
            }
            let position = Vector4::new(
                0.0,
                args.global_opts.camera_position[0],
                args.global_opts.camera_position[1],
                args.global_opts.camera_position[2],
            );

            let mut file =
                File::create(filename.clone()).map_err(RaytracerError::ConfigurationFileError)?;
            match config.geometry_type {
                GeometryType::Euclidean => {
                    println!("Rendering ray in Euclidean geometry");
                    render_euclidean_ray(
                        row,
                        col,
                        args.global_opts,
                        config,
                        Point::new_from_vector(position, CoordinateSystem::Cartesian),
                        &mut file,
                    )?;
                    println!("Saved integrated ray to {}", filename);
                }
                GeometryType::EuclideanSpherical => {
                    println!("Rendering ray in Euclidean spherical geometry");
                    render_euclidean_spherical_ray(
                        row,
                        col,
                        args.global_opts,
                        config,
                        Point::new_from_vector(position, CoordinateSystem::Cartesian),
                        &mut file,
                    )?;
                    println!("Saved integrated ray to {}", filename);
                }
                GeometryType::Schwarzschild {
                    radius,
                    horizon_epsilon,
                } => {
                    println!(
                        "Rendering ray in Schwarzschild geometry with radius: {}",
                        radius
                    );
                    render_schwarzschild_ray(
                        radius,
                        horizon_epsilon,
                        row,
                        col,
                        args.global_opts,
                        config,
                        Point::new_from_vector(position, CoordinateSystem::Cartesian),
                        &mut file,
                    )?;
                    println!("Saved integrated ray to {}", filename);
                }
            }
        }
        Action::RenderRayAt {
            position,
            direction,
            filename,
        } => {
            if position.len() != 3 {
                panic!("Position must be a vector of length 3");
            }
            if direction.len() != 3 {
                panic!("Direction must be a vector of length 3");
            }
            let mut file =
                File::create(filename.clone()).map_err(RaytracerError::ConfigurationFileError)?;
            let config_file = fs::read_to_string(args.config_file)
                .map_err(RaytracerError::ConfigurationFileError)?;
            let config: RenderConfig =
                toml::from_str(config_file.as_str()).map_err(RaytracerError::TomlError)?;
            match config.geometry_type {
                GeometryType::Euclidean => {
                    panic!("Rendering ray at not supported in Euclidean geometry yet");
                }
                GeometryType::EuclideanSpherical => {
                    panic!("Rendering ray at not supported in Euclidean spherical geometry yet");
                }
                GeometryType::Schwarzschild {
                    radius,
                    horizon_epsilon,
                } => {
                    render_schwarzschild_ray_at(
                        radius,
                        horizon_epsilon,
                        Point::new_cartesian(0.0, position[0], position[1], position[2]),
                        FourVector::new_cartesian(0.0, direction[0], direction[1], direction[2]),
                        args.global_opts,
                        &mut file,
                    )?;
                    println!("Saved integrated ray to {}", filename);
                }
            }
        }
    }

    let duration = start.elapsed();
    println!("Elapsed time: {:.2?}", duration);
    Ok(())
}

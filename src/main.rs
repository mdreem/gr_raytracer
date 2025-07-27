mod cli;
mod configuration;
mod geometry;
mod rendering;
mod scene_objects;

use crate::cli::cli::{Action, App};
use crate::cli::euclidean::{render_euclidean, render_euclidean_ray};
use crate::cli::euclidean_spherical::{render_euclidean_spherical, render_euclidean_spherical_ray};
use crate::cli::schwarzschild::{render_schwarzschild, render_schwarzschild_ray};
use crate::configuration::RenderConfig;
use crate::geometry::geometry::Geometry;
use clap::Parser;
use nalgebra::Vector4;
use std::fs;
use std::time::Instant;

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
        Action::RenderRay { row, col, filename } => {
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
                    println!("Rendering ray in Euclidean geometry");
                    render_euclidean_ray(row, col, args.global_opts, config, position, filename);
                }
                configuration::GeometryType::EuclideanSpherical => {
                    println!("Rendering ray in Euclidean spherical geometry");
                    render_euclidean_spherical_ray(
                        row,
                        col,
                        args.global_opts,
                        config,
                        position,
                        filename,
                    );
                }
                configuration::GeometryType::Schwarzschild { radius } => {
                    println!(
                        "Rendering ray in Schwarzschild geometry with radius: {}",
                        radius
                    );
                    render_schwarzschild_ray(
                        radius,
                        row,
                        col,
                        args.global_opts,
                        config,
                        position,
                        filename,
                    );
                }
            }
        }
    }

    let duration = start.elapsed();
    println!("Elapsed time: {:.2?}", duration);
}

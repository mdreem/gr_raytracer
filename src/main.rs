mod cli;
mod configuration;
mod geometry;
mod rendering;
mod scene_objects;

use crate::cli::cli::{Action, App};
use crate::configuration::RenderConfig;
use crate::geometry::four_vector::FourVector;
use crate::geometry::point::{CoordinateSystem, Point};
use crate::rendering::raytracer::RaytracerError;
use clap::Parser;
use log::info;
use nalgebra::Vector4;
use std::fs;
use std::fs::File;
use std::time::Instant;

use crate::cli::blackbody::print_blackbody_color;

fn main() {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();
    if let Err(e) = run() {
        eprintln!("Fatal error: {}", e);
        std::process::exit(1);
    }
}

fn run() -> Result<(), RaytracerError> {
    let args = App::parse();
    let start = Instant::now();

    if let Action::Blackbody {
        temperature,
        redshift,
    } = args.action
    {
        print_blackbody_color(temperature, redshift);
        return Ok(());
    }

    let config_path = args.config_file.as_ref().ok_or_else(|| {
        RaytracerError::InvalidConfiguration("Config file is required for this action".to_string())
    })?;
    if !std::path::Path::new(config_path).exists() {
        return Err(RaytracerError::InvalidConfiguration(format!(
            "Config file not found: {}",
            config_path
        )));
    }
    let config_content =
        fs::read_to_string(config_path).map_err(RaytracerError::ConfigurationFileError)?;
    let config: RenderConfig = toml::from_str(&config_content).map_err(RaytracerError::TomlError)?;

    match args.action {
        Action::Render {
            filename,
            from_row,
            from_col,
            to_row,
            to_col,
        } => {
            if args.global_opts.camera_position.len() != 3 {
                return Err(RaytracerError::InvalidConfiguration(
                    "Camera position must be a vector of length 3".to_string(),
                ));
            }
            let position = Vector4::new(
                0.0,
                args.global_opts.camera_position[0],
                args.global_opts.camera_position[1],
                args.global_opts.camera_position[2],
            );

            let geometry = config.geometry_type.get_renderable_geometry();
            info!(
                "Using coordinate system: {:?}",
                geometry.coordinate_system()
            );
            let position = Point::new_from_vector(position, CoordinateSystem::Cartesian);

            geometry.render(
                args.global_opts,
                config,
                position,
                filename,
                from_row,
                from_col,
                to_row,
                to_col,
            )?;
        }
        Action::RenderRay { row, col, filename } => {
            if args.global_opts.camera_position.len() != 3 {
                return Err(RaytracerError::InvalidConfiguration(
                    "Camera position must be a vector of length 3".to_string(),
                ));
            }
            let position = Vector4::new(
                0.0,
                args.global_opts.camera_position[0],
                args.global_opts.camera_position[1],
                args.global_opts.camera_position[2],
            );

            let mut file = File::create(filename.clone()).map_err(RaytracerError::IoError)?;

            config.geometry_type.get_renderable_geometry().render_ray(
                row,
                col,
                args.global_opts.clone(),
                config.clone(),
                Point::new_from_vector(position, CoordinateSystem::Cartesian),
                &mut file,
            )?;
            info!("Saved integrated ray to {}", filename);
        }
        Action::RenderRayAt {
            position,
            direction,
            filename,
        } => {
            if position.len() != 3 {
                return Err(RaytracerError::InvalidConfiguration(format!(
                    "Position must be a vector of length 3, got {:?}",
                    position
                )));
            }
            if direction.len() != 3 {
                return Err(RaytracerError::InvalidConfiguration(format!(
                    "Direction must be a vector of length 3, got {:?}",
                    direction
                )));
            }
            let mut file = File::create(filename.clone()).map_err(RaytracerError::IoError)?;

            config
                .geometry_type
                .get_renderable_geometry()
                .render_ray_at(
                    Point::new_cartesian(0.0, position[0], position[1], position[2]),
                    FourVector::new_cartesian(0.0, direction[0], direction[1], direction[2]),
                    args.global_opts,
                    &mut file,
                )?;
            info!("Saved integrated ray to {}", filename);
        }
        Action::Blackbody { .. } => unreachable!(),
    }

    let duration = start.elapsed();
    info!("Elapsed time: {:.2?}", duration);
    Ok(())
}

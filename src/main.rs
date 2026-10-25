mod cli;
mod configuration;
mod geometry;
mod rendering;
mod scene_objects;

use crate::cli::cli::{Action, App};
use crate::configuration::RenderConfig;
use crate::geometry::coordinates::cartesian::CartesianChart;
use crate::geometry::four_vector::FourVector;
use crate::geometry::point::Point;
use crate::rendering::error::RaytracerError;
use clap::Parser;
use log::info;
use nalgebra::Vector4;
use std::fs;
use std::fs::File;
use std::time::Instant;

use crate::cli::blackbody::{run_blackbody, run_blackbody_spectrum};

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

    match &args.action {
        Action::Blackbody {
            temperature,
            redshift,
        } => {
            run_blackbody(*temperature, *redshift);
            return Ok(());
        }
        Action::BlackbodySpectrum {
            min_temperature,
            max_temperature,
            min_redshift,
            max_redshift,
            width,
            height,
            filename,
        } => {
            run_blackbody_spectrum(
                *min_temperature,
                *max_temperature,
                *min_redshift,
                *max_redshift,
                *width,
                *height,
                filename.clone(),
                args.global_opts.tone_mapping,
            );
            return Ok(());
        }
        _ => {}
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
    let config: RenderConfig =
        toml::from_str(&config_content).map_err(RaytracerError::TomlError)?;

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
            info!("Using geometry: {:?}", config.geometry_type);
            let position = Point::new_cartesian(position[0], position[1], position[2], position[3]);

            let opts = args.global_opts.clone();
            let report = geometry.render(
                args.global_opts,
                config.clone(),
                position,
                filename.clone(),
                from_row,
                from_col,
                to_row,
                to_col,
            )?;
            let section = [from_row, from_col, to_row, to_col];
            write_render_report(&filename, &opts, section, &config, &report)?;
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
                Point::new_from_vector(position, CartesianChart::Cartesian),
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
                    FourVector::new(
                        0.0,
                        direction[0],
                        direction[1],
                        direction[2],
                        Point::new_cartesian(0.0, position[0], position[1], position[2]),
                    ),
                    args.global_opts,
                    &mut file,
                )?;
            info!("Saved integrated ray to {}", filename);
        }
        Action::Blackbody { .. } | Action::BlackbodySpectrum { .. } => unreachable!(),
    }

    let duration = start.elapsed();
    info!("Elapsed time: {:.2?}", duration);
    Ok(())
}

/// How an output was rendered, written next to it as `<filename>.render.toml`:
/// the command line, the integrator and camera settings, the accuracy check
/// if one ran, and the full scene. A `.hdr` cannot carry text metadata, so the
/// report travels beside the image instead.
#[derive(serde::Serialize)]
struct RenderInfo<'a> {
    version: &'static str,
    command: String,
    settings: RenderSettings,
    #[serde(skip_serializing_if = "Option::is_none")]
    accuracy_check: Option<&'a crate::rendering::raytracer::accuracy::AccuracyReport>,
    scene: &'a RenderConfig,
}

#[derive(serde::Serialize)]
struct RenderSettings {
    width: i64,
    height: i64,
    epsilon: f64,
    max_steps: usize,
    max_radius: f64,
    step_size: f64,
    camera_position: Vec<f64>,
    theta: f64,
    phi: f64,
    psi: f64,
    exposure: f64,
    tone_mapping: String,
    curved_star_membership: bool,
    /// from_row, from_col, to_row, to_col; omitted bounds are the full frame.
    section: Vec<i64>,
}

fn write_render_report(
    filename: &str,
    opts: &crate::cli::cli::GlobalOpts,
    section: [Option<u32>; 4],
    config: &RenderConfig,
    report: &crate::cli::shared::RenderReport,
) -> Result<(), RaytracerError> {
    let full = [0, 0, opts.height, opts.width];
    let info = RenderInfo {
        version: env!("CARGO_PKG_VERSION"),
        command: std::env::args().collect::<Vec<_>>().join(" "),
        settings: RenderSettings {
            width: opts.width,
            height: opts.height,
            epsilon: opts.epsilon,
            max_steps: opts.max_steps,
            max_radius: opts.max_radius,
            step_size: opts.step_size,
            camera_position: opts.camera_position.clone(),
            theta: opts.theta,
            phi: opts.phi,
            psi: opts.psi,
            exposure: opts.exposure,
            tone_mapping: format!("{:?}", opts.tone_mapping),
            curved_star_membership: opts.curved_star_membership,
            section: section
                .iter()
                .zip(full)
                .map(|(bound, default)| bound.map_or(default, i64::from))
                .collect(),
        },
        accuracy_check: report.accuracy.as_ref(),
        scene: config,
    };
    let text = toml::to_string(&info).map_err(|e| {
        RaytracerError::InvalidConfiguration(format!("could not serialise the render report: {e}"))
    })?;
    let path = format!("{filename}.render.toml");
    fs::write(&path, text).map_err(RaytracerError::IoError)?;
    info!("wrote render report to {}", path);
    Ok(())
}

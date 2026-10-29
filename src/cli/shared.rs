use crate::cli::cli::GlobalOpts;
use crate::configuration::{RenderConfig, TextureConfig};
use crate::geometry::atlases::Atlas;
use crate::geometry::four_vector::FourVector;
use crate::geometry::geometry::Geometry;
use crate::geometry::point::Point;
use crate::rendering::camera::Camera;
use crate::rendering::error::RaytracerError;
use crate::rendering::integration::integrator::IntegrationConfiguration;
use crate::rendering::light::color::{CIETristimulus, Color, ToneMappingMethod};
use crate::rendering::light::temperature::{ConstantTemperatureComputer, TemperatureComputer};
use crate::rendering::raytracer;
use crate::rendering::raytracer::accuracy::AccuracyReport;
use crate::rendering::scene::Scene;
use crate::rendering::stars::star_catalog::StarCatalog;
use crate::rendering::texture::{
    BlackBodyMapper, CheckerMapper, TextureData, TextureMapHandle, TextureMapperFactory,
};
use crate::scene_objects::objects::Objects;
use crate::{configuration, scene_objects};
use log::debug;
use nalgebra::Vector3;
use std::f64::consts::PI;
use std::io::Write;
use std::sync::Arc;

fn assert_directed<G: Geometry>(
    context: &str,
    geometry: &G,
    position: &Point<G::Chart>,
    four_vector: &FourVector<G::Chart>,
    expected_sign: f64,
) -> Result<(), RaytracerError> {
    let t = four_vector[0];
    let tetrad_t = geometry.get_tetrad_at(position).t;
    let orientation =
        geometry.signature()[0] * geometry.inner_product(position, &tetrad_t, four_vector);
    if !(t.is_finite() && expected_sign * orientation > 0.0) {
        return Err(RaytracerError::InvalidConfiguration(format!(
            "{} has wrong time orientation (orientation={} with t={}).",
            context, orientation, t
        )));
    }
    Ok(())
}

/// Resolve the camera four-velocity from the scene configuration.
///
/// `Explicit` is normalization-checked here so a bad vector fails as a clear
/// config error rather than an opaque downstream tetrad failure;
/// `StaticObserver`/`Zamo` are normalized by construction.
pub fn resolve_camera_velocity<G: Geometry>(
    geometry: &G,
    position: &Point<G::Chart>,
    config: &configuration::CameraVelocityConfig,
) -> Result<FourVector<G::Chart>, RaytracerError> {
    match config {
        configuration::CameraVelocityConfig::StaticObserver => {
            Ok(geometry.get_stationary_velocity_at(position))
        }
        configuration::CameraVelocityConfig::Zamo => Ok(geometry.get_zamo_velocity_at(position)),
        configuration::CameraVelocityConfig::Explicit { components } => {
            let velocity = FourVector::new(
                components[0],
                components[1],
                components[2],
                components[3],
                *position,
            );
            let norm = geometry.inner_product(position, &velocity, &velocity);
            let expected = geometry.signature()[0];
            if (norm - expected).abs() > 1e-6 {
                return Err(RaytracerError::InvalidConfiguration(format!(
                    "Explicit camera_velocity {:?} is not normalized: u.u = {} (expected {}).",
                    components, norm, expected
                )));
            }
            Ok(velocity)
        }
    }
}

pub fn assert_future_directed<G: Geometry>(
    context: &str,
    geometry: &G,
    position: &Point<G::Chart>,
    four_vector: &FourVector<G::Chart>,
) -> Result<(), RaytracerError> {
    assert_directed(context, geometry, position, four_vector, 1.0)
}

/// What a render produced besides the image, for the `<filename>.render.toml` report.
#[derive(Debug, Clone, Default)]
pub struct RenderReport {
    pub accuracy: Option<AccuracyReport>,
}

pub fn render<G: Geometry>(
    scene: Scene<G>,
    filename: String,
    tone_mapping: ToneMappingMethod,
    exposure: f64,
    verify_accuracy: usize,
    from_row: Option<u32>,
    from_col: Option<u32>,
    to_row: Option<u32>,
    to_col: Option<u32>,
) -> Result<RenderReport, RaytracerError> {
    if !(exposure.is_finite() && exposure > 0.0) {
        return Err(RaytracerError::InvalidConfiguration(format!(
            "exposure must be a positive finite number, got {exposure}"
        )));
    }
    let raytracer = raytracer::Raytracer::new(scene, tone_mapping, exposure);
    let (from_row, from_col) = (from_row.unwrap_or(0), from_col.unwrap_or(0));
    let to_row = to_row.unwrap_or(raytracer.scene.camera.rows as u32);
    let to_col = to_col.unwrap_or(raytracer.scene.camera.columns as u32);
    let luminance = raytracer.render_section(from_row, from_col, to_row, to_col, filename)?;
    let accuracy = (verify_accuracy > 0)
        .then(|| {
            raytracer.check_accuracy(
                from_row,
                from_col,
                to_row,
                to_col,
                &luminance,
                verify_accuracy,
            )
        })
        .transpose()?;
    Ok(RenderReport { accuracy })
}

/// The same point in the chart of its atlas best suited there, so that e.g. a camera on a
/// spherical chart's polar axis gets a regular tetrad.
pub fn in_best_chart<C: Atlas>(position: &Point<C>) -> Point<C> {
    position
        .chart
        .transition_point(position, position.chart.best_chart(position))
}

pub fn integrate_and_save_ray<G: Geometry>(
    geometry: &G,
    position: Point<G::Chart>,
    momentum: FourVector<G::Chart>,
    opts: GlobalOpts,
    write: &mut dyn Write,
) -> Result<(), RaytracerError> {
    let ray = crate::rendering::ray::Ray::new(0, 0, position, momentum);

    let integration_configuration = IntegrationConfiguration::new(
        opts.max_steps,
        opts.max_radius,
        opts.step_size,
        opts.epsilon,
    );
    let integrator =
        crate::rendering::integration::integrator::Integrator::new(geometry, integration_configuration);

    let (integrated_ray, stop_reason) = integrator.integrate(&ray)?;
    log::info!("Stop reason: {:?}", stop_reason);
    integrated_ray.save(write)?;
    Ok(())
}

pub fn create_scene<G: Geometry>(
    geometry: &G,
    camera_position: Point<G::Chart>,
    camera_momentum: FourVector<G::Chart>,
    opts: GlobalOpts,
    config: RenderConfig,
) -> Result<Scene<'_, G>, RaytracerError> {
    config
        .adaptive_sampling
        .validate()
        .map_err(RaytracerError::InvalidConfiguration)?;
    if let Some(star_catalog) = &config.star_catalog {
        star_catalog
            .validate()
            .map_err(RaytracerError::InvalidConfiguration)?;
    }
    let adaptive_sampling = config.adaptive_sampling.clone();
    let sampling_mask_color = opts
        .show_sampling_mask
        .then(|| CIETristimulus::from_color(&opts.sampling_mask_color));

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

    // Load the  Gaia star catalogue.
    let (
        star_catalog,
        star_flux_scale,
        winding_spread_threshold,
        max_subdivision_depth,
        curved_star_membership_config,
        magnification_cap_config,
    ) = match config.star_catalog {
        Some(catalog_config) => {
            let sky_rotation = match &catalog_config.rotation {
                Some(rotation) => Some(
                    rotation
                        .rotation()
                        .map_err(RaytracerError::InvalidConfiguration)?,
                ),
                None => None,
            };
            let catalog = StarCatalog::load_parquet_rotated(&catalog_config.path, sky_rotation)
                .map_err(|error| {
                    RaytracerError::InvalidConfiguration(format!(
                        "Failed to load star catalog {:?}: {}",
                        catalog_config.path, error
                    ))
                })?;
            debug!(
                "Loaded {} stars from {}{}",
                catalog.len(),
                catalog_config.path,
                if sky_rotation.is_some() {
                    " (celestial sphere rotated)"
                } else {
                    ""
                }
            );
            (
                Some(catalog),
                catalog_config.flux_scale,
                catalog_config.winding_spread_threshold,
                catalog_config.max_subdivision_depth,
                catalog_config.curved_star_membership,
                catalog_config.magnification_cap,
            )
        }
        None => (None, 1.0, std::f64::consts::PI, 6, false, None),
    };

    let mut objects = Objects::new(geometry);
    for object in config.objects {
        match object {
            configuration::ObjectsConfig::Sphere {
                radius,
                position,
                texture,
                temperature,
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
                    temperature,
                );
                objects.add_object(Box::new(sphere));
            }
            configuration::ObjectsConfig::Disc {
                inner_radius,
                outer_radius,
                texture,
                temperature,
                flux_scale,
            } => {
                debug!(
                    "Adding disc with inner radius: {}, outer radius: {}",
                    inner_radius, outer_radius
                );
                // Only a BlackBody-textured disc consumes the physical
                // Novikov-Thorne temperature profile; bitmap/checker discs
                // ignore temperature, and the NT calibration does not exist
                // past extremal spin (naked-singularity scenes, a > M).
                let physical_temperature = matches!(texture, TextureConfig::BlackBody { .. });
                let texture_mapper_disc = get_texture_mapper(&mut texture_mapper_factory, texture)?;
                let temperature_computer: Box<dyn TemperatureComputer> = if physical_temperature {
                    geometry.get_temperature_computer(temperature, inner_radius, outer_radius)?
                } else {
                    Box::new(ConstantTemperatureComputer::new(temperature))
                };
                let disc = scene_objects::disc::Disc::new(
                    inner_radius,
                    outer_radius,
                    texture_mapper_disc,
                    temperature_computer,
                    flux_scale,
                );
                objects.add_object(Box::new(disc));
            }
            configuration::ObjectsConfig::VolumetricDisc {
                inner_radius,
                outer_radius,
                texture,
                temperature,
                axis,
                num_octaves,
                perlin_seed,
                max_steps,
                step_size,
                thickness,
                density_multiplier,
                brightness_reference_temperature,
                absorption,
                scattering,
                noise_scale,
                noise_offset,
                flux_scale,
            } => {
                if outer_radius <= inner_radius {
                    return Err(RaytracerError::InvalidConfiguration(format!(
                        "VolumetricDisc requires outer_radius > inner_radius (got outer_radius={}, inner_radius={}).",
                        outer_radius, inner_radius
                    )));
                }
                if thickness <= 0.0 {
                    return Err(RaytracerError::InvalidConfiguration(format!(
                        "VolumetricDisc requires thickness > 0 (got thickness={}).",
                        thickness
                    )));
                }
                if max_steps == 0 {
                    return Err(RaytracerError::InvalidConfiguration(
                        "VolumetricDisc requires max_steps > 0 (got max_steps=0).".to_string(),
                    ));
                }
                if step_size <= 0.0 {
                    return Err(RaytracerError::InvalidConfiguration(format!(
                        "VolumetricDisc requires step_size > 0 (got step_size={}).",
                        step_size
                    )));
                }
                if brightness_reference_temperature <= 0.0 {
                    return Err(RaytracerError::InvalidConfiguration(format!(
                        "VolumetricDisc requires brightness_reference_temperature > 0 (got brightness_reference_temperature={}).",
                        brightness_reference_temperature
                    )));
                }
                if absorption < 0.0 {
                    return Err(RaytracerError::InvalidConfiguration(format!(
                        "VolumetricDisc requires absorption >= 0 (got absorption={}).",
                        absorption
                    )));
                }
                if scattering < 0.0 {
                    return Err(RaytracerError::InvalidConfiguration(format!(
                        "VolumetricDisc requires scattering >= 0 (got scattering={}).",
                        scattering
                    )));
                }

                debug!(
                    "Adding volumetric disc with inner radius: {}, outer radius: {}",
                    inner_radius, outer_radius
                );
                let physical_temperature = matches!(texture, TextureConfig::BlackBody { .. });
                let texture_mapper_disc = get_texture_mapper(&mut texture_mapper_factory, texture)?;
                let temperature_computer: Box<dyn TemperatureComputer> = if physical_temperature {
                    geometry.get_temperature_computer(temperature, inner_radius, outer_radius)?
                } else {
                    Box::new(ConstantTemperatureComputer::new(temperature))
                };
                let axis = axis
                    .map(|(x, y, z)| Vector3::new(x, y, z))
                    .unwrap_or(Vector3::new(0.0, 0.0, 1.0));
                let disc = scene_objects::volumetric_disc::VolumetricDisc::new(
                    inner_radius,
                    outer_radius,
                    texture_mapper_disc,
                    temperature_computer,
                    axis,
                    num_octaves,
                    perlin_seed.unwrap_or(1),
                    max_steps,
                    step_size,
                    thickness,
                    density_multiplier,
                    brightness_reference_temperature,
                    absorption,
                    scattering,
                    Vector3::new(noise_scale.0, noise_scale.1, noise_scale.2),
                    noise_offset,
                )
                .with_flux_scale(flux_scale);
                objects.add_object(Box::new(disc));
            }
        }
    }

    let mut scene = Scene::new(
        integration_configuration,
        objects,
        texture_data,
        geometry,
        camera,
        false,
        config.celestial_temperature,
    )
    .with_sampling_options(adaptive_sampling, sampling_mask_color)
    .with_star_catalog(
        star_catalog,
        star_flux_scale,
        winding_spread_threshold,
        max_subdivision_depth,
    );
    // Enabled by either the per-scene TOML flag or the CLI override.
    scene.curved_star_membership = opts.curved_star_membership || curved_star_membership_config;
    scene.star_magnification_cap = magnification_cap_config.unwrap_or(f64::INFINITY);
    Ok(scene)
}

fn get_texture_mapper(
    texture_mapper_factory: &mut TextureMapperFactory,
    texture: TextureConfig,
) -> Result<TextureMapHandle, RaytracerError> {
    let texture_mapper_sphere = match texture {
        TextureConfig::Bitmap {
            beaming_exponent,
            path,
        } => texture_mapper_factory.get_texture_mapper(beaming_exponent, path)?,
        TextureConfig::Checker {
            beaming_exponent,
            width,
            height,
            color1,
            color2,
        } => Arc::new(CheckerMapper::new(
            beaming_exponent,
            width,
            height,
            Color::new(color1.0, color1.1, color1.2, 255),
            Color::new(color2.0, color2.1, color2.2, 255),
        )),
        TextureConfig::BlackBody { beaming_exponent } => {
            Arc::new(BlackBodyMapper::new(beaming_exponent))
        }
    };
    Ok(texture_mapper_sphere)
}

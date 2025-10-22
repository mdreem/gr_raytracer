use crate::cli::cli::GlobalOpts;
use crate::cli::shared::{create_scene, render};
use crate::configuration::RenderConfig;
use crate::geometry::four_vector::FourVector;
use crate::geometry::geometry::{Geometry, InnerProduct};
use crate::geometry::kerr_analytic::KerrAnalytic;
use crate::geometry::point::Point;
use crate::geometry::spherical_coordinates_helper::cartesian_to_spherical;
use crate::rendering::integrator::{IntegrationConfiguration, Integrator};
use crate::rendering::ray::Ray;
use crate::rendering::raytracer;
use log::{debug, info};
use std::io::Write;

/// Render a full image using the analytic Kerr solution.
pub fn render_kerr_analytic(
    radius: f64,
    a: f64,
    horizon_epsilon: f64,
    opts: GlobalOpts,
    config: RenderConfig,
    camera_position: Point,
    filename: String,
) -> Result<(), raytracer::RaytracerError> {
    let camera_position = cartesian_to_spherical(&camera_position);

    let r = camera_position[1];
    let a_factor = 1.0 - radius / r;
    let momentum = FourVector::new_spherical(1.0 / a_factor, -(radius / r).sqrt(), 0.0, 0.0);
    debug!("momentum: {:?}", momentum);

    let geometry = KerrAnalytic::new(radius, a, horizon_epsilon);
    debug!(
        "m_s: {}",
        geometry.inner_product(&camera_position, &momentum, &momentum)
    );

    let scene = create_scene(&geometry, camera_position, momentum, opts, config.clone())?;

    render(scene, filename, config.color_normalization)
}

/// Render a single ray using the analytic Kerr solution.
pub fn render_kerr_analytic_ray(
    radius: f64,
    a: f64,
    horizon_epsilon: f64,
    row: i64,
    col: i64,
    opts: GlobalOpts,
    config: RenderConfig,
    camera_position: Point,
    write: &mut dyn Write,
) -> Result<(), raytracer::RaytracerError> {
    let camera_position = cartesian_to_spherical(&camera_position);
    let r = camera_position[1];
    let a_factor = 1.0 - radius / r;
    let momentum = FourVector::new_spherical(1.0 / a_factor, -(radius / r).sqrt(), 0.0, 0.0);
    let geometry = KerrAnalytic::new(radius, a, horizon_epsilon);

    let scene = create_scene(&geometry, camera_position, momentum, opts, config.clone())?;
    let raytracer = raytracer::Raytracer::new(scene, config.color_normalization);
    let (integrated_ray, stop_reason) = raytracer.integrate_ray_at_point(row, col)?;
    info!("Stop reason: {:?}", stop_reason);
    integrated_ray.save(write)?;
    Ok(())
}

/// Render a ray from a specific position and direction using the analytic Kerr solution.
pub fn render_kerr_analytic_ray_at(
    radius: f64,
    a: f64,
    horizon_epsilon: f64,
    position: Point,
    direction: FourVector,
    opts: GlobalOpts,
    write: &mut dyn Write,
) -> Result<(), raytracer::RaytracerError> {
    let position_spherical = cartesian_to_spherical(&position);
    let geometry = KerrAnalytic::new(radius, a, horizon_epsilon);

    let tetrad = geometry.get_tetrad_at(&position_spherical);
    info!("Tetrad at position {:?}: {}", position_spherical, tetrad);

    let momentum = tetrad.t * 1.0
        + tetrad.x * direction[1]
        + tetrad.y * direction[2]
        + tetrad.z * direction[3];

    let m_s = geometry.inner_product(&position_spherical, &momentum, &momentum);

    info!(
        "Momentum at position {:?}: {:?} with m_s={}",
        position_spherical, momentum, m_s
    );

    let ray = Ray::new(0, 0, position_spherical, momentum);

    let integration_configuration = IntegrationConfiguration::new(
        opts.max_steps,
        opts.max_radius,
        opts.step_size,
        opts.epsilon,
    );
    let integrator = Integrator::new(&geometry, integration_configuration);

    let (integrated_ray, stop_reason) = integrator.integrate(&ray)?;
    info!("Stop reason: {:?}", stop_reason);
    integrated_ray.save(write)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::cli::cli::GlobalOpts;
    use crate::cli::kerr_analytic::render_kerr_analytic_ray_at;
    use crate::geometry::four_vector::FourVector;
    use crate::geometry::point::Point;
    use std::f64::consts::PI;
    use std::io::BufWriter;

    #[test]
    fn test_render_kerr_analytic_ray_at() {
        let position = Point::new_cartesian(0.0, 0.0, 4.0, -18.0);
        let direction = FourVector::new_cartesian(0.0, 0.0, 1.0, 0.0);
        let radius = 1.0;
        let a = 0.5;
        let horizon_epsilon = 1e-5;
        let opts = GlobalOpts {
            width: 400,
            max_steps: 10,
            max_radius: 20.0,
            step_size: 0.01,
            epsilon: 1e-5,
            height: 400,
            phi: 0.0,
            theta: PI / 2.0,
            psi: PI / 2.0,
            camera_position: vec![],
        };
        let mut output_buffer = BufWriter::new(Vec::new());
        render_kerr_analytic_ray_at(
            radius,
            a,
            horizon_epsilon,
            position,
            direction,
            opts.clone(),
            &mut output_buffer,
        )
        .expect("Failed to render kerr analytic ray");
        let output = std::str::from_utf8(output_buffer.get_ref()).unwrap();
        let lines = output.split("\n").collect::<Vec<&str>>();

        assert_eq!(
            lines.len(),
            12,
            "Expected 12 lines in output, got {}",
            lines.len()
        );
        assert_eq!(
            lines[0], "i,t,tau,x,y,z",
            "First line does not match expected output"
        );
    }
}

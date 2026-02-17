use crate::cli::cli::GlobalOpts;
use crate::cli::shared::{create_scene, render};
use crate::configuration::RenderConfig;
use crate::geometry::four_vector::FourVector;
use crate::geometry::geometry::{Geometry, InnerProduct};
use crate::geometry::point::Point;
use crate::geometry::schwarzschild::Schwarzschild;
use crate::geometry::spherical_coordinates_helper::cartesian_to_spherical;
use crate::rendering::integrator::{IntegrationConfiguration, Integrator};
use crate::rendering::ray::Ray;
use crate::rendering::raytracer;
use crate::rendering::raytracer::RaytracerError;
use crate::rendering::scene::Scene;
use log::{debug, info};
use std::io::Write;

fn create_schwarzschild_scene<'a>(
    geometry: &'a Schwarzschild,
    radius: f64,
    opts: GlobalOpts,
    config: &RenderConfig,
    camera_position: Point,
) -> Result<Scene<'a, Schwarzschild>, RaytracerError> {
    let camera_position = cartesian_to_spherical(&camera_position);

    let r = camera_position[1];
    let a = 1.0 - radius / r;
    let momentum = FourVector::new_spherical(-1.0 / a.sqrt(), 0.0, 0.0, 0.0);
    debug!("momentum: {:?}", momentum);

    debug!(
        "m_s: {}",
        geometry.inner_product(&camera_position, &momentum, &momentum)
    );
    let scene = create_scene(geometry, camera_position, momentum, opts, config.clone())?;
    Ok(scene)
}

pub fn render_schwarzschild(
    radius: f64,
    horizon_epsilon: f64,
    opts: GlobalOpts,
    config: RenderConfig,
    camera_position: Point,
    filename: String,
    from_row: Option<u32>,
    from_col: Option<u32>,
    to_row: Option<u32>,
    to_col: Option<u32>,
) -> Result<(), RaytracerError> {
    let geometry = Schwarzschild::new(radius, horizon_epsilon);

    let scene = create_schwarzschild_scene(&geometry, radius, opts, &config, camera_position)?;

    render(
        scene,
        filename,
        config.color_normalization,
        from_row,
        from_col,
        to_row,
        to_col,
    )
}

pub fn render_schwarzschild_ray(
    radius: f64,
    horizon_epsilon: f64,
    row: i64,
    col: i64,
    opts: GlobalOpts,
    config: RenderConfig,
    camera_position: Point,
    write: &mut dyn Write,
) -> Result<(), RaytracerError> {
    let geometry = Schwarzschild::new(radius, horizon_epsilon);

    let scene = create_schwarzschild_scene(&geometry, radius, opts, &config, camera_position)?;

    let raytracer = raytracer::Raytracer::new(scene, config.color_normalization);
    let (integrated_ray, stop_reason) = raytracer.integrate_ray_at_point(row, col)?;
    info!("Stop reason: {:?}", stop_reason);
    integrated_ray.save(write)?;
    Ok(())
}

pub fn render_schwarzschild_ray_at(
    radius: f64,
    horizon_epsilon: f64,
    position: Point,
    direction: FourVector,
    opts: GlobalOpts,
    write: &mut dyn Write,
) -> Result<(), RaytracerError> {
    let position_spherical = cartesian_to_spherical(&position);
    let _r = position_spherical[1];
    let theta = position_spherical[2];
    let phi = position_spherical[3];

    let r_d = theta.sin() * phi.cos() * direction[1]
        + theta.sin() * phi.sin() * direction[2]
        + theta.cos() * direction[3];
    let theta_d = theta.cos() * phi.cos() * direction[1] + theta.cos() * phi.sin() * direction[2]
        - theta.sin() * direction[3];
    let phi_d = -phi.sin() * direction[1] + phi.cos() * direction[2];

    let geometry = Schwarzschild::new(radius, horizon_epsilon);

    let tetrad = geometry.get_tetrad_at(&position_spherical);
    info!("Tetrad at position {:?}: {}", position_spherical, tetrad);

    let momentum = tetrad.t * 1.0 + tetrad.x * (phi_d) + tetrad.y * (-theta_d) + tetrad.z * (-r_d);

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
    use crate::cli::schwarzschild::render_schwarzschild_ray_at;
    use crate::geometry::four_vector::FourVector;
    use crate::geometry::point::Point;
    use std::io::BufWriter;

    #[test]
    fn test_render_schwarzschild_ray_at() {
        let position = Point::new_cartesian(0.0, 0.0, 4.0, -18.0);
        let direction = FourVector::new_cartesian(0.0, 0.0, 1.0, 0.0);
        let radius = 1.0;
        let horizon_epsilon = 1e-5;
        let opts = GlobalOpts {
            width: 400,
            max_steps: 10,
            max_radius: 20.0,
            step_size: 0.01,
            epsilon: 1e-5,
            height: 400,
            phi: 0.0,
            theta: 0.0,
            psi: 0.0,
            camera_position: vec![],
        };
        let mut output_buffer = BufWriter::new(Vec::new());
        render_schwarzschild_ray_at(
            radius,
            horizon_epsilon,
            position,
            direction,
            opts.clone(),
            &mut output_buffer,
        )
        .expect("Failed to render schwarzschild ray");
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

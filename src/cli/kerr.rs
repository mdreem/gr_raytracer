use crate::cli::cli::GlobalOpts;
use crate::cli::shared::{create_scene, render};
use crate::configuration::RenderConfig;
use crate::geometry::four_vector::FourVector;
use crate::geometry::geometry::{Geometry, InnerProduct};
use crate::geometry::kerr::Kerr;
use crate::geometry::point::Point;
use crate::rendering::integrator::{IntegrationConfiguration, Integrator};
use crate::rendering::ray::Ray;
use crate::rendering::raytracer;
use crate::rendering::raytracer::RaytracerError;
use crate::rendering::scene::Scene;
use log::{info, trace};
use std::io::Write;

fn compute_f(geometry: &Kerr, camera_position: Point) -> f64 {
    let a = geometry.a;
    let radius = geometry.radius;
    let x = camera_position[1];
    let y = camera_position[2];
    let z = camera_position[3];
    let rho_sqr = x * x + y * y + z * z;
    let term = ((rho_sqr - a * a).powi(2) + 4.0 * a * a * z * z).sqrt();
    let r_sqr = 0.5 * (rho_sqr - a * a + term);
    let r = r_sqr.sqrt();
    let f = (r * r * r * radius) / (r * r * r * r + a * a * z * z);
    trace!("f: {}", f);
    f
}

fn create_kerr_scene<'a>(
    geometry: &'a Kerr,
    opts: GlobalOpts,
    config: &'a RenderConfig,
    camera_position: Point,
) -> Result<Scene<'a, Kerr>, RaytracerError> {
    let f = compute_f(geometry, camera_position);
    let momentum = FourVector::new_cartesian(1.0 / (1.0 - f).sqrt(), 0.0, 0.0, 0.0);
    let m_s = geometry.inner_product(&camera_position, &momentum, &momentum);
    info!(
        "Momentum at position {:?}: {:?} with m_s={}",
        camera_position, momentum, m_s
    );
    let scene = create_scene(geometry, camera_position, momentum, opts, config.clone())?;
    Ok(scene)
}

pub fn render_kerr(
    radius: f64,
    a: f64,
    horizon_epsilon: f64,
    opts: GlobalOpts,
    config: RenderConfig,
    camera_position: Point,
    filename: String,
) -> Result<(), RaytracerError> {
    let geometry = Kerr::new(radius, a, horizon_epsilon);
    let scene = create_kerr_scene(&geometry, opts, &config, camera_position)?;

    render(scene, filename, config.color_normalization)
}

pub fn render_kerr_ray(
    radius: f64,
    a: f64,
    horizon_epsilon: f64,
    row: i64,
    col: i64,
    opts: GlobalOpts,
    config: RenderConfig,
    camera_position: Point,
    write: &mut dyn Write,
) -> Result<(), RaytracerError> {
    let geometry = Kerr::new(radius, a, horizon_epsilon);
    let scene = create_kerr_scene(&geometry, opts, &config, camera_position)?;

    let raytracer = raytracer::Raytracer::new(scene, config.color_normalization);
    let (integrated_ray, stop_reason) = raytracer.integrate_ray_at_point(row, col)?;
    info!("Stop reason: {:?}", stop_reason);
    integrated_ray.save(write)?;
    Ok(())
}

pub fn render_kerr_ray_at(
    radius: f64,
    a: f64,
    horizon_epsilon: f64,
    position: Point,
    direction: FourVector,
    opts: GlobalOpts,
    write: &mut dyn Write,
) -> Result<(), RaytracerError> {
    let geometry = Kerr::new(radius, a, horizon_epsilon);

    let tetrad = geometry.get_tetrad_at(&position);
    info!("Tetrad at position {:?}: {}", position, tetrad);

    let space_part = tetrad.x * direction[1] + tetrad.y * direction[2] + tetrad.z * direction[3];
    let norm_space_part = geometry
        .inner_product(&position, &space_part, &space_part)
        .sqrt();

    let momentum = tetrad.t * 1.0
        + tetrad.x * direction[1] / norm_space_part
        + tetrad.y * direction[2] / norm_space_part
        + tetrad.z * direction[3] / norm_space_part;

    let m_s = geometry.inner_product(&position, &momentum, &momentum);

    info!(
        "Momentum at position {:?}: {:?} with m_s={}",
        position, momentum, m_s
    );

    let ray = Ray::new(0, 0, position, momentum);

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

    let final_step_position = integrated_ray.steps.last().unwrap().x;
    let final_step_momentum = integrated_ray.steps.last().unwrap().p;
    let final_m_s = geometry.inner_product(
        &final_step_position,
        &final_step_momentum,
        &final_step_momentum,
    );
    info!(
        "Final step momentum at position {:?}: {:?} with m_s={}",
        final_step_position, final_step_momentum, final_m_s
    );

    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::cli::cli::GlobalOpts;
    use crate::cli::kerr::render_kerr_ray_at;
    use crate::geometry::four_vector::FourVector;
    use crate::geometry::point::Point;
    use std::f64::consts::PI;
    use std::io::BufWriter;

    #[test]
    fn test_render_kerr_ray_at() {
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
        render_kerr_ray_at(
            radius,
            a,
            horizon_epsilon,
            position,
            direction,
            opts.clone(),
            &mut output_buffer,
        )
        .expect("Failed to render kerr ray");
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

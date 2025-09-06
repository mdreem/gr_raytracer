use crate::cli::cli::GlobalOpts;
use crate::cli::shared::{create_scene, render};
use crate::configuration::RenderConfig;
use crate::geometry::euclidean::EuclideanSpace;
use crate::geometry::euclidean_spherical::EuclideanSpaceSpherical;
use crate::geometry::four_vector::FourVector;
use crate::geometry::geometry::{Geometry, InnerProduct};
use crate::geometry::kerr::Kerr;
use crate::geometry::point::Point;
use crate::geometry::spherical_coordinates_helper::cartesian_to_spherical;
use crate::rendering::integrator::{IntegrationConfiguration, Integrator};
use crate::rendering::ray::Ray;
use crate::rendering::raytracer;
use log::info;
use std::io::Write;

pub fn render_kerr(
    radius: f64,
    a: f64,
    horizon_epsilon: f64,
    opts: GlobalOpts,
    config: RenderConfig,
    camera_position: Point,
    filename: String,
) -> Result<(), raytracer::RaytracerError> {
    let f = compute_f(radius, a, camera_position);
    let momentum = FourVector::new_cartesian(1.0, 0.0, 0.0, 0.0); // TODO: adapt to Kerr

    let geometry = Kerr::new(radius, a, horizon_epsilon);
    let scene = create_scene(&geometry, camera_position, momentum, opts, config.clone())?;

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
) -> Result<(), raytracer::RaytracerError> {
    let f = compute_f(radius, a, camera_position);
    let momentum = FourVector::new_spherical(1.0 / f.sqrt(), 0.0, 0.0, 0.0);
    let geometry = Kerr::new(radius, a, horizon_epsilon);

    let scene = create_scene(&geometry, camera_position, momentum, opts, config.clone())?;
    let raytracer = raytracer::Raytracer::new(scene, config.color_normalization);
    let (integrated_ray, stop_reason) = raytracer.integrate_ray_at_point(row, col)?;
    info!("Stop reason: {:?}", stop_reason);
    integrated_ray.save(write, &geometry)?;
    Ok(())
}

fn compute_f(radius: f64, a: f64, camera_position: Point) -> f64 {
    let x = camera_position[1];
    let y = camera_position[2];
    let z = camera_position[3];
    let rho_sqr = x * x + y * y + z * z;
    let r_sqr = 0.5
        * (rho_sqr - a * a + ((rho_sqr - a * a) * (rho_sqr - a * a) + 4.0 * a * a * z * z).sqrt());
    let r = r_sqr.sqrt();
    let f = (r * r * r * radius) / (r * r * r * r + a * a * z * z);
    f
}

pub fn render_kerr_ray_at(
    radius: f64,
    a: f64,
    horizon_epsilon: f64,
    position: Point,
    _direction: FourVector,
    opts: GlobalOpts,
    write: &mut dyn Write,
) -> Result<(), raytracer::RaytracerError> {
    let x = position[1];
    let y = position[2];
    let z = position[3];

    let bare_cartesian = EuclideanSpace::new(); // TODO: Use Kerr here.
    let geometry = Kerr::new(radius, a, horizon_epsilon);

    let tetrad = bare_cartesian.get_tetrad_at(&position);
    info!("Tetrad at position {:?}: {}", position, tetrad);

    let momentum = tetrad.t * 1.0 + tetrad.x * z + tetrad.y * y + tetrad.z * z;

    let m_s = geometry.inner_product(&position, &momentum, &momentum);

    info!(
        "Momentum at position {:?}: {:?} with m_s={}",
        position, momentum, m_s
    );

    let ray = Ray::new(0, 0, opts.width, opts.height, position, momentum);

    let integration_configuration = IntegrationConfiguration::new(
        opts.max_steps,
        opts.max_radius,
        opts.step_size,
        opts.epsilon,
    );
    let integrator = Integrator::new(&geometry, integration_configuration);

    let (integrated_ray, stop_reason) = integrator.integrate(&ray)?;
    info!("Stop reason: {:?}", stop_reason);
    integrated_ray.save(write, &geometry)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::cli::cli::GlobalOpts;
    use crate::cli::kerr::render_kerr_ray_at;
    use crate::geometry::four_vector::FourVector;
    use crate::geometry::point::Point;
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

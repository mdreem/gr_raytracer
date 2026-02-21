use crate::cli::cli::GlobalOpts;
use crate::cli::shared::{create_scene, integrate_and_save_ray, render};
use crate::configuration::RenderConfig;
use crate::geometry::four_vector::FourVector;
use crate::geometry::geometry::{Geometry, InnerProduct, RenderableGeometry};
use crate::geometry::kerr::Kerr;
use crate::geometry::point::Point;
use crate::rendering::raytracer;
use crate::rendering::raytracer::RaytracerError;
use crate::rendering::scene::Scene;
use log::info;
use std::io::Write;

fn compute_f(geometry: &Kerr, camera_position: &Point) -> f64 {
    let a = geometry.a;
    let radius = geometry.radius;
    let x = camera_position[1];
    let y = camera_position[2];
    let z = camera_position[3];
    let rho_sqr = x * x + y * y + z * z;
    let term = ((rho_sqr - a * a).powi(2) + 4.0 * a * a * z * z).sqrt();
    let r_sqr = 0.5 * (rho_sqr - a * a + term);
    let r = r_sqr.sqrt();
    (r * r * r * radius) / (r * r * r * r + a * a * z * z)
}

fn create_scene_internal<'a>(
    geometry: &'a Kerr,
    opts: GlobalOpts,
    config: &RenderConfig,
    camera_position: Point,
) -> Result<Scene<'a, Kerr>, RaytracerError> {
    let f = compute_f(geometry, &camera_position);
    let momentum = FourVector::new_cartesian(1.0 / (1.0 - f).sqrt(), 0.0, 0.0, 0.0);
    create_scene(geometry, camera_position, momentum, opts, config.clone())
}

impl RenderableGeometry for Kerr {
    fn render(
        &self,
        opts: GlobalOpts,
        config: RenderConfig,
        camera_position: Point,
        filename: String,
        from_row: Option<u32>,
        from_col: Option<u32>,
        to_row: Option<u32>,
        to_col: Option<u32>,
    ) -> Result<(), RaytracerError> {
        let scene = create_scene_internal(self, opts, &config, camera_position)?;

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

    fn render_ray(
        &self,
        row: i64,
        col: i64,
        opts: GlobalOpts,
        config: RenderConfig,
        camera_position: Point,
        write: &mut dyn Write,
    ) -> Result<(), RaytracerError> {
        let scene = create_scene_internal(self, opts, &config, camera_position)?;

        let raytracer = raytracer::Raytracer::new(scene, config.color_normalization);
        let (integrated_ray, stop_reason) = raytracer.integrate_ray_at_point(row, col)?;
        info!("Stop reason: {:?}", stop_reason);
        integrated_ray.save(write)?;
        Ok(())
    }

    fn render_ray_at(
        &self,
        position: Point,
        direction: FourVector,
        opts: GlobalOpts,
        write: &mut dyn Write,
    ) -> Result<(), RaytracerError> {
        let tetrad = self.get_tetrad_at(&position);
        info!("Tetrad at position {:?}: {}", position, tetrad);

        let space_part = tetrad.x * direction[1] + tetrad.y * direction[2] + tetrad.z * direction[3];
        let norm_space_part = self
            .inner_product(&position, &space_part, &space_part)
            .sqrt();

        let momentum = tetrad.t * 1.0
            + tetrad.x * direction[1] / norm_space_part
            + tetrad.y * direction[2] / norm_space_part
            + tetrad.z * direction[3] / norm_space_part;

        integrate_and_save_ray(self, position, momentum, opts, write)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::cli::GlobalOpts;
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
        let geometry = Kerr::new(radius, a, horizon_epsilon);
        let mut output_buffer = BufWriter::new(Vec::new());
        geometry
            .render_ray_at(position, direction, opts.clone(), &mut output_buffer)
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

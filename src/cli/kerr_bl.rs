use crate::cli::cli::GlobalOpts;
use crate::cli::shared::{assert_future_directed, create_scene, integrate_and_save_ray, render};
use crate::configuration::RenderConfig;
use crate::geometry::four_vector::FourVector;
use crate::geometry::geometry::{Geometry, InnerProduct, RenderableGeometry, SupportQuantities};
use crate::geometry::kerr_bl::KerrBL;
use crate::geometry::point::Point;
use crate::rendering::color::ToneMappingMethod;
use crate::rendering::raytracer;
use crate::rendering::raytracer::RaytracerError;
use crate::rendering::scene::Scene;
use log::info;
use std::io::Write;

fn create_scene_internal<'a>(
    geometry: &'a KerrBL,
    opts: GlobalOpts,
    config: &RenderConfig,
    camera_position: Point,
) -> Result<Scene<'a, KerrBL>, RaytracerError> {
    let bl_position = geometry.cartesian_to_bl(&camera_position);
    let momentum = geometry.get_stationary_velocity_at(&bl_position);
    assert_future_directed(
        "KerrBL camera four-velocity",
        geometry,
        &bl_position,
        &momentum,
    )?;
    create_scene(geometry, bl_position, momentum, opts, config.clone())
}

impl RenderableGeometry for KerrBL {
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
        let tone_mapping = opts.tone_mapping;
        let scene = create_scene_internal(self, opts, &config, camera_position)?;

        render(
            scene,
            filename,
            tone_mapping,
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

        let raytracer = raytracer::Raytracer::new(scene, ToneMappingMethod::default());
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
        let position_bl = self.cartesian_to_bl(&position);
        let tetrad = self.get_tetrad_at(&position_bl);
        info!("Tetrad at position {:?}: {}", position_bl, tetrad);

        let space_part =
            tetrad.x * direction[1] + tetrad.y * direction[2] + tetrad.z * direction[3];
        let norm_space_part = self
            .inner_product(&position_bl, &space_part, &space_part)
            .sqrt();

        let momentum = tetrad.t * 1.0
            + tetrad.x * direction[1] / norm_space_part
            + tetrad.y * direction[2] / norm_space_part
            + tetrad.z * direction[3] / norm_space_part;
        assert_future_directed(
            "KerrBL render_ray_at momentum",
            self,
            &position_bl,
            &momentum,
        )?;

        integrate_and_save_ray(self, position_bl, momentum, opts, write)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::cli::GlobalOpts;
    use crate::geometry::four_vector::FourVector;
    use crate::geometry::point::{CoordinateSystem, Point};
    use crate::rendering::color::ToneMappingMethod;
    use std::f64::consts::PI;
    use std::io::BufWriter;

    #[test]
    fn test_render_kerr_bl_ray_at() {
        let a = 0.5;
        let radius = 1.0;
        let horizon_epsilon = 1e-5;
        let geometry = KerrBL::new(radius, a, horizon_epsilon);

        // Provide a Cartesian position for render_ray_at (equatorial plane, r=18, phi=0)
        // cartesian_to_bl will convert this to BL (r=18, theta=PI/2, phi=0)
        let position = Point::new(0.0, 18.0, 0.0, 0.0, CoordinateSystem::Cartesian);
        // Direction: purely radial inward in BL frame
        let direction = FourVector::new(0.0, 1.0, 0.0, 0.0, CoordinateSystem::BoyerLindquist { a });

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
            tone_mapping: ToneMappingMethod::Reinhard,
        };

        let mut output_buffer = BufWriter::new(Vec::new());
        geometry
            .render_ray_at(position, direction, opts.clone(), &mut output_buffer)
            .expect("Failed to render KerrBL ray");
        let output = std::str::from_utf8(output_buffer.get_ref()).unwrap();
        let lines = output.split("\n").collect::<Vec<&str>>();

        assert_eq!(
            lines[0], "i,t,tau,x,y,z",
            "First line does not match expected output"
        );
    }
}

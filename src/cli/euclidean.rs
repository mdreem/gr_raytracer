use crate::cli::cli::GlobalOpts;
use crate::cli::shared::{create_scene, integrate_and_save_ray, render};
use crate::configuration::RenderConfig;
use crate::geometry::euclidean::EuclideanSpace;
use crate::geometry::four_vector::FourVector;
use crate::geometry::point::Point;
use crate::rendering::raytracer;
use log::debug;
use std::io::Write;
use crate::geometry::geometry::RenderableGeometry;
use crate::rendering::raytracer::RaytracerError;

impl RenderableGeometry for EuclideanSpace {
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
        let momentum = FourVector::new_cartesian(1.0, 0.0, 0.0, 0.0);
        let scene = create_scene(self, camera_position, momentum, opts, config.clone())?;

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
        let momentum = FourVector::new_cartesian(1.0, 0.0, 0.0, 0.0);

        let scene = create_scene(self, camera_position, momentum, opts, config.clone())?;
        let raytracer = raytracer::Raytracer::new(scene, config.color_normalization);
        let (integrated_ray, stop_reason) = raytracer.integrate_ray_at_point(row, col)?;
        debug!("Stop reason: {:?}", stop_reason);
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
        integrate_and_save_ray(self, position, direction, opts, write)
    }
}


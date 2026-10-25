use crate::cli::cli::GlobalOpts;
use crate::cli::shared::{
    assert_future_directed, create_scene, in_best_chart, integrate_and_save_ray, render,
    resolve_camera_velocity,
};
use crate::configuration::RenderConfig;
use crate::geometry::coordinates::euclidean::CartesianChart;
use crate::geometry::spacetimes::euclidean_spherical::EuclideanSpaceSpherical;
use crate::geometry::four_vector::FourVector;
use crate::geometry::geometry::RenderableGeometry;
use crate::geometry::point::Point;
use crate::geometry::spherical_coordinates_helper::cartesian_to_spherical;
use crate::rendering::light::color::ToneMappingMethod;
use crate::rendering::raytracer;
use crate::rendering::error::RaytracerError;
use log::debug;
use std::io::Write;

impl RenderableGeometry for EuclideanSpaceSpherical {
    fn render(
        &self,
        opts: GlobalOpts,
        config: RenderConfig,
        camera_position: Point<CartesianChart>,
        filename: String,
        from_row: Option<u32>,
        from_col: Option<u32>,
        to_row: Option<u32>,
        to_col: Option<u32>,
    ) -> Result<(), RaytracerError> {
        let camera_position = in_best_chart(&cartesian_to_spherical(&camera_position));
        let momentum = resolve_camera_velocity(self, &camera_position, &config.camera_velocity)?;
        assert_future_directed(
            "Euclidean-spherical camera four-velocity",
            self,
            &camera_position,
            &momentum,
        )?;
        let tone_mapping = opts.tone_mapping;
        let exposure = opts.exposure;
        let scene = create_scene(self, camera_position, momentum, opts, config.clone())?;

        render(
            scene,
            filename,
            tone_mapping,
            exposure,
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
        camera_position: Point<CartesianChart>,
        write: &mut dyn Write,
    ) -> Result<(), RaytracerError> {
        let camera_position = in_best_chart(&cartesian_to_spherical(&camera_position));
        let momentum = resolve_camera_velocity(self, &camera_position, &config.camera_velocity)?;
        assert_future_directed(
            "Euclidean-spherical camera four-velocity",
            self,
            &camera_position,
            &momentum,
        )?;

        let scene = create_scene(self, camera_position, momentum, opts, config.clone())?;
        let raytracer = raytracer::Raytracer::new(scene, ToneMappingMethod::default(), 1.0);
        let (integrated_ray, stop_reason) = raytracer.integrate_ray_at_point(row, col)?;
        debug!("Stop reason: {:?}", stop_reason);
        integrated_ray.save(write)?;
        Ok(())
    }

    fn render_ray_at(
        &self,
        position: Point<CartesianChart>,
        direction: FourVector<CartesianChart>,
        opts: GlobalOpts,
        write: &mut dyn Write,
    ) -> Result<(), RaytracerError> {
        let position_spherical = cartesian_to_spherical(&position);

        let r = position_spherical[1];
        let theta = position_spherical[2];
        let phi = position_spherical[3];
        if !(r.is_finite() && r > 0.0) {
            return Err(RaytracerError::InvalidConfiguration(format!(
                "Euclidean-spherical render_ray_at requires r > 0, got r={}.",
                r
            )));
        }
        let sin_theta = theta.sin();
        if !(sin_theta.is_finite() && sin_theta.abs() > 1e-12) {
            return Err(RaytracerError::InvalidConfiguration(format!(
                "Euclidean-spherical render_ray_at is undefined on the polar axis (theta={}).",
                theta
            )));
        }

        let dx = direction[1];
        let dy = direction[2];
        let dz = direction[3];
        let spatial_norm = (direction[1] * direction[1]
            + direction[2] * direction[2]
            + direction[3] * direction[3])
            .sqrt();
        if !(spatial_norm.is_finite() && spatial_norm > 0.0) {
            return Err(RaytracerError::InvalidConfiguration(
                "render_ray_at direction must have a non-zero finite spatial part.".to_string(),
            ));
        }

        // Convert Cartesian spatial direction into spherical coordinate components.
        let r_dot = sin_theta * phi.cos() * dx + sin_theta * phi.sin() * dy + theta.cos() * dz;
        let theta_unit_dot =
            theta.cos() * phi.cos() * dx + theta.cos() * phi.sin() * dy - sin_theta * dz;
        let phi_unit_dot = -phi.sin() * dx + phi.cos() * dy;
        let theta_dot = theta_unit_dot / r;
        let phi_dot = phi_unit_dot / (r * sin_theta);

        let momentum = FourVector::new(spatial_norm, r_dot, theta_dot, phi_dot, position_spherical);
        assert_future_directed(
            "Euclidean-spherical render_ray_at momentum",
            self,
            &position_spherical,
            &momentum,
        )?;

        integrate_and_save_ray(self, position_spherical, momentum, opts, write)
    }
}

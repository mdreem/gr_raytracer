use crate::cli::cli::GlobalOpts;
use crate::cli::shared::{
    assert_future_directed, create_scene, integrate_and_save_ray, render, resolve_camera_velocity,
};
use crate::configuration::RenderConfig;
use crate::geometry::coordinates::cartesian::CartesianChart;
use crate::geometry::four_vector::FourVector;
use crate::geometry::geometry::{Geometry, InnerProduct, RenderableGeometry};
use crate::geometry::point::Point;
use crate::geometry::spacetimes::kerr_bl::KerrBL;
use crate::rendering::error::RaytracerError;
use crate::rendering::light::color::ToneMappingMethod;
use crate::rendering::raytracer;
use crate::rendering::scene::Scene;
use log::info;
use std::io::Write;

fn create_scene_internal<'a>(
    geometry: &'a KerrBL,
    opts: GlobalOpts,
    config: &RenderConfig,
    camera_position: Point<CartesianChart>,
) -> Result<Scene<'a, KerrBL>, RaytracerError> {
    let bl_position = geometry.cartesian_to_bl(&camera_position);
    let momentum = resolve_camera_velocity(geometry, &bl_position, &config.camera_velocity)?;
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
        camera_position: Point<CartesianChart>,
        filename: String,
        from_row: Option<u32>,
        from_col: Option<u32>,
        to_row: Option<u32>,
        to_col: Option<u32>,
    ) -> Result<crate::cli::shared::RenderReport, RaytracerError> {
        let tone_mapping = opts.tone_mapping;
        let exposure = opts.exposure;
        let verify_accuracy = opts.verify_accuracy;
        let scene = create_scene_internal(self, opts, &config, camera_position)?;

        render(
            scene,
            filename,
            tone_mapping,
            exposure,
            verify_accuracy,
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
        let scene = create_scene_internal(self, opts, &config, camera_position)?;

        let raytracer = raytracer::Raytracer::new(scene, ToneMappingMethod::default(), 1.0);
        let (integrated_ray, stop_reason) = raytracer.integrate_ray_at_point(row, col)?;
        info!("Stop reason: {:?}", stop_reason);
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
    use crate::geometry::point::Point;
    use crate::rendering::light::color::ToneMappingMethod;
    use std::f64::consts::PI;
    use std::io::BufWriter;

    /// Global options for the README photon-ring crop: the hero camera at
    /// 7680x4320, so one pixel lands on the shadow's limb.
    fn photon_ring_opts(epsilon: f64) -> GlobalOpts {
        GlobalOpts {
            exposure: 1.0,
            width: 7680,
            height: 4320,
            max_steps: 5_000_000,
            max_radius: 15000.0,
            step_size: 0.01,
            epsilon,
            phi: 0.0,
            theta: -3.14159,
            psi: 0.0,
            camera_position: vec![-17.0, 0.0, 1.5],
            tone_mapping: ToneMappingMethod::Reinhard,
            show_sampling_mask: false,
            sampling_mask_color: crate::rendering::light::color::Color::new(255, 0, 255, 255),
            curved_star_membership: false,
            verify_accuracy: 0,
        }
    }

    /// Loops around the spin axis, escape direction and disc-plane crossing
    /// radii of the ray through `(row, col)`.
    fn photon_ring_ray(
        epsilon: f64,
        row: i64,
        col: i64,
    ) -> (f64, nalgebra::Vector3<f64>, Vec<f64>) {
        let geometry = KerrBL::new(1.0, 0.499, 1e-4);
        let config = RenderConfig {
            geometry_type: crate::configuration::GeometryType::KerrBL {
                radius: 1.0,
                a: 0.499,
                horizon_epsilon: 1e-4,
            },
            ..Default::default()
        };
        let camera = Point::new_cartesian(0.0, -17.0, 0.0, 1.5);
        let scene = create_scene_internal(&geometry, photon_ring_opts(epsilon), &config, camera)
            .expect("scene");
        let raytracer = raytracer::Raytracer::new(scene, ToneMappingMethod::default(), 1.0);
        let (ray, stop) = raytracer.integrate_ray_at_point(row, col).expect("ray");
        assert_eq!(
            stop,
            Some(crate::rendering::integration::integrator::StopReason::CelestialSphereReached)
        );

        let points: Vec<_> = ray
            .steps
            .iter()
            .map(|s| s.x.get_spatial_vector_cartesian())
            .collect();
        let mut loops = 0.0;
        let mut crossings = Vec::new();
        for pair in points.windows(2) {
            let (p, q) = (pair[0], pair[1]);
            let mut dphi = q.y.atan2(q.x) - p.y.atan2(p.x);
            if dphi > PI {
                dphi -= 2.0 * PI;
            } else if dphi < -PI {
                dphi += 2.0 * PI;
            }
            loops += dphi / (2.0 * PI);
            if p.z * q.z < 0.0 {
                crossings.push((p.x * p.x + p.y * p.y).sqrt());
            }
        }
        let escape = points.last().expect("steps").normalize();
        (loops.abs(), escape, crossings)
    }

    /// A ray on the shadow's limb loops around the photon orbit, where its
    /// deviation grows several hundredfold per loop. At the old default
    /// tolerance (1e-5) the integration error decided how many loops it made,
    /// which painted false nested arcs into the photon-ring close-up. The
    /// default must reproduce the converged path: loop count, disc crossings,
    /// and an escape direction within 0.01 rad. Pixel-exact escape at 7680x4320
    /// (~2e-4 rad) needs ~1e-11 for this ray; that is what --verify-accuracy
    /// is for.
    #[test]
    fn photon_ring_ray_converges_at_default_epsilon() {
        let (row, col) = (2046, 3674);
        let (ref_loops, ref_escape, ref_crossings) = photon_ring_ray(1e-12, row, col);
        let (loops, escape, crossings) =
            photon_ring_ray(crate::cli::cli::DEFAULT_EPSILON, row, col);

        assert!(
            (ref_loops - 3.03).abs() < 0.01,
            "reference loops {ref_loops}"
        );
        assert!(
            (loops - ref_loops).abs() < 0.01,
            "loops {loops} vs {ref_loops}"
        );
        assert_eq!(
            crossings.len(),
            ref_crossings.len(),
            "{crossings:?} vs {ref_crossings:?}"
        );
        for (r, r_ref) in crossings.iter().zip(&ref_crossings) {
            assert!((r - r_ref).abs() < 0.01, "crossing at {r} vs {r_ref}");
        }
        let deviation = escape.dot(&ref_escape).clamp(-1.0, 1.0).acos();
        assert!(deviation < 0.01, "escape direction off by {deviation} rad");

        // The test can fail: at the old default the same ray takes another path.
        let (old_loops, old_escape, _) = photon_ring_ray(1e-5, row, col);
        assert!(
            (old_loops - ref_loops).abs() > 0.5,
            "1e-5 loops {old_loops}"
        );
        assert!(old_escape.dot(&ref_escape).clamp(-1.0, 1.0).acos() > 0.5);
    }

    #[test]
    fn test_render_kerr_bl_ray_at() {
        let a = 0.5;
        let radius = 1.0;
        let horizon_epsilon = 1e-5;
        let geometry = KerrBL::new(radius, a, horizon_epsilon);

        // Provide a Cartesian position for render_ray_at (equatorial plane, r=18, phi=0)
        // cartesian_to_bl will convert this to BL (r=18, theta=PI/2, phi=0)
        let position = Point::new_cartesian(0.0, 18.0, 0.0, 0.0);
        // Direction: purely radial inward in BL frame
        let direction = FourVector::new(0.0, 1.0, 0.0, 0.0, position);

        let opts = GlobalOpts {
            exposure: 1.0,
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
            show_sampling_mask: false,
            sampling_mask_color: crate::rendering::light::color::Color::new(255, 0, 255, 255),
            curved_star_membership: false,
            verify_accuracy: 0,
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

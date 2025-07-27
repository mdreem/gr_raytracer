use crate::cli::cli::GlobalOpts;
use crate::cli::shared::{create_scene, render};
use crate::configuration::RenderConfig;
use crate::geometry::euclidean_spherical::EuclideanSpaceSpherical;
use crate::geometry::four_vector::FourVector;
use crate::geometry::spherical_coordinates_helper::cartesian_to_spherical;
use crate::rendering::raytracer;
use nalgebra::Vector4;
use std::io::Write;

pub fn render_euclidean_spherical(
    opts: GlobalOpts,
    config: RenderConfig,
    camera_position: Vector4<f64>,
    filename: String,
) {
    let camera_position = cartesian_to_spherical(&camera_position);
    let momentum = FourVector::new_spherical(1.0, 0.0, 0.0, 0.0);
    let geometry = EuclideanSpaceSpherical::new();
    let scene = create_scene(&geometry, camera_position, momentum, opts, config);

    render(scene, filename);
}

pub fn render_euclidean_spherical_ray(
    row: i64,
    col: i64,
    opts: GlobalOpts,
    config: RenderConfig,
    camera_position: Vector4<f64>,
    write: &mut dyn Write,
) {
    let camera_position = cartesian_to_spherical(&camera_position);
    let momentum = FourVector::new_spherical(1.0, 0.0, 0.0, 0.0);
    let geometry = EuclideanSpaceSpherical::new();

    let scene = create_scene(&geometry, camera_position, momentum, opts, config);
    let raytracer = raytracer::Raytracer::new(scene);
    let (integrated_ray, stop_reason) = raytracer.integrate_ray_at_point(row, col);
    println!("Stop reason: {:?}", stop_reason);
    integrated_ray.save(write, &geometry);
}

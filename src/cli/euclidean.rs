use crate::cli::cli::GlobalOpts;
use crate::cli::shared::{create_scene, render};
use crate::configuration::RenderConfig;
use crate::geometry::euclidean::EuclideanSpace;
use crate::geometry::four_vector::FourVector;
use crate::rendering::raytracer;
use nalgebra::Vector4;
use std::fs::File;
use std::io::Write;

pub fn render_euclidean(
    opts: GlobalOpts,
    config: RenderConfig,
    camera_position: Vector4<f64>,
    filename: String,
) {
    let geometry = EuclideanSpace::new();
    let momentum = FourVector::new_cartesian(1.0, 0.0, 0.0, 0.0);
    let scene = create_scene(&geometry, camera_position, momentum, opts, config);

    render(scene, filename);
}

pub fn render_euclidean_ray(
    row: i64,
    col: i64,
    opts: GlobalOpts,
    config: RenderConfig,
    camera_position: Vector4<f64>,
    write: &mut dyn Write,
) {
    let geometry = EuclideanSpace::new();
    let momentum = FourVector::new_cartesian(1.0, 0.0, 0.0, 0.0);

    let scene = create_scene(&geometry, camera_position, momentum, opts, config);
    let raytracer = raytracer::Raytracer::new(scene);
    let (integrated_ray, stop_reason) = raytracer.integrate_ray_at_point(row, col);
    println!("Stop reason: {:?}", stop_reason);
    integrated_ray.save(write, &geometry);
}

use crate::cli::cli::GlobalOpts;
use crate::cli::shared::{create_scene, render};
use crate::configuration::RenderConfig;
use crate::geometry::four_vector::FourVector;
use crate::geometry::schwarzschild::Schwarzschild;
use crate::geometry::spherical_coordinates_helper::cartesian_to_spherical;
use crate::rendering::raytracer;
use nalgebra::Vector4;

pub fn render_schwarzschild(
    radius: f64,
    opts: GlobalOpts,
    config: RenderConfig,
    camera_position: Vector4<f64>,
    filename: String,
) {
    let camera_position = cartesian_to_spherical(&camera_position);

    let r = camera_position[1];
    let a = 1.0 - radius / r;
    let momentum = FourVector::new_spherical(1.0 / a.sqrt(), 0.0, 0.0, 0.0);

    let geometry = Schwarzschild::new(radius);
    let scene = create_scene(&geometry, camera_position, momentum, opts, config);

    render(scene, filename);
}

pub fn render_schwarzschild_ray(
    radius: f64,
    row: i64,
    col: i64,
    opts: GlobalOpts,
    config: RenderConfig,
    camera_position: Vector4<f64>,
    filename: String,
) {
    let camera_position = cartesian_to_spherical(&camera_position);
    let r = camera_position[1];
    let a = 1.0 - radius / r;
    let momentum = FourVector::new_spherical(1.0 / a.sqrt(), 0.0, 0.0, 0.0);
    let geometry = Schwarzschild::new(radius);

    let scene = create_scene(&geometry, camera_position, momentum, opts, config);
    let raytracer = raytracer::Raytracer::new(scene);
    let (integrated_ray, stop_reason) = raytracer.integrate_ray_at_point(row, col);
    println!("Stop reason: {:?}", stop_reason);
    integrated_ray.save(filename.clone(), &geometry);
    println!("Saved integrated ray to {}", filename);
}

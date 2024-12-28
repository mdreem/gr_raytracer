use crate::scene::{Scene, TextureMapper};

mod camera;
mod four_vector;
mod raytracer;
mod runge_kutta;
mod scene;

use std::env;

fn main() {
    match env::current_dir() {
        Ok(path) => println!("Current directory: {}", path.display()),
        Err(e) => eprintln!("Error getting current directory: {}", e),
    }

    let texture_mapper_celestial = TextureMapper::new(String::from("./resources/celestial.png"));
    let texture_mapper_disk = TextureMapper::new(String::from("./resources/disk.png"));
    let texture_mapper_sphere = TextureMapper::new(String::from("./resources/sphere.png"));

    let scene = Scene::new(
        1000,
        0.01,
        2.0,
        3.0,
        5.0,
        texture_mapper_celestial,
        texture_mapper_disk,
        texture_mapper_sphere,
    );
    let raytracer = raytracer::Raytracer::new(500, 500, scene);
    raytracer.render();
}

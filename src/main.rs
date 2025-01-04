use crate::scene::{Scene, TextureMapper};
use nalgebra::Vector4;

mod camera;
mod euclidean;
mod euclidean_spherical;
mod four_vector;
mod geometry;
mod raytracer;
mod runge_kutta;
mod scene;
mod schwarzschild;
mod spherical_coordinates_helper;

use crate::euclidean::EuclideanSpace;
use crate::euclidean_spherical::EuclideanSpaceSpherical;
use crate::four_vector::FourVector;
use crate::schwarzschild::Schwarzschild;
use crate::spherical_coordinates_helper::cartesian_to_spherical;

fn render_euclidean() {
    let texture_mapper_celestial = TextureMapper::new(String::from("./resources/celestial.png"));
    let texture_mapper_disk = TextureMapper::new(String::from("./resources/disk.png"));
    let texture_mapper_sphere = TextureMapper::new(String::from("./resources/sphere.png"));

    let scene = Scene::new(
        10000,
        10.0,
        0.01,
        2.0,
        3.0,
        5.0,
        texture_mapper_celestial,
        texture_mapper_disk,
        texture_mapper_sphere,
        EuclideanSpace::new(),
        false,
    );
    let camera_position = Vector4::new(0.0, 0.0, 0.8, -7.0);
    let raytracer = raytracer::Raytracer::new(
        500,
        500,
        camera_position,
        FourVector::new_cartesian(1.0, 0.0, 0.0, 0.0),
        scene,
    );
    raytracer.render();
}

fn render_euclidean_spherical() {
    let texture_mapper_celestial = TextureMapper::new(String::from("./resources/celestial.png"));
    let texture_mapper_disk = TextureMapper::new(String::from("./resources/disk.png"));
    let texture_mapper_sphere = TextureMapper::new(String::from("./resources/sphere.png"));

    let scene = Scene::new(
        10000,
        10.0,
        0.01,
        2.0,
        3.0,
        5.0,
        texture_mapper_celestial,
        texture_mapper_disk,
        texture_mapper_sphere,
        EuclideanSpaceSpherical::new(),
        false,
    );
    let camera_position = cartesian_to_spherical(&Vector4::new(0.0, 0.0, 0.8, -7.0));
    let raytracer = raytracer::Raytracer::new(
        500,
        500,
        camera_position,
        FourVector::new_spherical(1.0, 0.0, 0.0, 0.0),
        scene,
    );
    raytracer.render();
}

fn render_schwarzschild() {
    let texture_mapper_celestial = TextureMapper::new(String::from("./resources/celestial.png"));
    let texture_mapper_disk = TextureMapper::new(String::from("./resources/disk.png"));
    let texture_mapper_sphere = TextureMapper::new(String::from("./resources/sphere.png"));

    let radius = 1.0;
    let scene = Scene::new(
        15000,
        20.0,
        0.01,
        0.01,
        3.0,
        5.0,
        texture_mapper_celestial,
        texture_mapper_disk,
        texture_mapper_sphere,
        Schwarzschild::new(radius),
        false,
    );
    let camera_position = cartesian_to_spherical(&Vector4::new(0.0, 0.0, 0.8, -10.0));

    let r = camera_position[1];
    let a = 1.0 - radius / r;
    let velocity = FourVector::new_spherical(1.0 / a, -(radius / r).sqrt(), 0.0, 0.0); // we have a freely falling observer here.

    let raytracer = raytracer::Raytracer::new(500, 500, camera_position, velocity, scene);
    raytracer.render();
}

fn main() {
    // render_euclidean();
    // render_euclidean_spherical();
    render_schwarzschild();
}

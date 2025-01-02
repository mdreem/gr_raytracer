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
    let raytracer = raytracer::Raytracer::new(500, 500, camera_position, scene);
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
    let camera_position_spatial = cartesian_to_spherical(&Vector4::new(0.0, 0.0, 0.8, -7.0));
    let camera_position = Vector4::new(
        0.0,
        camera_position_spatial[0],
        camera_position_spatial[1],
        camera_position_spatial[2],
    );
    let raytracer = raytracer::Raytracer::new(500, 500, camera_position, scene);
    raytracer.render();
}

fn render_schwarzschild() {
    let texture_mapper_celestial = TextureMapper::new(String::from("./resources/celestial.png"));
    let texture_mapper_disk = TextureMapper::new(String::from("./resources/disk.png"));
    let texture_mapper_sphere = TextureMapper::new(String::from("./resources/sphere.png"));

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
        Schwarzschild::new(1.0),
        false,
    );
    let camera_position_spatial = cartesian_to_spherical(&Vector4::new(0.0, 0.0, 0.8, -10.0));
    let camera_position = Vector4::new(
        0.0,
        camera_position_spatial[0],
        camera_position_spatial[1],
        camera_position_spatial[2],
    );
    let raytracer = raytracer::Raytracer::new(500, 500, camera_position, scene);
    raytracer.render();
}

fn main() {
    // render_euclidean();
    // render_euclidean_spherical();
    render_schwarzschild();
}

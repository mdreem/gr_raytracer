use crate::camera::Camera;
use crate::scene::{Scene};
use nalgebra::Vector4;
use std::f64::consts::PI;
use std::time::Instant;

mod camera;
mod color;
mod debug;
mod integrator;
mod raytracer;
mod redshift;
mod runge_kutta;
mod scene;
mod scene_objects;

mod geometry;
mod texture;

use crate::geometry::euclidean::EuclideanSpace;
use crate::geometry::euclidean_spherical::EuclideanSpaceSpherical;
use crate::geometry::four_vector::FourVector;
use crate::geometry::schwarzschild::Schwarzschild;
use crate::geometry::spherical_coordinates_helper::cartesian_to_spherical;
use crate::integrator::IntegrationConfiguration;
use crate::scene_objects::objects::Objects;
use crate::texture::{TextureData, TextureMapper};

fn render_euclidean() {
    let texture_mapper_celestial = TextureMapper::new(String::from("./resources/celestial.png"));
    let texture_mapper_disk = TextureMapper::new(String::from("./resources/disk.png"));
    let texture_mapper_sphere = TextureMapper::new(String::from("./resources/sphere.png"));
    let geometry = EuclideanSpace::new();
    let camera_position = Vector4::new(0.0, 0.0, 0.8, -15.0);

    let camera = Camera::new(
        camera_position,
        FourVector::new_cartesian(1.0, 0.0, 0.0, 0.0),
        PI / 4.0,
        500,
        500,
        &geometry,
    );

    let integration_configuration =
        IntegrationConfiguration::new(15000, 20.0, 0.01, 15000, 10000.0, 15.0);

    let texture_data = TextureData {
        celestial_map: texture_mapper_celestial,
        center_disk_map: texture_mapper_disk.clone(),
        center_sphere_map: texture_mapper_sphere.clone(),
    };

    let mut objects = Objects::new();
    let disc = scene_objects::disc::Disc::new(3.0, 5.0, texture_mapper_disk);
    let sphere = scene_objects::sphere::Sphere::new(2.0, texture_mapper_sphere);
    objects.add_object(Box::new(disc));
    objects.add_object(Box::new(sphere));

    let scene = Scene::new(
        integration_configuration,
        objects,
        texture_data,
        &geometry,
        camera,
        false,
    );

    let raytracer = raytracer::Raytracer::new(scene);
    raytracer.render();
}

fn render_euclidean_spherical() {
    let texture_mapper_celestial = TextureMapper::new(String::from("./resources/celestial.png"));
    let texture_mapper_disk = TextureMapper::new(String::from("./resources/disk.png"));
    let texture_mapper_sphere = TextureMapper::new(String::from("./resources/sphere.png"));
    let camera_position = cartesian_to_spherical(&Vector4::new(0.0, 0.0, 0.8, -7.0));
    let geometry = EuclideanSpaceSpherical::new();

    let camera = Camera::new(
        camera_position,
        FourVector::new_spherical(1.0, 0.0, 0.0, 0.0),
        PI / 4.0,
        500,
        500,
        &geometry,
    );

    let integration_configuration =
        IntegrationConfiguration::new(15000, 20.0, 0.01, 15000, 10000.0, 1.0);

    let texture_data = TextureData {
        celestial_map: texture_mapper_celestial,
        center_disk_map: texture_mapper_disk.clone(),
        center_sphere_map: texture_mapper_sphere.clone(),
    };

    let mut objects = Objects::new();
    let disc = scene_objects::disc::Disc::new(3.0, 5.0, texture_mapper_disk);
    let sphere = scene_objects::sphere::Sphere::new(2.0, texture_mapper_sphere);
    objects.add_object(Box::new(disc));
    objects.add_object(Box::new(sphere));

    let scene = Scene::new(
        integration_configuration,
        objects,
        texture_data,
        &geometry,
        camera,
        false,
    );

    let raytracer = raytracer::Raytracer::new(scene);
    raytracer.render();
}

fn render_schwarzschild() {
    let texture_mapper_celestial = TextureMapper::new(String::from("./resources/celestial.png"));
    let texture_mapper_disk = TextureMapper::new(String::from("./resources/disk.png"));
    let texture_mapper_sphere = TextureMapper::new(String::from("./resources/sphere.png"));

    let radius = 1.0;
    let camera_position = cartesian_to_spherical(&Vector4::new(0.0, 0.0, 0.8, -18.0));

    let r = camera_position[1];
    let a = 1.0 - radius / r;
    let velocity = FourVector::new_spherical(1.0 / a.sqrt(), 0.0, 0.0, 0.0); // we have a freely falling observer here.
    let geometry = Schwarzschild::new(radius);

    let camera = Camera::new(camera_position, velocity, PI / 4.0, 250, 250, &geometry);

    let integration_configuration =
        IntegrationConfiguration::new(15000, 25.0, 0.01, 15000, 10000.0, 15.0);

    let texture_data = TextureData {
        celestial_map: texture_mapper_celestial,
        center_disk_map: texture_mapper_disk.clone(),
        center_sphere_map: texture_mapper_sphere.clone(),
    };

    let mut objects = Objects::new();
    let disc = scene_objects::disc::Disc::new(3.0, 5.0, texture_mapper_disk);
    let sphere = scene_objects::sphere::Sphere::new(0.01, texture_mapper_sphere);
    objects.add_object(Box::new(disc));
    objects.add_object(Box::new(sphere));

    let scene = Scene::new(
        integration_configuration,
        objects,
        texture_data,
        &geometry,
        camera,
        false,
    );

    let raytracer = raytracer::Raytracer::new(scene);
    raytracer.render();
}

fn main() {
    let start = Instant::now();
    // render_euclidean();
    // render_euclidean_spherical();
    render_schwarzschild();

    let duration = start.elapsed();
    println!("Elapsed time: {:.2?}", duration);
}

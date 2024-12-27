use crate::scene::Scene;

mod camera;
mod four_vector;
mod raytracer;
mod runge_kutta;
mod scene;

fn main() {
    let scene = Scene::new(1000, 0.01, 2.0, 3.0, 5.0);
    let raytracer = raytracer::Raytracer::new(500, 500, scene);
    raytracer.render();
}

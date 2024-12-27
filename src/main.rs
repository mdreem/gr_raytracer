mod camera;
mod four_vector;
mod raytracer;
mod runge_kutta;
mod scene;

fn main() {
    let scene = scene::Scene::new();
    let raytracer = raytracer::Raytracer::new(500, 500, scene);
    raytracer.render();
}

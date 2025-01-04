use crate::camera::Camera;
use crate::four_vector::FourVector;
use crate::geometry::Geometry;
use crate::scene::{Scene, TextureMap};
use nalgebra::Vector4;
use rayon::iter::ParallelIterator;
use std::f64::consts::PI;
use std::sync::atomic::{AtomicUsize, Ordering};

pub struct Raytracer<T: TextureMap, G: Geometry> {
    camera: Camera<G>,
    scene: Scene<T, G>,
    image_width: i64,
    image_height: i64,
}

impl<T: TextureMap, G: Geometry> Raytracer<T, G> {
    pub fn new(
        image_width: i64,
        image_height: i64,
        camera_position: Vector4<f64>,
        velocity: FourVector,
        scene: Scene<T, G>,
    ) -> Self {
        let camera = Camera::new(
            camera_position,
            velocity,
            PI / 4.0,
            image_height,
            image_width,
            scene.geometry.clone(), // TODO see how geometry can be distributed to all needed places.
        );

        Self {
            camera,
            scene,
            image_width,
            image_height,
        }
    }

    pub fn render_ray_at(&self, row: i64, col: i64) {
        let ray = self.camera.get_ray_for(row, col);
        println!("ray: {:?}", ray);
        let color = self.scene.color_of_ray(&ray);
        println!("color: {:?}", color);
    }

    pub fn render(&self) {
        let mut imgbuf = image::ImageBuffer::new(self.image_width as u32, self.image_width as u32);
        // TODO: maybe change image width and height types to align through the codebase

        let count = AtomicUsize::new(0);
        let max_count = self.image_width * self.image_height;
        let last_progress = AtomicUsize::new(0);

        imgbuf.par_enumerate_pixels_mut().for_each(|(x, y, pixel)| {
            count.fetch_add(1, Ordering::SeqCst);

            let progress =
                (100.0 * (count.load(Ordering::Relaxed) as f64) / (max_count as f64)) as usize;
            if progress > last_progress.load(Ordering::Relaxed) {
                println!("progress: {}% ({}|{})", progress, x, y);
                last_progress.store(progress, Ordering::Relaxed);
            }

            let ray = self.camera.get_ray_for(y as i64, x as i64);
            let color = self.scene.color_of_ray(&ray);
            *pixel = image::Rgb(color.get_as_array());
        });

        imgbuf.save("render.png").unwrap();
        println!("saved image");
    }
}

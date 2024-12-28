use crate::camera::Camera;
use crate::scene::{Scene, TextureMap};
use nalgebra::Vector4;
use rayon::iter::ParallelIterator;
use std::sync::atomic::{AtomicUsize, Ordering};

pub struct Raytracer<T: TextureMap> {
    camera: Camera,
    scene: Scene<T>,
    image_width: i64,
    image_height: i64,
}

impl<T: TextureMap> Raytracer<T> {
    pub fn new(image_width: i64, image_height: i64, scene: Scene<T>) -> Self {
        let camera = Camera::new(
            Vector4::new(0.0, 0.0, 0.8, -7.0),
            std::f64::consts::PI / 4.0,
            image_height,
            image_width,
        );

        Self {
            camera,
            scene,
            image_width,
            image_height,
        }
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

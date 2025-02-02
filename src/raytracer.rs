use crate::geometry::Geometry;
use crate::scene::{Scene, TextureMap};
use rayon::iter::ParallelIterator;
use std::sync::atomic::{AtomicUsize, Ordering};

pub struct Raytracer<T: TextureMap, G: Geometry> {
    scene: Scene<T, G>,
}

impl<T: TextureMap, G: Geometry> Raytracer<T, G> {
    pub fn new(scene: Scene<T, G>) -> Self {
        Self { scene }
    }

    pub fn render_ray_at(&self, row: i64, col: i64) {
        let ray = self.scene.camera.get_ray_for(row, col);
        println!("ray: {:?}", ray);
        let color = self.scene.color_of_ray(&ray);
        println!("color: {:?}", color);
    }

    pub fn render(&self) {
        let mut imgbuf = image::ImageBuffer::new(
            self.scene.camera.columns as u32,
            self.scene.camera.rows as u32,
        );
        // TODO: maybe change image width and height types to align through the codebase

        let count = AtomicUsize::new(0);
        let max_count = self.scene.camera.rows * self.scene.camera.columns;
        let last_progress = AtomicUsize::new(0);

        imgbuf.par_enumerate_pixels_mut().for_each(|(x, y, pixel)| {
            count.fetch_add(1, Ordering::SeqCst);

            let progress =
                (100.0 * (count.load(Ordering::Relaxed) as f64) / (max_count as f64)) as usize;
            if progress > last_progress.load(Ordering::Relaxed) {
                println!("progress: {}% ({}|{})", progress, x, y);
                last_progress.store(progress, Ordering::Relaxed);
            }

            let ray = self.scene.camera.get_ray_for(y as i64, x as i64);
            let (color, _) = self.scene.color_of_ray(&ray);
            *pixel = image::Rgb(color.get_as_array());
        });

        imgbuf.save("render.png").unwrap();
        println!("saved image");
    }
}

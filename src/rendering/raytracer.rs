use crate::geometry::geometry::Geometry;
use crate::rendering::integrator::StopReason;
use crate::rendering::ray::IntegratedRay;
use crate::rendering::scene::Scene;
use crate::rendering::texture::TextureMap;
use rayon::iter::ParallelIterator;
use std::sync::atomic::{AtomicUsize, Ordering};

pub struct Raytracer<'a, T: TextureMap, G: Geometry> {
    scene: Scene<'a, T, G>,
}

impl<'a, T: TextureMap, G: Geometry> Raytracer<'a, T, G> {
    pub fn new(scene: Scene<'a, T, G>) -> Self {
        Self { scene }
    }

    pub fn render_ray_at(&self, row: i64, col: i64) {
        let ray = self.scene.camera.get_ray_for(row, col);
        println!("ray: {:?}", ray);
        let color = self.scene.color_of_ray(&ray);
        println!("color: {:?}", color);
    }

    pub fn render_section(
        &self,
        from_row: u32,
        from_col: u32,
        to_row: u32,
        to_col: u32,
        filename: String,
    ) {
        let mut imgbuf = image::ImageBuffer::new(to_col - from_col, to_row - from_row);

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

            let ray = self
                .scene
                .camera
                .get_ray_for((y + from_row) as i64, (x + from_col) as i64);
            let (color, _) = self.scene.color_of_ray(&ray);
            *pixel = image::Rgb(color.get_as_array());
        });

        imgbuf.save(&filename).unwrap();
        println!("saved image to {}", filename);
    }

    // TODO: maybe change image width and height types to align through the codebase
    pub fn render(&self, filename: String) {
        self.render_section(
            0,
            0,
            self.scene.camera.rows as u32,
            self.scene.camera.columns as u32,
            filename,
        );
    }

    pub fn integrate_ray_at_point(
        &self,
        row: i64,
        col: i64,
    ) -> (IntegratedRay, Option<StopReason>) {
        let ray = self.scene.camera.get_ray_for(row, col);
        self.scene.integrate_ray(&ray)
    }
}

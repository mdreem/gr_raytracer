use crate::geometry::geometry::Geometry;
use crate::rendering::color;
use crate::rendering::color::CIETristimulusNormalization::NoNormalization;
use crate::rendering::color::{xyz_to_srgb, CIETristimulusNormalization};
use crate::rendering::integrator::StopReason;
use crate::rendering::ray::IntegratedRay;
use crate::rendering::scene::Scene;
use image::{ImageBuffer, ImageFormat, Primitive, Rgb};
use rayon::iter::IndexedParallelIterator;
use rayon::iter::ParallelIterator;
use rayon::prelude::ParallelSliceMut;
use std::sync::atomic::{AtomicUsize, Ordering};

pub struct Raytracer<'a, G: Geometry> {
    scene: Scene<'a, G>,
    color_normalization: CIETristimulusNormalization,
}

impl<'a, G: Geometry> Raytracer<'a, G> {
    pub fn new(scene: Scene<'a, G>, color_normalization: CIETristimulusNormalization) -> Self {
        Self {
            scene,
            color_normalization,
        }
    }

    #[allow(dead_code)] // For testing
    pub fn render_ray_at(&self, row: i64, col: i64) {
        let ray = self.scene.camera.get_ray_for(row, col);
        println!("ray: {:?}", ray);
        let color = self.scene.color_of_ray(&ray);
        println!("color: {:?}", color);
    }

    pub fn render_section_to_buffer<F, B: Primitive + Sync + Send>(
        &self,
        from_row: u32,
        from_col: u32,
        to_row: u32,
        to_col: u32,
        buffer_transformer: F,
    ) -> Vec<B>
    where
        F: Fn(f32, f32, f32) -> (B, B, B) + Sync,
    {
        let count = AtomicUsize::new(0);
        let max_count = (to_row - from_row) * (to_col - from_col);
        let mut buffer: Vec<B> = vec![B::zero(); 3 * max_count as usize];

        use indicatif::{ProgressBar, ProgressStyle};
        let pb = ProgressBar::new(max_count as u64);
        pb.set_style(ProgressStyle::with_template("🎨 {spinner:.green} [{elapsed_precise}] [{wide_bar:.blue}] {pos}/{len} ({percent_precise}%, {eta})")
            .unwrap()
            .progress_chars("█▇▆▅▄▃▂▁  "));

        buffer.par_chunks_mut(3).enumerate().for_each(|(i, p)| {
            count.fetch_add(1, Ordering::SeqCst);
            pb.set_position(count.load(Ordering::Relaxed) as u64);

            let y = i as u32 / (to_col - from_col);
            let x = i as u32 % (to_col - from_col);

            let ray = self
                .scene
                .camera
                .get_ray_for((y + from_row) as i64, (x + from_col) as i64);
            let cie_tristimulus = self.scene.color_of_ray(&ray);
            (p[0], p[1], p[2]) = buffer_transformer(
                cie_tristimulus.x as f32,
                cie_tristimulus.y as f32,
                cie_tristimulus.z as f32,
            );
        });
        pb.finish();
        buffer
    }

    pub fn render_section(
        &self,
        from_row: u32,
        from_col: u32,
        to_row: u32,
        to_col: u32,
        filename: String,
    ) {
        if filename.ends_with(".hdr") {
            println!("Creating HDR image");
            let buffer =
                self.render_section_to_buffer(from_row, from_col, to_row, to_col, &|x, y, z| {
                    (x, y, z)
                });
            let imgbuf_hdr: ImageBuffer<Rgb<f32>, Vec<f32>> =
                image::ImageBuffer::from_vec(to_col - from_col, to_row - from_row, buffer).unwrap();
            imgbuf_hdr
                .save_with_format(&filename, ImageFormat::Hdr)
                .unwrap();
        } else {
            println!("Creating non-HDR image");
            let color_normalization = self.color_normalization;
            let buffer =
                self.render_section_to_buffer(from_row, from_col, to_row, to_col, |x, y, z| {
                    let cie_tristimulus = color::CIETristimulus {
                        x: x as f64,
                        y: y as f64,
                        z: z as f64,
                        alpha: 1.0,
                    };
                    let color = xyz_to_srgb(&cie_tristimulus.normalize(color_normalization), 1.0);
                    (color.r, color.g, color.b)
                });
            let imgbuf: ImageBuffer<Rgb<u8>, Vec<u8>> =
                image::ImageBuffer::from_vec(to_col - from_col, to_row - from_row, buffer).unwrap();
            imgbuf.save(&filename).unwrap();
        }

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
        println!("ray for {}-{} is: {:?}", row, col, ray);
        self.scene.integrate_ray(&ray)
    }
}

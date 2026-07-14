use crate::geometry::geometry::Geometry;
use crate::rendering::camera::CameraError;
use crate::rendering::color::{
    CIETristimulus, ToneMappingMethod, linear_srgb_to_srgb_buffer, xyz_to_linear_srgb_buffer,
};
use crate::rendering::integrator::{IntegrationError, StopReason};
use crate::rendering::ray::IntegratedRay;
use crate::rendering::scene::{RayClass, RaySample, Scene};
use crate::rendering::texture::TextureError;
use image::{ImageBuffer, ImageError, ImageFormat, Rgb};
use indicatif::style::TemplateError;
use log::{debug, error, info};
use rayon::iter::IndexedParallelIterator;
use rayon::iter::IntoParallelRefMutIterator;
use rayon::iter::ParallelIterator;
use std::io;
use std::sync::atomic::{AtomicUsize, Ordering};

#[derive(Debug, thiserror::Error)]
pub enum RaytracerError {
    #[error("Integration error: {0}")]
    IntegrationError(#[from] IntegrationError),
    #[error("I/O error: {0}")]
    IoError(#[from] io::Error),
    #[error("Configuration file error: {0}")]
    ConfigurationFileError(io::Error),
    #[error("TOML error: {0}")]
    TomlError(#[from] toml::de::Error),
    #[error("Texture error: {0}")]
    TextureError(#[from] TextureError),
    #[error("Image error: {0}")]
    ImageError(#[from] ImageError),
    #[error("Image buffer creation failed")]
    ImageBufferCreation,
    #[error("Progress bar template error: {0}")]
    ProgressBarTemplateError(#[from] TemplateError),
    #[error("Camera error: {0}")]
    CameraError(#[from] CameraError),
    #[error("Invalid configuration: {0}")]
    InvalidConfiguration(String),
    #[error("No circular orbit possible")]
    NoCircularOrbitPossible,
    #[error("Radius is below RISCO")]
    BelowRISCO,
    #[error("Radius is not finite")]
    NonFiniteRadius,
    #[error("Number is below zero")]
    NumberBelowZero,
    #[error("Denominator is close to zero")]
    DenominatorCloseToZero,
}

pub struct Raytracer<'a, G: Geometry> {
    pub scene: Scene<'a, G>,
    tone_mapping: ToneMappingMethod,
}

#[derive(PartialEq)]
struct PixelToSample {
    pub row: u32,
    pub col: u32,
    pub result: Option<CIETristimulus>,
}

impl<'a, G: Geometry> Raytracer<'a, G> {
    pub fn new(scene: Scene<'a, G>, tone_mapping: ToneMappingMethod) -> Self {
        Self {
            scene,
            tone_mapping,
        }
    }

    #[allow(dead_code)] // For testing
    pub fn render_ray_at(&self, row: i64, col: i64) {
        let ray = self.scene.camera.get_ray_for(row, col);
        debug!("ray: {:?}", ray);
        let sample = self.scene.color_of_ray(&ray);
        debug!("sample: {:?}", sample);
    }

    fn render_section_to_cie_buffer(
        &self,
        from_row: u32,
        from_col: u32,
        to_row: u32,
        to_col: u32,
    ) -> Result<Vec<CIETristimulus>, RaytracerError> {
        self.render_section_to_cie_buffer_supersampled(from_row, from_col, to_row, to_col)
    }

    fn render_section_to_cie_buffer_raw(
        &self,
        from_row: u32,
        from_col: u32,
        to_row: u32,
        to_col: u32,
    ) -> Result<Vec<RaySample>, RaytracerError> {
        let count = AtomicUsize::new(0);
        let max_count = (to_row - from_row) * (to_col - from_col);
        let mut buffer: Vec<RaySample> = vec![
            RaySample {
                color: CIETristimulus::new(0.0, 0.0, 0.0, 1.0),
                ray_class: RayClass::Escaped
            };
            max_count as usize
        ];

        use indicatif::{ProgressBar, ProgressStyle};
        let pb = ProgressBar::new(max_count as u64);
        pb.set_style(ProgressStyle::with_template("🎨 {spinner:.green} [{elapsed_precise}] [{wide_bar:.blue}] {pos}/{len} ({percent_precise}%, {eta})")
            .map_err(RaytracerError::ProgressBarTemplateError)?
            .progress_chars("█▇▆▅▄▃▂▁  "));

        buffer.par_iter_mut().enumerate().for_each(|(i, p)| {
            count.fetch_add(1, Ordering::SeqCst);
            pb.set_position(count.load(Ordering::Relaxed) as u64);

            let y = i as u32 / (to_col - from_col);
            let x = i as u32 % (to_col - from_col);

            let ray = self
                .scene
                .camera
                .get_ray_for((y + from_row) as i64, (x + from_col) as i64);

            match self.scene.color_of_ray(&ray) {
                Ok(sample) => *p = sample,
                Err(err) => {
                    error!(
                        "Unable to compute color for ray at pixel ({}, {}): {:?}",
                        x + from_col,
                        y + from_row,
                        err
                    );
                }
            }
        });
        pb.finish();
        Ok(buffer)
    }

    fn get_pixel_index(
        &self,
        row: u32,
        col: u32,
        width: u32,
        offset_row: u32,
        offset_col: u32,
    ) -> usize {
        ((row - offset_row) * width + (col - offset_col)) as usize
    }

    fn render_section_to_cie_buffer_supersampled(
        &self,
        from_row: u32,
        from_col: u32,
        to_row: u32,
        to_col: u32,
    ) -> Result<Vec<CIETristimulus>, RaytracerError> {
        info!(
            "Rendering section from ({}, {}) to ({}, {}) with supersampling",
            from_row, from_col, to_row, to_col
        );
        let buffer = self.render_section_to_cie_buffer_raw(from_row, from_col, to_row, to_col)?;
        let mut output_buffer: Vec<CIETristimulus> =
            vec![CIETristimulus::new(0.0, 0.0, 0.0, 1.0); buffer.len()];

        let debug_color_mask = false;
        let n_samples = 4;

        let mut pixels_to_sample: Vec<PixelToSample> = Vec::new();

        // Fetch candidate pixels for supersampling
        for row in from_row..to_row {
            for col in from_col..to_col {
                let pixel_index =
                    self.get_pixel_index(row, col, to_col - from_col, from_row, from_col);
                for (row_shift, col_shift) in [
                    (-1, -1),
                    (-1, 0),
                    (-1, 1),
                    (0, -1),
                    (0, 1),
                    (1, -1),
                    (1, 0),
                    (1, 1),
                ] {
                    let neighbor_row = row as i32 + row_shift;
                    let neighbor_col = col as i32 + col_shift;

                    if neighbor_row < from_row as i32
                        || neighbor_row >= to_row as i32
                        || neighbor_col < from_col as i32
                        || neighbor_col >= to_col as i32
                    {
                        continue;
                    }
                    let neighbor_index = self.get_pixel_index(
                        neighbor_row as u32,
                        neighbor_col as u32,
                        to_col - from_col,
                        from_row,
                        from_col,
                    );

                    if let (Some(pixel), Some(neighbor)) =
                        (buffer.get(pixel_index), buffer.get(neighbor_index))
                    {
                        if pixel.ray_class != neighbor.ray_class {
                            pixels_to_sample.push(PixelToSample {
                                row,
                                col,
                                result: None,
                            });
                            break; // one differing neighbor is enough to flag this pixel
                        }
                    }
                }
            }
        }

        output_buffer = buffer.into_iter().map(|s| s.color).collect();

        if debug_color_mask {
            for pixel in &pixels_to_sample {
                let pixel_index = self.get_pixel_index(
                    pixel.row,
                    pixel.col,
                    to_col - from_col,
                    from_row,
                    from_col,
                );
                output_buffer[pixel_index] = CIETristimulus::new(1.0, 0.0, 1.0, 1.0);
            }
        } else {
            info!("Supersampling {} pixels", pixels_to_sample.len());

            let count = AtomicUsize::new(0);
            use indicatif::{ProgressBar, ProgressStyle};
            let pb = ProgressBar::new(pixels_to_sample.len() as u64);
            pb.set_style(ProgressStyle::with_template("🎨 {spinner:.green} [{elapsed_precise}] [{wide_bar:.blue}] {pos}/{len} ({percent_precise}%, {eta})")
                .map_err(RaytracerError::ProgressBarTemplateError)?
                .progress_chars("█▇▆▅▄▃▂▁  "));

            pixels_to_sample
                .par_iter_mut()
                .for_each(|pixel| {
                    count.fetch_add(1, Ordering::SeqCst);
                    pb.set_position(count.load(Ordering::Relaxed) as u64);

                    let mut sample_color = CIETristimulus::new(0.0, 0.0, 0.0, 0.0);
                    let mut valid_samples = 0u32;
                    for s_row in 0..n_samples {
                        for s_col in 0..n_samples {
                            // Sample the centre of each stratum: (s + 0.5) / n.
                            // The 0.5 is where the per-pixel hash jitter will go later.
                            let ray = self.scene.camera.get_ray_for_offset(
                                pixel.row as i64,
                                pixel.col as i64,
                                (s_row as f64 + 0.5) / n_samples as f64,
                                (s_col as f64 + 0.5) / n_samples as f64,
                            );
                            match self.scene.color_of_ray(&ray) {
                                Ok(sample) => {
                                    sample_color = sample_color + sample.color;
                                    valid_samples += 1;
                                }
                                Err(err) => {
                                    error!(
                                        "Unable to compute color for ray at pixel ({}, {}): {:?}",
                                        pixel.col, pixel.row, err
                                    );
                                }
                            }
                        }
                    }
                    // Divide by the number of rays that actually returned a colour, so
                    // a failed sub-sample does not bias the pixel toward black. If all
                    // failed, leave result = None and keep the base 1-spp colour.
                    if valid_samples > 0 {
                        let inv = 1.0 / valid_samples as f64;
                        sample_color.x *= inv;
                        sample_color.y *= inv;
                        sample_color.z *= inv;
                        sample_color.alpha *= inv;
                        pixel.result = Some(sample_color);
                    }
                });
        }

        for pixel in pixels_to_sample {
            if let Some(sample_color) = pixel.result {
                let pixel_index = self.get_pixel_index(
                    pixel.row,
                    pixel.col,
                    to_col - from_col,
                    from_row,
                    from_col,
                );
                output_buffer[pixel_index] = sample_color;
            }
        }

        info!(
            "Finished rendering section from ({}, {}) to ({}, {})",
            from_row, from_col, to_row, to_col
        );
        Ok(output_buffer)
    }

    pub fn render_section(
        &self,
        from_row: u32,
        from_col: u32,
        to_row: u32,
        to_col: u32,
        filename: String,
    ) -> Result<(), RaytracerError> {
        if filename.ends_with(".hdr") {
            info!("Creating HDR image");
            let raw_cie = self.render_section_to_cie_buffer(from_row, from_col, to_row, to_col)?;
            let buffer: Vec<f32> = raw_cie
                .into_iter()
                .flat_map(|c| [c.x as f32, c.y as f32, c.z as f32])
                .collect();
            let imgbuf_hdr: ImageBuffer<Rgb<f32>, Vec<f32>> =
                image::ImageBuffer::from_vec(to_col - from_col, to_row - from_row, buffer)
                    .ok_or(RaytracerError::ImageBufferCreation)?;
            imgbuf_hdr
                .save_with_format(&filename, ImageFormat::Hdr)
                .map_err(RaytracerError::ImageError)?;
        } else {
            info!("Creating non-HDR image");
            info!("Tone mapping method: {:?}", self.tone_mapping);
            let cie_pixels =
                self.render_section_to_cie_buffer(from_row, from_col, to_row, to_col)?;
            let linear_srgb = xyz_to_linear_srgb_buffer(&cie_pixels);
            let colors = linear_srgb_to_srgb_buffer(&linear_srgb, 1.0, self.tone_mapping);
            let buffer: Vec<u8> = colors.iter().flat_map(|c| [c.r, c.g, c.b]).collect();
            let imgbuf: ImageBuffer<Rgb<u8>, Vec<u8>> =
                image::ImageBuffer::from_vec(to_col - from_col, to_row - from_row, buffer)
                    .ok_or(RaytracerError::ImageBufferCreation)?;
            imgbuf.save(&filename).map_err(RaytracerError::ImageError)?;
        }

        info!("saved image to {}", filename);
        Ok(())
    }

    pub fn integrate_ray_at_point(
        &self,
        row: i64,
        col: i64,
    ) -> Result<(IntegratedRay, Option<StopReason>), RaytracerError> {
        let ray = self.scene.camera.get_ray_for(row, col);
        info!("ray for {}-{} is: {:?}", row, col, ray);
        self.scene.integrate_ray(&ray)
    }
}

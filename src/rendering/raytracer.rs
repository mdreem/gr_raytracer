use crate::geometry::geometry::Geometry;
use crate::rendering::camera::CameraError;
use crate::rendering::color::{
    CIETristimulus, CIETristimulusNormalization, ToneMappingMethod,
    xyz_to_linear_srgb_buffer, linear_srgb_to_srgb_buffer,
};
use crate::rendering::integrator::{IntegrationError, StopReason};
use crate::rendering::ray::IntegratedRay;
use crate::rendering::scene::Scene;
use crate::rendering::texture::TextureError;
use image::{ImageBuffer, ImageError, ImageFormat, Rgb};
use indicatif::style::TemplateError;
use log::{debug, error, info};
use rayon::iter::IndexedParallelIterator;
use rayon::iter::ParallelIterator;
use rayon::iter::IntoParallelRefMutIterator;
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
    color_normalization: CIETristimulusNormalization,
    tone_mapping: ToneMappingMethod,
}

impl<'a, G: Geometry> Raytracer<'a, G> {
    pub fn new(
        scene: Scene<'a, G>,
        color_normalization: CIETristimulusNormalization,
        tone_mapping: ToneMappingMethod,
    ) -> Self {
        Self {
            scene,
            color_normalization,
            tone_mapping,
        }
    }

    #[allow(dead_code)] // For testing
    pub fn render_ray_at(&self, row: i64, col: i64) {
        let ray = self.scene.camera.get_ray_for(row, col);
        debug!("ray: {:?}", ray);
        let color = self.scene.color_of_ray(&ray);
        debug!("color: {:?}", color);
    }

    fn render_section_to_cie_buffer(
        &self,
        from_row: u32,
        from_col: u32,
        to_row: u32,
        to_col: u32,
    ) -> Result<Vec<CIETristimulus>, RaytracerError> {
        let count = AtomicUsize::new(0);
        let max_count = (to_row - from_row) * (to_col - from_col);
        let mut buffer: Vec<CIETristimulus> =
            vec![CIETristimulus::new(0.0, 0.0, 0.0, 1.0); max_count as usize];

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
                Ok(cie_tristimulus) => *p = cie_tristimulus,
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
            let raw_cie = self.render_section_to_cie_buffer(from_row, from_col, to_row, to_col)?;
            let cie_pixels: Vec<CIETristimulus> = raw_cie
                .into_iter()
                .map(|c| c.normalize(self.color_normalization))
                .collect();
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

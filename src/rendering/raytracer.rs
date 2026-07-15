use crate::configuration::AdaptiveSamplingConfig;
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

const MICHELSON_DENOMINATOR_EPSILON: f64 = 1e-4;

#[derive(PartialEq)]
struct PixelToSample {
    pub row: u32,
    pub col: u32,
    pub result: Option<CIETristimulus>,
}

/// Michelson (relative) luminance contrast; relative because HDR radiance is unbounded.
fn luminance_contrast(p: &CIETristimulus, q: &CIETristimulus) -> f64 {
    let l_p = p.y;
    let l_q = q.y;
    (l_p - l_q).abs() / (l_p + l_q + MICHELSON_DENOMINATOR_EPSILON)
}

/// Absolute opacity difference; absolute because alpha is already normalized to [0, 1].
fn opacity_contrast(p: &CIETristimulus, q: &CIETristimulus) -> f64 {
    (p.alpha - q.alpha).abs()
}

/// Faintness gate for the contrast triggers: true if the brighter of the pair
/// clears the floor (max keeps it symmetric). The floor is in linear CIE Y.
/// Scenes park the emissive disc far above unit luminance (Y > 1e4), while the
/// star-field background sits below ~1, so a floor of 1.0 keeps the disc fully
/// but stops the high-frequency background from flooding the contrast trigger
/// (it would otherwise flag ~all pixels). Class-change edges are not gated, so
/// silhouettes and the shadow rim are unaffected.
fn visible(p: &CIETristimulus, q: &CIETristimulus, minimum_luminance: f64) -> bool {
    p.y.max(q.y) > minimum_luminance
}

fn should_supersample_pair(
    pixel: &RaySample,
    neighbor: &RaySample,
    config: &AdaptiveSamplingConfig,
    minimum_luminance: f64,
) -> bool {
    if pixel.ray_class != neighbor.ray_class {
        return true;
    }

    if config.exclude_background_contrast && pixel.ray_class == RayClass::Escaped {
        return false;
    }

    visible(&pixel.color, &neighbor.color, minimum_luminance)
        && (luminance_contrast(&pixel.color, &neighbor.color) > config.luminance_contrast_threshold
            || opacity_contrast(&pixel.color, &neighbor.color) > config.opacity_contrast_threshold)
}

/// Fraction of the frame's 99th-percentile luminance used as the contrast
/// faintness floor when `minimum_luminance` is not set explicitly.
const RELATIVE_MINIMUM_LUMINANCE_FRACTION: f64 = 1e-3;

/// Resolve the contrast faintness floor: the configured absolute value, or a
/// scene-relative one. The relative floor is a small fraction of the
/// 99th-percentile luminance, which tracks the disc brightness while ignoring
/// firefly outliers (using `max` would be dominated by them).
fn resolve_minimum_luminance(config: &AdaptiveSamplingConfig, buffer: &[RaySample]) -> f64 {
    if let Some(value) = config.minimum_luminance {
        return value;
    }
    if buffer.is_empty() {
        return 0.0;
    }
    let mut luminances: Vec<f64> = buffer.iter().map(|sample| sample.color.y).collect();
    let index = (((luminances.len() - 1) as f64) * 0.99) as usize;
    luminances.select_nth_unstable_by(index, f64::total_cmp);
    RELATIVE_MINIMUM_LUMINANCE_FRACTION * luminances[index]
}

// https://rosettacode.org/wiki/Pseudo-random_numbers/Splitmix64
fn mix64(mut z: u64) -> u64 {
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
    z ^ (z >> 31)
}

fn hash_pixel_samples(row: i64, col: i64, k: usize) -> f64 {
    let z = mix64((row as u64).wrapping_add(mix64((col as u64).wrapping_add(mix64(k as u64)))));

    // Convert to a floating-point number in the range [0, 1)
    (z >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
}

fn stratified_sample_offset(
    row: i64,
    col: i64,
    stratum_row: usize,
    stratum_col: usize,
    samples_per_axis: usize,
) -> (f64, f64) {
    debug_assert!(samples_per_axis > 0);
    let sample_index = stratum_row * samples_per_axis + stratum_col;
    let dx = (stratum_col as f64 + hash_pixel_samples(row, col, 2 * sample_index))
        / samples_per_axis as f64;
    let dy = (stratum_row as f64 + hash_pixel_samples(row, col, 2 * sample_index + 1))
        / samples_per_axis as f64;
    (dx, dy)
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
        if self.scene.adaptive_sampling.enabled || self.scene.sampling_mask_color.is_some() {
            self.render_section_to_cie_buffer_supersampled(from_row, from_col, to_row, to_col)
        } else {
            Ok(self
                .render_section_to_cie_buffer_raw(from_row, from_col, to_row, to_col)?
                .into_iter()
                .map(|sample| sample.color)
                .collect())
        }
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
        let samples_per_axis = self.scene.adaptive_sampling.samples_per_axis;
        let minimum_luminance =
            resolve_minimum_luminance(&self.scene.adaptive_sampling, &buffer);

        let mut pixels_to_sample = self.collect_pixels_to_supersample(
            from_row,
            from_col,
            to_row,
            to_col,
            &buffer,
            minimum_luminance,
        );

        let mut output_buffer: Vec<CIETristimulus> =
            buffer.into_iter().map(|sample| sample.color).collect();

        if let Some(mask_color) = self.scene.sampling_mask_color {
            for pixel in &pixels_to_sample {
                let pixel_index = self.get_pixel_index(
                    pixel.row,
                    pixel.col,
                    to_col - from_col,
                    from_row,
                    from_col,
                );
                output_buffer[pixel_index] = mask_color;
            }
        } else {
            self.supersample(samples_per_axis, &mut pixels_to_sample)?;
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

    fn supersample(
        &self,
        samples_per_axis: usize,
        pixels_to_sample: &mut Vec<PixelToSample>,
    ) -> Result<(), RaytracerError> {
        info!("Supersampling {} pixels", pixels_to_sample.len());

        let count = AtomicUsize::new(0);
        use indicatif::{ProgressBar, ProgressStyle};
        let pb = ProgressBar::new(pixels_to_sample.len() as u64);
        pb.set_style(ProgressStyle::with_template("🎨 {spinner:.green} [{elapsed_precise}] [{wide_bar:.blue}] {pos}/{len} ({percent_precise}%, {eta})")
            .map_err(RaytracerError::ProgressBarTemplateError)?
            .progress_chars("█▇▆▅▄▃▂▁  "));

        pixels_to_sample.par_iter_mut().for_each(|pixel| {
            count.fetch_add(1, Ordering::SeqCst);
            pb.set_position(count.load(Ordering::Relaxed) as u64);

            let mut sample_color = CIETristimulus::new(0.0, 0.0, 0.0, 0.0);
            let mut valid_samples = 0u32;
            for stratum_row in 0..samples_per_axis {
                for stratum_col in 0..samples_per_axis {
                    let (dx, dy) = stratified_sample_offset(
                        pixel.row as i64,
                        pixel.col as i64,
                        stratum_row,
                        stratum_col,
                        samples_per_axis,
                    );
                    let ray = self.scene.camera.get_ray_for_offset(
                        pixel.row as i64,
                        pixel.col as i64,
                        dx,
                        dy,
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
        pb.finish();

        Ok(())
    }

    fn collect_pixels_to_supersample(
        &self,
        from_row: u32,
        from_col: u32,
        to_row: u32,
        to_col: u32,
        buffer: &[RaySample],
        minimum_luminance: f64,
    ) -> Vec<PixelToSample> {
        let mut pixels_to_sample: Vec<PixelToSample> = Vec::new();
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

                    // Neighbours outside this section are skipped. For a full
                    // render that is only the true image border, and for a
                    // standalone section crop (previewing a region) it just
                    // leaves the crop's 1px border compared within the crop,
                    // which is cosmetically negligible. It would only matter if
                    // separately-rendered sections were stitched into one image,
                    // where an edge crossing a seam could stay at 1 spp and show
                    // a line; making that seam-free would need a 1px selection
                    // halo around the section. Not needed for crop previews.
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
                        && should_supersample_pair(
                            pixel,
                            neighbor,
                            &self.scene.adaptive_sampling,
                            minimum_luminance,
                        )
                    {
                        pixels_to_sample.push(PixelToSample {
                            row,
                            col,
                            result: None,
                        });
                        break;
                    }
                }
            }
        }

        pixels_to_sample
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

#[cfg(test)]
mod tests {
    use super::{
        MICHELSON_DENOMINATOR_EPSILON, luminance_contrast, should_supersample_pair,
        stratified_sample_offset,
    };
    use crate::configuration::AdaptiveSamplingConfig;
    use crate::rendering::color::CIETristimulus;
    use crate::rendering::scene::{RayClass, RaySample};

    fn sample(y: f64, alpha: f64, ray_class: RayClass) -> RaySample {
        RaySample {
            color: CIETristimulus::new(0.0, y, 0.0, alpha),
            ray_class,
        }
    }

    #[test]
    fn stratified_offsets_stay_in_their_cells() {
        let samples_per_axis = 4;
        for stratum_row in 0..samples_per_axis {
            for stratum_col in 0..samples_per_axis {
                let (dx, dy) =
                    stratified_sample_offset(17, 23, stratum_row, stratum_col, samples_per_axis);
                let cell_size = 1.0 / samples_per_axis as f64;
                assert!(
                    (stratum_col as f64 * cell_size..(stratum_col + 1) as f64 * cell_size)
                        .contains(&dx)
                );
                assert!(
                    (stratum_row as f64 * cell_size..(stratum_row + 1) as f64 * cell_size)
                        .contains(&dy)
                );
                assert_eq!(
                    (dx, dy),
                    stratified_sample_offset(17, 23, stratum_row, stratum_col, samples_per_axis,)
                );
            }
        }
    }

    #[test]
    fn michelson_contrast_uses_the_named_epsilon() {
        let black = CIETristimulus::new(0.0, 0.0, 0.0, 1.0);
        let faint = CIETristimulus::new(0.0, MICHELSON_DENOMINATOR_EPSILON, 0.0, 1.0);
        assert_eq!(luminance_contrast(&black, &black), 0.0);
        assert_eq!(luminance_contrast(&black, &faint), 0.5);
    }

    #[test]
    fn class_boundaries_are_always_supersampled() {
        let config = AdaptiveSamplingConfig::default();
        let escaped = sample(0.0, 1.0, RayClass::Escaped);
        let captured = sample(0.0, 1.0, RayClass::Captured);

        assert!(should_supersample_pair(&escaped, &captured, &config, 100.0));
        assert!(should_supersample_pair(&captured, &escaped, &config, 100.0));
    }

    #[test]
    fn background_contrast_does_not_trigger_supersampling() {
        let config = AdaptiveSamplingConfig {
            luminance_contrast_threshold: 0.0,
            opacity_contrast_threshold: 0.0,
            ..Default::default()
        };
        let dark = sample(1.0, 0.0, RayClass::Escaped);
        let bright = sample(100.0, 1.0, RayClass::Escaped);

        assert!(!should_supersample_pair(&dark, &bright, &config, 0.0));
    }

    #[test]
    fn visible_object_contrast_triggers_supersampling() {
        let config = AdaptiveSamplingConfig {
            luminance_contrast_threshold: 0.2,
            opacity_contrast_threshold: 0.2,
            ..Default::default()
        };

        assert!(should_supersample_pair(
            &sample(2.0, 1.0, RayClass::Hit),
            &sample(1.0, 1.0, RayClass::Hit),
            &config,
            1.0,
        ));
        assert!(should_supersample_pair(
            &sample(2.0, 0.6, RayClass::Hit),
            &sample(2.0, 0.9, RayClass::Hit),
            &config,
            1.0,
        ));
    }

    #[test]
    fn faint_object_contrast_does_not_trigger_supersampling() {
        let config = AdaptiveSamplingConfig {
            luminance_contrast_threshold: 0.0,
            opacity_contrast_threshold: 0.0,
            ..Default::default()
        };

        assert!(!should_supersample_pair(
            &sample(1.0, 0.0, RayClass::Hit),
            &sample(0.0, 1.0, RayClass::Hit),
            &config,
            1.0,
        ));
    }
}

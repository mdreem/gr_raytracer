//! Initial ray samples, adaptive selection, and supersampling.

use super::progress::Progress;
use super::Raytracer;
use crate::configuration::AdaptiveSamplingConfig;
use crate::geometry::coordinates::chart::Chart;
use crate::geometry::four_vector::FourVector;
use crate::geometry::geometry::Geometry;
use crate::rendering::error::RaytracerError;
use crate::rendering::integration::integrator::StopReason;
use crate::rendering::light::radiance::{Radiance, RadianceMean};
use crate::rendering::ray::Ray;
use crate::rendering::scene::{EscapeInfo, RayClass, RaySample};
use log::{error, info};
use rayon::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering};

const MICHELSON_DENOMINATOR_EPSILON: f64 = 1e-4;

#[derive(PartialEq)]
pub(super) struct PixelToSample {
    pub row: u32,
    pub col: u32,
    pub result: Option<Radiance>,
}

/// Michelson (relative) luminance contrast; relative because HDR radiance is unbounded.
fn luminance_contrast(p: &Radiance, q: &Radiance) -> f64 {
    let l_p = p.y;
    let l_q = q.y;
    (l_p - l_q).abs() / (l_p + l_q + MICHELSON_DENOMINATOR_EPSILON)
}

/// Absolute opacity difference equals the absolute transmittance difference.
fn opacity_contrast(p: &Radiance, q: &Radiance) -> f64 {
    (p.transmittance - q.transmittance).abs()
}

/// Faintness gate for the contrast triggers: true if the brighter of the pair
/// clears the floor (max keeps it symmetric). The floor is in linear CIE Y.
/// Scenes park the emissive disc far above unit luminance (Y > 1e4), while the
/// star-field background sits below ~1, so a floor of 1.0 keeps the disc fully
/// but stops the high-frequency background from flooding the contrast trigger
/// (it would otherwise flag ~all pixels). Class-change edges are not gated, so
/// silhouettes and the shadow rim are unaffected.
fn visible(p: &Radiance, q: &Radiance, minimum_luminance: f64) -> bool {
    p.y.max(q.y) > minimum_luminance
}

fn should_supersample_pair<C: Chart>(
    pixel: &RaySample<C>,
    neighbor: &RaySample<C>,
    config: &AdaptiveSamplingConfig,
    minimum_luminance: f64,
) -> bool {
    if std::mem::discriminant(&pixel.ray_class) != std::mem::discriminant(&neighbor.ray_class) {
        return true;
    }

    if config.exclude_background_contrast && matches!(pixel.ray_class, RayClass::Escaped(_)) {
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
fn resolve_minimum_luminance<C: Chart>(
    config: &AdaptiveSamplingConfig,
    buffer: &[RaySample<C>],
) -> f64 {
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
    /// The section's pixel colours, row-major, and why each pixel's primary
    /// (centre) ray stopped.
    pub(super) fn render_section_to_radiance_buffer(
        &self,
        from_row: u32,
        from_col: u32,
        to_row: u32,
        to_col: u32,
    ) -> Result<(Vec<Radiance>, Vec<Option<StopReason>>), RaytracerError> {
        let buffer =
            self.render_section_to_radiance_buffer_raw(from_row, from_col, to_row, to_col)?;
        let stop_reasons = buffer.iter().map(|s| s.stop_reason).collect();
        let colors =
            if self.scene.adaptive_sampling.enabled || self.scene.sampling_mask_color.is_some() {
                self.render_section_to_radiance_buffer_supersampled(
                    from_row, from_col, to_row, to_col, &buffer,
                )?
            } else {
                let width = (to_col - from_col) as usize;
                let height = (to_row - from_row) as usize;
                let mut colors: Vec<Radiance> = buffer.iter().map(|s| s.color).collect();
                self.add_star_layer(&buffer, width, height, &mut colors)?;
                colors
            };
        Ok((colors, stop_reasons))
    }

    pub(super) fn render_section_to_radiance_buffer_raw(
        &self,
        from_row: u32,
        from_col: u32,
        to_row: u32,
        to_col: u32,
    ) -> Result<Vec<RaySample<G::Chart>>, RaytracerError> {
        let count = AtomicUsize::new(0);
        let max_count = (to_row - from_row) * (to_col - from_col);
        // Create some dummy rays. They will be replaced with the actual rays in the parallel loop below.
        let mut buffer: Vec<RaySample<G::Chart>> = vec![
            RaySample {
                color: Radiance::BLACK,
                ray_class: RayClass::Escaped(EscapeInfo {
                    x: 0.0,
                    y: 0.0,
                    z: 0.0,
                    redshift: 1.0,
                }),
                stop_reason: None,
                accumulated_angular_distance: 0.0,
                ray: Ray::new(
                    0,
                    0,
                    self.scene.camera.position,
                    FourVector::new(0.0, 0.0, 1.0, 0.0, self.scene.camera.position),
                ),
            };
            max_count as usize
        ];

        let progress = Progress::new(max_count as u64, "render")?;

        buffer.par_iter_mut().enumerate().for_each(|(i, p)| {
            let done = count.fetch_add(1, Ordering::SeqCst) + 1;
            progress.set(done as u64);

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
        progress.finish();
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

    fn render_section_to_radiance_buffer_supersampled(
        &self,
        from_row: u32,
        from_col: u32,
        to_row: u32,
        to_col: u32,
        buffer: &[RaySample<G::Chart>],
    ) -> Result<Vec<Radiance>, RaytracerError> {
        info!(
            "Rendering section from ({}, {}) to ({}, {}) with supersampling",
            from_row, from_col, to_row, to_col
        );
        let samples_per_axis = self.scene.adaptive_sampling.samples_per_axis;
        let minimum_luminance = resolve_minimum_luminance(&self.scene.adaptive_sampling, buffer);

        let mut pixels_to_sample = self.collect_pixels_to_supersample(
            from_row,
            from_col,
            to_row,
            to_col,
            buffer,
            minimum_luminance,
        );

        // Initialize the output buffer with the base 1-spp colors from the raw buffer.
        let mut output_buffer: Vec<Radiance> = buffer.iter().map(|sample| sample.color).collect();

        if let Some(mask_color) = self.scene.sampling_mask_color {
            // The mask is a fixed diagnostic color; pre-divide by the
            // exposure so the later exposure multiply restores it exactly
            // and the overlay looks identical at every --exposure.
            let mask_color = mask_color.mul_color_part(1.0 / self.exposure);
            for pixel in &pixels_to_sample {
                let pixel_index = self.get_pixel_index(
                    pixel.row,
                    pixel.col,
                    to_col - from_col,
                    from_row,
                    from_col,
                );
                output_buffer[pixel_index] = Radiance::from_straight(mask_color);
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

        let width = (to_col - from_col) as usize;
        let height = (to_row - from_row) as usize;
        self.add_star_layer(buffer, width, height, &mut output_buffer)?;

        info!(
            "Finished rendering section from ({}, {}) to ({}, {})",
            from_row, from_col, to_row, to_col
        );
        Ok(output_buffer)
    }

    pub(super) fn supersample(
        &self,
        samples_per_axis: usize,
        pixels_to_sample: &mut Vec<PixelToSample>,
    ) -> Result<(), RaytracerError> {
        info!("Supersampling {} pixels", pixels_to_sample.len());

        let count = AtomicUsize::new(0);
        let progress = Progress::new(pixels_to_sample.len() as u64, "supersample")?;

        pixels_to_sample.par_iter_mut().for_each(|pixel| {
            let done = count.fetch_add(1, Ordering::SeqCst) + 1;
            progress.set(done as u64);

            let mut sample_colors = RadianceMean::default();
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
                            sample_colors.add(sample.color, 1.0);
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
            pixel.result = sample_colors.mean();
        });
        progress.finish();

        Ok(())
    }

    pub(super) fn collect_pixels_to_supersample(
        &self,
        from_row: u32,
        from_col: u32,
        to_row: u32,
        to_col: u32,
        buffer: &[RaySample<G::Chart>],
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
}

#[cfg(test)]
mod tests {
    use super::{
        MICHELSON_DENOMINATOR_EPSILON, luminance_contrast, should_supersample_pair,
        stratified_sample_offset,
    };
    use crate::configuration::AdaptiveSamplingConfig;
    use crate::geometry::coordinates::cartesian::CartesianChart;
    use crate::geometry::four_vector::FourVector;
    use crate::geometry::point::Point;
    use crate::rendering::light::radiance::Radiance;
    use crate::rendering::ray::Ray;
    use crate::rendering::scene::{EscapeInfo, RayClass, RaySample};

    fn sample(y: f64, alpha: f64, ray_class: RayClass) -> RaySample<CartesianChart> {
        RaySample {
            color: Radiance::new(0.0, y, 0.0, 1.0 - alpha),
            ray_class,
            stop_reason: None,
            accumulated_angular_distance: 0.0,
            ray: Ray::new(
                0,
                0,
                Point::new_cartesian(0.0, 0.0, 0.0, 0.0),
                FourVector::new(0.0, 0.0, 1.0, 0.0, Point::new_cartesian(0.0, 0.0, 0.0, 0.0)),
            ),
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
        let black = Radiance::BLACK;
        let faint = Radiance::new(0.0, MICHELSON_DENOMINATOR_EPSILON, 0.0, 0.0);
        assert_eq!(luminance_contrast(&black, &black), 0.0);
        assert_eq!(luminance_contrast(&black, &faint), 0.5);
    }

    #[test]
    fn class_boundaries_are_always_supersampled() {
        let config = AdaptiveSamplingConfig::default();
        let escaped = sample(
            0.0,
            1.0,
            RayClass::Escaped(EscapeInfo {
                x: 0.0,
                y: 0.0,
                z: 1.0,
                redshift: 1.0,
            }),
        );
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
        let dark = sample(
            1.0,
            0.0,
            RayClass::Escaped(EscapeInfo {
                x: 0.0,
                y: 0.0,
                z: 1.0,
                redshift: 1.0,
            }),
        );
        let bright = sample(
            100.0,
            1.0,
            RayClass::Escaped(EscapeInfo {
                x: 1.0,
                y: 0.0,
                z: 0.0,
                redshift: 1.0,
            }),
        );

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
